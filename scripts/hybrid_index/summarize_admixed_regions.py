#!/usr/bin/env python3
"""
summarize_admixed_regions.py

Summarize admixed / introgressed genomic regions from sliding-window hybrid-index output.

Expected input: tab-delimited file with one row per (sample, chrom, window_start, window_end) and
a hybrid index column plus (optionally) CI bounds. Supported column names include:

  - hybrid index: mean_hybrid_index | hyb_index | hybrid_index | h
  - CI:           ci_low/ci_high    | lower_ci/upper_ci
  - extras:       mean_n_loci | n_aims | n_loci, and n_reps_used (optional)

Classification modes:
  * ci   (default): parent1 if ci_high <= t1; parent2 if ci_low >= t2; else admixed
  * mean:           parent1 if hi <= t1;     parent2 if hi >= t2;     else admixed
  * existing:       use the input "ancestry" column if present

Windows are merged within each (sample, chrom) when consecutive windows are within merge_distance
(i.e., next_start <= current_end + merge_distance).

Outputs:
  1) --output: per-sample merged regions.
  2) --cohort-output (optional): cohort-level merged regions with support counts across samples.

Hybrid-index orientation assumed:
  hi ~ 0 => Parent1
  hi ~ 1 => Parent2
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Optional, Tuple, List

import numpy as np
import pandas as pd


def _first_present(cols: List[str], df_cols: pd.Index) -> Optional[str]:
    for c in cols:
        if c in df_cols:
            return c
    return None


@dataclass(frozen=True)
class ColMap:
    sample: str
    chrom: str
    ws: str
    we: str
    hi: str
    ci_low: Optional[str]
    ci_high: Optional[str]
    n_loci: Optional[str]
    n_reps: Optional[str]


def infer_columns(df: pd.DataFrame) -> ColMap:
    dfc = df.columns

    sample = _first_present(["sample", "ind", "individual"], dfc)
    chrom  = _first_present(["chrom", "#CHROM", "CHROM"], dfc)
    ws     = _first_present(["window_start", "start", "win_start"], dfc)
    we     = _first_present(["window_end", "end", "win_end"], dfc)
    hi     = _first_present(["mean_hybrid_index", "hyb_index", "hybrid_index", "h", "HI"], dfc)

    if not all([sample, chrom, ws, we, hi]):
        missing = [k for k, v in dict(
            sample=sample, chrom=chrom, window_start=ws, window_end=we, hybrid_index=hi
        ).items() if v is None]
        raise ValueError(f"Missing required columns: {missing}")

    ci_low  = _first_present(["ci_low", "lower_ci", "ci_lower", "lowerCI"], dfc)
    ci_high = _first_present(["ci_high", "upper_ci", "ci_upper", "upperCI"], dfc)

    n_loci  = _first_present(["mean_n_loci", "n_aims", "n_loci", "mean_n_aims"], dfc)
    n_reps  = _first_present(["n_reps_used", "n_reps", "reps"], dfc)

    return ColMap(sample=sample, chrom=chrom, ws=ws, we=we, hi=hi,
                  ci_low=ci_low, ci_high=ci_high, n_loci=n_loci, n_reps=n_reps)


def classify_windows(df: pd.DataFrame, cols: ColMap, t1: float, t2: float, mode: str) -> pd.DataFrame:
    out = df.copy()

    out["hi"] = pd.to_numeric(out[cols.hi], errors="coerce")
    out["window_start"] = pd.to_numeric(out[cols.ws], errors="coerce").astype("Int64")
    out["window_end"]   = pd.to_numeric(out[cols.we], errors="coerce").astype("Int64")
    out["chrom"] = out[cols.chrom].astype(str)
    out["sample"] = out[cols.sample].astype(str)

    if cols.ci_low and cols.ci_high:
        out["ci_low"]  = pd.to_numeric(out[cols.ci_low], errors="coerce")
        out["ci_high"] = pd.to_numeric(out[cols.ci_high], errors="coerce")
    else:
        out["ci_low"]  = np.nan
        out["ci_high"] = np.nan

    if mode == "existing":
        if "ancestry" not in out.columns:
            raise ValueError("mode=existing requested but input has no 'ancestry' column.")
        out["call"] = out["ancestry"].astype(str)
        return out

    if mode == "ci":
        if not (cols.ci_low and cols.ci_high):
            raise ValueError("mode=ci requested but CI columns not found in input.")
        out["call"] = np.where(out["ci_high"] <= t1, "parent1",
                        np.where(out["ci_low"] >= t2, "parent2", "admixed"))
        return out

    if mode == "mean":
        out["call"] = np.where(out["hi"] <= t1, "parent1",
                        np.where(out["hi"] >= t2, "parent2", "admixed"))
        return out

    raise ValueError(f"Unknown mode: {mode}")


def merge_regions_per_sample(
    df_called: pd.DataFrame,
    cols: ColMap,
    target: Tuple[str, ...],
    merge_distance: int,
) -> pd.DataFrame:
    d = df_called[df_called["call"].isin(target)].copy()
    if d.empty:
        return pd.DataFrame(columns=[
            "sample","chrom","region_start","region_end","length",
            "n_windows","mean_hybrid_index","max_hybrid_index",
            "mean_n_loci","min_ci_low","max_ci_high"
        ])

    d = d.sort_values(["sample","chrom","window_start","window_end"])

    have_n_loci = cols.n_loci is not None and cols.n_loci in df_called.columns

    out_rows = []
    for (sample, chrom), g in d.groupby(["sample","chrom"], sort=False):
        region_start = None
        region_end = None

        n_windows = 0
        sum_hi = 0.0
        max_hi = -np.inf

        sum_nl = 0.0
        n_nl = 0

        min_ci_low = np.inf
        max_ci_high = -np.inf

        def flush():
            nonlocal region_start, region_end, n_windows, sum_hi, max_hi, sum_nl, n_nl, min_ci_low, max_ci_high
            if region_start is None:
                return
            out_rows.append({
                "sample": sample,
                "chrom": chrom,
                "region_start": int(region_start),
                "region_end": int(region_end),
                "length": int(region_end - region_start),
                "n_windows": int(n_windows),
                "mean_hybrid_index": (sum_hi / n_windows) if n_windows else np.nan,
                "max_hybrid_index": max_hi if n_windows else np.nan,
                "mean_n_loci": (sum_nl / n_nl) if n_nl else np.nan,
                "min_ci_low": min_ci_low if min_ci_low != np.inf else np.nan,
                "max_ci_high": max_ci_high if max_ci_high != -np.inf else np.nan,
            })

        for _, row in g.iterrows():
            ws = int(row["window_start"])
            we = int(row["window_end"])

            if region_start is None:
                region_start = ws
                region_end = we
            elif ws <= region_end + merge_distance:
                region_end = max(region_end, we)
            else:
                flush()
                region_start = ws
                region_end = we
                n_windows = 0
                sum_hi = 0.0
                max_hi = -np.inf
                sum_nl = 0.0
                n_nl = 0
                min_ci_low = np.inf
                max_ci_high = -np.inf

            if pd.notna(row["hi"]):
                n_windows += 1
                sum_hi += float(row["hi"])
                max_hi = max(max_hi, float(row["hi"]))

            if have_n_loci and pd.notna(row[cols.n_loci]):
                sum_nl += float(row[cols.n_loci])
                n_nl += 1

            if pd.notna(row["ci_low"]):
                min_ci_low = min(min_ci_low, float(row["ci_low"]))
            if pd.notna(row["ci_high"]):
                max_ci_high = max(max_ci_high, float(row["ci_high"]))

        flush()

    out = pd.DataFrame(out_rows).sort_values(["sample","chrom","region_start","region_end"])
    return out


def cohort_supported_regions(
    df_called: pd.DataFrame,
    min_samples: int,
    target: Tuple[str, ...],
    merge_distance: int,
) -> pd.DataFrame:
    win = (df_called[df_called["call"].isin(target)]
           .groupby(["chrom","window_start","window_end"], as_index=False)
           .agg(n_samples_support=("sample","nunique")))

    if win.empty:
        return pd.DataFrame(columns=[
            "chrom","region_start","region_end","length","n_windows",
            "mean_samples_support","max_samples_support"
        ])

    win = win[win["n_samples_support"] >= min_samples].copy()
    if win.empty:
        return pd.DataFrame(columns=[
            "chrom","region_start","region_end","length","n_windows",
            "mean_samples_support","max_samples_support"
        ])

    win = win.sort_values(["chrom","window_start","window_end"])

    out_rows = []
    for chrom, g in win.groupby("chrom", sort=False):
        region_start = None
        region_end = None
        n_windows = 0
        sum_support = 0.0
        max_support = 0

        def flush():
            nonlocal region_start, region_end, n_windows, sum_support, max_support
            if region_start is None:
                return
            out_rows.append({
                "chrom": chrom,
                "region_start": int(region_start),
                "region_end": int(region_end),
                "length": int(region_end - region_start),
                "n_windows": int(n_windows),
                "mean_samples_support": (sum_support / n_windows) if n_windows else np.nan,
                "max_samples_support": int(max_support),
            })

        for _, row in g.iterrows():
            ws = int(row["window_start"])
            we = int(row["window_end"])
            sup = int(row["n_samples_support"])

            if region_start is None:
                region_start, region_end = ws, we
                n_windows, sum_support, max_support = 0, 0.0, 0
            elif ws <= region_end + merge_distance:
                region_end = max(region_end, we)
            else:
                flush()
                region_start, region_end = ws, we
                n_windows, sum_support, max_support = 0, 0.0, 0

            n_windows += 1
            sum_support += sup
            max_support = max(max_support, sup)

        flush()

    return pd.DataFrame(out_rows).sort_values(["chrom","region_start","region_end"])


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--parent1-threshold", type=float, default=0.1)
    p.add_argument("--parent2-threshold", type=float, default=0.9)
    p.add_argument("--merge-distance", type=int, default=0)

    p.add_argument("--mode", choices=["ci","mean","existing"], default="ci")
    p.add_argument("--target", choices=["admixed","parent2","admixed_or_parent2"], default="admixed")

    p.add_argument("--cohort-output", default=None)
    p.add_argument("--cohort-min-samples", type=int, default=2)

    return p.parse_args()


def main():
    args = parse_args()

    df = pd.read_csv(args.input, sep="\t")
    cols = infer_columns(df)

    df_called = classify_windows(df, cols, args.parent1_threshold, args.parent2_threshold, args.mode)

    if args.target == "admixed":
        target_calls = ("admixed",)
    elif args.target == "parent2":
        target_calls = ("parent2",)
    else:
        target_calls = ("admixed", "parent2")

    regions = merge_regions_per_sample(df_called, cols, target_calls, args.merge_distance)
    regions.to_csv(args.output, sep="\t", index=False)

    if args.cohort_output:
        cohort = cohort_supported_regions(df_called, args.cohort_min_samples, target_calls, args.merge_distance)
        cohort.to_csv(args.cohort_output, sep="\t", index=False)

    print(f"Wrote {len(regions)} merged regions -> {args.output}")
    if args.cohort_output:
        print(f"Wrote {len(cohort)} cohort regions -> {args.cohort_output}")


if __name__ == "__main__":
    main()
