#!/usr/bin/env python3
"""
find_introgression_deserts.py

Identify long "introgression-resistant" (a.k.a. hybridization desert) regions where
windows are confidently assigned to a target ancestry across many samples.

Expected input columns (TSV):
  sample, chrom, window_start, window_end, mean_hybrid_index, ci_low, ci_high

Back-compat accepted:
  lower_ci/upper_ci instead of ci_low/ci_high
  hybrid_index/h instead of mean_hybrid_index

Default target = parent1 (S. catenatus), using CI-based classification:
  parent1 if ci_high <= parent1_threshold
  parent2 if ci_low  >= parent2_threshold
  else admixed
"""

import argparse
import sys
from dataclasses import dataclass
from typing import Dict, List
from urllib.parse import unquote

import pandas as pd


def _pick_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    cols = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols:
            return cols[cand.lower()]
    return None


def classify_ancestry(
    df: pd.DataFrame,
    parent1_threshold: float,
    parent2_threshold: float,
    use_ci: bool = True,
) -> pd.Series:
    """
    Returns: Series[str] in {'parent1','parent2','admixed','missing'}.
    """
    hi_col = _pick_col(df, ["mean_hybrid_index", "hybrid_index", "h"])
    lo_ci_col = _pick_col(df, ["ci_low", "lower_ci", "lci"])
    hi_ci_col = _pick_col(df, ["ci_high", "upper_ci", "uci"])

    if hi_col is None:
        raise ValueError(
            "Missing hybrid index column. Need one of: mean_hybrid_index, hybrid_index, h"
        )

    h = pd.to_numeric(df[hi_col], errors="coerce")

    if use_ci:
        if lo_ci_col is None or hi_ci_col is None:
            raise ValueError(
                "use_ci=True but CI columns not found. Need (ci_low, ci_high) or (lower_ci, upper_ci)."
            )
        lo = pd.to_numeric(df[lo_ci_col], errors="coerce")
        hi = pd.to_numeric(df[hi_ci_col], errors="coerce")

        out = pd.Series("missing", index=df.index, dtype="object")
        ok = lo.notna() & hi.notna()
        out.loc[ok & (hi <= parent1_threshold)] = "parent1"
        out.loc[ok & (lo >= parent2_threshold)] = "parent2"
        out.loc[ok & ~((hi <= parent1_threshold) | (lo >= parent2_threshold))] = "admixed"
        return out

    # mean-based fallback
    out = pd.Series("missing", index=df.index, dtype="object")
    ok = h.notna()
    out.loc[ok & (h <= parent1_threshold)] = "parent1"
    out.loc[ok & (h >= parent2_threshold)] = "parent2"
    out.loc[ok & ~((h <= parent1_threshold) | (h >= parent2_threshold))] = "admixed"
    return out


@dataclass
class Region:
    chrom: str
    start: int
    end: int
    n_windows: int
    mean_support: float
    min_support: float
    mean_hybrid_index: float
    n_samples_median: float


@dataclass(frozen=True)
class Gene:
    start: int
    end: int
    strand: str
    gene_id: str
    note: str


def parse_gff_genes(gff_path: str) -> Dict[str, List[Gene]]:
    genes: Dict[str, List[Gene]] = {}
    with open(gff_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _source, feature, start, end, _score, strand, _phase, attrs = parts
            if feature != "gene":
                continue

            start_i = int(start)
            end_i = int(end)

            d: Dict[str, str] = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    d[k] = unquote(v)

            gene_id = d.get("Name") or d.get("gene") or d.get("ID") or "unknown"
            note = d.get("product") or d.get("Note") or ""
            genes.setdefault(chrom, []).append(Gene(start=start_i, end=end_i, strand=strand, gene_id=gene_id, note=note))

    for chrom in genes:
        genes[chrom].sort(key=lambda g: (g.start, g.end))
    return genes


def find_genes_in_region(
    genes_by_chrom: Dict[str, List[Gene]],
    chrom: str,
    region_start: int,
    region_end: int,
) -> List[Gene]:
    """Find all genes that overlap with the region."""
    genes = genes_by_chrom.get(chrom, [])
    overlapping = [g for g in genes if g.end >= region_start and g.start <= region_end]
    return overlapping


def build_window_summary(
    df: pd.DataFrame,
    ancestry_col: str,
) -> pd.DataFrame:
    """
    Per-window summary across samples.
    """
    # Required window keys
    for c in ["sample", "chrom", "window_start", "window_end"]:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")

    hi_col = _pick_col(df, ["mean_hybrid_index", "hybrid_index", "h"])
    if hi_col is None:
        raise ValueError(
            "Missing hybrid index column. Need one of: mean_hybrid_index, hybrid_index, h"
        )

    dd = df[df[ancestry_col] != "missing"].copy()
    dd["__h__"] = pd.to_numeric(dd[hi_col], errors="coerce")

    gb = dd.groupby(["chrom", "window_start", "window_end"], sort=False)

    def _agg(g: pd.DataFrame) -> pd.Series:
        a = g[ancestry_col]
        return pd.Series(
            {
                "n_samples": g["sample"].nunique(),
                "n_parent1": int((a == "parent1").sum()),
                "n_parent2": int((a == "parent2").sum()),
                "n_admixed": int((a == "admixed").sum()),
                "parent1_frac": float((a == "parent1").mean()),
                "parent2_frac": float((a == "parent2").mean()),
                "admixed_frac": float((a == "admixed").mean()),
                "mean_hybrid_index": float(g["__h__"].mean()),
            }
        )

    out = gb.apply(_agg).reset_index()
    return out


def merge_windows_into_regions(
    win: pd.DataFrame,
    support_col: str,
    min_samples: int,
    min_support: float,
    merge_distance: int,
    min_length: int,
    min_windows: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns:
      (filtered_windows_used_for_merging, regions_df)
    """
    w = win.copy()
    w = w[(w["n_samples"] >= min_samples) & (w[support_col] >= min_support)].copy()
    w = w.sort_values(["chrom", "window_start", "window_end"], kind="mergesort")

    regions: list[Region] = []

    for chrom, g in w.groupby("chrom", sort=False):
        cur_start: int | None = None
        cur_end: int | None = None
        rows: list[pd.Series] = []

        for _, row in g.iterrows():
            s = int(row["window_start"])
            e = int(row["window_end"])
            if cur_start is None:
                cur_start, cur_end = s, e
                rows = [row]
                continue

            if s <= (cur_end + merge_distance):
                cur_end = max(cur_end, e)
                rows.append(row)
            else:
                rr = pd.DataFrame(rows)
                regions.append(
                    Region(
                        chrom=chrom,
                        start=int(cur_start),
                        end=int(cur_end),
                        n_windows=int(len(rows)),
                        mean_support=float(rr[support_col].mean()),
                        min_support=float(rr[support_col].min()),
                        mean_hybrid_index=float(rr["mean_hybrid_index"].mean()),
                        n_samples_median=float(rr["n_samples"].median()),
                    )
                )
                cur_start, cur_end = s, e
                rows = [row]

        if cur_start is not None and cur_end is not None and rows:
            rr = pd.DataFrame(rows)
            regions.append(
                Region(
                    chrom=chrom,
                    start=int(cur_start),
                    end=int(cur_end),
                    n_windows=int(len(rows)),
                    mean_support=float(rr[support_col].mean()),
                    min_support=float(rr[support_col].min()),
                    mean_hybrid_index=float(rr["mean_hybrid_index"].mean()),
                    n_samples_median=float(rr["n_samples"].median()),
                )
            )

    reg_df = pd.DataFrame(
        [
            {
                "chrom": r.chrom,
                "region_start": r.start,
                "region_end": r.end,
                "length": r.end - r.start,
                "n_windows": r.n_windows,
                "mean_support": r.mean_support,
                "min_support": r.min_support,
                "mean_hybrid_index": r.mean_hybrid_index,
                "n_samples_median": r.n_samples_median,
            }
            for r in regions
        ]
    )

    if not reg_df.empty:
        reg_df = reg_df[(reg_df["length"] >= min_length) & (reg_df["n_windows"] >= min_windows)].copy()
        reg_df = reg_df.sort_values(["chrom", "region_start"], kind="mergesort")

    return w, reg_df


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Find introgression-resistant (hybridization desert) regions from hybrid-index windows."
    )
    p.add_argument("--input", nargs="+", required=True, help="One or more hybrid_index_windows TSVs.")
    p.add_argument("--output", required=True, help="Output TSV of merged desert regions.")
    p.add_argument(
        "--output-windows",
        default=None,
        help="Optional: write per-window across-sample summary TSV (all windows).",
    )
    p.add_argument(
        "--output-windows-filtered",
        default=None,
        help="Optional: write per-window summary TSV after (min_samples/min_support) filtering.",
    )

    p.add_argument("--target", choices=["parent1", "parent2", "admixed"], default="parent1")
    p.add_argument("--parent1-threshold", type=float, default=0.1)
    p.add_argument("--parent2-threshold", type=float, default=0.9)
    p.add_argument("--no-ci", action="store_true", help="Use mean_hybrid_index thresholds (ignore CI).")

    p.add_argument("--min-samples", type=int, default=7, help="Min samples with data per window.")
    p.add_argument("--min-support", type=float, default=0.8, help="Min fraction of samples in target state.")
    p.add_argument("--merge-distance", type=int, default=0, help="Merge windows if gap <= this many bp.")
    p.add_argument("--min-length", type=int, default=200_000, help="Min merged region length (bp).")
    p.add_argument("--min-windows", type=int, default=3, help="Min windows per merged region.")
    p.add_argument("--chrom", default=None, help="Optional: restrict to a single chrom/contig ID.")
    
    p.add_argument("--gff", default=None, help="Optional: GFF file to annotate regions with overlapping genes.")

    return p.parse_args()


def main() -> None:
    args = parse_args()

    # Load inputs
    dfs = []
    for fp in args.input:
        d = pd.read_csv(fp, sep="\t", low_memory=False)
        dfs.append(d)
    df = pd.concat(dfs, ignore_index=True)

    # Optional chrom filter
    if args.chrom is not None:
        df = df[df["chrom"] == args.chrom].copy()

    # Classify
    df["__ancestry__"] = classify_ancestry(
        df,
        parent1_threshold=args.parent1_threshold,
        parent2_threshold=args.parent2_threshold,
        use_ci=(not args.no_ci),
    )

    # Window summary across samples
    win = build_window_summary(df, ancestry_col="__ancestry__")

    # Choose support column
    if args.target == "parent1":
        support_col = "parent1_frac"
    elif args.target == "parent2":
        support_col = "parent2_frac"
    else:
        support_col = "admixed_frac"

    win["support"] = win[support_col]

    if args.output_windows is not None:
        win.to_csv(args.output_windows, sep="\t", index=False)

    # Merge into regions
    win_filt, regions = merge_windows_into_regions(
        win=win,
        support_col="support",
        min_samples=args.min_samples,
        min_support=args.min_support,
        merge_distance=args.merge_distance,
        min_length=args.min_length,
        min_windows=args.min_windows,
    )

    if args.output_windows_filtered is not None:
        win_filt.to_csv(args.output_windows_filtered, sep="\t", index=False)

    regions.insert(0, "target", args.target)

    # Annotate with genes if GFF provided
    if args.gff is not None and not regions.empty:
        genes_by_chrom = parse_gff_genes(args.gff)
        
        gene_ids = []
        gene_counts = []
        
        for _, row in regions.iterrows():
            overlapping = find_genes_in_region(
                genes_by_chrom,
                row["chrom"],
                int(row["region_start"]),
                int(row["region_end"])
            )
            
            gene_counts.append(len(overlapping))
            if overlapping:
                gene_ids.append(",".join(g.gene_id for g in overlapping))
            else:
                gene_ids.append("")
        
        regions["n_genes"] = gene_counts
        regions["genes"] = gene_ids

    regions.to_csv(args.output, sep="\t", index=False)

    # Minimal stdout summary
    n_reg = 0 if regions is None or regions.empty else len(regions)
    print(f"[OK] target={args.target} regions={n_reg} written: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()