#!/usr/bin/env python3
"""
plot_hybrid_index.py

Publishable plotting of sliding-window hybrid index:
  - horizontal ancestry bands (Parent1 / Admixed / Parent2)
  - mean across individuals with ±1 SD and ±1 SE across individuals
  - top-N peak annotation with upstream/downstream gene context from a GFF

Input: tab-delimited windows file with one row per (sample, chrom, window_start, window_end).
Supported column names:
  hybrid index: mean_hybrid_index | hyb_index | hybrid_index | h
  CI:           ci_low/ci_high    | lower_ci/upper_ci (optional)

Outputs:
  --output-summary: single-page plot (PNG/PDF)
  --output-individual: multi-page PDF with one plot per sample (optional)
  --output-yaml / --output-peaks-tsv: peak table with gene context (optional)

Hybrid-index orientation assumed:
  hi ~ 0 => Parent1
  hi ~ 1 => Parent2
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from urllib.parse import unquote

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")  # safe on headless nodes
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

try:
    import yaml
except Exception:
    yaml = None


PARENT1_COLOR = "#117BE6"  # blue
ADMIXED_COLOR = "#FFF3E0"  # light orange
PARENT2_COLOR = "#B13447"  # pink


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

    return ColMap(sample=sample, chrom=chrom, ws=ws, we=we, hi=hi, ci_low=ci_low, ci_high=ci_high)


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


def _clean_product(note: str, max_len: int = 40) -> str:
    if not note:
        return ""
    m = re.search(r"Similar to\s+(.+?)(?:\s*\(|;|$)", note)
    s = m.group(1).strip() if m else note.strip()
    s = s.rstrip(":")
    if len(s) > max_len:
        s = s[: max_len - 3] + "..."
    return s


def gene_context(
    genes_by_chrom: Dict[str, List[Gene]],
    chrom: str,
    win_start: int,
    win_end: int,
    flank: int,
) -> List[Tuple[Gene, int]]:
    genes = genes_by_chrom.get(chrom, [])
    if not genes:
        return []

    s = win_start - flank
    e = win_end + flank
    near = [g for g in genes if g.end >= s and g.start <= e]

    def dist(g: Gene) -> int:
        if g.end < win_start:
            return win_start - g.end
        if g.start > win_end:
            return g.start - win_end
        return 0

    if near:
        near_sorted = sorted(near, key=lambda g: (dist(g), g.start))
        return [(g, dist(g)) for g in near_sorted]

    upstream = [g for g in genes if g.end < win_start]
    downstream = [g for g in genes if g.start > win_end]
    out: List[Tuple[Gene, int]] = []
    if upstream:
        g = max(upstream, key=lambda x: x.end)
        out.append((g, win_start - g.end))
    if downstream:
        g = min(downstream, key=lambda x: x.start)
        out.append((g, g.start - win_end))
    return out


def choose_peaks(summary: pd.DataFrame, top_n: int, min_distance_bp: int) -> pd.DataFrame:
    if top_n <= 0:
        return summary.iloc[0:0].copy()

    cand = summary.sort_values("mean_hi", ascending=False).copy()
    selected = []
    centers = []

    for _, row in cand.iterrows():
        if len(selected) >= top_n:
            break
        center = float(row["position_bp"])
        if min_distance_bp > 0 and any(abs(center - c) < min_distance_bp for c in centers):
            continue
        selected.append(row)
        centers.append(center)

    return pd.DataFrame(selected)


def plot_summary(
    summary: pd.DataFrame,
    chrom: str,
    n_samples: int,
    outpath: str,
    parent1_threshold: float,
    parent2_threshold: float,
    title: Optional[str],
    chrom_label: Optional[str],
    x_scale: float,
    peaks: pd.DataFrame,
    genes_by_chrom: Optional[Dict[str, List[Gene]]],
    flank: int,
    dpi: int,
    annotate: bool,
):
    fig, ax = plt.subplots(figsize=(16, 6))

    x = (summary["position_bp"] / x_scale).to_numpy()
    y = summary["mean_hi"].to_numpy()
    sd = summary["sd_hi"].to_numpy()
    se = summary["se_hi"].to_numpy()

    ax.axhspan(0, parent1_threshold, color=PARENT1_COLOR, alpha=0.20, lw=0)
    ax.axhspan(parent1_threshold, parent2_threshold, color=ADMIXED_COLOR, alpha=0.70, lw=0)
    ax.axhspan(parent2_threshold, 1, color=PARENT2_COLOR, alpha=0.25, lw=0)

    ax.fill_between(x, np.clip(y - sd, 0, 1), np.clip(y + sd, 0, 1),
                    color="gray", alpha=0.18, linewidth=0, zorder=1)
    ax.fill_between(x, np.clip(y - se, 0, 1), np.clip(y + se, 0, 1),
                    color="gray", alpha=0.30, linewidth=0, zorder=2)

    ax.plot(x, y, color="black", lw=2.0, zorder=3)

    ax.axhline(parent1_threshold, ls="--", lw=1.5, color="black", alpha=0.6)
    ax.axhline(parent2_threshold, ls="--", lw=1.5, color="black", alpha=0.6)

    if annotate:
        y_levels = [0.93, 0.86, 0.79, 0.72]
        for i, row in enumerate(peaks.itertuples(index=False), start=1):
            px = float(row.position_bp) / x_scale
            ax.axvline(px, color="red", ls=":", lw=2, alpha=0.65)

            label_lines = [f"Peak {i}"]
            if genes_by_chrom is not None:
                ctx = gene_context(genes_by_chrom, chrom, int(row.window_start), int(row.window_end), flank=flank)
                shown = 0
                for g, dist_bp in ctx[:2]:
                    prod = _clean_product(g.note)
                    if prod:
                        label_lines.append(f"{g.gene_id}: {prod}" + (f" ({dist_bp} bp)" if dist_bp else ""))
                    else:
                        label_lines.append(f"{g.gene_id}" + (f" ({dist_bp} bp)" if dist_bp else ""))
                    shown += 1
                extra = max(0, len(ctx) - shown)
                if extra:
                    label_lines.append(f"(+{extra} more)")

            ax.annotate(
                "\n".join(label_lines),
                xy=(px, float(row.mean_hi)),
                xytext=(px, y_levels[(i - 1) % len(y_levels)]),
                textcoords="data",
                ha="center",
                va="center",
                arrowprops=dict(arrowstyle="->", lw=1.5, color="red"),
                bbox=dict(boxstyle="round,pad=0.35", fc="yellow", ec="orange", alpha=0.95, lw=1.5),
                fontsize=8,
            )
    else:
        for row in peaks.itertuples(index=False):
            px = float(row.position_bp) / x_scale
            ax.axvline(px, color="red", ls=":", lw=2, alpha=0.65)

    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Hybrid Index")

    if annotate:
        ax.set_title(
            title if title else f"Mean Hybrid Index Across {n_samples} Samples (Top {len(peaks)} Peaks Annotated)",
            fontsize=14, fontweight="bold"
        )
    else:
        ax.set_title(
            title if title else f"Mean Hybrid Index Across {n_samples} Samples",
            fontsize=14, fontweight="bold"
        )

    if chrom_label:
        ax.text(0.5, -0.10, chrom_label, transform=ax.transAxes, ha="center", va="top", fontsize=11)

    legend_elements = [
        mpatches.Patch(facecolor=PARENT1_COLOR, alpha=0.40, label="Parent1"),
        mpatches.Patch(facecolor=ADMIXED_COLOR, alpha=0.60, label="Admixed"),
        mpatches.Patch(facecolor=PARENT2_COLOR, alpha=0.40, label="Parent2"),
        plt.Line2D([0], [0], color="gray", alpha=0.30, lw=8, label="± 1 SE"),
        plt.Line2D([0], [0], color="gray", alpha=0.18, lw=8, label="± 1 SD"),
        plt.Line2D([0], [0], color="black", lw=2, label="Mean"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(1.02, 1), 
              fontsize=10, framealpha=0.95, edgecolor="black")

    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_individuals(
    df: pd.DataFrame,
    cols: ColMap,
    outpath: str,
    parent1_threshold: float,
    parent2_threshold: float,
    chrom_label: Optional[str],
    x_scale: float,
    dpi: int,
):
    df = df.copy()
    df["hi"] = pd.to_numeric(df[cols.hi], errors="coerce")
    df["position_bp"] = (pd.to_numeric(df[cols.ws], errors="coerce") + pd.to_numeric(df[cols.we], errors="coerce")) / 2.0

    have_ci = cols.ci_low is not None and cols.ci_high is not None and cols.ci_low in df.columns and cols.ci_high in df.columns
    if have_ci:
        df["ci_low"] = pd.to_numeric(df[cols.ci_low], errors="coerce")
        df["ci_high"] = pd.to_numeric(df[cols.ci_high], errors="coerce")

    with PdfPages(outpath) as pdf:
        for sample, g in df.groupby(cols.sample, sort=True):
            g = g.sort_values("position_bp")
            x = (g["position_bp"] / x_scale).to_numpy()
            y = g["hi"].to_numpy()

            fig, ax = plt.subplots(figsize=(16, 4))
            ax.axhspan(0, parent1_threshold, color=PARENT1_COLOR, alpha=0.20, lw=0)
            ax.axhspan(parent1_threshold, parent2_threshold, color=ADMIXED_COLOR, alpha=0.70, lw=0)
            ax.axhspan(parent2_threshold, 1, color=PARENT2_COLOR, alpha=0.25, lw=0)

            if have_ci:
                ax.fill_between(x, g["ci_low"].to_numpy(), g["ci_high"].to_numpy(), color="gray", alpha=0.25, linewidth=0)

            ax.plot(x, y, color="black", lw=1.5)
            ax.axhline(parent1_threshold, ls="--", lw=1.2, color="black", alpha=0.6)
            ax.axhline(parent2_threshold, ls="--", lw=1.2, color="black", alpha=0.6)

            ax.set_ylim(-0.05, 1.05)
            ax.set_xlabel("Chromosome")
            ax.set_ylabel("Hybrid Index")
            ax.set_title(str(sample), fontsize=12, fontweight="bold")
            if chrom_label:
                ax.text(0.5, -0.10, chrom_label, transform=ax.transAxes, ha="center", va="top", fontsize=11)

            fig.tight_layout()
            pdf.savefig(fig, dpi=dpi, bbox_inches="tight")
            plt.close(fig)


def write_peaks_table(
    peaks: pd.DataFrame,
    chrom: str,
    genes_by_chrom: Optional[Dict[str, List[Gene]]],
    flank: int,
    out_tsv: Optional[str],
    out_yaml: Optional[str],
):
    if out_tsv is None and out_yaml is None:
        return

    rows = []
    for rank, row in enumerate(peaks.itertuples(index=False), start=1):
        entry = {
            "rank": rank,
            "chrom": chrom,
            "window_start": int(row.window_start),
            "window_end": int(row.window_end),
            "position_bp": float(row.position_bp),
            "mean_hi": float(row.mean_hi),
            "sd_hi": float(row.sd_hi),
            "se_hi": float(row.se_hi),
        }
        if genes_by_chrom is not None:
            ctx = gene_context(genes_by_chrom, chrom, int(row.window_start), int(row.window_end), flank=flank)
            entry["genes"] = [
                {
                    "gene_id": g.gene_id,
                    "start": g.start,
                    "end": g.end,
                    "strand": g.strand,
                    "distance_bp": int(dist),
                    "note": g.note,
                }
                for g, dist in ctx
            ]
        rows.append(entry)

    if out_tsv is not None:
        flat = []
        for r in rows:
            gene_str = ""
            if "genes" in r and r["genes"]:
                parts = []
                for g in r["genes"][:5]:
                    parts.append(g["gene_id"] if g["distance_bp"] == 0 else f'{g["gene_id"]}({g["distance_bp"]})')
                extra = max(0, len(r["genes"]) - len(parts))
                if extra:
                    parts.append(f"+{extra} more")
                gene_str = ",".join(parts)

            flat.append({
                "rank": r["rank"],
                "chrom": r["chrom"],
                "window_start": r["window_start"],
                "window_end": r["window_end"],
                "position_bp": int(r["position_bp"]),
                "mean_hi": r["mean_hi"],
                "sd_hi": r["sd_hi"],
                "se_hi": r["se_hi"],
                "gene_context": gene_str,
            })
        pd.DataFrame(flat).to_csv(out_tsv, sep="\t", index=False)

    if out_yaml is not None:
        if yaml is None:
            raise RuntimeError("PyYAML not available but --output-yaml was requested.")
        with open(out_yaml, "w") as f:
            yaml.safe_dump({"peaks": rows}, f, sort_keys=False)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output-summary", required=True)
    p.add_argument("--output-individual", default=None)

    p.add_argument("--gff", default=None)
    p.add_argument("--annotate", action="store_true", help="Annotate peaks with gene context (requires --gff)")
    p.add_argument("--output-yaml", default=None)
    p.add_argument("--output-peaks-tsv", default=None)

    p.add_argument("--top-peaks", type=int, default=10)
    p.add_argument("--peak-min-distance", type=int, default=0)
    p.add_argument("--flank", type=int, default=10000)

    p.add_argument("--parent1-threshold", type=float, default=0.1)
    p.add_argument("--parent2-threshold", type=float, default=0.9)

    p.add_argument("--title", default=None)
    p.add_argument("--chrom-label", default=None)

    p.add_argument("--x-scale", type=float, default=1e6)  # Mb by default
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def main():
    args = parse_args()

    if args.annotate and not args.gff:
        raise ValueError("--annotate requires --gff")

    df = pd.read_csv(args.input, sep="\t")
    cols = infer_columns(df)

    df["hi"] = pd.to_numeric(df[cols.hi], errors="coerce")
    df["window_start"] = pd.to_numeric(df[cols.ws], errors="coerce")
    df["window_end"] = pd.to_numeric(df[cols.we], errors="coerce")
    df["position_bp"] = (df["window_start"] + df["window_end"]) / 2.0
    df["chrom"] = df[cols.chrom].astype(str)
    df["sample"] = df[cols.sample].astype(str)

    chroms = df["chrom"].unique()
    if len(chroms) != 1:
        raise ValueError(f"Expected one chromosome per input file; found {len(chroms)}: {chroms}")
    chrom = chroms[0]

    n_samples = df["sample"].nunique()

    summary = (df.groupby(["chrom","window_start","window_end","position_bp"], as_index=False)
                 .agg(mean_hi=("hi","mean"),
                      sd_hi=("hi","std"),
                      n=("hi","count")))
    summary["sd_hi"] = summary["sd_hi"].fillna(0.0)
    summary["se_hi"] = summary["sd_hi"] / np.sqrt(summary["n"].clip(lower=1))
    summary = summary.sort_values("position_bp")

    peaks = choose_peaks(summary, top_n=args.top_peaks, min_distance_bp=args.peak_min_distance)

    genes_by_chrom = parse_gff_genes(args.gff) if args.gff else None

    plot_summary(
        summary=summary,
        chrom=chrom,
        n_samples=n_samples,
        outpath=args.output_summary,
        parent1_threshold=args.parent1_threshold,
        parent2_threshold=args.parent2_threshold,
        title=args.title,
        chrom_label=args.chrom_label or chrom,
        x_scale=args.x_scale,
        peaks=peaks,
        genes_by_chrom=genes_by_chrom,
        flank=args.flank,
        dpi=args.dpi,
        annotate=args.annotate,
    )

    if args.output_individual:
        plot_individuals(
            df=df,
            cols=cols,
            outpath=args.output_individual,
            parent1_threshold=args.parent1_threshold,
            parent2_threshold=args.parent2_threshold,
            chrom_label=args.chrom_label or chrom,
            x_scale=args.x_scale,
            dpi=args.dpi,
        )

    write_peaks_table(
        peaks=peaks,
        chrom=chrom,
        genes_by_chrom=genes_by_chrom,
        flank=args.flank,
        out_tsv=args.output_peaks_tsv,
        out_yaml=args.output_yaml,
    )

    print(f"Summary plot: {args.output_summary}")
    if args.output_individual:
        print(f"Individual plots: {args.output_individual}")
    if args.output_peaks_tsv:
        print(f"Peak table: {args.output_peaks_tsv}")
    if args.output_yaml:
        print(f"Peak YAML: {args.output_yaml}")


if __name__ == "__main__":
    main()