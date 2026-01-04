#!/usr/bin/env python3
"""
calculate_hybrid_index_ml.py

Sliding-window hybrid index estimator using diagnostic (fixed-difference) markers,
following the likelihood framework of Buerkle (2005).

Fixes implemented:
- No external precomputed frequency files (eliminates CHROM/POS row misalignment bugs).
- Skips partial / non-diploid genotypes (e.g. 0/., ./1, haploid), so fixed-model math is valid.
- Defensively skips multi-allelic and non-SNP records (but you should still filter upstream).
- Optional parent1 subsampling with permutations to mitigate parent sample-size mismatch:
  repeat hybrid index estimation using random subsets of parent1 and summarize by mean
  and percentile CI across replicates.

Expected input VCF:
- Contains *input samples* + *both parental sets*.
- Typically already subset to AIMs / candidate loci.

Output:
TSV: per-sample, per-window mean hybrid index and CI.
- If --n-reps > 1: CI = empirical percentiles across replicates.
- If --n-reps == 1: CI = support interval based on 2 log-likelihood units (Buerkle-style),
  assuming independent loci.
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


def eprint(*args, **kwargs) -> None:
    print(*args, file=sys.stderr, **kwargs)


def read_sample_list(path: str) -> List[str]:
    samples: List[str] = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            samples.append(s)
    # preserve order but drop duplicates
    seen = set()
    out = []
    for s in samples:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


def open_textmaybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_gt_to_n_alt(gt: str) -> Optional[int]:
    """
    Parse a GT string (e.g., '0/1', '1|1') and return ALT allele count (0,1,2).
    Returns None if missing, partially missing, non-diploid, or non-biallelic alleles.
    """
    if gt in (".", "./.", ".|."):
        return None
    sep = "/" if "/" in gt else ("|" if "|" in gt else None)
    if sep is None:
        return None
    a = gt.split(sep)
    if len(a) != 2:
        return None
    if a[0] == "." or a[1] == ".":
        return None
    if a[0] not in ("0", "1") or a[1] not in ("0", "1"):
        return None
    return (a[0] == "1") + (a[1] == "1")


def safe_log(x: float, eps: float = 1e-12) -> float:
    return math.log(min(max(x, eps), 1.0 - eps))


def classify_ancestry(h: float, low: float, high: float) -> str:
    if math.isnan(h):
        return "NA"
    if h < low:
        return "parent1"
    if h > high:
        return "parent2"
    return "admixed"


def support_interval_binomial(k: int, n: int, delta_ll: float = 2.0) -> Tuple[float, float]:
    """
    Likelihood support interval for h in Binomial(n, h) (constants omitted), i.e.
      ll(h) = k*log(h) + (n-k)*log(1-h)
    Return bounds where ll(h) >= ll(h_hat) - delta_ll.

    Handles edge cases k=0 or k=n analytically.
    """
    if n <= 0:
        return (float("nan"), float("nan"))

    if k == 0:
        upper = 1.0 - math.exp(-delta_ll / n)
        return (0.0, min(1.0, upper))
    if k == n:
        lower = math.exp(-delta_ll / n)
        return (max(0.0, lower), 1.0)

    h_hat = k / n
    ll_max = k * safe_log(h_hat) + (n - k) * safe_log(1.0 - h_hat)
    target = ll_max - delta_ll

    def ll(h: float) -> float:
        return k * safe_log(h) + (n - k) * safe_log(1.0 - h)

    # Lower bound in (0, h_hat)
    lo, hi = 0.0, h_hat
    if ll(lo) >= target:
        lower = lo
    else:
        for _ in range(60):
            mid = (lo + hi) / 2.0
            if ll(mid) >= target:
                hi = mid
            else:
                lo = mid
        lower = hi

    # Upper bound in (h_hat, 1)
    lo, hi = h_hat, 1.0
    if ll(hi) >= target:
        upper = hi
    else:
        for _ in range(60):
            mid = (lo + hi) / 2.0
            if ll(mid) >= target:
                lo = mid
            else:
                hi = mid
        upper = lo

    return (max(0.0, lower), min(1.0, upper))


@dataclass
class ChromData:
    chrom: str
    positions: np.ndarray
    g_input: np.ndarray
    g_p1: np.ndarray
    g_p2: np.ndarray


def vcf_chrom_iterator(
    vcf_path: str,
    input_indices: Sequence[int],
    p1_indices: Sequence[int],
    p2_indices: Sequence[int],
) -> Iterable[ChromData]:
    """Stream a VCF and yield per-chromosome genotype matrices for the three groups."""
    positions: List[int] = []
    g_in: List[List[int]] = []
    g_p1: List[List[int]] = []
    g_p2: List[List[int]] = []

    current_chrom: Optional[str] = None

    with open_textmaybe_gzip(vcf_path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue
            if not line or line[0] == "#":
                continue

            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            pos = int(parts[1])

            ref = parts[3]
            alt = parts[4]
            if "," in alt:
                continue  # multi-allelic
            if len(ref) != 1 or len(alt) != 1:
                continue  # not a SNP

            fmt = parts[8].split(":")
            try:
                gt_i = fmt.index("GT")
            except ValueError:
                gt_i = 0

            samples = parts[9:]

            def get_gt(idx: int) -> Optional[int]:
                fields = samples[idx].split(":")
                if gt_i >= len(fields):
                    return None
                return parse_gt_to_n_alt(fields[gt_i])

            if current_chrom is None:
                current_chrom = chrom
            if chrom != current_chrom:
                yield ChromData(
                    chrom=current_chrom,
                    positions=np.asarray(positions, dtype=np.int32),
                    g_input=np.asarray(g_in, dtype=np.int8),
                    g_p1=np.asarray(g_p1, dtype=np.int8),
                    g_p2=np.asarray(g_p2, dtype=np.int8),
                )
                current_chrom = chrom
                positions, g_in, g_p1, g_p2 = [], [], [], []

            positions.append(pos)

            row_in = [get_gt(i) for i in input_indices]
            row_p1 = [get_gt(i) for i in p1_indices]
            row_p2 = [get_gt(i) for i in p2_indices]

            g_in.append([x if x is not None else -1 for x in row_in])
            g_p1.append([x if x is not None else -1 for x in row_p1])
            g_p2.append([x if x is not None else -1 for x in row_p2])

    if current_chrom is not None and positions:
        yield ChromData(
            chrom=current_chrom,
            positions=np.asarray(positions, dtype=np.int32),
            g_input=np.asarray(g_in, dtype=np.int8),
            g_p1=np.asarray(g_p1, dtype=np.int8),
            g_p2=np.asarray(g_p2, dtype=np.int8),
        )


def compute_alt_freq_and_called(g: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute per-variant ALT allele frequency and number of called diploid genotypes.
    g: (V, n_samples) int8 with -1 for missing, else 0/1/2
    Returns:
      freq_alt: (V,) float64, NaN where no calls
      called_genotypes: (V,) int32
    """
    called_mask = (g >= 0)
    called = called_mask.sum(axis=1).astype(np.int32)
    alt_count = (g.astype(np.int32) * called_mask.astype(np.int32)).sum(axis=1).astype(np.int32)
    total_alleles = 2 * called

    freq = np.full(g.shape[0], np.nan, dtype=np.float64)
    nonzero = total_alleles > 0
    freq[nonzero] = alt_count[nonzero] / total_alleles[nonzero]
    return freq, called


def precompute_windows(positions: np.ndarray, window_size: int, step_size: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if positions.size == 0:
        return (np.zeros(0, dtype=np.int32), np.zeros(0, dtype=np.int32), np.zeros((0, 2), dtype=np.int32))

    min_pos = int(positions[0])
    max_pos = int(positions[-1])

    starts = np.arange(min_pos, max_pos + 1, step_size, dtype=np.int32)
    ends = starts + np.int32(window_size)

    i0 = np.searchsorted(positions, starts, side="left").astype(np.int32)
    i1 = np.searchsorted(positions, ends, side="left").astype(np.int32)
    bounds = np.stack([i0, i1], axis=1)
    return starts, ends, bounds


def process_chrom(
    cd: ChromData,
    input_names: Sequence[str],
    p2_freq: np.ndarray,
    p2_called: np.ndarray,
    p1_subsets: Sequence[np.ndarray],
    diagnostic_threshold: float,
    min_parent_called: int,
    window_size: int,
    step_size: int,
    min_aims: int,
    chunk_windows: int,
    ancestry_low: float,
    ancestry_high: float,
    out_handle,
    ci_method: str,
) -> None:
    V = cd.positions.size
    if V == 0:
        return

    starts, ends, bounds = precompute_windows(cd.positions, window_size, step_size)
    W = starts.size
    if W == 0:
        return

    n_reps = len(p1_subsets)
    n_input = cd.g_input.shape[1]

    for w0 in range(0, W, chunk_windows):
        w1 = min(W, w0 + chunk_windows)
        bounds_chunk = bounds[w0:w1, :]
        starts_chunk = starts[w0:w1]
        ends_chunk = ends[w0:w1]
        Wc = w1 - w0

        h_rep = np.full((n_reps, n_input, Wc), np.nan, dtype=np.float32)
        m_rep = np.full((n_reps, n_input, Wc), np.nan, dtype=np.float32)

        # For support CI (n_reps==1), store exact counts (avoid float rounding).
        k_rep0 = None
        c_rep0 = None
        if ci_method == "support" and n_reps == 1:
            k_rep0 = np.full((n_input, Wc), -1, dtype=np.int32)
            c_rep0 = np.full((n_input, Wc), -1, dtype=np.int32)

        for r, subcols in enumerate(p1_subsets):
            g_p1_sub = cd.g_p1[:, subcols] if subcols.size > 0 else cd.g_p1
            p1_freq, p1_called = compute_alt_freq_and_called(g_p1_sub)

            ok_parent = (p1_called >= min_parent_called) & (p2_called >= min_parent_called)
            diag_alt = ok_parent & (p1_freq < diagnostic_threshold) & (p2_freq > (1.0 - diagnostic_threshold))
            diag_ref = ok_parent & (p1_freq > (1.0 - diagnostic_threshold)) & (p2_freq < diagnostic_threshold)
            diag = diag_alt | diag_ref
            if not np.any(diag):
                continue
            p2_is_alt = diag_alt

            for j in range(n_input):
                n_alt = cd.g_input[:, j].astype(np.int16)
                valid = diag & (n_alt >= 0)
                if not np.any(valid):
                    continue

                n_p2 = np.where(p2_is_alt, n_alt, 2 - n_alt)
                n_p2 = np.where(valid, n_p2, 0).astype(np.int32)
                c = valid.astype(np.int32)

                pref_p2 = np.concatenate(([0], np.cumsum(n_p2, dtype=np.int64)))
                pref_c = np.concatenate(([0], np.cumsum(c, dtype=np.int64)))

                i0 = bounds_chunk[:, 0]
                i1 = bounds_chunk[:, 1]
                sum_p2 = pref_p2[i1] - pref_p2[i0]
                sum_c = pref_c[i1] - pref_c[i0]

                ok = sum_c >= min_aims
                if not np.any(ok):
                    continue

                h = np.full(Wc, np.nan, dtype=np.float32)
                h[ok] = (sum_p2[ok] / (2.0 * sum_c[ok])).astype(np.float32)

                h_rep[r, j, :] = h
                m_rep[r, j, :] = np.where(ok, sum_c.astype(np.float32), np.nan)

                if k_rep0 is not None and r == 0:
                    k_rep0[j, :] = np.where(ok, sum_p2.astype(np.int32), -1)
                    c_rep0[j, :] = np.where(ok, sum_c.astype(np.int32), -1)

        reps_used = np.sum(~np.isnan(h_rep), axis=0).astype(np.int32)
        mean_h = np.nanmean(h_rep, axis=0)
        mean_m = np.nanmean(m_rep, axis=0)

        if ci_method == "perm":
            ci_lo = np.nanquantile(h_rep, 0.025, axis=0)
            ci_hi = np.nanquantile(h_rep, 0.975, axis=0)
        else:
            ci_lo = np.full_like(mean_h, np.nan, dtype=np.float32)
            ci_hi = np.full_like(mean_h, np.nan, dtype=np.float32)
            if n_reps == 1 and k_rep0 is not None and c_rep0 is not None:
                for j in range(n_input):
                    for w in range(Wc):
                        k = int(k_rep0[j, w])
                        c = int(c_rep0[j, w])
                        if k < 0 or c < 0:
                            continue
                        lo, hi = support_interval_binomial(k, 2 * c, delta_ll=2.0)
                        ci_lo[j, w] = lo
                        ci_hi[j, w] = hi
            else:
                ci_lo = np.nanquantile(h_rep, 0.025, axis=0)
                ci_hi = np.nanquantile(h_rep, 0.975, axis=0)

        for j, sample in enumerate(input_names):
            for w in range(Wc):
                if reps_used[j, w] == 0:
                    continue
                h = float(mean_h[j, w])
                lo = float(ci_lo[j, w])
                hi = float(ci_hi[j, w])
                m = float(mean_m[j, w])
                anc = classify_ancestry(h, ancestry_low, ancestry_high)
                out_handle.write(
                    f"{sample}\t{cd.chrom}\t{int(starts_chunk[w])}\t{int(ends_chunk[w])}\t"
                    f"{m:.2f}\t{h:.6f}\t{lo:.6f}\t{hi:.6f}\t{anc}\t{int(reps_used[j, w])}\n"
                )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Sliding-window hybrid index (diagnostic markers) with optional parent1 subsampling permutations."
    )
    p.add_argument("--vcf", required=True, help="VCF/VCF.GZ containing input + parents (typically already subset to AIMs).")
    p.add_argument("--input-samples", required=True, help="File: one input sample ID per line (Iowa/hybrids).")
    p.add_argument("--parent1-samples", required=True, help="File: parent1 sample IDs per line (Scat).")
    p.add_argument("--parent2-samples", required=True, help="File: parent2 sample IDs per line (Ster).")

    p.add_argument("--window", type=int, default=50_000)
    p.add_argument("--step", type=int, default=10_000)
    p.add_argument("--min-aims", type=int, default=10)

    p.add_argument("--diagnostic-threshold", type=float, default=0.05)
    p.add_argument("--min-parent-called", type=int, default=1)

    p.add_argument("--n-reps", type=int, default=100)
    p.add_argument("--parent1-subsample", type=int, default=0,
                   help="0 => auto=min(n_parent1,n_parent2) when --n-reps>1; else use all when --n-reps==1.")
    p.add_argument("--seed", type=int, default=1)

    p.add_argument("--ancestry-low", type=float, default=0.1)
    p.add_argument("--ancestry-high", type=float, default=0.9)

    p.add_argument("--chunk-windows", type=int, default=5000)
    p.add_argument("--output", required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()

    input_names = read_sample_list(args.input_samples)
    p1_names = read_sample_list(args.parent1_samples)
    p2_names = read_sample_list(args.parent2_samples)

    if len(input_names) == 0:
        raise SystemExit("ERROR: --input-samples list is empty.")
    if len(p1_names) == 0:
        raise SystemExit("ERROR: --parent1-samples list is empty.")
    if len(p2_names) == 0:
        raise SystemExit("ERROR: --parent2-samples list is empty.")

    # Read VCF header to map sample indices
    with open_textmaybe_gzip(args.vcf) as f:
        header_samples: Optional[List[str]] = None
        for line in f:
            if line.startswith("#CHROM"):
                header_samples = line.rstrip("\n").split("\t")[9:]
                break
    if header_samples is None:
        raise SystemExit("ERROR: VCF header '#CHROM' line not found.")

    name_to_idx: Dict[str, int] = {s: i for i, s in enumerate(header_samples)}

    def map_indices(wanted: Sequence[str], label: str) -> List[int]:
        missing = [s for s in wanted if s not in name_to_idx]
        if missing:
            raise SystemExit(f"ERROR: {label}: {len(missing)} sample(s) not found in VCF header. Example: {missing[:5]}")
        return [name_to_idx[s] for s in wanted]

    input_idx = map_indices(input_names, "input")
    p1_idx = map_indices(p1_names, "parent1")
    p2_idx = map_indices(p2_names, "parent2")

    n_p1 = len(p1_idx)
    n_p2 = len(p2_idx)

    if args.n_reps <= 1:
        subsample = args.parent1_subsample if args.parent1_subsample > 0 else 0
        n_reps = 1
    else:
        subsample = args.parent1_subsample if args.parent1_subsample > 0 else min(n_p1, n_p2)
        n_reps = args.n_reps

    if subsample > 0 and subsample > n_p1:
        raise SystemExit(f"ERROR: parent1-subsample={subsample} > n_parent1={n_p1}.")
    if args.min_parent_called < 1:
        raise SystemExit("ERROR: --min-parent-called must be >= 1.")

    rng = np.random.default_rng(args.seed)

    p1_subsets: List[np.ndarray] = []
    if n_reps == 1 and subsample == 0:
        p1_subsets = [np.arange(n_p1, dtype=np.int32)]
        ci_method = "support"
    else:
        if subsample == 0:
            subsample = n_p1
        for _ in range(n_reps):
            cols = rng.choice(n_p1, size=subsample, replace=False).astype(np.int32)
            p1_subsets.append(cols)
        ci_method = "perm"

    eprint("=== Hybrid index settings ===")
    eprint(f"VCF: {args.vcf}")
    eprint(f"n_input={len(input_names)}  n_parent1={n_p1}  n_parent2={n_p2}")
    eprint(f"window={args.window}  step={args.step}  min_aims={args.min_aims}")
    eprint(f"diagnostic_threshold={args.diagnostic_threshold}")
    eprint(f"min_parent_called={args.min_parent_called}")
    eprint(f"n_reps={len(p1_subsets)}  parent1_subsample={len(p1_subsets[0]) if p1_subsets else 0}  seed={args.seed}")
    eprint(f"CI method: {'support (2 logL units)' if ci_method=='support' else 'permutation percentiles'}")
    eprint("=============================")

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as out:
        out.write(
            "sample\tchrom\twindow_start\twindow_end\tmean_n_loci\tmean_hybrid_index\tci_low\tci_high\tancestry\tn_reps_used\n"
        )

        for cd in vcf_chrom_iterator(args.vcf, input_idx, p1_idx, p2_idx):
            p2_freq, p2_called = compute_alt_freq_and_called(cd.g_p2)

            process_chrom(
                cd=cd,
                input_names=input_names,
                p2_freq=p2_freq,
                p2_called=p2_called,
                p1_subsets=p1_subsets,
                diagnostic_threshold=args.diagnostic_threshold,
                min_parent_called=args.min_parent_called,
                window_size=args.window,
                step_size=args.step,
                min_aims=args.min_aims,
                chunk_windows=args.chunk_windows,
                ancestry_low=args.ancestry_low,
                ancestry_high=args.ancestry_high,
                out_handle=out,
                ci_method=ci_method,
            )

    eprint(f"Done. Wrote: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
