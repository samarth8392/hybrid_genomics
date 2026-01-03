#!/usr/bin/env python3
"""
calculate_hybrid_index_ml.py

Maximum likelihood estimation of hybrid index following Buerkle (2005).
Supports both fixed differences (diagnostic markers) and shared alleles.
"""

import sys
import gzip
import argparse
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import comb
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description='ML estimation of hybrid index (Buerkle 2005)'
    )
    parser.add_argument('--vcf', required=True,
                       help='VCF with genotypes at AIMs')
    parser.add_argument('--parent1-freqs', required=True,
                       help='Parent1 allele frequencies (chrom, pos, ref, alt, alt_freq)')
    parser.add_argument('--parent2-freqs', required=True,
                       help='Parent2 allele frequencies (chrom, pos, ref, alt, alt_freq)')
    parser.add_argument('--window', type=int, default=50000,
                       help='Window size in bp [50000]')
    parser.add_argument('--step', type=int, default=10000,
                       help='Step size in bp [10000]')
    parser.add_argument('--min-aims', type=int, default=10,
                       help='Minimum AIMs per window [10]')
    parser.add_argument('--fixed', action='store_true',
                       help='Use fixed differences only (faster)')
    parser.add_argument('--diagnostic-threshold', type=float, default=0.05,
                       help='Threshold for diagnostic markers [0.05]')
    parser.add_argument('--output', required=True,
                       help='Output file')
    return parser.parse_args()


def load_allele_frequencies(freq_file):
    """Load allele frequencies from file."""
    freqs = {}
    with open(freq_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            chrom, pos, ref, alt, alt_freq = parts
            freqs[(chrom, int(pos))] = float(alt_freq)
    return freqs


def genotype_likelihood_fixed(n_p2_alleles, h):
    """
    Genotype likelihood for fixed differences (diagnostic markers).
    
    For diagnostic marker where p1=0, p2=1:
    P(n_p2 | h) = C(2, n_p2) * h^n_p2 * (1-h)^(2-n_p2)
    
    Args:
        n_p2_alleles: Number of parent2 alleles (0, 1, or 2)
        h: Hybrid index
    
    Returns:
        Log-likelihood
    """
    if n_p2_alleles == 0:
        return 2 * np.log(1 - h + 1e-10)
    elif n_p2_alleles == 1:
        return np.log(2) + np.log(h + 1e-10) + np.log(1 - h + 1e-10)
    else:  # n_p2_alleles == 2
        return 2 * np.log(h + 1e-10)


def genotype_likelihood_shared(n_alt_alleles, h, p1_alt, p2_alt):
    """
    Genotype likelihood for shared alleles.
    
    Expected ALT frequency in hybrid:
    p_hybrid = h * p2_alt + (1-h) * p1_alt
    
    P(n_alt | h) = C(2, n_alt) * p_hybrid^n_alt * (1-p_hybrid)^(2-n_alt)
    
    Args:
        n_alt_alleles: Number of ALT alleles observed (0, 1, or 2)
        h: Hybrid index
        p1_alt: ALT frequency in parent1
        p2_alt: ALT frequency in parent2
    
    Returns:
        Log-likelihood
    """
    # Expected ALT frequency in hybrid
    p_hybrid = h * p2_alt + (1 - h) * p1_alt
    p_hybrid = np.clip(p_hybrid, 1e-10, 1 - 1e-10)
    
    if n_alt_alleles == 0:
        return 2 * np.log(1 - p_hybrid)
    elif n_alt_alleles == 1:
        return np.log(2) + np.log(p_hybrid) + np.log(1 - p_hybrid)
    else:  # n_alt_alleles == 2
        return 2 * np.log(p_hybrid)


def parse_genotype(gt_string):
    """
    Parse genotype string to count ALT alleles.
    
    Returns:
        (n_alt_alleles, n_total_alleles) or (None, None) if missing
    """
    if gt_string.startswith('.'):
        return None, None
    
    gt = gt_string.split(':')[0]  # Extract GT field
    alleles = gt.replace('|', '/').split('/')
    
    n_alt = 0
    n_total = 0
    
    for allele in alleles:
        if allele == '.':
            continue
        n_total += 1
        if allele != '0':  # ALT allele
            n_alt += 1
    
    if n_total == 0:
        return None, None
    
    return n_alt, n_total


def negative_log_likelihood_fixed(h, genotype_data):
    """
    Negative log-likelihood for fixed differences approach.
    
    Args:
        h: Hybrid index (0 to 1)
        genotype_data: List of (n_p2_alleles,) tuples
    
    Returns:
        Negative log-likelihood
    """
    log_lik = 0
    for n_p2_alleles, in genotype_data:
        log_lik += genotype_likelihood_fixed(n_p2_alleles, h)
    
    return -log_lik


def negative_log_likelihood_shared(h, genotype_data):
    """
    Negative log-likelihood for shared alleles approach.
    
    Args:
        h: Hybrid index (0 to 1)
        genotype_data: List of (n_alt, p1_alt, p2_alt) tuples
    
    Returns:
        Negative log-likelihood
    """
    log_lik = 0
    for n_alt, p1_alt, p2_alt in genotype_data:
        log_lik += genotype_likelihood_shared(n_alt, h, p1_alt, p2_alt)
    
    return -log_lik


def estimate_h_ml(genotype_data, fixed=True):
    """
    Maximum likelihood estimation of hybrid index.
    
    Args:
        genotype_data: Genotype information for all loci
        fixed: Whether to use fixed differences model
    
    Returns:
        (h_ml, lower_ci, upper_ci)
    """
    if fixed:
        result = minimize_scalar(
            lambda h: negative_log_likelihood_fixed(h, genotype_data),
            bounds=(0, 1),
            method='bounded'
        )
    else:
        result = minimize_scalar(
            lambda h: negative_log_likelihood_shared(h, genotype_data),
            bounds=(0, 1),
            method='bounded'
        )
    
    h_ml = result.x
    max_log_lik = -result.fun
    
    # Calculate 95% CI (values within 2 log-likelihood units)
    threshold = max_log_lik - 2
    
    # Find lower bound
    lower_ci = 0
    for h_test in np.linspace(0, h_ml, 100):
        if fixed:
            ll = -negative_log_likelihood_fixed(h_test, genotype_data)
        else:
            ll = -negative_log_likelihood_shared(h_test, genotype_data)
        if ll >= threshold:
            lower_ci = h_test
            break
    
    # Find upper bound
    upper_ci = 1
    for h_test in np.linspace(h_ml, 1, 100):
        if fixed:
            ll = -negative_log_likelihood_fixed(h_test, genotype_data)
        else:
            ll = -negative_log_likelihood_shared(h_test, genotype_data)
        if ll < threshold:
            upper_ci = h_test
            break
    
    return h_ml, lower_ci, upper_ci


def classify_ancestry(h, threshold_low=0.1, threshold_high=0.9):
    """Classify ancestry based on hybrid index."""
    if h < threshold_low:
        return 'parent1'
    elif h > threshold_high:
        return 'parent2'
    else:
        return 'admixed'


def parse_vcf(vcf_file):
    """Parse VCF and extract genotypes by chromosome."""
    samples = []
    genotypes_by_chrom = defaultdict(list)
    
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'
    
    with open_func(vcf_file, mode) as f:
        for line in f:
            if line.startswith('##'):
                continue
            
            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
                continue
            
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            gts = parts[9:]
            
            genotypes_by_chrom[chrom].append((pos, gts))
    
    return samples, genotypes_by_chrom


def process_windows_fixed(samples, genotypes_by_chrom, p1_freqs, p2_freqs,
                          window_size, step_size, min_aims, diag_threshold, 
                          output_file):
    """
    Process windows using fixed differences (diagnostic markers only).
    This is equivalent to simple allele counting when markers are diagnostic.
    """
    with open(output_file, 'w') as out:
        out.write('sample\tchrom\twindow_start\twindow_end\t'
                 'n_aims\thyb_index\tlower_ci\tupper_ci\tancestry\n')
        
        for chrom in sorted(genotypes_by_chrom.keys()):
            variants = genotypes_by_chrom[chrom]
            
            if not variants:
                continue
            
            positions = [v[0] for v in variants]
            min_pos = min(positions)
            max_pos = max(positions)
            
            window_start = min_pos
            while window_start < max_pos:
                window_end = window_start + window_size
                
                window_variants = [
                    v for v in variants 
                    if window_start <= v[0] < window_end
                ]
                
                # Filter for diagnostic markers only
                diagnostic_variants = []
                for pos, gts in window_variants:
                    p1_freq = p1_freqs.get((chrom, pos))
                    p2_freq = p2_freqs.get((chrom, pos))
                    
                    if p1_freq is None or p2_freq is None:
                        continue
                    
                    # Check if diagnostic
                    if (p1_freq < diag_threshold and p2_freq > (1 - diag_threshold)) or \
                       (p1_freq > (1 - diag_threshold) and p2_freq < diag_threshold):
                        diagnostic_variants.append((pos, gts, p2_freq > 0.5))
                
                if len(diagnostic_variants) < min_aims:
                    window_start += step_size
                    continue
                
                # Process each sample
                for s_idx, sample in enumerate(samples):
                    genotype_data = []
                    
                    for pos, gts, p2_is_alt in diagnostic_variants:
                        n_alt, n_total = parse_genotype(gts[s_idx])
                        
                        if n_alt is None:
                            continue
                        
                        # Convert to parent2 allele count
                        if p2_is_alt:
                            n_p2 = n_alt
                        else:
                            n_p2 = n_total - n_alt
                        
                        genotype_data.append((n_p2,))
                    
                    if len(genotype_data) == 0:
                        continue
                    
                    # ML estimation
                    h_ml, lower_ci, upper_ci = estimate_h_ml(genotype_data, fixed=True)
                    ancestry = classify_ancestry(h_ml)
                    
                    out.write(f'{sample}\t{chrom}\t{window_start}\t'
                             f'{window_end}\t{len(genotype_data)}\t'
                             f'{h_ml:.4f}\t{lower_ci:.4f}\t{upper_ci:.4f}\t'
                             f'{ancestry}\n')
                
                window_start += step_size


def process_windows_shared(samples, genotypes_by_chrom, p1_freqs, p2_freqs,
                           window_size, step_size, min_aims, output_file):
    """
    Process windows using shared alleles approach (full ML).
    """
    with open(output_file, 'w') as out:
        out.write('sample\tchrom\twindow_start\twindow_end\t'
                 'n_aims\thyb_index\tlower_ci\tupper_ci\tancestry\n')
        
        for chrom in sorted(genotypes_by_chrom.keys()):
            variants = genotypes_by_chrom[chrom]
            
            if not variants:
                continue
            
            positions = [v[0] for v in variants]
            min_pos = min(positions)
            max_pos = max(positions)
            
            window_start = min_pos
            while window_start < max_pos:
                window_end = window_start + window_size
                
                window_variants = [
                    v for v in variants 
                    if window_start <= v[0] < window_end
                ]
                
                if len(window_variants) < min_aims:
                    window_start += step_size
                    continue
                
                # Process each sample
                for s_idx, sample in enumerate(samples):
                    genotype_data = []
                    
                    for pos, gts in window_variants:
                        p1_freq = p1_freqs.get((chrom, pos))
                        p2_freq = p2_freqs.get((chrom, pos))
                        
                        if p1_freq is None or p2_freq is None:
                            continue
                        
                        n_alt, n_total = parse_genotype(gts[s_idx])
                        
                        if n_alt is None or n_total != 2:
                            continue
                        
                        genotype_data.append((n_alt, p1_freq, p2_freq))
                    
                    if len(genotype_data) == 0:
                        continue
                    
                    # ML estimation with shared alleles
                    h_ml, lower_ci, upper_ci = estimate_h_ml(genotype_data, fixed=False)
                    ancestry = classify_ancestry(h_ml)
                    
                    out.write(f'{sample}\t{chrom}\t{window_start}\t'
                             f'{window_end}\t{len(genotype_data)}\t'
                             f'{h_ml:.4f}\t{lower_ci:.4f}\t{upper_ci:.4f}\t'
                             f'{ancestry}\n')
                
                window_start += step_size


def main():
    args = parse_args()
    
    print(f"Loading parent1 frequencies...", file=sys.stderr)
    p1_freqs = load_allele_frequencies(args.parent1_freqs)
    print(f"Loaded {len(p1_freqs)} sites", file=sys.stderr)
    
    print(f"Loading parent2 frequencies...", file=sys.stderr)
    p2_freqs = load_allele_frequencies(args.parent2_freqs)
    print(f"Loaded {len(p2_freqs)} sites", file=sys.stderr)
    
    print(f"Parsing VCF...", file=sys.stderr)
    samples, genotypes_by_chrom = parse_vcf(args.vcf)
    print(f"Found {len(samples)} samples", file=sys.stderr)
    
    print(f"Calculating ML hybrid index...", file=sys.stderr)
    
    if args.fixed:
        print("Using fixed differences (diagnostic markers only)", file=sys.stderr)
        process_windows_fixed(samples, genotypes_by_chrom, p1_freqs, p2_freqs,
                             args.window, args.step, args.min_aims, 
                             args.diagnostic_threshold, args.output)
    else:
        print("Using shared alleles (full ML)", file=sys.stderr)
        process_windows_shared(samples, genotypes_by_chrom, p1_freqs, p2_freqs,
                               args.window, args.step, args.min_aims, args.output)
    
    print(f"Results written to {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()