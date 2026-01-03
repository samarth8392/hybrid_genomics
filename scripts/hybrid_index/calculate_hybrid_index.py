#!/usr/bin/env python3
"""
calculate_hybrid_index.py

Calculate hybrid index in sliding windows from VCF with diagnostic AIMs.
"""

import sys
import gzip
import argparse
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description='Calculate hybrid index in sliding windows'
    )
    parser.add_argument('--vcf', required=True,
                       help='VCF file with genotypes at diagnostic AIMs')
    parser.add_argument('--diagnostic', required=True,
                       help='File with diagnostic alleles (chrom, pos, ref, alt, parent2_allele)')
    parser.add_argument('--window', type=int, default=50000,
                       help='Window size in bp [50000]')
    parser.add_argument('--step', type=int, default=10000,
                       help='Step size in bp [10000]')
    parser.add_argument('--min-aims', type=int, default=10,
                       help='Minimum AIMs per window [10]')
    parser.add_argument('--output', required=True,
                       help='Output file')
    return parser.parse_args()


def load_diagnostic_alleles(diagnostic_file):
    """Load diagnostic alleles from file."""
    diagnostic = {}
    with open(diagnostic_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            chrom, pos, ref, alt, parent2_allele = parts
            diagnostic[(chrom, int(pos))] = parent2_allele
    return diagnostic


def parse_vcf(vcf_file):
    """Parse VCF file and extract genotypes by chromosome."""
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
            ref = parts[3]
            alt = parts[4]
            gts = parts[9:]
            
            genotypes_by_chrom[chrom].append((pos, ref, alt, gts))
    
    return samples, genotypes_by_chrom


def calculate_hybrid_index(genotype, ref, alt, parent2_allele):
    """
    Calculate contribution of parent2 alleles for a single genotype.
    
    Returns: (parent2_count, total_alleles)
    """
    if genotype.startswith('.'):
        return 0, 0
    
    # Parse genotype
    alleles = genotype.replace('|', '/').split('/')
    
    p2_count = 0
    total = 0
    
    for allele_idx in alleles:
        if allele_idx == '.':
            continue
        
        total += 1
        
        if allele_idx == '0':
            observed_allele = ref
        elif allele_idx == '1':
            observed_allele = alt.split(',')[0]  # Handle multi-allelic
        else:
            continue
        
        if observed_allele == parent2_allele:
            p2_count += 1
    
    return p2_count, total


def classify_ancestry(hyb_index, threshold_low=0.1, threshold_high=0.9):
    """Classify ancestry based on hybrid index."""
    if hyb_index < threshold_low:
        return 'parent1'
    elif hyb_index > threshold_high:
        return 'parent2'
    else:
        return 'admixed'


def process_windows(samples, genotypes_by_chrom, diagnostic, 
                   window_size, step_size, min_aims, output_file):
    """Calculate hybrid index in sliding windows."""
    
    with open(output_file, 'w') as out:
        # Write header
        out.write('sample\tchrom\twindow_start\twindow_end\t'
                 'n_aims\thyb_index\tancestry\n')
        
        # Process each chromosome
        for chrom in sorted(genotypes_by_chrom.keys()):
            variants = genotypes_by_chrom[chrom]
            
            if not variants:
                continue
            
            positions = [v[0] for v in variants]
            min_pos = min(positions)
            max_pos = max(positions)
            
            # Sliding windows
            window_start = min_pos
            while window_start < max_pos:
                window_end = window_start + window_size
                
                # Get variants in this window
                window_variants = [
                    v for v in variants 
                    if window_start <= v[0] < window_end
                ]
                
                if len(window_variants) < min_aims:
                    window_start += step_size
                    continue
                
                # Calculate hybrid index for each sample
                for s_idx, sample in enumerate(samples):
                    p2_allele_count = 0
                    total_alleles = 0
                    
                    for pos, ref, alt, gts in window_variants:
                        gt = gts[s_idx]
                        
                        parent2_allele = diagnostic.get((chrom, pos))
                        if not parent2_allele:
                            continue
                        
                        p2_count, total = calculate_hybrid_index(
                            gt, ref, alt, parent2_allele
                        )
                        
                        p2_allele_count += p2_count
                        total_alleles += total
                    
                    if total_alleles == 0:
                        continue
                    
                    hyb_index = p2_allele_count / total_alleles
                    ancestry = classify_ancestry(hyb_index)
                    
                    out.write(f'{sample}\t{chrom}\t{window_start}\t'
                             f'{window_end}\t{len(window_variants)}\t'
                             f'{hyb_index:.4f}\t{ancestry}\n')
                
                window_start += step_size


def main():
    args = parse_args()
    
    print(f"Loading diagnostic alleles from {args.diagnostic}...", 
          file=sys.stderr)
    diagnostic = load_diagnostic_alleles(args.diagnostic)
    print(f"Loaded {len(diagnostic)} diagnostic markers", file=sys.stderr)
    
    print(f"Parsing VCF file {args.vcf}...", file=sys.stderr)
    samples, genotypes_by_chrom = parse_vcf(args.vcf)
    print(f"Found {len(samples)} samples", file=sys.stderr)
    print(f"Found {len(genotypes_by_chrom)} chromosomes", file=sys.stderr)
    
    print(f"Calculating hybrid index in {args.window}bp windows "
          f"(step={args.step}bp)...", file=sys.stderr)
    process_windows(samples, genotypes_by_chrom, diagnostic,
                   args.window, args.step, args.min_aims, args.output)
    
    print(f"Results written to {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()