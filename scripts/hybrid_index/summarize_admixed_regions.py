#!/usr/bin/env python3
"""
summarize_admixed_regions.py

Identify and summarize admixed/introgressed regions per sample.
"""

import argparse
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description='Summarize admixed regions per sample from hybrid index output'
    )
    parser.add_argument('--input', required=True,
                       help='hybrid_index_windows.txt file')
    parser.add_argument('--output', required=True,
                       help='Output file for admixed regions')
    parser.add_argument('--min-hybrid-index', type=float, default=0.1,
                       help='Minimum hybrid index for admixed [0.1]')
    parser.add_argument('--merge-distance', type=int, default=100000,
                       help='Merge regions within this distance [100000]')
    return parser.parse_args()


def load_windows(input_file, min_hybrid_index):
    """Load windows and filter for admixed/parent2 ancestry."""
    data_by_sample = defaultdict(lambda: defaultdict(list))
    
    with open(input_file) as f:
        header = f.readline()
        
        for line in f:
            parts = line.strip().split('\t')
            sample = parts[0]
            chrom = parts[1]
            window_start = int(parts[2])
            window_end = int(parts[3])
            n_aims = int(parts[4])
            hyb_index = float(parts[5])
            ancestry = parts[6]
            
            # Keep windows with hybrid index >= threshold
            if hyb_index >= min_hybrid_index:
                data_by_sample[sample][chrom].append({
                    'start': window_start,
                    'end': window_end,
                    'n_aims': n_aims,
                    'hyb_index': hyb_index,
                    'ancestry': ancestry
                })
    
    return data_by_sample


def merge_windows(windows, merge_distance):
    """Merge overlapping or nearby windows into regions."""
    if not windows:
        return []
    
    # Sort by start position
    windows = sorted(windows, key=lambda x: x['start'])
    
    regions = []
    current_region = {
        'start': windows[0]['start'],
        'end': windows[0]['end'],
        'n_windows': 1,
        'n_aims': windows[0]['n_aims'],
        'hyb_indices': [windows[0]['hyb_index']],
        'ancestries': [windows[0]['ancestry']]
    }
    
    for window in windows[1:]:
        # Check if window is within merge distance of current region
        if window['start'] <= current_region['end'] + merge_distance:
            # Extend current region
            current_region['end'] = max(current_region['end'], window['end'])
            current_region['n_windows'] += 1
            current_region['n_aims'] += window['n_aims']
            current_region['hyb_indices'].append(window['hyb_index'])
            current_region['ancestries'].append(window['ancestry'])
        else:
            # Save current region and start new one
            regions.append(current_region)
            current_region = {
                'start': window['start'],
                'end': window['end'],
                'n_windows': 1,
                'n_aims': window['n_aims'],
                'hyb_indices': [window['hyb_index']],
                'ancestries': [window['ancestry']]
            }
    
    # Add last region
    regions.append(current_region)
    
    return regions


def calculate_region_stats(region):
    """Calculate summary statistics for a region."""
    mean_hyb_index = sum(region['hyb_indices']) / len(region['hyb_indices'])
    min_hyb_index = min(region['hyb_indices'])
    max_hyb_index = max(region['hyb_indices'])
    
    # Determine dominant ancestry
    ancestry_counts = {}
    for anc in region['ancestries']:
        ancestry_counts[anc] = ancestry_counts.get(anc, 0) + 1
    
    dominant_ancestry = max(ancestry_counts.items(), key=lambda x: x[1])[0]
    
    # Check if region spans multiple ancestry types
    is_mixed = len(ancestry_counts) > 1
    
    return {
        'mean_hyb_index': mean_hyb_index,
        'min_hyb_index': min_hyb_index,
        'max_hyb_index': max_hyb_index,
        'dominant_ancestry': dominant_ancestry,
        'is_mixed': is_mixed,
        'ancestry_composition': ancestry_counts
    }


def main():
    args = parse_args()
    
    print(f"Loading windows from {args.input}...")
    data_by_sample = load_windows(args.input, args.min_hybrid_index)
    
    print(f"Found {len(data_by_sample)} samples with admixed regions")
    
    # Process each sample
    with open(args.output, 'w') as out:
        # Write header
        out.write('sample\tchrom\tregion_start\tregion_end\tregion_length\t'
                 'n_windows\tn_aims\tmean_hyb_index\tmin_hyb_index\tmax_hyb_index\t'
                 'dominant_ancestry\tis_mixed\n')
        
        total_regions = 0
        
        for sample in sorted(data_by_sample.keys()):
            sample_regions = 0
            
            for chrom in sorted(data_by_sample[sample].keys()):
                windows = data_by_sample[sample][chrom]
                regions = merge_windows(windows, args.merge_distance)
                
                for region in regions:
                    stats = calculate_region_stats(region)
                    
                    region_length = region['end'] - region['start']
                    
                    out.write(f"{sample}\t{chrom}\t{region['start']}\t{region['end']}\t"
                             f"{region_length}\t{region['n_windows']}\t{region['n_aims']}\t"
                             f"{stats['mean_hyb_index']:.4f}\t{stats['min_hyb_index']:.4f}\t"
                             f"{stats['max_hyb_index']:.4f}\t{stats['dominant_ancestry']}\t"
                             f"{stats['is_mixed']}\n")
                    
                    sample_regions += 1
                    total_regions += 1
            
            if sample_regions > 0:
                print(f"  {sample}: {sample_regions} admixed regions")
        
        print(f"\nTotal admixed regions identified: {total_regions}")
        print(f"Output written to {args.output}")
    
    # Generate per-sample summary
    summary_file = args.output.replace('.txt', '_summary.txt')
    with open(summary_file, 'w') as out:
        out.write('sample\tn_regions\ttotal_length\tmean_region_length\t'
                 'mean_hyb_index\tn_parent2_regions\tn_admixed_regions\n')
        
        for sample in sorted(data_by_sample.keys()):
            total_length = 0
            all_hyb_indices = []
            n_parent2 = 0
            n_admixed = 0
            n_regions = 0
            
            for chrom in data_by_sample[sample].keys():
                windows = data_by_sample[sample][chrom]
                regions = merge_windows(windows, args.merge_distance)
                
                for region in regions:
                    stats = calculate_region_stats(region)
                    region_length = region['end'] - region['start']
                    total_length += region_length
                    all_hyb_indices.extend(region['hyb_indices'])
                    
                    if stats['dominant_ancestry'] == 'parent2':
                        n_parent2 += 1
                    elif stats['dominant_ancestry'] == 'admixed':
                        n_admixed += 1
                    
                    n_regions += 1
            
            if n_regions > 0:
                mean_length = total_length / n_regions
                mean_hyb_index = sum(all_hyb_indices) / len(all_hyb_indices)
                
                out.write(f"{sample}\t{n_regions}\t{total_length}\t{mean_length:.0f}\t"
                         f"{mean_hyb_index:.4f}\t{n_parent2}\t{n_admixed}\n")
    
    print(f"Per-sample summary written to {summary_file}")


if __name__ == '__main__':
    main()