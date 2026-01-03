#!/usr/bin/env python3
"""
plot_hybrid_index_annotated.py

Visualize hybrid index across the genome with gene annotations at peaks:
1. Mean ± SD across all samples with peak annotations
2. Multi-page PDF with one sample per page with peak annotations
3. YAML summary of annotated peaks
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from collections import defaultdict
import re
import yaml


def parse_args():
    parser = argparse.ArgumentParser(
        description='Plot hybrid index across genome with peak annotations'
    )
    parser.add_argument('--input', required=True,
                       help='hybrid_index_windows.txt file')
    parser.add_argument('--output-summary', required=True,
                       help='Output summary plot (mean ± SD)')
    parser.add_argument('--output-individual', required=True,
                       help='Output multi-page PDF with individual samples')
    parser.add_argument('--output-yaml',
                       help='Output YAML file with peak annotations (optional, auto-generated if not specified)')
    parser.add_argument('--gff', 
                       help='GFF3 annotation file (optional)')
    parser.add_argument('--samples', nargs='+',
                       help='Specific samples to plot (default: all)')
    parser.add_argument('--chromosomes', nargs='+',
                       help='Specific chromosomes to plot (default: all)')
    parser.add_argument('--top-peaks', type=int, default=10,
                       help='Number of top peaks to annotate [10]')
    parser.add_argument('--flank', type=int, default=10000,
                       help='Distance to search for genes (bp) [10000]')
    parser.add_argument('--peak-method', 
                       choices=['highest_index', 'parent2_only', 'admixed_only'],
                       default='highest_index',
                       help='Method to identify peaks [highest_index]')
    parser.add_argument('--min-hybrid-index', type=float, default=0.5,
                       help='Minimum hybrid index for peak annotation [0.5]')
    parser.add_argument('--width', type=float, default=16,
                       help='Figure width in inches [16]')
    parser.add_argument('--height', type=float, default=6,
                       help='Figure height in inches [6]')
    parser.add_argument('--dpi', type=int, default=300,
                       help='DPI for output [300]')
    return parser.parse_args()


# Define high-contrast color scheme
PARENT1_COLOR = '#1E90FF'      # Dodger Blue
ADMIXED_COLOR = '#FF8C00'      # Dark Orange
PARENT2_COLOR = '#DC143C'      # Crimson


def load_data(input_file, samples=None, chromosomes=None):
    """Load hybrid index data."""
    df = pd.read_csv(input_file, sep='\t')
    
    if samples:
        df = df[df['sample'].isin(samples)]
    
    if chromosomes:
        df = df[df['chrom'].isin(chromosomes)]
    
    df['position'] = (df['window_start'] + df['window_end']) / 2
    
    return df


def parse_gff_attributes(attr_string):
    """Parse GFF attributes string into dictionary."""
    attr_dict = {}
    for attr in attr_string.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attr_dict[key.strip()] = value.strip()
    return attr_dict


def extract_note_description(note):
    """Extract clean description from Note field."""
    # Pattern: "Similar to <description> (<species> OX=<id>)"
    match = re.search(r'Similar to ([^(]+)\s*\([^)]+\)', note)
    if match:
        return match.group(1).strip()
    
    # If no match, return the whole note but clean it up
    return note.replace('Similar to ', '').split('(')[0].strip()


def parse_gff(gff_file):
    """Parse GFF file and extract genes with functional annotations."""
    print(f"Parsing GFF file: {gff_file}")
    
    features = defaultdict(list)
    
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            chrom = parts[0]
            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            attributes = parts[8]
            
            # Only keep genes
            if feature_type != 'gene':
                continue
            
            # Parse attributes
            attr_dict = parse_gff_attributes(attributes)
            
            # Extract gene information
            gene_id = attr_dict.get('Name', attr_dict.get('ID', 'unknown'))
            note = attr_dict.get('Note', '')
            
            # Extract functional description from Note
            if note:
                description = extract_note_description(note)
            else:
                description = ''
            
            features[chrom].append({
                'start': start,
                'end': end,
                'gene_id': gene_id,
                'description': description,
                'note': note,
                'midpoint': (start + end) / 2
            })
    
    total_genes = sum(len(features[chrom]) for chrom in features)
    print(f"Loaded {total_genes} genes from {len(features)} chromosomes")
    
    return features


def identify_peaks(df, method, top_n, min_hybrid_index):
    """Identify top peaks for annotation."""
    # Calculate mean hybrid index per window
    summary = df.groupby(['chrom', 'window_start', 'window_end', 'position'])['hyb_index'].agg(
        mean='mean',
        std='std'
    ).reset_index()
    
    # Filter by minimum
    summary = summary[summary['mean'] >= min_hybrid_index]
    
    if len(summary) == 0:
        return pd.DataFrame()
    
    if method == 'highest_index':
        peaks = summary.nlargest(top_n, 'mean')
    elif method == 'parent2_only':
        peaks = summary[summary['mean'] > 0.9].nlargest(top_n, 'mean')
    elif method == 'admixed_only':
        peaks = summary[(summary['mean'] >= 0.1) & (summary['mean'] <= 0.9)].nlargest(top_n, 'mean')
    
    return peaks.reset_index(drop=True)


def find_nearby_genes(peak_chrom, peak_start, peak_end, features, flank):
    """Find genes near a peak with their annotations."""
    if peak_chrom not in features:
        return []
    
    nearby = []
    for gene in features[peak_chrom]:
        gene_start = gene['start']
        gene_end = gene['end']
        
        if gene_end >= (peak_start - flank) and gene_start <= (peak_end + flank):
            # Calculate distance and position
            if gene_end < peak_start:
                distance = peak_start - gene_end
                position = 'upstream'
            elif gene_start > peak_end:
                distance = gene_start - peak_end
                position = 'downstream'
            else:
                distance = 0
                position = 'overlapping'
            
            nearby.append({
                'gene_id': gene['gene_id'],
                'description': gene['description'],
                'note': gene['note'],
                'start': gene_start,
                'end': gene_end,
                'distance': distance,
                'position': position
            })
    
    return nearby


def create_yaml_summary(peaks, features, flank, output_file):
    """Create YAML summary of annotated peaks."""
    print(f"\nCreating YAML summary: {output_file}")
    
    yaml_data = {
        'analysis_info': {
            'total_peaks': len(peaks),
            'flank_distance': flank,
            'timestamp': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')
        },
        'peaks': []
    }
    
    for idx, peak in peaks.iterrows():
        peak_info = {
            'rank': idx + 1,
            'location': {
                'chromosome': peak['chrom'],
                'start': int(peak['window_start']),
                'end': int(peak['window_end']),
                'length': int(peak['window_end'] - peak['window_start'])
            },
            'hybrid_index': {
                'mean': float(peak['mean']),
                'std': float(peak['std']) if not pd.isna(peak['std']) else 0.0
            },
            'genes': []
        }
        
        # Find nearby genes
        if features:
            genes = find_nearby_genes(peak['chrom'], peak['window_start'], 
                                     peak['window_end'], features, flank)
            
            for gene in genes:
                gene_info = {
                    'gene_id': gene['gene_id'],
                    'description': gene['description'] if gene['description'] else 'No description',
                    'full_note': gene['note'] if gene['note'] else 'N/A',
                    'location': {
                        'start': gene['start'],
                        'end': gene['end']
                    },
                    'distance_to_peak': gene['distance'],
                    'position_relative_to_peak': gene['position']
                }
                peak_info['genes'].append(gene_info)
        
        peak_info['gene_count'] = len(peak_info['genes'])
        yaml_data['peaks'].append(peak_info)
    
    # Write YAML file
    with open(output_file, 'w') as f:
        yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False, 
                 allow_unicode=True, width=100)
    
    print(f"Saved YAML summary with {len(peaks)} peaks")


def calculate_genome_positions(df):
    """Calculate cumulative genome positions for plotting."""
    chromosomes = sorted(df['chrom'].unique())
    
    chrom_offsets = {}
    cumulative_pos = 0
    chrom_boundaries = []
    
    for chrom in chromosomes:
        chrom_offsets[chrom] = cumulative_pos
        df_chrom = df[df['chrom'] == chrom]
        max_pos = df_chrom['position'].max()
        chrom_boundaries.append((chrom, cumulative_pos, cumulative_pos + max_pos))
        cumulative_pos += max_pos + 5e6
    
    df['genome_position'] = df.apply(
        lambda row: row['position'] + chrom_offsets[row['chrom']], axis=1
    )
    
    return df, chrom_boundaries, chrom_offsets


def format_gene_label(genes, max_genes=2):
    """Format gene information for plot labels."""
    if not genes:
        return None
    
    labels = []
    for gene in genes[:max_genes]:
        if gene['description']:
            # Shorten long descriptions
            desc = gene['description']
            if len(desc) > 30:
                desc = desc[:27] + '...'
            labels.append(f"{gene['gene_id']}: {desc}")
        else:
            labels.append(gene['gene_id'])
    
    result = '\n'.join(labels)
    if len(genes) > max_genes:
        result += f'\n(+{len(genes) - max_genes} more)'
    
    return result


def add_peak_annotations(ax, peaks, chrom_offsets, features, flank):
    """Add peak annotations to plot."""
    if peaks is None or len(peaks) == 0:
        return
    
    # Calculate genome positions for peaks
    peak_positions = []
    for idx, peak in peaks.iterrows():
        genome_pos = peak['position'] + chrom_offsets.get(peak['chrom'], 0)
        peak_positions.append({
            'genome_pos': genome_pos,
            'chrom': peak['chrom'],
            'start': peak['window_start'],
            'end': peak['window_end'],
            'hyb_index': peak['mean'],
            'rank': idx + 1
        })
    
    # Add vertical lines and annotations
    for i, peak_info in enumerate(peak_positions):
        genome_pos_mb = peak_info['genome_pos'] / 1e6
        
        # Add vertical line at peak
        ax.axvline(x=genome_pos_mb, color='red', linestyle=':', 
                  linewidth=1.5, alpha=0.7, zorder=4)
        
        # Find genes if features provided
        gene_label = f"Peak {peak_info['rank']}"
        if features:
            genes = find_nearby_genes(peak_info['chrom'], peak_info['start'], 
                                     peak_info['end'], features, flank)
            if genes:
                formatted_genes = format_gene_label(genes, max_genes=2)
                if formatted_genes:
                    gene_label = f"Peak {peak_info['rank']}\n{formatted_genes}"
        
        # Alternate annotation heights to avoid overlap
        y_pos = 0.95 if i % 2 == 0 else 0.80
        
        # Add annotation with arrow
        ax.annotate(gene_label,
                   xy=(genome_pos_mb, peak_info['hyb_index']),
                   xytext=(genome_pos_mb, y_pos),
                   fontsize=7,
                   ha='center',
                   va='bottom',
                   bbox=dict(boxstyle='round,pad=0.4', facecolor='yellow', 
                            alpha=0.8, edgecolor='red', linewidth=1.5),
                   arrowprops=dict(arrowstyle='->', color='red', lw=1.5, alpha=0.8),
                   zorder=5)


def plot_summary(df, chrom_boundaries, chrom_offsets, peaks, features, flank, output_file, args):
    """Plot mean ± SD hybrid index across all samples."""
    print("Creating summary plot (mean ± SD)...")
    
    # Calculate mean and SD per window
    summary = df.groupby(['chrom', 'window_start', 'window_end', 'genome_position'])['hyb_index'].agg(
        mean='mean',
        std='std',
        count='count'
    ).reset_index()
    
    summary['se'] = summary['std'] / np.sqrt(summary['count'])
    
    fig, ax = plt.subplots(figsize=(args.width, args.height))
    
    # Add shaded ancestry regions
    ax.axhspan(0, 0.1, alpha=0.15, color=PARENT1_COLOR, zorder=0)
    ax.axhspan(0.1, 0.9, alpha=0.15, color=ADMIXED_COLOR, zorder=0)
    ax.axhspan(0.9, 1.0, alpha=0.15, color=PARENT2_COLOR, zorder=0)
    
    # Sort by genome position
    summary = summary.sort_values('genome_position')
    
    # Plot mean line
    ax.plot(summary['genome_position'] / 1e6,
           summary['mean'],
           color='black',
           linewidth=2.5,
           label='Mean hybrid index',
           zorder=3)
    
    # Plot SD/SE
    ax.fill_between(summary['genome_position'] / 1e6,
                    summary['mean'] - summary['std'],
                    summary['mean'] + summary['std'],
                    alpha=0.25,
                    color='gray',
                    label='± 1 SD',
                    zorder=2)
    
    ax.fill_between(summary['genome_position'] / 1e6,
                    summary['mean'] - summary['se'],
                    summary['mean'] + summary['se'],
                    alpha=0.45,
                    color='dimgray',
                    label='± 1 SE',
                    zorder=2)
    
    # Add peak annotations
    if peaks is not None and len(peaks) > 0:
        add_peak_annotations(ax, peaks, chrom_offsets, features, flank)
    
    # Add chromosome boundaries
    for chrom, start, end in chrom_boundaries[:-1]:
        ax.axvline(x=end / 1e6, color='black', linestyle='-', linewidth=1.5, alpha=0.4)
    
    # Add chromosome labels
    for chrom, start, end in chrom_boundaries:
        mid_pos = (start + end) / 2 / 1e6
        chrom_label = chrom.replace('CM078', '').replace('.1', '')
        ax.text(mid_pos, -0.15, chrom_label,
               ha='center', va='top', fontsize=10, fontweight='bold')
    
    # Formatting
    ax.set_xlabel('Chromosome', fontsize=14, fontweight='bold')
    ax.set_ylabel('Hybrid Index', fontsize=14, fontweight='bold')
    
    title = f'Mean Hybrid Index Across {len(df["sample"].unique())} Samples'
    if peaks is not None and len(peaks) > 0:
        title += f' (Top {len(peaks)} Peaks Annotated)'
    ax.set_title(title, fontsize=16, fontweight='bold')
    
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    
    # Threshold lines
    ax.axhline(y=0.1, color='black', linestyle='--', linewidth=1.2, alpha=0.7)
    ax.axhline(y=0.9, color='black', linestyle='--', linewidth=1.2, alpha=0.7)
    
    # Legend
    legend_elements = [
        plt.Line2D([0], [0], color='black', linewidth=2.5, label='Mean'),
        mpatches.Patch(color='dimgray', alpha=0.45, label='± 1 SE'),
        mpatches.Patch(color='gray', alpha=0.25, label='± 1 SD'),
        mpatches.Patch(color=PARENT1_COLOR, alpha=0.4, label='Parent1 (S. catenatus)'),
        mpatches.Patch(color=ADMIXED_COLOR, alpha=0.4, label='Admixed'),
        mpatches.Patch(color=PARENT2_COLOR, alpha=0.4, label='Parent2 (S. tergeminus)')
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True, fontsize=10)
    
    plt.tight_layout()
    fig.savefig(output_file, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    
    print(f"Saved summary plot: {output_file}")


def plot_individual_sample(df_sample, sample, chrom_boundaries, chrom_offsets, 
                          peaks_sample, features, flank, ax, args):
    """Plot hybrid index for a single sample."""
    # Add shaded ancestry regions
    ax.axhspan(0, 0.1, alpha=0.15, color=PARENT1_COLOR, zorder=0)
    ax.axhspan(0.1, 0.9, alpha=0.15, color=ADMIXED_COLOR, zorder=0)
    ax.axhspan(0.9, 1.0, alpha=0.15, color=PARENT2_COLOR, zorder=0)
    
    # Sort by genome position
    df_sample = df_sample.sort_values('genome_position')
    
    # Determine line color
    mean_hyb_index = df_sample['hyb_index'].mean()
    if mean_hyb_index < 0.1:
        line_color = PARENT1_COLOR
        ancestry_label = 'Parent1-like'
    elif mean_hyb_index > 0.9:
        line_color = PARENT2_COLOR
        ancestry_label = 'Parent2-like'
    else:
        line_color = ADMIXED_COLOR
        ancestry_label = 'Admixed'
    
    # Plot hybrid index
    ax.plot(df_sample['genome_position'] / 1e6,
           df_sample['hyb_index'],
           color=line_color,
           linewidth=2.5,
           alpha=0.9,
           zorder=3)
    
    # Add peak annotations for this sample
    if peaks_sample is not None and len(peaks_sample) > 0:
        add_peak_annotations(ax, peaks_sample, chrom_offsets, features, flank)
    
    # Add chromosome boundaries
    for chrom, start, end in chrom_boundaries[:-1]:
        ax.axvline(x=end / 1e6, color='black', linestyle='-', linewidth=1.5, alpha=0.4)
    
    # Add chromosome labels
    for chrom, start, end in chrom_boundaries:
        mid_pos = (start + end) / 2 / 1e6
        chrom_label = chrom.replace('CM078', '').replace('.1', '')
        ax.text(mid_pos, -0.15, chrom_label,
               ha='center', va='top', fontsize=10, fontweight='bold')
    
    # Formatting
    ax.set_xlabel('Chromosome', fontsize=12, fontweight='bold')
    ax.set_ylabel('Hybrid Index', fontsize=12, fontweight='bold')
    
    title = f'{sample} (mean = {mean_hyb_index:.3f}, {ancestry_label})'
    if peaks_sample is not None and len(peaks_sample) > 0:
        title += f'\n(Top {len(peaks_sample)} Peaks Annotated)'
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    
    # Threshold lines
    ax.axhline(y=0.1, color='black', linestyle='--', linewidth=1.2, alpha=0.7)
    ax.axhline(y=0.9, color='black', linestyle='--', linewidth=1.2, alpha=0.7)
    
    # Legend
    legend_elements = [
        mpatches.Patch(color=PARENT1_COLOR, alpha=0.4, label='Parent1 (S. catenatus)'),
        mpatches.Patch(color=ADMIXED_COLOR, alpha=0.4, label='Admixed'),
        mpatches.Patch(color=PARENT2_COLOR, alpha=0.4, label='Parent2 (S. tergeminus)')
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True, fontsize=9)


def plot_individuals(df, chrom_boundaries, chrom_offsets, features, flank, output_file, args):
    """Create multi-page PDF with one sample per page."""
    print("Creating individual sample plots...")
    
    samples = sorted(df['sample'].unique())
    
    with PdfPages(output_file) as pdf:
        for idx, sample in enumerate(samples):
            print(f"  Plotting {sample} ({idx+1}/{len(samples)})...")
            
            df_sample = df[df['sample'] == sample]
            
            # Identify peaks for this sample
            peaks_sample = None
            if args.gff:
                df_sample_peaks = df_sample.copy()
                df_sample_peaks['mean'] = df_sample_peaks['hyb_index']
                df_sample_peaks['std'] = 0
                peaks_sample = identify_peaks(df_sample_peaks, args.peak_method, 
                                            args.top_peaks, args.min_hybrid_index)
            
            fig, ax = plt.subplots(figsize=(args.width, args.height))
            plot_individual_sample(df_sample, sample, chrom_boundaries, chrom_offsets,
                                 peaks_sample, features, flank, ax, args)
            plt.tight_layout()
            
            pdf.savefig(fig, dpi=args.dpi, bbox_inches='tight')
            plt.close(fig)
        
        # Add metadata
        d = pdf.infodict()
        d['Title'] = 'Hybrid Index Individual Samples'
        d['Author'] = 'hybrid_index_sliding_window'
        d['Subject'] = f'{len(samples)} samples'
    
    print(f"Saved individual plots: {output_file}")


def main():
    args = parse_args()
    
    print(f"Loading data from {args.input}...")
    df = load_data(args.input, args.samples, args.chromosomes)
    
    samples = sorted(df['sample'].unique())
    chromosomes = sorted(df['chrom'].unique())
    
    print(f"Processing {len(samples)} samples across {len(chromosomes)} chromosomes...")
    
    # Calculate genome-wide positions
    df, chrom_boundaries, chrom_offsets = calculate_genome_positions(df)
    
    # Parse GFF if provided
    features = None
    if args.gff:
        features = parse_gff(args.gff)
    
    # Identify peaks for summary plot
    peaks = None
    if args.gff:
        print(f"\nIdentifying top {args.top_peaks} peaks using {args.peak_method} method...")
        peaks = identify_peaks(df, args.peak_method, args.top_peaks, args.min_hybrid_index)
        if len(peaks) > 0:
            print(f"Identified {len(peaks)} peaks for annotation")
            
            # Create YAML summary
            yaml_output = args.output_yaml
            if not yaml_output:
                # Auto-generate YAML filename from summary plot name
                yaml_output = args.output_summary.replace('.pdf', '_peaks.yaml')
            create_yaml_summary(peaks, features, args.flank, yaml_output)
    
    # Create summary plot
    plot_summary(df, chrom_boundaries, chrom_offsets, peaks, features, 
                args.flank, args.output_summary, args)
    
    # Create individual plots
    plot_individuals(df, chrom_boundaries, chrom_offsets, features, 
                    args.flank, args.output_individual, args)
    
    print("\nDone!")
    print(f"  Summary plot: {args.output_summary}")
    print(f"  Individual plots: {args.output_individual}")
    if args.gff and peaks is not None and len(peaks) > 0:
        yaml_file = args.output_yaml if args.output_yaml else args.output_summary.replace('.pdf', '_peaks.yaml')
        print(f"  YAML annotation: {yaml_file}")


if __name__ == '__main__':
    main()