#!/usr/bin/env python3
"""
plot_hybrid_index_ci_annotated.py

CI-aware visualization of hybrid index across the genome with gene annotations.
Implements Burke et al. (2005) hybrid index interpretation.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from collections import defaultdict
import numpy as np
import yaml
import re

# ------------------ CONFIG ------------------

PARENT1_COLOR = '#1E90FF'
ADMIXED_COLOR = '#FFF3E0'
PARENT2_COLOR = '#FFB6C1'

# ------------------ ARGS ------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True)
    p.add_argument('--output-summary', required=True)
    p.add_argument('--output-individual', required=True)
    p.add_argument('--output-yaml')
    p.add_argument('--gff')
    p.add_argument('--top-peaks', type=int, default=10)
    p.add_argument('--flank', type=int, default=10000)
    p.add_argument('--parent1-threshold', type=float, default=0.1)
    p.add_argument('--parent2-threshold', type=float, default=0.9)
    p.add_argument('--dpi', type=int, default=300)
    return p.parse_args()

# ------------------ CORE LOGIC ------------------

def classify_ci(row, t1, t2):
    if row.upper_ci <= t1:
        return 'parent1'
    if row.lower_ci >= t2:
        return 'parent2'
    return 'admixed'

def load_data(path, t1, t2):
    df = pd.read_csv(path, sep='\t')
    df['position'] = (df.window_start + df.window_end) / 2
    df['ancestry'] = df.apply(classify_ci, axis=1, t1=t1, t2=t2)
    return df

# ------------------ GFF PARSING ------------------

def clean_gene_name(note):
    """Extract clean gene name from Note field."""
    if not note or note == '':
        return None
    
    match = re.search(r'Similar to (.+?)\s*\(', note)
    if match:
        name = match.group(1).strip()
        name = name.rstrip(':')
        return name
    
    return note

def parse_gff(gff):
    feats = defaultdict(list)
    with open(gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom, _, ftype, start, end = parts[0], parts[2], parts[2], int(parts[3]), int(parts[4])
            attrs = parts[8]
            if ftype != 'gene':
                continue
            attr = dict(x.split('=',1) for x in attrs.split(';') if '=' in x)
            gene_id = attr.get('Name', attr.get('ID', 'unknown'))
            note = attr.get('Note','')
            gene_name = clean_gene_name(note) or gene_id
            
            feats[chrom].append({
                'start': start,
                'end': end,
                'gene': gene_name,
                'gene_id': gene_id,
                'note': note
            })
    return feats

def nearby_genes(chrom, start, end, feats, flank):
    out = []
    for g in feats.get(chrom, []):
        if g['end'] >= start-flank and g['start'] <= end+flank:
            out.append(g)
    return out

# ------------------ PEAK LOGIC ------------------

def identify_peaks(df, top_n):
    summary = df.groupby(['chrom','window_start','window_end','position']).agg(
        mean=('hyb_index','mean'),
        std=('hyb_index','std'),
        lower_ci=('lower_ci','mean'),
        upper_ci=('upper_ci','mean')
    ).reset_index()
    summary = summary[summary.upper_ci > 0.1]
    return summary.nlargest(top_n, 'mean')

# ------------------ PLOTTING ------------------

def plot_summary(df, peaks, feats, args):
    fig, ax = plt.subplots(figsize=(16,6))
    
    # Background colored regions - fill entire plot height
    mean_df = df.groupby('position').agg({
        'hyb_index': 'mean',
        'lower_ci': 'mean',
        'upper_ci': 'mean'
    }).reset_index()
    
    # Sort by position
    mean_df = mean_df.sort_values('position')
    
    # Create continuous background
    prev_x = None
    for i, row in mean_df.iterrows():
        if row.upper_ci <= args.parent1_threshold:
            color = PARENT1_COLOR
        elif row.lower_ci >= args.parent2_threshold:
            color = PARENT2_COLOR
        else:
            color = ADMIXED_COLOR
        
        curr_x = row.position / 1e6
        
        if prev_x is None:
            x_start = curr_x
        else:
            x_start = prev_x
        
        if i < len(mean_df) - 1:
            x_end = mean_df.iloc[i+1].position / 1e6
            x_end = (curr_x + x_end) / 2
        else:
            x_end = curr_x
        
        ax.axvspan(x_start, x_end, ymin=0, ymax=1, 
                  alpha=0.25, color=color, linewidth=0, zorder=0)
        prev_x = x_end
    
    # Error band
    ax.fill_between(
        mean_df.position/1e6,
        mean_df.lower_ci,
        mean_df.upper_ci,
        color='gray',
        alpha=0.25,
        linewidth=0,
        zorder=2
    )
    
    # Mean line
    ax.plot(mean_df.position/1e6, mean_df.hyb_index, 
            color='black', lw=1.5, zorder=3)
    
    # Thresholds
    ax.axhline(args.parent1_threshold, ls='--', lw=1.5, color='#0D47A1', alpha=0.8, zorder=1)
    ax.axhline(args.parent2_threshold, ls='--', lw=1.5, color='#B71C1C', alpha=0.8, zorder=1)
    
    # Peak markers and annotations
    for i, p in peaks.iterrows():
        x = p.position/1e6
        ax.axvline(x, color='red', ls=':', lw=2, alpha=0.7, zorder=2)
        
        # Gene names only
        if feats:
            genes = nearby_genes(p.chrom, p.window_start, p.window_end, feats, args.flank)
            if genes:
                gene_list = [g['gene'] for g in genes[:2]]
                if len(genes) > 2:
                    gene_list.append(f"(+{len(genes)-2} more)")
                label = ", ".join(gene_list)
                
                # Stagger heights
                y_pos = 0.93 if i % 2 == 0 else 0.86
                
                ax.annotate(label, xy=(x, p['mean']), xytext=(x, y_pos),
                           textcoords='data',
                           arrowprops=dict(arrowstyle='->', lw=1.5, color='red'),
                           bbox=dict(boxstyle='round,pad=0.4', fc='yellow', 
                                    ec='orange', alpha=0.95, lw=1.5),
                           ha='center', fontsize=8, zorder=5)
    
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(mean_df.position.min()/1e6, mean_df.position.max()/1e6)
    ax.set_xlabel("Genomic position (Mb)", fontsize=12, fontweight='bold')
    ax.set_ylabel("Hybrid Index", fontsize=12, fontweight='bold')
    ax.set_title("Mean Hybrid Index Across 10 Samples (Top 10 Peaks Annotated)", 
                 fontsize=14, fontweight='bold')
    
    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=PARENT1_COLOR, label='Parent1 (S. catenatus)', alpha=0.5),
        mpatches.Patch(facecolor=ADMIXED_COLOR, label='Admixed', alpha=0.5),
        mpatches.Patch(facecolor=PARENT2_COLOR, label='Parent2 (S. tergeminus)', alpha=0.5),
        plt.Line2D([0], [0], color='gray', alpha=0.5, lw=8, label='± 1 SE'),
        plt.Line2D([0], [0], color='black', lw=2, label='Mean')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10, 
             framealpha=0.95, edgecolor='black')
    
    fig.tight_layout()
    fig.savefig(args.output_summary, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)

# ------------------ YAML ------------------

def write_yaml(peaks, feats, args):
    if not args.output_yaml:
        return
    out = {'peaks': []}
    for i, p in peaks.iterrows():
        entry = {
            'rank': i+1,
            'chrom': p.chrom,
            'start': int(p.window_start),
            'end': int(p.window_end),
            'mean': float(p['mean']),
            'genes': []
        }
        if feats:
            for g in nearby_genes(p.chrom, p.window_start, p.window_end, feats, args.flank):
                entry['genes'].append({
                    'gene': g['gene'],
                    'gene_id': g['gene_id'],
                    'note': g['note']
                })
        out['peaks'].append(entry)
    with open(args.output_yaml,'w') as f:
        yaml.dump(out, f, default_flow_style=False)

# ------------------ MAIN ------------------

def main():
    args = parse_args()
    df = load_data(args.input, args.parent1_threshold, args.parent2_threshold)
    feats = parse_gff(args.gff) if args.gff else None
    peaks = identify_peaks(df, args.top_peaks)
    
    plot_summary(df, peaks, feats, args)
    write_yaml(peaks, feats, args)
    
    print(f"Summary plot: {args.output_summary}")
    if args.output_yaml:
        print(f"Peak annotations: {args.output_yaml}")

if __name__ == '__main__':
    main()