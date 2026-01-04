#!/usr/bin/env python3
"""
summarize_admixed_regions.py

Identify admixed genomic regions using CI-based hybrid index classification
(Burke et al. 2005; Gompert et al. 2020).
"""

import argparse
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--parent1-threshold', type=float, default=0.1)
    p.add_argument('--parent2-threshold', type=float, default=0.9)
    p.add_argument('--merge-distance', type=int, default=0)
    return p.parse_args()

def classify(row, t1, t2):
    if row['upper_ci'] <= t1:
        return 'parent1'
    if row['lower_ci'] >= t2:
        return 'parent2'
    return 'admixed'

def main():
    args = parse_args()

    df = pd.read_csv(args.input, sep='\t')

    required = {'sample','chrom','window_start','window_end','lower_ci','upper_ci'}
    if not required.issubset(df.columns):
        raise ValueError(f"Missing required columns: {required - set(df.columns)}")

    df['ancestry'] = df.apply(
        classify, axis=1,
        t1=args.parent1_threshold,
        t2=args.parent2_threshold
    )

    admixed = df[df['ancestry'] == 'admixed']
    admixed = admixed.sort_values(['sample','chrom','window_start'])

    regions = []

    for (sample, chrom), g in admixed.groupby(['sample','chrom']):
        start = end = None

        for _, row in g.iterrows():
            if start is None:
                start = row.window_start
                end = row.window_end
            elif row.window_start <= end + args.merge_distance:
                end = max(end, row.window_end)
            else:
                regions.append((sample, chrom, start, end))
                start = row.window_start
                end = row.window_end

        if start is not None:
            regions.append((sample, chrom, start, end))

    out = pd.DataFrame(
        regions,
        columns=['sample','chrom','region_start','region_end']
    )
    out['length'] = out.region_end - out.region_start

    out.to_csv(args.output, sep='\t', index=False)
    print(f"Identified {len(out)} admixed regions")

if __name__ == '__main__':
    main()
