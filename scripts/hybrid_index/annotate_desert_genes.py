#!/usr/bin/env python3
"""
annotate_desert_genes.py

Expand introgression desert regions into per-gene annotations.
Reads output from find_introgression_deserts.py and creates one row per gene
with its annotation from GFF.

Input: TSV from find_introgression_deserts.py with 'genes' column
Output: TSV with one gene per row and its annotation
"""

import argparse
import sys
from typing import Dict
from urllib.parse import unquote

import pandas as pd


def parse_gff_annotations(gff_path: str) -> Dict[str, str]:
    """Parse GFF and return dict of gene_id -> annotation."""
    annotations: Dict[str, str] = {}
    
    with open(gff_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            
            _chrom, _source, feature, _start, _end, _score, _strand, _phase, attrs = parts
            if feature != "gene":
                continue
            
            d: Dict[str, str] = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    d[k] = unquote(v)
            
            gene_id = d.get("Name") or d.get("gene") or d.get("ID") or "unknown"
            note = d.get("product") or d.get("Note") or ""
            
            annotations[gene_id] = note
    
    return annotations


def expand_desert_genes(
    deserts_df: pd.DataFrame,
    annotations: Dict[str, str],
) -> pd.DataFrame:
    """Expand comma-separated genes into individual rows with annotations."""
    
    rows = []
    
    for _, region in deserts_df.iterrows():
        genes_str = str(region.get("genes", ""))
        
        if not genes_str or genes_str == "":
            continue
        
        gene_ids = [g.strip() for g in genes_str.split(",") if g.strip()]
        
        for gene_id in gene_ids:
            rows.append({
                "chrom": region["chrom"],
                "region_start": int(region["region_start"]),
                "region_end": int(region["region_end"]),
                "region_length": int(region["length"]),
                "n_windows": int(region["n_windows"]),
                "gene_id": gene_id,
                "annotation": annotations.get(gene_id, ""),
            })
    
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Annotate genes found in introgression desert regions."
    )
    p.add_argument("--input", required=True, help="TSV from find_introgression_deserts.py")
    p.add_argument("--gff", required=True, help="GFF file with gene annotations")
    p.add_argument("--output", required=True, help="Output TSV with per-gene annotations")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    
    # Load desert regions
    deserts = pd.read_csv(args.input, sep="\t")
    
    if "genes" not in deserts.columns:
        print("ERROR: Input must contain 'genes' column. Run find_introgression_deserts.py with --gff", file=sys.stderr)
        sys.exit(1)
    
    # Parse GFF annotations
    annotations = parse_gff_annotations(args.gff)
    
    # Expand genes
    genes_df = expand_desert_genes(deserts, annotations)
    
    # Write output
    genes_df.to_csv(args.output, sep="\t", index=False)
    
    print(f"[OK] {len(genes_df)} genes annotated from {len(deserts)} regions -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()