#!/bin/bash
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=hybrid_index_sliding_window
#SBATCH -e %x
#SBATCH -o %x

# Usage: /fs/ess/scratch/PAS1533/smathur/hybrid_project/code/scripts/hybrid_index/hybrid_index_sliding_window.sh --vcf <file> --input <list> --parent1 <list> --parent2 <list> [options]

# Required:
#   --vcf FILE           Input VCF file
#   --input FILE         Sample list for hybrid index estimation
#   --parent1 FILE       Parent population 1 (S. catenatus)
#   --parent2 FILE       Parent population 2 (S. tergeminus)

# Optional:
#   --window INT         Window size [50000]
#   --step INT           Step size [10000]
#   --min-aims INT       Minimum AIMs per window [10]
#   --fst-threshold FLOAT FST threshold for AIMs [0.5]
#   --min-maf FLOAT      Minimum MAF [0.05]
#   --max-missing FLOAT  Maximum missing data [0.2]
#   --threads INT        Number of threads [4]
#   --outdir DIR         Output directory [hybrid_index_output]

MAINDIR="/fs/ess/scratch/PAS1533/smathur/hybrid_project"
OUTDIR="$MAINDIR/hybrid_index_output"
CODEDIR="$MAINDIR/code/scripts"

module load python/3.12

$CODEDIR/hybrid_index/hybrid_index_sliding_window.sh \
    --vcf $MAINDIR/vcf/allsistrurus.CM078129.1.minDP4.noinf.norm.vcf.gz \
    --input $CODEDIR/metadata/iowa_samples.txt \
    --parent1 $CODEDIR/metadata/scat_samples.txt \
    --parent2 $CODEDIR/metadata/ster_samples.txt \
    --window 50000 \
    --step 10000 \
    --min-aims 10 \
    --fst-threshold 0.9 \
    --min-maf 0.05 \
    --max-missing 0.2 \
    --threads 4 \
    --outdir $MAINDIR/hybrid_index_output

python $CODEDIR/hybrid_index/summarize_admixed_regions.py \
    --input $MAINDIR/hybrid_index_output/hybrid_index_windows.txt \
    --output $MAINDIR/hybrid_index_output/admixed_regions.txt \
    --min-hybrid-index 0.1 \
    --merge-distance 100000

# Basic usage - top 10 peaks with highest hybrid index
python $CODEDIR/hybrid_index/plot_hybrid_index.py \
    --input $MAINDIR/hybrid_index_output/hybrid_index_windows.txt \
    --output-summary $MAINDIR/hybrid_index_output/hybrid_index_summary.pdf \
    --output-individual $MAINDIR/hybrid_index_output/hybrid_index_individuals.pdf \
    --gff $MAINDIR/ref/Scate_genbankLiftOff.gff \
    --top-peaks 10 \
    --flank 10000 \
    --min-hybrid-index 0.15