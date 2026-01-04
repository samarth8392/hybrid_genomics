#!/bin/bash


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
CODEDIR="$MAINDIR/code"

cat $CODEDIR/metadata/vcf_files.txt | while IFS=$'\t' read -r vcf
do
    basename=$(basename "$vcf" .vcf.gz)
    contig=${basename#*.}
    contig=${contig%%.*}

    cat > "$CODEDIR/scripts/per_sample/${vcf}.hybrid_index.sh" << EOF
#!/bin/bash
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 24:00:00
#SBATCH --job-name=hybrid_index_${contig}
#SBATCH -e %x
#SBATCH -o %x

module load python/3.12
module load vcftools/0.1.16

mkdir $MAINDIR/hybrid_index_output/${contig}

# $CODEDIR/scripts/hybrid_index/hybrid_index_sliding_window.sh \
#     --vcf $MAINDIR/vcf/${vcf} \
#     --input $CODEDIR/metadata/iowa_samples.txt \
#     --parent1 $CODEDIR/metadata/scat_samples.txt \
#     --parent2 $CODEDIR/metadata/ster_samples.txt \
#     --window 50000 \
#     --step 10000 \
#     --min-aims 10 \
#     --fst-threshold 0.5 \
#     --max-missing 0.2 \
#     --outdir $MAINDIR/hybrid_index_output/${contig}

python $CODEDIR/scripts/hybrid_index/summarize_admixed_regions.py \
    --input $MAINDIR/hybrid_index_output/${contig}/hybrid_index_windows.tsv \
    --output $MAINDIR/hybrid_index_output/${contig}/admixed_regions.tsv \
    --parent1-threshold 0.1 \
    --mode ci \
    --target admixed_or_parent2 \
    --merge-distance 100000 \
    --cohort-output admixed_regions_cohort.tsv \
    --cohort-min-samples 3
    
### Basic usage - top 10 peaks with highest hybrid index
python $CODEDIR/scripts/hybrid_index/plot_hybrid_index.py \
    --input $MAINDIR/hybrid_index_output/${contig}/hybrid_index_windows.tsv \
    --output-summary $MAINDIR/hybrid_index_output/${contig}/hybrid_index_summary.pdf \
    --output-individual $MAINDIR/hybrid_index_output/${contig}/hybrid_index_individuals.pdf \
    --output-yaml $MAINDIR/hybrid_index_output/${contig}/hybrid_peaks.yaml \
    --gff $MAINDIR/ref/Scate_genbankLiftOff.gff \
    --top-peaks 10 \
    --flank 10000 \
    --parent1-threshold 0.1 \
    --parent2-threshold 0.9

EOF
done

cat $CODEDIR/metadata/vcf_files.txt | while IFS=$'\t' read -r vcf
do
    cd /fs/ess/scratch/PAS1533/smathur/hybrid_project/errors
    sbatch $CODEDIR/scripts/per_sample/${vcf}.hybrid_index.sh
done


############

### Save final outputs
cat $CODEDIR/metadata/vcf_files.txt | while IFS=$'\t' read -r vcf
do
    basename=$(basename "$vcf" .vcf.gz)
    contig=${basename#*.}
    contig=${contig%%.*}
    
    cp $MAINDIR/hybrid_index_output/${vcf}/hybrid_index_individuals.pdf $MAINDIR/results/hybrid_index/$contig.hybrid_index_individuals.pdf
    cp $MAINDIR/hybrid_index_output/${vcf}/hybrid_index_summary.pdf $MAINDIR/results/hybrid_index/$contig.hybrid_index_summary.pdf
    cp $MAINDIR/hybrid_index_output/${vcf}/hybrid_peaks.yaml $MAINDIR/results/hybrid_index/$contig.hybrid_peaks.yaml
    cp $MAINDIR/hybrid_index_output/${vcf}/admixed_regions.txt $MAINDIR/results/hybrid_index/$contig.admixed_regions.txt
    cp $MAINDIR/hybrid_index_output/${vcf}/hybrid_index_windows.txt $MAINDIR/results/hybrid_index/$contig.hybrid_index_windows.txt
done