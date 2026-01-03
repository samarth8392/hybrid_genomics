#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=parse_vcf
#SBATCH -e %x
#SBATCH -o %x

MAINDIR="/fs/ess/scratch/PAS1533/smathur/hybrid_project/"

# Extract scaffolds from each VCF

cd $MAINDIR/vcf

module load vcftools/0.1.16

for vcf in allsistrurus.*.vcf.gz; do
    echo "=== $vcf ==="
    vcftools --gzvcf "$vcf" --freq --stdout 2>/dev/null | awk 'NR>1 {print $1}' | sort -u
    echo
done > scaffolds_per_vcf.txt

# Complete list of all unique scaffolds
for vcf in allsistrurus.*.vcf.gz; do
    vcftools --gzvcf "$vcf" --freq --stdout 2>/dev/null | awk 'NR>1 {print $1}' | sort -u
done | sort -u > all_scaffolds.txt
