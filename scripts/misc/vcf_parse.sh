#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 10:00:00
#SBATCH --mem=30G
#SBATCH --job-name=parse_vcf7
#SBATCH -e %x
#SBATCH -o %x

MAINDIR="/fs/ess/scratch/PAS1533/smathur/hybrid_project"

# Extract scaffolds from each VCF

cd $MAINDIR/vcf

# module load vcftools/0.1.16

# for vcf in allsistrurus.*.vcf.gz; do
#     echo "=== $vcf ==="
#     vcftools --gzvcf "$vcf" --freq --stdout 2>/dev/null | awk 'NR>1 {print $1}' | sort -u
#     echo
# done > scaffolds_per_vcf.txt

# # Complete list of all unique scaffolds
# for vcf in allsistrurus.*.vcf.gz; do
#     vcftools --gzvcf "$vcf" --freq --stdout 2>/dev/null | awk 'NR>1 {print $1}' | sort -u
# done | sort -u > all_scaffolds.txt

# Extract individual chromosomes/scaffolds and create new organized VCFs

# Extract chromosomes from MHC file
# vcftools --gzvcf allsistrurus.MHC.minDP4.noinf.norm.vcf.gz --chr CM078116.1 --recode --recode-INFO-all --stdout | gzip -c > allsistrurus.CM078116.1.minDP4.noinf.norm.vcf.gz

# vcftools --gzvcf allsistrurus.MHC.minDP4.noinf.norm.vcf.gz --chr JBAIFV010000021.1 --recode --recode-INFO-all --stdout | gzip -c > allsistrurus.JBAIFV010000021.1.minDP4.noinf.norm.vcf.gz

# vcftools --gzvcf allsistrurus.MHC.minDP4.noinf.norm.vcf.gz --chr JBAIFV010000077.1 --recode --recode-INFO-all --stdout | gzip -c > allsistrurus.JBAIFV010000077.1.minDP4.noinf.norm.vcf.gz

# # Extract chromosomes from venom file
# vcftools --gzvcf allsistrurus.venom.minDP4.noinf.norm.vcf.gz --chr CM078122.1 --recode --recode-INFO-all --stdout | gzip -c > allsistrurus.CM078122.1.minDP4.noinf.norm.vcf.gz

# vcftools --gzvcf allsistrurus.venom.minDP4.noinf.norm.vcf.gz --chr CM078123.1 --recode --recode-INFO-all --stdout | gzip -c > allsistrurus.CM078123.1.minDP4.noinf.norm.vcf.gz

# vcftools --gzvcf allsistrurus.venom.minDP4.noinf.norm.vcf.gz --chr CM078131.1 --recode --recode-INFO-all --stdout | gzip -c > allsistrurus.CM078131.1.minDP4.noinf.norm.vcf.gz

# vcftools --gzvcf allsistrurus.unplaced.minDP4.noinf.norm.vcf.gz --chr CM078157.1 --recode --recode-INFO-all --stdout | gzip -c > allsistrurus.CM078157.1.minDP4.noinf.norm.vcf.gz


module load miniconda3/24.1.2-py310
source activate bcftools_env
module load vcftools/0.1.16

vcftools --gzvcf allsistrurus.MHC.minDP4.noinf.norm.vcf.gz \
  --chr JBAIFV010000077.1 --recode --recode-INFO-all --stdout | \
  bgzip -c > allsistrurus.JBAIFV010000077.1.minDP4.noinf.norm.vcf.gz


bcftools index -t allsistrurus.JBAIFV010000077.1.minDP4.noinf.norm.vcf.gz
bcftools index -t allsistrurus.unplaced.minDP4.noinf.norm.vcf.gz

bcftools concat -a \
  allsistrurus.unplaced.minDP4.noinf.norm.vcf.gz \
  allsistrurus.JBAIFV010000077.1.minDP4.noinf.norm.vcf.gz \
  -O z -o allsistrurus.unplaced.merged.vcf.gz

bcftools sort allsistrurus.unplaced.merged.vcf.gz \
  -O z -o allsistrurus.unplaced.sorted.vcf.gz

# mv allsistrurus.unplaced.sorted.vcf.gz allsistrurus.unplaced.minDP4.noinf.norm.vcf.gz