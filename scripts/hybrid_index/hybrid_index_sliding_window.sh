#!/bin/bash

# hybrid_index_sliding_window.sh

set -e

# Defaults
WINDOW_SIZE=50000
STEP_SIZE=10000
MIN_AIMS=10
FST_THRESHOLD=0.5
MIN_MAF=0.05
MISSING_THRESHOLD=0.2
THREADS=4

usage() {
    cat << EOF
Usage: $0 --vcf <file> --input <list> --parent1 <list> --parent2 <list> [options]

Required:
  --vcf FILE           Input VCF file
  --input FILE         Sample list for hybrid index estimation
  --parent1 FILE       Parent population 1 (S. catenatus)
  --parent2 FILE       Parent population 2 (S. tergeminus)

Optional:
  --window INT         Window size [${WINDOW_SIZE}]
  --step INT           Step size [${STEP_SIZE}]
  --min-aims INT       Minimum AIMs per window [${MIN_AIMS}]
  --fst-threshold FLOAT FST threshold for AIMs [${FST_THRESHOLD}]
  --min-maf FLOAT      Minimum MAF [${MIN_MAF}]
  --max-missing FLOAT  Maximum missing data [${MISSING_THRESHOLD}]
  --threads INT        Number of threads [${THREADS}]
  --outdir DIR         Output directory [hybrid_index_output]

EOF
    exit 1
}

# Parse arguments
OUTDIR="hybrid_index_output"
while [[ $# -gt 0 ]]; do
    case $1 in
        --vcf) VCF="$2"; shift 2 ;;
        --input) INPUT="$2"; shift 2 ;;
        --parent1) PARENT1="$2"; shift 2 ;;
        --parent2) PARENT2="$2"; shift 2 ;;
        --window) WINDOW_SIZE="$2"; shift 2 ;;
        --step) STEP_SIZE="$2"; shift 2 ;;
        --min-aims) MIN_AIMS="$2"; shift 2 ;;
        --fst-threshold) FST_THRESHOLD="$2"; shift 2 ;;
        --min-maf) MIN_MAF="$2"; shift 2 ;;
        --max-missing) MISSING_THRESHOLD="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# Load modules

module load vcftools/0.1.16
module load python/3.12

# Check required arguments
if [[ -z "$VCF" || -z "$INPUT" || -z "$PARENT1" || -z "$PARENT2" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check dependencies
for cmd in vcftools awk python3; do
    if ! command -v $cmd &> /dev/null; then
        echo "Error: $cmd not found"
        exit 1
    fi
done

# Check for Python script
SCRIPT_DIR=$(dirname "$0")
if [[ ! -f "$SCRIPT_DIR/calculate_hybrid_index_ml.py" ]]; then
    echo "Error: calculate_hybrid_index.py not found in $SCRIPT_DIR"
    exit 1
fi



mkdir -p "$OUTDIR"

echo "[$(date)] Starting hybrid index analysis"
echo "VCF: $VCF"
echo "Input samples: $(wc -l < $INPUT)"
echo "Parent1 samples: $(wc -l < $PARENT1)"
echo "Parent2 samples: $(wc -l < $PARENT2)"

# Create combined sample list
cat $INPUT $PARENT1 $PARENT2 > "$OUTDIR/all_samples.txt"

# Step 1: Filter VCF
echo "[$(date)] Filtering VCF..."
MAX_MISSING_FRAC=$(echo "1 - $MISSING_THRESHOLD" | bc -l)

vcftools --gzvcf $VCF \
    --keep "$OUTDIR/all_samples.txt" \
    --maf $MIN_MAF \
    --max-missing $MAX_MISSING_FRAC \
    --remove-indels \
    --recode \
    --recode-INFO-all \
    --out "$OUTDIR/filtered"

# Compress filtered VCF
if command -v bgzip &> /dev/null; then
    bgzip -f "$OUTDIR/filtered.recode.vcf"
    mv "$OUTDIR/filtered.recode.vcf.gz" "$OUTDIR/filtered.vcf.gz"
else
    gzip -f "$OUTDIR/filtered.recode.vcf"
    mv "$OUTDIR/filtered.recode.vcf.gz" "$OUTDIR/filtered.vcf.gz"
fi

# Step 2: Calculate FST between parents
echo "[$(date)] Calculating FST between parental populations..."
vcftools --gzvcf "$OUTDIR/filtered.vcf.gz" \
    --weir-fst-pop $PARENT1 \
    --weir-fst-pop $PARENT2 \
    --out "$OUTDIR/parental_fst"

# Step 3: Identify AIMs
echo "[$(date)] Identifying ancestry-informative markers (FST > $FST_THRESHOLD)..."
awk -v thresh=$FST_THRESHOLD 'NR>1 && $3 != "-nan" && $3 >= thresh {print $1"\t"$2}' \
    "$OUTDIR/parental_fst.weir.fst" > "$OUTDIR/aims.txt"
N_AIMS=$(wc -l < "$OUTDIR/aims.txt")
echo "Identified $N_AIMS AIMs"

if [ $N_AIMS -eq 0 ]; then
    echo "Error: No AIMs identified. Lower --fst-threshold or check data."
    exit 1
fi

# Step 4: Extract parental alleles at AIMs
echo "[$(date)] Extracting parental allele frequencies..."

# Extract parent1 genotypes at AIMs
vcftools --gzvcf "$OUTDIR/filtered.vcf.gz" \
    --positions "$OUTDIR/aims.txt" \
    --keep $PARENT1 \
    --recode \
    --recode-INFO-all \
    --out "$OUTDIR/parent1_aims"

# Extract parent2 genotypes at AIMs  
vcftools --gzvcf "$OUTDIR/filtered.vcf.gz" \
    --positions "$OUTDIR/aims.txt" \
    --keep $PARENT2 \
    --recode \
    --recode-INFO-all \
    --out "$OUTDIR/parent2_aims"

# Calculate allele frequencies from VCF
awk '!/^#/ {
    chrom=$1; pos=$2; ref=$4; alt=$5;
    ref_count=0; alt_count=0; total=0;
    
    for(i=10; i<=NF; i++) {
        split($i, gt_fields, ":");
        gt = gt_fields[1];
        
        if(gt=="0/0" || gt=="0|0") {ref_count+=2; total+=2}
        else if(gt=="1/1" || gt=="1|1") {alt_count+=2; total+=2}
        else if(gt=="0/1" || gt=="1/0" || gt=="0|1" || gt=="1|0") {ref_count++; alt_count++; total+=2}
    }
    
    if(total>0) print chrom"\t"pos"\t"ref"\t"alt"\t"alt_count/total
}' "$OUTDIR/parent1_aims.recode.vcf" > "$OUTDIR/parent1_freqs.txt"

awk '!/^#/ {
    chrom=$1; pos=$2; ref=$4; alt=$5;
    ref_count=0; alt_count=0; total=0;
    
    for(i=10; i<=NF; i++) {
        split($i, gt_fields, ":");
        gt = gt_fields[1];
        
        if(gt=="0/0" || gt=="0|0") {ref_count+=2; total+=2}
        else if(gt=="1/1" || gt=="1|1") {alt_count+=2; total+=2}
        else if(gt=="0/1" || gt=="1/0" || gt=="0|1" || gt=="1|0") {ref_count++; alt_count++; total+=2}
    }
    
    if(total>0) print chrom"\t"pos"\t"ref"\t"alt"\t"alt_count/total
}' "$OUTDIR/parent2_aims.recode.vcf" > "$OUTDIR/parent2_freqs.txt"

# Step 5: Determine diagnostic alleles
echo "[$(date)] Determining diagnostic alleles for each AIM..."
paste "$OUTDIR/parent1_freqs.txt" "$OUTDIR/parent2_freqs.txt" | \
    awk '{
        chrom=$1; pos=$2; ref=$3; alt=$4; p1_freq=$5; p2_freq=$10;
        
        # Determine which allele is diagnostic for parent2
        if(p1_freq < 0.05 && p2_freq > 0.95) parent2_allele=alt
        else if(p1_freq > 0.95 && p2_freq < 0.05) parent2_allele=ref
        else parent2_allele="."
        
        if(parent2_allele != ".") print chrom"\t"pos"\t"ref"\t"alt"\t"parent2_allele
    }' > "$OUTDIR/diagnostic_alleles.txt"

N_DIAGNOSTIC=$(wc -l < "$OUTDIR/diagnostic_alleles.txt")
echo "Identified $N_DIAGNOSTIC diagnostic AIMs"

if [ $N_DIAGNOSTIC -eq 0 ]; then
    echo "Error: No diagnostic AIMs identified. Adjust thresholds or check data."
    exit 1
fi

# Step 6: Extract genotypes for input samples at diagnostic AIMs
echo "[$(date)] Extracting genotypes for input samples..."
awk '{print $1"\t"$2}' "$OUTDIR/diagnostic_alleles.txt" > "$OUTDIR/diagnostic_positions.txt"

vcftools --gzvcf "$OUTDIR/filtered.vcf.gz" \
    --positions "$OUTDIR/diagnostic_positions.txt" \
    --keep $INPUT \
    --recode \
    --recode-INFO-all \
    --out "$OUTDIR/input_aims"

# Compress if needed
if [ -f "$OUTDIR/input_aims.recode.vcf" ]; then
    if command -v bgzip &> /dev/null; then
        bgzip -f "$OUTDIR/input_aims.recode.vcf"
        mv "$OUTDIR/input_aims.recode.vcf.gz" "$OUTDIR/input_aims.vcf.gz"
    else
        gzip -f "$OUTDIR/input_aims.recode.vcf"
        mv "$OUTDIR/input_aims.recode.vcf.gz" "$OUTDIR/input_aims.vcf.gz"
    fi
fi

# Step 7: Calculate hybrid index in sliding windows
echo "[$(date)] Calculating hybrid index in sliding windows..."
python3 "$SCRIPT_DIR/calculate_hybrid_index_ml.py" \
    --vcf "$OUTDIR/input_aims.vcf.gz" \
    --parent1-freqs "$OUTDIR/parent1_freqs.txt" \
    --parent2-freqs "$OUTDIR/parent2_freqs.txt" \
    --window $WINDOW_SIZE \
    --step $STEP_SIZE \
    --min-aims $MIN_AIMS \
    --fixed \
    --output "$OUTDIR/hybrid_index_windows.txt"


echo "[$(date)] Analysis complete"
echo "Results: $OUTDIR/hybrid_index_windows.txt"

# Summary statistics
echo "[$(date)] Generating summary..."
awk 'NR>1 {ancestry[$7]++; total++} 
     END {
         print "Ancestry\tWindows\tProportion"
         for(a in ancestry) print a"\t"ancestry[a]"\t"ancestry[a]/total
     }' "$OUTDIR/hybrid_index_windows.txt" > "$OUTDIR/ancestry_summary.txt"

echo ""
echo "Ancestry distribution:"
cat "$OUTDIR/ancestry_summary.txt"

# Cleanup intermediate files
rm -f "$OUTDIR"/*.log
rm -f "$OUTDIR"/parent1_aims.recode.vcf "$OUTDIR"/parent2_aims.recode.vcf