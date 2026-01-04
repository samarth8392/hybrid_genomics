#!/usr/bin/env bash
# hybrid_index_sliding_window.sh
#
# Pipeline:
#  1) Filter to input+parents, enforce biallelic SNPs, apply missingness filter.
#  2) Compute parental FST (Weir-Cockerham).
#  3) Keep AIMs by FST threshold.
#  4) Subset to AIMs for all samples (input+parents).
#  5) Run calculate_hybrid_index_ml.py (diagnostic sites + sliding windows + optional parent1 subsampling perms).
#
# Critical fixes vs old script:
#  - Global MAF filter OFF by default (prevents dropping parent2-fixed alleles that are rare overall).
#  - Explicit biallelic SNP enforcement.
#  - No fragile paste/join of frequency files.
#  - Correct ancestry summary column.
#
set -euo pipefail

WINDOW_SIZE=50000
STEP_SIZE=10000
MIN_AIMS=10
FST_THRESHOLD=0.5
MISSING_THRESHOLD=0.2          # fraction missing allowed
GLOBAL_MIN_MAF=0.0             # default OFF

DIAG_THRESHOLD=0.05
MIN_PARENT_CALLED=1

N_REPS=100
P1_SUBSAMPLE=0                 # 0 => auto=min(n_parent1,n_parent2) when N_REPS>1

CHUNK_WINDOWS=5000
SEED=1

OUTDIR="hybrid_index_output"

usage() {
  cat << EOF
Usage:
  $0 --vcf FILE --input FILE --parent1 FILE --parent2 FILE [options]

Required:
  --vcf FILE
  --input FILE
  --parent1 FILE
  --parent2 FILE

Options:
  --window INT                [${WINDOW_SIZE}]
  --step INT                  [${STEP_SIZE}]
  --min-aims INT              [${MIN_AIMS}]
  --fst-threshold FLOAT       [${FST_THRESHOLD}]

  --max-missing FLOAT         fraction missing allowed (0.2 => require 0.8 called) [${MISSING_THRESHOLD}]
  --global-min-maf FLOAT      global MAF filter across ALL kept samples (OFF=0.0) [${GLOBAL_MIN_MAF}]

  --diagnostic-threshold FLOAT [${DIAG_THRESHOLD}]
  --min-parent-called INT      [${MIN_PARENT_CALLED}]

  --n-reps INT                [${N_REPS}]
  --parent1-subsample INT     (0=auto) [${P1_SUBSAMPLE}]
  --seed INT                  [${SEED}]

  --chunk-windows INT         [${CHUNK_WINDOWS}]
  --outdir DIR                [${OUTDIR}]
EOF
  exit 1
}

VCF=""
INPUT=""
PARENT1=""
PARENT2=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) VCF="$2"; shift 2 ;;
    --input) INPUT="$2"; shift 2 ;;
    --parent1) PARENT1="$2"; shift 2 ;;
    --parent2) PARENT2="$2"; shift 2 ;;

    --window) WINDOW_SIZE="$2"; shift 2 ;;
    --step) STEP_SIZE="$2"; shift 2 ;;
    --min-aims) MIN_AIMS="$2"; shift 2 ;;
    --fst-threshold) FST_THRESHOLD="$2"; shift 2 ;;

    --max-missing) MISSING_THRESHOLD="$2"; shift 2 ;;
    --global-min-maf) GLOBAL_MIN_MAF="$2"; shift 2 ;;

    --diagnostic-threshold) DIAG_THRESHOLD="$2"; shift 2 ;;
    --min-parent-called) MIN_PARENT_CALLED="$2"; shift 2 ;;

    --n-reps) N_REPS="$2"; shift 2 ;;
    --parent1-subsample) P1_SUBSAMPLE="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;

    --chunk-windows) CHUNK_WINDOWS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage ;;
  esac
done

if [[ -z "${VCF}" || -z "${INPUT}" || -z "${PARENT1}" || -z "${PARENT2}" ]]; then
  echo "ERROR: Missing required arguments." >&2
  usage
fi

# Optional module loads (HPC)
if command -v module >/dev/null 2>&1; then
  module load vcftools/0.1.16 >/dev/null 2>&1 || true
  module load python/3.12 >/dev/null 2>&1 || true
fi

for cmd in vcftools awk python3; do
  command -v "${cmd}" >/dev/null 2>&1 || { echo "ERROR: '${cmd}' not found in PATH." >&2; exit 1; }
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="${SCRIPT_DIR}/calculate_hybrid_index_ml.py"
[[ -f "${PY_SCRIPT}" ]] || { echo "ERROR: ${PY_SCRIPT} not found." >&2; exit 1; }

mkdir -p "${OUTDIR}"

VCFTOOLS_FLAG="--vcf"
[[ "${VCF}" == *.gz ]] && VCFTOOLS_FLAG="--gzvcf"

cat "${INPUT}" "${PARENT1}" "${PARENT2}" | awk 'NF{print $1}' | sort -u > "${OUTDIR}/all_samples.txt"

N_INPUT=$(wc -l < "${INPUT}" | awk '{print $1}')
N_P1=$(wc -l < "${PARENT1}" | awk '{print $1}')
N_P2=$(wc -l < "${PARENT2}" | awk '{print $1}')

echo "[$(date)] Starting hybrid index pipeline"
echo "VCF: ${VCF}"
echo "Input samples:   ${N_INPUT}"
echo "Parent1 samples: ${N_P1}"
echo "Parent2 samples: ${N_P2}"
echo "Output dir: ${OUTDIR}"

# 1) Filter VCF
echo "[$(date)] Filtering VCF (biallelic SNPs + missingness; global MAF=${GLOBAL_MIN_MAF})..."
MAX_MISSING_FRAC=$(awk -v m="${MISSING_THRESHOLD}" 'BEGIN{printf "%.6f", (1.0 - m)}')

FILTER_ARGS=(
  "${VCFTOOLS_FLAG}" "${VCF}"
  --keep "${OUTDIR}/all_samples.txt"
  --remove-indels
  --min-alleles 2
  --max-alleles 2
  --max-missing "${MAX_MISSING_FRAC}"
  --recode
  --recode-INFO-all
  --out "${OUTDIR}/filtered"
)

if awk -v x="${GLOBAL_MIN_MAF}" 'BEGIN{exit !(x>0)}'; then
  FILTER_ARGS+=( --maf "${GLOBAL_MIN_MAF}" )
fi

vcftools "${FILTER_ARGS[@]}"

if command -v bgzip >/dev/null 2>&1; then
  bgzip -f "${OUTDIR}/filtered.recode.vcf"
  mv "${OUTDIR}/filtered.recode.vcf.gz" "${OUTDIR}/filtered.vcf.gz"
  command -v tabix >/dev/null 2>&1 && tabix -f -p vcf "${OUTDIR}/filtered.vcf.gz" || true
else
  gzip -f "${OUTDIR}/filtered.recode.vcf"
  mv "${OUTDIR}/filtered.recode.vcf.gz" "${OUTDIR}/filtered.vcf.gz"
fi

# 2) FST
echo "[$(date)] Calculating Weir-Cockerham FST between parents..."
vcftools --gzvcf "${OUTDIR}/filtered.vcf.gz" \
  --weir-fst-pop "${PARENT1}" \
  --weir-fst-pop "${PARENT2}" \
  --out "${OUTDIR}/parental_fst"

# 3) AIMs
echo "[$(date)] Identifying AIMs (FST >= ${FST_THRESHOLD})..."
awk -v t="${FST_THRESHOLD}" 'NR>1 && $3 != "-nan" && $3 >= t {print $1"\t"$2}' \
  "${OUTDIR}/parental_fst.weir.fst" > "${OUTDIR}/aims.txt"

N_AIMS=$(wc -l < "${OUTDIR}/aims.txt" | awk '{print $1}')
echo "AIMs retained: ${N_AIMS}"
[[ "${N_AIMS}" -gt 0 ]] || { echo "ERROR: No AIMs found at this FST threshold." >&2; exit 1; }

# 4) Subset to AIMs for ALL samples
echo "[$(date)] Subsetting filtered VCF to AIMs for all samples..."
vcftools --gzvcf "${OUTDIR}/filtered.vcf.gz" \
  --positions "${OUTDIR}/aims.txt" \
  --keep "${OUTDIR}/all_samples.txt" \
  --recode \
  --recode-INFO-all \
  --out "${OUTDIR}/aims_all"

if command -v bgzip >/dev/null 2>&1; then
  bgzip -f "${OUTDIR}/aims_all.recode.vcf"
  mv "${OUTDIR}/aims_all.recode.vcf.gz" "${OUTDIR}/aims_all.vcf.gz"
  command -v tabix >/dev/null 2>&1 && tabix -f -p vcf "${OUTDIR}/aims_all.vcf.gz" || true
else
  gzip -f "${OUTDIR}/aims_all.recode.vcf"
  mv "${OUTDIR}/aims_all.recode.vcf.gz" "${OUTDIR}/aims_all.vcf.gz"
fi

# 5) Hybrid index windows
echo "[$(date)] Estimating hybrid index in sliding windows..."
python3 "${PY_SCRIPT}" \
  --vcf "${OUTDIR}/aims_all.vcf.gz" \
  --input-samples "${INPUT}" \
  --parent1-samples "${PARENT1}" \
  --parent2-samples "${PARENT2}" \
  --window "${WINDOW_SIZE}" \
  --step "${STEP_SIZE}" \
  --min-aims "${MIN_AIMS}" \
  --diagnostic-threshold "${DIAG_THRESHOLD}" \
  --min-parent-called "${MIN_PARENT_CALLED}" \
  --n-reps "${N_REPS}" \
  --parent1-subsample "${P1_SUBSAMPLE}" \
  --seed "${SEED}" \
  --chunk-windows "${CHUNK_WINDOWS}" \
  --output "${OUTDIR}/hybrid_index_windows.tsv"

echo "[$(date)] Done."
echo "Results: ${OUTDIR}/hybrid_index_windows.tsv"

# Summary (ancestry is column 9)
echo "[$(date)] Writing ancestry summary..."
awk 'NR>1 {a=$9; c[a]++; tot++} END{
  print "Ancestry\tWindows\tProportion"
  for(k in c){ printf "%s\t%d\t%.6f\n", k, c[k], c[k]/tot }
}' "${OUTDIR}/hybrid_index_windows.tsv" > "${OUTDIR}/ancestry_summary.tsv"

echo "Ancestry distribution written to: ${OUTDIR}/ancestry_summary.tsv"
