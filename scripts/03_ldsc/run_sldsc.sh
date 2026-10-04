#!/bin/bash

# run stratified ld score regression for selected gwas traits and annotations
#
# usage:
#   bash run_sldsc.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR TRAIT ANNOTATION
#
# example:
#   bash run_sldsc.sh \
#     ~/Scratch/gwas \
#     ~/Scratch/gwas/ldsc_annotations \
#     ~/Scratch/gwas/sldsc_results \
#     ~/Scratch/repos/gwas \
#     ANX \
#     ATAC_Overall

set -euo pipefail

if [ "$#" -ne 6 ]; then
    echo "Usage: bash run_sldsc.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR TRAIT ANNOTATION"
    exit 1
fi

DATA_DIR="$1"
ANNOT_DIR="$2"
OUTPUT_DIR="$3"
REPO_DIR="$4"
TRAIT="$5"
ANNOTATION="$6"

MANIFEST="${REPO_DIR}/manifests/GWAS_manifest.tsv"

LDSC="${DATA_DIR}/ldsc/ldsc.py"

BASELINE="${DATA_DIR}/reference/ldsc/baselineLD_v2.3/baselineLD."
WEIGHTS="${DATA_DIR}/reference/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
FRQ="${DATA_DIR}/reference/ldsc/1000G_Phase3_frq/1000G.EUR.QC."

# activate ldsc conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ldsc

# check required files and directories
if [ ! -f "${MANIFEST}" ]; then
    echo "ERROR: GWAS manifest not found: ${MANIFEST}"
    exit 1
fi

if [ ! -f "${LDSC}" ]; then
    echo "ERROR: ldsc.py not found: ${LDSC}"
    exit 1
fi

if [ ! -d "${ANNOT_DIR}" ]; then
    echo "ERROR: annotation directory not found: ${ANNOT_DIR}"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# get analysis id from manifest
ANALYSIS_ID=$(awk -F'\t' -v trait="${TRAIT}" '
    NR > 1 && $1 == trait {
        print $3
        exit
    }
' "${MANIFEST}")

if [ -z "${ANALYSIS_ID}" ]; then
    echo "ERROR: trait not found in GWAS manifest: ${TRAIT}"
    exit 1
fi

SUMSTATS="${DATA_DIR}/${ANALYSIS_ID}/${ANALYSIS_ID}.sumstats.gz"
ANNOT_PREFIX="${ANNOT_DIR}/${ANNOTATION}."
OUTPUT_PREFIX="${OUTPUT_DIR}/${ANALYSIS_ID}_${ANNOTATION}"

if [ ! -f "${SUMSTATS}" ]; then
    echo "ERROR: munged summary statistics not found: ${SUMSTATS}"
    exit 1
fi

# check annotation ld scores
for CHR in {1..22}; do
    if [ ! -f "${ANNOT_PREFIX}${CHR}.l2.ldscore.gz" ]; then
        echo "ERROR: missing annotation LD score: ${ANNOT_PREFIX}${CHR}.l2.ldscore.gz"
        exit 1
    fi

    if [ ! -f "${ANNOT_PREFIX}${CHR}.l2.M" ]; then
        echo "ERROR: missing annotation M file: ${ANNOT_PREFIX}${CHR}.l2.M"
        exit 1
    fi

    if [ ! -f "${ANNOT_PREFIX}${CHR}.l2.M_5_50" ]; then
        echo "ERROR: missing annotation M_5_50 file: ${ANNOT_PREFIX}${CHR}.l2.M_5_50"
        exit 1
    fi
done

echo "Running S-LDSC"
echo "Trait:      ${TRAIT}"
echo "Analysis:   ${ANALYSIS_ID}"
echo "Annotation: ${ANNOTATION}"
echo "Output:     ${OUTPUT_PREFIX}"
echo

python -u "${LDSC}" \
    --h2 "${SUMSTATS}" \
    --ref-ld-chr "${BASELINE},${ANNOT_PREFIX}" \
    --w-ld-chr "${WEIGHTS}" \
    --overlap-annot \
    --frqfile-chr "${FRQ}" \
    --print-coefficients \
    --out "${OUTPUT_PREFIX}"

if [ ! -f "${OUTPUT_PREFIX}.results" ]; then
    echo "ERROR: S-LDSC results file was not generated"
    exit 1
fi

echo
echo "S-LDSC complete"
echo "Results: ${OUTPUT_PREFIX}.results"
