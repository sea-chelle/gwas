#!/bin/bash

# run stratified ld score regression jointly for male and female annotations
#
# usage:
#   bash run_sldsc_sex_comparison.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR TRAIT MODALITY
#
# example:
#   bash run_sldsc_sex_comparison.sh \
#     ~/Scratch/gwas \
#     ~/Scratch/gwas/ldsc_annotations \
#     ~/Scratch/gwas/sldsc_results_sex \
#     ~/Scratch/repos/gwas \
#     MDD \
#     RNA

set -euo pipefail

if [ "$#" -ne 6 ]; then
    echo "Usage: bash run_sldsc_sex_comparison.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR TRAIT 
MODALITY"
    exit 1
fi

DATA_DIR="$1"
ANNOT_DIR="$2"
OUTPUT_DIR="$3"
REPO_DIR="$4"
TRAIT="$5"
MODALITY="$6"

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

# select male and female annotations
case "${MODALITY}" in
    RNA)
        MALE_ANNOTATION="RNA_Male_60"
        FEMALE_ANNOTATION="RNA_Female_60"
        ;;
    ATAC)
        MALE_ANNOTATION="ATAC_Male"
        FEMALE_ANNOTATION="ATAC_Female"
        ;;
    *)
        echo "ERROR: modality must be RNA or ATAC"
        exit 1
        ;;
esac

SUMSTATS="${DATA_DIR}/${ANALYSIS_ID}/${ANALYSIS_ID}.sumstats.gz"

MALE_PREFIX="${ANNOT_DIR}/${MALE_ANNOTATION}."
FEMALE_PREFIX="${ANNOT_DIR}/${FEMALE_ANNOTATION}."

OUTPUT_PREFIX="${OUTPUT_DIR}/${ANALYSIS_ID}_${MODALITY}_Male_Female"

if [ ! -f "${SUMSTATS}" ]; then
    echo "ERROR: munged summary statistics not found: ${SUMSTATS}"
    exit 1
fi

# check annotation ld scores
for CHR in {1..22}; do
    for ANNOT_PREFIX in "${MALE_PREFIX}" "${FEMALE_PREFIX}"; do

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
done

echo "Running S-LDSC sex comparison"
echo "Trait:             ${TRAIT}"
echo "Analysis:          ${ANALYSIS_ID}"
echo "Modality:          ${MODALITY}"
echo "Male annotation:   ${MALE_ANNOTATION}"
echo "Female annotation: ${FEMALE_ANNOTATION}"
echo "Output:            ${OUTPUT_PREFIX}"
echo

python -u "${LDSC}" \
    --h2 "${SUMSTATS}" \
    --ref-ld-chr "${BASELINE},${MALE_PREFIX},${FEMALE_PREFIX}" \
    --w-ld-chr "${WEIGHTS}" \
    --overlap-annot \
    --frqfile-chr "${FRQ}" \
    --print-coefficients \
    --print-cov \
    --out "${OUTPUT_PREFIX}"

if [ ! -f "${OUTPUT_PREFIX}.results" ]; then
    echo "ERROR: S-LDSC results file was not generated"
    exit 1
fi

echo
echo "S-LDSC sex comparison complete"
echo "Results: ${OUTPUT_PREFIX}.results"
