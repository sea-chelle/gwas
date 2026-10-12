
#!/bin/bash

# generate chromosome-specific S-LDSC annotations and LD scores
# from genomic BED files
#
# usage:
#   bash make_ldsc_annotations.sh DATA_DIR [ANNOTATION ...]
#
# examples:
#   bash make_ldsc_annotations.sh ~/Scratch/gwas
#   bash make_ldsc_annotations.sh ~/Scratch/gwas ATAC_Consensus
#   bash make_ldsc_annotations.sh ~/Scratch/gwas ATAC_Overall ATAC_Male
#
# automatically identifies all available annotations if none are specified.
# BED files must contain genomic intervals in hg19/GRCh37.
# mouse ATAC-seq coordinates must first be converted using liftOver.

set -euo pipefail

# check arguments
if [ "$#" -lt 1 ]; then
    echo "Usage: bash make_ldsc_annotations.sh DATA_DIR [ANNOTATION ...]"
    exit 1
fi

# define directories
DATA_DIR="$1"
shift

ANNOT_DIR="${ANNOTATION_DIR:-${DATA_DIR}/ldsc_annotations}"
LDSC="${DATA_DIR}/ldsc"
REFERENCE="${DATA_DIR}/reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC"
HM3="${DATA_DIR}/reference/ldsc/hm3_no_MHC.list.txt"

# activate ldsc conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ldsc

# check required files and directories
if [ ! -d "${ANNOT_DIR}" ]; then
    echo "ERROR: annotation directory not found: ${ANNOT_DIR}"
    exit 1
fi

if [ ! -f "${LDSC}/make_annot.py" ]; then
    echo "ERROR: make_annot.py not found"
    exit 1
fi

if [ ! -f "${LDSC}/ldsc.py" ]; then
    echo "ERROR: ldsc.py not found"
    exit 1
fi

if [ ! -f "${HM3}" ]; then
    echo "ERROR: HapMap3 SNP list not found: ${HM3}"
    exit 1
fi

# identify BED files when no annotations are specified
if [ "$#" -eq 0 ]; then

    mapfile -t BED_FILES < <(
        find "${ANNOT_DIR}" \
            -maxdepth 1 \
            -type f \
            \( -name "*_GRCh37_100kb.bed" -o -name "*_hg19.bed" \) |
        sort
    )

else

    BED_FILES=()

    for ANNOTATION in "$@"; do

        if [ -f "${ANNOT_DIR}/${ANNOTATION}_GRCh37_100kb.bed" ]; then
            BED_FILES+=("${ANNOT_DIR}/${ANNOTATION}_GRCh37_100kb.bed")

        elif [ -f "${ANNOT_DIR}/${ANNOTATION}_hg19.bed" ]; then
            BED_FILES+=("${ANNOT_DIR}/${ANNOTATION}_hg19.bed")

        else
            echo "ERROR: BED file not found for ${ANNOTATION}"
            exit 1
        fi

    done

fi

# check annotations were found
if [ "${#BED_FILES[@]}" -eq 0 ]; then
    echo "ERROR: no eligible BED files found"
    exit 1
fi

# report annotations
echo
echo "Found ${#BED_FILES[@]} annotations:"
printf '%s\n' "${BED_FILES[@]}"
echo

# process each BED file
for BED_FILE in "${BED_FILES[@]}"; do

    FILENAME="$(basename "${BED_FILE}")"

    ANNOTATION="${FILENAME%_GRCh37_100kb.bed}"
    ANNOTATION="${ANNOTATION%_hg19.bed}"

    echo
    echo "========================================"
    echo "Preparing ${ANNOTATION}"
    echo "========================================"

    # generate chromosome-specific annotations
    for chr in {1..22}; do

        ANNOT_FILE="${ANNOT_DIR}/${ANNOTATION}.${chr}.annot.gz"
        BIM_FILE="${REFERENCE}.${chr}.bim"

        if [ ! -f "${BIM_FILE}" ]; then
            echo "ERROR: missing reference BIM file: ${BIM_FILE}"
            exit 1
        fi

        if [ -f "${ANNOT_FILE}" ]; then
            echo "Skipping existing annotation: ${ANNOTATION} chr${chr}"
            continue
        fi

        echo "Generating annotation: ${ANNOTATION} chr${chr}"

        python "${LDSC}/make_annot.py" \
            --bed-file "${BED_FILE}" \
            --bimfile "${BIM_FILE}" \
            --annot-file "${ANNOT_FILE}"

    done

    # calculate chromosome-specific LD scores
    for chr in {1..22}; do

        ANNOT_FILE="${ANNOT_DIR}/${ANNOTATION}.${chr}.annot.gz"
        OUTPUT_PREFIX="${ANNOT_DIR}/${ANNOTATION}.${chr}"

        if [ -f "${OUTPUT_PREFIX}.l2.ldscore.gz" ] &&
           [ -f "${OUTPUT_PREFIX}.l2.M" ] &&
           [ -f "${OUTPUT_PREFIX}.l2.M_5_50" ]; then

            echo "Skipping existing LD scores: ${ANNOTATION} chr${chr}"
            continue
        fi

        echo "Calculating LD scores: ${ANNOTATION} chr${chr}"

        python "${LDSC}/ldsc.py" \
            --l2 \
            --bfile "${REFERENCE}.${chr}" \
            --annot "${ANNOT_FILE}" \
            --thin-annot \
            --ld-wind-cm 1 \
            --print-snps "${HM3}" \
            --out "${OUTPUT_PREFIX}"

    done

done

echo
echo "S-LDSC annotation preparation complete"
