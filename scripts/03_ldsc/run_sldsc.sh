
#!/bin/bash

# run stratified LD score regression for selected GWAS traits and annotations
#
# models:
#   unconditional: baseline-LD + annotation
#   conditional:   baseline-LD + ATAC consensus + ATAC DAR annotation
#   all:           run both models
#
# usage:
#   bash run_sldsc.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR TRAIT [ANNOTATION] [--model MODEL]
#
# examples:
#   bash run_sldsc.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR ANX ATAC_Overall
#   bash run_sldsc.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR ANX ATAC_Male --model conditional
#   bash run_sldsc.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR ANX --model all
#
# automatically identifies available annotations if none are specified

set -euo pipefail

# check arguments
if [ "$#" -lt 5 ]; then
    echo "Usage: bash run_sldsc.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR TRAIT [ANNOTATION] [--model MODEL]"
    exit 1
fi

DATA_DIR="$1"
ANNOT_DIR="$2"
OUTPUT_DIR="$3"
REPO_DIR="$4"
TRAIT="$5"
shift 5

ANNOTATION="all"
MODEL="unconditional"

# parse optional arguments
while [ "$#" -gt 0 ]; do
    case "$1" in
        --model)
            if [ "$#" -lt 2 ]; then
                echo "ERROR: --model requires a value"
                exit 1
            fi
            MODEL="$2"
            shift 2
            ;;
        --annotation)
            if [ "$#" -lt 2 ]; then
                echo "ERROR: --annotation requires a value"
                exit 1
            fi
            ANNOTATION="$2"
            shift 2
            ;;
        -*)
            echo "ERROR: unknown option: $1"
            exit 1
            ;;
        *)
            ANNOTATION="$1"
            shift
            ;;
    esac
done

if [[ "${MODEL}" != "unconditional" &&
      "${MODEL}" != "conditional" &&
      "${MODEL}" != "all" ]]; then
    echo "ERROR: model must be unconditional, conditional or all"
    exit 1
fi

MANIFEST="${REPO_DIR}/manifests/GWAS_manifest.tsv"
LDSC="${DATA_DIR}/ldsc/ldsc.py"

BASELINE="${DATA_DIR}/reference/ldsc/baselineLD_v2.3/baselineLD."
WEIGHTS="${DATA_DIR}/reference/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
FRQ="${DATA_DIR}/reference/ldsc/1000G_Phase3_frq/1000G.EUR.QC."

BACKGROUND="ATAC_Consensus"
BACKGROUND_PREFIX="${ANNOT_DIR}/${BACKGROUND}."

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

# get analysis ID from manifest
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

if [ ! -f "${SUMSTATS}" ]; then
    echo "ERROR: munged summary statistics not found: ${SUMSTATS}"
    exit 1
fi

# identify available annotations
if [ "${ANNOTATION}" = "all" ]; then

    mapfile -t ANNOTATIONS < <(
        find "${ANNOT_DIR}" \
            -maxdepth 1 \
            -type f \
            -name "*.1.l2.ldscore.gz" \
            -printf "%f\n" |
        sed 's/\.1\.l2\.ldscore\.gz$//' |
        sort
    )

else
    ANNOTATIONS=("${ANNOTATION}")
fi

if [ "${#ANNOTATIONS[@]}" -eq 0 ]; then
    echo "ERROR: no annotations found"
    exit 1
fi

# check chromosome-specific annotation LD scores
check_annotation() {

    local annotation="$1"
    local prefix="${ANNOT_DIR}/${annotation}."

    for chr in {1..22}; do

        for suffix in l2.ldscore.gz l2.M l2.M_5_50; do

            if [ ! -f "${prefix}${chr}.${suffix}" ]; then
                echo "ERROR: missing annotation file: ${prefix}${chr}.${suffix}"
                exit 1
            fi

        done

    done
}

# run S-LDSC for one annotation and model
run_model() {

    local annotation="$1"
    local model="$2"

    local annot_prefix="${ANNOT_DIR}/${annotation}."
    local ref_ld_chr="${BASELINE},${annot_prefix}"
    local output_prefix="${OUTPUT_DIR}/${ANALYSIS_ID}_${annotation}"

    # include consensus background for conditional model
    if [ "${model}" = "conditional" ]; then

        check_annotation "${BACKGROUND}"

        ref_ld_chr="${BASELINE},${BACKGROUND_PREFIX},${annot_prefix}"
        output_prefix="${OUTPUT_DIR}/${ANALYSIS_ID}_${annotation}_conditional_${BACKGROUND}"

    fi

    check_annotation "${annotation}"

    echo
    echo "========================================"
    echo "Running S-LDSC"
    echo "Trait:      ${TRAIT}"
    echo "Analysis:   ${ANALYSIS_ID}"
    echo "Annotation: ${annotation}"
    echo "Model:      ${model}"
    echo "Output:     ${output_prefix}"
    echo "========================================"
    echo

    python -u "${LDSC}" \
        --h2 "${SUMSTATS}" \
        --ref-ld-chr "${ref_ld_chr}" \
        --w-ld-chr "${WEIGHTS}" \
        --overlap-annot \
        --frqfile-chr "${FRQ}" \
        --print-coefficients \
        --out "${output_prefix}"

    if [ ! -f "${output_prefix}.results" ]; then
        echo "ERROR: S-LDSC results file was not generated"
        exit 1
    fi

    echo "S-LDSC complete: ${output_prefix}.results"
}

# process selected annotations
for annotation in "${ANNOTATIONS[@]}"; do

    if [ "${MODEL}" = "unconditional" ] || [ "${MODEL}" = "all" ]; then
        run_model "${annotation}" "unconditional"
    fi

    if [ "${MODEL}" = "conditional" ] || [ "${MODEL}" = "all" ]; then

        # conditional analysis applies only to ATAC DAR annotations
        if [[ "${annotation}" == ATAC_* &&
              "${annotation}" != "${BACKGROUND}" ]]; then

            run_model "${annotation}" "conditional"

        fi

    fi

done

echo
echo "All requested S-LDSC analyses complete"
