
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -l h_rt=2:30:00
#$ -l h_vmem=64G

# run S-LDSC trait x annotation x model combinations as an SGE array
#
# models:
#   unconditional: baseline-LD + individual annotation
#   conditional:   baseline-LD + ATAC consensus + ATAC DAR annotation
#   all:           run both models
#
# usage:
#   qsub -t 1-N run_sldsc_array.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR [--model MODEL]
#
# examples:
#   qsub -t 1-49 run_sldsc_array.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR
#   qsub -t 1-21 run_sldsc_array.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR --model conditional
#   qsub -t 1-70 run_sldsc_array.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR --model all

set -euo pipefail

# check arguments
if [ "$#" -lt 4 ]; then
    echo "Usage: qsub -t 1-N run_sldsc_array.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR [--model MODEL]"
    exit 1
fi

DATA_DIR="$1"
ANNOT_DIR="$2"
OUTPUT_DIR="$3"
REPO_DIR="$4"
shift 4

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
        *)
            echo "ERROR: unknown argument: $1"
            exit 1
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
RUN_SCRIPT="${REPO_DIR}/scripts/03_ldsc/run_sldsc.sh"

# check required files and directories
if [ ! -f "${MANIFEST}" ]; then
    echo "ERROR: GWAS manifest not found: ${MANIFEST}"
    exit 1
fi

if [ ! -f "${RUN_SCRIPT}" ]; then
    echo "ERROR: run_sldsc.sh not found: ${RUN_SCRIPT}"
    exit 1
fi

if [ ! -d "${ANNOT_DIR}" ]; then
    echo "ERROR: annotation directory not found: ${ANNOT_DIR}"
    exit 1
fi

# get traits from GWAS manifest
mapfile -t TRAITS < <(
    awk -F'\t' 'NR > 1 && $1 != "" {print $1}' "${MANIFEST}"
)

# identify available annotations from chromosome 1 LD scores
mapfile -t ANNOTATIONS < <(
    find "${ANNOT_DIR}" \
        -maxdepth 1 \
        -type f \
        -name "*.1.l2.ldscore.gz" \
        -printf "%f\n" |
    sed 's/\.1\.l2\.ldscore\.gz$//' |
    sort
)

# identify ATAC DAR annotations for conditional analyses
DAR_ANNOTATIONS=()

for annotation in "${ANNOTATIONS[@]}"; do
    if [[ "${annotation}" == ATAC_* &&
          "${annotation}" != "ATAC_Consensus" ]]; then
        DAR_ANNOTATIONS+=("${annotation}")
    fi
done

N_TRAITS=${#TRAITS[@]}
N_ANNOTATIONS=${#ANNOTATIONS[@]}
N_DARS=${#DAR_ANNOTATIONS[@]}

if [ "${N_TRAITS}" -eq 0 ]; then
    echo "ERROR: no GWAS traits found in manifest"
    exit 1
fi

# construct analysis combinations
JOB_ANNOTATIONS=()
JOB_MODELS=()

if [ "${MODEL}" = "unconditional" ] || [ "${MODEL}" = "all" ]; then

    for annotation in "${ANNOTATIONS[@]}"; do
        JOB_ANNOTATIONS+=("${annotation}")
        JOB_MODELS+=("unconditional")
    done

fi

if [ "${MODEL}" = "conditional" ] || [ "${MODEL}" = "all" ]; then

    for annotation in "${DAR_ANNOTATIONS[@]}"; do
        JOB_ANNOTATIONS+=("${annotation}")
        JOB_MODELS+=("conditional")
    done

fi

N_COMBINATIONS=${#JOB_ANNOTATIONS[@]}
N_JOBS=$((N_TRAITS * N_COMBINATIONS))

if [ "${N_COMBINATIONS}" -eq 0 ]; then
    echo "ERROR: no eligible annotation/model combinations found"
    exit 1
fi

# check SGE array task ID
if [ -z "${SGE_TASK_ID:-}" ]; then
    echo "ERROR: SGE_TASK_ID is not defined"
    echo "Submit this script using qsub -t 1-${N_JOBS}"
    exit 1
fi

if [ "${SGE_TASK_ID}" -lt 1 ] || [ "${SGE_TASK_ID}" -gt "${N_JOBS}" ]; then
    echo "ERROR: SGE_TASK_ID ${SGE_TASK_ID} is outside the valid range 1-${N_JOBS}"
    exit 1
fi

# convert array task ID to trait x annotation/model indices
INDEX=$((SGE_TASK_ID - 1))
TRAIT_INDEX=$((INDEX / N_COMBINATIONS))
COMBINATION_INDEX=$((INDEX % N_COMBINATIONS))

TRAIT="${TRAITS[$TRAIT_INDEX]}"
ANNOTATION="${JOB_ANNOTATIONS[$COMBINATION_INDEX]}"
TASK_MODEL="${JOB_MODELS[$COMBINATION_INDEX]}"

echo "S-LDSC array task ${SGE_TASK_ID}/${N_JOBS}"
echo "Model selection: ${MODEL}"
echo "Traits:          ${N_TRAITS}"
echo "Annotations:     ${N_ANNOTATIONS}"
echo "ATAC DARs:       ${N_DARS}"
echo "Trait:           ${TRAIT}"
echo "Annotation:      ${ANNOTATION}"
echo "Model:           ${TASK_MODEL}"
echo

# run selected S-LDSC analysis
bash "${RUN_SCRIPT}" \
    "${DATA_DIR}" \
    "${ANNOT_DIR}" \
    "${OUTPUT_DIR}" \
    "${REPO_DIR}" \
    "${TRAIT}" \
    "${ANNOTATION}" \
    --model "${TASK_MODEL}"
