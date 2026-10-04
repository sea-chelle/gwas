#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -l h_rt=2:30:00
#$ -l h_vmem=64G

# run s-ldsc trait x annotation combinations as an sge array

set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage: qsub -t 1-N run_sldsc_array.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR"
    exit 1
fi

DATA_DIR="$1"
ANNOT_DIR="$2"
OUTPUT_DIR="$3"
REPO_DIR="$4"

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

# get traits from gwas manifest
mapfile -t TRAITS < <(
    awk -F'\t' 'NR > 1 && $1 != "" {print $1}' "${MANIFEST}"
)

# get annotation names from chromosome 1 ld score files
mapfile -t ANNOTATIONS < <(
    find "${ANNOT_DIR}" \
        -maxdepth 1 \
        -type f \
        -name "*.1.l2.ldscore.gz" \
        -printf "%f\n" |
    sed 's/\.1\.l2\.ldscore\.gz$//' |
    sort
)

N_TRAITS=${#TRAITS[@]}
N_ANNOTATIONS=${#ANNOTATIONS[@]}
N_JOBS=$((N_TRAITS * N_ANNOTATIONS))

if [ "${N_TRAITS}" -eq 0 ]; then
    echo "ERROR: no GWAS traits found in manifest"
    exit 1
fi

if [ "${N_ANNOTATIONS}" -eq 0 ]; then
    echo "ERROR: no annotations found in ${ANNOT_DIR}"
    exit 1
fi

if [ -z "${SGE_TASK_ID:-}" ]; then
    echo "ERROR: SGE_TASK_ID is not defined"
    echo "Submit this script using qsub -t 1-${N_JOBS}"
    exit 1
fi

if [ "${SGE_TASK_ID}" -lt 1 ] || [ "${SGE_TASK_ID}" -gt "${N_JOBS}" ]; then
    echo "ERROR: SGE_TASK_ID ${SGE_TASK_ID} is outside the valid range 1-${N_JOBS}"
    exit 1
fi

# convert array task id to trait x annotation indices
INDEX=$((SGE_TASK_ID - 1))
TRAIT_INDEX=$((INDEX / N_ANNOTATIONS))
ANNOTATION_INDEX=$((INDEX % N_ANNOTATIONS))

TRAIT="${TRAITS[$TRAIT_INDEX]}"
ANNOTATION="${ANNOTATIONS[$ANNOTATION_INDEX]}"

echo "S-LDSC array task ${SGE_TASK_ID}/${N_JOBS}"
echo "Traits:      ${N_TRAITS}"
echo "Annotations: ${N_ANNOTATIONS}"
echo "Trait:       ${TRAIT}"
echo "Annotation:  ${ANNOTATION}"
echo

bash "${RUN_SCRIPT}" \
    "${DATA_DIR}" \
    "${ANNOT_DIR}" \
    "${OUTPUT_DIR}" \
    "${REPO_DIR}" \
    "${TRAIT}" \
    "${ANNOTATION}"
