#!/bin/bash

#$ -cwd
#$ -V
#$ -l h_rt=2:30:00
#$ -l h_vmem=64G

# run s-ldsc male-female comparison models as an sge array
#
# usage:
#   qsub -t 1-14 run_sldsc_sex_comparison_array.sh DATA_DIR ANNOT_DIR OUTPUT_DIR REPO_DIR
#
# example:
#   qsub -t 1-14 run_sldsc_sex_comparison_array.sh \
#     ~/Scratch/gwas \
#     ~/Scratch/gwas/ldsc_annotations \
#     ~/Scratch/gwas/sldsc_results_sex \
#     ~/Scratch/repos/gwas

set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage: qsub -t 1-14 run_sldsc_sex_comparison_array.sh DATA_DIR ANNOT_DIR OUTPUT_DIR 
REPO_DIR"
    exit 1
fi

DATA_DIR="$1"
ANNOT_DIR="$2"
OUTPUT_DIR="$3"
REPO_DIR="$4"

WORKER="${REPO_DIR}/scripts/03_ldsc/run_sldsc_sex_comparison.sh"

if [ ! -f "${WORKER}" ]; then
    echo "ERROR: worker script not found: ${WORKER}"
    exit 1
fi

TRAITS=(
    ASD
    ADHD
    ANX
    MDD
    PTSD
    SCZ
    BD
)

MODALITIES=(
    RNA
    ATAC
)

# generate all trait x modality combinations
COMBINATIONS=()

for TRAIT in "${TRAITS[@]}"; do
    for MODALITY in "${MODALITIES[@]}"; do
        COMBINATIONS+=("${TRAIT} ${MODALITY}")
    done
done

N=${#COMBINATIONS[@]}

if [ "${SGE_TASK_ID}" -lt 1 ] || [ "${SGE_TASK_ID}" -gt "${N}" ]; then
    echo "ERROR: SGE_TASK_ID ${SGE_TASK_ID} is outside valid range 1-${N}"
    exit 1
fi

# select combination for this array task
read -r TRAIT MODALITY <<< "${COMBINATIONS[$((SGE_TASK_ID - 1))]}"

echo "S-LDSC sex comparison array task"
echo "Task:     ${SGE_TASK_ID}/${N}"
echo "Trait:    ${TRAIT}"
echo "Modality: ${MODALITY}"
echo

bash "${WORKER}" \
    "${DATA_DIR}" \
    "${ANNOT_DIR}" \
    "${OUTPUT_DIR}" \
    "${REPO_DIR}" \
    "${TRAIT}" \
    "${MODALITY}"
