#!/bin/bash

# run MAGMA competitive gene-set association analysis
#
# usage:
#   bash gene_set_analysis.sh <data_directory> <trait> [trait ...]
#
# examples:
#   bash gene_set_analysis.sh ~/Scratch/gwas ASD
#   bash gene_set_analysis.sh ~/Scratch/gwas ASD ADHD ANX MDD PTSD SCZ BD

set -euo pipefail

# require a data directory and at least one trait
if [[ "$#" -lt 2 ]]; then
  echo "Usage: bash gene_set_analysis.sh <data_directory> <trait> [trait ...]"
  exit 1
fi

DATA_DIR="$1"
shift

# locate repository and manifest
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MANIFEST="${REPO_DIR}/manifests/GWAS_manifest.tsv"

# MAGMA executable and gene-set file
MAGMA="${DATA_DIR}/magma/magma"
GENE_SETS="${DATA_DIR}/gene_sets/stress_gene_sets.set.annot"

# check required files exist
if [[ ! -f "${MAGMA}" ]]; then
  echo "Error: MAGMA executable not found:"
  echo "${MAGMA}"
  exit 1
fi

if [[ ! -f "${GENE_SETS}" ]]; then
  echo "Error: gene-set annotation file not found:"
  echo "${GENE_SETS}"
  exit 1
fi

echo "Manifest: ${MANIFEST}"
echo "Data directory: ${DATA_DIR}"
echo "MAGMA: ${MAGMA}"
echo "Gene sets: ${GENE_SETS}"

# process each selected trait
for TRAIT in "$@"; do

  echo
  echo "========================================"
  echo "Running MAGMA gene-set analysis for ${TRAIT}"
  echo "========================================"

  # get manifest row for selected trait
  TRAIT_ROW=$(awk -F'\t' -v trait="${TRAIT}" '
    NR > 1 && $1 == trait {print; exit}
  ' "${MANIFEST}")

  if [[ -z "${TRAIT_ROW}" ]]; then
    echo "Error: trait '${TRAIT}' not found in manifest."
    exit 1
  fi

  # read fields required for file paths
  IFS=$'\t' read -r \
    TRAIT_NAME \
    STUDY \
    ANALYSIS_ID \
    FILE \
    REST <<< "${TRAIT_ROW}"

  ANALYSIS_DIR="${DATA_DIR}/$(dirname "${FILE}")"
  GENE_RESULTS="${ANALYSIS_DIR}/${ANALYSIS_ID}.genes.raw"
  OUTPUT_PREFIX="${ANALYSIS_DIR}/${ANALYSIS_ID}_gene_sets"

  # check gene-level results exist
  if [[ ! -f "${GENE_RESULTS}" ]]; then
    echo "Error: gene results file not found:"
    echo "${GENE_RESULTS}"
    exit 1
  fi

  echo "Gene results: ${GENE_RESULTS}"
  echo "Output prefix: ${OUTPUT_PREFIX}"

  "${MAGMA}" \
    --gene-results "${GENE_RESULTS}" \
    --set-annot "${GENE_SETS}" \
    --out "${OUTPUT_PREFIX}"

  echo
  echo "Gene-set analysis completed for ${TRAIT}."

done

echo
echo "Finished MAGMA gene-set analysis for selected GWAS."
