#!/bin/bash

# perform MAGMA SNP-to-gene annotation
#
# usage:
#   bash annotate_genes.sh <data_directory> <trait> [trait ...]
#
# example:
#   bash annotate_genes.sh ~/Scratch/gwas ASD ADHD ANX

set -euo pipefail

# require a data directory and at least one trait
if [[ "$#" -lt 2 ]]; then
  echo "Usage: bash annotate_genes.sh <data_directory> <trait> [trait ...]"
  exit 1
fi

DATA_DIR="$1"
shift

# locate repository and manifest
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MANIFEST="${REPO_DIR}/manifests/GWAS_manifest.tsv"

# MAGMA and reference files
MAGMA="${DATA_DIR}/magma/magma"
GENE_LOC="${DATA_DIR}/reference/magma/NCBI37.3/NCBI37.3.gene.loc"

# check required files exist
if [[ ! -f "${MAGMA}" ]]; then
  echo "Error: MAGMA executable not found:"
  echo "${MAGMA}"
  exit 1
fi

if [[ ! -f "${GENE_LOC}" ]]; then
  echo "Error: gene-location file not found:"
  echo "${GENE_LOC}"
  exit 1
fi

echo "Manifest: ${MANIFEST}"
echo "Data directory: ${DATA_DIR}"
echo "MAGMA: ${MAGMA}"
echo "Gene locations: ${GENE_LOC}"

# process each selected trait
for TRAIT in "$@"; do

  echo
  echo "========================================"
  echo "Annotating ${TRAIT} SNPs to genes"
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
  SNP_LOC="${ANALYSIS_DIR}/${ANALYSIS_ID}.magma.snploc"
  OUTPUT_PREFIX="${ANALYSIS_DIR}/${ANALYSIS_ID}"

  # check SNP-location file exists
  if [[ ! -f "${SNP_LOC}" ]]; then
    echo "Error: SNP-location file not found:"
    echo "${SNP_LOC}"
    exit 1
  fi

  echo "Input: ${SNP_LOC}"
  echo "Output prefix: ${OUTPUT_PREFIX}"

  "${MAGMA}" \
    --annotate \
    --snp-loc "${SNP_LOC}" \
    --gene-loc "${GENE_LOC}" \
    --out "${OUTPUT_PREFIX}"

  echo "Annotation completed for ${TRAIT}."

done

echo
echo "Finished annotating selected GWAS."
