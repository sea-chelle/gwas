#!/bin/bash

# run MAGMA gene-level association analysis
#
# usage:
#   bash gene_analysis.sh <data_directory> <trait> [trait ...]
#
# examples:
#   bash gene_analysis.sh ~/Scratch/gwas ASD
#   bash gene_analysis.sh ~/Scratch/gwas ASD ADHD ANX MDD PTSD SCZ BD

set -euo pipefail

# require a data directory and at least one trait
if [[ "$#" -lt 2 ]]; then
  echo "Usage: bash gene_analysis.sh <data_directory> <trait> [trait ...]"
  exit 1
fi

DATA_DIR="$1"
shift

# locate repository and manifest
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MANIFEST="${REPO_DIR}/manifests/GWAS_manifest.tsv"

# MAGMA executable and LD reference
MAGMA="${DATA_DIR}/magma/magma"
LD_REF="${DATA_DIR}/reference/magma/g1000_eur/g1000_eur"

# check required files exist
if [[ ! -f "${MAGMA}" ]]; then
  echo "Error: MAGMA executable not found:"
  echo "${MAGMA}"
  exit 1
fi

if [[ ! -f "${LD_REF}.bed" ||
      ! -f "${LD_REF}.bim" ||
      ! -f "${LD_REF}.fam" ]]; then
  echo "Error: LD reference files not found for prefix:"
  echo "${LD_REF}"
  exit 1
fi

echo "Manifest: ${MANIFEST}"
echo "Data directory: ${DATA_DIR}"
echo "MAGMA: ${MAGMA}"
echo "LD reference: ${LD_REF}"

# process each selected trait
for TRAIT in "$@"; do

  echo
  echo "========================================"
  echo "Running MAGMA gene analysis for ${TRAIT}"
  echo "========================================"

  # get manifest row for selected trait
  TRAIT_ROW=$(awk -F'\t' -v trait="${TRAIT}" '
    NR > 1 && $1 == trait {print; exit}
  ' "${MANIFEST}")

  if [[ -z "${TRAIT_ROW}" ]]; then
    echo "Error: trait '${TRAIT}' not found in manifest."
    exit 1
  fi

  # read manifest fields
  IFS=$'\t' read -r \
    TRAIT_NAME \
    STUDY \
    ANALYSIS_ID \
    FILE \
    ANCESTRY \
    GENOME_BUILD \
    EFFECT \
    SNP_COL \
    CHR_COL \
    BP_COL \
    P_COL \
    N_METHOD \
    N_SOURCE \
    N_TRANSFORM \
    SAMPLE_SIZE \
    NOTES <<< "${TRAIT_ROW}"

  # define analysis files
  ANALYSIS_DIR="${DATA_DIR}/$(dirname "${FILE}")"
  PVAL_FILE="${ANALYSIS_DIR}/${ANALYSIS_ID}.magma.pval"
  GENE_ANNOT="${ANALYSIS_DIR}/${ANALYSIS_ID}.genes.annot"
  OUTPUT_PREFIX="${ANALYSIS_DIR}/${ANALYSIS_ID}"

  # check required input files exist
  if [[ ! -f "${PVAL_FILE}" ]]; then
    echo "Error: P-value file not found:"
    echo "${PVAL_FILE}"
    exit 1
  fi

  if [[ ! -f "${GENE_ANNOT}" ]]; then
    echo "Error: gene annotation file not found:"
    echo "${GENE_ANNOT}"
    exit 1
  fi

  echo
  echo "Study: ${STUDY}"
  echo "Analysis ID: ${ANALYSIS_ID}"
  echo "P-value file: ${PVAL_FILE}"
  echo "Gene annotation: ${GENE_ANNOT}"
  echo "N method: ${N_METHOD}"

  # run MAGMA using fixed or per-SNP sample size
  if [[ "${N_METHOD}" == "fixed" ]]; then

    if [[ -z "${SAMPLE_SIZE}" || "${SAMPLE_SIZE}" == "NA" ]]; then
      echo "Error: fixed sample size not defined for ${TRAIT}."
      exit 1
    fi

    echo "Sample size: ${SAMPLE_SIZE}"

    "${MAGMA}" \
      --bfile "${LD_REF}" \
      --pval "${PVAL_FILE}" N="${SAMPLE_SIZE}" \
      --gene-annot "${GENE_ANNOT}" \
      --out "${OUTPUT_PREFIX}"

  elif [[ "${N_METHOD}" == "per_snp" ]]; then

    echo "Sample size: per-SNP N column"

    "${MAGMA}" \
      --bfile "${LD_REF}" \
      --pval "${PVAL_FILE}" ncol=N \
      --gene-annot "${GENE_ANNOT}" \
      --out "${OUTPUT_PREFIX}"

  else

    echo "Error: unsupported N method '${N_METHOD}' for ${TRAIT}." >&2
    exit 1

  fi

  echo
  echo "Gene analysis completed for ${TRAIT}."

done

echo
echo "Finished MAGMA gene analysis for selected GWAS."
