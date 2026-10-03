#!/bin/bash

# run MAGMA gene-level association analysis
#
# usage:
#   bash gene_analysis.sh <magma_executable> <ld_reference_prefix> \
#     <pvalue_file> <sample_size> <gene_annotation> <output_prefix>

set -euo pipefail

MAGMA="$1"
LD_REF="$2"
PVAL_FILE="$3"
SAMPLE_SIZE="$4"
GENE_ANNOT="$5"
OUTPUT_PREFIX="$6"

"${MAGMA}" \
  --bfile "${LD_REF}" \
  --pval "${PVAL_FILE}" N="${SAMPLE_SIZE}" \
  --gene-annot "${GENE_ANNOT}" \
  --out "${OUTPUT_PREFIX}"
