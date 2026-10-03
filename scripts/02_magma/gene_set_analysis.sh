#!/bin/bash

# run MAGMA competitive gene-set analysis
#
# usage:
#   bash gene_set_analysis.sh <magma_executable> <gene_results> \
#     <gene_set_annotation> <output_prefix>

set -euo pipefail

MAGMA="$1"
GENE_RESULTS="$2"
GENE_SET_ANNOT="$3"
OUTPUT_PREFIX="$4"

"${MAGMA}" \
  --gene-results "${GENE_RESULTS}" \
  --set-annot "${GENE_SET_ANNOT}" \
  --out "${OUTPUT_PREFIX}"
