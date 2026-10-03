#!/bin/bash

# perform MAGMA SNP-to-gene annotation
#
# usage:
#   bash annotate_genes.sh <magma_executable> <snp_location_file> <gene_location_file> <output_prefix>

set -euo pipefail

MAGMA="$1"
SNP_LOC="$2"
GENE_LOC="$3"
OUTPUT_PREFIX="$4"

"${MAGMA}" \
  --annotate \
  --snp-loc "${SNP_LOC}" \
  --gene-loc "${GENE_LOC}" \
  --out "${OUTPUT_PREFIX}"
