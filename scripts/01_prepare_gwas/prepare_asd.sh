#!/bin/bash

# prepare Grove et al. (2019) ASD GWAS summary statistics for MAGMA.
#
# usage:
#   bash prepare_asd.sh <input_gwas.gz> <output_directory>

set -euo pipefail

INPUT_GWAS="$1"
OUTPUT_DIR="$2"

mkdir -p "${OUTPUT_DIR}"

# extract SNP identifier and association P-value for MAGMA gene analysis
zcat "${INPUT_GWAS}" |
awk 'BEGIN {OFS="\t"} NR==1 {print "SNP","P"} NR>1 {print $2,$9}' \
> "${OUTPUT_DIR}/ASD_Grove_2019.magma.pval"

# extract SNP genomic coordinates for MAGMA SNP-to-gene annotation
zcat "${INPUT_GWAS}" |
awk 'BEGIN {OFS="\t"} NR==1 {print "SNP","CHR","BP"} NR>1 {print $2,$1,$3}' \
> "${OUTPUT_DIR}/ASD_Grove_2019.magma.snploc"
