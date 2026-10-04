#!/bin/bash

# lift mouse ATAC-seq DAR coordinates from mm39 to human hg19/GRCh37
#
# a minimum match ratio of 0.1 is used for cross-species liftOver,
# corresponding to the UCSC mouse-to-human liftOver setting.

# input directory and chain file
input_dir="ldsc_annotations"
chain_file="reference/ldsc/liftover/mm39ToHg19.over.chain.gz"

# lift ATAC-seq DARs to hg19
for group in Overall Male Female
do
  liftOver -minMatch=0.1 \
    "${input_dir}/ATAC_${group}_mm39.bed" \
    "${chain_file}" \
    "${input_dir}/ATAC_${group}_hg19.bed" \
    "${input_dir}/ATAC_${group}_unmapped.bed"
done

# report mapping success
for group in Overall Male Female
do
  echo "=== ${group} ==="
  echo "Input:"
  wc -l < "${input_dir}/ATAC_${group}_mm39.bed"
  echo "Mapped:"
  wc -l < "${input_dir}/ATAC_${group}_hg19.bed"
  echo "Unique mapped peaks:"
  cut -f4 "${input_dir}/ATAC_${group}_hg19.bed" | sort -u | wc -l
  echo "Unmapped:"
  grep -v '^#' "${input_dir}/ATAC_${group}_unmapped.bed" | wc -l
done
