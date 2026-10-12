
#!/bin/bash

# lift mouse ATAC-seq annotations from mm39 to human hg19/GRCh37
#
# includes stress-responsive DARs and consensus ATAC-seq peaks.
# automatically identifies all available ATAC-seq annotations.
#
# a minimum match ratio of 0.1 is used for cross-species liftOver,
# corresponding to the UCSC mouse-to-human liftOver setting.

set -euo pipefail

# input directory and chain file
input_dir="${ANNOTATION_DIR:-ldsc_annotations}"
chain_file="${CHAIN_FILE:-reference/ldsc/liftover/mm39ToHg19.over.chain.gz}"

# check required executable and chain file
if ! command -v liftOver >/dev/null 2>&1; then
  echo "ERROR: UCSC liftOver executable not found in PATH"
  exit 1
fi

if [ ! -f "${chain_file}" ]; then
  echo "ERROR: chain file not found: ${chain_file}"
  exit 1
fi

# identify available mouse ATAC-seq annotations
shopt -s nullglob
atac_files=("${input_dir}"/ATAC_*_mm39.bed)

if [ ${#atac_files[@]} -eq 0 ]; then
  echo "ERROR: No ATAC-seq BED files found in ${input_dir}"
  exit 1
fi

# report number of annotations
echo "Found ${#atac_files[@]} ATAC-seq annotations"
echo

# lift ATAC-seq annotations to hg19
for input_file in "${atac_files[@]}"
do
  annotation=$(basename "${input_file}" _mm39.bed)

  echo "=== ${annotation} ==="

  liftOver -minMatch=0.1 \
    "${input_file}" \
    "${chain_file}" \
    "${input_dir}/${annotation}_hg19.bed" \
    "${input_dir}/${annotation}_unmapped.bed"

  # report mapping success
  echo "Input:"
  wc -l < "${input_file}"

  echo "Mapped:"
  wc -l < "${input_dir}/${annotation}_hg19.bed"

  echo "Unique mapped peaks:"
  cut -f4 "${input_dir}/${annotation}_hg19.bed" | sort -u | wc -l

  echo "Unmapped:"
  grep -vc '^#' "${input_dir}/${annotation}_unmapped.bed" || true

  echo
done

echo "liftOver complete"

