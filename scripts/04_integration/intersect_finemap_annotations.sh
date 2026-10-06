#!/bin/bash
#$ -cwd
#$ -V
#$ -l h_rt=01:00:00
#$ -l h_vmem=4G
#$ -pe smp 1
#$ -N finemap_intersect
#$ -o logs/
#$ -e logs/

set -euo pipefail

# usage:
# qsub scripts/04_integration/intersect_finemap_annotations.sh \
#   ~/Scratch/gwas

DATA_DIR=${1:?Please provide DATA_DIR}

FINEMAP_DIR="${DATA_DIR}/locus_integration/finemap"
ANNOT_DIR="${DATA_DIR}/ldsc_annotations"
OUT_DIR="${DATA_DIR}/locus_integration/overlaps"

mkdir -p "${OUT_DIR}"

traits=(ASD ADHD ANX MDD PTSD SCZ BD)

annotations=(RNA_Overall_60 RNA_Male_60 RNA_Female_60 ATAC_Overall ATAC_Male ATAC_Female)

declare -A annotation_files=(
  [RNA_Overall_60]="RNA_Overall_60_GRCh37_100kb.bed"
  [RNA_Male_60]="RNA_Male_60_GRCh37_100kb.bed"
  [RNA_Female_60]="RNA_Female_60_GRCh37_100kb.bed"
  [ATAC_Overall]="ATAC_Overall_hg19.bed"
  [ATAC_Male]="ATAC_Male_hg19.bed"
  [ATAC_Female]="ATAC_Female_hg19.bed"
)

summary="${OUT_DIR}/finemap_overlap_summary.tsv"

printf "Trait\tAnnotation\tFineMappedVariants\tOverlapRecords\tUniqueOverlappingSNPs\n" \
  > "${summary}"

for trait in "${traits[@]}"; do

  finemap="${FINEMAP_DIR}/${trait}_finemap_GRCh37.bed"

  if [[ ! -f "${finemap}" ]]; then
    echo "ERROR: fine-mapping BED not found: ${finemap}"
    exit 1
  fi

  n_finemap=$(wc -l < "${finemap}")

  echo ""
  echo "${trait}: ${n_finemap} fine-mapped variants"

  for annotation in "${annotations[@]}"; do

    annot="${ANNOT_DIR}/${annotation_files[$annotation]}"
    outfile="${OUT_DIR}/${trait}_${annotation}_overlap.bed"

    if [[ ! -f "${annot}" ]]; then
      echo "ERROR: annotation BED not found: ${annot}"
      exit 1
    fi

    bedtools intersect \
      -a "${finemap}" \
      -b "${annot}" \
      -wa \
      -wb \
      > "${outfile}"

    n_records=$(wc -l < "${outfile}")

    if [[ "${n_records}" -gt 0 ]]; then
      n_snps=$(cut -f4 "${outfile}" | sort -u | wc -l)
    else
      n_snps=0
    fi

    printf "%s\t%s\t%s\t%s\t%s\n" \
      "${trait}" \
      "${annotation}" \
      "${n_finemap}" \
      "${n_records}" \
      "${n_snps}" \
      >> "${summary}"

    echo "  ${annotation}: ${n_snps} overlapping SNPs"

  done
done

echo ""
echo "Finished."
echo "Overlap files: ${OUT_DIR}"
echo "Summary: ${summary}"
