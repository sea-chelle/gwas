#!/bin/bash

# prepare GWAS summary statistics for MAGMA
#
# usage:
#   bash prepare_gwas.sh <data_directory> <trait> [trait ...]
#
# examples:
#   bash prepare_gwas.sh ~/Scratch/gwas ADHD
#   bash prepare_gwas.sh ~/Scratch/gwas ASD ADHD ANX

set -euo pipefail

# require a data directory and at least one trait
if [[ "$#" -lt 2 ]]; then
  echo "Usage: bash prepare_gwas.sh <data_directory> <trait> [trait ...]"
  exit 1
fi

DATA_DIR="$1"
shift

# locate repository and manifest
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MANIFEST="${REPO_DIR}/manifests/GWAS_manifest.tsv"

echo "Manifest: ${MANIFEST}"
echo "Data directory: ${DATA_DIR}"

# function to find the numerical index of a column from the GWAS header
get_col_index() {
  local column_name="$1"

  awk -v target="${column_name}" '
    BEGIN {FS="[[:space:]]+"}
    {
      for (i = 1; i <= NF; i++) {
        if ($i == target) {
          print i
          exit
        }
      }
    }
  ' <<< "${HEADER}"
}

# process each selected trait
for TRAIT in "$@"; do

  echo
  echo "========================================"
  echo "Preparing ${TRAIT} GWAS for MAGMA"
  echo "========================================"

  # get manifest row for selected trait
  TRAIT_ROW=$(awk -F'\t' -v trait="${TRAIT}" '
    NR > 1 && $1 == trait {print; exit}
  ' "${MANIFEST}")

  # stop if trait is not found
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
  A1_COL \
  A2_COL \
  INFO_COL \
  P_COL \
  N_METHOD \
  N_SOURCE \
  N_TRANSFORM \
  SAMPLE_SIZE \
  NOTES <<< "${TRAIT_ROW}"

  INPUT_FILE="${DATA_DIR}/${FILE}"

  # check input file exists
  if [[ ! -f "${INPUT_FILE}" ]]; then
    echo "Error: input file not found:"
    echo "${INPUT_FILE}"
    exit 1
  fi

  # show selected configuration
  echo
  echo "Study: ${STUDY}"
  echo "Analysis ID: ${ANALYSIS_ID}"
  echo "Input: ${INPUT_FILE}"
  echo "SNP column: ${SNP_COL}"
  echo "Chromosome column: ${CHR_COL}"
  echo "Position column: ${BP_COL}"
  echo "P-value column: ${P_COL}"
  echo "N method: ${N_METHOD}"
  echo "N source: ${N_SOURCE}"
  echo "N transform: ${N_TRANSFORM}"
  echo "Sample size: ${SAMPLE_SIZE}"

  # get first non-metadata line as the GWAS header
  HEADER=$(zgrep -m 1 -v '^##' "${INPUT_FILE}")

  echo
  echo "GWAS header:"
  echo "${HEADER}"

  # detect required column numbers
  SNP_IDX=$(get_col_index "${SNP_COL}")
  CHR_IDX=$(get_col_index "${CHR_COL}")
  BP_IDX=$(get_col_index "${BP_COL}")
  P_IDX=$(get_col_index "${P_COL}")

  echo
  echo "Detected column numbers:"
  echo "SNP: ${SNP_IDX}"
  echo "CHR: ${CHR_IDX}"
  echo "BP: ${BP_IDX}"
  echo "P: ${P_IDX}"

  # detect sample-size column(s) where required
  if [[ "${N_METHOD}" == "per_snp" ]]; then

    if [[ "${N_TRANSFORM}" == "case_control_effective_N" ]]; then

      IFS=',' read -r N_CASE_COL N_CONTROL_COL <<< "${N_SOURCE}"

      N_CASE_IDX=$(get_col_index "${N_CASE_COL}")
      N_CONTROL_IDX=$(get_col_index "${N_CONTROL_COL}")

      echo "N case: ${N_CASE_IDX}"
      echo "N control: ${N_CONTROL_IDX}"

    else

      N_IDX=$(get_col_index "${N_SOURCE}")
      echo "N: ${N_IDX}"

    fi
  fi

  # check required GWAS columns were found
  for idx in "${SNP_IDX}" "${CHR_IDX}" "${BP_IDX}" "${P_IDX}"; do
    if [[ -z "${idx}" ]]; then
      echo "Error: one or more required GWAS columns were not found."
      exit 1
    fi
  done

  # check required sample-size columns were found
  if [[ "${N_METHOD}" == "per_snp" ]]; then

    if [[ "${N_TRANSFORM}" == "case_control_effective_N" ]]; then

      if [[ -z "${N_CASE_IDX}" || -z "${N_CONTROL_IDX}" ]]; then
        echo "Error: case/control sample-size columns were not found."
        exit 1
      fi

    else

      if [[ -z "${N_IDX}" ]]; then
        echo "Error: sample-size column '${N_SOURCE}' was not found."
        exit 1
      fi

    fi
  fi

  # define output files
  OUTPUT_DIR="$(dirname "${INPUT_FILE}")"
  OUTPUT_PREFIX="${OUTPUT_DIR}/${ANALYSIS_ID}"

  SNPLOC_OUT="${OUTPUT_PREFIX}.magma.snploc"
  PVAL_OUT="${OUTPUT_PREFIX}.magma.pval"

  echo
  echo "Output SNP location file: ${SNPLOC_OUT}"
  echo "Output P-value file: ${PVAL_OUT}"

  # create MAGMA SNP location file
  echo
  echo "Creating SNP location file..."

  zcat "${INPUT_FILE}" | awk \
    -v snp_idx="${SNP_IDX}" \
    -v chr_idx="${CHR_IDX}" \
    -v bp_idx="${BP_IDX}" \
    '
    BEGIN {
      OFS="\t"
      print "SNP", "CHR", "BP"
    }

    /^##/ {next}

    !header_seen {
      header_seen=1
      next
    }

    {
      print $snp_idx, $chr_idx, $bp_idx
    }
    ' > "${SNPLOC_OUT}"

  echo "SNP location file created."

  # create MAGMA P-value file
  echo
  echo "Creating P-value file..."

  if [[ "${N_METHOD}" == "fixed" ]]; then

    zcat "${INPUT_FILE}" | awk \
      -v snp_idx="${SNP_IDX}" \
      -v p_idx="${P_IDX}" \
      '
      BEGIN {
        OFS="\t"
        print "SNP", "P"
      }

      /^##/ {next}

      !header_seen {
        header_seen=1
        next
      }

      {
        print $snp_idx, $p_idx
      }
      ' > "${PVAL_OUT}"

  elif [[ "${N_TRANSFORM}" == "none" ]]; then

    zcat "${INPUT_FILE}" | awk \
      -v snp_idx="${SNP_IDX}" \
      -v p_idx="${P_IDX}" \
      -v n_idx="${N_IDX}" \
      '
      BEGIN {
        OFS="\t"
        print "SNP", "P", "N"
      }

      /^##/ {next}

      !header_seen {
        header_seen=1
        next
      }

      {
        print $snp_idx, $p_idx, $n_idx
      }
      ' > "${PVAL_OUT}"

  elif [[ "${N_TRANSFORM}" == "multiply_by_2" ]]; then

    zcat "${INPUT_FILE}" | awk \
      -v snp_idx="${SNP_IDX}" \
      -v p_idx="${P_IDX}" \
      -v n_idx="${N_IDX}" \
      '
      BEGIN {
        OFS="\t"
        print "SNP", "P", "N"
      }

      /^##/ {next}

      !header_seen {
        header_seen=1
        next
      }

      {
        printf "%s\t%s\t%.2f\n", $snp_idx, $p_idx, 2 * $n_idx
      }
      ' > "${PVAL_OUT}"

  elif [[ "${N_TRANSFORM}" == "case_control_effective_N" ]]; then

    zcat "${INPUT_FILE}" | awk \
      -v snp_idx="${SNP_IDX}" \
      -v p_idx="${P_IDX}" \
      -v n_case_idx="${N_CASE_IDX}" \
      -v n_control_idx="${N_CONTROL_IDX}" \
      '
      BEGIN {
        OFS="\t"
        print "SNP", "P", "N"
      }

      /^##/ {next}

      !header_seen {
        header_seen=1
        next
      }

      {
        n_case = $n_case_idx
        n_control = $n_control_idx
        n_eff = 4 * n_case * n_control / (n_case + n_control)

        printf "%s\t%s\t%.2f\n", $snp_idx, $p_idx, n_eff
      }
      ' > "${PVAL_OUT}"

  else

    echo "Error: unsupported N transformation '${N_TRANSFORM}'." >&2
    exit 1

  fi

  echo "P-value file created."

  # QC generated MAGMA input files
  echo
  echo "Running QC..."

  # count variants, excluding header
  PVAL_N=$(awk 'END {print NR - 1}' "${PVAL_OUT}")
  SNPLOC_N=$(awk 'END {print NR - 1}' "${SNPLOC_OUT}")

  # check P-value file
  if [[ "${N_METHOD}" == "fixed" ]]; then

    read -r INVALID_PVAL_FIELDS INVALID_SNP INVALID_P < <(
      awk '
        NR == 1 {next}

        NF != 2 {bad_fields++}
        $1 == "" {bad_snp++}
        $2 !~ /^[0-9.eE+-]+$/ || $2 <= 0 || $2 > 1 {bad_p++}

        END {
          print bad_fields+0, bad_snp+0, bad_p+0
        }
      ' "${PVAL_OUT}"
    )

    INVALID_N="NA"

  else

    read -r INVALID_PVAL_FIELDS INVALID_SNP INVALID_P INVALID_N < <(
      awk '
        NR == 1 {next}

        NF != 3 {bad_fields++}
        $1 == "" {bad_snp++}
        $2 !~ /^[0-9.eE+-]+$/ || $2 <= 0 || $2 > 1 {bad_p++}
        $3 !~ /^[0-9.eE+-]+$/ || $3 <= 0 {bad_n++}

        END {
          print bad_fields+0, bad_snp+0, bad_p+0, bad_n+0
        }
      ' "${PVAL_OUT}"
    )

  fi

  # check SNP location file
  read -r INVALID_SNPLOC_FIELDS INVALID_SNPLOC_SNP INVALID_CHR INVALID_BP < <(
    awk '
      NR == 1 {next}

      NF != 3 {bad_fields++}
      $1 == "" {bad_snp++}
      $2 !~ /^([0-9]+|X|Y|XY|MT|M)$/ {bad_chr++}
      $3 !~ /^[0-9]+$/ || $3 <= 0 {bad_bp++}

      END {
        print bad_fields+0, bad_snp+0, bad_chr+0, bad_bp+0
      }
    ' "${SNPLOC_OUT}"
  )

  # report QC
  echo "  P-value variants: ${PVAL_N}"
  echo "  SNP-location variants: ${SNPLOC_N}"
  echo "  Invalid P-value rows: ${INVALID_PVAL_FIELDS}"
  echo "  Missing SNP IDs in P-value file: ${INVALID_SNP}"
  echo "  Invalid P-values: ${INVALID_P}"

  if [[ "${N_METHOD}" == "per_snp" ]]; then
    echo "  Invalid sample sizes: ${INVALID_N}"
  fi

  echo "  Invalid SNP-location rows: ${INVALID_SNPLOC_FIELDS}"
  echo "  Missing SNP IDs in SNP-location file: ${INVALID_SNPLOC_SNP}"
  echo "  Invalid chromosomes: ${INVALID_CHR}"
  echo "  Invalid positions: ${INVALID_BP}"

  # fail QC if problems are detected
  QC_FAILED=0

  if [[ "${PVAL_N}" -ne "${SNPLOC_N}" ]]; then
    echo "  Row counts match: NO"
    QC_FAILED=1
  else
    echo "  Row counts match: yes"
  fi

  if [[ "${INVALID_PVAL_FIELDS}" -gt 0 ||
        "${INVALID_SNP}" -gt 0 ||
        "${INVALID_P}" -gt 0 ||
        "${INVALID_SNPLOC_FIELDS}" -gt 0 ||
        "${INVALID_SNPLOC_SNP}" -gt 0 ||
        "${INVALID_CHR}" -gt 0 ||
        "${INVALID_BP}" -gt 0 ]]; then
    QC_FAILED=1
  fi

  if [[ "${N_METHOD}" == "per_snp" && "${INVALID_N}" -gt 0 ]]; then
    QC_FAILED=1
  fi

  if [[ "${QC_FAILED}" -eq 1 ]]; then
    echo "QC FAILED for ${TRAIT}." >&2
    exit 1
  fi

  echo "QC passed for ${TRAIT}."

done

echo
echo "Finished preparing selected GWAS."
