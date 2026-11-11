#!/bin/bash

# prepare GWAS summary statistics for S-LDSC
#
# usage:
#   bash prepare_gwas_ldsc.sh <data_directory> <trait> [trait ...]
#
# examples:
#   bash prepare_gwas_ldsc.sh ~/Scratch/gwas ADHD
#   bash prepare_gwas_ldsc.sh ~/Scratch/gwas ASD ADHD ANX

set -euo pipefail

# require a data directory and at least one trait
if [[ "$#" -lt 2 ]]; then
  echo "Usage: bash prepare_gwas_ldsc.sh <data_directory> <trait> [trait ...]"
  exit 1
fi

DATA_DIR="$1"
shift

# locate repository and manifest
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MANIFEST="${REPO_DIR}/manifests/GWAS_manifest.tsv"

# locate LDSC
LDSC_DIR="${DATA_DIR}/ldsc"
MUNGE_SUMSTATS="${LDSC_DIR}/munge_sumstats.py"

echo "Manifest: ${MANIFEST}"
echo "Data directory: ${DATA_DIR}"
echo "LDSC directory: ${LDSC_DIR}"

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
  echo "Preparing ${TRAIT} GWAS for S-LDSC"
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
  echo "Effect: ${EFFECT}"
  echo "SNP column: ${SNP_COL}"
  echo "A1 column: ${A1_COL}"
  echo "A2 column: ${A2_COL}"
  echo "P-value column: ${P_COL}"
  echo "INFO column: ${INFO_COL}"
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
  A1_IDX=$(get_col_index "${A1_COL}")
  A2_IDX=$(get_col_index "${A2_COL}")
  EFFECT_IDX=$(get_col_index "${EFFECT}")
  P_IDX=$(get_col_index "${P_COL}")

  echo
  echo "Detected column numbers:"
  echo "SNP: ${SNP_IDX}"
  echo "A1: ${A1_IDX}"
  echo "A2: ${A2_IDX}"
  echo "Effect: ${EFFECT_IDX}"
  echo "P: ${P_IDX}"

  # detect INFO column where available
  if [[ "${INFO_COL}" != "NA" ]]; then
    INFO_IDX=$(get_col_index "${INFO_COL}")
    echo "INFO: ${INFO_IDX}"
  fi

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
  for idx in "${SNP_IDX}" "${A1_IDX}" "${A2_IDX}" "${EFFECT_IDX}" "${P_IDX}"; do
    if [[ -z "${idx}" ]]; then
      echo "Error: one or more required GWAS columns were not found."
      exit 1
    fi
  done

  # check INFO column was found where required
  if [[ "${INFO_COL}" != "NA" && -z "${INFO_IDX}" ]]; then
    echo "Error: INFO column '${INFO_COL}' was not found."
    exit 1
  fi

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

  LDSC_INPUT="${OUTPUT_PREFIX}.ldsc.input"

  echo
  echo "Output LDSC input file: ${LDSC_INPUT}"
  echo "Output munged summary statistics: ${OUTPUT_PREFIX}.sumstats.gz"

  # create standardised LDSC input file
  echo
  echo "Creating LDSC input file..."

  if [[ "${N_METHOD}" == "fixed" ]]; then

    zcat "${INPUT_FILE}" | awk \
      -v snp_idx="${SNP_IDX}" \
      -v a1_idx="${A1_IDX}" \
      -v a2_idx="${A2_IDX}" \
      -v effect_idx="${EFFECT_IDX}" \
      -v p_idx="${P_IDX}" \
      -v info_idx="${INFO_IDX}" \
      -v n="${SAMPLE_SIZE}" \
      '
      BEGIN {
        OFS="\t"
        print "SNP", "A1", "A2", "EFFECT", "P", "INFO", "N"
      }

      /^##/ {next}

      !header_seen {
        header_seen=1
        next
      }

      {
        print $snp_idx, $a1_idx, $a2_idx, $effect_idx, $p_idx, $info_idx, n
      }
      ' > "${LDSC_INPUT}"

  elif [[ "${N_TRANSFORM}" == "none" ]]; then

    if [[ "${INFO_COL}" == "NA" ]]; then

      zcat "${INPUT_FILE}" | awk \
        -v snp_idx="${SNP_IDX}" \
        -v a1_idx="${A1_IDX}" \
        -v a2_idx="${A2_IDX}" \
        -v effect_idx="${EFFECT_IDX}" \
        -v p_idx="${P_IDX}" \
        -v n_idx="${N_IDX}" \
        '
        BEGIN {
          OFS="\t"
          print "SNP", "A1", "A2", "EFFECT", "P", "N"
        }

        /^##/ {next}

        !header_seen {
          header_seen=1
          next
        }

        {
          print $snp_idx, $a1_idx, $a2_idx, $effect_idx, $p_idx, $n_idx
        }
        ' > "${LDSC_INPUT}"

    else

      zcat "${INPUT_FILE}" | awk \
        -v snp_idx="${SNP_IDX}" \
        -v a1_idx="${A1_IDX}" \
        -v a2_idx="${A2_IDX}" \
        -v effect_idx="${EFFECT_IDX}" \
        -v p_idx="${P_IDX}" \
        -v info_idx="${INFO_IDX}" \
        -v n_idx="${N_IDX}" \
        '
        BEGIN {
          OFS="\t"
          print "SNP", "A1", "A2", "EFFECT", "P", "INFO", "N"
        }

        /^##/ {next}

        !header_seen {
          header_seen=1
          next
        }

        {
          print $snp_idx, $a1_idx, $a2_idx, $effect_idx, $p_idx, $info_idx, $n_idx
        }
        ' > "${LDSC_INPUT}"

    fi

  elif [[ "${N_TRANSFORM}" == "multiply_by_2" ]]; then

    zcat "${INPUT_FILE}" | awk \
      -v snp_idx="${SNP_IDX}" \
      -v a1_idx="${A1_IDX}" \
      -v a2_idx="${A2_IDX}" \
      -v effect_idx="${EFFECT_IDX}" \
      -v p_idx="${P_IDX}" \
      -v info_idx="${INFO_IDX}" \
      -v n_idx="${N_IDX}" \
      '
      BEGIN {
        OFS="\t"
        print "SNP", "A1", "A2", "EFFECT", "P", "INFO", "N"
      }

      /^##/ {next}

      !header_seen {
        header_seen=1
        next
      }

      {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%.2f\n",
          $snp_idx, $a1_idx, $a2_idx, $effect_idx, $p_idx, $info_idx, 2 * $n_idx
      }
      ' > "${LDSC_INPUT}"

  elif [[ "${N_TRANSFORM}" == "case_control_effective_N" ]]; then

    zcat "${INPUT_FILE}" | awk \
      -v snp_idx="${SNP_IDX}" \
      -v a1_idx="${A1_IDX}" \
      -v a2_idx="${A2_IDX}" \
      -v effect_idx="${EFFECT_IDX}" \
      -v p_idx="${P_IDX}" \
      -v info_idx="${INFO_IDX}" \
      -v n_case_idx="${N_CASE_IDX}" \
      -v n_control_idx="${N_CONTROL_IDX}" \
      '
      BEGIN {
        OFS="\t"
        print "SNP", "A1", "A2", "EFFECT", "P", "INFO", "N"
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

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%.2f\n",
          $snp_idx, $a1_idx, $a2_idx, $effect_idx, $p_idx, $info_idx, n_eff
      }
      ' > "${LDSC_INPUT}"

  else

    echo "Error: unsupported N transformation '${N_TRANSFORM}'." >&2
    exit 1

  fi

  echo "LDSC input file created."

  # define signed effect statistic
  if [[ "${EFFECT}" == "OR" ]]; then
    SIGNED_SUMSTATS="EFFECT,1"
  elif [[ "${EFFECT}" == "BETA" ]]; then
    SIGNED_SUMSTATS="EFFECT,0"
  else
    echo "Error: unsupported effect type '${EFFECT}'." >&2
    exit 1
  fi

  # munge summary statistics
  echo
  echo "Munging summary statistics..."

  if [[ "${INFO_COL}" == "NA" ]]; then

    python "${MUNGE_SUMSTATS}" \
      --sumstats "${LDSC_INPUT}" \
      --snp SNP \
      --a1 A1 \
      --a2 A2 \
      --p P \
      --N-col N \
      --signed-sumstats "${SIGNED_SUMSTATS}" \
      --out "${OUTPUT_PREFIX}"

  else

    python "${MUNGE_SUMSTATS}" \
      --sumstats "${LDSC_INPUT}" \
      --snp SNP \
      --a1 A1 \
      --a2 A2 \
      --p P \
      --info INFO \
      --N-col N \
      --signed-sumstats "${SIGNED_SUMSTATS}" \
      --out "${OUTPUT_PREFIX}"

  fi

  echo "Summary statistics munged."

  # check output files exist
  if [[ ! -f "${OUTPUT_PREFIX}.sumstats.gz" ]]; then
    echo "Error: munged summary statistics were not created." >&2
    exit 1
  fi

  if [[ ! -f "${OUTPUT_PREFIX}.log" ]]; then
    echo "Error: LDSC munging log was not created." >&2
    exit 1
  fi

  echo "QC passed for ${TRAIT}."

done

echo
echo "Finished preparing selected GWAS."
