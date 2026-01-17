#!/bin/bash
set -e

# ======================
# Project root
# ======================

PROJECT_ROOT="$(cd "$(dirname "$0")"/.. && pwd)"

# ======================
# Required inputs
# ======================

IDS_REL=""        # Relative path UNDER data/
TASKS=()          # Tasks to run (comments | features | ec | cross_references)

# ======================
# Optional inputs
# ======================

OUT_SUBDIR=""     # Relative subdirectory UNDER data/Sample/
BATCH_SIZE=100    # Number of IDs per batch (UniProt REST API max: 100 IDs/request)
OUTPUT_FORMAT="jsonl"   # Output format (jsonl | csv | tsv)
TEMPLATE_SUBDIR="Template"   # Relative template directory UNDER data/ (default: Template)

# ======================
# CLI parse
# ======================

usage() {
  echo "Usage: bash $0 --tasks <one_or_more_tasks...> --ids_rel <file_under_data> [--out_subdir <subdir_under_data_Sample>] [--template_subdir <dir_under_data>] [--batch_size <int>] [--output_format <jsonl|csv|tsv>]"
  echo
  echo "Tasks: comments | features | ec | cross_references"
  echo
  echo "Examples:"
  echo "  bash $0 --tasks comments --ids_rel demo_ids.txt"
  echo "  bash $0 --tasks comments cross_references ec features --ids_rel train_ids.txt --out_subdir train_samples_pt --template_subdir Template_pt --batch_size 100 --output_format tsv"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tasks)
      shift
      # Collect task names until the next flag or end of arguments
      while [[ $# -gt 0 && "$1" != --* ]]; do
        TASKS+=("$1")
        shift
      done
      ;;
    --ids_rel)
      IDS_REL="${2:-}"
      shift 2
      ;;
    --out_subdir)
      OUT_SUBDIR="${2:-}"
      shift 2
      ;;
    --template_subdir)
      TEMPLATE_SUBDIR="${2:-}"
      shift 2
      ;;
    --batch_size)
      BATCH_SIZE="${2:-}"
      shift 2
      ;;
    --output_format)
      OUTPUT_FORMAT="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1"
      usage
      exit 1
      ;;
  esac
done

# ======================
# Validate args
# ======================

if [[ ${#TASKS[@]} -eq 0 ]]; then
  echo "Error: --tasks is required (one or more of: comments features ec cross_references)"
  usage
  exit 1
fi

if [[ -z "${IDS_REL}" ]]; then
  echo "Error: --ids_rel is required (relative path under data/)"
  usage
  exit 1
fi

# Disallow absolute path & traversal
if [[ "${IDS_REL}" = /* ]] || [[ "${IDS_REL}" == *".."* ]]; then
  echo "Error: --ids_rel must be a relative path under data/ (no absolute path, no ..)"
  exit 1
fi

if [[ -n "${OUT_SUBDIR}" ]]; then
  if [[ "${OUT_SUBDIR}" = /* ]] || [[ "${OUT_SUBDIR}" == *".."* ]]; then
    echo "Error: --out_subdir must be a relative path under data/Sample/ (no absolute path, no ..)"
    exit 1
  fi
fi

if [[ "${TEMPLATE_SUBDIR}" = /* ]] || [[ "${TEMPLATE_SUBDIR}" == *".."* ]]; then
  echo "Error: --template_subdir must be a relative directory under data/ (no absolute path, no ..)"
  exit 1
fi

if [[ "${OUTPUT_FORMAT}" != "jsonl" && "${OUTPUT_FORMAT}" != "csv" && "${OUTPUT_FORMAT}" != "tsv" ]]; then
  echo "Error: --output_format must be one of: jsonl | csv | tsv"
  exit 1
fi

# Validate tasks values (keep in bash to fail fast)
for t in "${TASKS[@]}"; do
  case "${t}" in
    comments|features|ec|cross_references) ;;
    *)
      echo "Error: --tasks contains invalid value: ${t} (allowed: comments | features | ec | cross_references)"
      exit 1
      ;;
  esac
done

# ======================
# Compose final paths
# ======================

IDS_FILE="${PROJECT_ROOT}/data/${IDS_REL}"

if [[ -n "${OUT_SUBDIR}" ]]; then
  OUT_DIR="${PROJECT_ROOT}/data/Sample/${OUT_SUBDIR}"
else
  OUT_DIR="${PROJECT_ROOT}/data/Sample"
fi

# Early checks
if [[ ! -f "${IDS_FILE}" ]]; then
  echo "Error: IDS file not found: ${IDS_FILE}"
  exit 1
fi

mkdir -p "${OUT_DIR}"

# ======================
# Templates
# ======================

TEMPLATE_DIR="${PROJECT_ROOT}/data/${TEMPLATE_SUBDIR}"

# Validate tasks values (keep in bash to fail fast)
if [[ ! -d "${TEMPLATE_DIR}" ]]; then
  echo "Error: Template directory not found: ${TEMPLATE_DIR}"
  exit 1
fi

# comments
CA_A="${TEMPLATE_DIR}/catalytic_activity/catalytic_activity_templates_a.json"

CF_A="${TEMPLATE_DIR}/cofactor/cofactor_templates_a.json"
CF_B="${TEMPLATE_DIR}/cofactor/cofactor_templates_b.json"

SL_A="${TEMPLATE_DIR}/subcellular_location/subcellular_location_templates_a.json"
SL_B="${TEMPLATE_DIR}/subcellular_location/subcellular_location_templates_b.json"
SL_C="${TEMPLATE_DIR}/subcellular_location/subcellular_location_templates_c.json"
SL_D="${TEMPLATE_DIR}/subcellular_location/subcellular_location_templates_d.json"

PW_A="${TEMPLATE_DIR}/pathway/pathway_templates_a.json"
PW_B="${TEMPLATE_DIR}/pathway/pathway_templates_b.json"
PW_C="${TEMPLATE_DIR}/pathway/pathway_templates_c.json"
PW_D="${TEMPLATE_DIR}/pathway/pathway_templates_d.json"

# cross_references
GO_BP="${TEMPLATE_DIR}/gene_ontology/gene_ontology_templates_bp.json"
GO_MF="${TEMPLATE_DIR}/gene_ontology/gene_ontology_templates_mf.json"
GO_CC="${TEMPLATE_DIR}/gene_ontology/gene_ontology_templates_cc.json"

# ec
EC_A="${TEMPLATE_DIR}/ec/ec_numbers_templates_a.json"
EC_B="${TEMPLATE_DIR}/ec/ec_numbers_templates_b.json"
EC_C="${TEMPLATE_DIR}/ec/ec_numbers_templates_c.json"
EC_D="${TEMPLATE_DIR}/ec/ec_numbers_templates_d.json"
EC_E="${TEMPLATE_DIR}/ec/ec_numbers_templates_e.json"
EC_F="${TEMPLATE_DIR}/ec/ec_numbers_templates_f.json"
EC_G="${TEMPLATE_DIR}/ec/ec_numbers_templates_g.json"
EC_H="${TEMPLATE_DIR}/ec/ec_numbers_templates_h.json"
EC_I="${TEMPLATE_DIR}/ec/ec_numbers_templates_i.json"
EC_J="${TEMPLATE_DIR}/ec/ec_numbers_templates_j.json"
EC_K="${TEMPLATE_DIR}/ec/ec_numbers_templates_k.json"

# features
RG_A="${TEMPLATE_DIR}/regions/regions_templates_a.json"
RG_B="${TEMPLATE_DIR}/regions/regions_templates_b.json"

ST_A="${TEMPLATE_DIR}/sites/sites_templates_a.json"
ST_B="${TEMPLATE_DIR}/sites/sites_templates_b.json"
ST_C="${TEMPLATE_DIR}/sites/sites_templates_c.json"
ST_D="${TEMPLATE_DIR}/sites/sites_templates_d.json"

# ======================
# Helpers
# ======================

has_task() {
  local want="$1"
  for x in "${TASKS[@]}"; do
    if [[ "$x" == "$want" ]]; then
      return 0
    fi
  done
  return 1
}

# ======================
# Run
# ======================

cd "${PROJECT_ROOT}"

CMD=(python -m runners.run)
CMD+=(--tasks "${TASKS[@]}")
CMD+=(--ids "${IDS_FILE}")
CMD+=(--out_dir "${OUT_DIR}")
CMD+=(--batch_size "${BATCH_SIZE}")
CMD+=(--output_format "${OUTPUT_FORMAT}")

# Append only required template args for selected tasks
if has_task "comments"; then
  CMD+=(--ca_a "${CA_A}")
  CMD+=(--cf_a "${CF_A}" --cf_b "${CF_B}")
  CMD+=(--sl_a "${SL_A}" --sl_b "${SL_B}" --sl_c "${SL_C}" --sl_d "${SL_D}")
  CMD+=(--pw_a "${PW_A}" --pw_b "${PW_B}" --pw_c "${PW_C}" --pw_d "${PW_D}")
fi

if has_task "cross_references"; then
  CMD+=(--go_bp "${GO_BP}" --go_mf "${GO_MF}" --go_cc "${GO_CC}")
fi

if has_task "ec"; then
  CMD+=(--ec_a "${EC_A}" --ec_b "${EC_B}" --ec_c "${EC_C}" --ec_d "${EC_D}")
  CMD+=(--ec_e "${EC_E}" --ec_f "${EC_F}" --ec_g "${EC_G}" --ec_h "${EC_H}")
  CMD+=(--ec_i "${EC_I}" --ec_j "${EC_J}" --ec_k "${EC_K}")
fi

if has_task "features"; then
  CMD+=(--rg_a "${RG_A}" --rg_b "${RG_B}")
  CMD+=(--st_a "${ST_A}" --st_b "${ST_B}" --st_c "${ST_C}" --st_d "${ST_D}")
fi

"${CMD[@]}"