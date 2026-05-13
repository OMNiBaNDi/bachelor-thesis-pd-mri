#!/bin/bash
# Convert Forde BL/3Y/5Y DICOM staging folders to NIfTI via dcm2niix.
set -euo pipefail

DCM2NIIX=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/tools/dcm2niix/dcm2niix
SRC_BL=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_BL_DICOM
SRC_3Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_3Y_DICOM
SRC_5Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_5Y_DICOM
DEST_BL=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_BL_converted_to_nifti
DEST_3Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_3Y_converted_to_nifti
DEST_5Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_5Y_converted_to_nifti

PATIENT_IDS=(02 04 07 08 11 13 15 20 22 25 27 28 30)
MODE="all"
LIST_FILE=""

trim() {
  local s="$1"
  s="${s#${s%%[!$'\t\r\n ']*}}"
  s="${s%${s##*[!$'\t\r\n ']}}"
  printf '%s' "$s"
}

normalize() {
  local token="$(trim "$1")"
  [[ -n "$token" ]] || return 1
  token="${token//,/}"
  token="${token^^}"
  if [[ "$token" =~ ^F([0-9]+)$ ]]; then
    printf 'F%02d' "$((10#${BASH_REMATCH[1]}))"
    return 0
  fi
  if [[ "$token" =~ ^([0-9]+)$ ]]; then
    printf 'F%02d' "$((10#$token))"
    return 0
  fi
  printf '%s' "$token"
}

usage() {
  cat <<'EOF'
Usage: convert_forde_dicom.sh [--mode BL|3Y|5Y|all] [--patients-file FILE]
Default: convert BL, 3Y, and 5Y staging folders for the curated Forde IDs.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="$2"
      shift 2
      ;;
    --patients-file)
      LIST_FILE="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

IDS=()
if [[ -n "$LIST_FILE" ]]; then
  if [[ ! -f "$LIST_FILE" ]]; then
    echo "Patient list file not found: $LIST_FILE" >&2
    exit 1
  fi
  mapfile -t RAW_IDS < "$LIST_FILE"
  for entry in "${RAW_IDS[@]}"; do
    entry="$(trim "$entry")"
    [[ -z "$entry" ]] && continue
    IDS+=("$entry")
  done
else
  IDS=("${PATIENT_IDS[@]}")
fi

case "$MODE" in
  all|ALL) TPs=(BL 3Y 5Y) ;;
  BL|3Y|5Y) TPs=("$MODE") ;;
  *) echo "Invalid mode: $MODE" >&2; exit 1 ;;
esac

declare -A SRC=( [BL]=$SRC_BL [3Y]=$SRC_3Y [5Y]=$SRC_5Y )
declare -A DEST=( [BL]=$DEST_BL [3Y]=$DEST_3Y [5Y]=$DEST_5Y )

mkdir -p "$DEST_BL" "$DEST_3Y" "$DEST_5Y"

failures=0
completed=0

for tp in "${TPs[@]}"; do
  src_root="${SRC[$tp]}"
  dest_root="${DEST[$tp]}"
  if [[ ! -d "$src_root" ]]; then
    echo "WARN: source root missing for $tp ($src_root)"
    continue
  fi
  for raw_id in "${IDS[@]}"; do
    if ! pid="$(normalize "$raw_id")"; then
      echo "Skipping invalid ID: $raw_id"
      continue
    fi
    src_dir="$src_root/$pid"
    dest_dir="$dest_root/$pid"
    if [[ ! -d "$src_dir" ]]; then
      echo "Missing $src_dir"
      failures=$((failures + 1))
      continue
    fi
    rm -rf "$dest_dir"
    mkdir -p "$dest_dir"
    mapfile -t dicom_dirs < <(find "$src_dir" -type f ! -name '*.json' ! -name '*.txt' -printf '%h\n' | LC_ALL=C sort -u)
    if [[ ${#dicom_dirs[@]} -eq 0 ]]; then
      echo "No DICOM files found under $src_dir"
      failures=$((failures + 1))
      continue
    fi
    for dicom_dir in "${dicom_dirs[@]}"; do
      echo "[${tp}] Converting $dicom_dir → $dest_dir"
      if ! "$DCM2NIIX" -z y -o "$dest_dir" "$dicom_dir" >/dev/null 2>&1; then
        echo "dcm2niix failed for $dicom_dir ($pid $tp)"
        failures=$((failures + 1))
      else
        completed=$((completed + 1))
      fi
    done
  done
done

printf 'Conversion finished: %d successful directory conversions, %d issues.\n' "$completed" "$failures"
