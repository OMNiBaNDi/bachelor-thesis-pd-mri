#!/bin/bash
# Convert Forde healthy-control BL/3Y/5Y DICOM folders (FK01) to NIfTI using dcm2niix.
# Usage: convert_forde_healthy_dicom.sh [--mode bl|3y|5y|all]

set -euo pipefail

DCM2NIIX=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/tools/dcm2niix/dcm2niix

SRC_BL=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_healthy_BL_DICOM
DEST_BL=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_healthy_BL_converted_to_nifti
SRC_3Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_healthy_3Y_DICOM
DEST_3Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_healthy_3Y_converted_to_nifti
SRC_5Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_healthy_5Y_DICOM
DEST_5Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_healthy_5Y_converted_to_nifti

DEFAULT_PATIENT_IDS=(01)
MODE="all"

usage() {
  cat <<'EOF'
Usage: convert_forde_healthy_dicom.sh [--mode bl|3y|5y|all]
Defaults: convert BL, 3Y, and 5Y for FK01.
EOF
}

trim() {
  local s="$1"
  s="${s#${s%%[!$'\t\r\n ']*}}"
  s="${s%${s##*[!$'\t\r\n ']}}"
  printf '%s' "$s"
}

normalize_id() {
  local token="$(trim "$1")"
  [[ -n "$token" ]] || return 1
  token="${token^^}"
  token="${token//,/}"
  token="${token// /}"
  token="${token#FK}"
  if [[ "$token" =~ ^([0-9]+)$ ]]; then
    local num=$((10#$token))
    printf 'FK%02d' "$num"
    return 0
  fi
  printf '%s' "$token"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="$2"
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

case "$MODE" in
  bl|BL)
    TARGET_TPS=(BL)
    ;;
  3y|3Y|Y3)
    TARGET_TPS=(3Y)
    ;;
  5y|5Y|Y5)
    TARGET_TPS=(5Y)
    ;;
  all|ALL)
    TARGET_TPS=(BL 3Y 5Y)
    ;;
  *)
    echo "Invalid mode: $MODE" >&2
    usage >&2
    exit 1
    ;;
esac

declare -A SRC=( [BL]="$SRC_BL" [3Y]="$SRC_3Y" [5Y]="$SRC_5Y" )
declare -A DEST=( [BL]="$DEST_BL" [3Y]="$DEST_3Y" [5Y]="$DEST_5Y" )

convert_patient() {
  local id="$1" tp="$2"
  local src_root="${SRC[$tp]}"
  local dest_root="${DEST[$tp]}"
  local src_dir="$src_root/$id"
  local dest_dir="$dest_root/$id"

  if [[ ! -d "$src_dir" ]]; then
    printf 'Missing staging folder %s for %s (tp=%s)\n' "$src_dir" "$id" "$tp" >&2
    return 1
  fi
  mkdir -p "$dest_dir"
  printf 'Converting %s (%s) → %s\n' "$src_dir" "$tp" "$dest_dir"
  if ! "$DCM2NIIX" -z y -o "$dest_dir" "$src_dir"; then
    printf 'dcm2niix failed for %s (%s)\n' "$id" "$tp" >&2
    return 1
  fi
}

mkdir -p "$DEST_BL" "$DEST_3Y" "$DEST_5Y"

failures=0
for raw_id in "${DEFAULT_PATIENT_IDS[@]}"; do
  if ! norm="$(normalize_id "$raw_id")"; then
    echo "Skipping invalid ID: $raw_id" >&2
    failures=$((failures + 1))
    continue
  fi
  for tp in "${TARGET_TPS[@]}"; do
    if ! convert_patient "$norm" "$tp"; then
      failures=$((failures + 1))
    fi
  done
done

if [[ $failures -gt 0 ]]; then
  echo "Completed with $failures issue(s). Review the log above." >&2
  exit 1
else
  echo "Conversion finished for mode=$MODE."
fi
