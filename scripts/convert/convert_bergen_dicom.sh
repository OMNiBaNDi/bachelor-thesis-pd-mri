#!/bin/bash
# Convert Bergen BL and/or 3Y DICOM folders to NIfTI using dcm2niix.
# Run this on the UIS cluster where /nfs/br1_prosjekt/... is mounted.
set -euo pipefail

DCM2NIIX=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/tools/dcm2niix/dcm2niix
SRC_BL=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/BergenConversion/Bergen_BL_DICOM
DEST_BL=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/BergenConversion/Bergen_BL_converted_to_nifti
SRC_3Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/BergenConversion/Bergen_3Y_DICOM
DEST_3Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/BergenConversion/Bergen_3Y_converted_to_nifti
MODE="both"  # bl | 3y | both
LIST_FILE=""

DEFAULT_PATIENT_IDS=(
  05 14 15 17 20 21 24 25 26 28 29 33 34 35 36 38
  39 40 43 53 54 55 56 57 58 63 71 72 74 83 84 90
  92 94 95 101
)

usage() {
  cat <<'EOF'
Usage: convert_bergen_dicom.sh [--mode bl|3y|both] [--patients-file FILE]

Default: convert both BL and 3Y for the curated 36 Bergen patients.
  --mode bl     Convert only the BL folders.
  --mode 3y     Convert only the 3Y folders.
  --mode both   Convert both BL and 3Y (default).
  --patients-file FILE  Override the patient list (one ID per line, e.g., 05 or B005).
EOF
}

trim() {
  local s="$1"
  s="${s#${s%%[!$'\t\r\n ']*}}"
  s="${s%${s##*[!$'\t\r\n ']}}"
  printf '%s' "$s"
}

norm_id() {
  local token="$(trim "$1")"
  [[ -n "$token" ]] || return 1
  token="${token//,/}"
  if [[ "$token" =~ ^B[[:space:]]*([0-9]+)$ ]]; then
    token="${BASH_REMATCH[1]}"
  elif [[ "$token" =~ ^B([0-9]+)$ ]]; then
    token="${BASH_REMATCH[1]}"
  fi
  if [[ "$token" =~ ^[0-9]+$ ]]; then
    local num=$((10#$token))
    printf '%03d' "$num"
    return 0
  fi
  if [[ "$token" =~ ^[0-9]{3}$ ]]; then
    printf '%s' "$token"
    return 0
  fi
  return 1
}

resolve_bl_src() {
  local id3="$1"
  local num=$((10#$id3))
  local candidates=(
    "B $num"
    "B $(printf '%02d' "$num")"
    "B $(printf '%03d' "$num")"
    "B$num"
    "B$(printf '%02d' "$num")"
    "B$(printf '%03d' "$num")"
  )
  for cand in "${candidates[@]}"; do
    if [[ -d "$SRC_BL/$cand" ]]; then
      printf '%s' "$cand"
      return 0
    fi
  done
  return 1
}

resolve_3y_src() {
  local id3="$1"
  local num=$((10#$id3))
  local candidates=(
    "B$(printf '%03d' "$num")"
    "B$(printf '%02d' "$num")"
    "B$num"
    "B $(printf '%03d' "$num")"
  )
  for cand in "${candidates[@]}"; do
    if [[ -d "$SRC_3Y/$cand" ]]; then
      printf '%s' "$cand"
      return 0
    fi
  done
  return 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="$2"
      shift 2
      ;;
    -f|--patients-file)
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

case "$MODE" in
  bl|BL)
    MODE="bl"
    ;;
  3y|Y3|3Y)
    MODE="3y"
    ;;
  both|ALL|Both)
    MODE="both"
    ;;
  *)
    echo "Invalid mode: $MODE" >&2
    usage >&2
    exit 1
    ;;
esac

PATIENT_IDS=()
if [[ -n "$LIST_FILE" ]]; then
  if [[ ! -f "$LIST_FILE" ]]; then
    echo "Patient list file not found: $LIST_FILE" >&2
    exit 1
  fi
  mapfile -t RAW_IDS < "$LIST_FILE"
  for entry in "${RAW_IDS[@]}"; do
    entry="$(trim "$entry")"
    [[ -z "$entry" ]] && continue
    PATIENT_IDS+=("$entry")
  done
else
  PATIENT_IDS=("${DEFAULT_PATIENT_IDS[@]}")
fi

convert_patient() {
  local id3="$1" tp="$2"
  local src_dir dest_root src_root mode_label
  case "$tp" in
    BL)
      src_root="$SRC_BL"
      dest_root="$DEST_BL"
      src_dir="$(resolve_bl_src "$id3")" || {
        printf 'BL source missing for B%s\n' "$id3"
        return 1
      }
      mode_label="BL"
      ;;
    3Y)
      src_root="$SRC_3Y"
      dest_root="$DEST_3Y"
      src_dir="$(resolve_3y_src "$id3")" || {
        printf '3Y source missing for B%s\n' "$id3"
        return 1
      }
      mode_label="3Y"
      ;;
    *)
      echo "Unknown timepoint: $tp" >&2
      return 1
      ;;
  esac
  local dest_dir="B$id3"
  mkdir -p "$dest_root/$dest_dir"
  printf 'Converting %s (%s, patient B%s) → %s\n' "$src_root/$src_dir" "$mode_label" "$id3" "$dest_root/$dest_dir"
  if ! "$DCM2NIIX" -z y -o "$dest_root/$dest_dir" "$src_root/$src_dir"; then
    printf 'dcm2niix failed for %s (%s)\n' "$dest_dir" "$mode_label" >&2
    return 1
  fi
}

mkdir -p "$DEST_BL" "$DEST_3Y"

failures=0
for raw_id in "${PATIENT_IDS[@]}"; do
  if ! id3="$(norm_id "$raw_id")"; then
    printf 'Skipping %s (unable to normalize ID)\n' "$raw_id"
    failures=$((failures + 1))
    continue
  fi
  if [[ "$MODE" == "bl" || "$MODE" == "both" ]]; then
    if ! convert_patient "$id3" BL; then
      failures=$((failures + 1))
    fi
  fi
  if [[ "$MODE" == "3y" || "$MODE" == "both" ]]; then
    if ! convert_patient "$id3" 3Y; then
      failures=$((failures + 1))
    fi
  fi
done

if [[ $failures -gt 0 ]]; then
  printf 'Completed with %d issue(s). Review the log above for details.\n' "$failures" >&2
  exit 1
else
  printf 'Conversion finished for mode=%s.\n' "$MODE"
fi
