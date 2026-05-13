#!/bin/bash
# Copy SUS BL/3Y/5Y DICOM folders into the SUS conversion staging directories.
set -euo pipefail

SRC_BL=/nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_Baseline/Pasienter/Stavanger_p_anonym
SRC_3Y=/nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_3Y/Pasienter/Stavanger_p_anonym
SRC_5Y=/nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_5Y/Pasienter/PV_patients_SUS_5Y
DEST_BL=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_BL_DICOM
DEST_3Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_3Y_DICOM
DEST_5Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_5Y_DICOM

PATIENT_IDS=(02 05 06 08 09 10 11 13 14 18 22 26 28 30 32 34 38 41 42 43 44 45 46 49 52)
MODE="all"
LIST_FILE=""

usage() {
  cat <<'EOF'
Usage: copy_sus_dicoms.sh [--timepoint BL|3Y|5Y] [--patients-file FILE]
Defaults: copy all three timepoints for the 25 curated SUS IDs.
EOF
}

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
  if [[ "$token" =~ ^S([0-9]+)$ ]]; then
    printf 'S%02d' "$((10#${BASH_REMATCH[1]}))"
    return 0
  fi
  if [[ "$token" =~ ^([0-9]+)$ ]]; then
    printf 'S%02d' "$((10#$token))"
    return 0
  fi
  printf '%s' "$token"
}

resolve_src_dir() {
  local base="$1"
  local norm_id="$2"
  local digits="${norm_id#S}"
  local num=$((10#$digits))
  local candidates=(
    "$norm_id"
    "S${num}"
    "S$(printf '%02d' "$num")"
    "S$(printf '%03d' "$num")"
    "S $(printf '%02d' "$num")"
    "S $(printf '%03d' "$num")"
    "S $(printf '%d' "$num")"
  )
  for cand in "${candidates[@]}"; do
    if [[ -d "$base/$cand" ]]; then
      printf '%s' "$cand"
      return 0
    fi
  done
  local short=$(printf '%d' "$num")
  local found=$(find "$base" -maxdepth 1 -mindepth 1 -type d -iname "S*${short}" -print -quit)
  if [[ -n "$found" ]]; then
    basename "$found"
    return 0
  fi
  return 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --timepoint)
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
  mapfile -t RAW < "$LIST_FILE"
  for entry in "${RAW[@]}"; do
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
  *) echo "Invalid timepoint: $MODE" >&2; exit 1 ;;
esac

declare -A SRC=( [BL]=$SRC_BL [3Y]=$SRC_3Y [5Y]=$SRC_5Y )
declare -A DEST=( [BL]=$DEST_BL [3Y]=$DEST_3Y [5Y]=$DEST_5Y )

copy_count=0
skip_count=0

for tp in "${TPs[@]}"; do
  src_root="${SRC[$tp]}"
  dest_root="${DEST[$tp]}"
  if [[ ! -d "$src_root" ]]; then
    echo "WARN: source root missing for $tp ($src_root)"
    continue
  fi
  mkdir -p "$dest_root"
  for raw_id in "${IDS[@]}"; do
    if ! stage_id="$(normalize "$raw_id")"; then
      echo "Skipping invalid ID: $raw_id"
      continue
    fi
    if ! src_rel="$(resolve_src_dir "$src_root" "$stage_id")"; then
      echo "Missing $src_root/$stage_id (tp=$tp)"
      skip_count=$((skip_count + 1))
      continue
    fi
    src_dir="$src_root/$src_rel"
    dest_dir="$dest_root/$stage_id"
    printf 'Copying %s (%s) → %s\n' "$src_dir" "$tp" "$dest_dir"
    rsync -a --info=progress2 --human-readable "$src_dir/" "$dest_dir/"
    copy_count=$((copy_count + 1))
  done
done

echo "Done. Copied $copy_count folders (skipped $skip_count)."
