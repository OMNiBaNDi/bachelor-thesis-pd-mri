#!/bin/bash
# Copy SUS healthy-control 3Y/5Y DICOM folders into staging directories.

set -euo pipefail

SRC_3Y=/nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_3Y/Kontroller/Stavanger_k_anonym
SRC_5Y=/nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_5Y/Kontroller/PV_controls_SUS_5Y
DEST_3Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_healthy_3Y_DICOM
DEST_5Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_healthy_5Y_DICOM

PATIENT_IDS=(05 08 10 11 13 15 16 23 26 27 28 31 34 35 37 41 46 47)
MODE="all"
LIST_FILE=""

usage() {
  cat <<'EOF'
Usage: copy_sus_healthy_dicoms.sh [--timepoint 3Y|5Y] [--patients-file FILE]
Defaults: copy both 3Y and 5Y for the 18 SUS healthy IDs listed in the script.
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
  token="${token#SK}"
  if [[ "$token" =~ ^([0-9]+)$ ]]; then
    local num=$((10#$token))
    printf 'SK%02d' "$num"
    return 0
  fi
  printf '%s' "$token"
}

resolve_src_dir() {
  local base="$1"
  local norm_id="$2"
  local digits="${norm_id#SK}"
  local num=$((10#$digits))
  local candidates=(
    "$norm_id"
    "SK$(printf '%d' "$num")"
    "SK $(printf '%02d' "$num")"
    "SK $(printf '%d' "$num")"
  )
  for cand in "${candidates[@]}"; do
    if [[ -d "$base/$cand" ]]; then
      printf '%s' "$cand"
      return 0
    fi
  done
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
  all|ALL)
    SELECT_TPS=(3Y 5Y)
    ;;
  3Y|5Y)
    SELECT_TPS=("$MODE")
    ;;
  *)
    echo "Invalid timepoint: $MODE" >&2
    usage >&2
    exit 1
    ;;
esac

declare -A SRC=( [3Y]="$SRC_3Y" [5Y]="$SRC_5Y" )
declare -A DEST=( [3Y]="$DEST_3Y" [5Y]="$DEST_5Y" )

copy_count=0
skip_count=0

for tp in "${SELECT_TPS[@]}"; do
  src_root="${SRC[$tp]}"
  dest_root="${DEST[$tp]}"
  if [[ ! -d "$src_root" ]]; then
    echo "WARN: source root missing for $tp ($src_root)"
    continue
  fi
  mkdir -p "$dest_root"
  for raw_id in "${IDS[@]}"; do
    if ! norm_id="$(normalize_id "$raw_id")"; then
      echo "Skipping invalid ID: $raw_id" >&2
      skip_count=$((skip_count + 1))
      continue
    fi
    if ! src_rel="$(resolve_src_dir "$src_root" "$norm_id")"; then
      echo "Missing $src_root/$norm_id (tp=$tp)" >&2
      skip_count=$((skip_count + 1))
      continue
    fi
    src_dir="$src_root/$src_rel"
    dest_dir="$dest_root/$norm_id"
    printf 'Copying %s (%s) → %s\n' "$src_dir" "$tp" "$dest_dir"
    rsync -a --info=progress2 "$src_dir/" "$dest_dir/"
    copy_count=$((copy_count + 1))
  done
done

echo "Done. Copied $copy_count folders (skipped $skip_count)."
