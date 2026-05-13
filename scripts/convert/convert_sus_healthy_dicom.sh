#!/bin/bash
# Convert SUS healthy-control 3Y/5Y DICOM staging folders to NIfTI (handles nested subfolders).
# Usage: convert_sus_healthy_dicom.sh [--mode 3y|5y|all] [--patients-file FILE]

set -euo pipefail

DCM2NIIX=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/tools/dcm2niix/dcm2niix
SRC_3Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_healthy_3Y_DICOM
DEST_3Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_healthy_3Y_converted_to_nifti
SRC_5Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_healthy_5Y_DICOM
DEST_5Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_healthy_5Y_converted_to_nifti

DEFAULT_IDS=(05 08 10 11 13 15 16 23 26 27 28 31 34 35 37 41 46 47)
MODE="all"
LIST_FILE=""

usage() {
  cat <<'EOF'
Usage: convert_sus_healthy_dicom.sh [--mode 3y|5y|all] [--patients-file FILE]
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
  mapfile -t RAW < "$LIST_FILE"
  for entry in "${RAW[@]}"; do
    entry="$(trim "$entry")"
    [[ -z "$entry" ]] && continue
    IDS+=("$entry")
  done
else
  IDS=("${DEFAULT_IDS[@]}")
fi

case "$MODE" in
  3y|3Y|Y3)
    TPs=(3Y)
    ;;
  5y|5Y|Y5)
    TPs=(5Y)
    ;;
  all|ALL)
    TPs=(3Y 5Y)
    ;;
  *)
    echo "Invalid mode: $MODE" >&2
    usage >&2
    exit 1
    ;;
esac

declare -A SRC=( [3Y]="$SRC_3Y" [5Y]="$SRC_5Y" )
declare -A DEST=( [3Y]="$DEST_3Y" [5Y]="$DEST_5Y" )

declare -A TMPDIRS

delete_tmpdirs() {
  for tmp in "${TMPDIRS[@]}"; do
    [[ -n "$tmp" && -d "$tmp" ]] && rm -rf "$tmp"
  done
}

trap delete_tmpdirs EXIT

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
  rm -rf "$dest_dir"
  mkdir -p "$dest_dir"
  tmpdir=$(mktemp -d)
  TMPDIRS["$id$tp"]="$tmpdir"

  echo "[$tp] Flattening $src_dir → $tmpdir"
  find "$src_dir" -type f ! -name '*.json' ! -name '*.txt' -print0 | while IFS= read -r -d '' file; do
    base=$(basename "${file%.*}")
    cp "$file" "$tmpdir/${base}_$(date +%s%N)${file##*.}" 2>/dev/null || cp "$file" "$tmpdir/${base}_$(date +%s%N)"
  done

  if ! "$DCM2NIIX" -z y -o "$dest_dir" "$tmpdir"; then
    echo "dcm2niix failed for $src_dir ($id $tp)"
    return 1
  fi
}

failures=0
for raw_id in "${IDS[@]}"; do
  if ! norm="$(normalize_id "$raw_id")"; then
    echo "Skipping invalid ID: $raw_id" >&2
    continue
  fi
  for tp in "${TPs[@]}"; do
    if ! convert_patient "$norm" "$tp"; then
      failures=$((failures + 1))
    fi
  done
done

delete_tmpdirs

if [[ $failures -gt 0 ]]; then
  echo "Completed with $failures issue(s). Review the log above." >&2
  exit 1
else
  echo "Conversion finished for mode=$MODE."
fi
