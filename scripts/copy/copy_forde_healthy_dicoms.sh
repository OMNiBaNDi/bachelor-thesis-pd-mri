#!/bin/bash
# Copy Forde healthy-control BL/3Y/5Y DICOM folders (FK01) into dedicated staging directories.
# Usage: copy_forde_healthy_dicoms.sh [--timepoint BL|3Y|5Y] [--patients-file FILE]

set -euo pipefail

SRC_BL=/nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_Baseline/Kontroller/Bergen_Forde_k_anonym
SRC_3Y=/nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_3Y/Kontroller/Forde_k_anonym
SRC_5Y=/nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_5Y/Kontroller/PV_controls_Frde_5Y
DEST_BL=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_healthy_BL_DICOM
DEST_3Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_healthy_3Y_DICOM
DEST_5Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_healthy_5Y_DICOM

if [[ ! -d "$SRC_5Y" ]]; then
  alt=$(ls -d /nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_5Y/Kontroller/PV_controls_F*de_5Y 2>/dev/null | head -n 1)
  if [[ -n "$alt" ]]; then
    echo "Detected 5Y source directory: $alt"
    SRC_5Y="$alt"
  fi
fi

DEFAULT_PATIENT_IDS=(01)
MODE="all"
LIST_FILE=""

usage() {
  cat <<'EOF'
Usage: copy_forde_healthy_dicoms.sh [--timepoint BL|3Y|5Y] [--patients-file FILE]
Defaults: copy all three timepoints for FK01.
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

resolve_src_dir() {
  local base="$1"
  local norm_id="$2"  # e.g., FK01
  local tp="$3"
  local digits="${norm_id#FK}"
  local num=$((10#$digits))
  local candidates=()
  case "$tp" in
    BL)
      candidates+=("FK $(printf '%02d' "$num")")
      candidates+=("FK $(printf '%d' "$num")")
      candidates+=("FK$(printf '%02d' "$num")")
      ;;
    3Y|5Y)
      candidates+=("FK$(printf '%02d' "$num")")
      candidates+=("FK$(printf '%d' "$num")")
      candidates+=("FK $(printf '%02d' "$num")")
      ;;
    *)
      candidates+=("FK$(printf '%02d' "$num")")
      ;;
  esac
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
  IDS=("${DEFAULT_PATIENT_IDS[@]}")
fi

case "$MODE" in
  all|ALL)
    SELECT_TPS=(BL 3Y 5Y)
    ;;
  BL|3Y|5Y)
    SELECT_TPS=("$MODE")
    ;;
  *)
    echo "Invalid timepoint: $MODE" >&2
    usage >&2
    exit 1
    ;;
esac

declare -A SRC=( [BL]="$SRC_BL" [3Y]="$SRC_3Y" [5Y]="$SRC_5Y" )
declare -A DEST=( [BL]="$DEST_BL" [3Y]="$DEST_3Y" [5Y]="$DEST_5Y" )

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
    if ! src_rel="$(resolve_src_dir "$src_root" "$norm_id" "$tp")"; then
      echo "Missing $src_root/$norm_id (tp=$tp)" >&2
      skip_count=$((skip_count + 1))
      continue
    fi
    src_dir="$src_root/$src_rel"
    dest_dir="$dest_root/$norm_id"
    printf 'Copying %s (%s) → %s\n' "$src_dir" "$tp" "$dest_dir"
    rsync -a --info=progress2 --human-readable "$src_dir/" "$dest_dir/"
    copy_count=$((copy_count + 1))
  done
done

echo "Done. Copied $copy_count folders (skipped $skip_count)."
