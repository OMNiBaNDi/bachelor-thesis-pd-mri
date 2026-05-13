#!/bin/bash
# Copy Forde BL/3Y/5Y DICOM folders into the thesis staging directories.
set -euo pipefail

SRC_BL=/nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_Baseline/Pasienter/Bergen_Forde_p_anonym
SRC_3Y=/nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_3Y/Pasienter/Forde_p_anonym
SRC_5Y=/nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_5Y/Pasienter/PV_patients_Frde_5Y
DEST_BL=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_BL_DICOM
DEST_3Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_3Y_DICOM
DEST_5Y=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Forde_conversion/Forde_5Y_DICOM

declare -A SRC=( [BL]=$SRC_BL [3Y]=$SRC_3Y [5Y]=$SRC_5Y )
declare -A DEST=( [BL]=$DEST_BL [3Y]=$DEST_3Y [5Y]=$DEST_5Y )

if [[ ! -d "$SRC_5Y" ]]; then
  alt=$(ls -d /nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_5Y/Pasienter/PV_patients_F*de_5Y 2>/dev/null | head -n 1)
  if [[ -n "$alt" ]]; then
    echo "Detected 5Y source directory: $alt"
    SRC_5Y="$alt"
    SRC[5Y]="$SRC_5Y"
  fi
fi

default_ids=(02 04 07 08 11 13 15 20 22 25 27 28 30)
mode="all"
list_file=""

tusage() {
  cat <<'EOF'
Usage: copy_forde_dicoms.sh [--timepoint BL|3Y|5Y] [--patients-file FILE]
Defaults: copy BL, 3Y, and 5Y folders for the 13 curated Forde IDs listed in the script.
  --timepoint   Limit to a single timepoint (repeat the flag to run one at a time)
  --patients-file FILE  Provide a custom list of Forde IDs (one per line, e.g., 02 or F02)
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
  if [[ "$token" =~ ^F([0-9]+)$ ]]; then
    local num="${BASH_REMATCH[1]}"
    printf 'F%02d' "$((10#$num))"
    return 0
  fi
  if [[ "$token" =~ ^([0-9]+)$ ]]; then
    printf 'F%02d' "$((10#$token))"
    return 0
  fi
  printf '%s' "$token"
}

resolve_src_dir() {
  local base="$1"
  local norm_id="$2"  # e.g., F02
  local digits="${norm_id#F}"
  local num=$((10#$digits))
  local candidates=(
    "$norm_id"
    "F${num}"
    "F$(printf '%02d' "$num")"
    "F$(printf '%03d' "$num")"
    "F $(printf '%02d' "$num")"
    "F $(printf '%03d' "$num")"
    "F $(printf '%d' "$num")"
  )
  for cand in "${candidates[@]}"; do
    if [[ -d "$base/$cand" ]]; then
      printf '%s' "$cand"
      return 0
    fi
  done
  local short=$(printf '%d' "$num")
  local found=$(find "$base" -maxdepth 1 -mindepth 1 -type d -iname "F*${short}" -print -quit)
  if [[ -n "$found" ]]; then
    basename "$found"
    return 0
  fi
  return 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --timepoint)
      mode="$2"
      shift 2
      ;;
    --patients-file)
      list_file="$2"
      shift 2
      ;;
    -h|--help)
      tusage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      tusage >&2
      exit 1
      ;;
  esac
done

read_ids=()
if [[ -n "$list_file" ]]; then
  if [[ ! -f "$list_file" ]]; then
    echo "Patient list file not found: $list_file" >&2
    exit 1
  fi
  mapfile -t raw_ids < "$list_file"
  for entry in "${raw_ids[@]}"; do
    entry="$(trim "$entry")"
    [[ -z "$entry" ]] && continue
    read_ids+=("$entry")
  done
else
  read_ids=("${default_ids[@]}")
fi

select_timepoints=()
if [[ "$mode" == "all" ]]; then
  select_timepoints=(BL 3Y 5Y)
else
  select_timepoints=("$mode")
fi

copy_count=0
skip_count=0

for tp in "${select_timepoints[@]}"; do
  src_root="${SRC[$tp]}"
  dest_root="${DEST[$tp]}"
  if [[ -z "$src_root" || -z "$dest_root" ]]; then
    echo "WARN: undefined paths for timepoint $tp" >&2
    continue
  fi
  mkdir -p "$dest_root"
  for pid in "${read_ids[@]}"; do
    if ! norm_id="$(normalize "$pid")"; then
      echo "Skipping invalid ID: $pid" >&2
      skip_count=$((skip_count + 1))
      continue
    fi
    if ! src_rel="$(resolve_src_dir "$src_root" "$norm_id")"; then
      echo "Missing $src_root/$norm_id (tp=$tp)"
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
