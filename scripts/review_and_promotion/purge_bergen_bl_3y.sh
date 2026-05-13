#!/bin/bash
# Delete the existing BL and 3Y contents for the curated Bergen patients.
# Run this on the UIS cluster where /nfs/br1_prosjekt/... is mounted.
set -euo pipefail

DATA_ROOT=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Bergen
LIST_FILE=""
TIMEPOINTS=(BL 3Y)

DEFAULT_PATIENT_IDS=(
  05 14 15 17 20 21 24 25 26 28 29 33 34 35 36 38
  39 40 43 53 54 55 56 57 58 63 71 72 74 83 84 90
  92 94 95 101
)

usage() {
  cat <<'EOF'
Usage: purge_bergen_bl_3y.sh [--patients-file FILE]

Default (no flags): empties the BL and 3Y folders for the curated 36 Bergen patients.
  --patients-file FILE  Only clean the IDs listed in FILE (one per line).
EOF
}

trim() {
  local s="$1"
  s="${s#${s%%[!$'\t\r\n ']*}}"
  s="${s%${s##*[!$'\t\r\n ']}}"
  printf '%s' "$s"
}

normalize_patient_dir() {
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
    printf 'B%03d' "$num"
    return 0
  fi
  printf '%s' "$token"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
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

cleaned=0
skipped=0
for id in "${PATIENT_IDS[@]}"; do
  if ! patient_dir="$(normalize_patient_dir "$id")"; then
    printf 'Skipping %s (unable to normalize patient ID)\n' "$id"
    skipped=$((skipped + 1))
    continue
  fi
  for tp in "${TIMEPOINTS[@]}"; do
    target="$DATA_ROOT/$patient_dir/$tp"
    if [[ -d "$target" ]]; then
      if find "$target" -mindepth 1 -print -quit >/dev/null 2>&1; then
        printf 'Clearing %s/%s\n' "$patient_dir" "$tp"
        find "$target" -mindepth 1 -delete
      else
        printf 'Already empty: %s/%s\n' "$patient_dir" "$tp"
      fi
      cleaned=$((cleaned + 1))
    else
      printf 'Missing folder: %s/%s (skipping)\n' "$patient_dir" "$tp"
      skipped=$((skipped + 1))
    fi
  done
done

printf 'Done. Processed %d folders (skipped %d missing/unmapped entries).\n' "$cleaned" "$skipped"
