#!/bin/bash
# Remove contents of data/SUS/Sxxx/{BL,3Y,5Y} for the curated SUS ID list.
set -euo pipefail

DATA_ROOT=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS
PATIENT_IDS=(02 05 06 08 09 10 11 13 14 18 22 26 28 30 32 34 38 41 42 43 44 45 46 49 52)
TIMEPOINTS=(BL 3Y 5Y)
LIST_FILE=""

usage() {
  cat <<'EOF'
Usage: purge_sus_timepoints.sh [--patients-file FILE]
Clears data/SUS/Sxxx/{BL,3Y,5Y} for the 25 SUS patients.
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
    printf 'S%03d' "$((10#${BASH_REMATCH[1]}))"
    return 0
  fi
  if [[ "$token" =~ ^([0-9]+)$ ]]; then
    printf 'S%03d' "$((10#$token))"
    return 0
  fi
  printf '%s' "$token"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
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

cleaned=0
missing=0

for raw_id in "${IDS[@]}"; do
  if ! pid="$(normalize "$raw_id")"; then
    echo "Skipping invalid ID: $raw_id"
    continue
  fi
  for tp in "${TIMEPOINTS[@]}"; do
    target="$DATA_ROOT/$pid/$tp"
    if [[ -d "$target" ]]; then
      if find "$target" -mindepth 1 -print -quit >/dev/null 2>&1; then
        echo "Clearing $target"
        find "$target" -mindepth 1 -delete
      else
        echo "Already empty: $target"
      fi
      cleaned=$((cleaned + 1))
    else
      echo "Missing folder: $target"
      missing=$((missing + 1))
    fi
  done
done

printf 'Done. Processed %d folders (missing %d).\n' "$cleaned" "$missing"
