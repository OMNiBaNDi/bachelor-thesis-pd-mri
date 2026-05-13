#!/bin/bash
# Delete the contents of BL/3Y/5Y folders for the curated Bergen healthy-control IDs.
# Usage: purge_bergen_healthy_timepoints.sh [--patients-file FILE] [--timepoints "BL 3Y"]
# Defaults: clean BL, 3Y, and 5Y for the 17 standard BK healthy IDs listed below.

set -euo pipefail

DATA_ROOT=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Bergen_healthy
DEFAULT_PATIENT_IDS=(17 21 23 24 26 28 30 31 34 39 47 51 52 53 56 57 66)
DEFAULT_TIMEPOINTS=(BL 3Y 5Y)
LIST_FILE=""
TIMEPOINT_ARGS=()

usage() {
  cat <<'EOF'
Usage: purge_bergen_healthy_timepoints.sh [--patients-file FILE] [--timepoints "BL 3Y 5Y"]

Options:
  --patients-file FILE   limit to IDs listed in FILE (one per line: e.g., 17 or BK17)
  --timepoints "BL 5Y"   override which timepoints to purge (space-separated list)
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
  token="${token^^}"
  token="${token//,/}"
  token="${token// /}"
  token="${token#BK}"
  if [[ "$token" =~ ^([0-9]+)$ ]]; then
    local num=$((10#$token))
    printf 'BK%02d' "$num"
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
    --timepoints)
      shift
      while [[ $# -gt 0 && "$1" != --* ]]; do
        TIMEPOINT_ARGS+=("$1")
        shift
      done
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
  mapfile -t RAW < "$LIST_FILE"
  for entry in "${RAW[@]}"; do
    entry="$(trim "$entry")"
    [[ -z "$entry" ]] && continue
    PATIENT_IDS+=("$entry")
  done
else
  PATIENT_IDS=("${DEFAULT_PATIENT_IDS[@]}")
fi

TIMEPOINTS=()
if [[ ${#TIMEPOINT_ARGS[@]} -gt 0 ]]; then
  for tp in "${TIMEPOINT_ARGS[@]}"; do
    case "$tp" in
      BL|3Y|5Y)
        TIMEPOINTS+=("$tp")
        ;;
      *)
        echo "Invalid timepoint: $tp" >&2
        exit 1
        ;;
    esac
  done
else
  TIMEPOINTS=("${DEFAULT_TIMEPOINTS[@]}")
fi

cleaned=0
skipped=0

for raw_id in "${PATIENT_IDS[@]}"; do
  if ! patient_dir="$(normalize_patient_dir "$raw_id")"; then
    echo "Skipping invalid ID: $raw_id" >&2
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

echo "Done. Processed $cleaned folders (skipped $skipped)."
