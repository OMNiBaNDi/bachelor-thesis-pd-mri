#!/bin/bash
# Purge SUS healthy-control 3Y/5Y folders under data/SUS_healthy/SKxx/{3Y,5Y}.

set -euo pipefail

DATA_ROOT=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_healthy
DEFAULT_IDS=(05 08 10 11 13 15 16 23 26 27 28 31 34 35 37 41 46 47)
DEFAULT_TPS=(3Y 5Y)
LIST_FILE=""
TP_ARGS=()

usage() {
  cat <<'EOF'
Usage: purge_sus_healthy_timepoints.sh [--patients-file FILE] [--timepoints "3Y 5Y"]
Defaults: purge both 3Y and 5Y directories for the 18 SUS healthy IDs.
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
    --patients-file)
      LIST_FILE="$2"
      shift 2
      ;;
    --timepoints)
      shift
      while [[ $# -gt 0 && "$1" != --* ]]; do
        TP_ARGS+=("$1")
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

TIMEPOINTS=()
if [[ ${#TP_ARGS[@]} -gt 0 ]]; then
  for tp in "${TP_ARGS[@]}"; do
    case "$tp" in
      3Y|5Y) TIMEPOINTS+=("$tp") ;;
      *) echo "Invalid timepoint: $tp" >&2; exit 1 ;;
    esac
  done
else
  TIMEPOINTS=("${DEFAULT_TPS[@]}")
fi

cleaned=0
skipped=0
for raw_id in "${IDS[@]}"; do
  if ! norm_id="$(normalize_id "$raw_id")"; then
    echo "Skipping invalid ID: $raw_id" >&2
    skipped=$((skipped + 1))
    continue
  fi
  for tp in "${TIMEPOINTS[@]}"; do
    target="$DATA_ROOT/$norm_id/$tp"
    if [[ -d "$target" ]]; then
      if find "$target" -mindepth 1 -print -quit >/dev/null 2>&1; then
        printf 'Clearing %s/%s\n' "$norm_id" "$tp"
        find "$target" -mindepth 1 -delete
      else
        printf 'Already empty: %s/%s\n' "$norm_id" "$tp"
      fi
      cleaned=$((cleaned + 1))
    else
      printf 'Missing folder: %s/%s (skipping)\n' "$norm_id" "$tp"
      skipped=$((skipped + 1))
    fi
  done
done

echo "Done. Processed $cleaned folders (skipped $skipped)."
