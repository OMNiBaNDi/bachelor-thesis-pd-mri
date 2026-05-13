#!/bin/bash
# Copy the curated Bergen baseline DICOM folders into the thesis-local Bergen_BL_DICOM staging area.
# Run this on the UIS cluster where /nfs/br1_prosjekt/... is mounted.
set -euo pipefail

SRC=/nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_Baseline/Pasienter/Bergen_Forde_p_anonym
DEST=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/BergenConversion/Bergen_BL_DICOM
LIST_FILE=""
MODE="default"  # default=hardcoded list, file=list file, all=every B directory

DEFAULT_PATIENT_IDS=(
  05 14 15 17 20 21 24 25 26 28 29 33 34 35 36 38
  39 40 43 53 54 55 56 57 58 63 71 72 74 83 84 90
  92 94 95 101
)

usage() {
  cat <<'EOF'
Usage: copy_bergen_bl.sh [--patients-file FILE | --all]

Default (no flags): copies the curated 36 Bergen baseline patients required for
the longitudinal pipeline.
  --patients-file FILE  Copy only the IDs listed in FILE (one per line).
  --all                 Copy every directory matching "B *" under the source.
EOF
}

trim() {
  local s="$1"
  s="${s#${s%%[!$'\t\r\n ']*}}"
  s="${s%${s##*[!$'\t\r\n ']}}"
  printf '%s' "$s"
}

resolve_patient_dir() {
  local token="$(trim "$1")"
  [[ -n "$token" ]] || return 1
  if [[ "$token" != B* ]]; then
    token="B $token"
  fi
  local candidates=()
  candidates+=("$token")
  local nospace="${token// /}"
  [[ -n "$nospace" && "$nospace" != "$token" ]] && candidates+=("$nospace")
  if [[ "$nospace" =~ ^B([0-9]+)$ ]]; then
    local num="${BASH_REMATCH[1]}"
    num=$((10#$num))
    local fmt
    printf -v fmt 'B %d' "$num"; candidates+=("$fmt")
    printf -v fmt 'B %02d' "$num"; candidates+=("$fmt")
    printf -v fmt 'B %03d' "$num"; candidates+=("$fmt")
  fi
  local -A seen=()
  for cand in "${candidates[@]}"; do
    [[ -n "$cand" && -z "${seen[$cand]:-}" ]] || continue
    seen[$cand]=1
    if [[ -d "$SRC/$cand" ]]; then
      printf '%s' "$cand"
      return 0
    fi
  done
  return 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -f|--patients-file)
      LIST_FILE="$2"
      MODE="file"
      shift 2
      ;;
    --all)
      MODE="all"
      shift
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

mkdir -p "$DEST"
cd "$SRC"

PATIENT_TOKENS=()
case "$MODE" in
  default)
    for id in "${DEFAULT_PATIENT_IDS[@]}"; do
      PATIENT_TOKENS+=("B $id")
    done
    ;;
  file)
    if [[ ! -f "$LIST_FILE" ]]; then
      echo "Patient list file not found: $LIST_FILE" >&2
      exit 1
    fi
    mapfile -t RAW_PATIENTS < "$LIST_FILE"
    for entry in "${RAW_PATIENTS[@]}"; do
      entry="$(trim "$entry")"
      [[ -z "$entry" ]] && continue
      PATIENT_TOKENS+=("$entry")
    done
    ;;
  all)
    shopt -s nullglob
    PATIENT_TOKENS=(B\ *)
    ;;
  *)
    echo "Unknown mode: $MODE" >&2
    exit 1
    ;;
esac

count=0
skipped=0
for token in "${PATIENT_TOKENS[@]}"; do
  local_dir="$token"
  if [[ "$MODE" != "all" ]]; then
    if ! resolved="$(resolve_patient_dir "$token")"; then
      printf 'Skipping %s (no matching folder under %s)\n' "$token" "$SRC"
      skipped=$((skipped + 1))
      continue
    fi
    local_dir="$resolved"
  fi
  printf 'Copying %-8s → %s\n' "$local_dir" "$DEST/$local_dir"
  rsync -a --info=progress2 --human-readable \
    "$SRC/$local_dir/" \
    "$DEST/$local_dir/"
  count=$((count + 1))
done

printf 'Done. Copied %d patient folders' "$count"
if [[ "$MODE" != "all" ]]; then
  printf ' (skipped %d missing entries).' "$skipped"
fi
printf '\n'
