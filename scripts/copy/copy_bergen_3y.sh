#!/bin/bash
# Copy the curated Bergen 3Y DICOM folders into the thesis-local Bergen_3Y_DICOM staging area.
# Run this on the UIS cluster where /nfs/br1_prosjekt/... is mounted.
set -euo pipefail

SRC=/nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_3Y/Pasienter/Bergen_p_anonym
DEST=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/BergenConversion/Bergen_3Y_DICOM
LIST_FILE=""

DEFAULT_PATIENT_IDS=(
  05 14 15 17 20 21 24 25 26 28 29 33 34 35 36 38
  39 40 43 53 54 55 56 57 58 63 71 72 74 83 84 90
  92 94 95 101
)

usage() {
  cat <<'EOF'
Usage: copy_bergen_3y.sh [--patients-file FILE]

Default (no flags): copies the curated 36 Bergen 3Y patients required for the longitudinal pipeline.
  --patients-file FILE  Copy only the IDs listed in FILE (one per line).
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
  if [[ -d "$SRC/$token" ]]; then
    printf '%s' "$token"
    return 0
  fi
  return 1
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

mkdir -p "$DEST"
cd "$SRC"

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

count=0
skipped=0
for id in "${PATIENT_IDS[@]}"; do
  if ! folder="$(normalize_id "$id")"; then
    printf 'Skipping %s (unable to normalize to B### form)\n' "$id"
    skipped=$((skipped + 1))
    continue
  fi
  if [[ ! -d "$SRC/$folder" ]]; then
    printf 'Skipping %s (folder %s not found under %s)\n' "$id" "$folder" "$SRC"
    skipped=$((skipped + 1))
    continue
  fi
  printf 'Copying %-4s → %s/%s\n' "$folder" "$DEST" "$folder"
  rsync -a --info=progress2 --human-readable \
    "$SRC/$folder/" \
    "$DEST/$folder/"
  count=$((count + 1))
done

printf 'Done. Copied %d patient folders (skipped %d).\n' "$count" "$skipped"
