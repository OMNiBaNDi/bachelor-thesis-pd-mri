#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: run_fsqc_on_cohort.sh COHORT output_suffix" >&2
  echo "  COHORT = patients | controls" >&2
  echo "  output_suffix = label appended to outputs/fsqc/<cohort>_<suffix>" >&2
  exit 1
fi

COHORT="$1"
SUFFIX="$2"
BASE=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis
SUBJ_DIR="$BASE/outputs/fsqc_subjects/$COHORT"
LIST_FILE="$SUBJ_DIR/${COHORT}_subjects.txt"
OUT_DIR="$BASE/outputs/fsqc/${COHORT}_${SUFFIX}"
LOG_FILE="$OUT_DIR/fsqc_commands.log"

if [[ ! -d "$SUBJ_DIR" ]]; then
  echo "Missing subject directory: $SUBJ_DIR" >&2
  exit 1
fi
if [[ ! -f "$LIST_FILE" ]]; then
  echo "Missing subject list: $LIST_FILE" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

scripts/run_fsqc_chunk.sh \
  --subjects-dir "$SUBJ_DIR" \
  --subjects-file "$LIST_FILE" \
  --output-dir "$OUT_DIR" \
  --log-file "$LOG_FILE"
