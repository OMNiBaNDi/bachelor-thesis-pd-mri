#!/bin/bash
set -e

if [ $# -ne 1 ]; then
  echo "Usage: $0 /full/path/to/site_folder"
  exit 1
fi

INPUT_FOLDER="$1"

if [ ! -d "$INPUT_FOLDER" ]; then
  echo "ERROR: folder does not exist: $INPUT_FOLDER"
  exit 1
fi

SITE=$(basename "$INPUT_FOLDER")
SCRIPT=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/scripts/longitudinal_patient.slurm

echo "Running longitudinal submission on folder:"
echo "  $INPUT_FOLDER"
echo "Detected site:"
echo "  $SITE"
echo

submitted=0
skipped=0

for patient_dir in "$INPUT_FOLDER"/*; do
  [ -d "$patient_dir" ] || continue

  PATIENT=$(basename "$patient_dir")

  BL_FILE=$(find "$patient_dir/BL" -name "*.nii*" 2>/dev/null | head -n 1)
  Y3_FILE=$(find "$patient_dir/3Y" -name "*.nii*" 2>/dev/null | head -n 1)
  Y5_FILE=$(find "$patient_dir/5Y" -name "*.nii*" 2>/dev/null | head -n 1)

  if [ -n "$BL_FILE" ] && [ -n "$Y3_FILE" ] && [ -n "$Y5_FILE" ]; then
    echo "Submitting $SITE/$PATIENT"
    sbatch \
      --job-name=${SITE}_${PATIENT}_long \
      --export=ALL,SITE=$SITE,PATIENT=$PATIENT \
      "$SCRIPT"
    submitted=$((submitted + 1))
  else
    echo "Skipping $SITE/$PATIENT (missing BL/3Y/5Y nifti)"
    skipped=$((skipped + 1))
  fi
done

echo
echo "Done."
echo "Submitted: $submitted"
echo "Skipped:   $skipped"