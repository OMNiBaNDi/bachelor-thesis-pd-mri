#!/bin/bash
set -euo pipefail

PW=/nfs/br1_prosjekt/ParkWest
THESIS=$PW/user/2026vae/AmundEspen/pd_thesis
DCM2NIIX=$THESIS/tools/dcm2niix/dcm2niix
INPUT_ROOT=$THESIS/data/SUS_healthy_conversion/SUS_healthy_BL_DICOM
OUTPUT_ROOT=$THESIS/data/SUS_healthy_conversion/SUS_healthy_BL_converted_to_nifti

PATIENTS=(
  SK05
  SK08
  SK10
  SK11
  SK13
  SK15
  SK16
  SK23
  SK26
  SK27
  SK28
  SK31
  SK34
  SK35
  SK37
  SK41
  SK46
  SK47
)

if [ ! -x "$DCM2NIIX" ]; then
  echo "ERROR: dcm2niix not found or not executable: $DCM2NIIX"
  exit 1
fi

mkdir -p "$OUTPUT_ROOT"

echo "SUS healthy BL batch conversion"
echo "Input root:  $INPUT_ROOT"
echo "Output root: $OUTPUT_ROOT"
echo "dcm2niix:    $DCM2NIIX"
echo

converted=0
skipped=0
failed=0

for patient in "${PATIENTS[@]}"; do
  input_dir="$INPUT_ROOT/$patient"
  output_dir="$OUTPUT_ROOT/$patient"

  if [ ! -d "$input_dir" ]; then
    echo "SKIP $patient -- missing input folder: $input_dir"
    skipped=$((skipped + 1))
    continue
  fi

  mkdir -p "$output_dir"

  echo "=================================================="
  echo "Converting $patient"
  echo "  from: $input_dir"
  echo "  to:   $output_dir"

  if "$DCM2NIIX" -z y -b y -o "$output_dir" "$input_dir"; then
    echo "OK   $patient"
    converted=$((converted + 1))
  else
    echo "FAIL $patient"
    failed=$((failed + 1))
  fi
  echo

done

echo "Batch conversion done."
echo "Converted: $converted"
echo "Skipped:   $skipped"
echo "Failed:    $failed"
