#!/bin/bash
set -euo pipefail
shopt -s nullglob

THESIS=/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis
SRC_ROOT=$THESIS/data/SUS_healthy_conversion/SUS_healthy_BL_converted_to_nifti
DST_ROOT=$THESIS/data/SUS_healthy
BACKUP_TAG=before_BL_3D_FFE_replace_2026-04-19

PATIENTS=(
  SK05 SK08 SK10 SK11 SK13 SK15 SK16 SK23 SK26
  SK27 SK28 SK31 SK34 SK35 SK37 SK41 SK46 SK47
)

for p in "${PATIENTS[@]}"; do
  src_dir="$SRC_ROOT/$p"
  dst_dir="$DST_ROOT/$p/BL"
  backup_dir="$dst_dir/$BACKUP_TAG"

  src_nii=( "$src_dir"/"${p}"_PARKVEST_T1W_3D_FFE_*.nii.gz "$src_dir"/"${p}"_PARKVEST_T1W_3D_FFE_*.nii )
  src_json=( "$src_dir"/"${p}"_PARKVEST_T1W_3D_FFE_*.json )

  if [ "${#src_nii[@]}" -ne 1 ]; then
    echo "ERROR $p -- expected exactly 1 3D FFE NIfTI, found ${#src_nii[@]}"
    printf '  %s\n' "${src_nii[@]}"
    exit 1
  fi

  if [ ! -d "$dst_dir" ]; then
    echo "ERROR $p -- missing destination BL folder: $dst_dir"
    exit 1
  fi

  mkdir -p "$backup_dir"

  existing=( "$dst_dir"/*.nii "$dst_dir"/*.nii.gz "$dst_dir"/*.json )
  if [ "${#existing[@]}" -gt 0 ]; then
    mv "${existing[@]}" "$backup_dir"/
  fi

  cp -v "${src_nii[0]}" "$dst_dir"/
  if [ "${#src_json[@]}" -eq 1 ]; then
    cp -v "${src_json[0]}" "$dst_dir"/
  fi

  echo "DONE $p"
done

echo
echo "Replacement complete for ${#PATIENTS[@]} subjects."
echo "Backups stored under each BL folder as: $BACKUP_TAG"
