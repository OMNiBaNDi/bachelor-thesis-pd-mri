#!/usr/bin/env python3
"""Copy SELECTED SUS NIfTI+JSON files into data/SUS/Sxxx/{BL,3Y,5Y}."""

import argparse
import csv
import shutil
from pathlib import Path

LOG_PATH = Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/scripts/sus_selection_log.csv")
DEST_ROOT = Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS")
TIMEPOINT_DIRS = {"BL": "BL", "3Y": "3Y", "5Y": "5Y"}


def json_sidecar(path: Path) -> Path:
  if path.name.endswith(".nii.gz"):
    return path.with_suffix("").with_suffix(".json")
  return path.with_suffix(".json")


def copy_pair(src: Path, dest_dir: Path, force: bool, dry_run: bool) -> None:
  dest_dir.mkdir(parents=True, exist_ok=True)
  dest_nifti = dest_dir / src.name
  dest_json = dest_dir / json_sidecar(src).name

  if dest_nifti.exists() and not force:
    print(f"SKIP (exists): {dest_nifti}")
  else:
    print(f"COPY: {src} -> {dest_nifti}")
    if not dry_run:
      shutil.copy2(src, dest_nifti)

  src_json = json_sidecar(src)
  if src_json.exists():
    if dest_json.exists() and not force:
      print(f"SKIP (exists): {dest_json}")
    else:
      print(f"COPY: {src_json} -> {dest_json}")
      if not dry_run:
        shutil.copy2(src_json, dest_json)
  else:
    print(f"WARN: JSON sidecar missing for {src}")


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("--log-path", type=Path, default=LOG_PATH)
  parser.add_argument("--dest-root", type=Path, default=DEST_ROOT)
  parser.add_argument("--timepoint", choices=["BL", "3Y", "5Y"], help="Limit to a single timepoint")
  parser.add_argument("--dry-run", action="store_true")
  parser.add_argument("--force", action="store_true")
  args = parser.parse_args()

  if not args.log_path.exists():
    parser.error(f"Log file not found: {args.log_path}")

  copied = 0
  skipped = 0

  with args.log_path.open() as f:
    reader = csv.DictReader(f)
    for row in reader:
      tp = row["timepoint"].strip()
      if args.timepoint and tp != args.timepoint:
        continue
      if row["status"].strip().upper() != "SELECTED":
        continue
      dest_id = row["dest_id"].strip()
      selected_file = row["selected_file"].strip()
      if not selected_file:
        print(f"WARN: no selected file for {dest_id} {tp}")
        skipped += 1
        continue
      src = Path(selected_file)
      if not src.exists():
        print(f"WARN: source missing: {src}")
        skipped += 1
        continue
      dest_dir = args.dest_root / dest_id / TIMEPOINT_DIRS[tp]
      copy_pair(src, dest_dir, args.force, args.dry_run)
      copied += 1

  print(f"Done. Processed {copied} entries (skipped {skipped}).{' Dry run.' if args.dry_run else ''}")


if __name__ == "__main__":
  main()
