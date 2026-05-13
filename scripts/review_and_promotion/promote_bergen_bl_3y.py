#!/usr/bin/env python3
"""Copy the selected Bergen BL/3Y NIfTIs (and JSON) into data/Bergen/Bxxx/{BL,3Y}."""

import argparse
import csv
import shutil
from pathlib import Path

DEFAULT_LOG = Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/scripts/bergen_bl3y_selection_log.csv")
DEFAULT_DEST = Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Bergen")
TIMEPOINT_DIRS = {"BL": "BL", "3Y": "3Y"}


def json_sidecar(path: Path) -> Path:
    if path.name.endswith(".nii.gz"):
        return path.with_suffix("").with_suffix(".json")
    return path.with_suffix(".json")


def copy_with_json(src: Path, dest_dir: Path, dry_run: bool, force: bool) -> None:
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log-path", type=Path, default=DEFAULT_LOG)
    parser.add_argument("--dest-root", type=Path, default=DEFAULT_DEST)
    parser.add_argument("--timepoint", choices=["BL", "3Y"], help="Limit to one timepoint")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true", help="Overwrite existing files")
    args = parser.parse_args()

    if not args.log_path.exists():
        parser.error(f"Log file not found: {args.log_path}")

    copied = 0
    skipped = 0

    with args.log_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            timepoint = row["timepoint"].strip()
            if args.timepoint and timepoint != args.timepoint:
                continue
            if row["status"].strip().upper() != "SELECTED":
                continue
            padded_id = row["padded_id"].strip()
            selected_file = row["selected_file"].strip()
            if not selected_file:
                print(f"WARN: no selected file recorded for {padded_id} {timepoint}")
                skipped += 1
                continue
            src = Path(selected_file)
            if not src.exists():
                print(f"WARN: source missing {src}")
                skipped += 1
                continue
            tp_dir = TIMEPOINT_DIRS.get(timepoint)
            if not tp_dir:
                print(f"WARN: unknown timepoint {timepoint} for {padded_id}")
                skipped += 1
                continue
            dest_dir = args.dest_root / padded_id / tp_dir
            copy_with_json(src, dest_dir, args.dry_run, args.force)
            copied += 1

    print(f"Done. Processed {copied} entries (skipped {skipped}).{' Dry run only.' if args.dry_run else ''}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
