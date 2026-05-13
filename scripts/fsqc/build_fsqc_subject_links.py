#!/usr/bin/env python3
"""Create fsqc-friendly subject trees by symlinking FastSurfer outputs.

The script flattens the FastSurfer longitudinal directory layout into
cohort-specific folders (patients vs healthy controls) so that fsqc can iterate
through a single `subjects_dir`. Each symlink is named `<SITE>_<PATIENT>_<TP>`
and points to the corresponding FastSurfer cross-sectional folder (BL/3Y/5Y).
It also emits a text file listing all generated subject IDs for convenient use
with `run_fsqc --subjects-file`.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Iterable, List

DEFAULT_SOURCE = \
    "/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/outputs/fastsurfer_longitudinal"
DEFAULT_DEST = \
    "/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/outputs/fsqc_subjects"

COHORT_SITES = {
    "patients": ["Bergen", "SUS", "Forde"],
    "controls": ["Bergen_healthy", "SUS_healthy", "Forde_healthy"],
}

TIMEPOINTS = ["BL", "3Y", "5Y"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-root",
        default=DEFAULT_SOURCE,
        help="FastSurfer longitudinal root (default: %(default)s)",
    )
    parser.add_argument(
        "--dest-root",
        default=DEFAULT_DEST,
        help="Destination root for symlinks and subject lists (default: %(default)s)",
    )
    parser.add_argument(
        "--cohorts",
        choices=["patients", "controls", "all"],
        default=["patients", "controls"],
        nargs="+",
        help="Which cohort(s) to build. 'all' processes patients and controls.",
    )
    parser.add_argument(
        "--timepoints",
        choices=TIMEPOINTS,
        default=TIMEPOINTS,
        nargs="+",
        help="Timepoints to include (default: %(default)s)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing symlinks/lists if they already exist.",
    )
    return parser.parse_args()


def resolve_cohorts(raw: Iterable[str]) -> List[str]:
    cohorts: List[str] = []
    for entry in raw:
        if entry == "all":
            for name in ("patients", "controls"):
                if name not in cohorts:
                    cohorts.append(name)
        elif entry not in cohorts:
            cohorts.append(entry)
    return cohorts


def ensure_clean_path(path: Path, overwrite: bool) -> None:
    if path.is_symlink() or path.exists():
        if not overwrite:
            return
        if path.is_dir() and not path.is_symlink():
            # Symlink destinations should never be directories, but we guard anyway.
            raise RuntimeError(f"Refusing to remove existing directory: {path}")
        path.unlink()


def main() -> None:
    args = parse_args()
    source_root = Path(args.source_root).expanduser().resolve()
    dest_root = Path(args.dest_root).expanduser().resolve()
    timepoints = args.timepoints

    if not source_root.is_dir():
        raise SystemExit(f"Source root does not exist: {source_root}")

    dest_root.mkdir(parents=True, exist_ok=True)

    for cohort in resolve_cohorts(args.cohorts):
        allowed_sites = COHORT_SITES.get(cohort)
        if not allowed_sites:
            raise SystemExit(f"Unsupported cohort label: {cohort}")

        cohort_dir = dest_root / cohort
        cohort_dir.mkdir(parents=True, exist_ok=True)
        subjects: List[str] = []

        for site in sorted(allowed_sites):
            site_dir = source_root / site
            if not site_dir.is_dir():
                print(f"[warn] Missing site folder for {site} in {source_root}")
                continue
            for patient_dir in sorted(site_dir.iterdir()):
                if not patient_dir.is_dir():
                    continue
                patient = patient_dir.name
                for tp in timepoints:
                    fs_dir = patient_dir / tp
                    if not fs_dir.is_dir():
                        continue
                    subject_id = f"{site}_{patient}_{tp}"
                    link_path = cohort_dir / subject_id
                    ensure_clean_path(link_path, overwrite=args.overwrite)
                    if link_path.exists():
                        # Already pointing somewhere; skip unless overwrite was specified.
                        print(f"[skip] {link_path} already exists")
                        continue
                    os.symlink(fs_dir, link_path)
                    subjects.append(subject_id)

        subjects.sort()
        list_path = cohort_dir / f"{cohort}_subjects.txt"
        if list_path.exists() and not args.overwrite:
            print(f"[info] Subject list exists (use --overwrite to refresh): {list_path}")
        else:
            with list_path.open("w", encoding="utf-8") as fh:
                fh.write("\n".join(subjects) + ("\n" if subjects else ""))
        print(f"[done] {cohort}: linked {len(subjects)} subjects (list: {list_path})")


if __name__ == "__main__":
    main()
