#!/usr/bin/env python3
"""List JSON metadata for every candidate NIfTI in a Bergen healthy staging folder."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, List, Tuple


def find_candidates(staging_root: Path, padded_id: str) -> List[Tuple[Path, Path]]:
    base = staging_root / padded_id
    if not base.exists():
        return []
    candidates: List[Tuple[Path, Path]] = []
    for path in sorted(base.rglob("*.nii*")):
        if not path.is_file():
            continue
        json_path = path.with_suffix("").with_suffix(".json") if path.name.endswith(".nii.gz") else path.with_suffix(".json")
        candidates.append((path, json_path))
    return candidates


def print_metadata(nifti_path: Path, json_path: Path) -> None:
    print(f"\n=== {nifti_path.name} ===")
    print(f"NIfTI: {nifti_path}")
    if json_path.exists():
        try:
            with json_path.open("r", encoding="utf-8") as f:
                meta = json.load(f)
        except Exception as exc:  # pragma: no cover
            print(f"[WARN] Failed to parse {json_path}: {exc}")
            return
        for key in (
            "SeriesDescription",
            "ProtocolName",
            "SequenceName",
            "ImageType",
            "MRAcquisitionType",
            "EchoTime",
            "RepetitionTime",
            "FlipAngle",
        ):
            value = meta.get(key)
            if value is not None:
                print(f"{key}: {value}")
    else:
        print(f"[WARN] Missing JSON sidecar: {json_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staging-root", required=True, help="Path to the staging folder (converted NIfTIs)")
    parser.add_argument("patient", help="Patient ID (e.g., 39 or BK39)")
    args = parser.parse_args()

    staging_root = Path(args.staging_root).expanduser().resolve()
    if not staging_root.exists():
        parser.error(f"Staging root not found: {staging_root}")

    patient = args.patient.upper()
    if patient.startswith("BK"):
        padded_id = patient
    else:
        padded_id = f"BK{int(patient):02d}"

    candidates = find_candidates(staging_root, padded_id)
    if not candidates:
        print(f"No candidates found under {staging_root}/{padded_id}")
        return 0

    print(f"Found {len(candidates)} candidate(s) for {padded_id} in {staging_root}")
    for nifti_path, json_path in candidates:
        print_metadata(nifti_path, json_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
