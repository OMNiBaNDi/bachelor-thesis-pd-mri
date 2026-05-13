#!/usr/bin/env python3
"""List SUS healthy 3Y/5Y converted files with JSON metadata and write to a CSV."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import List

DEFAULT_ROOTS = {
    "3Y": Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_healthy_3Y_converted_to_nifti"),
    "5Y": Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_healthy_5Y_converted_to_nifti"),
}
DEFAULT_IDS = [
    "05", "08", "10", "11", "13", "15", "16", "23", "26",
    "27", "28", "31", "34", "35", "37", "41", "46", "47",
]

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT = SCRIPT_DIR / "sus_healthy_converted_series.csv"


def find_files(root: Path, patient: str) -> List[Path]:
    patient_dir = root / patient
    if not patient_dir.exists():
        return []
    return sorted(p for p in patient_dir.glob("*.nii*") if p.is_file())


def json_for_nifti(path: Path) -> Path:
    if path.name.endswith(".nii.gz"):
        return path.with_suffix("").with_suffix(".json")
    return path.with_suffix(".json")


def load_metadata(json_path: Path) -> dict:
    try:
        with json_path.open("r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timepoints", nargs="+", choices=["3Y", "5Y"], default=["3Y", "5Y"])
    parser.add_argument("--patients", nargs="+", help="Optional subset of IDs (e.g., 05 08)")
    parser.add_argument("--root-3y", type=Path, default=DEFAULT_ROOTS["3Y"])
    parser.add_argument("--root-5y", type=Path, default=DEFAULT_ROOTS["5Y"])
    args = parser.parse_args()

    patient_ids = args.patients if args.patients else DEFAULT_IDS
    root_map = {"3Y": args.root_3y, "5Y": args.root_5y}

    rows: List[List[str]] = []

    for tp in args.timepoints:
        root = root_map[tp]
        for raw_id in patient_ids:
            padded = f"SK{int(raw_id):02d}"
            nifti_files = find_files(root, padded)
            if not nifti_files:
                rows.append([tp, raw_id, padded, "", "", "", "", "", ""])
                continue
            for nifti in nifti_files:
                meta = load_metadata(json_for_nifti(nifti))
                rows.append([
                    tp,
                    raw_id,
                    padded,
                    nifti.name,
                    str(nifti),
                    meta.get("SeriesDescription", ""),
                    meta.get("ProtocolName", ""),
                    meta.get("SequenceName", ""),
                    meta.get("ImageType", ""),
                ])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "timepoint",
            "raw_id",
            "padded_id",
            "filename",
            "path",
            "SeriesDescription",
            "ProtocolName",
            "SequenceName",
            "ImageType",
        ])
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
