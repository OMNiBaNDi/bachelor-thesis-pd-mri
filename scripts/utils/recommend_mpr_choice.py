#!/usr/bin/env python3
"""Recommend which t1_mpr JSON/NIfTI pair to keep when multiple candidates exist.

Usage:
  python3 scripts/recommend_mpr_choice.py --patients B70 B92

The script scans each patient's staging folder under
`/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Bergen_5Y_converted_to_nifti/<rawID>`
for `.json` files, extracts key metadata, and prints the top-ranked recommendation
along with all candidate details so you can double-check the decision.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence, Tuple

STAGING_ROOT = Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Bergen_5Y_converted_to_nifti")


@dataclass
class Candidate:
    json_path: Path
    nii_path: Path
    series_description: str
    protocol_name: str
    sequence_name: str
    series_number: int
    acquisition_time: str
    series_time: str
    image_type: List[str]
    is_t1_mpr: bool
    is_original: bool

    @property
    def score(self) -> Tuple[int, int, int, str]:
        return (
            1 if self.is_t1_mpr else 0,
            1 if self.is_original else 0,
            self.series_number,
            self.acquisition_time or self.series_time,
        )


def detect_t1_mpr(name: str, series_desc: str, protocol: str, sequence: str) -> bool:
    text = " ".join(filter(None, [name, series_desc, protocol, sequence])).lower()
    return ("t1" in text) and any(token in text for token in ("mpr", "mprage", "tfl"))


def load_candidates(patient: str) -> List[Candidate]:
    patient_dir = STAGING_ROOT / patient
    candidates: List[Candidate] = []
    if not patient_dir.exists():
        return candidates
    for json_path in sorted(patient_dir.glob("*.json")):
        nii_path = json_path.with_suffix("")
        if not nii_path.exists() and not nii_path.with_suffix(".nii.gz").exists():
            continue
        if not nii_path.exists():
            nii_path = nii_path.with_suffix(".nii.gz")
        with json_path.open() as f:
            meta = json.load(f)
        series_desc = meta.get("SeriesDescription", "")
        protocol = meta.get("ProtocolName", "")
        sequence = meta.get("SequenceName", "")
        series_number = int(meta.get("SeriesNumber", 0) or 0)
        acquisition_time = str(meta.get("AcquisitionTime", ""))
        series_time = str(meta.get("SeriesTime", ""))
        image_type = meta.get("ImageType", [])
        if isinstance(image_type, str):
            image_type = [image_type]
        is_original = any("ORIGINAL" in item.upper() for item in image_type)
        name = json_path.stem
        is_t1_mpr = detect_t1_mpr(name, series_desc, protocol, sequence)
        candidates.append(
            Candidate(
                json_path=json_path,
                nii_path=nii_path,
                series_description=series_desc,
                protocol_name=protocol,
                sequence_name=sequence,
                series_number=series_number,
                acquisition_time=acquisition_time,
                series_time=series_time,
                image_type=image_type,
                is_t1_mpr=is_t1_mpr,
                is_original=is_original,
            )
        )
    return [c for c in candidates if c.is_t1_mpr]


def recommend_for_patient(patient: str) -> None:
    candidates = load_candidates(patient)
    print(f"\nPatient {patient}")
    if not candidates:
        print("  No t1_mpr-like candidates found.")
        return
    candidates.sort(key=lambda c: c.score, reverse=True)
    best = candidates[0]
    print("  Recommendation:")
    print(f"    {best.nii_path.name} (SeriesNumber={best.series_number}, AcquisitionTime={best.acquisition_time or best.series_time})")
    print("  Candidates:")
    for cand in candidates:
        warn = []
        if not cand.is_original:
            warn.append("not ORIGINAL")
        if warn:
            warn_text = " | warnings: " + ", ".join(warn)
        else:
            warn_text = ""
        print(
            f"    - {cand.nii_path.name} (SeriesNumber={cand.series_number}, "
            f"AcqTime={cand.acquisition_time or cand.series_time}, ImageType={';'.join(cand.image_type)}){warn_text}"
        )


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patients", nargs="*", default=["B70", "B92"])
    args = parser.parse_args(argv)

    for patient in args.patients:
        recommend_for_patient(patient)


if __name__ == "__main__":
    main()
