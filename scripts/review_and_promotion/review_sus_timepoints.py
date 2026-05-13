#!/usr/bin/env python3
"""Generate the SUS BL/3Y/5Y NIfTI selection log (best T1 per patient/timepoint)."""

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

if sys.version_info < (3, 8):
    sys.exit("Requires Python 3.8+. Load `uenv python-3.11.1` before running.")

DEFAULT_IDS = ["02", "05", "06", "08", "09", "10", "11", "13", "14", "18", "22", "26", "28", "30", "32", "34", "38", "41", "42", "43", "44", "45", "46", "49", "52"]

STAGING_ROOTS = {
    "BL": Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_BL_converted_to_nifti"),
    "3Y": Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_3Y_converted_to_nifti"),
    "5Y": Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_5Y_converted_to_nifti"),
}

OUTPUT_LOG = Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/scripts/sus_selection_log.csv")

PREFERRED = {
    "BL": {"shapes": [(256, 192, 256), (256, 224, 256)], "vox": [(0.9765625, 1.0, 0.9765625)]},
    "3Y": {"shapes": [(160, 256, 256)], "vox": [(1.0, 1.0, 1.0)]},
    "5Y": {"shapes": [(160, 256, 256)], "vox": [(1.0, 1.0, 1.0)]},
}


def approx_equal(lhs: Sequence[float], rhs: Sequence[float], tol: float = 1e-3) -> bool:
    return all(abs(a - b) <= tol for a, b in zip(lhs, rhs))


def is_nifti(path: Path) -> bool:
    return path.is_file() and path.name.endswith((".nii", ".nii.gz"))


def json_for_nifti(path: Path) -> Path:
    if path.name.endswith(".nii.gz"):
        return path.with_suffix("").with_suffix(".json")
    return path.with_suffix(".json")


def read_nifti_shape(path: Path) -> Tuple[Tuple[int, int, int], Tuple[float, float, float]]:
    import gzip
    import struct

    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rb") as f:  # type: ignore[arg-type]
        header = f.read(348)
    if len(header) < 348:
        raise ValueError(f"{path} header too short")
    dims = struct.unpack("<8h", header[40:56])
    pixdim = struct.unpack("<8f", header[76:108])
    shape = (int(dims[1]), int(dims[2]), int(dims[3]))
    vox = (float(pixdim[1]), float(pixdim[2]), float(pixdim[3]))
    return shape, vox


def normalize_text_list(value) -> List[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(v) for v in value]
    if isinstance(value, str):
        return [v.strip() for v in re.split(r"[,\\]", value) if v.strip()]
    return [str(value)]


@dataclass
class Candidate:
    path: Path
    json_path: Path
    dim: Tuple[int, int, int]
    vox: Tuple[float, float, float]
    series_description: str
    protocol_name: str
    sequence_name: str
    bids_guess: List[str]
    image_type: List[str]
    acquisition_type: str
    has_t1: bool
    has_mpr: bool
    is_original: bool
    fail_reasons: List[str] = field(default_factory=list)
    warn_reasons: List[str] = field(default_factory=list)
    preferred_shape: bool = False
    preferred_vox: bool = False


def evaluate_candidate(path: Path, preferred_shapes, preferred_voxels) -> Candidate:
    json_path = json_for_nifti(path)
    series_description = ""
    protocol_name = ""
    sequence_name = ""
    bids_guess: List[str] = []
    image_type: List[str] = []
    acquisition_type = ""
    fail_reasons: List[str] = []
    warn_reasons: List[str] = []

    text_sources = [path.name]

    if not json_path.exists():
        fail_reasons.append("missing JSON sidecar")
        meta = {}
    else:
        with json_path.open() as f:
            meta = json.load(f)
        series_description = meta.get("SeriesDescription", "")
        protocol_name = meta.get("ProtocolName", "")
        sequence_name = meta.get("SequenceName", "")
        bids_guess = normalize_text_list(meta.get("BidsGuess"))
        acquisition_type = str(meta.get("MRAcquisitionType", ""))
        image_type = normalize_text_list(meta.get("ImageType"))
        text_sources.extend([series_description, protocol_name, sequence_name, " ".join(bids_guess)])

    text = " ".join(ts for ts in text_sources if ts).lower()
    has_t1 = "t1" in text or any(token in path.name.lower() for token in ("t1w", "t1w_3d", "3d_t1", "t1_se"))
    has_mpr = any(token in text for token in ("mpr", "mprage", "tfl", "ffe", "tfe"))
    if not has_t1:
        fail_reasons.append("metadata/filename lacks t1 keyword")
    if not has_mpr:
        fail_reasons.append("metadata/filename lacks mpr/tfl/ffe/tfe signal")
    is_philips_tfe = any(token in text for token in ("3d_tfe", "t1w_3d_ffe", "st1w_3d", "mpr"))

    is_original = any("original" in item.lower() for item in image_type)
    if image_type and any("derived" in item.lower() for item in image_type):
        warn_reasons.append("ImageType reports DERIVED")
    if image_type and not is_original:
        warn_reasons.append("ImageType missing ORIGINAL flag")

    if acquisition_type and acquisition_type.upper() != "3D":
        warn_reasons.append(f"acquisition type {acquisition_type} != 3D")

    dim, vox = read_nifti_shape(path)
    if dim[2] <= 1 or dim[1] <= 1 or dim[0] <= 1:
        fail_reasons.append("single-slice volume")

    preferred_shape = dim in preferred_shapes
    preferred_vox = any(approx_equal(vox, target) for target in preferred_voxels)

    return Candidate(
        path=path,
        json_path=json_path,
        dim=dim,
        vox=vox,
        series_description=series_description,
        protocol_name=protocol_name,
        sequence_name=sequence_name,
        bids_guess=bids_guess,
        image_type=image_type,
        acquisition_type=acquisition_type,
        has_t1=has_t1,
        has_mpr=has_mpr,
        is_original=is_original,
        fail_reasons=fail_reasons,
        warn_reasons=warn_reasons,
        preferred_shape=preferred_shape,
        preferred_vox=preferred_vox,
    )


def candidate_score(candidate: Candidate) -> Tuple[int, int, int]:
    shape_score = 1 if candidate.preferred_shape else 0
    vox_score = 1 if candidate.preferred_vox else 0
    warning_penalty = len(candidate.warn_reasons)
    return (shape_score, vox_score, -warning_penalty)


def pick_candidate(candidates: List[Candidate]) -> Tuple[str, Optional[Candidate], str]:
    if not candidates:
        return "NO_CANDIDATES", None, "no NIfTI files found"
    clean = [c for c in candidates if not c.fail_reasons]
    if not clean:
        reasons = "; ".join(["/".join(sorted(set(c.fail_reasons))) for c in candidates])
        return "MANUAL_REVIEW", None, f"all candidates flagged: {reasons[:400]}"
    # Prefer Philips 3D TFE/MPR series if multiple candidates remain (especially for 5Y)
    clean.sort(key=candidate_score, reverse=True)
    best = clean[0]
    alt = [c for c in clean if c is not best and any(token in c.series_description.lower() for token in ("tfe", "ffe", "mpr"))]
    if alt:
        best = alt[0]
    clean = [best]
    top = clean[0]
    if len(clean) > 1 and candidate_score(clean[0]) == candidate_score(clean[1]):
        return "MANUAL_REVIEW", None, "multiple equally strong candidates"
    detail = " | ".join([
        f"shape={top.dim}",
        f"vox={','.join(f'{v:.3f}' for v in top.vox)}",
        f"series='{top.series_description}'",
    ])
    if top.warn_reasons:
        detail += " | warnings=" + "/".join(top.warn_reasons)
    return "SELECTED", top, detail


def iter_niftis(root: Path, stage_id: str) -> Iterable[Path]:
    patient_dir = root / stage_id
    if not patient_dir.exists():
        return []
    return sorted(p for p in patient_dir.rglob("*.nii*") if is_nifti(p))


def pad_stage(raw_id: str) -> str:
    return f"S{int(raw_id):02d}"


def pad_dest(raw_id: str) -> str:
    return f"S{int(raw_id):03d}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patients", nargs="*", help="Optional subset of SUS IDs (e.g., 02 05 06)")
    parser.add_argument("--timepoints", nargs="*", choices=["BL", "3Y", "5Y"], help="Optional subset of timepoints")
    args = parser.parse_args()

    ids = args.patients if args.patients else DEFAULT_IDS
    timepoints = args.timepoints if args.timepoints else ["BL", "3Y", "5Y"]

    rows: List[List[str]] = []
    summary: List[str] = []

    for raw_id in ids:
        stage_id = pad_stage(raw_id)
        dest_id = pad_dest(raw_id)
        for tp in timepoints:
            root = STAGING_ROOTS[tp]
            cfg = PREFERRED[tp]
            print(f"Reviewing {dest_id} {tp}...", flush=True)
            if not root.exists():
                status = "MISSING_STAGING_ROOT"
                detail = f"{root} not found"
                rows.append([tp, raw_id, stage_id, dest_id, status, detail, "", "", "", "", "", ""])
                summary.append(f"{dest_id} {tp}: {status} ({detail})")
                continue
            nifti_files = list(iter_niftis(root, stage_id))
            candidates = [evaluate_candidate(p, cfg["shapes"], cfg["vox"]) for p in nifti_files]
            status, selection, detail = pick_candidate(candidates)
            selected_path = selection.path if selection else None
            rows.append([
                tp,
                raw_id,
                stage_id,
                dest_id,
                status,
                detail,
                str(selected_path) if selected_path else "",
                "x".join(str(d) for d in (selection.dim if selection else ("", "", ""))),
                "x".join(f"{v:.3f}" for v in (selection.vox if selection else (0.0, 0.0, 0.0))),
                selection.series_description if selection else "",
                selection.protocol_name if selection else "",
                selection.sequence_name if selection else "",
            ])
            summary.append(f"{dest_id} {tp}: {status} → {selected_path if selected_path else detail}")

    OUTPUT_LOG.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_LOG.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "timepoint",
            "raw_id",
            "stage_id",
            "dest_id",
            "status",
            "detail",
            "selected_file",
            "shape",
            "voxel",
            "series_description",
            "protocol_name",
            "sequence_name",
        ])
        writer.writerows(rows)

    print(f"Wrote log to {OUTPUT_LOG} ({len(rows)} rows)")
    print("Summary:")
    for line in summary:
        print(" -", line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
