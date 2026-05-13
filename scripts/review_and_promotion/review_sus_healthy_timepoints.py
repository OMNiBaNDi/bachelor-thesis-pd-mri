#!/usr/bin/env python3
"""Review SUS healthy-control 3Y/5Y conversions (SKxx) and log selected T1 volumes."""

from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

DEFAULT_IDS = [
    "05", "08", "10", "11", "13", "15", "16", "23", "26",
    "27", "28", "31", "34", "35", "37", "41", "46", "47",
]

DEFAULT_3Y_STAGING = Path(
    "/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_healthy_3Y_converted_to_nifti"
)
DEFAULT_5Y_STAGING = Path(
    "/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_healthy_5Y_converted_to_nifti"
)

TIMEPOINT_PREFS: Dict[str, Dict[str, Sequence[Tuple[int, int, int]] | Sequence[Tuple[float, float, float]]]] = {
    "3Y": {
        "preferred_shapes": [(160, 256, 256)],
        "preferred_voxels": [(1.0, 1.0, 1.0)],
    },
    "5Y": {
        "preferred_shapes": [(160, 256, 256)],
        "preferred_voxels": [(1.0, 1.0, 1.0)],
    },
}

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_LOG = SCRIPT_DIR / "sus_healthy_selection_log.csv"
PREFERRED_KEYWORDS = (
    "st1w", "s_t1w", "3d_tfe", "t1w", "t1/tfe", "t1w/tfe", "mpr", "mpr_ns",
)


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
    is_roi: bool
    is_single_slice_aux: bool
    is_preferred_series: bool
    fail_reasons: List[str] = field(default_factory=list)
    warn_reasons: List[str] = field(default_factory=list)


def evaluate_candidate(
    path: Path,
    preferred_shapes: Sequence[Tuple[int, int, int]],
    preferred_voxels: Sequence[Tuple[float, float, float]],
) -> Candidate:
    name_lower = path.name.lower()
    json_path = json_for_nifti(path)
    series_description = ""
    protocol_name = ""
    sequence_name = ""
    bids_guess: List[str] = []
    image_type: List[str] = []
    acquisition_type = ""
    fail_reasons: List[str] = []
    warn_reasons: List[str] = []

    text_sources: List[str] = [path.name]

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
    has_t1 = any(token in text for token in ("t1", "t1w", "st1", "s_t1"))
    has_mpr = any(keyword in text for keyword in PREFERRED_KEYWORDS)
    if not has_t1:
        fail_reasons.append("metadata/filename lacks t1 keyword")
    if not has_mpr:
        fail_reasons.append("metadata/filename lacks mpr/tfl signal")

    is_original = any("original" in item.lower() for item in image_type)
    if image_type and any("derived" in item.lower() for item in image_type):
        warn_reasons.append("ImageType reports DERIVED")
    if image_type and not is_original:
        warn_reasons.append("ImageType missing ORIGINAL flag")

    if acquisition_type and acquisition_type.upper() != "3D":
        warn_reasons.append(f"acquisition type {acquisition_type} != 3D")

    dim, vox = read_nifti_shape(path)
    if dim[2] <= 1:
        fail_reasons.append("single-slice volume")

    is_roi = "roi" in name_lower or "roi" in series_description.lower()
    if is_roi:
        fail_reasons.append("ROI/derived series")

    is_single_slice_aux = bool(re.search(r"_a(?:\.nii)?(?:\.gz)?$", name_lower))
    if is_single_slice_aux:
        fail_reasons.append("auxiliary *_a file")

    is_preferred_series = any(keyword in text for keyword in PREFERRED_KEYWORDS)

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
        is_roi=is_roi,
        is_single_slice_aux=is_single_slice_aux,
        is_preferred_series=is_preferred_series,
        fail_reasons=fail_reasons,
        warn_reasons=warn_reasons,
    )


def pad_id(raw: str) -> str:
    return f"SK{int(raw):02d}"


def iter_nifti_files(root: Path, patient: str) -> Iterable[Path]:
    patient_dir = root / patient
    if not patient_dir.exists():
        return []
    return sorted(patient_dir.glob("*.nii*"))


def pick_candidate(tp: str, candidates: List[Candidate]) -> Tuple[str, Optional[Candidate], str]:
    if not candidates:
        return "NO_CANDIDATES", None, "no NIfTI files found"
    clean = [c for c in candidates if not c.fail_reasons]
    if not clean:
        reasons = "; ".join(["/".join(sorted(set(c.fail_reasons))) for c in candidates])
        return "MANUAL_REVIEW", None, f"all candidates flagged: {reasons[:400]}"

    prefs = TIMEPOINT_PREFS[tp]
    pref_shapes: Sequence[Tuple[int, int, int]] = prefs["preferred_shapes"]  # type: ignore[assignment]
    pref_voxels: Sequence[Tuple[float, float, float]] = prefs["preferred_voxels"]  # type: ignore[assignment]

    def score(c: Candidate) -> Tuple[int, int, int, int]:
        preferred_series_score = 1 if c.is_preferred_series else 0
        shape_score = 1 if c.dim in pref_shapes else 0
        vox_score = 1 if any(approx_equal(c.vox, target) for target in pref_voxels) else 0
        warning_penalty = len(c.warn_reasons)
        return (preferred_series_score, shape_score, vox_score, -warning_penalty)

    clean.sort(key=score, reverse=True)
    top = clean[0]
    if len(clean) > 1 and score(clean[0]) == score(clean[1]):
        return "MANUAL_REVIEW", None, "multiple equally strong candidates"

    detail = " | ".join([
        f"shape={top.dim}",
        f"vox={','.join(f'{v:.3f}' for v in top.vox)}",
        f"series='{top.series_description}'",
    ])
    if top.warn_reasons:
        detail += " | warnings=" + "/".join(top.warn_reasons)
    return "SELECTED", top, detail


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--timepoints",
        nargs="+",
        choices=["3Y", "5Y"],
        default=["3Y", "5Y"],
        help="Timepoints to inspect (default: 3Y 5Y)",
    )
    parser.add_argument("--patients", nargs="+", help="Optional subset of IDs (e.g., 05 08)")
    parser.add_argument("--log-path", type=Path, default=DEFAULT_LOG)
    parser.add_argument("--staging-3y", type=Path, default=DEFAULT_3Y_STAGING)
    parser.add_argument("--staging-5y", type=Path, default=DEFAULT_5Y_STAGING)
    args = parser.parse_args()

    patient_ids = args.patients if args.patients else DEFAULT_IDS
    staging_map = {"3Y": args.staging_3y, "5Y": args.staging_5y}

    rows: List[List[str]] = []
    summary: List[str] = []

    for raw_id in patient_ids:
        padded = pad_id(raw_id)
        for tp in args.timepoints:
            staging_root = staging_map[tp]
            if not staging_root.exists():
                status = "MISSING_STAGING"
                detail = f"{staging_root} not found"
                rows.append([tp, raw_id, padded, status, detail, "", "", "", "", "", ""])
                summary.append(f"{padded} {tp}: {status} ({detail})")
                continue

            nifti_files = list(iter_nifti_files(staging_root, padded))
            prefs = TIMEPOINT_PREFS[tp]
            candidates = [evaluate_candidate(p, prefs["preferred_shapes"], prefs["preferred_voxels"]) for p in nifti_files]  # type: ignore[arg-type]
            status, selection, detail = pick_candidate(tp, candidates)

            selected_path = selection.path if selection else None
            rows.append([
                tp,
                raw_id,
                padded,
                status,
                detail,
                str(selected_path) if selected_path else "",
                "x".join(str(d) for d in (selection.dim if selection else ("", "", ""))),
                "x".join(f"{v:.3f}" for v in (selection.vox if selection else (0.0, 0.0, 0.0))),
                selection.series_description if selection else "",
                selection.protocol_name if selection else "",
                selection.sequence_name if selection else "",
            ])
            summary.append(f"{padded} {tp}: {status} → {selected_path if selected_path else detail}")

    args.log_path.parent.mkdir(parents=True, exist_ok=True)
    with args.log_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "timepoint",
            "raw_id",
            "padded_id",
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

    print(f"Wrote log to {args.log_path} ({len(rows)} rows)")
    print("\nSummary:")
    for line in summary:
        print(" -", line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
