#!/usr/bin/env python3
"""Inspect converted Bergen BL/3Y NIfTI files and recommend the best T1 per patient."""

import sys

if sys.version_info < (3, 8):
    sys.exit("This script requires Python 3.8+. Please load `uenv python-3.11` before running it.")

import argparse
import csv
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

DEFAULT_PATIENT_IDS = [
    "05", "14", "15", "17", "20", "21", "24", "25", "26", "28", "29", "33",
    "34", "35", "36", "38", "39", "40", "43", "53", "54", "55", "56", "57",
    "58", "63", "71", "72", "74", "83", "84", "90", "92", "94", "95", "101",
]

DEFAULT_BL_STAGING = Path(
    "/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/BergenConversion/Bergen_BL_converted_to_nifti"
)
DEFAULT_3Y_STAGING = Path(
    "/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/BergenConversion/Bergen_3Y_converted_to_nifti"
)

TIMEPOINT_SHAPES: Dict[str, Dict[str, object]] = {
    "BL": {
        "preferred_shapes": [(256, 192, 256), (256, 224, 256)],
        "preferred_voxels": [(0.9765625, 1.0, 0.9765625)],
    },
    "3Y": {
        "preferred_shapes": [(160, 256, 256)],
        "preferred_voxels": [(1.0, 1.0, 1.0)],
    },
}

DEFAULT_LOG = Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/outputs/bergen_bl3y_selection_log.csv")


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
    fail_reasons: List[str] = field(default_factory=list)
    warn_reasons: List[str] = field(default_factory=list)
    preferred_shape: bool = False
    preferred_vox: bool = False


def evaluate_candidate(path: Path, preferred_shapes: Sequence[Tuple[int, int, int]], preferred_voxels: Sequence[Tuple[float, float, float]]) -> Candidate:
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
        text_sources.extend([
            series_description,
            protocol_name,
            sequence_name,
            " ".join(bids_guess),
        ])

    text = " ".join(ts for ts in text_sources if ts).lower()
    has_t1 = "t1" in text
    has_mpr = any(token in text for token in ("mpr", "mprage", "tfl"))
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
        is_roi=is_roi,
        is_single_slice_aux=is_single_slice_aux,
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
    clean.sort(key=candidate_score, reverse=True)
    top = clean[0]
    if len(clean) > 1 and candidate_score(clean[0]) == candidate_score(clean[1]):
        return "MANUAL_REVIEW", None, "multiple equally strong candidates"
    detail = " | ".join(
        [
            f"shape={top.dim}",
            f"vox={','.join(f'{v:.3f}' for v in top.vox)}",
            f"series='{top.series_description}'",
        ]
    )
    if top.warn_reasons:
        detail += " | warnings=" + "/".join(top.warn_reasons)
    return "SELECTED", top, detail


def pad_id(raw: str) -> str:
    digits = int(raw)
    return f"B{digits:03d}"


def iter_nifti_files(root: Path, patient: str) -> Iterable[Path]:
    patient_dir = root / patient
    if not patient_dir.exists():
        return []
    return sorted(p for p in patient_dir.rglob("*.nii*") if is_nifti(p))


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["bl", "3y", "both"], default="both")
    parser.add_argument("--patients", nargs="*", help="Optional list of patient IDs (e.g., 05 or B005)")
    parser.add_argument("--log-path", type=Path, default=DEFAULT_LOG)
    parser.add_argument("--bl-staging", type=Path, default=DEFAULT_BL_STAGING, help="Path to the BL staging directory")
    parser.add_argument("--y3-staging", type=Path, default=DEFAULT_3Y_STAGING, help="Path to the 3Y staging directory")
    args = parser.parse_args(argv)

    if args.patients:
        patient_ids = args.patients
    else:
        patient_ids = DEFAULT_PATIENT_IDS

    requested_tps = ["BL", "3Y"] if args.mode == "both" else ["BL" if args.mode == "bl" else "3Y"]

    timepoint_configs: Dict[str, Dict[str, object]] = {}
    for tp in requested_tps:
        if tp == "BL":
            staging = args.bl_staging
        else:
            staging = args.y3_staging
        cfg = TIMEPOINT_SHAPES[tp]
        timepoint_configs[tp] = {
            "staging": staging,
            "preferred_shapes": list(cfg["preferred_shapes"]),  # type: ignore[index]
            "preferred_voxels": list(cfg["preferred_voxels"]),  # type: ignore[index]
        }

    rows: List[List[str]] = []
    summary: List[str] = []

    for raw_id in patient_ids:
        normalized = raw_id.upper().replace("B", "").strip()
        padded = pad_id(normalized)
        for tp in requested_tps:
            cfg = timepoint_configs[tp]
            staging_root: Path = cfg["staging"]  # type: ignore[assignment]
            preferred_shapes = cfg["preferred_shapes"]  # type: ignore[assignment]
            preferred_voxels = cfg["preferred_voxels"]  # type: ignore[assignment]

            patient_dir = staging_root / padded
            if not patient_dir.exists():
                status = "MISSING_STAGING"
                detail = f"{patient_dir} not found"
                rows.append([tp, raw_id, padded, status, detail, "", "", "", "", "", ""])
                summary.append(f"{padded} {tp}: {status} ({detail})")
                continue

            nifti_files = list(iter_nifti_files(staging_root, padded))
            candidates = [evaluate_candidate(p, preferred_shapes, preferred_voxels) for p in nifti_files]
            status, selection, detail = pick_candidate(candidates)

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
    sys.exit(main())
