#!/usr/bin/env python3
"""Select and copy the correct Bergen 5Y T1w volume into the curated dataset.

The script inspects all converted NIfTI files inside the staging directory
(`Bergen_5Y_converted_to_nifti/<rawID>`), applies the selection rule discussed in
`research/dicom_to_nifti.md`, copies the chosen NIfTI+JSON pair into
`data/Bergen/<paddedID>/5Y`, and emits a CSV log describing every decision.

Key heuristics (configurable via CLI):
- prefer filenames/metadata containing both `t1` and `mpr`
- skip ROI / derived / one-slice ("*_a") files
- require `ImageType` to include ORIGINAL when present
- prefer the canonical `(160,256,256)` shape at 1.0 mm isotropic
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

PREFERRED_SHAPE = (160, 256, 256)
PREFERRED_VOXEL = (1.0, 1.0, 1.0)
DEFAULT_STAGING = "/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Bergen_5Y_converted_to_nifti"
DEFAULT_FINAL = "/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/Bergen"
DEFAULT_LOG = "/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/outputs/bergen_5y_selection_log.csv"


def parse_tuple(text: str, cast=float) -> Tuple:
    parts = [p.strip() for p in text.replace("x", ",").split(",") if p.strip()]
    return tuple(cast(p) for p in parts)


def is_nifti(path: Path) -> bool:
    return path.is_file() and (path.suffix in {".nii", ".gz"} and path.name.endswith((".nii", ".nii.gz")))


def json_for_nifti(path: Path) -> Path:
    if path.name.endswith(".nii.gz"):
        return path.with_suffix("").with_suffix(".json")
    return path.with_suffix(".json")


def approx_equal(lhs: Sequence[float], rhs: Sequence[float], tol: float = 1e-3) -> bool:
    return all(abs(a - b) <= tol for a, b in zip(lhs, rhs))


def read_nifti_shape(path: Path) -> Tuple[Tuple[int, int, int], Tuple[float, float, float]]:
    import gzip
    import struct

    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rb") as f:  # type: ignore[arg-type]
        header = f.read(348)
    if len(header) < 348:
        raise ValueError(f"{path} header too short")
    dims = struct.unpack("<8h", header[40:40 + 16])
    pixdim = struct.unpack("<8f", header[76:76 + 32])
    shape = (int(dims[1]), int(dims[2]), int(dims[3]))
    vox = (float(pixdim[1]), float(pixdim[2]), float(pixdim[3]))
    return shape, vox


def normalize_list(value) -> List[str]:
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

    @property
    def preferred_shape(self) -> bool:
        return self.dim == preferred_shape

    @property
    def preferred_vox(self) -> bool:
        return approx_equal(self.vox, preferred_voxel)


# globals set at runtime for dataclass properties
preferred_shape = PREFERRED_SHAPE
preferred_voxel = PREFERRED_VOXEL


def evaluate_candidate(path: Path, preferred_shape_cfg, preferred_voxel_cfg) -> Candidate:
    global preferred_shape, preferred_voxel
    preferred_shape = preferred_shape_cfg
    preferred_voxel = preferred_voxel_cfg

    name_lower = path.name.lower()
    json_path = json_for_nifti(path)
    series_description = ""
    protocol_name = ""
    sequence_name = ""
    bids_guess: List[str] = []
    image_type: List[str] = []
    acquisition_type = ""
    has_t1 = False
    has_mpr = False
    is_original = False
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
        bids_guess = normalize_list(meta.get("BidsGuess"))
        acquisition_type = str(meta.get("MRAcquisitionType", ""))
        image_type = normalize_list(meta.get("ImageType"))
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
    )


def candidate_score(candidate: Candidate) -> Tuple[int, int, int]:
    shape_score = 1 if candidate.dim == preferred_shape else 0
    vox_score = 1 if approx_equal(candidate.vox, preferred_voxel) else 0
    warning_penalty = len(candidate.warn_reasons)
    return (shape_score, vox_score, -warning_penalty)


def pick_candidate(candidates: List[Candidate]) -> Tuple[str, Optional[Candidate], str]:
    clean = [c for c in candidates if not c.fail_reasons]
    if not candidates:
        return "NO_CANDIDATES", None, "no NIfTI files found"
    if not clean:
        reasons = "; ".join(["/".join(sorted(set(c.fail_reasons))) for c in candidates])
        return "MANUAL_REVIEW", None, f"all candidates flagged: {reasons[:400]}"
    clean.sort(key=candidate_score, reverse=True)
    best = clean[0]
    if len(clean) > 1 and candidate_score(clean[0]) == candidate_score(clean[1]):
        return "MANUAL_REVIEW", None, "multiple equally strong candidates"
    detail_parts = [
        f"shape={best.dim}",
        f"vox={','.join(f'{v:.3f}' for v in best.vox)}",
        f"series='{best.series_description}'",
    ]
    if best.warn_reasons:
        detail_parts.append("warnings=" + "/".join(best.warn_reasons))
    return "SELECTED", best, " | ".join(detail_parts)


def pad_patient_id(raw_id: str) -> str:
    m = re.match(r"([A-Za-z]+)(\d+)", raw_id)
    if not m:
        return raw_id
    prefix, digits = m.groups()
    return f"{prefix.upper()}{int(digits):03d}"


def iter_patients(staging_root: Path) -> Iterable[Path]:
    for entry in sorted(staging_root.iterdir()):
        if entry.is_dir():
            yield entry


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_log(log_path: Path, rows: List[List[str]]) -> None:
    ensure_dir(log_path.parent)
    with log_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "raw_id",
            "padded_id",
            "status",
            "detail",
            "selected_file",
            "dest_file",
            "shape",
            "voxel",
            "series_description",
            "protocol_name",
            "sequence_name",
        ])
        writer.writerows(rows)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staging-root", default=DEFAULT_STAGING, type=Path)
    parser.add_argument("--final-root", default=DEFAULT_FINAL, type=Path)
    parser.add_argument("--log-path", default=DEFAULT_LOG, type=Path)
    parser.add_argument("--preferred-shape", default="160,256,256")
    parser.add_argument("--preferred-voxel", default="1.0,1.0,1.0")
    parser.add_argument("--patients", nargs="*", help="limit to these raw or padded IDs")
    parser.add_argument("--dry-run", action="store_true", help="compute selections without copying")
    parser.add_argument("--force", action="store_true", help="overwrite existing files in the final dataset")
    args = parser.parse_args(argv)

    preferred_shape_cfg = tuple(int(x) for x in parse_tuple(args.preferred_shape, int))
    preferred_voxel_cfg = tuple(float(x) for x in parse_tuple(args.preferred_voxel, float))

    staging_root: Path = args.staging_root
    final_root: Path = args.final_root
    log_rows: List[List[str]] = []

    if not staging_root.exists():
        parser.error(f"Staging root {staging_root} does not exist")

    patient_filter = {p.upper() for p in args.patients} if args.patients else None

    for patient_dir in iter_patients(staging_root):
        raw_id = patient_dir.name
        padded_id = pad_patient_id(raw_id)
        if patient_filter and raw_id.upper() not in patient_filter and padded_id.upper() not in patient_filter:
            continue

        nifti_files = [p for p in patient_dir.rglob("*.nii*") if is_nifti(p)]
        candidates = [evaluate_candidate(p, preferred_shape_cfg, preferred_voxel_cfg) for p in nifti_files]
        status, candidate, detail = pick_candidate(candidates)

        selected_file = str(candidate.path) if candidate else ""
        dest_file = ""

        if status == "SELECTED" and candidate:
            dest_dir = final_root / padded_id / "5Y"
            ensure_dir(dest_dir)
            dest_path = dest_dir / candidate.path.name
            dest_json = dest_dir / candidate.json_path.name
            if dest_path.exists() and not args.force:
                status = "ALREADY_PRESENT"
                detail = f"dest file exists: {dest_path}"
            else:
                dest_file = str(dest_path)
                if not args.dry_run:
                    shutil.copy2(candidate.path, dest_path)
                    if candidate.json_path.exists():
                        shutil.copy2(candidate.json_path, dest_json)
            if args.dry_run:
                status = "WOULD_COPY"
                dest_file = str(dest_path)

        log_rows.append([
            raw_id,
            padded_id,
            status,
            detail,
            selected_file,
            dest_file,
            "x".join(str(d) for d in (candidate.dim if candidate else ("", "", ""))),
            "x".join(f"{v:.3f}" for v in (candidate.vox if candidate else (0.0, 0.0, 0.0))),
            candidate.series_description if candidate else "",
            candidate.protocol_name if candidate else "",
            candidate.sequence_name if candidate else "",
        ])

    write_log(args.log_path, log_rows)
    print(f"Wrote log to {args.log_path} ({len(log_rows)} patients)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
