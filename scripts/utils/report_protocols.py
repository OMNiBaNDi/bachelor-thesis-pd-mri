#!/usr/bin/env python3
"""Scan thesis data folders and report the scanner protocol for each BL/3Y/5Y NIfTI.

If no paths are provided, the script automatically scans the standard thesis data
sites (Bergen, Bergen_healthy, Førde, Førde_healthy, SUS, SUS_healthy).
"""

import argparse
import json
import logging
import re
from pathlib import Path
from typing import Iterable, List, Optional

try:
    import nibabel as nib  # type: ignore
except ModuleNotFoundError:  # pragma: no cover
    nib = None

LOG = logging.getLogger(__name__)
PROTOCOL_KEYS = ("ProtocolName", "SeriesDescription", "SequenceName")
DEFAULT_SITE_DIRS = [
    "Bergen",
    "Bergen_healthy",
    "Forde",
    "Forde_healthy",
    "SUS",
    "SUS_healthy",
]


def default_paths() -> List[str]:
    script_dir = Path(__file__).resolve().parent
    data_root = script_dir.parent / "data"
    return [str((data_root / site).resolve()) for site in DEFAULT_SITE_DIRS]


def find_nifti_files(root: Path) -> Iterable[Path]:
    """Yield the first NIfTI inside every BL/3Y/5Y folder under the root."""
    if not root.exists():
        LOG.warning("%s does not exist", root)
        return []

    for folder in sorted(root.rglob("*")):
        if not folder.is_dir():
            continue
        if folder.name.upper() not in {"BL", "3Y", "5Y"}:
            continue
        candidates = sorted(folder.glob("*.nii")) + sorted(folder.glob("*.nii.gz"))
        if not candidates:
            LOG.warning("No NIfTI found in %s", folder)
            continue
        yield candidates[0]


def json_sidecar(nifti_path: Path) -> Optional[Path]:
    if nifti_path.suffix == ".gz" and nifti_path.stem.endswith(".nii"):
        json_path = nifti_path.with_suffix("").with_suffix(".json")
    else:
        json_path = nifti_path.with_suffix(".json")
    return json_path if json_path.exists() else None


def protocol_from_json(json_path: Path) -> Optional[str]:
    try:
        with json_path.open("r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception as exc:  # pragma: no cover
        LOG.warning("Failed to parse %s: %s", json_path, exc)
        return None
    return extract_protocol_from_dict(data)


def extract_protocol_from_dict(data: dict) -> Optional[str]:
    for key in PROTOCOL_KEYS:
        value = data.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return None


def protocol_from_nifti(nifti_path: Path) -> Optional[str]:
    """Use nibabel to read JSON extensions or descrip text."""
    if nib is None:
        LOG.debug("nibabel not installed; cannot inspect %s", nifti_path)
        return None
    try:
        img = nib.load(str(nifti_path))
    except Exception as exc:  # pragma: no cover
        LOG.warning("Failed to load %s: %s", nifti_path, exc)
        return None

    for ext in img.header.extensions:
        try:
            payload = ext.get_content()
        except Exception:
            continue
        if isinstance(payload, bytes):
            payload = payload.decode("utf-8", errors="ignore")
        if not isinstance(payload, str):
            continue
        text = payload.strip()
        if not text.startswith("{"):
            continue
        try:
            data = json.loads(text)
        except json.JSONDecodeError:
            continue
        protocol = extract_protocol_from_dict(data)
        if protocol:
            return protocol

    descrip = img.header.get("descrip")
    if descrip is not None:
        if hasattr(descrip, "tobytes"):
            descrip = descrip.tobytes()
        if isinstance(descrip, bytes):
            descrip = descrip.decode("utf-8", errors="ignore")
        if isinstance(descrip, str):
            descrip = descrip.strip()
            if descrip:
                match = re.search(r"ProtocolName=([^;\\n]+)", descrip)
                if match:
                    return match.group(1).strip()
                return descrip
    return None


def report_protocols(paths: Iterable[str]) -> None:
    for raw in paths:
        root = Path(raw).expanduser()
        print(f"\n# {root}")
        if not root.exists():
            print("(path not found)")
            continue
        for nifti in find_nifti_files(root):
            protocol = None
            json_path = json_sidecar(nifti)
            if json_path:
                protocol = protocol_from_json(json_path)
            if protocol is None:
                protocol = protocol_from_nifti(nifti)
            rel = nifti.relative_to(root)
            print(f"{rel}\t{protocol or '<protocol not found>'}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "paths",
        nargs="*",
        help="Paths to scan (default: the standard thesis data cohorts)",
    )
    parser.add_argument("--log", default="WARNING", help="Logging level")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    logging.basicConfig(level=args.log.upper())
    paths = args.paths if args.paths else default_paths()
    report_protocols(paths)


if __name__ == "__main__":
    main()
