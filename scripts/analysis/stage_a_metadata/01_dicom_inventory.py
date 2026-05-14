#!/usr/bin/env python3
"""Walk the ParkWest DICOM tree and write a per-series metadata CSV.

Reads one representative DICOM per series across all timepoints, sites,
and groups, and writes dicom_inventory_all.csv (one row per series)
for 02_scanner_metadata.py to consume.

Inconsistent on-disk layouts handled:
  Baseline/
    Pasienter/
      Stavanger_p_anonym/S05/1/1/*.dcm              # SUS: numeric series
      Bergen_Forde_p_anonym/B 03/t1_mpr_n/*.dcm     # combined Bergen/Forde
      Bergen_Forde_p_anonym/F01/t1_mpr_n/*.dcm
    Kontroller/
      Stavanger_k_anonym/SK01/...
      Bergen_Forde_k_anonym/BK 17/...

  3Y/
    Pasienter/
      Stavanger_p_anonym/S05/1/1/*.dcm
      Bergen_p_anonym/B005/1/1/*.dcm                # now 3-digit
      Forde_p_anonym/F01/1/1/*.dcm
    Kontroller/
      Stavanger_k_anonym/ / Bergen_k_anonym/ / Forde_k_anonym/

  5Y/
    Pasienter/
      PV_patients_SUS_5Y/S05/DICOM/<hash>/...
      PV_patients_Bergen_5Y/B10/DICOM/<numeric>/...
      PV_patients_Frde_5Y/F04/DICOM/<hash>/...      # Forde encoding varies
    Kontroller/
      PV_controls_SUS_5Y/ / PV_controls_Bergen_5Y/ / PV_controls_Frde_5Y/

Subject IDs are run through pipeline_lib.ids.normalize_subject_id.

Usage:
    python scripts/stage_a_metadata/01_dicom_inventory.py \\
        --output-dir outputs/stage_a_metadata/

This is slow (~30-60 min) but only runs once; the output is reused by
every subsequent 02_scanner_metadata.py run.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import warnings
from pathlib import Path
from typing import Optional

import pandas as pd

from pipeline_lib.ids import normalize_subject_id


# Silence pydicom's "value length exceeds VR maximum" warning. Some
# SUS 5Y DICOMs have a SoftwareVersions string at 65 chars (one over
# the VR LO limit of 64). Cosmetic, not a parsing problem.
warnings.filterwarnings("ignore", category=UserWarning, module="pydicom.valuerep")


# Read-only NFS location of the ParkWest imaging archive.
PARKWEST_ROOT = Path("/nfs/br1_prosjekt/ParkWest/ImageData")


def resolve_path(parent: Path, *name_candidates: str) -> Optional[Path]:
    """Find a subdirectory of `parent` matching any candidate name.

    Tries an exact match first, then a normalized match (lowercased,
    non-ASCII stripped) so different on-disk encodings of 'ø' end up
    at the same folder.
    """
    # Try exact match first.
    for n in name_candidates:
        try:
            p = parent / n
            if p.is_dir():
                return p
        except (UnicodeEncodeError, OSError):
            pass

    if not parent.is_dir():
        return None

    def normalize(s: str) -> str:
        low = s.lower()
        return re.sub(r"[^a-z0-9_\-]", "", low)

    norm_candidates = [normalize(n) for n in name_candidates]

    try:
        entries = list(parent.iterdir())
    except OSError:
        return None

    for entry in entries:
        if not entry.is_dir():
            continue
        norm_entry = normalize(entry.name)
        for nc in norm_candidates:
            if nc and (nc == norm_entry or nc in norm_entry or norm_entry in nc):
                return entry

    # Couldn't match; dump what's actually there.
    print(f"  resolve_path failed for candidates {name_candidates}")
    print(f"  entries under {parent}:")
    for entry in entries[:20]:
        try:
            raw_bytes = entry.name.encode("utf-8", errors="backslashreplace")
        except Exception:
            raw_bytes = b"<encoding error>"
        print(f"    name={entry.name!r}  bytes={raw_bytes!r}")
    return None


def find_forde_folder(parent: Path, prefix: str) -> Optional[Path]:
    """Find a subfolder starting with `prefix` and containing 'rde'.

    Matches Forde / Førde / Frde variants regardless of how 'ø' is
    encoded on disk.
    """
    if not parent.is_dir():
        print(f"  WARN: parent {parent} not a directory")
        return None
    prefix_low = prefix.lower()
    try:
        entries = list(parent.iterdir())
    except OSError as e:
        print(f"  WARN: cannot list {parent}: {e}")
        return None
    matches = []
    for entry in entries:
        if not entry.is_dir():
            continue
        name_low = entry.name.lower()
        if name_low.startswith(prefix_low) and "rde" in name_low:
            matches.append(entry)
    if len(matches) == 1:
        print(f"  Found Forde folder: {matches[0]}")
        return matches[0]
    if len(matches) > 1:
        print(f"  WARN: multiple matches for prefix {prefix!r}, using first: {matches}")
        return matches[0]
    print(f"  WARN: no Forde folder matching prefix {prefix!r} found under {parent}")
    print(f"  entries in {parent}:")
    for entry in entries[:30]:
        print(f"    {entry.name!r}  (bytes: {os.fsencode(entry.name)!r})")
    return None


def _build_roots():
    """Build the {(Site, Group, Timepoint): root_dir} map.

    Forde folders are resolved with fuzzy matching because the on-disk
    encoding of 'ø' varies.
    """
    BL = PARKWEST_ROOT / "ParkVest_Baseline"
    Y3 = PARKWEST_ROOT / "ParkVest_3Y"
    Y5 = PARKWEST_ROOT / "ParkVest_5Y"
    r = {
        # Baseline
        ("SUS",    "Pasienter",  "BL"): BL / "Pasienter" / "Stavanger_p_anonym",
        ("BerFor", "Pasienter",  "BL"): BL / "Pasienter" / "Bergen_Forde_p_anonym",
        ("SUS",    "Kontroller", "BL"): BL / "Kontroller" / "Stavanger_k_anonym",
        ("BerFor", "Kontroller", "BL"): BL / "Kontroller" / "Bergen_Forde_k_anonym",
        # 3Y
        ("SUS",    "Pasienter",  "3Y"): Y3 / "Pasienter" / "Stavanger_p_anonym",
        ("Bergen", "Pasienter",  "3Y"): Y3 / "Pasienter" / "Bergen_p_anonym",
        ("Forde",  "Pasienter",  "3Y"): Y3 / "Pasienter" / "Forde_p_anonym",
        ("SUS",    "Kontroller", "3Y"): Y3 / "Kontroller" / "Stavanger_k_anonym",
        ("Bergen", "Kontroller", "3Y"): Y3 / "Kontroller" / "Bergen_k_anonym",
        ("Forde",  "Kontroller", "3Y"): Y3 / "Kontroller" / "Forde_k_anonym",
        # 5Y. Forde paths use varying encodings of 'ø', so resolve them
        # via find_forde_folder.
        ("SUS",    "Pasienter",  "5Y"): Y5 / "Pasienter" / "PV_patients_SUS_5Y",
        ("Bergen", "Pasienter",  "5Y"): Y5 / "Pasienter" / "PV_patients_Bergen_5Y",
        ("Forde",  "Pasienter",  "5Y"): find_forde_folder(
            Y5 / "Pasienter", "PV_patients_F"
        ),
        ("SUS",    "Kontroller", "5Y"): Y5 / "Kontroller" / "PV_controls_SUS_5Y",
        ("Bergen", "Kontroller", "5Y"): Y5 / "Kontroller" / "PV_controls_Bergen_5Y",
        ("Forde",  "Kontroller", "5Y"): find_forde_folder(
            Y5 / "Kontroller", "PV_controls_F"
        ),
    }
    # Drop None entries so the scanner doesn't try to walk them.
    return {k: v for k, v in r.items() if v is not None}


# Built lazily in main() to avoid import-time disk I/O (which would
# break testing and --help on machines without NFS access).
ROOTS: dict = {}

# DICOM tags we want from every series.
TAGS = [
    # Scanner / protocol
    "Manufacturer", "ManufacturerModelName", "DeviceSerialNumber",
    "StationName", "SoftwareVersions", "MagneticFieldStrength",
    "RepetitionTime", "EchoTime", "InversionTime", "FlipAngle",
    "PixelSpacing", "SliceThickness", "AcquisitionMatrix",
    "SeriesDescription", "ProtocolName", "SeriesNumber",
    "InstitutionName", "StudyDate", "AcquisitionDate",
    "ReceiveCoilName", "ScanningSequence", "SequenceVariant",
    "Modality",
    # Demographics (may be blank when fully anonymized).
    "PatientSex", "PatientAge", "PatientBirthDate", "PatientID",
    "PatientWeight", "PatientSize",
]


def scalar(v):
    """Flatten a DICOM value to a plain Python scalar."""
    if v is None:
        return None
    if hasattr(v, "__iter__") and not isinstance(v, str):
        return str(list(v))
    return str(v)


def infer_site(site_key: str, subject_id_raw: str) -> str:
    """Map a ROOTS key's Site token to the actual site name.

    'BerFor' roots combine Bergen and Forde; split by the first letter
    of the subject ID (B vs F) to recover the true site.
    """
    if site_key != "BerFor":
        return site_key
    s = subject_id_raw.strip().replace(" ", "")
    if s.startswith("B"):
        return "Bergen"
    elif s.startswith("F"):
        return "Forde"
    return "Unknown"


def walk_subject(subject_dir: Path):
    """Yield one representative DICOM per SeriesInstanceUID under subject_dir.

    Some subjects have folder-per-series (1/, 2/, ...) and some have
    everything dumped flat in one folder (e.g. Forde 5Y). For each
    folder, if it has more than FLAT_THRESHOLD files we treat it as a
    flat dump and probe several files to find each series UID;
    otherwise we read just the first file.

    Yields (dataset, series_folder, n_files_in_series).
    """
    import pydicom

    folders: "dict[Path, list[Path]]" = {}
    for f in subject_dir.rglob("*"):
        if not f.is_file() or f.name in {"VERSION", "DICOMDIR"}:
            continue
        folders.setdefault(f.parent, []).append(f)

    FLAT_THRESHOLD = 50           # files per folder above which we suspect a flat dump
    MAX_PROBES_PER_FOLDER = 40    # cap on the number of DICOMs we sniff

    seen_series_uids: set = set()

    for folder, files in folders.items():
        files = sorted(files)
        if len(files) <= FLAT_THRESHOLD:
            # Folder per series: read the first file.
            for f in files:
                try:
                    d = pydicom.dcmread(str(f), stop_before_pixels=True, force=True)
                except Exception:
                    continue
                modality = getattr(d, "Modality", None)
                if modality is not None and modality != "MR":
                    continue
                series_uid = getattr(d, "SeriesInstanceUID", None) or f"folder:{folder}"
                if series_uid in seen_series_uids:
                    break
                seen_series_uids.add(series_uid)
                yield d, folder, len(files)
                break
        else:
            # Flat dump: probe enough files to cover every series.
            step = max(1, len(files) // MAX_PROBES_PER_FOLDER)
            series_first_seen: "dict[str, tuple]" = {}
            series_count: "dict[str, int]" = {}
            for f in files[::step]:
                try:
                    d = pydicom.dcmread(str(f), stop_before_pixels=True, force=True)
                except Exception:
                    continue
                modality = getattr(d, "Modality", None)
                if modality is not None and modality != "MR":
                    continue
                uid = getattr(d, "SeriesInstanceUID", None)
                if uid is None:
                    continue
                if uid not in series_first_seen:
                    series_first_seen[uid] = (d, f)
                series_count[uid] = series_count.get(uid, 0) + 1
            # Per-series file count is extrapolated from the sample.
            for uid, (d, _) in series_first_seen.items():
                if uid in seen_series_uids:
                    continue
                seen_series_uids.add(uid)
                estimated_n = series_count[uid] * step
                yield d, folder, estimated_n


def scan_one_root(site_key: str, group: str, tp: str, root: Path):
    rows = []
    if not root.is_dir():
        print(f"  SKIP: {root} does not exist")
        return rows
    print(f"  Scanning {site_key}/{group}/{tp}: {root}")
    for subj_dir in sorted(root.iterdir()):
        if not subj_dir.is_dir():
            continue
        raw_id = subj_dir.name
        subj_id = normalize_subject_id(raw_id)
        actual_site = infer_site(site_key, raw_id)
        n_series_in_subj = 0
        for d, folder, n_files in walk_subject(subj_dir):
            row = {
                "Site":         actual_site,
                "Group":        group,
                "Subject":      subj_id,
                "Subject_raw":  raw_id,
                "Timepoint":    tp,
                "SeriesFolder": str(folder.relative_to(subj_dir)),
                "NFiles":       n_files,
            }
            for t in TAGS:
                row[t] = scalar(getattr(d, t, None))
            rows.append(row)
            n_series_in_subj += 1
        if n_series_in_subj == 0:
            print(f"    WARN: no MR series for {raw_id}")
    return rows


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--output-dir", type=Path,
        default=Path("outputs/stage_a_metadata"),
        help="Directory for dicom_inventory_all.csv "
             "(default: outputs/stage_a_metadata).",
    )
    args = p.parse_args()

    # Late import so --help works without pydicom installed.
    try:
        import pydicom  # noqa: F401
    except ImportError:
        print("ERROR: pydicom not installed. Run: pip install pydicom")
        return 1

    global ROOTS
    ROOTS = _build_roots()

    all_rows = []
    print(f"Scanning {len(ROOTS)} roots under {PARKWEST_ROOT}")
    for (site_key, group, tp), root in ROOTS.items():
        all_rows.extend(scan_one_root(site_key, group, tp, root))

    if not all_rows:
        print("No DICOMs found. Check ROOTS paths.")
        return 1

    df = pd.DataFrame(all_rows)

    # Some DICOM tags (notably SoftwareVersions) stringify as lists with
    # commas and single-quotes, so we quote every non-numeric field.
    # Also scrub embedded null bytes seen in some corrupted strings.
    def scrub(v):
        if isinstance(v, str):
            return v.replace("\x00", "").replace("\r", " ").replace("\n", " ")
        return v

    df = df.apply(lambda col: col.map(scrub) if col.dtype == object else col)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    inv_path = args.output_dir / "dicom_inventory_all.csv"
    df.to_csv(inv_path, index=False, quoting=csv.QUOTE_NONNUMERIC)
    print(f"\nWrote {inv_path} ({len(df)} rows)")

    print("\nSubjects per Site x Group x Timepoint:")
    print(
        df.groupby(["Site", "Group", "Timepoint"])["Subject"].nunique()
          .unstack("Timepoint", fill_value=0)
    )

    print("\nDemographics availability:")
    for col in ["PatientSex", "PatientAge", "PatientBirthDate", "PatientWeight"]:
        filled = df[col].notna().sum()
        total  = len(df)
        print(f"  {col}: {filled}/{total} non-null ({100*filled/total:.1f}%)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
