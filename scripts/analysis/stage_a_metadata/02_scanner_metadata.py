#!/usr/bin/env python3
"""Filter the DICOM inventory to 3D T1 series and assign batch IDs.

One row per (Site x Group x Subject x Timepoint), restricted to the 3D
T1 series, with these columns added:

    ScannerID   fine-grained protocol identifier (~10 values).
    BatchID     collapsed batch label (6 values) passed to longCombat.
                Functionally-equivalent protocols are collapsed so
                batch sizes stay large.
    in_cohort   True if this subject is in cohort.csv (all three
                timepoints complete). Everything downstream filters
                on this flag.

Usage:
    python scripts/stage_a_metadata/02_scanner_metadata.py \\
        --inventory outputs/stage_a_metadata/dicom_inventory_all.csv \\
        --cohort    outputs/stage_a_metadata/cohort.csv \\
        --output    outputs/stage_a_metadata/scanner_metadata.csv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


# --- 3D T1 structural series filter ---

T1_3D_KEYWORDS = [
    # Philips (SUS)
    r"3d/ffe", r"3d_ffe", r"ffe\b", r"tfe\b", r"3d/tfe",
    # Siemens (Bergen, Forde): MPRAGE variants on the syngo console
    # like t1_mpr_ns_cor, t1_mpr_ns_sag, mprage, mpr_ns.
    r"mpr_ns", r"mprage", r"t1_mpr",
    # Generic
    r"bravo", r"3dt1",
]
T1_3D_RE = re.compile("|".join(T1_3D_KEYWORDS), re.IGNORECASE)

# Exclude these even if the keyword filter matched (localizers, 2D SE).
T1_3D_EXCLUDE = re.compile(
    r"\bscout\b|\blocalizer\b|\btse\b|t1_se|t1_tse",
    re.IGNORECASE,
)


def is_t1_3d_series(row: pd.Series) -> bool:
    """True if this row is a 3D structural T1 suitable for FastSurfer."""
    desc = str(row.get("SeriesDescription", "") or "")
    proto = str(row.get("ProtocolName", "") or "")
    combined = f"{desc} {proto}"
    if not T1_3D_RE.search(combined):
        return False
    # Exclude localizers / 2D SE sequences with 'T1' in the name.
    if T1_3D_EXCLUDE.search(combined):
        return False
    # Enough slices for 3D.
    n = row.get("NFiles")
    try:
        if pd.isna(n) or float(n) < 120:
            return False
    except (TypeError, ValueError):
        return False
    # Sub-2.5mm slice thickness (exclude anisotropic 2D T1s).
    st = row.get("SliceThickness")
    try:
        if pd.isna(st) or float(st) > 2.5:
            return False
    except (TypeError, ValueError):
        return False
    return True


# --- Protocol descriptor helpers ---

def normalize_desc(s: str) -> str:
    """Strip 'PARKVEST 08' and punctuation; return a snake_case token."""
    s = str(s or "").lower()
    s = re.sub(r"parkvest\s*0?\d*\s*", "", s)
    s = re.sub(r"[^a-z0-9]+", "_", s)
    return s.strip("_")


def short_software_version(sw_raw) -> str:
    """Reduce the multi-valued SoftwareVersions tag to one short token.

      "['NT 10.3.1', 'PIIM V2.1.4.1 MIMIT MCS']"         -> 'NT 10.3.1'
      "['10.6.2', '10.6.2.5', 'Gyroscan PMS/DICOM...']"  -> '10.6.2'
      "syngo MR 2004A 4VA25A"                            -> 'syngo 2004A'
      "syngo MR A30 4VA30A"                              -> 'syngo A30'
    """
    if sw_raw is None or (isinstance(sw_raw, float) and pd.isna(sw_raw)):
        return "unknown"
    s = str(sw_raw).strip(" []'\"")
    first = s.split(",")[0].strip(" '\"")
    if not first:
        return "unknown"
    # For Siemens 'syngo MR 2004A 4VA25A' pick the date-ish token.
    sy = re.match(r"syngo\s*MR\s*(\S+)", first, re.IGNORECASE)
    if sy:
        return f"syngo {sy.group(1)}"
    return first


def is_modern_protocol(row: pd.Series) -> bool:
    """High-res 3D T1 (sub-mm slice, low flip); used as a tie-breaker
    when a subject has multiple candidate T1s at the same timepoint."""
    try:
        st = float(row["SliceThickness"])
        fa = float(row["FlipAngle"])
        return st <= 1.5 and fa < 15
    except (TypeError, ValueError):
        return False


def round_or_na(x, dp=2) -> str:
    try:
        return str(round(float(x), dp))
    except (TypeError, ValueError):
        return "NA"


# --- Batch ID construction ---

def full_scanner_id(row: pd.Series) -> str:
    """Fine-grained protocol identifier: manufacturer|model|field|
    sequence|slice|flip|SW."""
    return "|".join([
        str(row.get("Manufacturer", "NA") or "NA"),
        str(row.get("ManufacturerModelName", "NA") or "NA"),
        round_or_na(row.get("MagneticFieldStrength"), 1),
        normalize_desc(row.get("SeriesDescription")),
        round_or_na(row.get("SliceThickness"), 2),
        round_or_na(row.get("FlipAngle"), 1),
        short_software_version(row.get("SoftwareVersions")),
    ])


def collapsed_batch_id(row: pd.Series) -> str:
    """Collapsed batch label for longCombat.

      SUS_BL_FFE      Philips 3D-FFE at BL (2mm, flip-30)
      SUS_TFE         Philips sT1W/3D/TFE at 3Y/5Y (1mm, flip-8 MPRAGE).
                      Collapses NT10 vs 10.6.2 SW variants.
      Bergen_MPR_cor  Siemens t1_mpr_ns_cor (coronal MPRAGE, TR~2130, flip 15)
      Bergen_MPR_sag  Siemens t1_mpr_ns_sag (sagittal MPRAGE, TR~1950, flip 8)
      Forde_MPR_cor   Siemens t1_mpr_ns_cor at Forde
      Forde_MPR_sag   Siemens t1_mpr_ns_sag at Forde

    Bergen and Forde stay separate: same protocol and model but
    different physical scanners.
    """
    site = str(row.get("Site", ""))
    desc = str(row.get("SeriesDescription", "") or "").lower()
    if site == "SUS":
        if "ffe" in desc:
            return "SUS_BL_FFE"
        if "tfe" in desc:
            return "SUS_TFE"
        return f"SUS_other_{normalize_desc(desc)}"
    if site in {"Bergen", "Forde"}:
        if "cor" in desc:
            return f"{site}_MPR_cor"
        if "sag" in desc:
            return f"{site}_MPR_sag"
        return f"{site}_other_{normalize_desc(desc)}"
    return f"{site}_unknown"


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--inventory", type=Path, required=True,
        help="dicom_inventory_all.csv from 01_dicom_inventory.py",
    )
    p.add_argument(
        "--cohort", type=Path, required=True,
        help="cohort.csv from 00. Used to tag in_cohort=True for "
             "subjects in the FastSurfer cohort.",
    )
    p.add_argument(
        "--output", type=Path, required=True,
        help="Path to write scanner_metadata.csv",
    )
    args = p.parse_args()

    inv = pd.read_csv(args.inventory)
    inv["NFiles"] = pd.to_numeric(inv["NFiles"], errors="coerce")

    cohort = pd.read_csv(args.cohort)
    if "Subject" not in cohort.columns:
        print(
            f"ERROR: --cohort CSV must have a 'Subject' column. "
            f"Columns found: {list(cohort.columns)}",
            file=sys.stderr,
        )
        return 2
    cohort_subjects = set(cohort["Subject"].astype(str))
    print(f"Cohort: {len(cohort_subjects)} subjects from {args.cohort}")

    # Filter to 3D T1 structural series.
    mask = inv.apply(is_t1_3d_series, axis=1)
    t1 = inv[mask].copy()
    print(f"[1] {len(t1)} candidate 3D T1 series "
          f"covering {t1['Subject'].nunique()} subjects across "
          f"{t1['Timepoint'].nunique()} timepoints")
    print(f"    Series descriptions kept:")
    for desc, n in t1["SeriesDescription"].value_counts().items():
        print(f"      n={n:3d}  {desc}")

    # Resolve duplicate scans within a timepoint: prefer modern (1mm iso,
    # low flip), then highest slice count. Some subjects (e.g. S026)
    # have both old-FFE and new-TFE at the same session.
    t1["is_modern"] = t1.apply(is_modern_protocol, axis=1)
    t1_one = (
        t1.sort_values(["is_modern", "NFiles"], ascending=[False, False])
          .drop_duplicates(["Site", "Group", "Subject", "Timepoint"],
                           keep="first")
          .sort_values(["Site", "Group", "Subject", "Timepoint"])
          .reset_index(drop=True)
    )
    print(f"[2] After de-duplication: {len(t1_one)} rows "
          f"(one per Site x Group x Subject x Timepoint)")

    # Build ScannerID (full) and BatchID (collapsed).
    t1_one["SoftwareShort"] = t1_one["SoftwareVersions"].apply(
        short_software_version
    )
    t1_one["SeriesDescNormalized"] = t1_one["SeriesDescription"].apply(
        normalize_desc
    )
    t1_one["ScannerID"] = t1_one.apply(full_scanner_id, axis=1)
    t1_one["BatchID"]   = t1_one.apply(collapsed_batch_id, axis=1)

    # Tag in_cohort from cohort.csv.
    t1_one["in_cohort"] = t1_one["Subject"].isin(cohort_subjects)

    n_rows_in_cohort = t1_one["in_cohort"].sum()
    n_subs_in_cohort = t1_one.loc[t1_one["in_cohort"], "Subject"].nunique()
    print(f"[3] Tagged {n_rows_in_cohort} rows as in_cohort=True "
          f"({n_subs_in_cohort} / {len(cohort_subjects)} cohort subjects "
          f"matched)")

    unmatched = cohort_subjects - set(t1_one["Subject"].unique())
    if unmatched:
        print(f"    WARN: {len(unmatched)} cohort subject(s) not "
              f"matched in DICOM inventory.")
        print(f"    First 10 unmatched: "
              f"{sorted(unmatched)[:10]}")

    # Reorder output columns.
    out_cols = [
        "Site", "Group", "Subject", "Subject_raw", "Timepoint",
        "BatchID", "ScannerID", "in_cohort",
        "Manufacturer", "ManufacturerModelName", "DeviceSerialNumber",
        "StationName", "SoftwareShort", "SoftwareVersions",
        "MagneticFieldStrength", "SeriesDescription", "SeriesDescNormalized",
        "ProtocolName", "SliceThickness", "PixelSpacing", "FlipAngle",
        "RepetitionTime", "EchoTime", "InversionTime", "AcquisitionMatrix",
        "ReceiveCoilName", "ScanningSequence", "SequenceVariant",
        "InstitutionName", "AcquisitionDate", "StudyDate",
        "PatientSex", "PatientAge",
        "SeriesFolder", "NFiles",
    ]
    out_cols = [c for c in out_cols if c in t1_one.columns]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    t1_one[out_cols].to_csv(args.output, index=False)
    print(f"[4] Wrote {args.output} ({len(t1_one)} rows)")

    # Diagnostics.
    print("\nBatchID x Site x Timepoint distribution (all subjects):")
    tbl = (
        t1_one.groupby(["Site", "Timepoint", "BatchID"])
              .size()
              .unstack("BatchID", fill_value=0)
    )
    print(tbl.to_string())

    print("\nSame, restricted to in_cohort=True subjects:")
    cohort_rows = t1_one[t1_one["in_cohort"]]
    if len(cohort_rows):
        tbl_c = (
            cohort_rows.groupby(["Site", "Timepoint", "BatchID"])
                       .size()
                       .unstack("BatchID", fill_value=0)
        )
        print(tbl_c.to_string())
    else:
        print("  (no in_cohort rows)")

    print("\nSubjects with N distinct BatchIDs across their timepoints:")
    per = t1_one.groupby("Subject")["BatchID"].nunique()
    print(per.value_counts().sort_index().to_string())

    print("\nBatch sizes (longCombat adequacy check):")
    bc = t1_one["BatchID"].value_counts()
    for bid, n in bc.items():
        print(f"  n={n:3d}  {bid}")
    print(f"  min batch size: {bc.min()}")
    if (bc < 5).any():
        print(f"  Batches with n<5 are unreliable for ComBat:")
        print(f"    {bc[bc<5].to_dict()}")
    if (bc < 3).any():
        print(f"  Batches with n<3 will be rejected by longCombat.")

    print("\nFine-grained ScannerIDs (for reference):")
    for sid, n in t1_one["ScannerID"].value_counts().items():
        print(f"  n={n:3d}  {sid}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
