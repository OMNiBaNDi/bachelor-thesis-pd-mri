#!/usr/bin/env python3
"""Extract per-subject clinical covariates from the ParkWest SPSS file.

Pulls age at baseline, sex, and group for the subjects in the analysis
cohort. The output (clinical_covariates.csv) is consumed by
04_merge_metadata.py.

Run from the pd_thesis venv (pyreadstat isn't in the rpy2 venv used
elsewhere in the pipeline). 04 renames the output columns (sex,
group_label) to the pipeline-wide names on merge.

Usage:
    python scripts/stage_a_metadata/03_clinical_covariates.py \\
        --sav     /nfs/br1_prosjekt/ParkWest/ClinicalData/ParkVest_V02-V21_191219.sav \\
        --cohort  outputs/stage_a_metadata/cohort.csv \\
        --output  outputs/stage_a_metadata/clinical_covariates.csv

    # Exploration: list columns and value labels:
    python scripts/stage_a_metadata/03_clinical_covariates.py \\
        --sav     .../ParkVest_V02-V21_191219.sav \\
        --cohort  outputs/stage_a_metadata/cohort.csv \\
        --output  /dev/null \\
        --explore
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

from pipeline_lib.ids import normalize_subject_id

try:
    import pyreadstat
except ImportError:
    print(
        "ERROR: pyreadstat not installed. This script must run in the "
        "pd_thesis venv.",
        file=sys.stderr,
    )
    print(
        "If clinical_covariates.csv already exists, skip to "
        "04_merge_metadata.py.",
        file=sys.stderr,
    )
    sys.exit(1)


# Candidate column names per variable. First match wins; case-insensitive.
CANDIDATE_IMAGING_ID_COLS = [
    "BL_MRI_ID", "BL_MRI_I", "IDcoderForMRI", "MRI_ID", "IDcoder",
    "imagingID", "ImagingID",
]
CANDIDATE_AGE_COLS = [
    "BL_MRI_AGE", "BL_MRI_A", "BL_Age", "AGE_BL", "age_bl",
    "alder_bl", "AGE_MRI", "age_MRI",
]
CANDIDATE_SEX_COLS = [
    "BL_Sex", "SEX", "sex", "KJONN", "kjonn", "GENDER", "gender",
]
CANDIDATE_GROUP_COLS = [
    "BL_Type", "GROUP", "group", "GRUPPE", "gruppe", "CASE", "case",
    "STATUS", "DIAGNOSIS", "diagnosis", "BL_TYPE",
]
CANDIDATE_CLINICAL_ID_COLS = [
    "BL_CASE", "CASE_ID", "id", "subjid", "lopenr", "PatientID",
    "SUBJECT", "subject",
]
CANDIDATE_CENTER_COLS = [
    "BL_Center", "BL_centre", "Center", "centre", "Site", "site",
]
CANDIDATE_BL_MRI_DATE_COLS = [
    # Primary names in the ParkVest SPSS export. Used as the BL date
    # fallback for SUS (DICOM was de-identified).
    "BL_MRI_date", "BL_MRI_Date",
    # Fallbacks in case a different export uses another name.
    "BL_DateVisit", "BL_VisitDate", "BL_Date",
]


def find_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first candidate that exists in df (case-insensitive)."""
    lower_to_real = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_to_real:
            return lower_to_real[cand.lower()]
    return None


def fuzzy_search_column(df: pd.DataFrame, keywords: list[str]) -> list[str]:
    """Return columns whose names contain any of the keywords."""
    hits: list[str] = []
    for c in df.columns:
        cl = c.lower()
        for kw in keywords:
            if kw.lower() in cl:
                hits.append(c)
                break
    return hits


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--sav", type=Path, required=True,
        help="ParkVest .sav file from ClinicalData/",
    )
    p.add_argument(
        "--fam", type=Path, default=None,
        help="Optional: PLINK .fam file for sex cross-check",
    )
    p.add_argument(
        "--cohort", type=Path, required=True,
        help="cohort.csv from 00 (only the Subject column is read).",
    )
    p.add_argument(
        "--output", type=Path, required=True,
        help="Where to write clinical_covariates.csv",
    )
    p.add_argument(
        "--explore", action="store_true",
        help="Print all columns + labels, don't extract anything",
    )
    args = p.parse_args()

    # Load SPSS.
    print(f"Loading SPSS file: {args.sav}")
    clin, meta = pyreadstat.read_sav(str(args.sav))
    print(f"  Shape: {clin.shape[0]} rows x {clin.shape[1]} columns")

    # Explore mode.
    if args.explore:
        print("\nAll variables:")
        for i, col in enumerate(clin.columns):
            label = meta.column_names_to_labels.get(col, "") or ""
            val_map = meta.variable_value_labels.get(col, {})
            vm_str = f"  [values: {val_map}]" if val_map else ""
            print(f"  {i:4d}  {col:30s}  {label}{vm_str}")
        print("\nFirst 3 rows (selected plausible columns):")
        plausible = []
        for cands in [CANDIDATE_CLINICAL_ID_COLS, CANDIDATE_IMAGING_ID_COLS,
                      CANDIDATE_AGE_COLS, CANDIDATE_SEX_COLS,
                      CANDIDATE_GROUP_COLS]:
            c = find_column(clin, cands)
            if c and c not in plausible:
                plausible.append(c)
        for extra in fuzzy_search_column(clin, ["mri", "age", "alder",
                                                "sex", "kjonn", "dob",
                                                "birth", "_id"]):
            if extra not in plausible:
                plausible.append(extra)
        print(clin[plausible].head(3).to_string())
        return 0

    # Locate key columns.
    imaging_id_col = find_column(clin, CANDIDATE_IMAGING_ID_COLS)
    age_col        = find_column(clin, CANDIDATE_AGE_COLS)
    sex_col        = find_column(clin, CANDIDATE_SEX_COLS)
    group_col      = find_column(clin, CANDIDATE_GROUP_COLS)
    clin_id_col    = find_column(clin, CANDIDATE_CLINICAL_ID_COLS)
    center_col     = find_column(clin, CANDIDATE_CENTER_COLS)
    bl_date_col    = find_column(clin, CANDIDATE_BL_MRI_DATE_COLS)

    print("\nColumn detection:")
    print(f"  Imaging ID column (BL_MRI_ID):     {imaging_id_col}")
    print(f"  Age at BL column (BL_MRI_AGE):     {age_col}")
    print(f"  Sex column (BL_Sex):               {sex_col}")
    print(f"  Group / diagnosis column:          {group_col}")
    print(f"  Clinical ID column (BL_CASE):      {clin_id_col}")
    print(f"  Center column (BL_Center):         {center_col}")
    print(f"  BL MRI date column:                {bl_date_col}")

    if imaging_id_col is None:
        print("\nWARN: could not find an imaging-ID link column automatically.",
              file=sys.stderr)
        print(f"  Candidates tried: {CANDIDATE_IMAGING_ID_COLS}", file=sys.stderr)
        print("\n  Columns in the SPSS file that mention 'MRI' or 'ID':",
              file=sys.stderr)
        for c in fuzzy_search_column(clin, ["mri", "id", "coder"]):
            lbl = meta.column_names_to_labels.get(c, "")
            print(f"    {c:30s}  {lbl}", file=sys.stderr)
        return 1

    if age_col is None or sex_col is None:
        print("\nWARN: missing age or sex column. Re-run with --explore.",
              file=sys.stderr)
        return 1

    # ID prefix distribution.
    def prefix_of(s):
        if not isinstance(s, str):
            return None
        m = re.match(r"^([A-Za-z]+)", s)
        return m.group(1) if m else None

    prefix_counts = clin[imaging_id_col].apply(prefix_of).value_counts(dropna=False)
    print(f"\nID prefix distribution in {imaging_id_col}:")
    print(prefix_counts.head(20).to_string())

    if center_col is not None:
        center_vmap = meta.variable_value_labels.get(center_col, {})
        if center_vmap:
            print(f"\nSubjects per center ({center_col}):")
            center_counts = clin[center_col].map(center_vmap).value_counts(dropna=False)
            print(center_counts.to_string())

    # Load cohort subject IDs.
    cohort = pd.read_csv(args.cohort)
    if "Subject" not in cohort.columns:
        print(f"ERROR: --cohort CSV must have a 'Subject' column. "
              f"Columns found: {list(cohort.columns)}",
              file=sys.stderr)
        return 2
    imaging_subs = set(cohort["Subject"].astype(str).unique())
    print(f"\nCohort: {len(imaging_subs)} distinct subjects from {args.cohort}")

    # Normalize SPSS imaging IDs.
    # The SPSS file uses short IDs (B03, BK07, S02, SK01); cohort.csv
    # uses 3-digit IDs (B003, BK007, ...). Normalize so they join.
    def _normalize_spss_id(x):
        if pd.isna(x):
            return None
        s = str(x).strip()
        return normalize_subject_id(s) if s else None

    clin["imaging_id_normalized"] = clin[imaging_id_col].apply(_normalize_spss_id)

    # Match cohort subjects to clinical rows.
    matched = clin[clin["imaging_id_normalized"].isin(imaging_subs)]
    print(f"Matched {matched['imaging_id_normalized'].nunique()} of "
          f"{len(imaging_subs)} cohort subjects to clinical records.")
    unmatched = imaging_subs - set(matched["imaging_id_normalized"])
    if unmatched:
        print(f"Unmatched cohort subjects: {sorted(unmatched)}")

    # Extract the fields we want.
    keep = [imaging_id_col, "imaging_id_normalized", age_col, sex_col]
    if group_col:
        keep.append(group_col)
    if clin_id_col:
        keep.append(clin_id_col)
    if center_col:
        keep.append(center_col)
    if bl_date_col:
        keep.append(bl_date_col)
    out = matched[keep].copy()

    # 04 renames sex_clinical / group_clinical to the pipeline-wide form.
    out = out.rename(columns={
        "imaging_id_normalized": "Subject",
        age_col: "age_at_BL",
        sex_col: "sex_clinical",
    })
    if group_col:
        out = out.rename(columns={group_col: "group_clinical"})
    if clin_id_col:
        out = out.rename(columns={clin_id_col: "clinical_id"})
    if center_col:
        out = out.rename(columns={center_col: "center_clinical"})
    if bl_date_col:
        out = out.rename(columns={bl_date_col: "bl_mri_date"})
    out = out.rename(columns={imaging_id_col: "imaging_id_raw"})

    # Collapse multi-visit rows to one row per subject.
    out = (out.sort_values("Subject")
              .groupby("Subject", as_index=False)
              .agg({c: "first" for c in out.columns if c != "Subject"}))
    print(f"\nAfter dedup: {len(out)} rows ({out['Subject'].nunique()} subjects)")

    # Decode sex to 'M' / 'F' to match the DICOM PatientSex convention.
    # SPSS BL_Sex is 0.0='Man', 1.0='Woman'.
    sex_vmap = meta.variable_value_labels.get(sex_col, {})
    print(f"\nSex value map from SPSS ({sex_col}): {sex_vmap}")

    def decode_sex(v):
        if pd.isna(v):
            return None
        label = sex_vmap.get(v, v) if sex_vmap else v
        if isinstance(label, str):
            low = label.lower().strip()
            if low in ("man", "m", "male", "mann", "herre"):
                return "M"
            if low in ("woman", "w", "f", "female", "kvinne", "kvinn", "dame"):
                return "F"
            if "female" in low or "woman" in low or "kvinn" in low:
                return "F"
            if "male" in low and "female" not in low:
                return "M"
        try:
            vf = float(v)
            if vf == 0.0: return "M"
            if vf == 1.0:
                # 1.0 is 'Woman' in BL_Sex but 'Male' in PLINK. If no
                # value map matched, assume PLINK.
                return "M" if not sex_vmap else "F"
            if vf == 2.0: return "F"  # PLINK-style
        except (TypeError, ValueError):
            pass
        return None

    out["sex"] = out["sex_clinical"].apply(decode_sex)

    # Normalize bl_mri_date to "YYYY-MM-DD". SAV returns datetime.date;
    # DICOM dates elsewhere are YYYYMMDD.
    if "bl_mri_date" in out.columns:
        def _format_date(v):
            if pd.isna(v):
                return None
            if hasattr(v, "strftime"):  # datetime.date / Timestamp
                return v.strftime("%Y-%m-%d")
            s = str(v).strip()
            if not s or s.lower() in ("nan", "nat", "none"):
                return None
            # Try parsing already-formatted strings; pass through on
            # parse failure so the bad value stays visible.
            try:
                return pd.to_datetime(s).strftime("%Y-%m-%d")
            except (ValueError, TypeError):
                return s
        out["bl_mri_date"] = out["bl_mri_date"].apply(_format_date)
        n_filled = out["bl_mri_date"].notna().sum()
        print(f"\nbl_mri_date populated: {n_filled} / {len(out)} subjects")

    # Decode group and center labels for readability.
    if group_col:
        group_vmap = meta.variable_value_labels.get(group_col, {})
        if group_vmap:
            out["group_label"] = out["group_clinical"].map(group_vmap)
    if center_col:
        center_vmap = meta.variable_value_labels.get(center_col, {})
        if center_vmap:
            out["center_label"] = out["center_clinical"].map(center_vmap)

    # Optional .fam cross-check.
    if args.fam is not None and args.fam.exists():
        fam = pd.read_csv(args.fam, sep=r"\s+", header=None,
                          names=["fid", "iid", "pid", "mid",
                                 "sex_fam", "pheno"])
        fam["sex_fam_decoded"] = fam["sex_fam"].map({1: "M", 2: "F", 0: None})
        if "clinical_id" in out.columns:
            fam_ids = set(fam["iid"].astype(str))
            out_ids = set(out["clinical_id"].astype(str))
            overlap = out_ids & fam_ids
            print(f"\nFam cross-check: {len(overlap)} of {len(out_ids)} "
                  f"clinical IDs found in .fam file")
            if overlap:
                fam_sex = fam.set_index("iid")["sex_fam_decoded"]
                out["sex_fam"] = out["clinical_id"].astype(str).map(fam_sex)
                mismatch = out[(out["sex_fam"].notna()) &
                               (out["sex"].notna()) &
                               (out["sex_fam"] != out["sex"])]
                if len(mismatch):
                    print(f"  WARN: {len(mismatch)} sex mismatches between "
                          f".sav and .fam:")
                    print(mismatch[["Subject", "clinical_id",
                                    "sex", "sex_fam"]].to_string(index=False))
                else:
                    print(f"  All {len(out[out['sex_fam'].notna()])} "
                          f"overlapping subjects have consistent sex.")

    # Reorder and write.
    out_cols = [c for c in ["Subject", "imaging_id_raw", "clinical_id",
                            "age_at_BL", "sex", "sex_clinical",
                            "sex_fam", "group_clinical", "group_label",
                            "center_clinical", "center_label",
                            "bl_mri_date"]
                if c in out.columns]
    out = out[out_cols].sort_values("Subject").reset_index(drop=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(f"\nWrote {args.output}")

    print(f"\nSummary:")
    print(f"Subjects with age_at_BL: {out['age_at_BL'].notna().sum()}")
    print(f"Subjects with sex:       {out['sex'].notna().sum()}")
    print(f"Sex distribution:")
    print(out["sex"].value_counts(dropna=False).to_string())
    if "center_label" in out.columns:
        print(f"\nCenter distribution:")
        print(out["center_label"].value_counts(dropna=False).to_string())
    if "group_label" in out.columns:
        print(f"\nGroup distribution:")
        print(out["group_label"].value_counts(dropna=False).to_string())
    if out["age_at_BL"].notna().any():
        print(f"\nAge at BL: mean={out['age_at_BL'].mean():.1f}, "
              f"sd={out['age_at_BL'].std():.1f}, "
              f"range=[{out['age_at_BL'].min():.0f},"
              f" {out['age_at_BL'].max():.0f}]")

    print(f"\nFirst 10 rows:")
    print(out.head(10).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
