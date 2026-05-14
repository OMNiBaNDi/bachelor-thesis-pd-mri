#!/usr/bin/env python3
"""Merge per-subject clinical covariates onto per-scan scanner metadata.

The output is one row per scan with these columns added:

    age_at_BL              age at baseline MRI
    PatientSex_clinical    'M' or 'F' from SPSS BL_Sex
    Group_clinical         'Pasient' or 'Control'
    bl_mri_date            'YYYY-MM-DD' from SPSS BL_MRI_date. Optional;
                           used as the BL date fallback for SUS, where
                           DICOM headers were de-identified.

03 uses different names (sex, group_label) so those get renamed here.

Usage:
    python scripts/stage_a_metadata/04_merge_metadata.py \\
        --scanner-metadata outputs/stage_a_metadata/scanner_metadata.csv \\
        --clinical         outputs/stage_a_metadata/clinical_covariates.csv \\
        --output           outputs/stage_a_metadata/scanner_metadata_with_covariates.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--scanner-metadata", type=Path, required=True,
        help="scanner_metadata.csv from 02 (one row per scan).",
    )
    p.add_argument(
        "--clinical", type=Path, required=True,
        help="clinical_covariates.csv from 03 (one row per subject).",
    )
    p.add_argument(
        "--output", type=Path, required=True,
        help="Output CSV.",
    )
    args = p.parse_args()

    meta = pd.read_csv(args.scanner_metadata)
    clin = pd.read_csv(args.clinical)
    print(f"scanner_metadata.csv: {len(meta):5d} rows, "
          f"{meta['Subject'].nunique()} unique subjects")
    print(f"clinical_covariates:  {len(clin):5d} rows, "
          f"{clin['Subject'].nunique()} unique subjects")

    # Required clinical columns. bl_mri_date is optional.
    needed_clin = ["Subject", "age_at_BL", "sex", "group_label"]
    missing = [c for c in needed_clin if c not in clin.columns]
    if missing:
        print(f"ERROR: --clinical CSV is missing required column(s): "
              f"{missing}. Columns found: {list(clin.columns)}",
              file=sys.stderr)
        return 2

    optional_clin = ["bl_mri_date"]
    optional_present = [c for c in optional_clin if c in clin.columns]
    optional_missing = [c for c in optional_clin if c not in clin.columns]
    if optional_missing:
        print(f"  --clinical CSV doesn't have optional column(s) "
              f"{optional_missing}.")

    if "in_cohort" not in meta.columns:
        print(f"ERROR: --scanner-metadata CSV is missing the 'in_cohort' "
              f"column. Columns found: {list(meta.columns)}",
              file=sys.stderr)
        return 2

    # Pick the columns we want and rename to the pipeline-wide names.
    clin_slim = clin[needed_clin + optional_present].rename(columns={
        "sex":         "PatientSex_clinical",
        "group_label": "Group_clinical",
    })

    merged = meta.merge(clin_slim, on="Subject", how="left")
    print(f"\nAfter merge: {len(merged)} rows")

    n_with_age = merged["age_at_BL"].notna().sum()
    n_with_sex = merged["PatientSex_clinical"].notna().sum()
    print(f"  Rows with age_at_BL:           {n_with_age} / {len(merged)}")
    print(f"  Rows with PatientSex_clinical: {n_with_sex} / {len(merged)}")

    missing_clin = merged[merged["age_at_BL"].isna()]
    if len(missing_clin):
        missing_subs = sorted(missing_clin["Subject"].unique())
        print(f"\n  Subjects missing clinical data ({len(missing_subs)}):")
        for s in missing_subs[:10]:
            n_scans = (missing_clin["Subject"] == s).sum()
            in_coh = merged.loc[merged["Subject"] == s, "in_cohort"].iloc[0]
            print(f"    {s}  ({n_scans} scans, in_cohort={in_coh})")
        if len(missing_subs) > 10:
            print(f"    ... and {len(missing_subs) - 10} more")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, index=False)
    print(f"\nWrote {args.output} ({len(merged)} rows, "
          f"{len(merged.columns)} columns)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
