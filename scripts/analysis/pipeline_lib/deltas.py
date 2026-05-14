"""Per-subject deltas over the 12-ROI panel, plus elapsed-time computation.

Used by 05_extract.py for the unharmonized deltas and by
06_harmonize.py for the harmonized deltas and interval_summary.csv.
The functions work on either harmonized or unharmonized inputs since
the schemas are preserved through harmonization.

The ROI set comes from pipeline_lib.constants.ROI_PANEL_THESIS.

subject_roi_deltas.csv columns: Subject, Site, Group, ROI, Tissue,
Cluster, Measure, Delta_Window, Value_t1, Value_t2, Delta_abs, Delta_pct.

interval_summary.csv columns: Subject, Delta_Window, bl_date, t2_date,
interval_years, interval_source. interval_source is 'dicom',
'sav_bl_mri_date', or 'fallback_nominal'.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from pipeline_lib.constants import (
    DELTA_WINDOWS,
    ROI_PANEL_THESIS,
    TIMEPOINT_YEARS,
    TIMEPOINTS,
)


# --- Cortical bilateral averaging ---

def _prepare_cortical_for_analysis(
    cortical_df: pd.DataFrame,
) -> pd.DataFrame:
    """Collapse the cortical table to one row per (Subject, Timepoint, StructName).

    ThickAvg_bilateral is the mean of left and right hemispheres; if
    only one hemisphere is present, the mean reduces to that single
    value.

    Returns columns: Subject, Timepoint, StructName, ThickAvg_bilateral.
    """
    if cortical_df.empty:
        return pd.DataFrame(columns=[
            "Subject", "Timepoint", "StructName", "ThickAvg_bilateral",
        ])

    df = cortical_df[
        ["Subject", "Timepoint", "Hemisphere", "StructName", "ThickAvg"]
    ].copy()

    bilateral = (
        df.groupby(["Subject", "Timepoint", "StructName"], as_index=False)
          ["ThickAvg"]
          .mean()
          .rename(columns={"ThickAvg": "ThickAvg_bilateral"})
    )
    return bilateral


# --- Per-subject deltas ---

def compute_analysis_deltas(
    cohort_wide: pd.DataFrame,
    cohort_cortical: pd.DataFrame,
) -> pd.DataFrame:
    """Compute subject-level deltas for every ROI x window.

    Only subjects with Data_Present=True at all three timepoints
    (BL, 3Y, 5Y) are included.

        Delta_abs = value(t2) - value(t1)
        Delta_pct = 100 * (value(t2) - value(t1)) / value(t1)   (NaN if v1==0)

    Returns the long-format DataFrame described in the module docstring,
    or an empty DataFrame if no subject has all three timepoints.
    """
    rows: List[Dict] = []

    # Identify complete subjects (BL + 3Y + 5Y all present).
    complete = (
        cohort_wide[cohort_wide["Data_Present"]]
        .groupby(["Site", "Group", "Subject"])["Timepoint"]
        .apply(set)
    )
    complete_triples = [
        (site, group, subj)
        for (site, group, subj), tps in complete.items()
        if {"BL", "3Y", "5Y"}.issubset(tps)
    ]
    incomplete_count = len(complete) - len(complete_triples)
    if incomplete_count:
        print(
            f"  {incomplete_count} subject(s) missing one or more of "
            f"(BL, 3Y, 5Y); excluded from deltas."
        )

    if not complete_triples:
        print("  No subjects have all 3 timepoints; deltas empty.")
        return pd.DataFrame()

    complete_ids = {subj for _, _, subj in complete_triples}
    site_lookup  = {subj: site  for site, _,     subj in complete_triples}
    group_lookup = {subj: group for _,    group, subj in complete_triples}

    # Subset to complete subjects.
    sub_wide = cohort_wide[
        cohort_wide["Subject"].isin(complete_ids)
        & cohort_wide["Data_Present"]
    ].copy()

    cortical_bilateral = _prepare_cortical_for_analysis(
        cohort_cortical[cohort_cortical["Subject"].isin(complete_ids)]
        if not cohort_cortical.empty else pd.DataFrame()
    )

    # Compute deltas per panel entry.
    n_windows_skipped_nan = 0
    n_subcortical_col_missing = 0
    n_cortical_region_missing = 0

    for entry in ROI_PANEL_THESIS:
        tissue       = entry["tissue"]
        roi_display  = entry["roi"]
        cluster      = entry["cluster"]

        if tissue == "subcortical":
            col = entry["analysis_col"]   # e.g. 'Putamen_Mean_norm'
            if col not in sub_wide.columns:
                n_subcortical_col_missing += 1
                continue

            measure = "volume_norm"

            for subj in complete_ids:
                subj_data = sub_wide[sub_wide["Subject"] == subj]
                if subj_data.empty:
                    continue
                tp_vals = subj_data.set_index("Timepoint")[col]

                for tp1, tp2 in DELTA_WINDOWS:
                    if tp1 not in tp_vals.index or tp2 not in tp_vals.index:
                        continue
                    v1 = tp_vals[tp1]
                    v2 = tp_vals[tp2]
                    if pd.isna(v1) or pd.isna(v2):
                        n_windows_skipped_nan += 1
                        continue

                    delta_abs = v2 - v1
                    delta_pct = (100.0 * (v2 - v1) / v1
                                 if v1 != 0 else np.nan)

                    rows.append({
                        "Subject":      subj,
                        "Site":         site_lookup[subj],
                        "Group":        group_lookup[subj],
                        "ROI":          roi_display,
                        "Tissue":       "subcortical",
                        "Cluster":      cluster,
                        "Measure":      measure,
                        "Delta_Window": f"{tp1}→{tp2}",
                        "Value_t1":     v1,
                        "Value_t2":     v2,
                        "Delta_abs":    delta_abs,
                        "Delta_pct":    delta_pct,
                    })

        else:  # cortical
            # Cortical reads from the bilateralized long-format table.
            fs_region = entry["fs_region"]   # e.g. 'entorhinal'
            region_data = cortical_bilateral[
                cortical_bilateral["StructName"] == fs_region
            ]
            if region_data.empty:
                n_cortical_region_missing += 1
                continue

            measure = "thickness_mm"

            for subj in complete_ids:
                subj_data = region_data[region_data["Subject"] == subj]
                if subj_data.empty:
                    continue
                tp_vals = subj_data.set_index("Timepoint")["ThickAvg_bilateral"]

                for tp1, tp2 in DELTA_WINDOWS:
                    if tp1 not in tp_vals.index or tp2 not in tp_vals.index:
                        continue
                    v1 = tp_vals[tp1]
                    v2 = tp_vals[tp2]
                    if pd.isna(v1) or pd.isna(v2):
                        n_windows_skipped_nan += 1
                        continue

                    delta_abs = v2 - v1
                    delta_pct = (100.0 * (v2 - v1) / v1
                                 if v1 != 0 else np.nan)

                    rows.append({
                        "Subject":      subj,
                        "Site":         site_lookup[subj],
                        "Group":        group_lookup[subj],
                        "ROI":          roi_display,
                        "Tissue":       "cortical",
                        "Cluster":      cluster,
                        "Measure":      measure,
                        "Delta_Window": f"{tp1}→{tp2}",
                        "Value_t1":     v1,
                        "Value_t2":     v2,
                        "Delta_abs":    delta_abs,
                        "Delta_pct":    delta_pct,
                    })



    return pd.DataFrame(rows)


# --- Long-format builder for LME analysis ---

def build_cohort_long(
    cohort_wide: pd.DataFrame,
    cohort_cortical: pd.DataFrame,
    scanner_meta: pd.DataFrame,
    interval_summary: pd.DataFrame,
) -> pd.DataFrame:
    """Build the long-format dataset for the LME analysis.

    One row per (Subject, Timepoint, ROI), with the regional value plus
    the covariates the LME needs. Works on either harmonized or
    unharmonized inputs; the typical use is on harmonized inputs to
    produce cohort_long_harmonized.csv. The unharmonized variant is
    built for the Risk-A sensitivity check.

    Years_from_BL is the actual elapsed time from
    compute_interval_years, not the nominal 0/3/5 stored in cohort_wide.
    The LME needs real elapsed time so subjects with off-nominal
    follow-up intervals contribute the right slope.

    Required columns:
        cohort_wide:      Subject, Timepoint, Years_from_BL, Data_Present,
                          plus the analysis_col for each subcortical
                          entry of ROI_PANEL_THESIS.
        cohort_cortical:  Subject, Timepoint, Hemisphere, StructName,
                          ThickAvg.
        scanner_meta:     Subject, Timepoint, BatchID, age_at_BL,
                          PatientSex_clinical, Group_clinical.
        interval_summary: Subject, Delta_Window, interval_years.

    Subjects missing one or more timepoints, and rows with NaN values
    or missing covariates, are dropped with a count printed.

    Returns columns: Subject, Timepoint, Years_from_BL, ROI, Tissue,
    Measure, Value, Group_clinical, age_at_BL, PatientSex_clinical,
    BatchID.
    """
    rows: List[Dict] = []

    # Validate input schemas.
    required_wide = {"Subject", "Timepoint", "Data_Present"}
    missing_wide = required_wide - set(cohort_wide.columns)
    if missing_wide:
        raise ValueError(
            f"build_cohort_long: cohort_wide is missing "
            f"required columns: {sorted(missing_wide)}"
        )

    required_cort = {"Subject", "Timepoint", "Hemisphere",
                     "StructName", "ThickAvg"}
    missing_cort = required_cort - set(cohort_cortical.columns)
    if missing_cort and not cohort_cortical.empty:
        raise ValueError(
            f"build_cohort_long: cohort_cortical is "
            f"missing required columns: {sorted(missing_cort)}"
        )

    required_meta = {"Subject", "Timepoint", "BatchID", "age_at_BL",
                     "PatientSex_clinical", "Group_clinical"}
    missing_meta = required_meta - set(scanner_meta.columns)
    if missing_meta:
        raise ValueError(
            f"build_cohort_long: scanner_meta is missing "
            f"required columns: {sorted(missing_meta)}"
        )

    required_iv = {"Subject", "Delta_Window", "interval_years"}
    missing_iv = required_iv - set(interval_summary.columns)
    if missing_iv:
        raise ValueError(
            f"build_cohort_long: interval_summary is "
            f"missing required columns: {sorted(missing_iv)}"
        )

    # Identify complete subjects.
    complete = (
        cohort_wide[cohort_wide["Data_Present"]]
        .groupby("Subject")["Timepoint"]
        .apply(set)
    )
    complete_ids = {
        subj for subj, tps in complete.items()
        if set(TIMEPOINTS).issubset(tps)
    }
    incomplete_count = len(complete) - len(complete_ids)
    if incomplete_count:
        print(
            f"  {incomplete_count} subject(s) missing one or more of "
            f"{TIMEPOINTS}; excluded from long format."
        )

    if not complete_ids:
        print("  No subjects have all 3 timepoints; long-format empty.")
        return pd.DataFrame()

    # Subset and prepare.
    sub_wide = cohort_wide[
        cohort_wide["Subject"].isin(complete_ids)
        & cohort_wide["Data_Present"]
    ].copy()

    cortical_bilateral = _prepare_cortical_for_analysis(
        cohort_cortical[cohort_cortical["Subject"].isin(complete_ids)]
        if not cohort_cortical.empty else pd.DataFrame()
    )

    # (Subject, Timepoint) -> covariates lookup. BatchID is per-scan;
    # the other covariates are subject-invariant. drop_duplicates is
    # defensive.
    meta_slim = (
        scanner_meta[
            ["Subject", "Timepoint", "BatchID",
             "age_at_BL", "PatientSex_clinical", "Group_clinical"]
        ]
        .drop_duplicates(subset=["Subject", "Timepoint"], keep="first")
        .set_index(["Subject", "Timepoint"])
    )

    # (Subject, Timepoint) -> actual elapsed years. BL is 0; follow-up
    # timepoints come from interval_summary's second TP in Delta_Window.
    years_lookup: Dict[Tuple[str, str], float] = {}
    for subj in complete_ids:
        years_lookup[(subj, "BL")] = 0.0
    for _, row in interval_summary.iterrows():
        subj = row["Subject"]
        if subj not in complete_ids:
            continue
        try:
            tp2 = row["Delta_Window"].split("→")[1]
        except (AttributeError, IndexError):
            continue
        years_lookup[(subj, tp2)] = float(row["interval_years"])

    # Iterate over (panel, subject, timepoint).
    n_value_nan = 0
    n_subcortical_col_missing = 0
    n_cortical_region_missing = 0
    n_years_missing = 0
    n_meta_missing = 0

    for entry in ROI_PANEL_THESIS:
        roi_display = entry["roi"]
        tissue = entry["tissue"]

        if tissue == "subcortical":
            col = entry["analysis_col"]   # e.g. "Hippocampus_Mean_norm"
            if col not in sub_wide.columns:
                n_subcortical_col_missing += 1
                continue
            measure = "volume_norm"
            value_source = sub_wide[
                ["Subject", "Timepoint", col]
            ].rename(columns={col: "Value"})
        elif tissue == "cortical":
            fs_region = entry["fs_region"]
            region_data = cortical_bilateral[
                cortical_bilateral["StructName"] == fs_region
            ]
            if region_data.empty:
                n_cortical_region_missing += 1
                continue
            measure = "thickness_mm"
            value_source = region_data[
                ["Subject", "Timepoint", "ThickAvg_bilateral"]
            ].rename(columns={"ThickAvg_bilateral": "Value"})
        else:
            raise ValueError(
                f"Unknown tissue type {tissue!r} in ROI_PANEL_THESIS "
                f"entry for {roi_display!r}"
            )

        for subj in complete_ids:
            subj_data = value_source[value_source["Subject"] == subj]
            if subj_data.empty:
                continue
            tp_vals = subj_data.set_index("Timepoint")["Value"]

            for tp in TIMEPOINTS:
                if tp not in tp_vals.index:
                    continue
                value = tp_vals[tp]
                if pd.isna(value):
                    n_value_nan += 1
                    continue

                # Actual elapsed time per subject.
                years = years_lookup.get((subj, tp))
                if years is None:
                    n_years_missing += 1
                    continue

                # Per-subject covariates.
                if (subj, tp) not in meta_slim.index:
                    n_meta_missing += 1
                    continue
                meta_row = meta_slim.loc[(subj, tp)]
                # If duplicates remain, .loc returns a DataFrame.
                if isinstance(meta_row, pd.DataFrame):
                    meta_row = meta_row.iloc[0]

                rows.append({
                    "Subject":              subj,
                    "Timepoint":            tp,
                    "Years_from_BL":        float(years),
                    "ROI":                  roi_display,
                    "Tissue":               tissue,
                    "Measure":              measure,
                    "Value":                float(value),
                    "Group_clinical":       meta_row["Group_clinical"],
                    "age_at_BL":            float(meta_row["age_at_BL"]),
                    "PatientSex_clinical":  meta_row["PatientSex_clinical"],
                    "BatchID":              meta_row["BatchID"],
                })

    # Diagnostics.
    if n_subcortical_col_missing:
        print(f"  {n_subcortical_col_missing} subcortical ROI(s) "
              f"missing from cohort_wide.")
    if n_cortical_region_missing:
        print(f"  {n_cortical_region_missing} cortical ROI(s) missing "
              f"from cohort_cortical_regions.")
    if n_value_nan:
        print(f"  {n_value_nan} obs dropped (NaN value).")
    if n_years_missing:
        print(f"  {n_years_missing} obs dropped (no DICOM/SAV anchor).")
    if n_meta_missing:
        print(f"  {n_meta_missing} obs dropped (missing covariates).")

    out = pd.DataFrame(rows)

    if out.empty:
        print("  Long-format output is empty.")
        return out

    # Enforce no NaN in Value (already filtered; sanity check).
    n_nan_values = out["Value"].isna().sum()
    if n_nan_values:
        raise ValueError(
            f"build_cohort_long: {n_nan_values} NaN values "
            f"in Value column."
        )

    # Stable column order.
    col_order = [
        "Subject", "Timepoint", "Years_from_BL",
        "ROI", "Tissue", "Measure", "Value",
        "Group_clinical", "age_at_BL", "PatientSex_clinical",
        "BatchID",
    ]
    return out[col_order].sort_values(
        ["Subject", "Timepoint", "Tissue", "ROI"]
    ).reset_index(drop=True)


# --- ROI-level aggregation ---

def compute_roi_level_summary(deltas_df: pd.DataFrame) -> pd.DataFrame:
    """ROI-level mean/SD summary of subject-level deltas.

    One row per (ROI, Tissue, Cluster, Measure, Delta_Window), with
    mean and SD of Delta_abs and Delta_pct plus n distinct subjects.
    Stage B QC only; the inferential analysis is in 07_analysis.py.

    Returns an empty DataFrame if the input is empty.
    """
    if deltas_df.empty:
        return pd.DataFrame()

    summary = (
        deltas_df
        .groupby(["ROI", "Tissue", "Cluster", "Measure", "Delta_Window"])
        .agg(
            mean_delta_abs=("Delta_abs", "mean"),
            sd_delta_abs=("Delta_abs",   "std"),
            mean_delta_pct=("Delta_pct", "mean"),
            sd_delta_pct=("Delta_pct",   "std"),
            n=("Subject", "nunique"),
        )
        .reset_index()
    )

    return summary


# --- Per-subject elapsed time ---

# Dates before this are treated as sentinels. ParkWest started enrolling
# around 2004.
_PLAUSIBLE_DATE_MIN = pd.Timestamp("2003-01-01")

# Tolerance (years) on either side of a window's nominal length for
# accepting a parsed interval. 1.5y covers the real cohort range
# (BL-3Y: ~2.5-3.7y, BL-5Y: ~3.8-5.4y).
_PLAUSIBLE_INTERVAL_TOLERANCE_YEARS = 1.5


def _parse_study_date(v) -> Optional[pd.Timestamp]:
    """Parse a DICOM YYYYMMDD StudyDate to a Timestamp.

    Accepts string or float input (pandas auto-casts "20081112" to
    20081112.0). Returns None for NaN/None, dates before
    _PLAUSIBLE_DATE_MIN (catches the 1900-01-01 sentinel), future dates,
    and anything that doesn't fit the 8-digit pattern.
    """
    if pd.isna(v):
        return None
    try:
        s = str(int(float(v)))
    except (ValueError, TypeError):
        return None
    if len(s) != 8:
        return None
    ts = pd.to_datetime(s, format="%Y%m%d", errors="coerce")
    if pd.isna(ts):
        return None
    if ts < _PLAUSIBLE_DATE_MIN or ts > pd.Timestamp.today():
        return None
    return ts


def _parse_iso_date(v) -> Optional[pd.Timestamp]:
    """Parse a YYYY-MM-DD date to a Timestamp.

    Used for bl_mri_date, which comes from the SPSS file via
    03_clinical_covariates.py. Returns None for missing or unparseable
    values and for dates outside the plausible range.
    """
    if pd.isna(v):
        return None
    s = str(v).strip()
    if not s or s.lower() in ("nan", "nat", "none"):
        return None
    ts = pd.to_datetime(s, format="%Y-%m-%d", errors="coerce")
    if pd.isna(ts):
        return None
    if ts < _PLAUSIBLE_DATE_MIN or ts > pd.Timestamp.today():
        return None
    return ts


def compute_interval_years(
    scanner_meta: pd.DataFrame,
    delta_windows: Optional[List[Tuple[str, str]]] = None,
) -> pd.DataFrame:
    """Per-subject elapsed time in years between BL and follow-up timepoints.

    scanner_meta needs Subject, Timepoint, and StudyDate; usually read
    from scanner_metadata_with_covariates.csv. delta_windows defaults
    to DELTA_WINDOWS from constants.

    interval_source records where the dates came from:
        'dicom'             both dates from real DICOM headers.
        'sav_bl_mri_date'   BL DICOM is missing but the SPSS bl_mri_date
                            column has it. Expected for SUS BL, where
                            the DICOM headers were de-identified.
        'fallback_nominal'  no usable date; falls back to TIMEPOINT_YEARS.

    Returns one row per (Subject, Delta_Window) with columns Subject,
    Delta_Window, bl_date, t2_date, interval_years, interval_source.
    Intervals outside [nominal +/- _PLAUSIBLE_INTERVAL_TOLERANCE_YEARS]
    are warned about.
    """
    if delta_windows is None:
        delta_windows = DELTA_WINDOWS

    needed_cols = ["Subject", "Timepoint", "StudyDate"]
    missing = [c for c in needed_cols if c not in scanner_meta.columns]
    if missing:
        raise ValueError(
            f"compute_interval_years: scanner_meta is missing required "
            f"column(s): {missing}. Expected {needed_cols}."
        )

    # bl_mri_date is optional; used as a BL fallback if DICOM is sentinel.
    has_bl_mri_date = "bl_mri_date" in scanner_meta.columns
    if not has_bl_mri_date:
        print("  scanner_meta has no 'bl_mri_date'; BL rows with "
              "sentinel DICOM dates will use nominal intervals.")

    # Pivot StudyDate for (Subject, Timepoint) lookups. drop_duplicates
    # is defensive; 02_scanner_metadata.py only emits one row per pair.
    pivot = (
        scanner_meta[needed_cols]
        .drop_duplicates(subset=["Subject", "Timepoint"], keep="first")
        .set_index(["Subject", "Timepoint"])["StudyDate"]
        .unstack("Timepoint")
    )

    # Per-subject BL date from SAV (one value per subject).
    if has_bl_mri_date:
        sav_bl_per_subject = (
            scanner_meta[["Subject", "bl_mri_date"]]
            .drop_duplicates(subset=["Subject"], keep="first")
            .set_index("Subject")["bl_mri_date"]
            .to_dict()
        )
    else:
        sav_bl_per_subject = {}

    rows: List[Dict] = []
    plausibility_warnings: List[str] = []
    source_counts: Dict[str, int] = {
        "dicom": 0, "sav_bl_mri_date": 0, "fallback_nominal": 0,
    }

    for subject in pivot.index:
        # Parse all available DICOM dates once per subject.
        dicom_dates = {tp: _parse_study_date(pivot.loc[subject].get(tp))
                       for tp in pivot.columns}
        sav_bl = _parse_iso_date(sav_bl_per_subject.get(subject))

        for tp1, tp2 in delta_windows:
            d1 = dicom_dates.get(tp1)
            d2 = dicom_dates.get(tp2)

            if d1 is not None and d2 is not None:
                # Both DICOM dates present.
                bl_date = d1
                t2_date = d2
                interval = (t2_date - bl_date).days / 365.25
                source = "dicom"
            elif (d1 is None and d2 is not None
                    and tp1 == "BL" and sav_bl is not None):
                # BL DICOM missing but SAV has it (SUS BL case).
                bl_date = sav_bl
                t2_date = d2
                interval = (t2_date - bl_date).days / 365.25
                source = "sav_bl_mri_date"
            else:
                # No usable anchor.
                bl_date = None
                t2_date = None
                interval = TIMEPOINT_YEARS[tp2] - TIMEPOINT_YEARS[tp1]
                source = "fallback_nominal"

            source_counts[source] += 1

            # Window-aware plausibility check.
            nominal = TIMEPOINT_YEARS[tp2] - TIMEPOINT_YEARS[tp1]
            plaus_min = nominal - _PLAUSIBLE_INTERVAL_TOLERANCE_YEARS
            plaus_max = nominal + _PLAUSIBLE_INTERVAL_TOLERANCE_YEARS
            if not (plaus_min <= interval <= plaus_max):
                plausibility_warnings.append(
                    f"  {subject}  {tp1}→{tp2}  "
                    f"interval={interval:.3f}y  source={source}  "
                    f"(plausible range for this window: "
                    f"[{plaus_min:.1f}, {plaus_max:.1f}])"
                )

            rows.append({
                "Subject":         subject,
                "Delta_Window":    f"{tp1}→{tp2}",
                "bl_date":         (bl_date.strftime("%Y-%m-%d")
                                    if bl_date is not None else None),
                "t2_date":         (t2_date.strftime("%Y-%m-%d")
                                    if t2_date is not None else None),
                "interval_years":  float(interval),
                "interval_source": source,
            })

    out = pd.DataFrame(rows)

    print(f"  compute_interval_years: {len(out)} rows from "
          f"{pivot.shape[0]} subjects.")
    print(f"    sources: dicom={source_counts['dicom']}, "
          f"sav_bl_mri_date={source_counts['sav_bl_mri_date']}, "
          f"fallback_nominal={source_counts['fallback_nominal']}.")
    if plausibility_warnings:
        print(
            f"  {len(plausibility_warnings)} interval(s) outside "
            f"+/-{_PLAUSIBLE_INTERVAL_TOLERANCE_YEARS}y of nominal:"
        )
        for w in plausibility_warnings:
            print(w)

    return out
