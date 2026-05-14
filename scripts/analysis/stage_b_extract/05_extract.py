#!/usr/bin/env python3
"""Parse FastSurfer stats files and compute per-subject analysis deltas.

Walks the longitudinal output tree, reads each timepoint's *.stats
files, and writes:

    cohort_wide.csv              one row per (Subject, Timepoint)
    cohort_cortical_regions.csv  one row per (Subject, Timepoint, Hemi, Region)
    subject_roi_deltas.csv       per-subject deltas over the panel
    roi_level_summary.csv        ROI-level QC summary
    skip_log.txt                 rejected/warned subjects (if any)

With --full-output it also writes cohort_long.csv and per-site-group /
per-subject splits, which downstream scripts don't read but help for
manual QC.

Usage:
    python scripts/stage_b_extract/05_extract.py \\
        --input-dir  outputs/fastsurfer_longitudinal \\
        --output-dir outputs/stage_b_extract/

    # With the extra splits and cohort_long.csv:
    ... --full-output

    # Dry run, no files written:
    ... --no-save
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from pipeline_lib.constants import (
    ASEG_DKT_STATS_FILENAME,
    ASEG_STATS_FILENAME,
    CEREBELLUM_STATS_FILENAME,
    LH_APARC_FILENAME,
    ORDERED_METRIC_COLS,
    RH_APARC_FILENAME,
    SITE_GROUP_DIRS,
    TIMEPOINTS,
    TIMEPOINT_YEARS,
    VOLUME_COLS_FOR_NORM,
)
from pipeline_lib.deltas import (
    compute_analysis_deltas,
    compute_roi_level_summary,
)
from pipeline_lib.ids import normalize_subject_id


# --- MetricType classification (for long-format output) ---

def _build_metric_type_map() -> Dict[str, str]:
    """Build MetricType lookup for the long-format 'MetricType' column."""
    m: Dict[str, str] = {}

    for col in [
        "MaskVol", "eTIV",
        "BrainSegVol", "BrainSegVolNotVent", "VentricleChoroidVol",
        "lhCortexVol", "rhCortexVol", "CortexVol",
        "lhCerebralWhiteMatterVol", "rhCerebralWhiteMatterVol",
        "CerebralWhiteMatterVol",
        "SubCortGrayVol", "TotalGrayVol",
        "SupraTentorialVol", "SupraTentorialVolNotVent",
        "Left-Thalamus", "Right-Thalamus", "Thalamus_Mean",
        "Left-Caudate",  "Right-Caudate",  "Caudate_Mean",
        "Left-Putamen",  "Right-Putamen",  "Putamen_Mean",
        "Left-Pallidum", "Right-Pallidum", "Pallidum_Mean",
        "Left-Hippocampus", "Right-Hippocampus", "Hippocampus_Mean",
        "Left-Amygdala",    "Right-Amygdala",    "Amygdala_Mean",
        "Left-Accumbens-area", "Right-Accumbens-area", "Accumbens_Mean",
        "Left-VentralDC",      "Right-VentralDC",
        "Brain-Stem", "WM-hypointensities",
        "Left-Lateral-Ventricle", "Right-Lateral-Ventricle",
        "Left-Inf-Lat-Vent", "Right-Inf-Lat-Vent",
        "3rd-Ventricle", "4th-Ventricle",
        "LateralVentricles_Total", "InfLatVentricles_Total",
        "Ventricles_Total_Main",
        "lh_RegionalGrayVol_Sum", "rh_RegionalGrayVol_Sum",
        "CortexRegionalGrayVol_Sum",
        "Cbm_Left_I_IV", "Cbm_Right_I_IV",
        "Cbm_Left_V", "Cbm_Right_V",
        "Cbm_Anterior_Total",
        "Cbm_Left_CortexVol", "Cbm_Right_CortexVol", "Cbm_CortexVol_Total",
        "Cbm_Left_WhiteMatter", "Cbm_Right_WhiteMatter",
        "Cbm_WhiteMatter_Total",
        "Cbm_Vermis_VI", "Cbm_Vermis_VII", "Cbm_Vermis_VIII",
        "Cbm_Vermis_IX", "Cbm_Vermis_X", "Cbm_Vermis_Total",
    ]:
        m[col] = "raw"

    for col in [
        "lh_NumVert", "rh_NumVert",
        "lh_WhiteSurfArea", "rh_WhiteSurfArea", "Total_WhiteSurfArea",
        "lh_MeanThickness", "rh_MeanThickness", "MeanThickness_Weighted",
    ]:
        m[col] = "surface"

    for col in ["BrainSegVol_to_eTIV", "MaskVol_to_eTIV"]:
        m[col] = "ratio"

    for col in [
        "lh_RegionalGrayVol_vs_lhCortexVol_pctdiff",
        "rh_RegionalGrayVol_vs_rhCortexVol_pctdiff",
        "CortexRegionalGrayVol_vs_CortexVol_pctdiff",
    ]:
        m[col] = "consistency"

    return m


METRIC_TYPE_MAP: Dict[str, str] = _build_metric_type_map()


def metric_type(col_name: str) -> str:
    """Return the MetricType label for a column name."""
    if col_name in METRIC_TYPE_MAP:
        return METRIC_TYPE_MAP[col_name]
    if col_name.endswith("_norm"):
        return "normalized"
    if col_name.startswith("AI_") and col_name.endswith("_pct"):
        return "asymmetry_index"
    return "other"


# --- Math helpers ---

def asymmetry_index(left: float, right: float) -> float:
    """AI (%) = 100 * (L - R) / ((L + R) / 2). Positive -> left larger."""
    if pd.isna(left) or pd.isna(right) or (left + right) == 0:
        return np.nan
    return 100.0 * (left - right) / ((left + right) / 2.0)


def safe_float(x) -> float:
    try:
        return float(x)
    except Exception:
        return np.nan


def nansum_or_nan(values: List[float]) -> float:
    """nansum that returns NaN only when ALL inputs are NaN."""
    arr = np.array(values, dtype=float)
    if np.isnan(arr).all():
        return np.nan
    return float(np.nansum(arr))


def nanmean_or_nan(values: List[float]) -> float:
    """nanmean that returns NaN only when ALL inputs are NaN.

    Used for bilateral aggregation of subcortical volumes (the ROI
    panel uses means, not sums). With one hemisphere missing, the
    result is just the available value.
    """
    arr = np.array(values, dtype=float)
    if np.isnan(arr).all():
        return np.nan
    return float(np.nanmean(arr))


def percent_difference(reference: float, observed: float) -> float:
    """100 * (observed - reference) / reference"""
    if pd.isna(reference) or pd.isna(observed) or reference == 0:
        return np.nan
    return 100.0 * (observed - reference) / reference


def pick_first(d: Dict[str, float], candidates: List[str]) -> float:
    for key in candidates:
        if key in d:
            return d[key]
    return np.nan


# --- File parsers ---

def parse_measure_lines(text: str) -> Dict[str, float]:
    """Parse '# Measure ..., short_name, ..., value, unit' header lines."""
    measures: Dict[str, float] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line.startswith("# Measure"):
            continue
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 5:
            continue
        value = safe_float(parts[3])
        if not pd.isna(value):
            measures[parts[1]] = value
    return measures


def _detect_col_positions(text: str) -> Optional[Tuple[int, int]]:
    """Find Volume_mm3 and StructName column indices from a '# ColHeaders' line."""
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("# ColHeaders"):
            headers = line.split()[2:]
            try:
                return headers.index("Volume_mm3"), headers.index("StructName")
            except ValueError:
                return None
    return None


def parse_structure_table(text: str) -> Dict[str, float]:
    """Parse aseg/aseg+DKT data rows -> {StructName: volume_mm3}."""
    pos = _detect_col_positions(text)
    vol_col, name_col = pos if pos else (3, 4)
    structures: Dict[str, float] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) <= max(vol_col, name_col):
            continue
        volume = safe_float(parts[vol_col])
        if not pd.isna(volume):
            structures[parts[name_col]] = volume
    return structures


def parse_aparc_region_table(text: str) -> pd.DataFrame:
    """Parse lh/rh aparc stats region rows into a DataFrame."""
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 10 or parts[0].lower() == "structname":
            continue
        rows.append({
            "StructName": parts[0],
            "NumVert":  safe_float(parts[1]),
            "SurfArea": safe_float(parts[2]),
            "GrayVol":  safe_float(parts[3]),
            "ThickAvg": safe_float(parts[4]),
            "ThickStd": safe_float(parts[5]),
            "MeanCurv": safe_float(parts[6]),
            "GausCurv": safe_float(parts[7]),
            "FoldInd":  safe_float(parts[8]),
            "CurvInd":  safe_float(parts[9]),
        })
    cols = ["StructName", "NumVert", "SurfArea", "GrayVol", "ThickAvg",
            "ThickStd", "MeanCurv", "GausCurv", "FoldInd", "CurvInd"]
    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=cols)


def load_stats_file(path: Path) -> Tuple[Dict[str, float], Dict[str, float]]:
    text = path.read_text(errors="ignore")
    return parse_measure_lines(text), parse_structure_table(text)


def load_aparc_file(path: Path) -> Tuple[Dict[str, float], pd.DataFrame]:
    text = path.read_text(errors="ignore")
    return parse_measure_lines(text), parse_aparc_region_table(text)


def load_cerebellum_file(path: Path) -> Dict[str, float]:
    """Load cerebellum.CerebNet.stats and return {structure: volume}."""
    text = path.read_text(errors="ignore")
    return parse_structure_table(text)


# --- Metric extractors ---

def extract_global_metrics(m: Dict[str, float]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    out["MaskVol"]                  = pick_first(m, ["MaskVol"])
    out["eTIV"]                     = pick_first(m, ["eTIV", "EstimatedTotalIntraCranialVol"])
    out["BrainSegVol"]              = pick_first(m, ["BrainSegVol", "BrainSeg"])
    out["BrainSegVolNotVent"]       = pick_first(m, ["BrainSegVolNotVent", "BrainSegNotVent"])
    out["VentricleChoroidVol"]      = pick_first(m, ["VentricleChoroidVol"])
    out["lhCortexVol"]              = pick_first(m, ["lhCortexVol"])
    out["rhCortexVol"]              = pick_first(m, ["rhCortexVol"])
    out["CortexVol"]                = pick_first(m, ["CortexVol"])
    out["lhCerebralWhiteMatterVol"] = pick_first(m, ["lhCerebralWhiteMatterVol"])
    out["rhCerebralWhiteMatterVol"] = pick_first(m, ["rhCerebralWhiteMatterVol"])
    out["CerebralWhiteMatterVol"]   = pick_first(m, ["CerebralWhiteMatterVol"])
    out["SubCortGrayVol"]           = pick_first(m, ["SubCortGrayVol"])
    out["TotalGrayVol"]             = pick_first(m, ["TotalGrayVol"])
    out["SupraTentorialVol"]        = pick_first(m, ["SupraTentorialVol", "SupraTentorial"])
    out["SupraTentorialVolNotVent"] = pick_first(m, ["SupraTentorialVolNotVent"])
    out["BrainSegVol_to_eTIV"]      = pick_first(m, ["BrainSegVol-to-eTIV"])
    out["MaskVol_to_eTIV"]          = pick_first(m, ["MaskVol-to-eTIV"])
    out["AI_CortexVol_pct"]         = asymmetry_index(out["lhCortexVol"], out["rhCortexVol"])
    return out


def extract_subcortical_metrics(s: Dict[str, float]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for prefix in ["Left", "Right"]:
        out[f"{prefix}-Thalamus"]          = pick_first(s, [f"{prefix}-Thalamus", f"{prefix}-Thalamus-Proper"])
        out[f"{prefix}-Caudate"]           = pick_first(s, [f"{prefix}-Caudate"])
        out[f"{prefix}-Putamen"]           = pick_first(s, [f"{prefix}-Putamen"])
        out[f"{prefix}-Pallidum"]          = pick_first(s, [f"{prefix}-Pallidum"])
        out[f"{prefix}-Hippocampus"]       = pick_first(s, [f"{prefix}-Hippocampus"])
        out[f"{prefix}-Amygdala"]          = pick_first(s, [f"{prefix}-Amygdala"])
        out[f"{prefix}-Accumbens-area"]    = pick_first(s, [f"{prefix}-Accumbens-area"])
        out[f"{prefix}-VentralDC"]         = pick_first(s, [f"{prefix}-VentralDC"])
        out[f"{prefix}-Lateral-Ventricle"] = pick_first(s, [f"{prefix}-Lateral-Ventricle"])
        out[f"{prefix}-Inf-Lat-Vent"]      = pick_first(s, [f"{prefix}-Inf-Lat-Vent"])

    for name in ["Thalamus", "Caudate", "Putamen", "Pallidum", "Hippocampus", "Amygdala"]:
        out[f"{name}_Mean"] = nanmean_or_nan([out[f"Left-{name}"], out[f"Right-{name}"]])

    out["Accumbens_Mean"] = nanmean_or_nan([
        out["Left-Accumbens-area"], out["Right-Accumbens-area"],
    ])

    out["3rd-Ventricle"]      = pick_first(s, ["3rd-Ventricle"])
    out["4th-Ventricle"]      = pick_first(s, ["4th-Ventricle"])
    out["Brain-Stem"]         = pick_first(s, ["Brain-Stem"])
    out["WM-hypointensities"] = pick_first(s, ["WM-hypointensities"])

    out["LateralVentricles_Total"] = nansum_or_nan([
        out["Left-Lateral-Ventricle"], out["Right-Lateral-Ventricle"],
    ])
    out["InfLatVentricles_Total"] = nansum_or_nan([
        out["Left-Inf-Lat-Vent"], out["Right-Inf-Lat-Vent"],
    ])
    out["Ventricles_Total_Main"] = nansum_or_nan([
        out["Left-Lateral-Ventricle"], out["Right-Lateral-Ventricle"],
        out["Left-Inf-Lat-Vent"],      out["Right-Inf-Lat-Vent"],
        out["3rd-Ventricle"],          out["4th-Ventricle"],
    ])

    for name in ["Thalamus", "Caudate", "Putamen", "Pallidum", "Hippocampus", "Amygdala"]:
        out[f"AI_{name}_pct"] = asymmetry_index(out[f"Left-{name}"], out[f"Right-{name}"])

    return out


def extract_cortical_summary(lh_m: Dict[str, float], rh_m: Dict[str, float]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    out["lh_NumVert"]          = pick_first(lh_m, ["NumVert"])
    out["rh_NumVert"]          = pick_first(rh_m, ["NumVert"])
    out["lh_WhiteSurfArea"]    = pick_first(lh_m, ["WhiteSurfArea"])
    out["rh_WhiteSurfArea"]    = pick_first(rh_m, ["WhiteSurfArea"])
    out["Total_WhiteSurfArea"] = nansum_or_nan([out["lh_WhiteSurfArea"], out["rh_WhiteSurfArea"]])
    out["lh_MeanThickness"]    = pick_first(lh_m, ["MeanThickness"])
    out["rh_MeanThickness"]    = pick_first(rh_m, ["MeanThickness"])

    lh_t, rh_t = out["lh_MeanThickness"], out["rh_MeanThickness"]
    lh_a, rh_a = out["lh_WhiteSurfArea"],  out["rh_WhiteSurfArea"]

    if all(pd.notna(v) for v in [lh_t, rh_t, lh_a, rh_a]) and (lh_a + rh_a) > 0:
        out["MeanThickness_Weighted"] = (lh_t * lh_a + rh_t * rh_a) / (lh_a + rh_a)
    else:
        out["MeanThickness_Weighted"] = np.nan

    out["AI_WhiteSurfArea_pct"] = asymmetry_index(lh_a, rh_a)
    out["AI_MeanThickness_pct"] = asymmetry_index(lh_t, rh_t)
    return out


def extract_cortical_regional_consistency(
    lh_df: pd.DataFrame, rh_df: pd.DataFrame,
    lh_cortex_vol: float, rh_cortex_vol: float,
) -> Dict[str, float]:
    out: Dict[str, float] = {}

    # sum(min_count=1) so an all-NaN GrayVol column returns NaN, not 0.0.
    out["lh_RegionalGrayVol_Sum"] = (
        float(lh_df["GrayVol"].sum(min_count=1)) if not lh_df.empty else np.nan
    )
    out["rh_RegionalGrayVol_Sum"] = (
        float(rh_df["GrayVol"].sum(min_count=1)) if not rh_df.empty else np.nan
    )
    out["CortexRegionalGrayVol_Sum"] = nansum_or_nan([
        out["lh_RegionalGrayVol_Sum"], out["rh_RegionalGrayVol_Sum"],
    ])
    out["AI_RegionalGrayVol_pct"] = asymmetry_index(
        out["lh_RegionalGrayVol_Sum"], out["rh_RegionalGrayVol_Sum"],
    )
    out["lh_RegionalGrayVol_vs_lhCortexVol_pctdiff"] = percent_difference(
        lh_cortex_vol, out["lh_RegionalGrayVol_Sum"],
    )
    out["rh_RegionalGrayVol_vs_rhCortexVol_pctdiff"] = percent_difference(
        rh_cortex_vol, out["rh_RegionalGrayVol_Sum"],
    )
    out["CortexRegionalGrayVol_vs_CortexVol_pctdiff"] = percent_difference(
        nansum_or_nan([lh_cortex_vol, rh_cortex_vol]),
        out["CortexRegionalGrayVol_Sum"],
    )
    return out


def extract_cerebellar_metrics(s: Dict[str, float]) -> Dict[str, float]:
    """Extract cerebellar sub-region volumes from cerebellum.CerebNet.stats.

    Pulls anterior lobe (I-IV + V) per hemisphere and bilateral total,
    vermis sub-regions and total, and total cerebellar cortex and
    white matter per hemisphere.
    """
    out: Dict[str, float] = {}

    # Anterior lobe = lobules I-IV + V.
    out["Cbm_Left_I_IV"]  = pick_first(s, ["Cbm_Left_I_IV"])
    out["Cbm_Right_I_IV"] = pick_first(s, ["Cbm_Right_I_IV"])
    out["Cbm_Left_V"]     = pick_first(s, ["Cbm_Left_V"])
    out["Cbm_Right_V"]    = pick_first(s, ["Cbm_Right_V"])
    out["Cbm_Anterior_Total"] = nansum_or_nan([
        out["Cbm_Left_I_IV"], out["Cbm_Right_I_IV"],
        out["Cbm_Left_V"],    out["Cbm_Right_V"],
    ])

    out["Cbm_Left_CortexVol"]  = pick_first(s, ["Left-Cerebellum-Cortex"])
    out["Cbm_Right_CortexVol"] = pick_first(s, ["Right-Cerebellum-Cortex"])
    out["Cbm_CortexVol_Total"] = nansum_or_nan([
        out["Cbm_Left_CortexVol"], out["Cbm_Right_CortexVol"],
    ])

    out["Cbm_Left_WhiteMatter"]  = pick_first(s, ["Left-Cerebellum-White-Matter"])
    out["Cbm_Right_WhiteMatter"] = pick_first(s, ["Right-Cerebellum-White-Matter"])
    out["Cbm_WhiteMatter_Total"] = nansum_or_nan([
        out["Cbm_Left_WhiteMatter"], out["Cbm_Right_WhiteMatter"],
    ])

    for lobule in ["VI", "VII", "VIII", "IX", "X"]:
        out[f"Cbm_Vermis_{lobule}"] = pick_first(s, [f"Cbm_Vermis_{lobule}"])

    vermis_summary = pick_first(s, ["Cbm_Vermis"])
    if pd.notna(vermis_summary):
        out["Cbm_Vermis_Total"] = vermis_summary
    else:
        out["Cbm_Vermis_Total"] = nansum_or_nan([
            out[f"Cbm_Vermis_{lb}"] for lb in ["VI", "VII", "VIII", "IX", "X"]
        ])

    return out


# --- Subject discovery ---

def stats_files_present(tp_dir: Path) -> bool:
    """True only if all four required stats files exist."""
    stats = tp_dir / "stats"
    return all((stats / f).exists() for f in [
        ASEG_STATS_FILENAME, ASEG_DKT_STATS_FILENAME,
        LH_APARC_FILENAME, RH_APARC_FILENAME,
    ])


def discover_subjects(
    input_root: Path,
) -> Tuple[List[Tuple[str, str, str, Path]], List[str]]:
    """Walk input_root for every site/group in SITE_GROUP_DIRS.

    Accepts subjects with at least one complete timepoint. Stricter
    cohort-level acceptance lives in cohort.csv (stage A, step 00).

    Returns (valid, skipped):
      valid    [(site, group, subject_id, subject_path), ...]
      skipped  human-readable messages for rejected or warned subjects.
    """
    valid: List[Tuple[str, str, str, Path]] = []
    skipped: List[str] = []

    for dir_name, (site, group) in SITE_GROUP_DIRS.items():
        site_dir = input_root / dir_name
        if not site_dir.is_dir():
            skipped.append(f"SITE/GROUP NOT FOUND: {site_dir}")
            continue

        for subj_dir in sorted(p for p in site_dir.iterdir() if p.is_dir()):
            if "_template" in subj_dir.name.lower():
                continue

            subject_id = normalize_subject_id(subj_dir.name)
            available  = [
                tp for tp in TIMEPOINTS
                if stats_files_present(subj_dir / tp)
            ]

            if not available:
                skipped.append(
                    f"SKIP  {site}/{group}/{subject_id}: no complete "
                    f"timepoint folder found (checked {list(TIMEPOINTS)})"
                )
                continue

            missing = [tp for tp in TIMEPOINTS if tp not in available]
            if missing:
                skipped.append(
                    f"WARN  {site}/{group}/{subject_id}: missing "
                    f"timepoints {missing} -- NaN rows produced."
                )

            valid.append((site, group, subject_id, subj_dir))

    return valid, skipped


# --- Per-subject extraction ---

def extract_subject_timepoint(
    site: str,
    group: str,
    subject_id: str,
    subj_dir: Path,
    timepoint: str,
) -> Tuple[Optional[Dict], Optional[pd.DataFrame]]:
    """Extract all metrics for one subject x timepoint.

    Returns (wide_row, cortical_region_df), or (None, None) if files
    are missing.
    """
    stats_dir = subj_dir / timepoint / "stats"

    aseg_path     = stats_dir / ASEG_STATS_FILENAME
    aseg_dkt_path = stats_dir / ASEG_DKT_STATS_FILENAME
    lh_path       = stats_dir / LH_APARC_FILENAME
    rh_path       = stats_dir / RH_APARC_FILENAME

    if not all(p.exists() for p in [aseg_path, aseg_dkt_path, lh_path, rh_path]):
        return None, None

    aseg_measures, _          = load_stats_file(aseg_path)
    _, aseg_dkt_structures    = load_stats_file(aseg_dkt_path)
    lh_measures, lh_region_df = load_aparc_file(lh_path)
    rh_measures, rh_region_df = load_aparc_file(rh_path)

    global_m      = extract_global_metrics(aseg_measures)
    subcortical_m = extract_subcortical_metrics(aseg_dkt_structures)
    cortical_m    = extract_cortical_summary(lh_measures, rh_measures)
    consistency_m = extract_cortical_regional_consistency(
        lh_region_df, rh_region_df,
        global_m["lhCortexVol"], global_m["rhCortexVol"],
    )

    # Cerebellum is optional.
    cerebellum_path = stats_dir / CEREBELLUM_STATS_FILENAME
    if cerebellum_path.exists():
        cbm_structures = load_cerebellum_file(cerebellum_path)
        cerebellar_m   = extract_cerebellar_metrics(cbm_structures)
    else:
        cerebellar_m = {k: np.nan for k in [
            "Cbm_Left_I_IV", "Cbm_Right_I_IV", "Cbm_Left_V", "Cbm_Right_V",
            "Cbm_Anterior_Total",
            "Cbm_Left_CortexVol", "Cbm_Right_CortexVol", "Cbm_CortexVol_Total",
            "Cbm_Left_WhiteMatter", "Cbm_Right_WhiteMatter", "Cbm_WhiteMatter_Total",
            "Cbm_Vermis_VI", "Cbm_Vermis_VII", "Cbm_Vermis_VIII",
            "Cbm_Vermis_IX", "Cbm_Vermis_X", "Cbm_Vermis_Total",
        ]}

    row: Dict = {
        "Site":          site,
        "Group":         group,
        "Subject":       subject_id,
        "Timepoint":     timepoint,
        "Years_from_BL": TIMEPOINT_YEARS[timepoint],
        "Data_Present":  True,
    }
    row.update(global_m)
    row.update(subcortical_m)
    row.update(cortical_m)
    row.update(consistency_m)
    row.update(cerebellar_m)

    etiv = global_m.get("eTIV", np.nan)
    for col in VOLUME_COLS_FOR_NORM:
        raw_val = row.get(col, np.nan)
        row[f"{col}_norm"] = (
            raw_val / etiv
            if pd.notna(raw_val) and pd.notna(etiv) and etiv > 0
            else np.nan
        )

    region_frames = []
    for hemi, hemi_df in [("lh", lh_region_df), ("rh", rh_region_df)]:
        if not hemi_df.empty:
            hemi_copy = hemi_df.copy()
            hemi_copy["GrayVol_norm"] = (
                hemi_copy["GrayVol"] / etiv
                if pd.notna(etiv) and etiv > 0
                else np.nan
            )
            hemi_copy.insert(0, "Hemisphere",    hemi)
            hemi_copy.insert(0, "Years_from_BL", TIMEPOINT_YEARS[timepoint])
            hemi_copy.insert(0, "Timepoint",     timepoint)
            hemi_copy.insert(0, "Subject",       subject_id)
            hemi_copy.insert(0, "Group",         group)
            hemi_copy.insert(0, "Site",          site)
            region_frames.append(hemi_copy)

    cortical_df = (
        pd.concat(region_frames, ignore_index=True)
        if region_frames else pd.DataFrame()
    )
    return row, cortical_df


def extract_subject(
    site: str,
    group: str,
    subject_id: str,
    subj_dir: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Extract all timepoints for one subject.

    Missing timepoints get a Data_Present=False sentinel row.
    """
    wide_rows: List[Dict] = []
    region_frames: List[pd.DataFrame] = []

    for tp in TIMEPOINTS:
        row, cortical_df = extract_subject_timepoint(site, group, subject_id, subj_dir, tp)

        if row is None:
            row = {
                "Site":          site,
                "Group":         group,
                "Subject":       subject_id,
                "Timepoint":     tp,
                "Years_from_BL": TIMEPOINT_YEARS[tp],
                "Data_Present":  False,
            }

        wide_rows.append(row)
        if cortical_df is not None and not cortical_df.empty:
            region_frames.append(cortical_df)

    wide_df = pd.DataFrame(wide_rows)

    id_cols   = ["Site", "Group", "Subject", "Timepoint",
                 "Years_from_BL", "Data_Present"]
    norm_cols = [c for c in wide_df.columns if c.endswith("_norm")]
    ordered   = id_cols + [c for c in ORDERED_METRIC_COLS if c in wide_df.columns] + norm_cols
    extra     = [c for c in wide_df.columns if c not in ordered]
    wide_df   = wide_df[ordered + extra]

    cortical_df = (
        pd.concat(region_frames, ignore_index=True)
        if region_frames else pd.DataFrame()
    )
    return wide_df, cortical_df


# --- Long format conversion ---

def wide_to_long(wide_df: pd.DataFrame) -> pd.DataFrame:
    """Melt the wide DataFrame to long format."""
    id_vars     = ["Site", "Group", "Subject", "Timepoint",
                   "Years_from_BL", "Data_Present"]
    metric_cols = [c for c in wide_df.columns if c not in id_vars]

    long_df = wide_df.melt(
        id_vars=id_vars,
        value_vars=metric_cols,
        var_name="Metric",
        value_name="Value",
    )

    long_df["MetricType"] = long_df["Metric"].map(metric_type)

    long_df = long_df.sort_values(
        ["Site", "Group", "Subject", "MetricType", "Metric", "Years_from_BL"]
    ).reset_index(drop=True)

    col_order = ["Site", "Group", "Subject", "Timepoint", "Years_from_BL",
                 "Data_Present", "Metric", "MetricType", "Value"]
    return long_df[col_order]


# --- Output helpers ---

def save_csv(df: pd.DataFrame, path: Path, index: bool = False) -> None:
    # 10dp preserves precision across the stage-B / stage-C boundary.
    # Subcortical *_norm values are ~1e-3 to 1e-2; 6dp introduced
    # ~1e-6 roundoff in harmonized deltas.
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=index, float_format="%.10f")


def write_skip_log(skipped: List[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(f"Skip/warning log -- {len(skipped)} entries\n")
        f.write("=" * 70 + "\n")
        for msg in skipped:
            f.write(msg + "\n")


def print_summary(wide_df: pd.DataFrame, skipped: List[str]) -> None:
    n_subjects = wide_df.groupby(["Site", "Group", "Subject"]).ngroups
    n_skip = len([s for s in skipped if s.startswith("SKIP")])
    n_warn = len([s for s in skipped if s.startswith("WARN")])

    completeness = (
        wide_df[wide_df["Data_Present"]]
        .groupby(["Site", "Group", "Subject"])["Timepoint"]
        .count()
    )
    n_full    = int((completeness == 3).sum())
    n_partial = int((completeness < 3).sum())

    print(
        f"\nSummary: {n_subjects} subjects processed "
        f"({n_full} complete, {n_partial} partial); "
        f"{n_skip} skipped, {n_warn} warned."
    )

    site_group = (
        wide_df.drop_duplicates(["Site", "Group", "Subject"])
               .groupby(["Site", "Group"]).size()
               .unstack(fill_value=0)
    )
    print(site_group.to_string())
    print()


def print_analysis_summary(
    roi_summary: pd.DataFrame,
    n_complete: int,
) -> None:
    """Print a one-line summary of the delta computation.

    Detailed per-ROI/per-domain summaries are in 07_analysis.py; stage
    B just confirms the deltas were produced.
    """
    n_rois = roi_summary["ROI"].nunique() if not roi_summary.empty else 0
    n_windows = (roi_summary["Delta_Window"].nunique()
                 if not roi_summary.empty else 0)
    print(f"  Computed deltas for {n_complete} subjects across "
          f"{n_rois} ROIs x {n_windows} windows.")


# --- Output writers (split by mode) ---

def write_default_outputs(
    output_root: Path,
    cohort_wide: pd.DataFrame,
    cohort_cortical: pd.DataFrame,
    analysis_deltas: pd.DataFrame,
    roi_summary: pd.DataFrame,
    skipped: List[str],
) -> int:
    """Write the default output set. Returns file count."""
    if skipped:
        write_skip_log(skipped, output_root / "skip_log.txt")

    save_csv(cohort_wide, output_root / "cohort_wide.csv")
    if not cohort_cortical.empty:
        save_csv(cohort_cortical, output_root / "cohort_cortical_regions.csv")

    if not analysis_deltas.empty:
        save_csv(analysis_deltas, output_root / "subject_roi_deltas.csv")
    if not roi_summary.empty:
        save_csv(roi_summary, output_root / "roi_level_summary.csv")

    return sum(1 for _ in output_root.rglob("*.csv")) + (
        1 if (output_root / "skip_log.txt").exists() else 0
    )


def write_full_outputs(
    output_root: Path,
    cohort_wide: pd.DataFrame,
    cohort_long: pd.DataFrame,
    cohort_cortical: pd.DataFrame,
) -> None:
    """Write the --full-output extras: cohort_long.csv, by_site_group/,
    by_subject/."""
    save_csv(cohort_long, output_root / "cohort_long.csv")

    for (site_name, group_name), _ in cohort_wide.groupby(["Site", "Group"]):
        site_dir = output_root / "by_site_group"
        tag = f"{site_name}_{group_name}"
        site_mask = (
            (cohort_wide["Site"] == site_name)
            & (cohort_wide["Group"] == group_name)
        )
        save_csv(cohort_wide[site_mask], site_dir / f"{tag}_wide.csv")
        save_csv(
            cohort_long[
                (cohort_long["Site"] == site_name)
                & (cohort_long["Group"] == group_name)
            ],
            site_dir / f"{tag}_long.csv",
        )
        if not cohort_cortical.empty:
            sub_cortical = cohort_cortical[
                (cohort_cortical["Site"] == site_name)
                & (cohort_cortical["Group"] == group_name)
            ]
            if not sub_cortical.empty:
                save_csv(sub_cortical, site_dir / f"{tag}_cortical_regions.csv")

    for (site_name, group_name, subject), grp in cohort_wide.groupby(
        ["Site", "Group", "Subject"]
    ):
        subj_dir = output_root / "by_subject"
        tag      = f"{site_name}_{group_name}_{subject}"
        save_csv(grp.reset_index(drop=True), subj_dir / f"{tag}_wide.csv")

        if not cohort_cortical.empty:
            subj_cortical = cohort_cortical[
                (cohort_cortical["Site"] == site_name)
                & (cohort_cortical["Group"] == group_name)
                & (cohort_cortical["Subject"] == subject)
            ]
            if not subj_cortical.empty:
                save_csv(subj_cortical, subj_dir / f"{tag}_cortical_regions.csv")


# --- Main ---

def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-dir", type=Path, required=True,
        help="Root directory containing the FastSurfer site folders.",
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Where to write the output CSVs.",
    )
    parser.add_argument(
        "--full-output", action="store_true",
        help="Also write cohort_long.csv plus by_site_group/ and "
             "by_subject/ splits. Off by default; downstream scripts "
             "don't read these but they help for manual QC.",
    )
    parser.add_argument(
        "--no-save", action="store_true",
        help="Print summary only; don't write any files.",
    )
    args = parser.parse_args()

    input_root  = args.input_dir
    output_root = args.output_dir

    if not input_root.is_dir():
        print(f"ERROR: --input-dir does not exist: {input_root}",
              file=sys.stderr)
        return 1

    # Discover subjects.
    print(f"Scanning {input_root} ...")
    valid_subjects, skipped = discover_subjects(input_root)
    print(f"  Found {len(valid_subjects)} valid subjects, "
          f"{len(skipped)} skipped/warned.")

    if not valid_subjects:
        print("No valid subjects found. Check --input-dir and SITE_GROUP_DIRS.")
        return 1

    # Extract all subjects.
    all_wide:     List[pd.DataFrame] = []
    all_cortical: List[pd.DataFrame] = []

    for i, (site, group, subject_id, subj_dir) in enumerate(valid_subjects, 1):
        print(f"  [{i:3d}/{len(valid_subjects)}] {site}/{group}/{subject_id}",
              end="\r", flush=True)
        wide_df, cortical_df = extract_subject(site, group, subject_id, subj_dir)
        all_wide.append(wide_df)
        if not cortical_df.empty:
            all_cortical.append(cortical_df)
    print()

    cohort_wide     = pd.concat(all_wide, ignore_index=True)
    cohort_cortical = (
        pd.concat(all_cortical, ignore_index=True)
        if all_cortical else pd.DataFrame()
    )

    # Summary.
    print_summary(cohort_wide, skipped)

    # Per-subject deltas over the ROI panel.
    print("Computing per-subject deltas over the ROI panel ...")
    analysis_deltas = compute_analysis_deltas(cohort_wide, cohort_cortical)

    if not analysis_deltas.empty:
        roi_summary = compute_roi_level_summary(analysis_deltas)
        n_complete  = analysis_deltas["Subject"].nunique()
        print_analysis_summary(roi_summary, n_complete)
    else:
        roi_summary = pd.DataFrame()
        print("  No complete subjects found; analysis skipped.")

    if args.no_save:
        print("(--no-save: file output skipped)")
        return 0

    # Write outputs.
    print(f"Writing output to {output_root} ...")
    n_files = write_default_outputs(
        output_root,
        cohort_wide, cohort_cortical,
        analysis_deltas, roi_summary,
        skipped,
    )
    if args.full_output:
        print("  --full-output: writing extra splits and cohort_long.csv ...")
        cohort_long = wide_to_long(cohort_wide)
        write_full_outputs(
            output_root, cohort_wide, cohort_long, cohort_cortical,
        )
        n_files = sum(1 for _ in output_root.rglob("*.csv")) + (
            1 if (output_root / "skip_log.txt").exists() else 0
        )

    print(f"Done. {n_files} files written to {output_root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
