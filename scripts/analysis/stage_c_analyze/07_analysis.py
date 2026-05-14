"""LME-based stage C analysis for the ParkWest ROI panel.

Inputs:
    cohort_long_harmonized.csv      from 06_harmonize.py
    cohort_long_unharmonized.csv    sibling file, used for Risk-A

Outputs:
    stats_lme_between_group.csv     12 rows; primary inferential CSV
    stats_lme_within_group.csv      24 rows; descriptive Test A CSV
    roi_display_order.csv           canonical |d_adj_5y| ordering
    lme_sanity_check.txt            Risk-A all-ROI summary
    risk_a_all_rois.csv             Risk-A per-ROI comparison
    lme_diagnostics_test_b.csv      per-ROI Test B model diagnostics
    demographics.csv                cohort summary (9 rows x 3 cols)
    figure1_forest_d_adj_5Y.png             forest plot of d_adj_5y
    figure2_trajectory_fdrsig[_clean].png   FDR-sig trajectories
    figure3_trajectory_subcortical_supp[_clean].png
    figure4_trajectory_cortical_supp[_clean].png
    figure5_heatmap_per_group_atrophy.png   within-group %-change heatmap
    figure_risk_a_scatter_supp.png          harm. vs unharm. scatter
    figure_lme_diagnostics_test_b_supp.png  residual diagnostics

Test B (between-group, 12 fits):
    Value ~ C(Group_clinical, Treatment(reference='Control'))
            * Years_from_BL + age_at_BL + C(PatientSex_clinical)
            + (1 | Subject)
REML; Wald z on the Group x Years_from_BL interaction (no
Kenward-Roger or Satterthwaite, so p-values are approximate).

Test A (within-group, 24 fits):
    Value ~ Years_from_BL + age_at_BL + C(PatientSex_clinical)
            + (1 | Subject)
Descriptive: tests each group's slope per ROI, not PD vs Control.

Effect sizes from Test B:
    predicted_diff_5y_pct = 100 * 5*b / mean_BL_Control
    d_adj_5y = (5*b / sqrt(residual_var)) * J   (Hedges' J)

FDR: Benjamini-Hochberg per tissue (6 subcortical, 6 cortical tests
as separate families).

Risk-A refits Test B on unharmonized data and compares sign and
magnitude vs harmonized. Descriptive, not a hard gate.
"""
from __future__ import annotations

import argparse
import math
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FormatStrFormatter, FuncFormatter, MaxNLocator
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests

# Convergence is checked via mdf.converged; suppress statsmodels' noise.
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


# --- Constants ---

# Patsy's exact label for the Group x Years_from_BL interaction under
# C(Group_clinical, Treatment(reference='Control')) with Group values
# "Control" / "Pasient". Determined by single-fit sandbox.
INTERACTION_COEF_LABEL = (
    "C(Group_clinical, Treatment(reference='Control'))[T.Pasient]"
    ":Years_from_BL"
)

# Group main-effect label (PD vs Control at baseline).
GROUP_MAIN_COEF_LABEL = (
    "C(Group_clinical, Treatment(reference='Control'))[T.Pasient]"
)

# Time main-effect = Control's slope (Control is the reference).
TIME_MAIN_COEF_LABEL = "Years_from_BL"

REFERENCE_D_VALUES = (0.2, 0.5, 0.8)

# Star annotation thresholds (FDR-corrected p-value cutoffs).
STAR_THRESHOLDS = [
    (0.001, "***"),
    (0.010, "**"),
    (0.050, "*"),
]

CSV_FLOAT_FORMAT = "%.12g"


# --- Figure-style constants ---


COLOR_PD = "#D55E00"          # vermillion / orange (Patient)
COLOR_CONTROL = "#0072B2"     # blue (Control)
COLOR_NEUTRAL = "#2C2C2A"     # near-black (forest markers, neutral elements)
COLOR_GRAY_DARK = "#555555"   # axis labels, n.s. text
COLOR_GRAY_LIGHT = "#BBBBBB"  # reference lines, gridlines
COLOR_GRAY_BG = "#E8E8E8"     # subtle backgrounds

# FreeSurferColorLUT-inspired ROI colors for the supplementary Risk-A
# scatter. The figure is split by tissue, so color encodes ROI identity
# within each panel. Dark marker edges keep light colors legible.
RISK_A_SUBCORTICAL_COLORS = {
    "Amygdala": "#67FFFF",      # FreeSurfer LUT: 103, 255, 255
    "Caudate": "#7ABADC",       # FreeSurfer LUT: 122, 186, 220
    "Putamen": "#EC0DB0",       # FreeSurfer LUT: 236, 13, 176
    "Hippocampus": "#DCD814",   # FreeSurfer LUT: 220, 216, 20
    "Thalamus": "#00760E",      # FreeSurfer LUT: 0, 118, 14
    "Accumbens": "#FFA500",     # FreeSurfer LUT: 255, 165, 0
}

RISK_A_CORTICAL_COLORS = {
    "Parahippocampal gyrus": "#14DC3C",        # ctx-*h-parahippocampal
    "Precuneus": "#A08CB4",                    # ctx-*h-precuneus
    "Lingual gyrus": "#E18C8C",                # ctx-*h-lingual
    "Caudal anterior cingulate": "#7D64A0",    # ctx-*h-caudalanteriorcingulate
    "Superior frontal": "#14DCA0",             # ctx-*h-superiorfrontal
    "Entorhinal cortex": "#DC140A",            # ctx-*h-entorhinal
}

RISK_A_SUBCORTICAL_ORDER = [
    "Amygdala", "Caudate", "Putamen",
    "Hippocampus", "Thalamus", "Accumbens",
]

RISK_A_CORTICAL_ORDER = [
    "Parahippocampal gyrus", "Precuneus", "Lingual gyrus",
    "Caudal anterior cingulate", "Superior frontal",
    "Entorhinal cortex",
]

# Figure-level constants.
FIG_DPI = 200
TEXTWIDTH_INCHES = 6.3       

# Per-figure dimensions.
FIG1_FOREST_DIMS = (TEXTWIDTH_INCHES, 5.0)
FIG2_TRAJ_FDRSIG_DIMS = (TEXTWIDTH_INCHES, 5.5)   

# Font sizes.
FONTSIZE_TITLE = 12
FONTSIZE_SUBTITLE = 9
FONTSIZE_LABEL = 11
FONTSIZE_TICK = 10
FONTSIZE_ANNOTATION = 9
FONTSIZE_CAPTION = 8


# --- Loaders ---

def load_long_harmonized(path: Path) -> pd.DataFrame:
    """Load cohort_long_harmonized.csv and validate its schema.

    Checks: required columns present, (Subject, Timepoint, ROI)
    uniqueness, no NaN Value, at least one BL row per subject.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"cohort_long_harmonized.csv not found at {path}"
        )
    df = pd.read_csv(path)
    print(f"Loaded {path.name}: {df.shape}")

    required = {
        "Subject", "Timepoint", "Years_from_BL", "ROI", "Tissue",
        "Measure", "Value", "Group_clinical", "age_at_BL",
        "PatientSex_clinical", "BatchID",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"cohort_long is missing columns: {sorted(missing)}"
        )

    n_dup = df.duplicated(subset=["Subject", "Timepoint", "ROI"]).sum()
    if n_dup:
        raise ValueError(
            f"cohort_long has {n_dup} duplicate (Subject, Timepoint, "
            f"ROI) rows"
        )

    n_nan = df["Value"].isna().sum()
    if n_nan:
        raise ValueError(
            f"cohort_long has {n_nan} NaN Value entries"
        )

    n_subjects = df["Subject"].nunique()
    n_rois = df["ROI"].nunique()
    print(f"  {n_subjects} subjects, {n_rois} ROIs, "
          f"{df['Timepoint'].nunique()} timepoints, "
          f"{df['Group_clinical'].value_counts().to_dict()} per group")

    return df


# --- LME fitting helpers ---

def fit_between_group_lme(roi_data: pd.DataFrame) -> Dict[str, Any]:
    """Fit the between-group LME for one ROI.

    Model:
        Value ~ C(Group_clinical, Treatment(reference='Control'))
                * Years_from_BL
                + age_at_BL
                + C(PatientSex_clinical)
                + (1 | Subject)

    REML; Wald z on the Group x Years_from_BL interaction.

    Returns a dict carrying all fitted coefficients/SEs/p-values
    plus n_obs, n_subjects per group, residual_var, df_resid,
    convergence flag, and per-group BL means (for predicted_diff_5y).
    Status is "ok", "singular_fit", or "convergence_failed:<reason>".
    """
    formula = (
        "Value ~ C(Group_clinical, Treatment(reference='Control'))"
        " * Years_from_BL"
        " + age_at_BL"
        " + C(PatientSex_clinical)"
    )

    n_obs = len(roi_data)
    n_subjects = roi_data["Subject"].nunique()
    n_subjects_pd = roi_data[
        roi_data["Group_clinical"] == "Pasient"
    ]["Subject"].nunique()
    n_subjects_co = roi_data[
        roi_data["Group_clinical"] == "Control"
    ]["Subject"].nunique()

    # Mean BL value per group (for the predicted_diff_5y denominator).
    bl_data = roi_data[roi_data["Timepoint"] == "BL"]
    mean_bl_control = float(
        bl_data[bl_data["Group_clinical"] == "Control"]["Value"].mean()
    )
    mean_bl_pd = float(
        bl_data[bl_data["Group_clinical"] == "Pasient"]["Value"].mean()
    )

    nan_row: Dict[str, Any] = {
        "status": "convergence_failed",
        "n_obs": n_obs,
        "n_subjects": n_subjects,
        "n_subjects_pd": n_subjects_pd,
        "n_subjects_co": n_subjects_co,
        "intercept": np.nan, "intercept_se": np.nan,
        "group_main_coef": np.nan, "group_main_se": np.nan,
        "group_main_p": np.nan,
        "time_main_coef": np.nan, "time_main_se": np.nan,
        "time_main_p": np.nan,
        "interaction_coef": np.nan, "interaction_se": np.nan,
        "interaction_z": np.nan, "interaction_p": np.nan,
        "interaction_ci_low": np.nan, "interaction_ci_high": np.nan,
        "age_coef": np.nan, "age_se": np.nan,
        "sex_coef": np.nan, "sex_se": np.nan,
        "residual_var": np.nan, "random_intercept_var": np.nan,
        "df_resid": 0, "converged": False,
        "mean_bl_control": mean_bl_control,
        "mean_bl_pd": mean_bl_pd,
    }

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            md = smf.mixedlm(formula, data=roi_data,
                             groups=roi_data["Subject"])
            mdf = md.fit(reml=True, method="lbfgs")
    except (np.linalg.LinAlgError, ValueError, OverflowError) as exc:
        nan_row["status"] = f"convergence_failed: {type(exc).__name__}"
        return nan_row

    if not mdf.converged:
        nan_row["status"] = "convergence_failed"
        nan_row["converged"] = False
        return nan_row

    # Detect singular random-effect fit (random-intercept variance
    # at the boundary).
    rand_var = float(mdf.cov_re.iloc[0, 0])
    if rand_var < 1e-20:
        # Effectively a fixed-effects-only fit; flag but don't fail.
        status = "singular_fit"
    else:
        status = "ok"

    # Extract parameters by literal label for robustness across
    # statsmodels versions.
    fe_params = mdf.fe_params
    fe_se = mdf.bse_fe
    fe_pvalues = mdf.pvalues.loc[fe_params.index]
    fe_tvalues = mdf.tvalues.loc[fe_params.index]
    ci = mdf.conf_int().loc[fe_params.index]

    return {
        "status": status,
        "n_obs": n_obs,
        "n_subjects": n_subjects,
        "n_subjects_pd": n_subjects_pd,
        "n_subjects_co": n_subjects_co,
        "intercept": float(fe_params["Intercept"]),
        "intercept_se": float(fe_se["Intercept"]),
        "group_main_coef": float(fe_params[GROUP_MAIN_COEF_LABEL]),
        "group_main_se": float(fe_se[GROUP_MAIN_COEF_LABEL]),
        "group_main_p": float(fe_pvalues[GROUP_MAIN_COEF_LABEL]),
        "time_main_coef": float(fe_params[TIME_MAIN_COEF_LABEL]),
        "time_main_se": float(fe_se[TIME_MAIN_COEF_LABEL]),
        "time_main_p": float(fe_pvalues[TIME_MAIN_COEF_LABEL]),
        "interaction_coef": float(fe_params[INTERACTION_COEF_LABEL]),
        "interaction_se": float(fe_se[INTERACTION_COEF_LABEL]),
        "interaction_z": float(fe_tvalues[INTERACTION_COEF_LABEL]),
        "interaction_p": float(fe_pvalues[INTERACTION_COEF_LABEL]),
        "interaction_ci_low": float(ci.loc[INTERACTION_COEF_LABEL, 0]),
        "interaction_ci_high": float(ci.loc[INTERACTION_COEF_LABEL, 1]),
        "age_coef": float(fe_params["age_at_BL"]),
        "age_se": float(fe_se["age_at_BL"]),
        "sex_coef": float(fe_params["C(PatientSex_clinical)[T.M]"]),
        "sex_se": float(fe_se["C(PatientSex_clinical)[T.M]"]),
        "residual_var": float(mdf.scale),
        "random_intercept_var": rand_var,
        "df_resid": int(n_obs - len(fe_params)),
        "converged": bool(mdf.converged),
        "mean_bl_control": mean_bl_control,
        "mean_bl_pd": mean_bl_pd,
    }


def fit_within_group_lme(
    roi_data: pd.DataFrame, group: str
) -> Dict[str, Any]:
    """Fit the within-group LME for one ROI on a single group's data.

    Model:
        Value ~ Years_from_BL + age_at_BL
              + C(PatientSex_clinical) + (1 | Subject)
    """
    sub = roi_data[roi_data["Group_clinical"] == group]
    n_obs = len(sub)
    n_subjects = sub["Subject"].nunique()

    nan_row: Dict[str, Any] = {
        "status": "convergence_failed",
        "n_obs": n_obs, "n_subjects": n_subjects,
        "intercept": np.nan, "intercept_se": np.nan,
        "slope": np.nan, "slope_se": np.nan,
        "slope_z": np.nan, "slope_p": np.nan,
        "slope_ci_low": np.nan, "slope_ci_high": np.nan,
        "residual_var": np.nan, "random_intercept_var": np.nan,
        "converged": False,
    }

    formula = (
        "Value ~ Years_from_BL + age_at_BL"
        " + C(PatientSex_clinical)"
    )

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            md = smf.mixedlm(formula, data=sub, groups=sub["Subject"])
            mdf = md.fit(reml=True, method="lbfgs")
    except (np.linalg.LinAlgError, ValueError, OverflowError) as exc:
        nan_row["status"] = f"convergence_failed: {type(exc).__name__}"
        return nan_row

    if not mdf.converged:
        return nan_row

    rand_var = float(mdf.cov_re.iloc[0, 0])
    status = "singular_fit" if rand_var < 1e-20 else "ok"

    fe_params = mdf.fe_params
    fe_se = mdf.bse_fe
    fe_pvalues = mdf.pvalues.loc[fe_params.index]
    fe_tvalues = mdf.tvalues.loc[fe_params.index]
    ci = mdf.conf_int().loc[fe_params.index]

    return {
        "status": status,
        "n_obs": n_obs,
        "n_subjects": n_subjects,
        "intercept": float(fe_params["Intercept"]),
        "intercept_se": float(fe_se["Intercept"]),
        "slope": float(fe_params[TIME_MAIN_COEF_LABEL]),
        "slope_se": float(fe_se[TIME_MAIN_COEF_LABEL]),
        "slope_z": float(fe_tvalues[TIME_MAIN_COEF_LABEL]),
        "slope_p": float(fe_pvalues[TIME_MAIN_COEF_LABEL]),
        "slope_ci_low": float(ci.loc[TIME_MAIN_COEF_LABEL, 0]),
        "slope_ci_high": float(ci.loc[TIME_MAIN_COEF_LABEL, 1]),
        "residual_var": float(mdf.scale),
        "random_intercept_var": rand_var,
        "converged": bool(mdf.converged),
    }


# --- Effect-size derivations ---

def compute_predicted_diff_at_horizon(
    b: float, se: float, control_bl_mean: float, t: float
) -> Tuple[float, float, float]:
    """Convert the LME interaction coefficient to a predicted t-year
    PD-vs-Control difference, as % of the Control baseline mean.

        predicted_diff_pct = 100 * (t * b) / control_bl_mean

    The same LME fit supports any horizon; only the multiplier
    changes. The Wald z and p-value on b are horizon-invariant.

    Returns (point, ci_low, ci_high). CI uses 1.96 * SE on b.
    """
    if control_bl_mean == 0 or np.isnan(control_bl_mean):
        return (np.nan, np.nan, np.nan)
    point = 100.0 * (t * b) / control_bl_mean
    ci_low = 100.0 * (t * (b - 1.96 * se)) / control_bl_mean
    ci_high = 100.0 * (t * (b + 1.96 * se)) / control_bl_mean
    return (point, ci_low, ci_high)


def compute_predicted_diff_5y(
    b: float, se: float, control_bl_mean: float
) -> Tuple[float, float, float]:
    """Predicted PD-vs-Control difference at t=5."""
    return compute_predicted_diff_at_horizon(b, se, control_bl_mean, 5.0)


def compute_d_adj_5y(
    b: float, se: float, residual_var: float, df_resid: int
) -> Tuple[float, float, float]:
    """Standardized 5-year effect size with Hedges' small-sample J.

        d_uncorrected = (5 * b) / sqrt(residual_var)
        J = 1 - 3 / (4 * df_resid - 1)
        d_adj_5y = d_uncorrected * J

    Returns (d_adj, ci_low, ci_high). CI from 1.96 * SE_d.
    """
    if residual_var <= 0 or np.isnan(residual_var) or df_resid <= 1:
        return (np.nan, np.nan, np.nan)
    sigma = np.sqrt(residual_var)
    J = 1.0 - 3.0 / (4.0 * df_resid - 1.0)
    d = (5.0 * b) / sigma * J
    se_d = (5.0 * se) / sigma * J
    return (d, d - 1.96 * se_d, d + 1.96 * se_d)


def compute_within_group_pct_change(
    slope: float, slope_se: float, mean_bl: float, t: float
) -> Tuple[float, float, float]:
    """Convert a within-group LME slope to a predicted t-year percent
    change for that group, relative to the group's own BL mean.

        pct_t = 100 * (t * slope) / mean_bl

    Returns (point, ci_low, ci_high). NaN if mean_bl is NaN/zero or
    inputs are NaN.
    """
    if (
        pd.isna(slope) or pd.isna(slope_se)
        or pd.isna(mean_bl) or mean_bl == 0
    ):
        return (np.nan, np.nan, np.nan)
    point = 100.0 * (t * slope) / mean_bl
    ci_low = 100.0 * (t * (slope - 1.96 * slope_se)) / mean_bl
    ci_high = 100.0 * (t * (slope + 1.96 * slope_se)) / mean_bl
    return (point, ci_low, ci_high)


# --- Test B and Test A ---

def run_test_b(cohort_long: pd.DataFrame) -> pd.DataFrame:
    """Between-group LME for every ROI in cohort_long.

    Returns one row per ROI with all columns from fit_between_group_lme,
    plus derived effect sizes and per-tissue FDR.
    """
    rois = sorted(cohort_long["ROI"].unique())
    rows: List[Dict[str, Any]] = []

    for roi in rois:
        roi_data = cohort_long[cohort_long["ROI"] == roi]
        tissue = roi_data["Tissue"].iloc[0]

        result = fit_between_group_lme(roi_data)
        result["ROI"] = roi
        result["Tissue"] = tissue

        # Derived effect sizes (skipped if status != ok).
        if result["status"] in ("ok", "singular_fit"):
            # 5-year horizon (headline).
            pdif5 = compute_predicted_diff_at_horizon(
                result["interaction_coef"], result["interaction_se"],
                result["mean_bl_control"], 5.0
            )
            result["predicted_diff_5y_pct"] = pdif5[0]
            result["predicted_diff_5y_ci_low"] = pdif5[1]
            result["predicted_diff_5y_ci_high"] = pdif5[2]

            # 3-year horizon (CSV-only; mathematically (3/5) of the
            # 5y values from the same LME fit).
            pdif3 = compute_predicted_diff_at_horizon(
                result["interaction_coef"], result["interaction_se"],
                result["mean_bl_control"], 3.0
            )
            result["predicted_diff_3y_pct"] = pdif3[0]
            result["predicted_diff_3y_ci_low"] = pdif3[1]
            result["predicted_diff_3y_ci_high"] = pdif3[2]

            d = compute_d_adj_5y(
                result["interaction_coef"], result["interaction_se"],
                result["residual_var"], result["df_resid"]
            )
            result["d_adj_5y"] = d[0]
            result["d_adj_5y_ci_low"] = d[1]
            result["d_adj_5y_ci_high"] = d[2]
        else:
            result["predicted_diff_5y_pct"] = np.nan
            result["predicted_diff_5y_ci_low"] = np.nan
            result["predicted_diff_5y_ci_high"] = np.nan
            result["predicted_diff_3y_pct"] = np.nan
            result["predicted_diff_3y_ci_low"] = np.nan
            result["predicted_diff_3y_ci_high"] = np.nan
            result["d_adj_5y"] = np.nan
            result["d_adj_5y_ci_low"] = np.nan
            result["d_adj_5y_ci_high"] = np.nan

        rows.append(result)

    df = pd.DataFrame(rows)

    # Apply per-tissue FDR.
    df = apply_tissue_fdr(df, p_col="interaction_p")

    return df


def run_test_a(cohort_long: pd.DataFrame) -> pd.DataFrame:
    """Within-group LME for every (ROI, Group) combination.

    Returns 24 rows (12 ROIs x 2 groups) with each group's slope,
    significance, BL mean, and predicted %-change at 3y/5y horizons
    (with 95% CIs). Descriptive: tests each group's slope, not
    PD vs Control.
    """
    rois = sorted(cohort_long["ROI"].unique())
    groups_in_order = ["Pasient", "Control"]
    rows: List[Dict[str, Any]] = []

    for roi in rois:
        roi_data = cohort_long[cohort_long["ROI"] == roi]
        tissue = roi_data["Tissue"].iloc[0]

        for group in groups_in_order:
            result = fit_within_group_lme(roi_data, group)
            result["ROI"] = roi
            result["Tissue"] = tissue
            result["Group_clinical"] = group

            # Per-group BL mean (used for within-group %-change).
            bl_group = roi_data[
                (roi_data["Group_clinical"] == group)
                & (roi_data["Timepoint"] == "BL")
            ]
            if len(bl_group) > 0:
                mean_bl_group = float(bl_group["Value"].mean())
            else:
                mean_bl_group = np.nan
            result["mean_bl_group"] = mean_bl_group

            # Within-group predicted % change at 5y and 3y, relative
            # to this group's own BL mean.
            if result["status"] in ("ok", "singular_fit"):
                p5 = compute_within_group_pct_change(
                    result["slope"], result["slope_se"],
                    mean_bl_group, 5.0,
                )
                result["pct_5y"] = p5[0]
                result["pct_5y_ci_low"] = p5[1]
                result["pct_5y_ci_high"] = p5[2]

                p3 = compute_within_group_pct_change(
                    result["slope"], result["slope_se"],
                    mean_bl_group, 3.0,
                )
                result["pct_3y"] = p3[0]
                result["pct_3y_ci_low"] = p3[1]
                result["pct_3y_ci_high"] = p3[2]
            else:
                result["pct_5y"] = np.nan
                result["pct_5y_ci_low"] = np.nan
                result["pct_5y_ci_high"] = np.nan
                result["pct_3y"] = np.nan
                result["pct_3y_ci_low"] = np.nan
                result["pct_3y_ci_high"] = np.nan

            rows.append(result)

    df = pd.DataFrame(rows)

    # Per-tissue FDR (per Group x Tissue).
    df["p_fdr_tissue"] = np.nan
    df["sig_star_tissue"] = ""
    for grp in groups_in_order:
        for tissue in ["subcortical", "cortical"]:
            mask = (df["Group_clinical"] == grp) & (df["Tissue"] == tissue)
            valid = mask & df["slope_p"].notna()
            if valid.sum() == 0:
                continue
            p_vals = df.loc[valid, "slope_p"].to_numpy()
            _, p_fdr, _, _ = multipletests(p_vals, method="fdr_bh")
            df.loc[valid, "p_fdr_tissue"] = p_fdr

    df["sig_star_tissue"] = df["p_fdr_tissue"].apply(_pvalue_to_star)

    return df


def apply_tissue_fdr(
    df: pd.DataFrame, p_col: str
) -> pd.DataFrame:
    """BH-FDR correction per tissue family.

    Subcortical and cortical (6 tests each) treated as separate
    inferential families. Adds p_fdr_tissue and sig_star_tissue.
    """
    df = df.copy()

    df["p_fdr_tissue"] = np.nan
    for tissue in ["subcortical", "cortical"]:
        mask = df["Tissue"] == tissue
        valid = mask & df[p_col].notna()
        if valid.sum() == 0:
            continue
        p_vals = df.loc[valid, p_col].to_numpy()
        _, p_fdr, _, _ = multipletests(p_vals, method="fdr_bh")
        df.loc[valid, "p_fdr_tissue"] = p_fdr
    df["sig_star_tissue"] = df["p_fdr_tissue"].apply(_pvalue_to_star)

    return df


def _pvalue_to_star(p: float) -> str:
    """Convert a p-value to its significance star annotation."""
    if pd.isna(p):
        return "n.s."
    for thresh, star in STAR_THRESHOLDS:
        if p < thresh:
            return star
    return "n.s."


# --- ROI ordering ---

def compute_roi_ordering(test_b_results: pd.DataFrame) -> pd.DataFrame:
    """Canonical ROI display order: descending |d_adj_5y| within each
    tissue, subcortical block first.

    Returns columns ROI, Tissue, ordering_rank.
    """
    rows = []
    for tissue in ["subcortical", "cortical"]:
        block = test_b_results[
            test_b_results["Tissue"] == tissue
        ].copy()
        block["abs_d"] = block["d_adj_5y"].abs()
        block_sorted = block.sort_values(
            "abs_d", ascending=False
        ).reset_index(drop=True)
        for rank, row in block_sorted.iterrows():
            rows.append({
                "ROI": row["ROI"],
                "Tissue": tissue,
                "ordering_rank": rank,
            })
    return pd.DataFrame(rows)


# --- Risk-A harmonization check ---

def _safe_abs_ratio(a: float, b: float) -> float:
    """abs(a)/abs(b) with NaN/zero guards. NaN otherwise."""
    if pd.isna(a) or pd.isna(b) or b == 0:
        return np.nan
    return abs(a) / abs(b)


def run_risk_a_all_rois(
    test_b_harmonized: pd.DataFrame,
    cohort_long_unharm: pd.DataFrame,
) -> pd.DataFrame:
    """Per-ROI harmonized-vs-unharmonized comparison of Test B's
    Group x Years_from_BL interaction.

    For every ROI in test_b_harmonized, refit the same Test B LME on
    the unharmonized cohort_long subset, then compute sign-match and
    magnitude fields. If the unharmonized fit fails, those fields are
    NaN.

    Descriptive: shows whether direction and approximate magnitude
    of the primary inferential term survive harmonization. Not an
    inferential gate.
    """
    rows: List[Dict[str, Any]] = []

    test_b_idx = test_b_harmonized.set_index("ROI")
    rois = list(test_b_harmonized["ROI"])

    for roi in rois:
        harm_row = test_b_idx.loc[roi]
        tissue = harm_row["Tissue"]

        # Harmonized values come straight from Test B output.
        beta_h = float(harm_row["interaction_coef"])
        p_h = float(harm_row["interaction_p"])
        d_h = float(harm_row["d_adj_5y"])
        p_fdr_t_h = float(harm_row["p_fdr_tissue"])
        sig_t_h = harm_row["sig_star_tissue"]

        # Refit Test B on the unharmonized data for this ROI.
        unharm_roi = cohort_long_unharm[cohort_long_unharm["ROI"] == roi]
        if unharm_roi.empty:
            beta_u = np.nan
            p_u = np.nan
            d_u = np.nan
        else:
            unharm_fit = fit_between_group_lme(unharm_roi)
            if unharm_fit["status"] in ("ok", "singular_fit"):
                beta_u = float(unharm_fit["interaction_coef"])
                p_u = float(unharm_fit["interaction_p"])
                d_tuple = compute_d_adj_5y(
                    unharm_fit["interaction_coef"],
                    unharm_fit["interaction_se"],
                    unharm_fit["residual_var"],
                    unharm_fit["df_resid"],
                )
                d_u = float(d_tuple[0]) if not np.isnan(d_tuple[0]) else np.nan
            else:
                beta_u = np.nan
                p_u = np.nan
                d_u = np.nan

        # Sign-match guarded for NaN.
        if pd.isna(beta_h) or pd.isna(beta_u):
            sign_match: Any = np.nan
        else:
            sign_match = bool(np.sign(beta_h) == np.sign(beta_u))

        rows.append({
            "ROI": roi,
            "Tissue": tissue,
            "beta_harmonized": beta_h,
            "beta_unharmonized": beta_u,
            "beta_abs_ratio": _safe_abs_ratio(beta_h, beta_u),
            "beta_difference": (
                beta_h - beta_u
                if not (pd.isna(beta_h) or pd.isna(beta_u))
                else np.nan
            ),
            "sign_match": sign_match,
            "p_harmonized": p_h,
            "p_unharmonized": p_u,
            "d_harmonized": d_h,
            "d_unharmonized": d_u,
            "d_abs_ratio": _safe_abs_ratio(d_h, d_u),
            "p_fdr_tissue_harmonized": p_fdr_t_h,
            "sig_star_tissue_harmonized": sig_t_h,
        })

    return pd.DataFrame(rows)


def write_risk_a_csv(df: pd.DataFrame, path: Path) -> None:
    """Write risk_a_all_rois.csv with stable column order."""
    col_order = [
        "ROI", "Tissue",
        "beta_harmonized", "beta_unharmonized",
        "beta_abs_ratio", "beta_difference", "sign_match",
        "p_harmonized", "p_unharmonized",
        "d_harmonized", "d_unharmonized", "d_abs_ratio",
        "p_fdr_tissue_harmonized", "sig_star_tissue_harmonized",
    ]
    out = df[col_order].copy()
    out.to_csv(path, index=False, float_format=CSV_FLOAT_FORMAT)
    print(f"  Wrote: {path} ({len(out)} rows, {len(out.columns)} cols)")


def write_risk_a_summary_text(
    risk_a: pd.DataFrame,
    summary_path: Path,
    test_b_results: pd.DataFrame,
    cohort_long_unharm_present: bool,
) -> None:
    """Write lme_sanity_check.txt: an all-ROI Risk-A summary.

    Aggregate stats (sign-match counts, Pearson r of harmonized vs
    unharmonized beta_interaction, FDR-sig sign-preservation), then a
    per-ROI table and the FDR-significant findings on harmonized data.
    """
    lines: List[str] = ["LME SANITY CHECK -- Risk A: harmonization preservation"]
    lines.append("")

    if not cohort_long_unharm_present:
        lines.append("STATUS: skipped (cohort_long_unharmonized.csv not found)")
        summary_path.write_text("\n".join(lines) + "\n")
        print(f"  Wrote: {summary_path} (sanity check skipped)")
        return

    # Aggregate stats across ROIs.
    n_total = len(risk_a)
    sign_match_col = risk_a["sign_match"]
    n_same_sign = int(sign_match_col.eq(True).sum())
    n_diff_sign = int(sign_match_col.eq(False).sum())
    n_unknown_sign = int(sign_match_col.isna().sum())

    # Pearson correlation between harmonized and unharmonized beta.
    bh = risk_a["beta_harmonized"].to_numpy(dtype=float)
    bu = risk_a["beta_unharmonized"].to_numpy(dtype=float)
    mask = np.isfinite(bh) & np.isfinite(bu)
    if mask.sum() >= 3:
        r_pearson = float(np.corrcoef(bh[mask], bu[mask])[0, 1])
    else:
        r_pearson = np.nan

    # Median magnitude ratio.
    ratio_col = risk_a["beta_abs_ratio"]
    med_ratio = float(ratio_col.median(skipna=True))
    p25_ratio = float(ratio_col.quantile(0.25, interpolation="linear"))
    p75_ratio = float(ratio_col.quantile(0.75, interpolation="linear"))

    # FDR-significant harmonized subset.
    sig_mask = risk_a["sig_star_tissue_harmonized"].isin(["*", "**", "***"])
    n_sig_harm = int(sig_mask.sum())
    sig_same_sign = int(
        (sign_match_col.eq(True) & sig_mask).sum()
    )

    lines.append("ALL-ROI SUMMARY")
    lines.append("")
    lines.append(
        f"  Total ROIs compared:              {n_total}"
    )
    lines.append(
        f"  ROIs with same beta_interaction sign: {n_same_sign} / {n_total}"
    )
    if n_diff_sign:
        lines.append(
            f"  ROIs with opposite sign:           {n_diff_sign} / {n_total}"
        )
    if n_unknown_sign:
        lines.append(
            f"  ROIs with undetermined sign:       {n_unknown_sign} / {n_total}"
        )
    if not np.isnan(r_pearson):
        lines.append(
            f"  Pearson r (b_harm vs b_unharm):    {r_pearson:.3f}"
        )
    if not np.isnan(med_ratio):
        lines.append(
            f"  |b_harm| / |b_unharm| median:      "
            f"{med_ratio:.3f}  (IQR {p25_ratio:.3f}-{p75_ratio:.3f})"
        )
    lines.append(
        f"  Harmonized FDR-tissue sig. ROIs:   {n_sig_harm} / {n_total}"
    )
    lines.append(
        f"  ...of those, same sign in unharm.:  {sig_same_sign} / {n_sig_harm}"
        if n_sig_harm > 0
        else "  ...of those, same sign in unharm.:  N/A (no FDR-sig ROIs)"
    )
    lines.append("")

    # Per-ROI table (sorted by tissue, then by harmonized |d|).
    lines.append("PER-ROI COMPARISON (sorted by Tissue, then |d_harm| desc.)")
    lines.append("")
    lines.append(
        f"{'ROI':<26}{'Tissue':<12}"
        f"{'b_harm':>12}{'b_unharm':>12}"
        f"{'|H/U|':>8}"
        f"{'sign?':>7}"
        f"{'sig':>6}"
    )
    lines.append("-" * 72)

    table_rows = risk_a.copy()
    table_rows["abs_d_harm"] = table_rows["d_harmonized"].abs()
    tissue_order = ["subcortical", "cortical"]
    for tissue in tissue_order:
        block = table_rows[table_rows["Tissue"] == tissue].sort_values(
            "abs_d_harm", ascending=False
        )
        for _, r in block.iterrows():
            ratio_str = (
                f"{r['beta_abs_ratio']:.2f}"
                if not pd.isna(r["beta_abs_ratio"]) else " - "
            )
            sign_str = (
                "Y" if r["sign_match"] is True
                else ("N" if r["sign_match"] is False else "?")
            )
            sig = r["sig_star_tissue_harmonized"]
            sig_str = sig if sig in ("*", "**", "***") else ""
            lines.append(
                f"{r['ROI']:<26}{r['Tissue']:<12}"
                f"{r['beta_harmonized']:>12.3e}"
                f"{r['beta_unharmonized']:>12.3e}"
                f"{ratio_str:>8}"
                f"{sign_str:>7}"
                f"{sig_str:>6}"
            )

    lines.append("")

    # FDR-significant findings on harmonized data.
    lines.append("FULL-PANEL FDR-SIGNIFICANT FINDINGS (harmonized data)")
    lines.append("")
    sig = test_b_results[
        test_b_results["sig_star_tissue"].isin(["*", "**", "***"])
    ].copy()
    if sig.empty:
        lines.append("(none -- no ROI surviving FDR within tissue)")
    else:
        lines.append(
            f"{'ROI':<25} {'Tissue':<13} "
            f"{'d_adj_5y':>10} {'p_fdr_tissue':>14} {'sig':>5}"
        )
        lines.append("-" * 72)
        for _, r in sig.sort_values("d_adj_5y").iterrows():
            lines.append(
                f"{r['ROI']:<25} {r['Tissue']:<13} "
                f"{r['d_adj_5y']:>10.3f} {r['p_fdr_tissue']:>14.4f}"
                f" {r['sig_star_tissue']:>5}"
            )

    summary_path.write_text("\n".join(lines) + "\n")
    print(f"  Wrote: {summary_path}")


# --- Figure helpers ---

def _apply_figure_style() -> None:
    """Set matplotlib rcParams for thesis-quality figures."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": FONTSIZE_TICK,
        "axes.titlesize": FONTSIZE_TITLE,
        "axes.labelsize": FONTSIZE_LABEL,
        "xtick.labelsize": FONTSIZE_TICK,
        "ytick.labelsize": FONTSIZE_TICK,
        "legend.fontsize": FONTSIZE_ANNOTATION,
        "figure.titlesize": FONTSIZE_TITLE,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "savefig.bbox": "tight",
    })


def _pvalue_to_label(p: float) -> str:
    """Star label or 'n.s.' for figure annotations."""
    if pd.isna(p):
        return "n.s."
    for thresh, star in STAR_THRESHOLDS:
        if p < thresh:
            return star
    return "n.s."


def _resolve_auto_xlim(values: np.ndarray, floor: float) -> float:
    """Symmetric xlim half-width: max(|values|) rounded up to the next
    0.5, but never less than `floor`."""
    if len(values) == 0 or np.all(np.isnan(values)):
        return floor
    m = float(np.nanmax(np.abs(values)))
    rounded = math.ceil(m / 0.5) * 0.5
    return max(rounded, floor)


def _format_roi_label(roi: str, tissue: str) -> str:
    """ROI label string for figures.

    matplotlib handles italics via fontstyle on the Axes call rather
    than the string; the caller picks the style.
    """
    return roi


def _build_lme_predictions(
    fit_result: Dict[str, Any],
    roi_data: pd.DataFrame,
    n_grid: int = 50,
) -> Dict[str, Any]:
    """Per-group fitted-mean trajectories with 95% CI ribbons.

    For each group g in (Control, Pasient):
        Y_hat(t) = b0 + b1*[g=PD] + b2*t + b3*[g=PD]*t
                 + b4*age_mean + b5*sex_mean
    with t swept across the cohort's Years_from_BL range. Variance
    of Y_hat(t) via delta method using the LME's fixed-effects
    covariance matrix.

    Refits the LME internally to access the full Cov(b); fit_result
    only carries diagonal SEs.

    Returns a dict with t_grid, pd_mean, pd_ci_low, pd_ci_high,
    co_mean, co_ci_low, co_ci_high, age_mean, sex_mean_M.
    """
    formula = (
        "Value ~ C(Group_clinical, Treatment(reference='Control'))"
        " * Years_from_BL"
        " + age_at_BL"
        " + C(PatientSex_clinical)"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        md = smf.mixedlm(formula, data=roi_data,
                         groups=roi_data["Subject"])
        mdf = md.fit(reml=True, method="lbfgs")

    # Cohort-mean covariates (one row per subject so subjects with
    # multiple observations aren't double-weighted).
    per_subj = roi_data.drop_duplicates(subset="Subject")
    age_mean = float(per_subj["age_at_BL"].mean())
    sex_mean_M = float(
        (per_subj["PatientSex_clinical"] == "M").mean()
    )

    # Prediction grid.
    t_min = float(roi_data["Years_from_BL"].min())
    t_max = float(roi_data["Years_from_BL"].max())
    t_grid = np.linspace(t_min, t_max, n_grid)

    # Build design rows. fe_params index order is:
    # Intercept, Group[T.Pasient], Sex[T.M], Years_from_BL,
    # Group[T.Pasient]:Years_from_BL, age_at_BL.
    param_order = list(mdf.fe_params.index)
    cov_fe = mdf.cov_params().loc[param_order, param_order].to_numpy()

    def design_row(group: str, t: float) -> np.ndarray:
        is_pd = 1.0 if group == "Pasient" else 0.0
        row = np.zeros(len(param_order))
        for i, p in enumerate(param_order):
            if p == "Intercept":
                row[i] = 1.0
            elif p == GROUP_MAIN_COEF_LABEL:
                row[i] = is_pd
            elif p == "C(PatientSex_clinical)[T.M]":
                row[i] = sex_mean_M
            elif p == TIME_MAIN_COEF_LABEL:
                row[i] = t
            elif p == INTERACTION_COEF_LABEL:
                row[i] = is_pd * t
            elif p == "age_at_BL":
                row[i] = age_mean
            else:
                # Unknown parameter; leave zero.
                pass
        return row

    beta = mdf.fe_params.loc[param_order].to_numpy()

    pd_pred = np.empty(n_grid)
    pd_se = np.empty(n_grid)
    co_pred = np.empty(n_grid)
    co_se = np.empty(n_grid)

    for i, t in enumerate(t_grid):
        x_pd = design_row("Pasient", t)
        x_co = design_row("Control", t)
        pd_pred[i] = float(x_pd @ beta)
        co_pred[i] = float(x_co @ beta)
        pd_se[i] = float(np.sqrt(x_pd @ cov_fe @ x_pd))
        co_se[i] = float(np.sqrt(x_co @ cov_fe @ x_co))

    return {
        "t_grid": t_grid,
        "pd_mean": pd_pred,
        "pd_ci_low": pd_pred - 1.96 * pd_se,
        "pd_ci_high": pd_pred + 1.96 * pd_se,
        "co_mean": co_pred,
        "co_ci_low": co_pred - 1.96 * co_se,
        "co_ci_high": co_pred + 1.96 * co_se,
        "age_mean": age_mean,
        "sex_mean_M": sex_mean_M,
    }


# --- Output writers ---

def write_test_b_csv(df: pd.DataFrame, path: Path) -> None:
    """Write stats_lme_between_group.csv with stable column order."""
    col_order = [
        "ROI", "Tissue",
        "n_obs", "n_subjects", "n_subjects_pd", "n_subjects_co",
        "status", "converged",
        "intercept", "intercept_se",
        "group_main_coef", "group_main_se", "group_main_p",
        "time_main_coef", "time_main_se", "time_main_p",
        "interaction_coef", "interaction_se",
        "interaction_z", "interaction_p",
        "interaction_ci_low", "interaction_ci_high",
        "p_fdr_tissue", "sig_star_tissue",
        "predicted_diff_5y_pct",
        "predicted_diff_5y_ci_low", "predicted_diff_5y_ci_high",
        "predicted_diff_3y_pct",
        "predicted_diff_3y_ci_low", "predicted_diff_3y_ci_high",
        "d_adj_5y", "d_adj_5y_ci_low", "d_adj_5y_ci_high",
        "age_coef", "age_se",
        "sex_coef", "sex_se",
        "residual_var", "random_intercept_var",
        "df_resid",
        "mean_bl_control", "mean_bl_pd",
    ]
    out = df[col_order].sort_values(["Tissue", "ROI"]).reset_index(drop=True)
    out.to_csv(path, index=False, float_format=CSV_FLOAT_FORMAT)
    print(f"  Wrote: {path} ({len(out)} rows, {len(out.columns)} cols)")


def write_test_a_csv(df: pd.DataFrame, path: Path) -> None:
    """Write stats_lme_within_group.csv."""
    col_order = [
        "Group_clinical", "ROI", "Tissue",
        "n_obs", "n_subjects",
        "status", "converged",
        "intercept", "intercept_se",
        "slope", "slope_se", "slope_z", "slope_p",
        "slope_ci_low", "slope_ci_high",
        "mean_bl_group",
        "pct_5y", "pct_5y_ci_low", "pct_5y_ci_high",
        "pct_3y", "pct_3y_ci_low", "pct_3y_ci_high",
        "p_fdr_tissue", "sig_star_tissue",
        "residual_var", "random_intercept_var",
    ]
    out = df[col_order].sort_values(
        ["Group_clinical", "Tissue", "ROI"]
    ).reset_index(drop=True)
    out.to_csv(path, index=False, float_format=CSV_FLOAT_FORMAT)
    print(f"  Wrote: {path} ({len(out)} rows, {len(out.columns)} cols)")


def write_roi_ordering_csv(df: pd.DataFrame, path: Path) -> None:
    """Write roi_display_order.csv."""
    df.to_csv(path, index=False)
    print(f"  Wrote: {path} ({len(df)} rows)")


# --- Figures ---

# Figure 1: forest plot of d_adj_5y across all 12 ROIs.

def plot_figure1_forest(
    test_b: pd.DataFrame,
    ordering: pd.DataFrame,
    out_path: Path,
    xlim: Optional[float] = None,
) -> None:
    """Headline inferential figure: d_adj_5y per ROI with 95% CI.

    Subcortical block on top, cortical block below; within each
    block, ROIs ordered by descending |d_adj_5y|. Filled markers
    indicate FDR-significance (p_fdr_tissue); hollow markers are
    n.s. Dashed reference lines at d=+/-0.2/0.5/0.8 are conventional
    benchmarks for interpretation only.
    """
    _apply_figure_style()

    # Row order: subcortical first (sorted), then cortical.
    sub_block = ordering[ordering["Tissue"] == "subcortical"].sort_values(
        "ordering_rank"
    )
    cor_block = ordering[ordering["Tissue"] == "cortical"].sort_values(
        "ordering_rank"
    )
    rows = list(sub_block["ROI"]) + list(cor_block["ROI"])
    n_rows = len(rows)
    n_subcortical = len(sub_block)

    # Effect-size data per ROI.
    test_b_idx = test_b.set_index("ROI")
    d_vals = np.array([test_b_idx.loc[r, "d_adj_5y"] for r in rows])
    d_low = np.array([test_b_idx.loc[r, "d_adj_5y_ci_low"] for r in rows])
    d_high = np.array([test_b_idx.loc[r, "d_adj_5y_ci_high"] for r in rows])
    p_fdr = np.array([test_b_idx.loc[r, "p_fdr_tissue"] for r in rows])
    tissues = [test_b_idx.loc[r, "Tissue"] for r in rows]

    # x-limit includes CI bars.
    if xlim is None:
        all_d_vals = np.concatenate([d_vals, d_low, d_high])
        xlim = _resolve_auto_xlim(all_d_vals, floor=1.5)

    fig, ax = plt.subplots(figsize=FIG1_FOREST_DIMS)

    # Highest |d| at top; first entry in `rows` lands at top.
    y_positions = np.arange(n_rows)[::-1]

    # Vertical reference lines at d=0 and at +/-0.2/0.5/0.8.
    ax.axvline(0.0, color=COLOR_GRAY_DARK, linewidth=0.8,
               alpha=0.6, zorder=1)
    for ref in REFERENCE_D_VALUES:
        for sign in (-1, 1):
            ax.axvline(sign * ref, color=COLOR_GRAY_LIGHT,
                       linewidth=0.6, linestyle="--",
                       alpha=0.8, zorder=1)

    # Per-ROI marker + CI.
    for i, (y, d, lo, hi, p) in enumerate(zip(
        y_positions, d_vals, d_low, d_high, p_fdr
    )):
        is_sig = (p < 0.05) and not pd.isna(p)
        # CI bars
        ax.plot(
            [lo, hi], [y, y],
            color=COLOR_NEUTRAL if is_sig else COLOR_GRAY_DARK,
            linewidth=1.4 if is_sig else 1.0,
            alpha=1.0 if is_sig else 0.7, zorder=3,
        )
        # End caps
        cap_h = 0.18
        for x in (lo, hi):
            ax.plot(
                [x, x], [y - cap_h, y + cap_h],
                color=COLOR_NEUTRAL if is_sig else COLOR_GRAY_DARK,
                linewidth=1.4 if is_sig else 1.0,
                alpha=1.0 if is_sig else 0.7, zorder=3,
            )
        # Marker
        marker_face = COLOR_NEUTRAL if is_sig else "white"
        ax.plot(
            d, y, marker="o", markersize=7,
            markerfacecolor=marker_face,
            markeredgecolor=COLOR_NEUTRAL,
            markeredgewidth=1.4, zorder=4,
        )

    # Significance annotation column (right of the plot area).
    annot_x = xlim * 1.04
    for y, p in zip(y_positions, p_fdr):
        label = _pvalue_to_label(p)
        is_sig = label != "n.s."
        ax.text(
            annot_x, y, label,
            ha="left", va="center",
            fontsize=FONTSIZE_ANNOTATION,
            color=COLOR_NEUTRAL if is_sig else COLOR_GRAY_DARK,
            fontweight=600 if is_sig else "normal",
        )

    # Tissue divider between subcortical and cortical blocks.
    if n_subcortical > 0 and n_subcortical < n_rows:
        sep_y = (n_rows - n_subcortical) - 0.5
        ax.axhline(
            y=sep_y, xmin=0.0, xmax=1.0,
            color=COLOR_GRAY_DARK, linewidth=0.8,
            alpha=0.6, zorder=2,
        )

    # Y-axis: ROI labels, italic for subcortical.
    ax.set_yticks(y_positions)
    ax.set_yticklabels(rows)
    for tick_label, tissue in zip(ax.get_yticklabels(), tissues):
        if tissue == "subcortical":
            tick_label.set_style("italic")
        else:
            tick_label.set_weight(500)

    # X-axis.
    ax.set_xlim(-xlim, xlim)
    ax.set_xlabel(
        "Standardized 5-year slope-difference effect ($d_{adj,5Y}$)"
    )

    # Y-axis bounds (avoid clipping markers).
    ax.set_ylim(-0.6, n_rows - 0.4)

    # Right-margin tissue labels (rotated 90deg).
    if n_subcortical > 0:
        sub_y_center = np.mean(y_positions[:n_subcortical])
        sub_y_frac = (sub_y_center - (-0.6)) / ((n_rows - 0.4) - (-0.6))
        ax.text(
            1.10, sub_y_frac, "Subcortical volume",
            rotation=270, ha="center", va="center",
            fontsize=FONTSIZE_ANNOTATION, color=COLOR_GRAY_DARK,
            fontweight=500, transform=ax.transAxes,
        )
    if n_subcortical < n_rows:
        cor_y_center = np.mean(y_positions[n_subcortical:])
        cor_y_frac = (cor_y_center - (-0.6)) / ((n_rows - 0.4) - (-0.6))
        ax.text(
            1.10, cor_y_frac, "Cortical thickness",
            rotation=270, ha="center", va="center",
            fontsize=FONTSIZE_ANNOTATION, color=COLOR_GRAY_DARK,
            fontweight=500, transform=ax.transAxes,
        )

    # Title and subtitle with explicit top padding.
    fig.suptitle(
        "Standardized atrophy difference (PD - Control) over 5 years",
        fontsize=FONTSIZE_TITLE, y=0.98,
    )
    fig.text(
        0.5, 0.93,
        "LME-derived $d_{adj,5Y}$ with 95% CI; FDR within tissue",
        ha="center", va="center",
        fontsize=FONTSIZE_SUBTITLE,
        color=COLOR_GRAY_DARK, style="italic",
    )
    fig.subplots_adjust(top=0.88, bottom=0.12, left=0.18, right=0.88)


    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=FIG_DPI)
    plt.close(fig)
    print(f"  Wrote: {out_path}")


# Per-panel trajectory renderer (shared by Figures 2 and 3).

def _plot_trajectory_panel(
    ax: "plt.Axes",
    roi_data: pd.DataFrame,
    roi: str,
    test_b_row: pd.Series,
    show_d_annotation: bool = True,
    line_alpha: float = 0.20,
    line_width: float = 0.6,
    show_xlabel: bool = True,
    show_ylabel: bool = True,
    title_fontsize: int = FONTSIZE_LABEL,
    show_individual_lines: bool = True,
) -> None:
    """Render one trajectory panel: per-subject thin lines (optional)
    plus bold LME-fitted group means and 95% CI ribbons.

    Used by Figure 2 (FDR-sig) and Figures 3/4 (subcortical / cortical
    splits of the full panel).

    Y-axis: per-panel auto-scaled raw units (Volume x 10^3 / eTIV for
    subcortical; thickness in mm for cortical). The d_adj annotation
    handles cross-panel slope comparison.

    X-axis: continuous Years_from_BL capped at [0, 5.7]; ticks at 0,
    3, 5. Titles show "ROI ***" or "ROI (n.s.)".

    show_individual_lines=False omits the per-subject spaghetti
    (clean variant; used for *_clean.png).
    """
    tissue = test_b_row["Tissue"]
    d_adj = test_b_row["d_adj_5y"]
    sig = test_b_row["sig_star_tissue"]

    # Subcortical: x 1000 for readability.
    y_scale = 1000.0 if tissue == "subcortical" else 1.0

    # Vertical reference lines at nominal timepoints (zorder=0).
    for t in (0.0, 3.0, 5.0):
        ax.axvline(t, color=COLOR_GRAY_LIGHT, linewidth=0.5,
                   alpha=0.5, zorder=0)

    # Per-subject lines (optional).
    if show_individual_lines:
        for subj, group_data in roi_data.groupby("Subject"):
            gd = group_data.sort_values("Years_from_BL")
            group = gd["Group_clinical"].iloc[0]
            color = COLOR_PD if group == "Pasient" else COLOR_CONTROL
            ax.plot(
                gd["Years_from_BL"], gd["Value"] * y_scale,
                color=color, alpha=line_alpha, linewidth=line_width,
                zorder=1,
            )

    # LME predictions; refits internally for full Cov(b).
    fit_result_dummy = {
        "residual_var": float(roi_data["Value"].var())
    }
    preds = _build_lme_predictions(
        fit_result_dummy, roi_data, n_grid=50
    )

    # CI ribbons
    ax.fill_between(
        preds["t_grid"],
        preds["pd_ci_low"] * y_scale,
        preds["pd_ci_high"] * y_scale,
        color=COLOR_PD, alpha=0.20, linewidth=0, zorder=2,
    )
    ax.fill_between(
        preds["t_grid"],
        preds["co_ci_low"] * y_scale,
        preds["co_ci_high"] * y_scale,
        color=COLOR_CONTROL, alpha=0.20, linewidth=0, zorder=2,
    )

    # Bold group-mean lines
    ax.plot(
        preds["t_grid"], preds["pd_mean"] * y_scale,
        color=COLOR_PD, linewidth=2.2, zorder=3,
    )
    ax.plot(
        preds["t_grid"], preds["co_mean"] * y_scale,
        color=COLOR_CONTROL, linewidth=2.2, zorder=3,
    )

    # Title with FDR-significance suffix.
    if sig and sig != "n.s.":
        title_str = f"{roi} {sig}"
    else:
        title_str = f"{roi} (n.s.)"
    if tissue == "subcortical":
        ax.set_title(title_str, fontsize=title_fontsize,
                     fontstyle="italic", pad=4)
    else:
        ax.set_title(title_str, fontsize=title_fontsize, pad=4)

    # d_adj annotation in upper-right of panel.
    if show_d_annotation and not pd.isna(d_adj):
        ax.text(
            0.97, 0.95,
            f"$d_{{adj}}$ = {d_adj:.2f}",
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=FONTSIZE_CAPTION,
            color=COLOR_GRAY_DARK,
        )

    # X-axis cap.
    ax.set_xlim(0.0, 5.7)
    ax.set_xticks([0, 3, 5])
    if show_xlabel:
        ax.set_xlabel("Years from baseline", fontsize=FONTSIZE_TICK)
    if show_ylabel:
        if tissue == "subcortical":
            ax.set_ylabel("Volume x 10^3 / eTIV",
                          fontsize=FONTSIZE_TICK)
        else:
            ax.set_ylabel("Cortical thickness (mm)",
                          fontsize=FONTSIZE_TICK)


# Figure 2: FDR-significant trajectories.

def plot_figure2_trajectory_fdrsig(
    cohort_long: pd.DataFrame,
    test_b: pd.DataFrame,
    ordering: pd.DataFrame,
    out_path: Path,
    show_individual_lines: bool = True,
) -> None:
    """Trajectory panel for the FDR-significant ROIs.

    Layout: 2x2 grid for 4 FDR-sig ROIs. Layout adapts if the FDR-sig
    set is a different size (1xN for <=3, 2xceil for 4-6, 3xceil for
    7+). Wong palette; thin per-subject lines under bold LME-fitted
    group means with 95% CI ribbons.
    """
    _apply_figure_style()

    # FDR-sig set.
    sig_mask = test_b["sig_star_tissue"].isin(["*", "**", "***"])
    sig_rois = test_b[sig_mask].copy()
    if sig_rois.empty:
        print("  Skipping Figure 2: no FDR-significant ROIs to plot.")
        return

    # Order by canonical |d_adj_5y| ranking from `ordering`.
    ordering_idx = ordering.set_index("ROI")
    sig_rois["ordering_rank"] = sig_rois["ROI"].map(
        ordering_idx["ordering_rank"]
    )
    # Tissue-block primary, ordering_rank secondary.
    sig_rois["tissue_priority"] = sig_rois["Tissue"].map(
        {"subcortical": 0, "cortical": 1}
    )
    sig_rois = sig_rois.sort_values(
        ["tissue_priority", "ordering_rank"]
    ).reset_index(drop=True)

    n_sig = len(sig_rois)
    # Layout by count.
    if n_sig <= 3:
        n_rows, n_cols = 1, n_sig
    elif n_sig <= 6:
        n_rows, n_cols = 2, math.ceil(n_sig / 2)
    else:
        n_rows, n_cols = 3, math.ceil(n_sig / 3)

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=FIG2_TRAJ_FDRSIG_DIMS,
        sharex=True, squeeze=False,
    )

    # Render each panel.
    for i, sig_row in sig_rois.iterrows():
        r, c = divmod(i, n_cols)
        ax = axes[r, c]
        roi = sig_row["ROI"]
        roi_data = cohort_long[cohort_long["ROI"] == roi]
        _plot_trajectory_panel(
            ax=ax, roi_data=roi_data, roi=roi,
            test_b_row=sig_row,
            show_d_annotation=True,
            line_alpha=0.20, line_width=0.7,
            show_xlabel=(r == n_rows - 1),  # bottom row
            show_ylabel=(c == 0),            # left column
            show_individual_lines=show_individual_lines,
        )

    # Hide unused axes.
    for i in range(n_sig, n_rows * n_cols):
        r, c = divmod(i, n_cols)
        axes[r, c].set_visible(False)

    # Figure-level title and subtitle (subtitle adapts to clean variant).
    fig.suptitle(
        "PD vs. Control trajectories: FDR-significant ROIs",
        fontsize=FONTSIZE_TITLE, y=0.99,
    )
    if show_individual_lines:
        subtitle_text = (
            "Thin lines = individual subjects; bold = LME-fitted "
            "group mean; ribbon = 95% CI. "
            "Y-axis scales differ by panel."
        )
    else:
        subtitle_text = (
            "Bold = LME-fitted group mean; ribbon = 95% CI "
            "(individual-subject trajectories suppressed for clarity). "
            "Y-axis scales differ by panel."
        )
    fig.text(
        0.5, 0.945,
        subtitle_text,
        ha="center", va="center",
        fontsize=FONTSIZE_SUBTITLE,
        color=COLOR_GRAY_DARK, style="italic",
    )

    # Legend below subtitle so it doesn't collide with panel titles.
    legend_handles = [
        Line2D([0], [0], color=COLOR_PD, linewidth=2.5,
               label="PD (n=67)"),
        Line2D([0], [0], color=COLOR_CONTROL, linewidth=2.5,
               label="Control (n=33)"),
    ]
    fig.legend(
        handles=legend_handles, loc="upper center",
        bbox_to_anchor=(0.5, 0.92),
        frameon=False, fontsize=FONTSIZE_ANNOTATION,
        ncol=2, columnspacing=2.0,
    )

    fig.subplots_adjust(top=0.84, bottom=0.10, left=0.10,
                        right=0.97, hspace=0.35, wspace=0.32)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=FIG_DPI)
    plt.close(fig)
    print(f"  Wrote: {out_path}")


# Figure 3 / 4: tissue trajectory grids.

def _plot_tissue_trajectory_grid(
    cohort_long: pd.DataFrame,
    test_b: pd.DataFrame,
    ordering: pd.DataFrame,
    out_path: Path,
    tissue: str,
    title: str,
    show_individual_lines: bool = True,
) -> None:
    """Render a 2x3 trajectory grid for the 6 ROIs of one tissue.

    Shared by plot_figure3_trajectory_subcortical and
    plot_figure4_trajectory_cortical. Each panel uses raw units with
    per-panel auto-scaled y-axes; the d_adj annotation handles
    cross-panel slope comparison.
    """
    _apply_figure_style()

    # ROI order: 6 ROIs of this tissue, ranked by |d_adj_5y|.
    block = ordering[ordering["Tissue"] == tissue].sort_values(
        "ordering_rank"
    )
    rois = list(block["ROI"])
    test_b_idx = test_b.set_index("ROI")

    n_rows, n_cols = 2, 3
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(7.5, 5.0),
        sharex=True, squeeze=False,
    )

    for i, roi in enumerate(rois):
        r, c = divmod(i, n_cols)
        ax = axes[r, c]
        roi_data = cohort_long[cohort_long["ROI"] == roi]
        sig_row = test_b_idx.loc[roi]
        _plot_trajectory_panel(
            ax=ax, roi_data=roi_data, roi=roi,
            test_b_row=sig_row,
            show_d_annotation=True,
            line_alpha=0.20,
            line_width=0.6,
            show_xlabel=(r == n_rows - 1),
            show_ylabel=(c == 0),
            title_fontsize=FONTSIZE_LABEL,
            show_individual_lines=show_individual_lines,
        )

    # Hide unused axes (shouldn't happen with 6 ROIs in 2x3).
    for i in range(len(rois), n_rows * n_cols):
        r, c = divmod(i, n_cols)
        axes[r, c].set_visible(False)

    # Figure-level title and subtitle (adapts to clean variant).
    fig.suptitle(title, fontsize=FONTSIZE_TITLE, y=0.985)
    if show_individual_lines:
        subtitle_text = (
            "Thin lines = individual subjects (real elapsed time); "
            "bold = LME-fitted group mean; ribbon = 95% CI. "
            "Y-axis scales differ by panel."
        )
    else:
        subtitle_text = (
            "Bold = LME-fitted group mean; ribbon = 95% CI "
            "(individual-subject trajectories suppressed for clarity). "
            "Y-axis scales differ by panel."
        )
    fig.text(
        0.5, 0.943,
        subtitle_text,
        ha="center", va="center",
        fontsize=FONTSIZE_SUBTITLE,
        color=COLOR_GRAY_DARK, style="italic",
    )

    # Legend below subtitle.
    legend_handles = [
        Line2D([0], [0], color=COLOR_PD, linewidth=2.2,
               label="PD (n=67)"),
        Line2D([0], [0], color=COLOR_CONTROL, linewidth=2.2,
               label="Control (n=33)"),
    ]
    fig.legend(
        handles=legend_handles, loc="upper center",
        bbox_to_anchor=(0.5, 0.915),
        frameon=False, fontsize=FONTSIZE_ANNOTATION,
        ncol=2, columnspacing=2.0,
    )

    fig.subplots_adjust(top=0.83, bottom=0.10, left=0.10,
                        right=0.97, hspace=0.40, wspace=0.32)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=FIG_DPI)
    plt.close(fig)
    print(f"  Wrote: {out_path}")


def plot_figure3_trajectory_subcortical(
    cohort_long: pd.DataFrame,
    test_b: pd.DataFrame,
    ordering: pd.DataFrame,
    out_path: Path,
    show_individual_lines: bool = True,
) -> None:
    """Supplementary trajectory grid: 6 subcortical ROIs in a 2x3
    grid, ranked by |d_adj_5y|. Italic ROI titles. Y-axis in raw
    eTIV-normalized volume units (x 10^3), per-panel auto-scaled.

    show_individual_lines=False produces the clean variant.
    """
    _plot_tissue_trajectory_grid(
        cohort_long=cohort_long,
        test_b=test_b,
        ordering=ordering,
        out_path=out_path,
        tissue="subcortical",
        title="PD vs. Control trajectories: subcortical volumes",
        show_individual_lines=show_individual_lines,
    )


def plot_figure4_trajectory_cortical(
    cohort_long: pd.DataFrame,
    test_b: pd.DataFrame,
    ordering: pd.DataFrame,
    out_path: Path,
    show_individual_lines: bool = True,
) -> None:
    """Supplementary trajectory grid: 6 cortical ROIs in a 2x3 grid,
    ranked by |d_adj_5y|. Regular-weight (non-italic) ROI titles.
    Y-axis in cortical thickness (mm), per-panel auto-scaled.

    show_individual_lines=False produces the clean variant.
    """
    _plot_tissue_trajectory_grid(
        cohort_long=cohort_long,
        test_b=test_b,
        ordering=ordering,
        out_path=out_path,
        tissue="cortical",
        title="PD vs. Control trajectories: cortical thickness",
        show_individual_lines=show_individual_lines,
    )


# Figure 5: per-group atrophy heatmap.

def plot_figure5_heatmap_per_group_atrophy(
    cohort_long: pd.DataFrame,
    test_a: pd.DataFrame,
    test_b: pd.DataFrame,
    out_path: Path,
) -> None:
    """Heatmap of LME-predicted within-group % change from baseline at
    the 3-year and 5-year horizons, separately for Control and PD.

    Layout (12 rows x 4 cols). Rows: 12 ROIs in two stacked tissue
    blocks (subcortical then cortical), each sorted by Patient's
    |predicted_5y_pct| descending. Cols: BL->3Y Control, BL->3Y
    Patient, BL->5Y Control, BL->5Y Patient.

    Cell value: pct(t) = 100 * (slope * t) / mean_bl_group, where
    slope is from Test A's within-group LME and mean_bl_group is the
    simple mean of that group's BL observations.

    Color: diverging colormap centered at 0%, symmetric vmin/vmax
    from the global max |%|. Cell text shows "%value sig" where sig
    is Test A's per-group p_fdr_tissue star (non-sig -> "n.s.").
    """
    _apply_figure_style()

    # Per-group BL means for the % normalization denominator.
    bl_means = (
        cohort_long[cohort_long["Timepoint"] == "BL"]
        .groupby(["ROI", "Group_clinical"])["Value"]
        .mean()
        .reset_index()
        .rename(columns={"Value": "mean_bl"})
    )

    # Join test_a with the BL means.
    ta = test_a.merge(
        bl_means, on=["ROI", "Group_clinical"], how="left"
    )

    # Per-cell predicted % atrophy at t=3 and t=5.
    horizons = (3.0, 5.0)
    cells: List[Dict[str, Any]] = []
    for _, r in ta.iterrows():
        if pd.isna(r["mean_bl"]) or r["mean_bl"] == 0:
            continue
        for t in horizons:
            pct = 100.0 * (r["slope"] * t) / r["mean_bl"]
            cells.append({
                "ROI": r["ROI"],
                "Tissue": r["Tissue"],
                "Group": r["Group_clinical"],   # Pasient or Control
                "Horizon": int(t),
                "pct": pct,
                "p_fdr": r["p_fdr_tissue"],
                "sig": r["sig_star_tissue"],
            })
    cell_df = pd.DataFrame(cells)

    # Row order: subcortical block first, each sorted by Patient's
    # |predicted_5y_pct| from test_b (descending). Derive directly
    # from the per-group t=5 Pasient cell value:
    sort_key = (
        test_b[["ROI", "Tissue", "predicted_5y_pct_pd"]]
        if "predicted_5y_pct_pd" in test_b.columns
        else None
    )
    pat_t5 = cell_df[
        (cell_df["Group"] == "Pasient") & (cell_df["Horizon"] == 5)
    ][["ROI", "Tissue", "pct"]].rename(columns={"pct": "patient_5y_pct"})
    pat_t5["abs_pct"] = pat_t5["patient_5y_pct"].abs()

    sub_rois = list(
        pat_t5[pat_t5["Tissue"] == "subcortical"]
        .sort_values("abs_pct", ascending=False)["ROI"]
    )
    cor_rois = list(
        pat_t5[pat_t5["Tissue"] == "cortical"]
        .sort_values("abs_pct", ascending=False)["ROI"]
    )
    roi_order = sub_rois + cor_rois
    n_sub = len(sub_rois)
    n_total = len(roi_order)

    # Cell matrices: rows = ROIs, cols = (Co_3y, PD_3y, Co_5y, PD_5y).
    col_specs = [
        ("Control", 3), ("Pasient", 3),
        ("Control", 5), ("Pasient", 5),
    ]
    col_labels = [
        "BL→3Y\nControl", "BL→3Y\nPatient",
        "BL→5Y\nControl", "BL→5Y\nPatient",
    ]
    n_cols = len(col_specs)

    pct_matrix = np.full((n_total, n_cols), np.nan)
    sig_matrix: List[List[str]] = [
        ["n.s." for _ in range(n_cols)] for _ in range(n_total)
    ]
    for i, roi in enumerate(roi_order):
        for j, (group, horizon) in enumerate(col_specs):
            cell = cell_df[
                (cell_df["ROI"] == roi)
                & (cell_df["Group"] == group)
                & (cell_df["Horizon"] == horizon)
            ]
            if cell.empty:
                continue
            pct_matrix[i, j] = float(cell["pct"].iloc[0])
            sig_matrix[i][j] = str(cell["sig"].iloc[0])

    # Symmetric vmin/vmax based on max |%| (global).
    vmax_abs = float(np.nanmax(np.abs(pct_matrix)))
    if not np.isfinite(vmax_abs) or vmax_abs == 0:
        vmax_abs = 5.0
    # Pad slightly so cell text isn't lost in the most-saturated cells.
    vmax_abs = vmax_abs * 1.1

    # Figure size: row height ~0.45in plus margin.
    fig_h = max(5.5, 0.45 * n_total + 2.0)
    fig, ax = plt.subplots(figsize=(8.0, fig_h))

    # Heatmap.
    im = ax.imshow(
        pct_matrix, cmap="RdBu",
        vmin=-vmax_abs, vmax=vmax_abs,
        aspect="auto",
    )

    # Cell text annotations.
    for i in range(n_total):
        for j in range(n_cols):
            val = pct_matrix[i, j]
            if not np.isfinite(val):
                continue
            sig = sig_matrix[i][j]
            sig_text = sig if sig != "n.s." else "n.s."
            label = f"{val:+.2f}%\n{sig_text}"
            # Dark text on light bg, light text on dark bg.
            normalized = abs(val) / vmax_abs
            text_color = "white" if normalized > 0.55 else "black"
            ax.text(
                j, i, label,
                ha="center", va="center",
                fontsize=FONTSIZE_CAPTION,
                color=text_color, linespacing=1.1,
            )

    # Tissue divider between subcortical and cortical blocks.
    if 0 < n_sub < n_total:
        ax.axhline(
            y=n_sub - 0.5, color=COLOR_GRAY_DARK,
            linewidth=1.2, alpha=0.7, zorder=4,
        )

    # Axis labels: ROI on rows, group x horizon on cols.
    ax.set_yticks(np.arange(n_total))
    ax.set_yticklabels(roi_order, fontsize=FONTSIZE_TICK)
    # Italic for subcortical.
    for tick_label, roi in zip(ax.get_yticklabels(), roi_order):
        tissue = (
            "subcortical" if roi in sub_rois else "cortical"
        )
        if tissue == "subcortical":
            tick_label.set_style("italic")

    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(col_labels, fontsize=FONTSIZE_TICK)
    ax.tick_params(axis="x", which="both", length=0, pad=4)
    ax.tick_params(axis="y", which="both", length=0)

    # Light vertical separator between BL->3Y and BL->5Y blocks.
    ax.axvline(
        x=1.5, color=COLOR_GRAY_DARK,
        linewidth=0.8, alpha=0.4, zorder=4,
    )

    # Title and colorbar.
    ax.set_title(
        "LME-predicted within-group % change from baseline",
        fontsize=FONTSIZE_TITLE, pad=14,
    )

    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cbar.set_label("% change from baseline",
                   fontsize=FONTSIZE_TICK, rotation=270, labelpad=14)
    cbar.ax.tick_params(labelsize=FONTSIZE_CAPTION)

    fig.subplots_adjust(top=0.92, bottom=0.10, left=0.20, right=0.95)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=FIG_DPI)
    plt.close(fig)
    print(f"  Wrote: {out_path}")


# Risk-A scatter (supplementary).

def plot_figure_risk_a_scatter(
    risk_a: pd.DataFrame,
    out_path: Path,
) -> None:
    """Supplementary Risk-A scatter: harmonized vs unharmonized
    beta_interaction for all 12 ROIs, split into tissue-specific panels.

    x-axis: beta_unharmonized (per year); y-axis: beta_harmonized.
    One point per ROI; marker shape encodes tissue, color encodes
    ROI within tissue. Axes are scaled separately per tissue.

    Harmonization sensitivity check, not primary inference.
    """
    _apply_figure_style()

    fig, axes = plt.subplots(
        1, 2,
        figsize=(10.2, 5.6),
        sharex=False,
        sharey=False,
        squeeze=False,
    )
    axes = axes[0]

    panel_specs = [
        {
            "tissue": "subcortical",
            "title": "Subcortical volume",
            "marker": "o",
            "color_map": RISK_A_SUBCORTICAL_COLORS,
            "roi_order": RISK_A_SUBCORTICAL_ORDER,
            "ticks": [-2e-5, -1e-5, 0.0, 1e-5, 2e-5],
            "formatter": None,
        },
        {
            "tissue": "cortical",
            "title": "Cortical thickness",
            "marker": "s",
            "color_map": RISK_A_CORTICAL_COLORS,
            "roi_order": RISK_A_CORTICAL_ORDER,
            "ticks": [-0.010, -0.005, 0.000, 0.005, 0.010],
            "formatter": FormatStrFormatter("%.3f"),
        },
    ]

    for ax, spec in zip(axes, panel_specs):
        tissue = spec["tissue"]
        marker = spec["marker"]
        color_map = spec["color_map"]
        available = set(risk_a.loc[risk_a["Tissue"] == tissue, "ROI"])
        roi_order = [roi for roi in spec["roi_order"] if roi in available]
        # Include any unexpected future ROI at the end rather than dropping it.
        present_rois = list(risk_a.loc[risk_a["Tissue"] == tissue, "ROI"])
        roi_order.extend([roi for roi in present_rois if roi not in roi_order])

        sub = risk_a[risk_a["Tissue"] == tissue].copy()
        x = sub["beta_unharmonized"].to_numpy(dtype=float)
        y = sub["beta_harmonized"].to_numpy(dtype=float)

        finite = np.isfinite(x) & np.isfinite(y)
        if finite.any():
            bound = np.nanmax(np.abs(np.concatenate([x[finite], y[finite]])))
            if not np.isfinite(bound) or bound == 0:
                bound = 1e-5
            bound *= 1.20
        else:
            bound = 1e-5

        # Keep fixed ticks visible even if future values are smaller.
        ticks = spec["ticks"]
        if ticks:
            bound = max(bound, max(abs(t) for t in ticks) * 1.05)

        # Reference lines.
        ax.axhline(
            0.0, color=COLOR_GRAY_DARK,
            linewidth=0.6, alpha=0.6, zorder=1,
        )
        ax.axvline(
            0.0, color=COLOR_GRAY_DARK,
            linewidth=0.6, alpha=0.6, zorder=1,
        )

        # Diagonal y = x.
        diag = np.linspace(-bound, bound, 100)
        ax.plot(
            diag, diag,
            color=COLOR_GRAY_LIGHT,
            linewidth=0.8,
            linestyle="--",
            alpha=0.8,
            zorder=1,
        )

        # Plot one ROI at a time so color encodes ROI identity and the
        # legend identifies points without inline annotations.
        legend_handles = []
        sub_idx = sub.set_index("ROI")
        for roi in roi_order:
            if roi not in sub_idx.index:
                continue
            r = sub_idx.loc[roi]
            if isinstance(r, pd.DataFrame):
                r = r.iloc[0]
            if pd.isna(r["beta_harmonized"]) or pd.isna(r["beta_unharmonized"]):
                continue

            colour = color_map.get(roi, COLOR_NEUTRAL)
            sig = r["sig_star_tissue_harmonized"]
            is_sig = sig in ("*", "**", "***")
            legend_label = f"{roi} {sig}" if is_sig else roi

            ax.scatter(
                r["beta_unharmonized"],
                r["beta_harmonized"],
                s=78,
                c=colour,
                edgecolors=COLOR_GRAY_DARK,
                linewidths=0.8,
                marker=marker,
                zorder=4,
            )
            legend_handles.append(
                Line2D(
                    [0], [0], marker=marker, linestyle="",
                    markerfacecolor=colour,
                    markeredgecolor=COLOR_GRAY_DARK,
                    markeredgewidth=0.8,
                    markersize=6.5,
                    label=legend_label,
                )
            )

        ax.set_xlim(-bound, bound)
        ax.set_ylim(-bound, bound)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(spec["title"], fontsize=FONTSIZE_LABEL, pad=8)

        if tissue == "subcortical":
            ax.set_xlabel(
                r"$\beta_{\mathrm{interaction}}$ "
                r"(unharmonized, $\times 10^{-5}$ yr$^{-1}$)"
            )
            ax.set_ylabel(
                r"$\beta_{\mathrm{interaction}}$ "
                r"(harmonized, $\times 10^{-5}$ yr$^{-1}$)"
            )
        else:
            ax.set_xlabel(
                r"$\beta_{\mathrm{interaction}}$ "
                r"(unharmonized, yr$^{-1}$)"
            )
            ax.set_ylabel(
                r"$\beta_{\mathrm{interaction}}$ "
                r"(harmonized, yr$^{-1}$)"
            )

        # Controlled tick density prevents cortical x-axis label overlap.
        # Subcortical ticks are displayed as multiples of 10^-5 to avoid
        # Matplotlib's offset text colliding with the x-axis label.
        if ticks:
            ax.xaxis.set_major_locator(FixedLocator(ticks))
            ax.yaxis.set_major_locator(FixedLocator(ticks))
        formatter = spec["formatter"]
        if tissue == "subcortical":
            sub_formatter = FuncFormatter(lambda val, _pos: f"{val / 1e-5:.0f}")
            ax.xaxis.set_major_formatter(sub_formatter)
            ax.yaxis.set_major_formatter(sub_formatter)
        elif formatter is not None:
            ax.xaxis.set_major_formatter(formatter)
            ax.yaxis.set_major_formatter(formatter)

        ax.tick_params(axis="both", labelsize=FONTSIZE_CAPTION)

        ax.legend(
            handles=legend_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.23),
            frameon=False,
            fontsize=FONTSIZE_CAPTION,
            ncol=2,
            columnspacing=1.0,
            handletextpad=0.45,
            borderaxespad=0.0,
        )

    fig.suptitle(
        "Risk-A: harmonized vs unharmonized Group x Time interaction",
        fontsize=FONTSIZE_TITLE,
        y=0.985,
    )


    fig.subplots_adjust(
        top=0.83,
        bottom=0.31,
        left=0.08,
        right=0.98,
        wspace=0.32,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Wrote: {out_path}")


# --- LME diagnostics (Test B) ---

def _safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson correlation guarded against degenerate inputs."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return np.nan
    if np.nanstd(x[mask]) == 0 or np.nanstd(y[mask]) == 0:
        return np.nan
    return float(np.corrcoef(x[mask], y[mask])[0, 1])


def compute_lme_diagnostics_test_b(
    cohort_long: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Per-ROI Test B model diagnostics.

    Refits the same Test B LME formula and extracts convergence
    status, variance components, ICC, max |standardized residual|
    (and count > 3), residual mean/sd/skew/excess kurtosis, and
    correlations of residual with fitted/years. Diagnostic flags
    are concatenated into a semicolon-separated string.

    Returns (diagnostics_df, residuals_df). residuals_df has one
    row per (ROI, scan) for the residual-vs-fitted plot.
    """
    formula = (
        "Value ~ C(Group_clinical, Treatment(reference='Control'))"
        " * Years_from_BL"
        " + age_at_BL"
        " + C(PatientSex_clinical)"
    )

    diag_rows: List[Dict[str, Any]] = []
    resid_rows: List[Dict[str, Any]] = []

    rois = sorted(cohort_long["ROI"].unique())

    for roi in rois:
        roi_data = cohort_long[cohort_long["ROI"] == roi].copy()
        tissue = roi_data["Tissue"].iloc[0]

        n_obs = len(roi_data)
        n_subjects = roi_data["Subject"].nunique()

        # Fit
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                md = smf.mixedlm(
                    formula, data=roi_data,
                    groups=roi_data["Subject"],
                )
                mdf = md.fit(reml=True, method="lbfgs")
        except (np.linalg.LinAlgError, ValueError, OverflowError) as exc:
            diag_rows.append({
                "ROI": roi, "Tissue": tissue,
                "n_obs": n_obs, "n_subjects": n_subjects,
                "status": f"diagnostic_fit_failed:{type(exc).__name__}",
                "converged": False,
                "residual_var": np.nan,
                "random_intercept_var": np.nan,
                "icc": np.nan,
                "max_abs_standardized_residual": np.nan,
                "n_abs_standardized_residual_gt_3": np.nan,
                "residual_mean": np.nan,
                "residual_sd": np.nan,
                "residual_skew": np.nan,
                "residual_kurtosis_excess": np.nan,
                "corr_resid_fitted": np.nan,
                "corr_abs_resid_fitted": np.nan,
                "corr_resid_years": np.nan,
                "diagnostic_flags": (
                    f"diagnostic_fit_failed:{type(exc).__name__}"
                ),
            })
            continue

        converged = bool(mdf.converged)
        status = "ok" if converged else "convergence_failed"

        # Variance components and ICC.
        residual_var = float(mdf.scale)
        residual_sd_total = float(np.sqrt(residual_var)) if residual_var > 0 else np.nan
        try:
            random_intercept_var = float(mdf.cov_re.iloc[0, 0])
        except Exception:
            random_intercept_var = np.nan
        if (
            not np.isnan(random_intercept_var)
            and not np.isnan(residual_var)
            and (random_intercept_var + residual_var) > 0
        ):
            icc = random_intercept_var / (random_intercept_var + residual_var)
            if random_intercept_var < 1e-20:
                status = "singular_fit" if status == "ok" else status
        else:
            icc = np.nan

        # Residuals and fitted values.
        fitted = np.asarray(mdf.fittedvalues, dtype=float)
        observed = np.asarray(roi_data["Value"], dtype=float)
        residual = observed - fitted

        if not np.isnan(residual_sd_total) and residual_sd_total > 0:
            std_residual = residual / residual_sd_total
        else:
            std_residual = np.full_like(residual, np.nan)

        # Residual-distribution stats.
        residual_mean = float(np.nanmean(residual))
        residual_sd = float(np.nanstd(residual, ddof=1)) if residual.size > 1 else np.nan
        try:
            residual_skew = float(scipy_stats.skew(residual, bias=False))
        except Exception:
            residual_skew = np.nan
        try:
            residual_kurt = float(
                scipy_stats.kurtosis(residual, fisher=True, bias=False)
            )
        except Exception:
            residual_kurt = np.nan

        # Standardized-residual extremes.
        if np.all(np.isnan(std_residual)):
            max_abs_std = np.nan
            n_gt_3 = np.nan
        else:
            max_abs_std = float(np.nanmax(np.abs(std_residual)))
            n_gt_3 = int(np.nansum(np.abs(std_residual) > 3.0))

        # Correlations.
        years = np.asarray(roi_data["Years_from_BL"], dtype=float)
        corr_rf = _safe_corr(residual, fitted)
        corr_arf = _safe_corr(np.abs(std_residual), fitted)
        corr_ry = _safe_corr(residual, years)

        # Flags.
        flags: List[str] = []
        if status != "ok":
            flags.append(status)
        if not np.isnan(max_abs_std) and max_abs_std > 4.0:
            flags.append("max_std_resid_gt_4")
        if not np.isnan(n_gt_3) and n_gt_3 > 3:
            flags.append("n_std_resid_gt_3_more_than_3")
        if not np.isnan(residual_skew) and abs(residual_skew) > 1.0:
            flags.append("residual_skew_gt_1")
        if not np.isnan(residual_kurt) and residual_kurt > 4.0:
            flags.append("residual_kurtosis_excess_gt_4")
        if not np.isnan(corr_rf) and abs(corr_rf) > 0.30:
            flags.append("corr_resid_fitted_gt_0.30")
        if not np.isnan(corr_arf) and abs(corr_arf) > 0.30:
            flags.append("corr_abs_resid_fitted_gt_0.30")
        if not np.isnan(corr_ry) and abs(corr_ry) > 0.30:
            flags.append("corr_resid_years_gt_0.30")

        diag_rows.append({
            "ROI": roi, "Tissue": tissue,
            "n_obs": n_obs, "n_subjects": n_subjects,
            "status": status, "converged": converged,
            "residual_var": residual_var,
            "random_intercept_var": random_intercept_var,
            "icc": icc,
            "max_abs_standardized_residual": max_abs_std,
            "n_abs_standardized_residual_gt_3": n_gt_3,
            "residual_mean": residual_mean,
            "residual_sd": residual_sd,
            "residual_skew": residual_skew,
            "residual_kurtosis_excess": residual_kurt,
            "corr_resid_fitted": corr_rf,
            "corr_abs_resid_fitted": corr_arf,
            "corr_resid_years": corr_ry,
            "diagnostic_flags": ";".join(flags) if flags else "none",
        })

        # Per-row residuals for the plot.
        for i in range(len(roi_data)):
            row = roi_data.iloc[i]
            resid_rows.append({
                "ROI": roi, "Tissue": tissue,
                "Subject": row["Subject"],
                "Timepoint": row["Timepoint"],
                "Group_clinical": row["Group_clinical"],
                "Years_from_BL": float(row["Years_from_BL"]),
                "Value": float(row["Value"]),
                "fitted": float(fitted[i]),
                "residual": float(residual[i]),
                "std_residual": (
                    float(std_residual[i])
                    if not np.isnan(std_residual[i]) else np.nan
                ),
            })

    diag_df = pd.DataFrame(diag_rows)
    resid_df = pd.DataFrame(resid_rows)
    return diag_df, resid_df


def write_lme_diagnostics_csv(df: pd.DataFrame, path: Path) -> None:
    """Write lme_diagnostics_test_b.csv with stable column order."""
    col_order = [
        "ROI", "Tissue",
        "n_obs", "n_subjects",
        "status", "converged",
        "residual_var", "random_intercept_var", "icc",
        "max_abs_standardized_residual",
        "n_abs_standardized_residual_gt_3",
        "residual_mean", "residual_sd",
        "residual_skew", "residual_kurtosis_excess",
        "corr_resid_fitted", "corr_abs_resid_fitted",
        "corr_resid_years",
        "diagnostic_flags",
    ]
    out = df[col_order].copy()
    out.to_csv(path, index=False, float_format=CSV_FLOAT_FORMAT)
    print(f"  Wrote: {path} ({len(out)} rows, {len(out.columns)} cols)")


def plot_lme_diagnostics_test_b(
    residuals_df: pd.DataFrame,
    out_path: Path,
    ordering: Optional[pd.DataFrame] = None,
) -> None:
    """Supplementary diagnostic figure: standardized residuals vs
    fitted values, one panel per Test B ROI.

    Layout: 3x4 grid. Reference lines at 0 and +/-3 standardized
    residuals. Subcortical titles in italic.

    If `ordering` is provided, panels follow the canonical |d_adj_5y|
    order; otherwise sorted by tissue then ROI.
    """
    _apply_figure_style()

    if ordering is not None and not ordering.empty:
        sub = ordering[ordering["Tissue"] == "subcortical"].sort_values(
            "ordering_rank"
        )
        cor = ordering[ordering["Tissue"] == "cortical"].sort_values(
            "ordering_rank"
        )
        roi_order = list(sub["ROI"]) + list(cor["ROI"])
    else:
        roi_order = sorted(residuals_df["ROI"].unique(),
                           key=lambda r: (
                               residuals_df.loc[residuals_df["ROI"] == r,
                                                "Tissue"].iloc[0],
                               r,
                           ))

    n_rois = len(roi_order)
    n_cols = 4
    n_rows = int(np.ceil(n_rois / n_cols))

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(11.0, 2.5 * n_rows + 0.5),
        sharex=False, sharey=True,
    )
    axes = np.atleast_2d(axes)

    for idx, roi in enumerate(roi_order):
        r = idx // n_cols
        c = idx % n_cols
        ax = axes[r, c]
        sub = residuals_df[residuals_df["ROI"] == roi]

        if sub.empty:
            ax.text(0.5, 0.5, "no data",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=FONTSIZE_ANNOTATION,
                    color=COLOR_GRAY_DARK)
            ax.set_xticks([]); ax.set_yticks([])
            continue

        tissue = sub["Tissue"].iloc[0]

        # Reference lines.
        ax.axhline(0.0, color=COLOR_GRAY_DARK,
                   linewidth=0.7, alpha=0.7, zorder=1)
        for ref_y in (-3.0, 3.0):
            ax.axhline(ref_y, color=COLOR_GRAY_LIGHT,
                       linewidth=0.6, linestyle="--",
                       alpha=0.8, zorder=1)

        # Color by group.
        for grp, color in [
            ("Pasient", COLOR_PD), ("Control", COLOR_CONTROL)
        ]:
            mask = sub["Group_clinical"] == grp
            ax.scatter(
                sub.loc[mask, "fitted"],
                sub.loc[mask, "std_residual"],
                s=14, c=color, alpha=0.55,
                edgecolors="none", zorder=3,
            )

        max_abs_resid = sub["std_residual"].abs().max(skipna=True)
        if pd.notna(max_abs_resid):
            ax.text(
                0.98,
                0.04,
                rf"max $|r_z|$={max_abs_resid:.2f}",
                ha="right",
                va="bottom",
                transform=ax.transAxes,
                fontsize=FONTSIZE_CAPTION,
                color=COLOR_GRAY_DARK,
                alpha=0.9,
            )

        title = roi
        if tissue == "subcortical":
            ax.set_title(title, fontsize=FONTSIZE_SUBTITLE, fontstyle="italic")
        else:
            ax.set_title(title, fontsize=FONTSIZE_SUBTITLE)

        if c == 0:
            ax.set_ylabel("Std. residual",
                          fontsize=FONTSIZE_ANNOTATION)
        if r == n_rows - 1:
            ax.set_xlabel("Fitted",
                          fontsize=FONTSIZE_ANNOTATION)

        ax.tick_params(axis="both", labelsize=FONTSIZE_CAPTION)

        # Cap tick density so x-axis labels stay readable in the
        # compact 3x4 layout (Caudate etc. otherwise get long
        # decimal labels from the default locator).
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))

        # Subcortical fitted values are small normalized volumes;
        # scientific notation keeps labels short. Cortical thickness
        # panels stay in plain decimal.
        fitted_abs_max = sub["fitted"].abs().max(skipna=True)
        if pd.notna(fitted_abs_max) and fitted_abs_max < 0.01:
            ax.ticklabel_format(
                axis="x",
                style="sci",
                scilimits=(0, 0),
                useOffset=False,
            )
        else:
            ax.ticklabel_format(
                axis="x",
                style="plain",
                useOffset=False,
            )
        ax.xaxis.get_offset_text().set_fontsize(FONTSIZE_CAPTION)

    # Hide unused subplots.
    for idx in range(n_rois, n_rows * n_cols):
        r = idx // n_cols
        c = idx % n_cols
        axes[r, c].axis("off")

    # Legend on the figure.
    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="",
               markerfacecolor=COLOR_PD, markeredgecolor="none",
               markersize=7, label="PD"),
        Line2D([0], [0], marker="o", linestyle="",
               markerfacecolor=COLOR_CONTROL, markeredgecolor="none",
               markersize=7, label="Control"),
    ]

    fig.suptitle(
        "Test B residual diagnostics: std. residual vs fitted",
        fontsize=FONTSIZE_TITLE, y=0.995,
    )
    fig.legend(
        handles=legend_handles,
        loc="upper right", frameon=False,
        fontsize=FONTSIZE_ANNOTATION,
        bbox_to_anchor=(0.99, 0.985),
    )


    fig.tight_layout(rect=(0, 0.04, 1, 0.96))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=FIG_DPI)
    plt.close(fig)
    print(f"  Wrote: {out_path}")


# --- Demographics ---

def compute_demographics(cohort_long: pd.DataFrame) -> pd.DataFrame:
    """Cohort demographics summary.

    Reduces cohort_long to one row per subject (BL of one ROI) and
    computes descriptives per group. Returns columns variable, PD,
    Control with 9 rows: n, age_at_BL_mean, age_at_BL_sd, sex_female,
    sex_male, follow_up_3Y_mean, follow_up_3Y_sd, follow_up_5Y_mean,
    follow_up_5Y_sd.
    """
    any_roi = cohort_long["ROI"].iloc[0]
    bl_slice = cohort_long[
        (cohort_long["Timepoint"] == "BL")
        & (cohort_long["ROI"] == any_roi)
    ].copy()
    pd_subj = bl_slice[bl_slice["Group_clinical"] == "Pasient"]
    co_subj = bl_slice[bl_slice["Group_clinical"] == "Control"]

    pd_age = pd_subj["age_at_BL"].dropna()
    co_age = co_subj["age_at_BL"].dropna()

    pd_f = int((pd_subj["PatientSex_clinical"] == "F").sum())
    pd_m = int((pd_subj["PatientSex_clinical"] == "M").sum())
    co_f = int((co_subj["PatientSex_clinical"] == "F").sum())
    co_m = int((co_subj["PatientSex_clinical"] == "M").sum())

    fu3 = cohort_long[
        (cohort_long["Timepoint"] == "3Y")
        & (cohort_long["ROI"] == any_roi)
    ]
    fu5 = cohort_long[
        (cohort_long["Timepoint"] == "5Y")
        & (cohort_long["ROI"] == any_roi)
    ]
    pd_fu3 = fu3[fu3["Group_clinical"] == "Pasient"]["Years_from_BL"]
    co_fu3 = fu3[fu3["Group_clinical"] == "Control"]["Years_from_BL"]
    pd_fu5 = fu5[fu5["Group_clinical"] == "Pasient"]["Years_from_BL"]
    co_fu5 = fu5[fu5["Group_clinical"] == "Control"]["Years_from_BL"]

    rows = [
        ("n", str(len(pd_subj)), str(len(co_subj))),
        ("age_at_BL_mean", f"{pd_age.mean():.1f}", f"{co_age.mean():.1f}"),
        ("age_at_BL_sd", f"{pd_age.std():.1f}", f"{co_age.std():.1f}"),
        ("sex_female", str(pd_f), str(co_f)),
        ("sex_male", str(pd_m), str(co_m)),
        ("follow_up_3Y_mean", f"{pd_fu3.mean():.2f}", f"{co_fu3.mean():.2f}"),
        ("follow_up_3Y_sd", f"{pd_fu3.std():.2f}", f"{co_fu3.std():.2f}"),
        ("follow_up_5Y_mean", f"{pd_fu5.mean():.2f}", f"{co_fu5.mean():.2f}"),
        ("follow_up_5Y_sd", f"{pd_fu5.std():.2f}", f"{co_fu5.std():.2f}"),
    ]
    return pd.DataFrame(rows, columns=["variable", "PD", "Control"])


def write_demographics_csv(df: pd.DataFrame, path: Path) -> None:
    """Write demographics.csv."""
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    print(f"  Wrote: {path} ({len(df)} rows)")


# --- Main ---

def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--cohort-long", type=Path, required=True,
                   help="Path to cohort_long_harmonized.csv. The "
                        "Risk-A sensitivity check expects "
                        "cohort_long_unharmonized.csv as a sibling "
                        "file in the same directory.")
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    cohort_long = load_long_harmonized(args.cohort_long)

    # Test B.
    print("\nFitting Test B (between-group, 12 ROIs)...")
    test_b = run_test_b(cohort_long)
    n_ok = (test_b["status"] == "ok").sum()
    n_sing = (test_b["status"] == "singular_fit").sum()
    n_failed = (
        test_b["status"].str.startswith("convergence_failed").sum()
    )
    print(f"  Fitted: {n_ok} ok, {n_sing} singular, {n_failed} failed")

    n_sig_tissue = (test_b["p_fdr_tissue"] < 0.05).sum()
    print(f"  FDR within tissue: {n_sig_tissue} of 12 ROIs < 0.05")

    if n_ok > 0:
        ordered = test_b.copy()
        ordered["abs_d"] = ordered["d_adj_5y"].abs()
        ordered = ordered.sort_values("abs_d", ascending=False)
        print("\n  All 12 ROIs sorted by |d_adj_5y|:")
        for _, r in ordered.iterrows():
            print(f"    {r['ROI']:<25} d={r['d_adj_5y']:>+.3f}  "
                  f"p_fdr={r['p_fdr_tissue']:.4g} {r['sig_star_tissue']}")

    # Test A.
    print("\nFitting Test A (within-group, 24 fits)...")
    test_a = run_test_a(cohort_long)
    n_ok_a = (test_a["status"] == "ok").sum()
    n_sing_a = (test_a["status"] == "singular_fit").sum()
    print(f"  Fitted: {n_ok_a} ok, {n_sing_a} singular, "
          f"{len(test_a) - n_ok_a - n_sing_a} failed")
    n_sig_a = (test_a["p_fdr_tissue"] < 0.05).sum()
    print(f"  {n_sig_a} of {len(test_a)} (Group x ROI) cells "
          f"FDR-significant within (Group x Tissue)")

    # ROI ordering.
    print("\nComputing canonical ROI display order...")
    ordering = compute_roi_ordering(test_b)
    print("  Subcortical block (descending |d_adj_5y|):")
    for _, r in ordering[ordering["Tissue"] == "subcortical"].iterrows():
        print(f"    {r['ordering_rank']+1}. {r['ROI']}")
    print("  Cortical block (descending |d_adj_5y|):")
    for _, r in ordering[ordering["Tissue"] == "cortical"].iterrows():
        print(f"    {r['ordering_rank']+1}. {r['ROI']}")

    # Write CSVs.
    print("\nWriting output CSVs...")
    write_test_b_csv(
        test_b, args.output_dir / "stats_lme_between_group.csv"
    )
    write_test_a_csv(
        test_a, args.output_dir / "stats_lme_within_group.csv"
    )
    write_roi_ordering_csv(
        ordering, args.output_dir / "roi_display_order.csv"
    )

    # Risk-A.
    print("\nRisk-A harmonization check...")
    cohort_long_unharm: Optional[pd.DataFrame] = None
    unharm_path = args.cohort_long.parent / "cohort_long_unharmonized.csv"
    if unharm_path.exists():
        cohort_long_unharm = pd.read_csv(unharm_path)
        print(f"  Loaded {unharm_path.name}: {cohort_long_unharm.shape}")
    else:
        print(f"  WARNING: {unharm_path} not found, skipping Risk-A",
              file=sys.stderr)

    risk_a: Optional[pd.DataFrame] = None
    if cohort_long_unharm is not None:
        risk_a = run_risk_a_all_rois(
            test_b_harmonized=test_b,
            cohort_long_unharm=cohort_long_unharm,
        )
        write_risk_a_csv(
            risk_a, args.output_dir / "risk_a_all_rois.csv"
        )

    write_risk_a_summary_text(
        risk_a=(
            risk_a
            if risk_a is not None
            else pd.DataFrame(columns=[
                "ROI", "Tissue", "beta_harmonized", "beta_unharmonized",
                "beta_abs_ratio", "beta_difference", "sign_match",
                "p_harmonized", "p_unharmonized",
                "d_harmonized", "d_unharmonized", "d_abs_ratio",
                "p_fdr_tissue_harmonized", "sig_star_tissue_harmonized",
            ])
        ),
        summary_path=args.output_dir / "lme_sanity_check.txt",
        test_b_results=test_b,
        cohort_long_unharm_present=(cohort_long_unharm is not None),
    )

    # LME diagnostics for the 12 Test B fits.
    print("\nTest B LME diagnostics...")
    diag_b, resid_b = compute_lme_diagnostics_test_b(cohort_long)
    write_lme_diagnostics_csv(
        diag_b, args.output_dir / "lme_diagnostics_test_b.csv"
    )
    plot_lme_diagnostics_test_b(
        residuals_df=resid_b,
        out_path=args.output_dir / "figure_lme_diagnostics_test_b_supp.png",
        ordering=ordering,
    )

    # Demographics.
    print("\nWriting demographics...")
    demo_df = compute_demographics(cohort_long)
    write_demographics_csv(demo_df, args.output_dir / "demographics.csv")

    # Main-body figures.
    print("\nRendering main-body figures...")
    plot_figure1_forest(
        test_b, ordering,
        args.output_dir / "figure1_forest_d_adj_5Y.png",
    )
    plot_figure2_trajectory_fdrsig(
        cohort_long, test_b, ordering,
        args.output_dir / "figure2_trajectory_fdrsig.png",
    )
    plot_figure2_trajectory_fdrsig(
        cohort_long, test_b, ordering,
        args.output_dir / "figure2_trajectory_fdrsig_clean.png",
        show_individual_lines=False,
    )

    plot_figure5_heatmap_per_group_atrophy(
        cohort_long, test_a, test_b,
        args.output_dir / "figure5_heatmap_per_group_atrophy.png",
    )

    # Supplementary figures.
    print("\nRendering supplementary figures...")
    plot_figure3_trajectory_subcortical(
        cohort_long, test_b, ordering,
        args.output_dir / "figure3_trajectory_subcortical_supp.png",
    )
    plot_figure3_trajectory_subcortical(
        cohort_long, test_b, ordering,
        args.output_dir / "figure3_trajectory_subcortical_supp_clean.png",
        show_individual_lines=False,
    )
    plot_figure4_trajectory_cortical(
        cohort_long, test_b, ordering,
        args.output_dir / "figure4_trajectory_cortical_supp.png",
    )
    plot_figure4_trajectory_cortical(
        cohort_long, test_b, ordering,
        args.output_dir / "figure4_trajectory_cortical_supp_clean.png",
        show_individual_lines=False,
    )

    if risk_a is not None:
        plot_figure_risk_a_scatter(
            risk_a,
            args.output_dir / "figure_risk_a_scatter_supp.png",
        )
    else:
        print(
            "  Skipping figure_risk_a_scatter_supp.png "
            "(no unharmonized cohort_long available)."
        )

    print("\nDone.")
    print(f"All outputs in: {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
