#!/usr/bin/env python3
"""Harmonize FastSurfer measurements with longitudinal ComBat.

Runs longCombat on the 12 panel ROIs (subcortical volumes and cortical
thicknesses), rebuilds bilateral means and eTIV-normalized values from
the harmonized halves, and recomputes per-subject deltas. Outputs:

    cohort_wide_harmonized.csv               31 cols, panel only
    cohort_cortical_regions_harmonized.csv   12 panel pairs
    subject_roi_deltas_harmonized.csv
    cohort_long_harmonized.csv               LME input
    cohort_long_unharmonized.csv             Risk-A reference
    interval_summary.csv
    longcombat_diagnostics.csv
    longcombat_summary.txt

Each tissue block is harmonized in one longCombat call (V=12 features)
so the empirical-Bayes shrinkage pools across all features. Hemispheres
are harmonized separately and bilateral means are rebuilt afterwards
as (left + right) / 2. The preserve formula matches the downstream
LME (Years_from_BL * Group_clinical + age_at_BL + PatientSex_clinical),
per Beer 2020 section 6. eTIV itself isn't harmonized; *_norm columns
are recomputed from harmonized raw / original eTIV.

Environment: needs the rpy2_clean venv plus R 4.5.2 and longCombat.
Source the activation wrapper first:

    source ~/rpy2_env.sh

Usage:
    python scripts/stage_c_analyze/06_harmonize.py \\
        --cohort-wide       outputs/stage_b_extract/cohort_wide.csv \\
        --cortical-regions  outputs/stage_b_extract/cohort_cortical_regions.csv \\
        --scanner-metadata  outputs/stage_a_metadata/scanner_metadata_with_covariates.csv \\
        --output-dir        outputs/stage_c_analyze/ \\
        --cohort-only

References:
    Beer JC et al. (2020). NeuroImage 220:117129.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from pipeline_lib.constants import (
    ALL_SUBCORTICAL_HARMONIZE_COLS,
    ALL_SUBCORTICAL_MEAN_COLS,
    BILATERAL_SUBCORTICAL_MEANS,
    COHORT_CORTICAL_HARMONIZED_OUTPUT_STRUCTS,
    COHORT_WIDE_HARMONIZED_OUTPUT_COLS,
    MIN_BATCH_SIZE_FAIL,
    MIN_BATCH_SIZE_WARN,
    ROI_PANEL_CORTICAL_STRUCTS,
)
from pipeline_lib.deltas import (
    build_cohort_long,
    compute_analysis_deltas,
    compute_interval_years,
)


# --- rpy2 bridge ---

_R_HANDLES: Optional[dict] = None


def get_r_handles() -> dict:
    """Import longCombat + rpy2 once and cache the handles."""
    global _R_HANDLES
    if _R_HANDLES is not None:
        return _R_HANDLES
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        from rpy2.robjects import default_converter, globalenv
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter

        longcombat = importr("longCombat")
        importr("nlme")
        importr("invgamma")

        converter = default_converter + pandas2ri.converter

        _R_HANDLES = {
            "ro": ro,
            "longcombat": longcombat,
            "pandas2ri": pandas2ri,
            "globalenv": globalenv,
            "converter": converter,
            "localconverter": localconverter,
        }
        return _R_HANDLES
    except ImportError as exc:
        print(f"ERROR: rpy2 / longCombat not available: {exc}",
              file=sys.stderr)
        print("Source ~/rpy2_env.sh first.", file=sys.stderr)
        sys.exit(1)


# --- Variance-explained-by-batch helper ---

def _variance_explained_by_batch(values: np.ndarray, batch: np.ndarray) -> float:
    """Fraction of variance in `values` explained by `batch` (one-way
    ANOVA R-squared). Harmonization should drive this close to zero."""
    values = np.asarray(values, dtype=float)
    batch = np.asarray(batch)
    finite = np.isfinite(values)
    if finite.sum() < 3:
        return np.nan
    v = values[finite]
    b = batch[finite]
    grand_mean = v.mean()
    ss_total = ((v - grand_mean) ** 2).sum()
    if ss_total <= 0:
        return 0.0
    ss_between = 0.0
    for lvl in np.unique(b):
        mask = b == lvl
        if mask.sum() == 0:
            continue
        ss_between += mask.sum() * (v[mask].mean() - grand_mean) ** 2
    return ss_between / ss_total


# --- Multi-feature longCombat call ---

def _r_safe_name(col: str) -> str:
    """Return an R-safe identifier for a Python column name.

    R's formula parser treats "-", ".", ":" and spaces specially.
    Replacing all four with underscores is enough; the inputs here
    always start with a letter, so leading-underscore isn't an issue.
    Cortical keys like 'lh:entorhinal' need the colon replacement
    specifically.
    """
    return (col.replace("-", "_")
               .replace(" ", "_")
               .replace(".", "_")
               .replace(":", "_"))


def _build_formula(
    usable_preserve: List[str],
    time_col: str,
    interaction_with: str,
) -> str:
    """Build the R formula string.

    time_col * interaction_with expands to main effects plus the
    interaction; the rest of preserve is added as additive terms.
    """
    other_terms = [
        c for c in usable_preserve
        if c not in (time_col, interaction_with)
    ]
    formula = f"{time_col} * {interaction_with}"
    if other_terms:
        formula += " + " + " + ".join(other_terms)
    return formula


def run_longcombat_block(
    long_df: pd.DataFrame,
    value_cols: List[str],
    batch_col: str,
    subject_col: str,
    time_col: str,
    preserve_cols: List[str],
    interaction_with: str,
    block_label: str = "",
) -> Tuple[Optional[pd.DataFrame], List[Dict[str, object]]]:
    """Run longCombat on a block of related features in a single call.

    The empirical-Bayes shrinkage pools across all features in
    value_cols. Must have at least 2 features (longCombat doesn't do
    single-feature mode), and they should be on a similar measurement
    scale.

    `interaction_with` is the preserve covariate expanded into
    time x interaction (R's `*` syntax). It must be in preserve_cols
    and multi-level after the single-level auto-drop step; otherwise
    the block returns None with an
    'error_interaction_covariate_single_level' status.

    Returns (harmonized_df, status_list):
      harmonized_df  same index as long_df, one column per value_col
                     with harmonized values (NaN where the row was
                     dropped). None if longCombat failed.
      status_list    one dict per feature with status, n_used,
                     var_explained_by_batch_before/after, etc.

    Doesn't raise; the caller decides what to do on failure.
    """
    H = get_r_handles()
    ro = H["ro"]
    globalenv = H["globalenv"]
    converter = H["converter"]
    localconverter = H["localconverter"]

    n_input = len(long_df)

    # Initialize per-feature status records.
    statuses: List[Dict[str, object]] = [
        {
            "feature":            f,
            "block":              block_label,
            "n_input":            n_input,
            "n_used":             0,
            "n_dropped_missing":  0,
            "status":             "pending",
        }
        for f in value_cols
    ]

    # Hard precondition: longCombat requires V >= 2.
    if len(value_cols) < 2:
        for s in statuses:
            s["status"] = (
                f"error_block_too_small:V={len(value_cols)} "
                f"(longCombat requires >=2 features)"
            )
        return None, statuses

    # Drop rows missing any value or covariate.
    required_cols = (
        [subject_col, batch_col, time_col]
        + list(preserve_cols)
        + list(value_cols)
    )
    # Deduplicate while preserving order.
    seen = set()
    required_cols = [c for c in required_cols
                     if not (c in seen or seen.add(c))]
    df = long_df.dropna(subset=required_cols).copy()
    n_used = len(df)
    n_dropped = n_input - n_used
    for s in statuses:
        s["n_used"] = n_used
        s["n_dropped_missing"] = n_dropped

    if n_used == 0:
        for s in statuses:
            s["status"] = "error_no_complete_rows"
        return None, statuses

    # Per-feature variance check. Zero variance in a pre-specified
    # panel feature is a data-extraction/merge bug, so abandon the
    # whole block rather than emit a partial output.
    zero_var_features: List[str] = []
    for f, s in zip(value_cols, statuses):
        if df[f].nunique() < 2:
            s["status"] = "error_zero_variance_feature"
            zero_var_features.append(f)

    if zero_var_features:
        for s in statuses:
            if s["status"] == "pending":
                s["status"] = (
                    f"error_block_abandoned_due_to_zero_var:"
                    f"{zero_var_features}"
                )
        return None, statuses

    usable_value_cols = list(value_cols)

    # Batch-size precondition.
    batch_counts = df[batch_col].value_counts()
    tiny_batches = batch_counts[batch_counts < MIN_BATCH_SIZE_FAIL]
    if len(tiny_batches) > 0:
        for f, s in zip(value_cols, statuses):
            if s["status"] == "pending":
                s["status"] = (
                    f"error_tiny_batch:{tiny_batches.to_dict()}"
                )
        return None, statuses

    # Preserve original row labels for alignment.
    # After dropna(), df.index is a subset of long_df.index. We use
    # label-based alignment to put harmonized values back, which works
    # for any index type.
    preserved_labels = df.index.copy()

    # Reset to a clean 0..M-1 RangeIndex for R. _orig_idx is kept as a
    # column so we can verify R didn't reorder the rows.
    df = df.reset_index(drop=True)
    df["_orig_idx"] = np.arange(len(df))

    # R-safe rename for value columns. Python keeps the original name;
    # R uses the safe form. The original is restored on read-back.
    rename_map: Dict[str, str] = {
        f: _r_safe_name(f) for f in usable_value_cols
        if _r_safe_name(f) != f
    }
    if rename_map:
        df = df.rename(columns=rename_map)
    safe_value_cols = [rename_map.get(f, f) for f in usable_value_cols]

    # Single-level preserve auto-drop.
    # A one-level preserve covariate has zero df for contrasts and
    # crashes R's `contrasts<-`. Time variable is always kept.
    usable_preserve: List[str] = []
    dropped_single_level: List[Tuple[str, object]] = []
    for c in preserve_cols:
        if c == time_col:
            usable_preserve.append(c)
            continue
        if c not in df.columns:
            continue
        n_levels = df[c].nunique(dropna=True)
        if n_levels < 2:
            only_val = df[c].iloc[0] if len(df) else None
            dropped_single_level.append((c, only_val))
        else:
            usable_preserve.append(c)

    # Build formula. If the interaction covariate didn't survive
    # auto-drop, fail the block; a partial-effect formula would
    # silently change what's modeled.
    if interaction_with not in usable_preserve:
        for s in statuses:
            if s["status"] == "pending":
                s["status"] = (
                    f"error_interaction_covariate_single_level:"
                    f"{interaction_with}"
                )
        return None, statuses

    formula_str = _build_formula(usable_preserve, time_col, interaction_with)

    # Convert string columns to factors in R.
    cols_r_will_see = ([subject_col, batch_col, time_col]
                       + usable_preserve)
    seen = set()
    cols_r_will_see = [c for c in cols_r_will_see
                       if not (c in seen or seen.add(c))]

    # Identifier/covariate names are passed through verbatim; value
    # columns get renamed by _r_safe_name above. Check that no id/
    # covariate has chars that would break R's formula parser.
    unsafe_id_cols = [
        c for c in cols_r_will_see
        if _r_safe_name(c) != c
    ]
    if unsafe_id_cols:
        for s in statuses:
            if s["status"] == "pending":
                s["status"] = (
                    f"error_unsafe_r_names:{unsafe_id_cols}. "
                    f"Identifier/covariate names containing '-', ' ', "
                    f"or '.' would crash R's formula parser. Rename "
                    f"these columns upstream of run_longcombat_block."
                )
        return None, statuses

    string_cols_to_factor = [
        c for c in cols_r_will_see
        if c in df.columns and df[c].dtype == object
    ]
    factor_conversion = "\n".join(
        f'    in_df${c} <- as.factor(in_df${c})'
        for c in string_cols_to_factor
    )

    # Build R script.
    features_r = "c(" + ", ".join(f'"{c}"' for c in safe_value_cols) + ")"
    r_script = f"""
    {factor_conversion}
    res <- longCombat::longCombat(
        idvar    = "{subject_col}",
        timevar  = "{time_col}",
        batchvar = "{batch_col}",
        features = {features_r},
        formula  = "{formula_str}",
        ranef    = "(1|{subject_col})",
        data     = in_df,
        verbose  = FALSE
    )
    harm <- res$data_combat
    """

    # Run R.
    try:
        with localconverter(converter):
            r_df = ro.conversion.get_conversion().py2rpy(df)
            globalenv["in_df"] = r_df
            ro.r(r_script)
            harm_df_r = ro.r("as.data.frame(harm)")
            harm_df = ro.conversion.get_conversion().rpy2py(harm_df_r)
    except Exception as exc:
        for s in statuses:
            if s["status"] == "pending":
                s["status"] = (
                    f"error_r:{type(exc).__name__}:{str(exc)[:200]}"
                )
        return None, statuses

    # Validate output columns.
    expected_combat_cols = [f"{c}.combat" for c in safe_value_cols]
    missing = [c for c in expected_combat_cols
               if c not in harm_df.columns]
    if missing:
        for s in statuses:
            if s["status"] == "pending":
                s["status"] = f"error_missing_output_cols:{missing}"
        return None, statuses

    # Validate row order and align harmonized values. longCombat
    # preserves input order; we verify the row count and _orig_idx
    # round-trip, then align using preserved_labels.
    if len(harm_df) != len(df):
        for s in statuses:
            if s["status"] == "pending":
                s["status"] = (
                    f"error_row_count_mismatch:"
                    f"input={len(df)}, output={len(harm_df)}"
                )
        return None, statuses

    if "_orig_idx" in harm_df.columns:
        round_tripped = np.asarray(harm_df["_orig_idx"].values, dtype=int)
        expected = np.arange(len(df))
        if not np.array_equal(round_tripped, expected):
            for s in statuses:
                if s["status"] == "pending":
                    s["status"] = (
                        "error_row_order_changed_in_R"
                    )
            return None, statuses

    harmonized_df = pd.DataFrame(
        np.nan, index=long_df.index, columns=value_cols, dtype=float
    )
    for orig_col, safe_col in zip(usable_value_cols, safe_value_cols):
        combat_col = f"{safe_col}.combat"
        harm_vals = np.asarray(harm_df[combat_col].values, dtype=float)
        # preserved_labels[i] is the original long_df index for row i.
        harmonized_df.loc[preserved_labels, orig_col] = harm_vals

    # Per-feature batch R² before/after.
    # Use a feature -> status lookup so reordering doesn't break the
    # R² assignment.
    status_by_feature: Dict[str, Dict[str, object]] = {
        s["feature"]: s for s in statuses
    }
    batch_arr = df[batch_col].values
    for orig_col, safe_col in zip(usable_value_cols, safe_value_cols):
        s = status_by_feature[orig_col]
        before = df[safe_col].values
        after = np.asarray(harm_df[f"{safe_col}.combat"].values, dtype=float)
        try:
            ve_before = _variance_explained_by_batch(before, batch_arr)
            ve_after = _variance_explained_by_batch(after, batch_arr)
        except Exception:
            ve_before, ve_after = np.nan, np.nan
        s["var_explained_by_batch_before"] = round(float(ve_before), 4)
        s["var_explained_by_batch_after"] = round(float(ve_after), 4)
        s["status"] = "ok"

    # Annotate statuses with block-level metadata.
    for s in statuses:
        s["formula_used"] = formula_str
        s["interaction_with"] = interaction_with
        s["features_in_block"] = len(usable_value_cols)
        s["block_n_used"] = n_used
        if dropped_single_level:
            s["dropped_single_level"] = dropped_single_level

    return harmonized_df, statuses


# --- Input preparation ---

def _apply_actual_years_from_bl(
    df: pd.DataFrame,
    interval_summary: pd.DataFrame,
) -> pd.DataFrame:
    """Override the nominal Years_from_BL with actual elapsed time.

    In Stage B writes Years_from_BL as the nominal map (BL=0, 3Y=3, 5Y=5),
    but the harmonization model should match the analysis model (Beer
    2020 section 6), so we swap nominal with actual before longCombat
    fits its preserve formula. BL rows get 0.0; follow-up rows get
    interval_years from compute_interval_years.

    Raises RuntimeError if any row has no matching entry in
    interval_summary.
    """
    required_in_iv = {"Subject", "Delta_Window", "interval_years"}
    missing_iv = required_in_iv - set(interval_summary.columns)
    if missing_iv:
        raise ValueError(
            f"_apply_actual_years_from_bl: interval_summary is "
            f"missing required columns: {sorted(missing_iv)}"
        )
    required_in_df = {"Subject", "Timepoint", "Years_from_BL"}
    missing_df = required_in_df - set(df.columns)
    if missing_df:
        raise ValueError(
            f"_apply_actual_years_from_bl: input df is missing "
            f"required columns: {sorted(missing_df)}"
        )

    # (Subject, Timepoint) -> actual_years lookup.
    years_lookup: Dict[Tuple[str, str], float] = {}
    for subj in df["Subject"].dropna().unique():
        years_lookup[(subj, "BL")] = 0.0
    for _, row in interval_summary.iterrows():
        try:
            tp2 = row["Delta_Window"].split("→")[1]
        except (AttributeError, IndexError):
            # Malformed Delta_Window; downstream miss will hard-fail.
            continue
        years_lookup[(row["Subject"], tp2)] = float(row["interval_years"])

    # Apply override.
    out = df.copy()
    keys = list(zip(out["Subject"], out["Timepoint"]))
    actual = pd.Series(
        [years_lookup.get(k, np.nan) for k in keys],
        index=out.index,
        dtype=float,
    )
    n_missing = int(actual.isna().sum())
    if n_missing:
        missing_examples = [
            k for k, v in zip(keys, actual)
            if pd.isna(v)
        ][:5]
        raise RuntimeError(
            f"_apply_actual_years_from_bl: {n_missing} row(s) have "
            f"no matching entry in interval_summary "
            f"({len(keys) - n_missing} of {len(keys)} matched). "
            f"First 5 missing (Subject, Timepoint) keys: "
            f"{missing_examples}."
        )

    out["Years_from_BL"] = actual
    return out


def merge_scanner_metadata(
    cohort_wide: pd.DataFrame,
    scanner_meta: pd.DataFrame,
    preserve_cols: List[str],
    cohort_only: bool,
) -> pd.DataFrame:
    """Merge BatchID and preserve covariates onto cohort_wide.

    Returns the same rows as cohort_wide (or filtered to cohort-only)
    with the added columns. Raises if any row is missing BatchID or
    a preserve covariate after merge.
    """
    needed = ["Subject", "Timepoint", "BatchID"] + preserve_cols
    missing_in_meta = [c for c in needed if c not in scanner_meta.columns]
    if missing_in_meta:
        raise ValueError(
            f"scanner_metadata is missing required columns: {missing_in_meta}"
        )

    meta_slim = scanner_meta[needed].drop_duplicates(["Subject", "Timepoint"])
    merged = cohort_wide.merge(
        meta_slim, on=["Subject", "Timepoint"], how="left"
    )

    if cohort_only and "in_cohort" in scanner_meta.columns:
        cohort_subs = set(
            scanner_meta.loc[scanner_meta["in_cohort"], "Subject"]
        )
        merged = merged[merged["Subject"].isin(cohort_subs)].copy()

    missing_batch = merged["BatchID"].isna().sum()
    if missing_batch:
        bad = (merged.loc[merged["BatchID"].isna(),
                          ["Subject", "Timepoint"]]
                     .drop_duplicates().values.tolist())
        raise ValueError(
            f"{missing_batch} rows have no BatchID after merge. "
            f"Examples: {bad[:5]}."
        )
    for col in preserve_cols:
        n_miss = merged[col].isna().sum()
        if n_miss:
            bad = (merged.loc[merged[col].isna(),
                              ["Subject", "Timepoint"]]
                         .drop_duplicates().values.tolist())
            raise ValueError(
                f"{n_miss} rows are missing preserve covariate "
                f"'{col}'. Examples: {bad[:5]}."
            )
    return merged


def merge_scanner_metadata_cortical(
    cohort_cortical: pd.DataFrame,
    scanner_meta: pd.DataFrame,
    preserve_cols: List[str],
    cohort_only: bool,
) -> pd.DataFrame:
    """Same as merge_scanner_metadata but for the cortical long-format table."""
    needed = ["Subject", "Timepoint", "BatchID"] + preserve_cols
    meta_slim = scanner_meta[needed].drop_duplicates(["Subject", "Timepoint"])
    merged = cohort_cortical.merge(
        meta_slim, on=["Subject", "Timepoint"], how="left"
    )
    if cohort_only and "in_cohort" in scanner_meta.columns:
        cohort_subs = set(
            scanner_meta.loc[scanner_meta["in_cohort"], "Subject"]
        )
        merged = merged[merged["Subject"].isin(cohort_subs)].copy()

    missing_batch = merged["BatchID"].isna().sum()
    if missing_batch:
        raise ValueError(
            f"{missing_batch} cortical rows have no BatchID after merge."
        )
    return merged


# --- Harmonization drivers ---

def _print_block_preconditions(
    df: pd.DataFrame,
    batch_col: str,
    preserve_cols: List[str],
    time_col: str,
) -> None:
    """Log batch sizes and single-level preserve covariates before the
    block call. Informational; the real checks are in
    run_longcombat_block."""
    batch_sizes = df[batch_col].value_counts().sort_index()
    print("\n  Batch sizes in harmonization input:")
    for bid, n in batch_sizes.items():
        flag = ""
        if n < MIN_BATCH_SIZE_FAIL:
            flag = "  <- TOO SMALL (call will fail)"
        elif n < MIN_BATCH_SIZE_WARN:
            flag = "  <- below n=10 guideline (unreliable)"
        print(f"    {bid:20s}  n={n}{flag}")

    single_level = []
    for c in preserve_cols:
        if c == time_col:
            continue
        if c in df.columns:
            n_levels = df[c].nunique(dropna=True)
            if n_levels < 2:
                only_val = df[c].dropna().iloc[0] if df[c].notna().any() else None
                single_level.append((c, only_val))
    if single_level:
        print("\n  These preserve covariates have only one level and "
              "will be auto-dropped from the longCombat formula:")
        for col, val in single_level:
            print(f"    {col} = {val!r} (constant; cannot estimate contrasts)")


def _write_failure_diagnostics(
    reports: List[Dict[str, object]],
    output_dir: Path,
    section: str,
) -> Path:
    """Write a partial diagnostics CSV before raising on failure.

    The filename is suffixed _FAILED_<section> so a successful re-run
    doesn't overwrite the failure record.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"longcombat_diagnostics_FAILED_{section}.csv"
    pd.DataFrame(reports).to_csv(path, index=False)
    print(f"  Wrote failure diagnostics: {path}", file=sys.stderr)
    return path


def harmonize_cohort_wide(
    wide_merged: pd.DataFrame,
    preserve_cols: List[str],
    interaction_with: str,
    output_dir: Path,
) -> Tuple[pd.DataFrame, List[Dict[str, object]]]:
    """Harmonize the 12 panel hemispheric subcortical volumes.

    Single multi-feature longCombat call; bilateral means and *_norm
    columns are rebuilt afterwards. Other columns pass through
    unchanged in memory and get dropped by the slim output writer.

    Raises RuntimeError on longCombat failure; a failure-diagnostics
    CSV is written to output_dir first.
    """
    out = wide_merged.copy()

    missing_cols = [c for c in ALL_SUBCORTICAL_HARMONIZE_COLS
                    if c not in out.columns]
    if missing_cols:
        raise ValueError(
            f"cohort_wide is missing expected harmonization columns: "
            f"{missing_cols}"
        )

    _print_block_preconditions(
        out, batch_col="BatchID",
        preserve_cols=preserve_cols, time_col="Years_from_BL",
    )

    print(f"\n  Harmonizing {len(ALL_SUBCORTICAL_HARMONIZE_COLS)} "
          f"subcortical volume features (V={len(ALL_SUBCORTICAL_HARMONIZE_COLS)}, "
          f"single multi-feature call)...")

    harmonized_df, reports = run_longcombat_block(
        long_df=out,
        value_cols=list(ALL_SUBCORTICAL_HARMONIZE_COLS),
        batch_col="BatchID",
        subject_col="Subject",
        time_col="Years_from_BL",
        preserve_cols=list(preserve_cols),
        interaction_with=interaction_with,
        block_label="subcortical",
    )

    if harmonized_df is None:
        first_status = reports[0].get("status", "unknown") if reports else "unknown"
        for r in reports:
            print(f"    FAIL {r['feature']:30s}  {r.get('status', 'unknown')}",
                  file=sys.stderr)
        fail_path = _write_failure_diagnostics(
            reports, output_dir, "subcortical"
        )
        raise RuntimeError(
            f"Subcortical longCombat block failed: {first_status}. "
            f"Per-feature details: {fail_path}"
        )

    # Apply harmonized values to the output DataFrame and log per-feature R².
    for feat in ALL_SUBCORTICAL_HARMONIZE_COLS:
        out[feat] = harmonized_df[feat]
        s = next((r for r in reports if r["feature"] == feat), None)
        if s and s.get("status") == "ok":
            vb = s.get("var_explained_by_batch_before", np.nan)
            va = s.get("var_explained_by_batch_after", np.nan)
            print(f"    ok  {feat:30s}  batch R² {vb:.3f} -> {va:.3f}")
        else:
            print(f"    SKIP {feat:30s}  "
                  f"{s.get('status') if s else 'no report'}")

    if reports:
        formula = reports[0].get("formula_used", "?")
        print(f"\n  Formula used: {formula!r}")

    # Rebuild bilateral mean columns from harmonized halves.
    for mean_col, (lc, rc) in BILATERAL_SUBCORTICAL_MEANS.items():
        out[mean_col] = (out[lc] + out[rc]) / 2.0

    # Recompute *_norm columns from harmonized raws / original eTIV.
    # eTIV itself isn't harmonized, which avoids partial cancellation
    # between the harmonization shift and the normalization step.
    if "eTIV" not in out.columns:
        raise ValueError(
            "cohort_wide has no eTIV column; cannot recompute _norm."
        )
    etiv = out["eTIV"]
    cols_to_renorm = (
        list(ALL_SUBCORTICAL_HARMONIZE_COLS)
        + list(ALL_SUBCORTICAL_MEAN_COLS)
    )
    n_renormed = 0
    for col in cols_to_renorm:
        norm_col = f"{col}_norm"
        if norm_col in out.columns:
            out[norm_col] = np.where(
                (etiv > 0) & np.isfinite(etiv),
                out[col] / etiv,
                np.nan,
            )
            n_renormed += 1
    print(f"\n  Recomputed {n_renormed} *_norm columns from harmonized "
          f"raws / original eTIV")

    return out, reports


def harmonize_cortical(
    cort_merged: pd.DataFrame,
    preserve_cols: List[str],
    interaction_with: str,
    output_dir: Path,
) -> Tuple[pd.DataFrame, List[Dict[str, object]]]:
    """Harmonize ThickAvg for the 12 panel (Hemisphere, StructName) pairs.

    Single multi-feature longCombat call. Returns only the in-scope
    panel rows; out-of-scope regions (postcentral, fusiform, etc.)
    are dropped so the *_harmonized.csv file only holds values that
    were actually harmonized.

    Raises RuntimeError on longCombat failure; a failure-diagnostics
    CSV is written to output_dir first.
    """
    in_scope_mask = cort_merged["StructName"].isin(ROI_PANEL_CORTICAL_STRUCTS)
    in_scope = cort_merged[in_scope_mask].copy()
    n_in = len(in_scope)
    n_out = len(cort_merged) - n_in
    print(f"\n  Cortical rows total: {len(cort_merged)}  "
          f"(in panel scope: {n_in}, dropped from harmonized output: "
          f"{n_out})")

    # Build the (Hemisphere, StructName) feature key column.
    in_scope["_feature_key"] = (
        in_scope["Hemisphere"].astype(str)
        + ":"
        + in_scope["StructName"].astype(str)
    )

    # Pivot to wide: one row per (Subject, Timepoint, covariates),
    # one column per feature key.
    id_cols = (
        ["Subject", "Timepoint", "Years_from_BL", "BatchID"]
        + [c for c in preserve_cols
           if c != "Years_from_BL" and c in in_scope.columns]
    )
    seen = set()
    id_cols = [c for c in id_cols if not (c in seen or seen.add(c))]

    missing_id = [c for c in id_cols if c not in in_scope.columns]
    if missing_id:
        raise ValueError(
            f"cortical input missing id columns: {missing_id}"
        )

    # aggfunc='first'; the key (Subject, Timepoint, Hemisphere,
    # StructName) should be unique anyway.
    wide = in_scope.pivot_table(
        index=id_cols,
        columns="_feature_key",
        values="ThickAvg",
        aggfunc="first",
    ).reset_index()

    # Canonical feature ordering: lh and rh for each panel ROI.
    expected_features = [
        f"{h}:{s}" for h in ("lh", "rh")
        for s in ROI_PANEL_CORTICAL_STRUCTS
    ]
    available_features = set(wide.columns) - set(id_cols)
    feature_cols = [f for f in expected_features if f in available_features]
    missing_features = [f for f in expected_features if f not in available_features]
    if missing_features:
        print(f"  WARNING: expected feature(s) absent from input: "
              f"{missing_features}", file=sys.stderr)

    print(f"  Harmonizing {len(feature_cols)} (hemisphere x region) "
          f"features (V={len(feature_cols)}, single multi-feature call)...")

    if len(feature_cols) < 2:
        raise RuntimeError(
            f"Cortical block has only {len(feature_cols)} features "
            f"after panel filter; longCombat requires >=2."
        )

    harmonized_wide, reports = run_longcombat_block(
        long_df=wide,
        value_cols=feature_cols,
        batch_col="BatchID",
        subject_col="Subject",
        time_col="Years_from_BL",
        preserve_cols=list(preserve_cols),
        interaction_with=interaction_with,
        block_label="cortical",
    )

    if harmonized_wide is None:
        first_status = reports[0].get("status", "unknown") if reports else "unknown"
        for r in reports:
            print(f"    FAIL {r['feature']:30s}  {r.get('status', 'unknown')}",
                  file=sys.stderr)
        fail_path = _write_failure_diagnostics(
            reports, output_dir, "cortical"
        )
        raise RuntimeError(
            f"Cortical longCombat block failed: {first_status}. "
            f"Per-feature details: {fail_path}"
        )

    # Log per-feature R² values.
    for feat in feature_cols:
        s = next((r for r in reports if r["feature"] == feat), None)
        if s and s.get("status") == "ok":
            vb = s.get("var_explained_by_batch_before", np.nan)
            va = s.get("var_explained_by_batch_after", np.nan)
            print(f"    ok  {feat:30s}  batch R² {vb:.3f} -> {va:.3f}")
        else:
            print(f"    SKIP {feat:30s}  "
                  f"{s.get('status') if s else 'no report'}")

    if reports:
        formula = reports[0].get("formula_used", "?")
        print(f"\n  Formula used: {formula!r}")

    # Project harmonized values back into long format. In-scope
    # long-format with ThickAvg replaced; out-of-scope rows are not
    # carried through.

    # Melt feature_cols back into rows.
    harm_wide_with_id = wide[id_cols].copy()
    for feat in feature_cols:
        harm_wide_with_id[feat] = harmonized_wide[feat].values

    harm_long = harm_wide_with_id.melt(
        id_vars=id_cols,
        value_vars=feature_cols,
        var_name="_feature_key",
        value_name="ThickAvg_harmonized",
    )
    # Split _feature_key back into Hemisphere and StructName.
    split_keys = harm_long["_feature_key"].str.split(":", n=1, expand=True)
    harm_long["Hemisphere"] = split_keys[0]
    harm_long["StructName"] = split_keys[1]

    # Merge harmonized values onto the in-scope long-format frame.
    # Check both sides for duplicate keys first.
    in_scope_dups = in_scope.duplicated(
        subset=["Subject", "Timepoint", "Hemisphere", "StructName"]
    ).sum()
    harm_dups = harm_long.duplicated(
        subset=["Subject", "Timepoint", "Hemisphere", "StructName"]
    ).sum()
    if in_scope_dups or harm_dups:
        raise RuntimeError(
            f"Duplicate (Subject, Timepoint, Hemisphere, StructName) "
            f"keys detected: in_scope={in_scope_dups}, "
            f"harmonized={harm_dups}."
        )

    in_scope = in_scope.merge(
        harm_long[["Subject", "Timepoint", "Hemisphere", "StructName",
                   "ThickAvg_harmonized"]],
        on=["Subject", "Timepoint", "Hemisphere", "StructName"],
        how="left",
    )

    # Fail loudly on partial-harmonized output: a row with non-NaN
    # original ThickAvg but no harmonized replacement would be
    # misleading in a file named *_harmonized.csv.
    missing_harm = (
        in_scope["ThickAvg"].notna()
        & in_scope["ThickAvg_harmonized"].isna()
    )
    if missing_harm.any():
        bad_examples = (
            in_scope.loc[
                missing_harm,
                ["Subject", "Timepoint", "Hemisphere", "StructName"]
            ]
            .head(10)
            .to_dict("records")
        )
        raise RuntimeError(
            f"{int(missing_harm.sum())} cortical panel rows had a "
            f"non-NaN original ThickAvg but no harmonized replacement. "
            f"First 10 examples: {bad_examples}"
        )

    in_scope["ThickAvg"] = in_scope["ThickAvg_harmonized"]
    in_scope = in_scope.drop(columns=["ThickAvg_harmonized", "_feature_key"])

    return in_scope, reports


# --- Delta computation and covariate attach ---

CLINICAL_COVARIATE_COLS = ["age_at_BL", "PatientSex_clinical", "Group_clinical"]


def attach_clinical_covariates(
    deltas: pd.DataFrame,
    scanner_meta: pd.DataFrame,
) -> pd.DataFrame:
    """Merge per-subject clinical covariates onto the deltas table.

    Covariates are subject-invariant, so scanner_meta is deduplicated
    to one row per subject before merging.
    """
    missing_in_meta = [c for c in CLINICAL_COVARIATE_COLS
                       if c not in scanner_meta.columns]
    if missing_in_meta:
        raise ValueError(
            f"scanner_metadata is missing required covariate columns: "
            f"{missing_in_meta}."
        )

    meta_per_subject = (
        scanner_meta[["Subject"] + CLINICAL_COVARIATE_COLS]
        .drop_duplicates(subset="Subject", keep="first")
    )
    return deltas.merge(meta_per_subject, on="Subject", how="left")


# --- Output writing ---

def write_harmonized_outputs(
    wide_h: pd.DataFrame,
    cort_h: pd.DataFrame,
    output_dir: Path,
    wide_reports: List[dict],
    cort_reports: List[dict],
) -> Tuple[Path, Path]:
    """Write the slim harmonized outputs plus diagnostics.

    cohort_wide_harmonized.csv (31 cols),
    cohort_cortical_regions_harmonized.csv (panel pairs only),
    longcombat_diagnostics.csv, and longcombat_summary.txt.

    Returns (wide_path, cort_path); main re-reads these to feed
    compute_analysis_deltas.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Slim wide schema.
    missing_wide = [c for c in COHORT_WIDE_HARMONIZED_OUTPUT_COLS
                    if c not in wide_h.columns]
    if missing_wide:
        raise RuntimeError(
            f"Harmonized wide DataFrame is missing "
            f"{len(missing_wide)} expected column(s): {missing_wide}."
        )
    wide_out_cols = list(COHORT_WIDE_HARMONIZED_OUTPUT_COLS)

    wide_path = output_dir / "cohort_wide_harmonized.csv"
    # 10dp matches Stage B's float_format so cohort_wide.csv and
    # cohort_wide_harmonized.csv can be diffed cleanly.
    wide_h[wide_out_cols].to_csv(wide_path, index=False,
                                 float_format="%.10f")
    print(f"\n  Wrote: {wide_path}  "
          f"({len(wide_h)} rows x {len(wide_out_cols)} columns)")

    # Slim cortical schema. cort_h is already filtered by
    # harmonize_cortical; this re-filter and the schema asserts are
    # defensive.
    cort_in_panel = cort_h[
        cort_h["StructName"].isin(COHORT_CORTICAL_HARMONIZED_OUTPUT_STRUCTS)
    ].copy()

    expected_structs = set(COHORT_CORTICAL_HARMONIZED_OUTPUT_STRUCTS)
    actual_structs = set(cort_in_panel["StructName"].dropna().unique())
    if actual_structs != expected_structs:
        missing_s = expected_structs - actual_structs
        extra_s = actual_structs - expected_structs
        raise RuntimeError(
            f"Cortical StructName set mismatch in harmonized output. "
            f"Missing: {sorted(missing_s)}. Extra: {sorted(extra_s)}."
        )

    expected_pairs = 2 * len(expected_structs)
    actual_pairs = (
        cort_in_panel[["Hemisphere", "StructName"]]
        .drop_duplicates()
        .shape[0]
    )
    if actual_pairs != expected_pairs:
        raise RuntimeError(
            f"Expected {expected_pairs} unique "
            f"(Hemisphere, StructName) pairs, got {actual_pairs}."
        )

    cort_path = output_dir / "cohort_cortical_regions_harmonized.csv"
    cort_in_panel.to_csv(cort_path, index=False, float_format="%.10f")
    print(f"  Wrote: {cort_path}  "
          f"({len(cort_in_panel)} rows x {len(cort_in_panel.columns)} columns)")

    # longcombat_diagnostics.csv.
    report_path = output_dir / "longcombat_diagnostics.csv"
    rep_rows = []
    for section, reports in [("cohort_wide", wide_reports),
                             ("cortical", cort_reports)]:
        for r in reports:
            rep_rows.append({"section": section, **r})
    pd.DataFrame(rep_rows).to_csv(report_path, index=False)
    print(f"  Wrote: {report_path}")

    # longcombat_summary.txt.
    summary_path = output_dir / "longcombat_summary.txt"
    with summary_path.open("w") as f:
        for label, reports in [("subcortical", wide_reports),
                               ("cortical", cort_reports)]:
            ok = [r for r in reports if r.get("status") == "ok"]
            if not ok:
                f.write(f"{label}: 0 features harmonized\n")
                continue
            avg_b = np.mean([r["var_explained_by_batch_before"] for r in ok])
            avg_a = np.mean([r["var_explained_by_batch_after"] for r in ok])
            f.write(f"{label}: {len(ok)} features harmonized, "
                    f"batch R^2 {avg_b:.4f} -> {avg_a:.4f}, "
                    f"formula {ok[0].get('formula_used', '?')!r}\n")
    print(f"  Wrote: {summary_path}")
    return wide_path, cort_path


# --- Main ---

def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--cohort-wide", type=Path, required=True)
    p.add_argument("--cortical-regions", type=Path, required=True)
    p.add_argument("--scanner-metadata", type=Path, required=True,
                   help="scanner_metadata_with_covariates.csv from "
                        "04_merge_metadata.py.")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument(
        "--preserve", nargs="+",
        default=["Years_from_BL", "age_at_BL",
                 "PatientSex_clinical", "Group_clinical"],
        help="Biological covariates whose variance longCombat should "
             "preserve. Years_from_BL is the time axis and must be "
             "present.",
    )
    p.add_argument(
        "--cohort-only", action="store_true",
        help="Only harmonize rows where scanner_metadata has "
             "in_cohort=True. Recommended.",
    )
    args = p.parse_args()

    # Validate preserve list.
    if "Years_from_BL" not in args.preserve:
        print("ERROR: --preserve must include Years_from_BL",
              file=sys.stderr)
        return 1
    if "Group_clinical" not in args.preserve:
        print("ERROR: --preserve must include Group_clinical "
              "(needed for the Years_from_BL x Group_clinical "
              "interaction).", file=sys.stderr)
        return 1

    interaction_with = "Group_clinical"

    # Load inputs.
    print(f"Loading {args.cohort_wide}...")
    wide = pd.read_csv(args.cohort_wide)
    print(f"  {len(wide)} rows, {wide['Subject'].nunique()} subjects")

    print(f"Loading {args.cortical_regions}...")
    cort = pd.read_csv(args.cortical_regions)
    print(f"  {len(cort)} rows, {cort['Subject'].nunique()} subjects, "
          f"{cort['StructName'].nunique()} regions x "
          f"{cort['Hemisphere'].nunique()} hemispheres")

    print(f"Loading {args.scanner_metadata}...")
    scan = pd.read_csv(args.scanner_metadata)
    print(f"  {len(scan)} rows, {scan['Subject'].nunique()} subjects")

    for df, name in [(wide, "cohort_wide"), (cort, "cortical")]:
        if "Years_from_BL" not in df.columns:
            print(f"ERROR: {name} is missing Years_from_BL",
                  file=sys.stderr)
            return 1

    preserve_from_meta = [c for c in args.preserve if c != "Years_from_BL"]

    # Merge batch + clinical.
    print("\nMerging scanner metadata + clinical covariates into cohort_wide...")
    wide_merged = merge_scanner_metadata(
        wide, scan, preserve_cols=preserve_from_meta,
        cohort_only=args.cohort_only,
    )
    print(f"  {len(wide_merged)} rows after merge "
          f"(cohort_only={args.cohort_only})")

    print("\nMerging scanner metadata + clinical covariates into cortical...")
    cort_merged = merge_scanner_metadata_cortical(
        cort, scan, preserve_cols=preserve_from_meta,
        cohort_only=args.cohort_only,
    )
    print(f"  {len(cort_merged)} rows after merge")

    # Compute actual elapsed time and override Years_from_BL. Stage B
    # writes the nominal map; we swap to actual elapsed time so the
    # harmonization model matches the analysis model.
    print("\nComputing per-subject interval_years from DICOM dates "
          "(pre-harmonization)...")
    cohort_subjects = set(wide_merged["Subject"].unique())
    scan_cohort = scan[scan["Subject"].isin(cohort_subjects)].copy()
    n_excluded = scan["Subject"].nunique() - len(cohort_subjects)
    if n_excluded:
        print(f"  Filtering scanner_meta to {len(cohort_subjects)} cohort "
              f"subjects ({n_excluded} non-cohort excluded from "
              f"interval_summary).")
    intervals = compute_interval_years(scan_cohort)
    interval_summary_path = args.output_dir / "interval_summary.csv"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    intervals.to_csv(interval_summary_path, index=False,
                     float_format="%.10f")
    print(f"  Wrote: {interval_summary_path} ({len(intervals)} rows)")

    print("\n  Overriding Years_from_BL on wide_merged and cort_merged "
          "with actual elapsed time...")
    wide_merged = _apply_actual_years_from_bl(wide_merged, intervals)
    cort_merged = _apply_actual_years_from_bl(cort_merged, intervals)
    # Range of actual times by timepoint.
    for label, df in [("wide_merged", wide_merged),
                      ("cort_merged", cort_merged)]:
        per_tp = df.groupby("Timepoint")["Years_from_BL"].agg(
            ["min", "mean", "max"]
        ).round(3)
        print(f"  {label} Years_from_BL by timepoint:")
        for tp, row in per_tp.iterrows():
            print(f"    {tp:4s}  min={row['min']:.3f}  "
                  f"mean={row['mean']:.3f}  max={row['max']:.3f}")

    # Harmonize.
    print("\nHarmonizing cohort_wide (subcortical volumes)...")
    wide_h, wide_reports = harmonize_cohort_wide(
        wide_merged,
        preserve_cols=args.preserve,
        interaction_with=interaction_with,
        output_dir=args.output_dir,
    )

    print("\nHarmonizing cohort_cortical_regions (ThickAvg)...")
    cort_h, cort_reports = harmonize_cortical(
        cort_merged,
        preserve_cols=args.preserve,
        interaction_with=interaction_with,
        output_dir=args.output_dir,
    )

    # Write harmonized cohort tables.
    print("\nWriting harmonized cohort tables (slim schemas)...")
    wide_path, cort_path = write_harmonized_outputs(
        wide_h=wide_h, cort_h=cort_h,
        output_dir=args.output_dir,
        wide_reports=wide_reports,
        cort_reports=cort_reports,
    )

    # Compute harmonized deltas. Read the just-written CSVs back so
    # compute_analysis_deltas sees the same values a downstream
    # consumer would.
    print("\nComputing per-subject analysis deltas (harmonized)...")
    harm_wide = pd.read_csv(wide_path)
    harm_cort = pd.read_csv(cort_path)

    deltas = compute_analysis_deltas(harm_wide, harm_cort)

    if deltas.empty:
        raise RuntimeError(
            "compute_analysis_deltas returned 0 rows after harmonization."
        )

    print(f"  {len(deltas)} delta rows "
          f"({deltas['Subject'].nunique()} subjects, "
          f"{deltas['ROI'].nunique()} ROIs, "
          f"{deltas['Delta_Window'].nunique()} delta windows)")
    print(f"  Tissue breakdown:  {dict(deltas['Tissue'].value_counts())}")
    print(f"  Cluster breakdown: {dict(deltas['Cluster'].value_counts())}")

    print("  Attaching clinical covariates to deltas...")
    deltas = attach_clinical_covariates(deltas, scan)
    n_with_age = deltas["age_at_BL"].notna().sum()
    print(f"  {n_with_age} / {len(deltas)} rows have clinical data")

    # Merge interval_years onto deltas.
    print("\nMerging interval_years onto deltas...")
    n_rows_before = len(deltas)
    deltas = deltas.merge(
        intervals[["Subject", "Delta_Window",
                   "interval_years", "interval_source"]],
        on=["Subject", "Delta_Window"], how="left",
    )
    assert len(deltas) == n_rows_before, "interval merge changed row count"
    n_with_interval = deltas["interval_years"].notna().sum()
    print(f"  {n_with_interval} / {len(deltas)} delta rows have "
          f"interval_years populated")
    if n_with_interval < len(deltas):
        missing = deltas[deltas["interval_years"].isna()]
        raise RuntimeError(
            f"{len(missing)} delta rows have no interval_years. "
            f"First 5 affected subjects: "
            f"{sorted(missing['Subject'].unique())[:5]}."
        )

    deltas_path = args.output_dir / "subject_roi_deltas_harmonized.csv"
    # Match Stage B's float_format so harmonized and unharmonized
    # deltas can be diffed cleanly.
    deltas.to_csv(deltas_path, index=False, float_format="%.10f")
    print(f"  Wrote: {deltas_path} ({len(deltas)} rows, "
          f"{len(deltas.columns)} columns)")

    # Build long-format dataset for LME analysis.
    # cohort_long_harmonized.csv is the main input to 07_analysis.py.
    # One row per (Subject, Timepoint, ROI) with covariates.
    print("\nBuilding long-format dataset for LME analysis...")
    cohort_long = build_cohort_long(
        cohort_wide=harm_wide,
        cohort_cortical=harm_cort,
        scanner_meta=scan_cohort,
        interval_summary=intervals,
    )
    if cohort_long.empty:
        print("  WARNING: cohort_long is empty; LME analysis cannot run.",
              file=sys.stderr)
    else:
        n_subjects = cohort_long["Subject"].nunique()
        n_rois = cohort_long["ROI"].nunique()
        n_timepoints = cohort_long["Timepoint"].nunique()
        print(f"  {len(cohort_long)} rows "
              f"({n_subjects} subjects x {n_timepoints} timepoints x "
              f"{n_rois} ROIs)")
        print(f"  Tissue breakdown:  "
              f"{dict(cohort_long['Tissue'].value_counts())}")

    cohort_long_path = args.output_dir / "cohort_long_harmonized.csv"
    cohort_long.to_csv(cohort_long_path, index=False, float_format="%.10f")
    print(f"  Wrote: {cohort_long_path} ({len(cohort_long)} rows, "
          f"{len(cohort_long.columns)} columns)")

    # Build cohort_long_unharmonized.csv for the Risk-A sanity check
    # in 07, which verifies the Group x Time interaction sign and
    # magnitude survive harmonization.
    print("\nBuilding cohort_long_unharmonized.csv (Risk-A artifact)...")
    cohort_long_unharm = build_cohort_long(
        cohort_wide=wide_merged,
        cohort_cortical=cort_merged,
        scanner_meta=scan_cohort,
        interval_summary=intervals,
    )
    if cohort_long_unharm.empty:
        print("  WARNING: cohort_long_unharmonized is empty.",
              file=sys.stderr)
    else:
        n_subjects = cohort_long_unharm["Subject"].nunique()
        n_rois = cohort_long_unharm["ROI"].nunique()
        n_timepoints = cohort_long_unharm["Timepoint"].nunique()
        print(f"  {len(cohort_long_unharm)} rows "
              f"({n_subjects} subjects x {n_timepoints} timepoints x "
              f"{n_rois} ROIs)")

    cohort_long_unharm_path = (
        args.output_dir / "cohort_long_unharmonized.csv"
    )
    cohort_long_unharm.to_csv(
        cohort_long_unharm_path, index=False, float_format="%.10f"
    )
    print(f"  Wrote: {cohort_long_unharm_path} "
          f"({len(cohort_long_unharm)} rows, "
          f"{len(cohort_long_unharm.columns)} columns)")

    # Final batch R² summary.
    ok_wide = [r for r in wide_reports if r.get("status") == "ok"]
    ok_cort = [r for r in cort_reports if r.get("status") == "ok"]
    all_ok = ok_wide + ok_cort
    if all_ok:
        avg_b = np.mean([r["var_explained_by_batch_before"] for r in all_ok])
        avg_a = np.mean([r["var_explained_by_batch_after"]  for r in all_ok])
        print(f"\nMean batch R² across {len(all_ok)} harmonized features:")
        print(f"  Before: {avg_b:.4f}")
        print(f"  After:  {avg_a:.4f}  "
              f"(should be << before if harmonization worked)")
        print(f"  Formula: {all_ok[0].get('formula_used')!r}")

    print("\nNext: run 07_analysis.py on "
          f"{cohort_long_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
