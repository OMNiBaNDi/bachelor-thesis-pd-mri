#!/usr/bin/env python3
"""Aggregate FSQC metrics across cohorts and apply threshold-based flags.

Usage example:
  python3 scripts/fsqc_analyze.py \
    --patients-metrics /nfs/.../outputs/fsqc/patients_20260428/metrics/fsqc_metrics.csv \
    --controls-metrics /nfs/.../outputs/fsqc/controls_20260428/metrics/fsqc_metrics.csv \
    --config projects/bachelor-thesis/config/fsqc_thresholds.yaml \
    --output-parquet /nfs/.../outputs/fsqc/data/fsqc_20260428.parquet \
    --decision-log /nfs/.../outputs/fsqc/qc/fsqc_decision_log.csv
"""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import yaml

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = WORKSPACE_ROOT / "config/fsqc_thresholds.yaml"


def _load_metrics(path_or_dir: Path, cohort: str) -> pd.DataFrame:
    """Return DataFrame of metrics for a cohort (aggregated file or per-subject tree)."""
    if path_or_dir.is_file():
        aggregated = path_or_dir
    else:
        for candidate in ("fsqc_metrics.csv", "fsqc_results.csv"):
            aggregated_candidate = path_or_dir / candidate
            if aggregated_candidate.is_file():
                aggregated = aggregated_candidate
                break
        else:
            aggregated = None

    if aggregated and aggregated.is_file():
        df = pd.read_csv(aggregated)
        if "subject" not in df.columns:
            raise ValueError(f"Missing 'subject' column in {aggregated}")
    else:
        root = path_or_dir if path_or_dir.is_dir() else path_or_dir.parent
        metrics_root = root / "metrics"
        rows = []
        if metrics_root.is_dir():
            for subject_dir in sorted(metrics_root.iterdir()):
                metrics_file = subject_dir / "metrics.csv"
                if not metrics_file.is_file():
                    continue
                df_sub = pd.read_csv(metrics_file)
                df_sub["subject"] = subject_dir.name
                rows.append(df_sub)
        if not rows:
            raise FileNotFoundError(
                f"Could not locate fsqc_metrics.csv or metrics/*.csv under {path_or_dir}"
            )
        df = pd.concat(rows, ignore_index=True)
    df["cohort"] = cohort
    df["subject"] = df["subject"].astype(str)
    parts = df["subject"].str.split("_", expand=False)
    df["site"] = [p[0] if p else "" for p in parts]
    df["timepoint"] = [p[-1] if p else "" for p in parts]
    df["participant"] = ["_".join(p[1:-1]) if len(p) > 2 else (p[1] if len(p) > 1 else "") for p in parts]
    return df


def _zscore(series: pd.Series, group: pd.Series) -> pd.Series:
    grouped = series.groupby(group)
    mean = grouped.transform("mean")
    std = grouped.transform("std").replace({0: pd.NA})
    z = (series - mean) / std
    return z.fillna(0)


def _evaluate_rules(df: pd.DataFrame, config: Dict) -> pd.DataFrame:
    reasons: List[List[str]] = [[] for _ in range(len(df))]
    reason_order = config.get("reporting", {}).get("flag_priority", [])

    def add_reason(mask: pd.Series, reason: str):
        for idx in mask.index[mask]:
            reasons[idx].append(reason)

    # Metric rules
    for rule in config.get("metric_rules", []):
        columns = rule.get("columns", [])
        reason = rule.get("reason", "metric_flag")
        min_abs = rule.get("min_absolute")
        max_abs = rule.get("max_absolute")
        min_z = rule.get("min_zscore")
        max_z = rule.get("max_zscore")
        for col in columns:
            if col not in df.columns:
                continue
            col_series = df[col]
            if min_abs is not None:
                add_reason(col_series < min_abs, reason)
            if max_abs is not None:
                add_reason(col_series > max_abs, reason)
            if min_z is not None or max_z is not None:
                z = _zscore(col_series, df["cohort"])
                if min_z is not None:
                    add_reason(z < min_z, reason)
                if max_z is not None:
                    add_reason(z > max_z, reason)

    # Outlier rules
    for rule in config.get("outlier_rules", []):
        col = rule.get("column")
        reason = rule.get("reason", "outlier_flag")
        if col not in df.columns:
            continue
        max_abs = rule.get("max_absolute")
        min_abs = rule.get("min_absolute")
        if max_abs is not None:
            add_reason(df[col] > max_abs, reason)
        if min_abs is not None:
            add_reason(df[col] < min_abs, reason)

    df = df.copy()
    df["flag_reasons"] = [sorted(set(r), key=lambda x: reason_order.index(x) if x in reason_order else len(reason_order))
                           for r in reasons]
    ok_label = config.get("status_mapping", {}).get("ok_label", "OK")
    exclude_label = config.get("status_mapping", {}).get("exclude_label", "EXCLUDE")
    df["auto_status"] = [exclude_label if r else ok_label for r in df["flag_reasons"]]
    df["auto_reason_text"] = [",".join(r) for r in df["flag_reasons"]]
    df["review_status"] = ""
    df["review_notes"] = ""
    df["reviewer"] = ""
    df["review_date"] = ""
    return df


def _infer_tag(cohort_config: Dict[str, Dict[str, str]]) -> Tuple[str, Path, Dict[str, Path]]:
    cohort_dirs = {}
    tags = set()
    base_dir = None
    for cohort, cfg in cohort_config.items():
        out_dir = Path(cfg.get("output_dir", "")).expanduser()
        cohort_dirs[cohort] = out_dir
        if out_dir.name and "_" in out_dir.name:
            tags.add(out_dir.name.split("_", 1)[1])
        elif out_dir.name:
            tags.add(out_dir.name)
        if base_dir is None and out_dir.parent:
            base_dir = out_dir.parent
    tag = tags.pop() if len(tags) == 1 else dt.datetime.now().strftime("%Y%m%d")
    return tag, base_dir or Path.cwd(), cohort_dirs


def _default_paths(tag: str, base_dir: Path) -> Dict[str, Path]:
    return {
        "parquet": base_dir / "data" / f"fsqc_{tag}.parquet",
        "csv": base_dir / "data" / f"fsqc_{tag}.csv",
        "decision": base_dir / "qc" / f"fsqc_decision_log_{tag}.csv",
    }


def _write_decision_log(df: pd.DataFrame, path: Path, overwrite: bool) -> None:
    cols = [
        "subject",
        "cohort",
        "site",
        "participant",
        "timepoint",
        "auto_status",
        "auto_reason_text",
        "review_status",
        "review_notes",
        "reviewer",
        "review_date",
    ]
    subset = df[cols]
    if path.exists() and not overwrite:
        print(f"[info] Decision log exists (skip): {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    subset.to_csv(path, index=False)
    print(f"[info] Wrote decision log template → {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze FSQC metrics and apply thresholds")
    parser.add_argument("--patients-metrics")
    parser.add_argument("--controls-metrics")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--output-parquet")
    parser.add_argument("--output-csv", default="")
    parser.add_argument("--decision-log", default="")
    parser.add_argument("--overwrite-decision-log", action="store_true")
    parser.add_argument("--tag", help="Override inferred tag for output filenames")
    args = parser.parse_args()

    config = yaml.safe_load(Path(args.config).read_text())
    tag, base_dir, cohort_dirs = _infer_tag(config.get("cohorts", {}))
    tag = args.tag or tag
    defaults = _default_paths(tag, base_dir)

    if args.patients_metrics:
        patients_source = Path(args.patients_metrics)
    else:
        patients_source = cohort_dirs.get("patients")
        if not patients_source:
            raise ValueError("Patients output_dir missing in config; specify --patients-metrics")
    if args.controls_metrics:
        controls_source = Path(args.controls_metrics)
    else:
        controls_source = cohort_dirs.get("controls")
        if not controls_source:
            raise ValueError("Controls output_dir missing in config; specify --controls-metrics")

    out_parquet = Path(args.output_parquet or defaults["parquet"])
    out_csv = Path(args.output_csv or defaults["csv"]) if (args.output_csv or defaults["csv"]) else None
    decision_log = Path(args.decision_log or defaults["decision"])

    df_pat = _load_metrics(patients_source, cohort="patients")
    df_ctrl = _load_metrics(controls_source, cohort="controls")
    df = pd.concat([df_pat, df_ctrl], ignore_index=True)

    df = _evaluate_rules(df, config)

    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_parquet, index=False)
    print(f"[info] Wrote {len(df)} rows → {out_parquet}")

    if out_csv:
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)
        print(f"[info] Wrote CSV → {out_csv}")

    _write_decision_log(df, decision_log, args.overwrite_decision_log)


if __name__ == "__main__":
    main()
