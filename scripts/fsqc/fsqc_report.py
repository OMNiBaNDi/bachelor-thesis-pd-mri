#!/usr/bin/env python3
"""Generate FSQC summary report + partner export."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import yaml

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = WORKSPACE_ROOT / "config/fsqc_thresholds.yaml"


def _load_status(df: pd.DataFrame, decision_log: Path | None) -> pd.DataFrame:
    if decision_log and decision_log.exists():
        log_df = pd.read_csv(decision_log)
        merge_cols = ["subject", "cohort", "site", "participant", "timepoint"]
        df = df.merge(log_df, on=merge_cols, how="left", suffixes=("", "_log"))
        final_status = []
        final_reason = []
        for _, row in df.iterrows():
            if isinstance(row.get("review_status"), str) and row["review_status"].strip():
                final_status.append(row["review_status"].strip())
                note = row.get("review_notes", "")
                final_reason.append(note if isinstance(note, str) else "")
            else:
                final_status.append(row["auto_status"])
                final_reason.append(row.get("auto_reason_text", ""))
        df["status"] = final_status
        df["status_reason"] = final_reason
    else:
        df["status"] = df["auto_status"]
        df["status_reason"] = df.get("auto_reason_text", "")
    return df


def _format_metric_table(df: pd.DataFrame, metrics: List[str]) -> str:
    rows = ["| Cohort | Metric | Mean | SD | Min | Max |", "| --- | --- | --- | --- | --- | --- |"]
    for cohort, cohort_df in df.groupby("cohort"):
        for metric in metrics:
            if metric not in cohort_df.columns:
                continue
            series = cohort_df[metric]
            rows.append(
                f"| {cohort} | {metric} | {series.mean():.2f} | {series.std():.2f} | {series.min():.2f} | {series.max():.2f} |"
            )
    return "\n".join(rows)


def _status_counts_table(df: pd.DataFrame) -> str:
    rows = ["| Cohort | Status | Count |", "| --- | --- | --- |"]
    counts = df.groupby(["cohort", "status"]).size().reset_index(name="count")
    for _, row in counts.iterrows():
        rows.append(f"| {row['cohort']} | {row['status']} | {int(row['count'])} |")
    return "\n".join(rows)


def _exclusion_table(df: pd.DataFrame) -> str:
    excluded = df[df["status"].str.upper() == "EXCLUDE"]
    if excluded.empty:
        return "All subjects cleared FSQC (no exclusions)."
    rows = ["| Cohort | Subject | Timepoint | Reason |", "| --- | --- | --- | --- |"]
    for _, row in excluded.sort_values(["cohort", "subject", "timepoint"]).iterrows():
        rows.append(
            f"| {row['cohort']} | {row['subject']} | {row['timepoint']} | {row['status_reason'] or row['auto_reason_text']} |"
        )
    return "\n".join(rows)


def _infer_paths(config: Dict, tag_override: str | None = None) -> Tuple[str, Path, Path, Path, Path]:
    cohort_cfg = config.get("cohorts", {})
    out_dirs = [Path(cfg.get("output_dir", "")).expanduser() for cfg in cohort_cfg.values() if cfg.get("output_dir")]
    base_dir = out_dirs[0].parent if out_dirs else Path.cwd()
    names = {p.name for p in out_dirs}
    tag = tag_override
    if not tag:
        suffixes = {name.split("_", 1)[1] for name in names if "_" in name}
        tag = suffixes.pop() if len(suffixes) == 1 else dt.datetime.now().strftime("%Y%m%d")
    data_parquet = base_dir / "data" / f"fsqc_{tag}.parquet"
    decision_log = base_dir / "qc" / f"fsqc_decision_log_{tag}.csv"
    report_md = base_dir / "reports" / f"fsqc_{tag}.md"
    status_csv = base_dir / "exports" / f"fsqc_status_{tag}.csv"
    return tag, data_parquet, decision_log, report_md, status_csv


def main() -> None:
    parser = argparse.ArgumentParser(description="Build FSQC report + status export")
    parser.add_argument("--data-parquet")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--decision-log", default="")
    parser.add_argument("--report-md")
    parser.add_argument("--status-export")
    parser.add_argument("--tag", help="Override inferred tag")
    args = parser.parse_args()

    config = yaml.safe_load(Path(args.config).read_text())
    tag, default_parquet, default_decision, default_report, default_export = _infer_paths(config, args.tag)

    data_parquet = Path(args.data_parquet or default_parquet)
    decision_log = Path(args.decision_log or default_decision)
    report_md = Path(args.report_md or default_report)
    status_export = Path(args.status_export or default_export)

    df = pd.read_parquet(data_parquet)
    df = _load_status(df, decision_log if decision_log.exists() else None)

    timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    metrics = config.get("reporting", {}).get("cohort_summary_metrics", [])

    report_lines = []
    report_lines.append("# FSQC Summary")
    report_lines.append("")
    report_lines.append(f"- Generated: {timestamp}")
    report_lines.append(f"- Total subjects: {len(df)}")
    report_lines.append(f"- Threshold config: {Path(args.config)} (version {config.get('version', 'n/a')})")
    report_lines.append("")

    report_lines.append("## Cohort metrics")
    report_lines.append(_format_metric_table(df, metrics))
    report_lines.append("")

    report_lines.append("## Status counts")
    report_lines.append(_status_counts_table(df))
    report_lines.append("")

    report_lines.append("## Exclusions")
    report_lines.append(_exclusion_table(df))
    report_lines.append("")

    report_lines.append("## Notes")
    report_lines.append("- Status labels originate from fsqc_analyze auto flags unless overridden in fsqc_decision_log.csv.")
    report_lines.append("- Screenshot evidence lives under the cohort output directories in outputs/fsqc/.")

    report_md.parent.mkdir(parents=True, exist_ok=True)
    report_md.write_text("\n".join(report_lines), encoding="utf-8")
    print(f"[info] Report written → {report_md}")

    export_cols = ["subject", "cohort", "site", "participant", "timepoint", "status", "status_reason"]
    export_df = df[export_cols].copy()
    export_df.rename(columns={"status": "fsqc_status", "status_reason": "fsqc_reason"}, inplace=True)
    status_export.parent.mkdir(parents=True, exist_ok=True)
    export_df.to_csv(status_export, index=False)
    print(f"[info] Status export written → {status_export}")


if __name__ == "__main__":
    main()
