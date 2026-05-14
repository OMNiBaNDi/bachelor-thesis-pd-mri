"""Find cohort subjects in the FastSurfer output tree.

A subject is in the cohort if it has all three timepoints (BL, 3Y, 5Y)
and each timepoint has the required stats files. Everything else is
skipped.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple

import pandas as pd

from pipeline_lib.constants import (
    REQUIRED_STATS_FILES,
    SITE_GROUP_DIRS,
    TEMPLATE_SUFFIX,
    TIMEPOINTS,
)
from pipeline_lib.ids import normalize_subject_id


@dataclass(frozen=True)
class SkippedSubject:
    """A subject folder rejected during discovery."""
    dir_name: str           # parent site/group dir
    subject_raw: str        # raw (unnormalized) folder name
    reason: str
    missing_timepoints: Tuple[str, ...] = ()


def discover_cohort(
    fs_root: Path,
    timepoints: Iterable[str] = TIMEPOINTS,
) -> Tuple[pd.DataFrame, List[SkippedSubject]]:
    """Walk the FastSurfer output tree and return the cohort.

    A subject folder is kept if its name doesn't contain TEMPLATE_SUFFIX
    and every timepoint has a populated stats/ folder. Returns
    (cohort_df, skipped): cohort_df has columns Site, Group, Subject_raw,
    Subject, FS_dir; skipped lists folders that were missing one or more
    timepoints.
    """
    timepoints = tuple(timepoints)
    rows: List[dict] = []
    skipped: List[SkippedSubject] = []

    if not fs_root.is_dir():
        raise FileNotFoundError(
            f"FastSurfer root directory does not exist: {fs_root}"
        )

    for dir_name, (site, group) in SITE_GROUP_DIRS.items():
        site_dir = fs_root / dir_name
        if not site_dir.is_dir():
            # Some site/group combos don't exist (e.g. a site without controls).
            continue

        for subj_dir in sorted(site_dir.iterdir()):
            if not subj_dir.is_dir():
                continue
            if TEMPLATE_SUFFIX in subj_dir.name.lower():
                # Skip FastSurfer longitudinal template folders.
                continue

            missing = tuple(
                tp for tp in timepoints
                if not _timepoint_is_populated(subj_dir / tp)
            )
            if missing:
                skipped.append(SkippedSubject(
                    dir_name=dir_name,
                    subject_raw=subj_dir.name,
                    reason=f"missing timepoint(s): {', '.join(missing)}",
                    missing_timepoints=missing,
                ))
                continue

            rows.append({
                "Site":        site,
                "Group":       group,
                "Subject_raw": subj_dir.name,
                "Subject":     normalize_subject_id(subj_dir.name),
                "FS_dir":      str(subj_dir.resolve()),
            })

    cohort = (
        pd.DataFrame(rows, columns=["Site", "Group", "Subject_raw",
                                     "Subject", "FS_dir"])
        .sort_values(["Site", "Group", "Subject"])
        .reset_index(drop=True)
    )
    return cohort, skipped


def _timepoint_is_populated(tp_dir: Path) -> bool:
    """True if tp_dir/stats/ contains every required stats file."""
    if not tp_dir.is_dir():
        return False
    stats_dir = tp_dir / "stats"
    if not stats_dir.is_dir():
        return False
    return all((stats_dir / f).exists() for f in REQUIRED_STATS_FILES)
