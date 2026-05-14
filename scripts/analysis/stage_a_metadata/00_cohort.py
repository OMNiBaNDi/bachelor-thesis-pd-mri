v#!/usr/bin/env python3
"""Define the cohort by walking the FastSurfer output tree.

Subjects with all three timepoints (BL, 3Y, 5Y) complete are accepted;
everything else is skipped.

Usage:
    python scripts/stage_a_metadata/00_cohort.py \\
        --fs-root outputs/fastsurfer_longitudinal \\
        --output  outputs/stage_a_metadata/cohort.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pipeline_lib.cohort import discover_cohort


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--fs-root", type=Path, required=True,
        help="FastSurfer longitudinal output root.",
    )
    p.add_argument(
        "--output", type=Path, required=True,
        help="Where to write cohort.csv.",
    )
    args = p.parse_args()

    cohort, skipped = discover_cohort(args.fs_root)

    print(f"FastSurfer root: {args.fs_root}")
    print(f"Discovered {len(cohort)} accepted subjects.\n")

    if len(cohort) == 0:
        print(
            "ERROR: no subjects accepted. Expected layout is "
            "<fs-root>/<SiteGroupDir>/<Subject>/{BL,3Y,5Y}/stats/.",
            file=sys.stderr,
        )
        return 2

    print("Per-site/group breakdown:")
    counts = (
        cohort.groupby(["Site", "Group"])
              .size()
              .rename("n")
              .reset_index()
              .sort_values(["Site", "Group"])
    )
    for _, r in counts.iterrows():
        print(f"  {r['Site']:<8s} {r['Group']:<12s} {r['n']:>3d}")
    print(f"  {'TOTAL':<21s} {len(cohort):>3d}\n")

    if skipped:
        print(f"Rejected {len(skipped)} folder(s):")
        for s in skipped:
            print(f"  [{s.dir_name}] {s.subject_raw}: {s.reason}")
        print()
    else:
        print("No folders rejected.\n")

    # Preview normalized vs raw IDs.
    print("First 10 rows of cohort.csv:")
    preview_cols = ["Site", "Group", "Subject_raw", "Subject"]
    print(cohort[preview_cols].head(10).to_string(index=False))
    print()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    cohort.to_csv(args.output, index=False)
    print(f"Wrote {args.output} ({len(cohort)} rows, "
          f"{len(cohort.columns)} columns)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
