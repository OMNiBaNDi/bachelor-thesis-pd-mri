"""
pipeline_lib
============

Shared utilities for the ParkWest longitudinal MRI harmonization pipeline.

This package holds code that more than one numbered CLI script needs:

  - ids.py:       normalize_subject_id(), the canonical subject-ID rule.
  - constants.py: ROI lists, domain mapping, batch thresholds, site/group
                  directory names. One source of truth, imported by every
                  stage that touches these values.
  - cohort.py:    FastSurfer-tree cohort discovery (used by stage A
                  step 00 and, in the future, stage B step 05).

Numbered scripts under scripts/ are CLI entry points only; they import
DOWN into pipeline_lib but never sideways into each other. This keeps
the dependency graph between scripts linear and matches the three-stage
pipeline diagram in the thesis methods chapter.
"""
__version__ = "0.1.0"
