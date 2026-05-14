# bachelor-thesis-pd-mri

Technical repository for the bachelor thesis workflow on **longitudinal
structural MRI analysis in Parkinson's disease**, covering the full
technical pipeline from raw DICOMs through final longitudinal statistical
analysis.

## Scope

This repository documents the complete technical lane used in the thesis,
in two stages:

**Pre-analysis lane** (UiS Slurm cluster):

- DICOM copy and staging
- DICOM-to-NIfTI conversion
- Candidate-series review and promotion
- Longitudinal FastSurfer execution
- FSQC quality-control gate (thresholding, decision logging,
  representative artifacts)

**Post-FSQC analysis lane** (Gorina1 server):

- Stage A — cohort definition + scanner/clinical metadata
- Stage B — feature extraction from FastSurfer outputs and 12-ROI panel
  delta computation
- Stage C-1 — multi-feature longCombat harmonization (per tissue)
- Stage C-2 — longitudinal linear mixed-effects (LME) modelling,
  per-tissue BH-FDR, effect-size estimation, Risk-A harmonization
  sensitivity check, figures and result tables

### Included

- All scripts that drive the pipeline end-to-end
- Selection logs, conversion reports, and the FSQC decision log
- Representative QC montages for the FSQC gate
- Threshold configuration for FSQC
- `pipeline_lib/` — the shared Python library used by the post-FSQC
  stages (subject-ID normalization, cohort discovery, ROI panel
  definitions, delta computation)

### Not included

- Raw DICOM data
- Converted NIfTI volumes
- Full FastSurfer output trees
- Full FSQC screenshot / dashboard directories
- Stage A / Stage B / Stage C derived output tables and figures
- Subject-level data

## Thesis boundary

This repository corresponds to the full thesis pipeline:

1. Curation of one validated T1-weighted MRI per subject per timepoint
2. Longitudinal FastSurfer processing
3. FSQC-based quality control and decision logging
4. Cohort definition + scanner / clinical metadata extraction (Stage A)
5. Feature extraction from FastSurfer outputs (Stage B)
6. Multi-feature longCombat harmonization (Stage C-1)
7. Longitudinal LME modelling, FDR-corrected inference, and effect-size
   estimation (Stage C-2)

Raw imaging data, derived outputs, and cluster-local paths are not in
the repo.

## Cohort scope

The technical workflow here covers the ParkWest structural MRI cohorts
used in the thesis:

- Parkinson's disease subjects from Bergen, SUS / Stavanger, and Førde
- Healthy controls from Bergen, SUS / Stavanger, and Førde

## Repository structure

```text
config/
  fsqc_thresholds.yaml

metadata/
  conversion_reports/
  qc/
  selection_logs/

artifacts/
  fsqc_good_samples/
  fsqc_bad_samples/

scripts/
  copy/
  convert/
  review_and_promotion/
  longitudinal/
  fsqc/
  analysis/                                # post-FSQC pipeline
    pipeline_lib/
      ids.py
      cohort.py
      constants.py
      deltas.py
    stage_a_metadata/
      00_cohort.py
      01_dicom_inventory.py
      02_scanner_metadata.py
      03_clinical_covariates.py
      04_merge_metadata.py
    stage_b_extract/
      05_extract.py
    stage_c_analyze/
      06_harmonize.py
      07_analysis.py
  utils/

runbook_bachelor_thesis.md
README.md
```

`scripts/analysis/` preserves the on-disk layout used on Gorina1
(`pipeline_lib/` alongside `scripts/`) so imports of the form
`from pipeline_lib.constants import ...` resolve correctly.

## Main workflow

### 1. Data curation

The curation lane follows this high-level sequence:

1. Copy DICOMs into thesis-local staging folders
2. Convert DICOMs to NIfTI with `dcm2niix`
3. Review candidate T1 series
4. Log the selected series
5. Promote selected files into the curated longitudinal dataset structure
6. Purge / rerun when necessary

### 2. Curated dataset rule

The curated longitudinal dataset uses the structure:

```text
<data>/<Site>/<PatientID>/BL
<data>/<Site>/<PatientID>/3Y
<data>/<Site>/<PatientID>/5Y
```

Each timepoint folder contains exactly **one** validated T1-weighted
NIfTI (`.nii` or `.nii.gz`) for downstream longitudinal processing.

### 3. Longitudinal FastSurfer execution

The main execution files are:

- `scripts/longitudinal/longitudinal_patient.slurm`
- `scripts/longitudinal/run_longitudinal_on_folder.sh`

These define how complete BL / 3Y / 5Y subject series were submitted and
processed on the UiS Slurm cluster.

### 4. FSQC gate

The main FSQC files are:

- `scripts/fsqc/build_fsqc_subject_links.py`
- `scripts/fsqc/run_fsqc_on_cohort.sh`
- `scripts/fsqc/run_fsqc_on_cohort.slurm`
- `scripts/fsqc/fsqc_analyze.py`
- `scripts/fsqc/fsqc_report.py`
- `config/fsqc_thresholds.yaml`

These define the QC input assembly, cohort-level execution, threshold
application, and report-generation steps used for the FSQC gate.

### 5. Cohort and metadata foundation (Stage A)

Once the FSQC-passing subjects are known, Stage A defines the analysis
cohort and assembles per-scan metadata:

- `scripts/analysis/stage_a_metadata/00_cohort.py` — walks the FastSurfer
  longitudinal output tree and writes `cohort.csv` listing every subject
  with all three timepoints complete.
- `scripts/analysis/stage_a_metadata/01_dicom_inventory.py` — walks the
  DICOM tree and produces an inventory of all series.
- `scripts/analysis/stage_a_metadata/02_scanner_metadata.py` — filters
  the inventory to the structural 3D T1 series FastSurfer consumed and
  assigns the 6-level `BatchID` per (Site × Protocol).
- `scripts/analysis/stage_a_metadata/03_clinical_covariates.py` — reads
  the ParkVest SPSS register and extracts per-subject age at baseline,
  sex, and group. Runs in a separate Python environment (see runbook
  §13.2).
- `scripts/analysis/stage_a_metadata/04_merge_metadata.py` — merges
  scanner metadata with clinical covariates.

### 6. Feature extraction (Stage B)

- `scripts/analysis/stage_b_extract/05_extract.py` — walks every
  FastSurfer subject's `*.stats` files, parses volumetric and cortical
  thickness measurements, and computes per-subject deltas over the
  12-ROI thesis panel.

The 12-ROI panel is defined in
`scripts/analysis/pipeline_lib/constants.py` (`ROI_PANEL_THESIS`): six
subcortical volumes (Hippocampus, Amygdala, Caudate, Putamen, Thalamus,
Accumbens; eTIV-normalized) and six cortical thicknesses (Entorhinal,
Parahippocampal, Precuneus, Lingual, Caudal Anterior Cingulate, Superior
Frontal).

### 7. Harmonization and longitudinal analysis (Stage C)

- `scripts/analysis/stage_c_analyze/06_harmonize.py` — runs longCombat
  per tissue block (V = 12 multi-feature call) via rpy2. Bilateral
  aggregation runs after harmonization. Builds the long-format LME
  inputs.
- `scripts/analysis/stage_c_analyze/07_analysis.py` — fits the LME
  models (Test A: within-group trajectories; Test B: between-group
  Years × Group interaction), applies per-tissue Benjamini-Hochberg FDR
  correction, computes age- and sex-adjusted Cohen's d at the 5-year
  window, runs the Risk-A harmonization sensitivity check, computes LME
  diagnostics, and writes the figures and result tables.

Harmonization formula:
`Years_from_BL * Group_clinical + age_at_BL + PatientSex_clinical`.
FDR is applied per tissue (subcortical and cortical treated as
separate families of 6 tests each). Effect sizes are reported as
age- and sex-adjusted Cohen's d at 5 years post-baseline.

### 8. Result summaries and figures

Stage C-2 produces the CSV tables and figures on Gorina1. Outputs are
not pushed to the repo. See runbook §15 for the full output list.

## Why the script naming is uneven

Some scripts are specific to individual sites, timepoints, or rerun
batches. This is intentional.

For the pre-analysis lane, the underlying archive structure, naming
conventions, and later rerun/recovery needs differed across BL, 3Y, and
5Y. As a result, the curation script set reflects the **actual
operational history** of the work rather than an idealized symmetric
pipeline.

For the post-FSQC lane, the script names follow a `NN_purpose.py`
convention (`00_cohort.py` through `07_analysis.py`) reflecting the
stage-A / stage-B / stage-C structure.

## Key provenance files

### Selection logs

These document curation decisions made before longitudinal processing:

- `metadata/selection_logs/bergen_bl3y_selection_log.csv`
- `metadata/selection_logs/bergen_healthy_selection_log.csv`
- `metadata/selection_logs/forde_selection_log.csv`
- `metadata/selection_logs/sus_selection_log.csv`

### Conversion reports

- `metadata/conversion_reports/sus_converted_series.csv`
- `metadata/conversion_reports/sus_healthy_converted_series.csv`

### QC decision log

- `metadata/qc/fsqc_decision_log_20260428.csv`

This file records FSQC decisions tied to the QC gate used in the thesis
workflow.

## Representative artifacts

Representative QC figures are stored under:

- `artifacts/fsqc_good_samples/`
- `artifacts/fsqc_bad_samples/`

These are compact examples of approved and excluded outputs, not a full
export of all QC screenshots.

## Environment assumptions

This workflow was designed around two different compute environments:

**Pre-analysis lane** — UiS Unix / Slurm cluster:

- Slurm-based execution for FastSurfer jobs
- Singularity-based FastSurfer container execution
- Project-local FreeSurfer license
- Thesis-local ParkWest directory structure under
  `/nfs/br1_prosjekt/ParkWest/...`

**Post-FSQC lane** — Gorina1 server:

- Two Python virtualenvs (`rpy2_clean` for stages A/B/C, `pd_thesis`
  for `03_clinical_covariates.py`)
- R 4.5.2 with `longCombat`, `nlme`, and `invgamma` packages
- Thesis-local working directory under
  `/nfs/br1_prosjekt/ParkWest/user/...`

See `runbook_bachelor_thesis.md` for environment, path, and workflow
notes covering both lanes.

## Notes on reproducibility

This repository is intended to make the **technical preprocessing,
QC, and analysis lane** understandable and auditable.

Because patient imaging data, cluster-local paths, and large derived
outputs cannot be published here, the repo should be read as a
structured technical record of:

- what was run
- how it was organized
- which scripts governed curation, QC, and analysis
- which provenance tables support the workflow narrative

## Main reference document

For the detailed operational description of the technical workflow, see:

- `runbook_bachelor_thesis.md`

That file is the best starting point for understanding the practical
execution of both the pre-analysis and post-FSQC lanes.