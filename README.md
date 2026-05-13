# bachelor-thesis-pd-mri

Technical repository for the bachelor thesis workflow on **longitudinal structural MRI analysis in Parkinson's disease**, with scope limited to the imaging pipeline **up to and including the FSQC gate**.

## Scope

This repository documents the technical lane used to prepare and quality-control the MRI data before downstream statistical analysis.

### Included
- DICOM copy / staging workflows
- DICOM-to-NIfTI conversion workflows
- candidate-series review and promotion workflows
- longitudinal FastSurfer execution scripts
- FSQC setup, thresholding, and report-generation scripts
- selected provenance tables and QC decision logs
- representative QC image artifacts

### Not included
- raw DICOM data
- converted NIfTI data
- full FastSurfer output trees
- full FSQC screenshot/dashboard directories
- post-FSQC harmonization and statistical-analysis materials
- partner-owned downstream result scripts and outputs

## Thesis boundary

This repo corresponds to the thesis pipeline through:

1. curation of one validated T1-weighted MRI per subject per timepoint
2. longitudinal FastSurfer processing
3. FSQC-based quality control and decision logging

It does **not** attempt to be a full end-to-end public reproduction package for the final inferential results.

## Cohort scope

The technical workflow here covers the ParkWest structural MRI cohorts used in the thesis:
- Parkinson's disease subjects from Bergen, SUS/Stavanger, and Førde
- healthy controls from Bergen, SUS/Stavanger, and Førde

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
  longitudinal/
  fsqc/
  review_and_promotion/
  utils/

runbook_bachelor_thesis.md
README.md
```

## Main workflow

### 1. Data curation
The curation lane follows this high-level sequence:
1. copy DICOMs into thesis-local staging folders
2. convert DICOMs to NIfTI with `dcm2niix`
3. review candidate T1 series
4. log the selected series
5. promote selected files into the curated longitudinal dataset structure
6. purge / rerun when necessary

### 2. Curated dataset rule
The curated longitudinal dataset uses the structure:

```text
<data>/<Site>/<PatientID>/BL
<data>/<Site>/<PatientID>/3Y
<data>/<Site>/<PatientID>/5Y
```

Each timepoint folder is expected to contain exactly **one** validated T1-weighted NIfTI (`.nii` or `.nii.gz`) for downstream longitudinal processing.

### 3. Longitudinal FastSurfer execution
The main execution files are:
- `scripts/longitudinal/longitudinal_patient.slurm`
- `scripts/longitudinal/run_longitudinal_on_folder.sh`

These define how complete BL/3Y/5Y subject series were submitted and processed on the UiS Slurm cluster.

### 4. FSQC gate
The main FSQC files are:
- `scripts/fsqc/build_fsqc_subject_links.py`
- `scripts/fsqc/run_fsqc_on_cohort.sh`
- `scripts/fsqc/run_fsqc_on_cohort.slurm`
- `scripts/fsqc/fsqc_analyze.py`
- `scripts/fsqc/fsqc_report.py`
- `config/fsqc_thresholds.yaml`

These files define the QC input assembly, cohort-level execution, threshold application, and report-generation steps used for the thesis FSQC gate.

## Why the script naming is uneven

Some scripts are specific to individual sites, timepoints, or rerun batches. This is intentional.

The underlying archive structure, naming conventions, and later rerun/recovery needs differed across BL, 3Y, and 5Y. As a result, the script set reflects the **actual operational history** of the curation work rather than an idealized symmetric pipeline.

## Key provenance files

### Selection logs
These document curation decisions made before longitudinal processing, for example:
- `metadata/selection_logs/bergen_bl3y_selection_log.csv`
- `metadata/selection_logs/bergen_healthy_selection_log.csv`
- `metadata/selection_logs/forde_selection_log.csv`
- `metadata/selection_logs/sus_selection_log.csv`

### Conversion reports
These summarize selected conversion-related inventory/provenance for the SUS cohorts:
- `metadata/conversion_reports/sus_converted_series.csv`
- `metadata/conversion_reports/sus_healthy_converted_series.csv`

### QC decision log
- `metadata/qc/fsqc_decision_log_20260428.csv`

This file records FSQC decisions tied to the QC gate used in the thesis workflow.

## Representative artifacts

Representative QC figures are stored under:
- `artifacts/fsqc_good_samples/`
- `artifacts/fsqc_bad_samples/`

These are intended as compact examples of approved and excluded outputs, not as a full export of all QC screenshots.

## Environment assumptions

This workflow was designed around the UiS Unix / Slurm cluster environment.

Important assumptions include:
- Slurm-based execution for FastSurfer jobs
- Singularity-based FastSurfer container execution
- a project-local FreeSurfer license
- thesis-local ParkWest directory structure under `/nfs/br1_prosjekt/ParkWest/...`

See `runbook_bachelor_thesis.md` for the detailed environment, path, and workflow notes.

## Notes on reproducibility

This repository is intended to make the **technical preprocessing and QC lane** understandable and auditable.

Because patient imaging data, cluster-local paths, and large derived outputs cannot be published here, the repo should be read as a structured technical record of:
- what was run
- how it was organized
- which scripts governed curation and QC
- which provenance tables support the workflow narrative

## Main reference document

For the detailed operational description of the technical workflow, see:
- `runbook_bachelor_thesis.md`

That file is the best starting point for understanding the practical execution of the pipeline.