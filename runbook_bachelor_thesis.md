# Bachelor Thesis Technical Runbook

## Purpose

This runbook documents the **technical workflow** for the bachelor thesis imaging pipeline up to and including the **FSQC gate**.

---

## 1. Execution environment

### UIS cluster overview

- UiS Unix / Slurm cluster is the supported execution site.
- Cluster login / Slurm account used for this work: `u216087`.
- FastSurfer processing was run through Slurm rather than interactive Singularity sessions.
- Proven GPU settings:
  - `--partition=gpu`
  - `--gres=gpu:1`
  - `--cpus-per-task=8`
  - `--mem=32G`
  - `--time=12:00:00`
- `uenv verbose singularity-4.2.2` must be invoked before running `long_fastsurfer.sh` inside Slurm scripts.

### Core directories

- Project root:
  - `/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis`
- Useful path variables:
  - `PW=/nfs/br1_prosjekt/ParkWest`
  - `DATA=${PW}/user/2026vae/AmundEspen/pd_thesis/data`
  - `SCRIPTS=${PW}/user/2026vae/AmundEspen/pd_thesis/scripts`
  - `OUT=${PW}/user/2026vae/AmundEspen/pd_thesis/outputs`
- FastSurfer container:
  - `${PW}/user/2026vae/AmundEspen/pd_thesis/containers/fastsurfer_sm120_sandbox`
- FreeSurfer license:
  - `${PW}/user/2026vae/AmundEspen/pd_thesis/fastsurfer/license/license.txt`
- `dcm2niix` binary:
  - `${PW}/user/2026vae/AmundEspen/pd_thesis/tools/dcm2niix/dcm2niix`

---

## 2. Dataset layout

All curated subjects were organized in the same folder structure:

```text
.../pd_thesis/data/<Site>/<PatientID>/BL
.../pd_thesis/data/<Site>/<PatientID>/3Y
.../pd_thesis/data/<Site>/<PatientID>/5Y
```

Rules:

- Each timepoint folder must contain exactly **one** validated T1 `.nii` or `.nii.gz` for longitudinal processing.
- The final curated dataset only promotes timepoints that satisfy this one-volume rule.

Examples:

- `data/Forde/F02/{BL,3Y,5Y}`
- `data/SUS/S02/{BL,3Y,5Y}`
- `data/Bergen/B005/{BL,3Y,5Y}`

---

## 3. Cohort scope for this repo/runbook

### Included cohorts

- Core ParkWest Parkinson cohort:
  - Bergen
  - SUS / Stavanger
  - Førde
- ParkWest healthy controls:
  - `SUS_healthy`
  - `Bergen_healthy`
  - `Forde_healthy`

---

## 4. Data curation workflow

### Goal

Prepare one validated T1-weighted MRI input per subject per timepoint for BL / 3Y / 5Y, in a uniform folder structure suitable for FastSurfer LONG.

### High-level workflow

1. Copy DICOMs into staging folders.
2. Convert DICOM to NIfTI with `dcm2niix`.
3. Review candidate series.
4. Log the selected series.
5. Promote selected NIfTIs into the final curated dataset structure.
6. Purge / rerun when necessary.

### Main script families

Timepoint-specific scripts exist because the underlying archive structure, naming conventions, and later rerun/recovery needs differed across BL, 3Y, and 5Y; these scripts document the actual operational history of curation rather than an idealized symmetric pipeline.

#### Copy / staging

- `scripts/copy_bergen_bl.sh`
- `scripts/copy_bergen_3y.sh`
- `scripts/copy_bergen_healthy_dicoms.sh`
- `scripts/copy_forde_dicoms.sh`
- `scripts/copy_forde_healthy_dicoms.sh`
- `scripts/copy_sus_dicoms.sh`
- `scripts/copy_sus_healthy_dicoms.sh`

#### Conversion

- `scripts/convert_bergen_dicom.sh`
- `scripts/convert_bergen_healthy_dicom.sh`
- `scripts/convert_forde_dicom.sh`
- `scripts/convert_forde_healthy_dicom.sh`
- `scripts/convert_sus_dicom.sh`
- `scripts/convert_sus_healthy_dicom.sh`
- `scripts/convert_sus_healthy_3y_batch.sh`
- `scripts/convert_sus_healthy_bl_batch.sh`

#### Review / promotion

- `scripts/review_bergen_bl_3y.py`
- `scripts/promote_bergen_bl_3y.py`
- `scripts/review_bergen_healthy_timepoints.py`
- `scripts/promote_bergen_healthy_timepoints.py`
- `scripts/review_forde_timepoints.py`
- `scripts/promote_forde_timepoints.py`
- `scripts/review_forde_healthy_timepoints.py`
- `scripts/promote_forde_healthy_timepoints.py`
- `scripts/review_sus_timepoints.py`
- `scripts/promote_sus_timepoints.py`
- `scripts/review_sus_healthy_timepoints.py`
- `scripts/promote_sus_healthy_timepoints.py`

#### Purge / reset helpers

- `scripts/purge_bergen_bl_3y.sh`
- `scripts/purge_bergen_healthy_timepoints.sh`
- `scripts/purge_forde_timepoints.sh`
- `scripts/purge_sus_timepoints.sh`
- `scripts/purge_sus_healthy_timepoints.sh`

---

## 5. Important curation notes

### Bergen 5Y selection rules

- Raw DICOM source:
  - `/nfs/br1_prosjekt/ParkWest/ImageData/ParkVest_5Y/Pasienter/PV_patients_Bergen_5Y`
- Preferred T1 candidate characteristics:
  - `SeriesDescription` contains `t1` and `mpr`
  - voxel size close to `1.0 x 1.0 x 1.0`
  - shape near `(160, 256, 256)`
- Exclude:
  - ROI helper series
  - PD/T2 sequences
  - one-slice auxiliary outputs
- Raw IDs like `B05` must be normalized to the final dataset convention, e.g. `B005`.

### Bergen healthy-control naming quirk

- Bergen BL healthy-control filenames use `BK_XX`.
- Bergen 3Y / 5Y use `BKXX`.
- Baseline lookup must be normalized before intersecting or copying IDs.

### SUS healthy baseline rerun note

For the 18-subject `SUS_healthy` BL rerun/debug set:

- Prefer the converted `PARKVEST_T1W_3D_FFE` series.
- Do **not** prefer the thick-slice `T1W_TSE_TRAN` series.
- The candidate-ranking scripts can mis-rank this because the good 3D FFE series may be tagged `SECONDARY`.

---

## 6. Longitudinal FastSurfer execution

### Core execution scripts

- `scripts/longitudinal_patient.slurm`
- `scripts/run_longitudinal_on_folder.sh`

### `longitudinal_patient.slurm`

Purpose:

- run one subject’s BL / 3Y / 5Y time series through `long_fastsurfer.sh`

Key behavior:

- expects `SITE` and `PATIENT` environment variables
- reads the first `*.nii*` inside each BL / 3Y / 5Y directory
- writes outputs to:
  - `outputs/fastsurfer_longitudinal/<Site>/<Patient>`
- binds the ParkWest project tree into the container

### `run_longitudinal_on_folder.sh`

Purpose:

- iterate through all subjects in a site folder and submit one Slurm job per complete subject

Usage examples:

```bash
scripts/run_longitudinal_on_folder.sh ${PW}/user/2026vae/AmundEspen/pd_thesis/data/Forde
scripts/run_longitudinal_on_folder.sh ${PW}/user/2026vae/AmundEspen/pd_thesis/data/SUS
scripts/run_longitudinal_on_folder.sh ${PW}/user/2026vae/AmundEspen/pd_thesis/data/Bergen
```

Behavior:

- detects site from the folder name
- checks for BL / 3Y / 5Y NIfTI presence
- skips incomplete subjects
- submits complete subjects via `sbatch`

### Example manual submission

```bash
sbatch --job-name=${SITE}_${PATIENT}_long \
  --export=ALL,SITE=$SITE,PATIENT=$PATIENT \
  scripts/longitudinal_patient.slurm
```

---

## 7. Job monitoring and verification

### Queue / job inspection

```bash
squeue -u $USER -o "%.10i %.20j %.9P %.8T %.10M %.6D %R"
scontrol show job <JOBID>
sacct -u $USER --format=JobID,JobName%25,Partition,State,ExitCode,Elapsed,Start,End
scancel <JOBID>
```

### Log-based success check

```bash
cd ${PW}/user/2026vae/AmundEspen/pd_thesis/outputs/fastsurfer_longitudinal/<Site>

for p in */*_template/scripts/long_fastsurfer.log; do
  patient=$(basename "$(dirname "$(dirname "$(dirname "$p")")")")
  if grep -q "Full longitudinal processing" "$p"; then
    echo "$patient SUCCESS"
  else
    echo "$patient FAILED"
  fi
done
```

---

## 8. Common FastSurfer execution pitfalls

1. **Wrong partition name**
   - `gpuA100` fails.
   - Use `gpu`.

2. **Missing Slurm log directory**
   - Ensure `outputs/slurm_logs` exists before first submission.

3. **DOS line endings**
   - Windows line endings can break Slurm scripts.
   - Fix with `dos2unix` or `sed -i 's/\r$//' script.slurm`.

4. **Do not pass `-s` to `long_fastsurfer.sh`**
   - FastSurfer errors with “Flag '-s' unrecognized”.

5. **Parallel settings**
   - The tested-safe settings are:
     - `--parallel 4`
     - `--parallel_surf 4`
     - `--parallel_seg 1`

6. **Interactive Singularity failures**
   - Login-node interactive runs can fail with GLIBC version errors.
   - Prefer batch jobs.

7. **GPU memory failures**
   - If `FastSurferCNN/run_prediction.py` reports insufficient GPU memory, rerun with CPU view aggregation (`--viewagg_device cpu`) and document the change.

---

## 9. FSQC gate

### Goal

Apply a reproducible quality-control gate to the longitudinal FastSurfer outputs before any downstream extraction or interpretation.

### Environment

- Run from the UIS cluster.
- FSQC environment used in this project:
  - `micromamba activate fsqc`
- FreeSurfer sourcing is not required unless optional FSQC modules are reintroduced.

### Core FSQC scripts

#### Subject assembly and execution

- `scripts/build_fsqc_subject_links.py`
- `scripts/run_fsqc_chunk.sh`
- `scripts/run_fsqc_on_cohort.sh`
- `scripts/run_fsqc_on_cohort.slurm`

#### Thresholding and reporting

- `scripts/fsqc_analyze.py`
- `scripts/fsqc_report.py`
- `config/fsqc_thresholds.yaml`

#### QC figure / review helpers

- `scripts/fsqc_collect_review_bundle.py`
- `scripts/fsqc_collect_good_samples.py`
- `scripts/fsqc_build_subject_montages.py`
- `scripts/fsqc_make_subject_timepoint_montage.py`
- `scripts/fsqc_make_good_montage.py`
- `scripts/fsqc_metric_snapshot.py`
- `scripts/run_fsqc_collect_good_samples.sh`
- `scripts/run_fsqc_low_snr_figure.sh`
- `scripts/run_fsqc_subject_montages.sh`

### Inputs

- Source tree:
  - `outputs/fastsurfer_longitudinal/<Site>/<Patient>/{BL,3Y,5Y}`
- Symlink tree + subject lists:
  - `outputs/fsqc_subjects/{patients,controls}/`
- Threshold config:
  - `config/fsqc_thresholds.yaml`

### Repeatable FSQC procedure

1. Rebuild FSQC subject links and subject lists.
2. Smoke-test 1–2 subjects with `run_fsqc_chunk.sh`.
3. Submit patient and control cohort jobs via `run_fsqc_on_cohort.slurm`.
4. Regenerate dashboards if needed with `run_fsqc --group-only`.
5. Run `fsqc_analyze.py` and `fsqc_report.py` to apply thresholds and export decision summaries.

### Example cohort submission

```bash
sbatch scripts/run_fsqc_on_cohort.slurm patients 20260428
sbatch scripts/run_fsqc_on_cohort.slurm controls 20260428
```

### Example subject-link build

```bash
python3 scripts/build_fsqc_subject_links.py \
  --source-root ${PW}/user/2026vae/AmundEspen/pd_thesis/outputs/fastsurfer_longitudinal \
  --dest-root   ${PW}/user/2026vae/AmundEspen/pd_thesis/outputs/fsqc_subjects \
  --cohorts patients controls --overwrite
```

---

## 10. FSQC thresholds and decision logic

The threshold policy is defined in:

- `config/fsqc_thresholds.yaml`

Main classes of checks:

- white-matter SNR
- gray-matter SNR
- contrast SNR
- Talairach rotation
- topology / holes / defects
- corpus callosum / fornix outlier behavior
- anatomical outliers

Decision rule:

- threshold failures are automatically flagged
- manual review can confirm exclusion or reinstate a borderline case
- the final QC state is recorded in the decision log

---

## 11. FSQC outputs and artifacts

### Main outputs

- Cohort-level FSQC directories:
  - `outputs/fsqc/patients_20260428/`
  - `outputs/fsqc/controls_20260428/`
- Decision log:
  - `outputs/fsqc/qc/fsqc_decision_log.csv`

### Representative thesis artifacts

- Approved montage:
  - `artifacts/fsqc_good_samples/good_samples_montage.png`
- Excluded montage:
  - `artifacts/fsqc_bad_samples/bad_samples_montage.png`

### Additional provenance files

Useful QC summaries generated during the project include:

- `scripts/fsqc_decision_log_20260428.csv`
- `scripts/fsqc_good_subjects.csv`
- `scripts/fsqc_bad_subjects.csv`
- `scripts/fsqc_low_snr_subjects.csv`
- `scripts/fsqc_subject_montages.json`

---
