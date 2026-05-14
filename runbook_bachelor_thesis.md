# Bachelor Thesis Technical Runbook

## Purpose

This runbook documents the **technical workflow** for the bachelor
thesis imaging pipeline, covering the full lane from raw DICOMs through
final longitudinal statistical analysis:

- Curation of DICOM/NIfTI per timepoint
- Longitudinal FastSurfer processing
- FSQC quality-control gate
- Cohort definition and metadata foundation (Stage A)
- FSQC-approved feature extraction (Stage B)
- Multi-feature longCombat harmonization (Stage C-1)
- Longitudinal LME modelling, FDR-corrected inference, and effect-size
  estimation (Stage C-2)

The pre-analysis lane (curation → FastSurfer → FSQC) runs on the UiS
Slurm cluster; the post-FSQC analysis lane (Stages A / B / C) runs on
the Gorina1 server. Environment differences between the two are
documented in §1 and §12.

---

## 1. Execution environment

### Two compute environments

The thesis pipeline spans two different servers:

| Lane | Server | Account |
|---|---|---|
| Curation → FastSurfer → FSQC (§§2–11) | UiS Slurm cluster | `u216087` |
| Stages A / B / C (§§12–15) | Gorina1 | `u258907` |

Outputs from the UiS lane are written to the ParkWest NFS share at
`/nfs/br1_prosjekt/ParkWest/...` and read from there by the Gorina1
lane.

### UIS cluster overview

- UiS Unix / Slurm cluster is the supported execution site for §§2–11.
- Cluster login / Slurm account used for this work: `u216087`.
- FastSurfer processing was run through Slurm rather than interactive
  Singularity sessions.
- Proven GPU settings:
  - `--partition=gpu`
  - `--gres=gpu:1`
  - `--cpus-per-task=8`
  - `--mem=32G`
  - `--time=12:00:00`
- `uenv verbose singularity-4.2.2` must be invoked before running
  `long_fastsurfer.sh` inside Slurm scripts.

### Core directories (UiS lane)

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

The Gorina1 directories used by Stages A / B / C are documented in §12.

---

## 2. Dataset layout

All curated subjects were organized in the same folder structure:

```text
.../pd_thesis/data/<Site>/<PatientID>/BL
.../pd_thesis/data/<Site>/<PatientID>/3Y
.../pd_thesis/data/<Site>/<PatientID>/5Y
```

Rules:

- Each timepoint folder must contain exactly **one** validated T1 `.nii`
  or `.nii.gz` for longitudinal processing.
- The final curated dataset only promotes timepoints that satisfy this
  one-volume rule.

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

Prepare one validated T1-weighted MRI input per subject per timepoint
for BL / 3Y / 5Y, in a uniform folder structure suitable for FastSurfer
LONG.

### High-level workflow

1. Copy DICOMs into staging folders.
2. Convert DICOM to NIfTI with `dcm2niix`.
3. Review candidate series.
4. Log the selected series.
5. Promote selected NIfTIs into the final curated dataset structure.
6. Purge / rerun when necessary.

### Main script families

Timepoint-specific scripts exist because the underlying archive
structure, naming conventions, and later rerun/recovery needs differed
across BL, 3Y, and 5Y; these scripts document the actual operational
history of curation rather than an idealized symmetric pipeline.

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
- Raw IDs like `B05` must be normalized to the final dataset
  convention, e.g. `B005`.

### Bergen healthy-control naming quirk

- Bergen BL healthy-control filenames use `BK_XX`.
- Bergen 3Y / 5Y use `BKXX`.
- Baseline lookup must be normalized before intersecting or copying IDs.

### SUS healthy baseline rerun note

For the 18-subject `SUS_healthy` BL rerun/debug set:

- Prefer the converted `PARKVEST_T1W_3D_FFE` series.
- Do **not** prefer the thick-slice `T1W_TSE_TRAN` series.
- The candidate-ranking scripts can mis-rank this because the good 3D
  FFE series may be tagged `SECONDARY`.

---

## 6. Longitudinal FastSurfer execution

### Core execution scripts

- `scripts/longitudinal_patient.slurm`
- `scripts/run_longitudinal_on_folder.sh`

### `longitudinal_patient.slurm`

Purpose:

- run one subject's BL / 3Y / 5Y time series through `long_fastsurfer.sh`

Key behavior:

- expects `SITE` and `PATIENT` environment variables
- reads the first `*.nii*` inside each BL / 3Y / 5Y directory
- writes outputs to:
  - `outputs/fastsurfer_longitudinal/<Site>/<Patient>`
- binds the ParkWest project tree into the container

### `run_longitudinal_on_folder.sh`

Purpose:

- iterate through all subjects in a site folder and submit one Slurm job
  per complete subject

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
   - FastSurfer errors with "Flag '-s' unrecognized".

5. **Parallel settings**
   - The tested-safe settings are:
     - `--parallel 4`
     - `--parallel_surf 4`
     - `--parallel_seg 1`

6. **Interactive Singularity failures**
   - Login-node interactive runs can fail with GLIBC version errors.
   - Prefer batch jobs.

7. **GPU memory failures**
   - If `FastSurferCNN/run_prediction.py` reports insufficient GPU
     memory, rerun with CPU view aggregation
     (`--viewagg_device cpu`) and document the change.

---

## 9. FSQC gate

### Goal

Apply a reproducible quality-control gate to the longitudinal FastSurfer
outputs before any downstream extraction or interpretation.

### Environment

- Run from the UIS cluster.
- FSQC environment used in this project:
  - `micromamba activate fsqc`
- FreeSurfer sourcing is not required unless optional FSQC modules are
  reintroduced.

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
3. Submit patient and control cohort jobs via
   `run_fsqc_on_cohort.slurm`.
4. Regenerate dashboards if needed with `run_fsqc --group-only`.
5. Run `fsqc_analyze.py` and `fsqc_report.py` to apply thresholds and
   export decision summaries.

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

## 12. Post-FSQC environment (Gorina1)

The post-FSQC pipeline runs on the Gorina1 server. Inputs come from the
FastSurfer longitudinal output tree on the ParkWest NFS share.

### Server and directory layout

**Pipeline code** (scripts + shared library):

```text
/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pipeline/
└── parkwest_pipeline/
    ├── pipeline_lib/                          # shared modules
    └── scripts/
        ├── stage_a_metadata/                  # 00–04
        ├── stage_b_extract/                   # 05
        └── stage_c_analyze/                   # 06, 07
```

**FastSurfer input tree** (read-only; produced by the UiS lane and
written to the shared NFS):

```text
/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/outputs/fastsurfer_longitudinal/
├── Bergen/                  # Bergen Pasienter
├── Bergen_healthy/          # Bergen Kontroller
├── SUS/                     # SUS Pasienter
├── SUS_healthy/             # SUS Kontroller
├── Forde/                   # Forde Pasienter
└── Forde_healthy/           # Forde Kontroller
```

**Pipeline outputs** (written by Stages A / B / C):

```text
/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pipeline/outputs/
├── stage_a_metadata/
├── stage_b_extract/
└── stage_c_analyze/
```

### Two Python virtualenvs

The post-FSQC pipeline uses two Python environments. The split exists
because `rpy2` and `pyreadstat` have incompatible Python-base
requirements:

| Virtualenv | Python | Activation | Used by | Purpose |
|---|---|---|---|---|
| `~/.venv/rpy2_clean` | 3.9.25 | `source ~/.venv/rpy2_clean/bin/activate` | Every script except 03 | `rpy2 3.6.7` + scientific stack (`pandas 2.3.3`, `numpy 2.0.2`, `scipy 1.13.1`, `matplotlib 3.9.4`, `statsmodels 0.14.6`, `pydicom 2.4.4`). |
| `~/.venv/pd_thesis` | 3.10.11 | `uenv miniconda3-py310` + `source ~/.venv/pd_thesis/bin/activate` | `03_clinical_covariates.py` only | `pyreadstat 1.3.4`, `pandas 2.3.3`. |

The `pd_thesis` venv produces `clinical_covariates.csv` from the SPSS
register.

### rpy2 + R activation

The rpy2 venv needs Python and R on the PATH simultaneously. The helper
script `~/rpy2_env.sh` handles this:

- activates `~/.venv/rpy2_clean`
- pins R 4.5.2 to the front of PATH
- exports `R_HOME` and `LD_LIBRARY_PATH` so rpy2 can find R's shared
  libraries
- sets `RPY2_CFFI_MODE=ABI`
- creates a writable `$TMPDIR=$HOME/tmp`

### Working-directory activation wrapper

A second helper script in the working directory chains `~/rpy2_env.sh`
with the absolute-path environment variables used in §13–§15. It is
named `activate_refactor.sh` and does the following:

- sources `~/rpy2_env.sh`
- exports `$WORK` = pipeline root on Gorina1
- exports `$OLD` = directory holding the upstream FastSurfer outputs
  and the existing DICOM inventory
- exports `$FS_ROOT` = FastSurfer longitudinal tree
- exports `$CLINICAL_SAV` = ParkVest SPSS register
- exports `$DICOM_INVENTORY` = path to an existing
  `dicom_inventory_all.csv`
- exports `$OUT_A`, `$OUT_B` = Stage A / Stage B output directories
- prepends `$WORK/parkwest_pipeline` to `PYTHONPATH`

Source it once per shell session:

```bash
source activate_refactor.sh
```

### R installation

R 4.5.2 is at `/opt/R-4.5.2_centos9/`. The CRAN packages needed by
Stage C-1 are `longCombat`, `nlme`, and `invgamma`. Install once
inside the rpy2_clean environment:

```r
install.packages(c("longCombat", "nlme", "invgamma"))
```

---

## 13. Cohort and metadata foundation (Stage A)

### Goal

Define the analysis cohort, build a per-scan metadata table with batch
labels for harmonization, and extract per-subject clinical covariates
(age at baseline, sex, group) from the ParkVest SPSS register.

### Stage-A scripts

- `scripts/analysis/stage_a_metadata/00_cohort.py`
- `scripts/analysis/stage_a_metadata/01_dicom_inventory.py`
- `scripts/analysis/stage_a_metadata/02_scanner_metadata.py`
- `scripts/analysis/stage_a_metadata/03_clinical_covariates.py`
- `scripts/analysis/stage_a_metadata/04_merge_metadata.py`

### 13.1 Cohort definition — `00_cohort.py`

Walks the FastSurfer longitudinal output tree and accepts every subject
whose three timepoints (BL, 3Y, 5Y) each contain a complete `stats/`
subdirectory with the four required FastSurfer stats files.

```bash
python $WORK/parkwest_pipeline/scripts/stage_a_metadata/00_cohort.py \
    --fs-root $FS_ROOT \
    --output  $OUT_A/cohort.csv
```

Output: `cohort.csv` (Site, Group, Subject_raw, Subject, FS_dir).

### 13.2 DICOM inventory — `01_dicom_inventory.py`

Walks the ParkWest DICOM tree and produces a per-(Subject × Timepoint ×
Series) inventory. Runtime is approximately 30–60 minutes.

```bash
python $WORK/parkwest_pipeline/scripts/stage_a_metadata/01_dicom_inventory.py \
    --output-dir $OUT_A/
```

Output: `dicom_inventory_all.csv`.

### 13.3 Scanner metadata — `02_scanner_metadata.py`

Filters the inventory to the structural 3D T1 series FastSurfer
consumed, joins to the cohort, and assigns the 6-level `BatchID` per
(Site × Protocol).

```bash
python $WORK/parkwest_pipeline/scripts/stage_a_metadata/02_scanner_metadata.py \
    --inventory $OUT_A/dicom_inventory_all.csv \
    --cohort    $OUT_A/cohort.csv \
    --output    $OUT_A/scanner_metadata.csv
```

Output: `scanner_metadata.csv` (one row per scan with `BatchID`,
`ScannerID`, and `in_cohort` flag).

### 13.4 Clinical covariates — `03_clinical_covariates.py`

Reads the ParkVest SPSS register and extracts per-subject age at
baseline, sex, and group. Runs in the `pd_thesis` venv.

```bash
uenv miniconda3-py310
source ~/.venv/pd_thesis/bin/activate
python $WORK/parkwest_pipeline/scripts/stage_a_metadata/03_clinical_covariates.py \
    --sav     $CLINICAL_SAV \
    --cohort  $OUT_A/cohort.csv \
    --output  $OUT_A/clinical_covariates.csv
deactivate
source activate_refactor.sh
```

Optional flags:

- `--fam <plink.fam>` — cross-check sex against a PLINK family file.
- `--explore` — dump the SPSS file's columns and value labels.

Output: `clinical_covariates.csv` (one row per subject).

### 13.5 Master metadata merge — `04_merge_metadata.py`

Merges scanner metadata with clinical covariates. Renames
`sex_clinical` → `PatientSex_clinical` and `group_clinical` →
`Group_clinical`.

```bash
python $WORK/parkwest_pipeline/scripts/stage_a_metadata/04_merge_metadata.py \
    --scanner-metadata $OUT_A/scanner_metadata.csv \
    --clinical         $OUT_A/clinical_covariates.csv \
    --output           $OUT_A/scanner_metadata_with_covariates.csv
```

Output: `scanner_metadata_with_covariates.csv`.

---

## 14. FSQC-approved feature extraction (Stage B)

### Goal

Walk every cohort subject's `*.stats` files, parse the volumetric and
cortical-thickness measurements, and compute per-subject deltas over
the 12-ROI thesis panel.

### Stage-B script

- `scripts/analysis/stage_b_extract/05_extract.py`

### The 12-ROI thesis panel

Defined in `scripts/analysis/pipeline_lib/constants.py`
(`ROI_PANEL_THESIS`):

- **Subcortical (6, volume, eTIV-normalized):** Hippocampus, Amygdala,
  Caudate, Putamen, Thalamus, Accumbens.
- **Cortical (6, thickness, no ICV correction):** Entorhinal,
  Parahippocampal, Precuneus, Lingual, Caudal Anterior Cingulate,
  Superior Frontal.

### Canonical Stage B command

```bash
python $WORK/parkwest_pipeline/scripts/stage_b_extract/05_extract.py \
    --input-dir  $FS_ROOT \
    --output-dir $OUT_B/ \
    --full-output
```

Optional flag:

- `--full-output` — emits per-site-group and per-subject splits plus
  `cohort_long.csv`.

### Outputs

- `cohort_wide.csv` — one row per (Subject × Timepoint) with all
  FreeSurfer-derived metrics.
- `cohort_cortical_regions.csv` — one row per
  (Subject × Timepoint × Hemisphere × StructName), all DKT atlas regions.
- `subject_roi_deltas.csv` — one row per (Subject × ROI × Delta_Window)
  for the 12-ROI panel across the BL→3Y and BL→5Y windows.
- `roi_level_summary.csv` — aggregated ROI-level statistics.
- `skip_log.txt` — subjects rejected or warned during extraction.

With `--full-output`: `cohort_long.csv` plus per-site-group and
per-subject splits.

---

## 15. Harmonization and longitudinal analysis (Stage C)

### Goal

Remove scanner / protocol batch effects via multi-feature longCombat
harmonization (per tissue), then fit longitudinal linear mixed-effects
models to estimate per-ROI atrophy trajectories and the PD-vs-Control
difference at the 5-year follow-up.

### Stage-C scripts

- `scripts/analysis/stage_c_analyze/06_harmonize.py` (Stage C-1)
- `scripts/analysis/stage_c_analyze/07_analysis.py` (Stage C-2)

### 15.1 Multi-feature longCombat harmonization — `06_harmonize.py`

Runs longCombat per tissue block (V = 12 multi-feature call) via rpy2.
Harmonization formula:
`Years_from_BL * Group_clinical + age_at_BL + PatientSex_clinical`.
Bilateral aggregation runs after harmonization.

```bash
python $WORK/parkwest_pipeline/scripts/stage_c_analyze/06_harmonize.py \
    --cohort-wide       $OUT_B/cohort_wide.csv \
    --cortical-regions  $OUT_B/cohort_cortical_regions.csv \
    --scanner-metadata  $OUT_A/scanner_metadata_with_covariates.csv \
    --output-dir        $WORK/outputs/stage_c_analyze/ \
    --cohort-only
```

Flag notes:

- `--cohort-only` filters to `in_cohort=True` rows from scanner metadata.
- `--preserve` defaults to
  `Years_from_BL age_at_BL PatientSex_clinical Group_clinical`.

### Stage C-1 outputs

Written to `$WORK/outputs/stage_c_analyze/`:

- `cohort_wide_harmonized.csv` — 31 cols, panel only.
- `cohort_cortical_regions_harmonized.csv` — 12 panel hemisphere pairs.
- `subject_roi_deltas_harmonized.csv` — post-harmonization per-subject
  ROI deltas.
- `cohort_long_harmonized.csv` — long format, joined with scanner
  metadata.
- `cohort_long_unharmonized.csv` — Risk-A reference.
- `interval_summary.csv` — per-subject elapsed-time table.
- `longcombat_diagnostics.csv` — per-feature longCombat status,
  batch-R² before/after, feature lists.
- `longcombat_summary.txt` — compact text summary.

### 15.2 Longitudinal LME analysis — `07_analysis.py`

Fits the two LME models, applies per-tissue BH-FDR, computes age- and
sex-adjusted Cohen's d at the 5-year window, runs the Risk-A
harmonization sensitivity check, computes LME diagnostics, and writes
the figures and result tables.

```bash
python $WORK/parkwest_pipeline/scripts/stage_c_analyze/07_analysis.py \
    --cohort-long $WORK/outputs/stage_c_analyze/cohort_long_harmonized.csv \
    --output-dir  $WORK/outputs/stage_c_analyze/
```

`cohort_long_unharmonized.csv` is auto-discovered as a sibling of the
harmonized input.

#### LME design

- **Test A — within-group trajectories.** For each group (Pasient,
  Control) separately: `Value ~ Years_from_BL + age_at_BL + PatientSex_clinical
  + (1 + Years_from_BL | Subject) + (1 | BatchID)`. Produces 24 rows
  (12 ROIs × 2 groups).
- **Test B — between-group trajectories.** Full cohort:
  `Value ~ Years_from_BL * Group_clinical + age_at_BL + PatientSex_clinical
  + (1 + Years_from_BL | Subject) + (1 | BatchID)`. The
  `Years_from_BL:Group_clinical` interaction term is the inferential
  effect. Produces 12 rows.

#### FDR scheme

Per-tissue Benjamini-Hochberg. Subcortical (6 tests) and cortical
(6 tests) are treated as separate inferential families. Output columns:
`p_fdr_tissue`, `sig_star_tissue`.

#### Effect-size derivation

For each ROI: model-predicted PD − Control difference at 5 years
(`pred_diff_5y`) divided by the pooled within-group SD at baseline,
adjusted for age and sex. Reported as `d_adj_5y` (95 % CI in
`d_adj_5y_lo`, `d_adj_5y_hi`).

### Stage C-2 outputs

CSV tables:

- `stats_lme_between_group.csv` — 12 rows, 42 columns; Test B results.
- `stats_lme_within_group.csv` — 24 rows, 26 columns; Test A results.
- `roi_display_order.csv` — `|d_adj_5y|` ordering used for figure axes.
- `risk_a_all_rois.csv` — per-ROI harmonized-vs-unharmonized
  comparison (sign agreement, β ratio, FDR-significance preservation).
- `lme_sanity_check.txt` — Risk-A all-ROI summary.
- `lme_diagnostics_test_b.csv` — per-ROI Test B model diagnostics
  (residual normality, random-effect variances, fit warnings).
- `demographics.csv` — 9 rows × 3 cols (variable, PD, Control) covering
  n, age_at_BL (mean / SD), sex (female / male counts), and follow-up
  time (3Y / 5Y means and SDs).

Figures:

- `figure1_forest_d_adj_5Y.png` — forest plot of effect sizes by ROI.
- `figure2_trajectory_fdrsig.png` — predicted-trajectory panel for the
  FDR-significant ROIs.
- `figure2_trajectory_fdrsig_clean.png` — clean variant.
- `figure3_trajectory_subcortical_supp.png` — all six subcortical
  trajectories.
- `figure3_trajectory_subcortical_supp_clean.png` — clean variant.
- `figure4_trajectory_cortical_supp.png` — all six cortical trajectories.
- `figure4_trajectory_cortical_supp_clean.png` — clean variant.
- `figure5_heatmap_per_group_atrophy.png` — per-group atrophy heatmap.
- `figure_risk_a_scatter_supp.png` — Risk-A scatter (harmonized vs
  unharmonized β).
- `figure_lme_diagnostics_test_b_supp.png` — LME diagnostics panel.