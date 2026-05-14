"""Pipeline-wide constants: ROI panel, file conventions, batch settings.

The 12 analysis ROIs live in ROI_PANEL_THESIS.
"""
from __future__ import annotations

from typing import Dict, List, Tuple


# --- FastSurfer output tree layout ---
#
# Subjects live in per-site folders; controls go in a sibling folder
# with the "_healthy" suffix:
#
#     outputs/fastsurfer_longitudinal/
#         Bergen/           <- Pasienter
#         Bergen_healthy/   <- Kontroller
#         SUS/              <- Pasienter
#         SUS_healthy/      <- Kontroller
#         Forde/            <- Pasienter
#         Forde_healthy/    <- Kontroller

SITE_GROUP_DIRS: Dict[str, Tuple[str, str]] = {
    "Bergen":         ("Bergen", "Pasienter"),
    "Bergen_healthy": ("Bergen", "Kontroller"),
    "SUS":            ("SUS",    "Pasienter"),
    "SUS_healthy":    ("SUS",    "Kontroller"),
    "Forde":          ("Forde",  "Pasienter"),
    "Forde_healthy":  ("Forde",  "Kontroller"),
}

# Timepoints required for a subject to be in the cohort.
TIMEPOINTS: Tuple[str, ...] = ("BL", "3Y", "5Y")

# Nominal years since baseline. Real elapsed time per subject comes
# from DICOM dates via pipeline_lib.deltas.compute_interval_years;
# this is the fallback.
TIMEPOINT_YEARS: Dict[str, float] = {"BL": 0.0, "3Y": 3.0, "5Y": 5.0}

# Suffix on FastSurfer longitudinal template folders, which we ignore.
TEMPLATE_SUFFIX: str = "_template"

# Stats filenames. Named individually because the stage B extractor
# opens each one by name.
ASEG_STATS_FILENAME:       str = "aseg.stats"
ASEG_DKT_STATS_FILENAME:   str = "aseg+DKT.stats"
LH_APARC_FILENAME:         str = "lh.aparc.DKTatlas.mapped.stats"
RH_APARC_FILENAME:         str = "rh.aparc.DKTatlas.mapped.stats"

# Optional. If missing, the QC script NaN-fills the cerebellar columns.
CEREBELLUM_STATS_FILENAME: str = "cerebellum.CerebNet.stats"

# A timepoint counts as complete when its stats/ folder has all four
# of these files.
REQUIRED_STATS_FILES: Tuple[str, ...] = (
    ASEG_STATS_FILENAME,
    ASEG_DKT_STATS_FILENAME,
    LH_APARC_FILENAME,
    RH_APARC_FILENAME,
)


# --- Subcortical structures harmonized per hemisphere ---
#
# Each bilateral structure is harmonized as two longCombat features
# (one per hemisphere). Bilateral means are rebuilt from the harmonized
# halves afterwards.
#
# 06_harmonize.py renames the hyphenated names to underscore form for
# R's formula parser, then restores the hyphens on output.

BILATERAL_SUBCORTICAL_MEANS: Dict[str, Tuple[str, str]] = {
    "Putamen_Mean":     ("Left-Putamen",        "Right-Putamen"),
    "Caudate_Mean":     ("Left-Caudate",        "Right-Caudate"),
    "Hippocampus_Mean": ("Left-Hippocampus",    "Right-Hippocampus"),
    "Amygdala_Mean":    ("Left-Amygdala",       "Right-Amygdala"),
    "Accumbens_Mean":   ("Left-Accumbens-area", "Right-Accumbens-area"),
    "Thalamus_Mean":    ("Left-Thalamus",       "Right-Thalamus"),
}

# Midline structures (no hemisphere split). Not in the panel and not
# harmonized; listed here for reference. The _Total suffix means
# component sum, not bilateral mean.
MIDLINE_SUBCORTICAL: List[str] = [
    "Brain-Stem",
    "Cbm_Vermis_Total",
]

# Cerebellar anterior lobules. Not in the panel and not harmonized.
# Cbm_Anterior_Total is a component sum, not a bilateral mean.
CBM_ANTERIOR_COMPONENTS: List[str] = [
    "Cbm_Left_I_IV",  "Cbm_Right_I_IV",
    "Cbm_Left_V",     "Cbm_Right_V",
]
CBM_ANTERIOR_TOTAL_COL: str = "Cbm_Anterior_Total"

# Columns passed to longCombat: 12 hemispheric volumes (6 panel
# structures x 2 hemispheres).
ALL_SUBCORTICAL_HARMONIZE_COLS: List[str] = [
    c for pair in BILATERAL_SUBCORTICAL_MEANS.values() for c in pair
]

# Bilateral means rebuilt from harmonized halves.
ALL_SUBCORTICAL_MEAN_COLS: List[str] = list(BILATERAL_SUBCORTICAL_MEANS.keys())


# --- Stage C output schemas ---
#
# Columns in cohort_wide_harmonized.csv and
# cohort_cortical_regions_harmonized.csv. Both files only contain
# harmonized values plus identifiers; unharmonized values for non-panel
# structures stay in the stage B files.
#
# Merged scanner metadata (BatchID, age_at_BL, etc.) isn't in the slim
# wide schema; build_cohort_long merges it in later.

COHORT_WIDE_HARMONIZED_OUTPUT_COLS: List[str] = [
    # Identifiers
    "Site", "Group", "Subject", "Timepoint", "Years_from_BL", "Data_Present",
    # Needed downstream
    "eTIV",
    # Per-hemisphere harmonized raws (12 panel halves).
    "Left-Hippocampus", "Right-Hippocampus",
    "Left-Amygdala", "Right-Amygdala",
    "Left-Caudate", "Right-Caudate",
    "Left-Putamen", "Right-Putamen",
    "Left-Thalamus", "Right-Thalamus",
    "Left-Accumbens-area", "Right-Accumbens-area",
    # Bilateral means rebuilt from harmonized halves.
    "Hippocampus_Mean", "Amygdala_Mean", "Caudate_Mean",
    "Putamen_Mean", "Thalamus_Mean", "Accumbens_Mean",
    # eTIV-normalized bilateral means. These are the Value column for
    # subcortical ROIs in cohort_long_harmonized.csv; computed as
    # harmonized mean / original (unharmonized) eTIV.
    "Hippocampus_Mean_norm", "Amygdala_Mean_norm", "Caudate_Mean_norm",
    "Putamen_Mean_norm", "Thalamus_Mean_norm", "Accumbens_Mean_norm",
]

# Cortical StructNames kept in the harmonized output. Out-of-panel
# regions are dropped so the file only holds harmonized values.
COHORT_CORTICAL_HARMONIZED_OUTPUT_STRUCTS: Tuple[str, ...] = (
    "entorhinal",
    "parahippocampal",
    "precuneus",
    "lingual",
    "caudalanteriorcingulate",
    "superiorfrontal",
)


# --- ROI panel (the 12 analysis ROIs) ---
#
# 6 subcortical (volume, eTIV-normalized) + 6 cortical (thickness, no
# ICV correction). Mainly based on ENIGMA-PD (Laansma et al. 2021).
#
# Subcortical entries have left_source/right_source/bilateral_col
# because the data sits in flat columns of cohort_wide.csv. Cortical
# entries only have fs_region because cortical data is in long format
# keyed by (Hemisphere, StructName).
#
# Fields:
#   roi          display name for tables/figures
#   tissue       'subcortical' or 'cortical'
#   measure      'volume' or 'thickness'
#   unit         'norm_volume' or 'mm'
#   cluster      figure ordering only
#   harmonize    whether the ROI goes through longCombat (always True)
#   icv_correct  divide by eTIV (subcortical only)
#   analysis_col column the analysis reads from
#   left_source, right_source, bilateral_col   subcortical only
#   fs_region    DKT atlas label, cortical only

ROI_PANEL_THESIS: List[Dict] = [
    # Subcortical (volume, eTIV-normalized)
    {
        "roi":           "Hippocampus",
        "tissue":        "subcortical",
        "measure":       "volume",
        "unit":          "norm_volume",
        "cluster":       "limbic_mtl",
        "left_source":   "Left-Hippocampus",
        "right_source":  "Right-Hippocampus",
        "bilateral_col": "Hippocampus_Mean",
        "analysis_col":  "Hippocampus_Mean_norm",
        "harmonize":     True,
        "icv_correct":   True,
    },
    {
        "roi":           "Amygdala",
        "tissue":        "subcortical",
        "measure":       "volume",
        "unit":          "norm_volume",
        "cluster":       "limbic_mtl",
        "left_source":   "Left-Amygdala",
        "right_source":  "Right-Amygdala",
        "bilateral_col": "Amygdala_Mean",
        "analysis_col":  "Amygdala_Mean_norm",
        "harmonize":     True,
        "icv_correct":   True,
    },
    {
        "roi":           "Caudate",
        "tissue":        "subcortical",
        "measure":       "volume",
        "unit":          "norm_volume",
        "cluster":       "basal_ganglia",
        "left_source":   "Left-Caudate",
        "right_source":  "Right-Caudate",
        "bilateral_col": "Caudate_Mean",
        "analysis_col":  "Caudate_Mean_norm",
        "harmonize":     True,
        "icv_correct":   True,
    },
    {
        "roi":           "Putamen",
        "tissue":        "subcortical",
        "measure":       "volume",
        "unit":          "norm_volume",
        "cluster":       "basal_ganglia",
        "left_source":   "Left-Putamen",
        "right_source":  "Right-Putamen",
        "bilateral_col": "Putamen_Mean",
        "analysis_col":  "Putamen_Mean_norm",
        "harmonize":     True,
        "icv_correct":   True,
    },
    {
        "roi":           "Thalamus",
        "tissue":        "subcortical",
        "measure":       "volume",
        "unit":          "norm_volume",
        "cluster":       "basal_ganglia",
        "left_source":   "Left-Thalamus",
        "right_source":  "Right-Thalamus",
        "bilateral_col": "Thalamus_Mean",
        "analysis_col":  "Thalamus_Mean_norm",
        "harmonize":     True,
        "icv_correct":   True,
    },
    {
        "roi":           "Accumbens",
        "tissue":        "subcortical",
        "measure":       "volume",
        "unit":          "norm_volume",
        "cluster":       "basal_ganglia",
        "left_source":   "Left-Accumbens-area",
        "right_source":  "Right-Accumbens-area",
        "bilateral_col": "Accumbens_Mean",
        "analysis_col":  "Accumbens_Mean_norm",
        "harmonize":     True,
        "icv_correct":   True,
    },

    # Cortical (thickness, no ICV correction)
    {
        "roi":           "Entorhinal cortex",
        "tissue":        "cortical",
        "measure":       "thickness",
        "unit":          "mm",
        "cluster":       "limbic_mtl",
        "fs_region":     "entorhinal",
        "analysis_col":  "Entorhinal_MeanThickness",
        "harmonize":     True,
        "icv_correct":   False,
    },
    {
        "roi":           "Parahippocampal gyrus",
        "tissue":        "cortical",
        "measure":       "thickness",
        "unit":          "mm",
        "cluster":       "limbic_mtl",
        "fs_region":     "parahippocampal",
        "analysis_col":  "Parahippocampal_MeanThickness",
        "harmonize":     True,
        "icv_correct":   False,
    },
    {
        "roi":           "Precuneus",
        "tissue":        "cortical",
        "measure":       "thickness",
        "unit":          "mm",
        "cluster":       "posterior_cortex",
        "fs_region":     "precuneus",
        "analysis_col":  "Precuneus_MeanThickness",
        "harmonize":     True,
        "icv_correct":   False,
    },
    {
        "roi":           "Lingual gyrus",
        "tissue":        "cortical",
        "measure":       "thickness",
        "unit":          "mm",
        "cluster":       "posterior_cortex",
        "fs_region":     "lingual",
        "analysis_col":  "Lingual_MeanThickness",
        "harmonize":     True,
        "icv_correct":   False,
    },
    {
        "roi":           "Caudal anterior cingulate",
        "tissue":        "cortical",
        "measure":       "thickness",
        "unit":          "mm",
        "cluster":       "frontal",
        "fs_region":     "caudalanteriorcingulate",
        "analysis_col":  "CaudalACC_MeanThickness",
        "harmonize":     True,
        "icv_correct":   False,
    },
    {
        "roi":           "Superior frontal",
        "tissue":        "cortical",
        "measure":       "thickness",
        "unit":          "mm",
        "cluster":       "frontal",
        "fs_region":     "superiorfrontal",
        "analysis_col":  "SuperiorFrontal_MeanThickness",
        "harmonize":     True,
        "icv_correct":   False,
    },
]

# Cluster order on the forest plot (top to bottom) and heatmap rows.
ROI_PANEL_CLUSTER_ORDER: Tuple[str, ...] = (
    "limbic_mtl",
    "basal_ganglia",
    "posterior_cortex",
    "frontal",
)


# --- Derived from ROI_PANEL_THESIS ---
# Flat lists for consumers that don't need the full panel metadata.

# Subcortical analysis columns. Used by 06 and 07.
ROI_PANEL_SUBCORTICAL_MEAN_NORM_COLS: List[str] = [
    e["analysis_col"] for e in ROI_PANEL_THESIS
    if e["tissue"] == "subcortical"
]

# Cortical StructNames in the panel. 06 uses these to scope the
# cortical longCombat run; 07 bilateralizes and analyzes them.
ROI_PANEL_CORTICAL_STRUCTS: List[str] = [
    e["fs_region"] for e in ROI_PANEL_THESIS
    if e["tissue"] == "cortical"
]


# Timepoint pairs for delta computation. BL-5Y is the primary window;
# BL-3Y is supplementary.
DELTA_WINDOWS: List[Tuple[str, str]] = [("BL", "3Y"), ("BL", "5Y")]


# --- Batch-size guardrails ---

# Warn when a batch has fewer than this many observations. ComBat's
# empirical-Bayes shrinkage gets unreliable below ~10.
MIN_BATCH_SIZE_WARN: int = 10

# Reject features whose smallest batch is below this. longCombat tends
# to fail or produce nonsense.
MIN_BATCH_SIZE_FAIL: int = 5


# --- Stage B column schemas ---
# Shape cohort_wide.csv (produced by 05_extract.py).

# Raw volume columns. Each gets a paired *_norm (raw / eTIV) column.
VOLUME_COLS_FOR_NORM: List[str] = [
    "BrainSegVol", "BrainSegVolNotVent",
    "CortexVol", "lhCortexVol", "rhCortexVol",
    "CerebralWhiteMatterVol",
    "lhCerebralWhiteMatterVol", "rhCerebralWhiteMatterVol",
    "SubCortGrayVol", "TotalGrayVol",
    "SupraTentorialVol", "SupraTentorialVolNotVent",
    "VentricleChoroidVol",
    "Left-Thalamus",    "Right-Thalamus",    "Thalamus_Mean",
    "Left-Caudate",     "Right-Caudate",     "Caudate_Mean",
    "Left-Putamen",     "Right-Putamen",     "Putamen_Mean",
    "Left-Pallidum",    "Right-Pallidum",    "Pallidum_Mean",
    "Left-Hippocampus", "Right-Hippocampus", "Hippocampus_Mean",
    "Left-Amygdala",    "Right-Amygdala",    "Amygdala_Mean",
    "Left-Accumbens-area", "Right-Accumbens-area", "Accumbens_Mean",
    "Left-VentralDC",   "Right-VentralDC",
    "Brain-Stem", "WM-hypointensities",
    "lh_RegionalGrayVol_Sum", "rh_RegionalGrayVol_Sum",
    "CortexRegionalGrayVol_Sum",
    "Left-Lateral-Ventricle", "Right-Lateral-Ventricle",
    "Left-Inf-Lat-Vent",      "Right-Inf-Lat-Vent",
    "3rd-Ventricle", "4th-Ventricle",
    "LateralVentricles_Total", "InfLatVentricles_Total",
    "Ventricles_Total_Main",
    # Cerebellar volumes from CerebNet. _Total here means component
    # sum, not bilateral mean.
    "Cbm_Left_I_IV", "Cbm_Right_I_IV",
    "Cbm_Left_V",    "Cbm_Right_V",
    "Cbm_Anterior_Total",
    "Cbm_Left_CortexVol", "Cbm_Right_CortexVol", "Cbm_CortexVol_Total",
    "Cbm_Left_WhiteMatter", "Cbm_Right_WhiteMatter",
    "Cbm_WhiteMatter_Total",
    "Cbm_Vermis_VI", "Cbm_Vermis_VII", "Cbm_Vermis_VIII",
    "Cbm_Vermis_IX", "Cbm_Vermis_X", "Cbm_Vermis_Total",
]

# Column order in cohort_wide.csv. Columns produced by the extractor
# but not listed here are appended at the end in insertion order.
ORDERED_METRIC_COLS: List[str] = [
    # Global
    "MaskVol", "eTIV",
    "BrainSegVol", "BrainSegVolNotVent", "VentricleChoroidVol",
    "lhCortexVol", "rhCortexVol", "CortexVol", "AI_CortexVol_pct",
    "lhCerebralWhiteMatterVol", "rhCerebralWhiteMatterVol",
    "CerebralWhiteMatterVol",
    "SubCortGrayVol", "TotalGrayVol",
    "SupraTentorialVol", "SupraTentorialVolNotVent",
    "BrainSegVol_to_eTIV", "MaskVol_to_eTIV",
    # Subcortical ROIs
    "Left-Thalamus",    "Right-Thalamus",    "Thalamus_Mean",
    "Left-Caudate",     "Right-Caudate",     "Caudate_Mean",
    "Left-Putamen",     "Right-Putamen",     "Putamen_Mean",
    "Left-Pallidum",    "Right-Pallidum",    "Pallidum_Mean",
    "Left-Hippocampus", "Right-Hippocampus", "Hippocampus_Mean",
    "Left-Amygdala",    "Right-Amygdala",    "Amygdala_Mean",
    "Left-Accumbens-area", "Right-Accumbens-area", "Accumbens_Mean",
    "Left-VentralDC",   "Right-VentralDC",
    "Brain-Stem", "WM-hypointensities",
    # Ventricles (component sums; _Total kept)
    "Left-Lateral-Ventricle", "Right-Lateral-Ventricle",
    "Left-Inf-Lat-Vent",      "Right-Inf-Lat-Vent",
    "3rd-Ventricle", "4th-Ventricle",
    "LateralVentricles_Total", "InfLatVentricles_Total",
    "Ventricles_Total_Main",
    # Cerebellar volumes (CerebNet). _Total = component sum.
    "Cbm_Left_I_IV", "Cbm_Right_I_IV",
    "Cbm_Left_V",    "Cbm_Right_V",
    "Cbm_Anterior_Total",
    "Cbm_Left_CortexVol", "Cbm_Right_CortexVol", "Cbm_CortexVol_Total",
    "Cbm_Left_WhiteMatter", "Cbm_Right_WhiteMatter",
    "Cbm_WhiteMatter_Total",
    "Cbm_Vermis_VI", "Cbm_Vermis_VII", "Cbm_Vermis_VIII",
    "Cbm_Vermis_IX", "Cbm_Vermis_X", "Cbm_Vermis_Total",
    # Asymmetry indices
    "AI_Thalamus_pct", "AI_Caudate_pct", "AI_Putamen_pct",
    "AI_Pallidum_pct", "AI_Hippocampus_pct", "AI_Amygdala_pct",
    # Cortical surface summary
    "lh_NumVert", "rh_NumVert",
    "lh_WhiteSurfArea", "rh_WhiteSurfArea", "Total_WhiteSurfArea",
    "lh_MeanThickness", "rh_MeanThickness", "MeanThickness_Weighted",
    "AI_WhiteSurfArea_pct", "AI_MeanThickness_pct",
    # Cortical regional consistency
    "lh_RegionalGrayVol_Sum", "rh_RegionalGrayVol_Sum",
    "CortexRegionalGrayVol_Sum",
    "AI_RegionalGrayVol_pct",
    "lh_RegionalGrayVol_vs_lhCortexVol_pctdiff",
    "rh_RegionalGrayVol_vs_rhCortexVol_pctdiff",
    "CortexRegionalGrayVol_vs_CortexVol_pctdiff",
]
