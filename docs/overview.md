# MESO-FC-Controls — Repository Overview & Data Flow

This document explains how the repository is organised and how data flows through
the pipeline, then gives a concise, file-by-file reference for every script.

The pipeline replicates Schmidt et al. (2025), which studies **mesoscale
resting-state functional connectivity (rs-FC)** in the human primary visual cortex
(V1) from high-resolution 7T fMRI. The core question is how rs-FC changes with
**cortical distance**, **ocular polarity** (same- vs. different-eye preference of the connected vertex pair), **cortical depth** (layers), and **V1 subregion**, and whether the rs-FC spatial pattern mirrors the **ocular dominance column (ODC)** maps.

---

## 1. Pipeline stages

The repository is organised into five stages, executed in order:

| # | Folder | Language | Purpose |
|---|--------|----------|---------|
| 1 | `preprocessing/` | MATLAB (+ FreeSurfer) | Turn raw fMRI into clean, surface-sampled, layer-resolved time series. |
| 2 | `processing/` | MATLAB | Compute rs-FC/selectivity per subject, layer, ROI, and ocular-polarity group; save `.mat` results. |
| 3 | `stats/` | MATLAB | Statistical modeling (repeated-measures ANOVA, mixed-effects models) of the saved results. |
| 4 | `plotting/` | Python | Load the `.mat` results and produce publication figures. |
| 5 | `control_analyses/` | MATLAB + Python | Supplementary/control analyses and figures. |

**Script-naming convention.** Within each stage, scripts are prefixed to indicate the analysis they belong to:

- `A` — Analysis A: effect of **distance** (and ocular polarity) on rs-FC → Fig 2.
- `B` — Analysis B: effect of **ocular-preference strength (beta)**, **layer**, and
  **V1 subregion** (intra- and inter-hemispheric) on rs-FC / selectivity → Figs 3–7.
- `C` — Analysis C: **robustness/replicability** of the rs-FC measure across scanning
  sessions → Fig S3.
- `Z` — Analysis Z: **correlation of rs-FC patterns with the ODC differential map** →
  Fig 8.
- `X` — Control / supplementary analyses → Figs 1, S1–S4 and various sanity checks.
- `A001a` / `B001a` / `B001b` / `Z001a` / `Z001b` — the per-subject *processing
  functions* called by the corresponding *Miner* scripts (`A01`, `B01`, `Z01`).

---

## 2. Data flow

### 2.1 Inputs (per subject)

The analysis relies on several per-subject datasets:

- **Raw resting-state fMRI runs** (magnitude data, `.nii`), one or more runs.
- **Upsampled anatomical surfaces** (e.g. `*_anat_upsample_B1justsub`) with the
  inflated and flattened occipital patch surfaces.
- **Retinotopy labels** defining V1 and its subregions (dorsal/ventral,
  posterior/central, anterior/peripheral), e.g. `V1_Upsampled_Rtopy.label`.
- **ODC differential beta maps** (`beta.nii`) from monocular-stimulation analyses —
  the signed value encodes ocular preference (negative ≈ left eye, positive ≈ right
  eye) and its magnitude encodes preference strength.

### 2.2 Stage 1 — Preprocessing (`preprocessing/`)

1. **NORDIC denoising** (`A01`) removes thermal noise from the magnitude time series,
   producing `data_denoised10` volumes for each run.
2. **FSFAST + layer preprocessing** (`A02`) runs the FreeSurfer FSFAST stream:
   - slice-timing correction (`stc-sess`),
   - upsampling (`mri_convert --upsample 2`),
   - motion/registration via `bbregister`,
   - `preproc-sess` for the remainder of the preprocessing (detrending, motion
     correction),
   - **layer sampling**: `mri_vol2surf` samples the signal to 9 intermediate cortical
     surfaces (layers `00`–`10`),
   - **intracortical smoothing**: `mris_smooth_intracortical` radially smooths into
     layer-averaged targets (e.g. rad `0-2`, `4-6`, `8-10`, `0-10`),
   - **covariate generation**: `fcseed-sess` produces the white-matter (`wm.dat`) and
     motion (`mcprextreg`) regressors.

The result for each run is a set of surface overlay files named like
`<hemi>.fmcpr.sm0.self.midgray.00.nb1_rad<layers>.mgz` plus the covariate files.

### 2.3 Stage 2 — Processing (`processing/`)

A **Miner** script (`A01_FC_Miner`, `B01_Selectivity_Miner`, `Z01_FC_ODC_corr_Miner`) defines the analysis parameters (ROI, layers, subjects, detrending/high-pass flags, quantile counts) and loops over layers and subjects, calling a per-subject function.

Each per-subject function typically does the following:

1. **Load ODC maps + surfaces/labels** and crop the occipital patch to the ROI
   (V1 or a subregion), building binary vertex masks.
2. **Split V1 vertices by ocular polarity** into negative (eye 1) and positive
   (eye 2) groups, binned into `nQuant = 10` beta quantiles each.
3. **Load resting-state data** for the ROI at the requested layer target, then apply
   optional **detrending** (`detrend`) and **high-pass filtering** (`highpass`).
4. **Compute partial correlations** (`partialcorr`) between vertex pairs, regressing
   out white-matter (`wm.dat`) and motion (`mcprextreg`) covariates.
   - Intra-hemispheric analyses correlate vertices within each hemisphere.
   - Inter-hemispheric analyses (`B001b`) correlate the left-hemisphere ROI with the right-hemisphere ROI.
5. **Quantile grouping & subsampling**: distances between vertex pairs are grouped
   into 10 distance quantiles. To compare *alike* vs *unalike* pairs fairly,
   vertex-pair **subsampling** matches the distance distribution of the two groups.
6. **Save results** — the mean rs-FC per distance/beta/type/layer is written to
   `CorrelationMtx_*.mat` per subject, together with an `AnalysisParam` struct
   (parameters, mask sizes, beta medians, distance medians, and the selected-pair
   indices used for subsampling).

The stored `Data_Combined` array dimensions encode:
`[distance-quantile, type, hemisphere, beta-quantile-1, beta-quantile-2, z-or-r]`.

### 2.4 Stage 3 — Statistics (`stats/`)

Statistics scripts load the per-subject `CorrelationMtx_*.mat` files, assemble
subject × condition matrices, and run:

- **repeated-measures ANOVA** with the Statistics Toolbox (`fitrm` + `ranova`) to test main effects and interactions of distance, type (alike/unalike), beta, layer, ROI, and session;
- **linear mixed-effects models** (`fitlme`) with a random subject intercept to account for the hierarchical (subject × vertex-pair) structure.

Results are displayed and optionally written to `.xlsx` tables.

### 2.5 Stage 4 — Plotting (`plotting/`)

Python scripts read the saved `.mat` results and generate the manuscript figures
(violin/box plots, beta×beta correlation matrices per layer, rs-FC↔ODC radius curves, etc.). Figures are saved as high-resolution `.tiff` files.

### 2.6 Stage 5 — Control analyses (`control_analyses/`)

Supplementary analyses validating robustness and sanity-checking assumptions: fMRI
timeseries and percent-signal-change plots, cortical-depth statistics, subsampling-iteration comparison, vertex-number and beta-strength comparisons, ODC distance distributions, and an example rs-FC-map renderer.

---

## 3. File-by-file reference (concise)

### `preprocessing/`

| File | What it does |
|------|--------------|
| `A01_preproc_data_NORDIC.m` | Runs NORDIC thermal-noise denoising on raw fMRI magnitude data for each subject/run. |
| `A02_preproc_data_FSFAST.m` | Full FSFAST pipeline: denoising, slice-timing correction, upsampling, registration, `preproc-sess`, layer (`mri_vol2surf`) sampling, intracortical smoothing, wm/motion covariate generation. |

### `processing/`

| File | What it does |
|------|--------------|
| `A01_FC_Miner.m` | Analysis A driver: loops layers/subjects, calls `A001a` to compute rs-FC vs distance & ocular polarity (Fig 2). |
| `A001a_FC_Proc_Data_subsample_beta.m` | Per-subject: loads ODC + rs data, detrends/high-pass filters, partial-correlation, beta/distance quantiles, subsampling; saves `CorrelationMtx_FC*` + subsample indices. |
| `B01_Selectivity_Miner.m` | Analysis B driver: loops layers/subjects; routes intra- vs inter-hemispheric to `B001a`/`B001b` (Figs 3–7). |
| `B001a_Selectivity_Proc_Data.m` | Per-subject intra-hemispheric selectivity (alike/eye1/eye2/unalike rs-FC) per layer, ROI, beta quantiles; reuses A001a subsample indices; saves `CorrelationMtx_Selectivity*`. |
| `B001b_Selectivity_Proc_Data_Interhemi.m` | Per-subject inter-hemispheric selectivity: partial-correlation between left and right hemisphere ROIs; saves `CorrelationMtx_Selectivity.mat`. |
| `Z01_FC_ODC_corr_Miner.m` | Analysis Z driver: loops layers/subjects; routes to `Z001a` or `Z001b` depending on chosen null hypothesis (Fig 8). |
| `Z001a_FC_ODC_corr_Proc_Data.m` | Per-subject: correlates 2D rs-FC ring maps (1000 random seeds, several radii) with a spatially shifted ODC map; saves `meanAbsR`. |
| `Z001b_FC_ODC_corr_Proc_Data_rotated.m` | Per-subject: like `Z001a` but the null uses a 180°-rotated ODC map; saves `meanAbsR_..._rotated`. |

### `stats/`

| File | What it does |
|------|--------------|
| `A02_FC_Stats.m` | rm-ANOVA for distance × type (alike/unalike) effects on rs-FC, V1 layer 0-10 (Fig 2). |
| `B02a_Selectivity_Stats.m` | rm-ANOVA + LME for beta × layer × type on rs-FC, V1 intra-hemispheric (Fig 3). |
| `B02b_Selectivity_Stats_DorsVSVent.m` | rm-ANOVA comparing V1_Dorsal vs V1_Ventral (Fig 4). |
| `B02c_Selectivity_Stats_CentVSPeri.m` | rm-ANOVA comparing V1_Anterior vs V1_Posterior (central/peripheral) (Fig 5). |
| `B02d_Selectivity_Stats_interhemi_DorsVSVent.m` | Inter-hemispheric rm-ANOVA, V1_Dorsal vs V1_Ventral (Fig 6). |
| `B02e_Selectivity_Stats_interhemi_CentVSPeri.m` | Inter-hemispheric rm-ANOVA, V1_Anterior vs V1_Posterior (Fig 7). |
| `C01_FC_Replicability_Stats.m` | rm-ANOVA + between-session correlations testing rs-FC replicability (Fig S3). |
| `Z02_FC_ODC_corr_Stats.m` | rm-ANOVA on rs-FC↔ODC correlation vs radius, layer, and null type (shifted/rotated) (Fig 8). |

### `plotting/`

| File | What it does |
|------|--------------|
| `A03a_FC_Plots_boxplot.py` | Split-violin plots of rs-FC vs distance for alike/unalike + selectivity, V1 (Fig 2). |
| `A03b_FC_Plots_boxplot_replicability1.py` | Violin plot of rs-FC vs distance for the old control group (left Fig S3). |
| `A03c_FC_Plots_boxplot_replicability2.py` | Violin plot of rs-FC vs distance for the new/control-2 group (right Fig S3). |
| `B03a_Selectivity_Plots.py` | Beta×beta rs-FC correlation matrices for alike/unalike/selectivity per layer (Fig 3). |
| `B03b_Selectivity_Plots_Subregions.py` | rs-FC and selectivity matrices per layer for a V1 subregion, intra-hemispheric (Figs 4–5). |
| `B04c_Selectivity_Plots_interhemi.py` | rs-FC and selectivity matrices per layer for a V1 subregion, inter-hemispheric (Figs 6–7). |
| `Z02_FC_ODC_corr_Plot.py` | Plots rs-FC↔ODC mean-|r| vs ring radius for H1, shifted-H0, rotated-H0 (+p-values) (Fig 8). |

### `control_analyses/`

| File | What it does |
|------|--------------|
| `X01_FMRI_Timeseries_Plots_psc_2.py` | Plots vertex/mean resting-state time series as % signal change across layers (Fig 1). |
| `X02_FMRI_Depth_Stats_psc.m` | rm-ANOVAs + bar plots of %-signal-change across cortical depth for V1 and subregions. |
| `X03_ocular_polarity_dist_comparison_Stats.m` | Tests whether alike/unalike vertex-pair distance distributions match after subsampling. |
| `X04_vertex_number_comparison_Stats.m` | Compares number of vertices assigned to dominant vs non-dominant eye (Fig S3). |
| `X05_beta_comparison_Stats.m` | Compares ocular-preference strength (beta) between dominant and non-dominant eye vertices. |
| `X06_distribution_distance_quantile_medians.py` | Plots median beta within each beta quantile across subjects. |
| `X07_Selectivity_Stats_subsample_iter_comp.m` | Tests whether results are stable across random subsampling iterations. |
| `X08_beta_quantile_medians.m` | Plots median beta for each beta quantile from saved params. |
| `X09_rsFC_map.m` | Renders an example 2D rs-FC map for a seed vertex, saved as `.mgz`. |
| `X10_distribution_distance_quantile_medians.py` | Plots median distance within each distance quantile across subjects (Fig S1). |

---

## 4. Key concepts glossary

- **rs-FC** — resting-state functional connectivity: the (partial) correlation between
  vertex time series.
- **Ocular polarity** — the sign of the ODC beta at a vertex (≈ left/right eye
  preference). Vertex pairs are **alike** (both prefer the same eye) or **unalike**
  (prefer different eyes).
- **Selectivity** — the difference between alike and unalike rs-FC.
- **Beta quantile** — a binning of vertices by ocular-preference *strength*
  (`nQuant = 10` bins).
- **Distance quantile** — a binning of vertex-pair Euclidean distances (on the
  flattened patch, `distThresh = 3 mm` minimum).
- **Subsampling** — randomly dropping vertex pairs so alike/unalike groups have
  matched distance distributions before comparison.
- **Intra- vs inter-hemispheric** — correlations within one hemisphere vs across the left/right hemispheres of the same ROI.
