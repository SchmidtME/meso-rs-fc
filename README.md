# Mesoscale functional connectivity of human V1 — analysis pipeline

This repository contains the complete, end-to-end code used to replicate the study
presented in

**Schmidt et al. (2025).** *Unraveling the mesoscale functional connectivity of the
human primary visual cortex using high-resolution functional MRI.* bioRxiv.
DOI: <https://doi.org/10.1101/2025.03.27.645795>

The pipeline analyses high-resolution 7T resting-state fMRI data to investigate how
**resting-state functional connectivity (rs-FC)** within and between the primary
visual cortices (V1) depends on:

- **Cortical distance** between vertices,
- **Ocular polarity** — whether connected vertex pairs prefer the same eye
  (*alike*) or different eyes (*unalike*),
- **Cortical depth** (layers, e.g. superficial / middle / deep),
- **V1 subregions** (dorsal/ventral, central/peripheral),
- and whether rs-FC spatial patterns relate to the **ocular dominance column (ODC)**
  maps derived from population receptive-field / monocular stimulation data.

The repository is organised into five sequential stages. Code is a mix of MATLAB
(preprocessing, processing, statistics) and Python (plotting + some control analyses).

## Repository layout

| Folder               | Stage                                                  | Manuscript figures |
|----------------------|--------------------------------------------------------|--------------------|
| `preprocessing/`     | NORDIC denoising + FSFAST surface preprocessing       | —                  |
| `processing/`        | Compute rs-FC / selectivity and save per-subject `.mat`| 2–8                |
| `stats/`             | Repeated-measures ANOVA + mixed-effects models         | accompanying stats |
| `plotting/`          | Publication-quality figures from saved `.mat` results  | 2–8, S3            |
| `control_analyses/`  | Supplementary / control analyses & figures             | 1, S1–S4           |

See [`docs/overview.md`](docs/overview.md) for a detailed description of the data
flow and a concise, file-by-file reference.

## Quick start

The scripts reference **absolute cluster paths** (e.g. `/space/ardebil_001/...` and
`/autofs/space/...`) matching the environment where the analysis was originally run.
These paths do **not** exist outside that cluster. Before running anything you must
replace them with the locations of your own data (see *Configuring paths* below).

The recommended execution order follows the data flow:

```text
preprocessing  ->  processing  ->  stats  ->  plotting
                      (control_analyses explored as needed)
```

### 0. Configuring paths

Every script hard-codes its input/output directories at the top. To run the pipeline
you should:

1. Copy the relevant data (raw fMRI runs, upsampled anatomical surfaces, retinotopy /
   ODC beta maps) into your own directory structure, or
2. Edit the directory variables at the top of each script (`dir_orig_data`,
   `dir_data`, `Root`, `rsFolder`, `anatFolder`, `odcFolder`, `labelFolder`, …) to
   point at your data, and
3. Remove/replace any missing `addpath(...)` lines that point to cluster-specific
   software (e.g. the NORDIC and FreeSurfer paths in `preprocessing/`).

A clean approach for a sandbox is to centralise these paths as environment variables
or a shared config snippet and reference them from each script.

### 1. Preprocessing (MATLAB + FreeSurfer)

Denoise and bring data onto cortical surfaces before connectivity analysis.

```text
preprocessing/A01_preproc_data_NORDIC.m     % NORDIC denoising of raw fMRI magnitude data
preprocessing/A02_preproc_data_FSFAST.m     % slice-timing, upsample, registration,
                                            % layer sampling, smoothing, covariates
```

> **Note:** `A02` flag-gates its individual steps with `denoise`, `stc`, `upsample`,
> `do_preproc_*`, `gen_layers`, `vert_smooth`, `gen_covariates`. Set these to `1` as
> needed (they are disabled by default).

### 2. Processing (MATLAB)

Compute rs-FC / selectivity for each subject, layer, ROI, and ocular-polarity
combination and save `CorrelationMtx_*.mat` result files.

```text
processing/A01_FC_Miner.m                  % Analysis A  -> Fig 2
processing/B01_Selectivity_Miner.m         % Analysis B  -> Figs 3-7
processing/Z01_FC_ODC_corr_Miner.m         % Analysis Z  -> Fig 8
```

Each *Miner* script loops over layers/subjects and calls a per-subject function
(`A001a_FC_Proc_Data_subsample_beta.m`, `B001a/B001b_Selectivity_Proc_Data*.m`,
`Z001a/Z001b_FC_ODC_corr_Proc_Data*.m`) that loads the ODC maps and resting-state
data, applies detrending/high-pass filtering, computes partial correlations, and
groups vertex pairs into distance and ocular-preference (beta) quantiles.

### 3. Statistics (MATLAB)

Run repeated-measures ANOVAs (`fitrm`/`ranova`) and mixed-effects models (`fitlme`)
on the saved `.mat` results.

```text
stats/A02_FC_Stats.m                                  % Fig 2        (distance x type)
stats/B02a_Selectivity_Stats.m                        % Fig 3        (beta x layer x type)
stats/B02b/B02c_Selectivity_Stats_*                  % Figs 4-5     (V1 subregions, intra)
stats/B02d/B02e_Selectivity_Stats_interhemi_*        % Figs 6-7     (V1 subregions, inter)
stats/C01_FC_Replicability_Stats.m                   % Fig S3       (session replicability)
stats/Z02_FC_ODC_corr_Stats.m                        % Fig 8        (FC <-> ODC correlation)
```

### 4. Plotting (Python)

Load the saved `.mat` results and render the publication figures.

```text
plotting/A03a_FC_Plots_boxplot.py             % Fig 2
plotting/A03b/A03c_FC_Plots_boxplot_replicability*.py   % Fig S3
plotting/B03a_Selectivity_Plots.py            % Fig 3
plotting/B03b_Selectivity_Plots_Subregions.py % Figs 4-5
plotting/B04c_Selectivity_Plots_interhemi.py  % Figs 6-7
plotting/Z02_FC_ODC_corr_Plot.py              % Fig 8
```

### 5. Control / supplementary analyses

Exploratory and supplementary analyses (timeseries, depth statistics, subsampling
robustness, vertex/beta comparisons, rs-FC map rendering). See
[`docs/overview.md`](docs/overview.md) for per-file details.

## Dependencies

The pipeline uses several external software packages and libraries:

- **FreeSurfer (v7.4.0)** with the **FSFAST** analysis stream (`preproc-sess`,
  `bbregister`, `mri_vol2surf`, `mris_smooth_intracortical`, `fcseed-sess`, …).
- **MATLAB (R2023b)** with the Statistics and Machine Learning Toolbox
  (`fitrm`, `ranova`, `fitlme`).
- **Python (v3.12)** with `numpy`, `scipy`, `matplotlib`, `seaborn`, `pandas`,
  `nibabel`, and `h5py`.
- **NORDIC_Raw** (Moeller et al., 2021) for thermal-noise removal of the fMRI data.
- In-house I/O helpers for FreeSurfer surface/NIfTI data (`read_surf`, `read_patch`,
  `read_ROIlabel`, `load_nifti`, `load_mgh`, `highpass`, …) located in a separate
  utility package (referenced via `addpath` / `sys.path.append` in the scripts).

## License

This repository is distributed under the **GNU General Public License v3.0**. See the
[`LICENSE`](LICENSE) file for details.

## Attribution

**Authors / Maintainers:**

- Marianna Elisa Schmidt — marianna.schmidt@maxplanckschools.de
- Iman Aganj — iaganj@mgh.harvard.edu
- Shahin Nasr — shahin.nasr@mgh.harvard.edu
