#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script plots fMRI time series (percent signal change) of V1 resting-state
data for a single subject, broken down by cortical depth level (deep, middle,
superficial). It loads the resting-state data masked to the V1 patch, optionally
detrends and high-pass filters the signal, and produces two figures: the time
series of one selected vertex and the mean time series across all vertices in V1.
These correspond to panels E and D of Figure 1 of the manuscript. If requested,
the underlying time series are also written to CSV files for further plotting.

@author: Marianna Elisa Schmidt (marianna.schmidt@maxplanckschools.de)
"""

from datetime import datetime
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.signal import detrend
import sys
sys.path.append('/autofs/space/ardebil_001/users/Others/Marianna/Code/MesoVision/Utils')
from load_data import read_surf, read_patch, read_ROIlabel, load_nifti, highpass

#%% Specification of parameters

save_figures = 0
if save_figures:
    date_stamp = datetime.now().strftime('%Y-%m-%d')
    save_dir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Figure_1/{date_stamp}'
    os.makedirs(save_dir, exist_ok=True)

# for the manuscript we used rh of ylri

# Define parameters
subName = 'aroo'

AnalysisParam = {
    'AlsoCombineSessions': True,
    'withCovariates': True,
    'StrtPnt': 0,
    'TmSeries_Length': 128,
    'nPC': 0,
    'prc': [0, 0],
    'nQuant': 5,
    'distThresh': 3,
    'stepLine': 0.5,
    'denoised': True,
    'stc': True,
    'detrending': True,
    'hpf': True,
    'pca_denoise': False,
    'pca_comps': 0,
    'pca_Thresh': [0],
    'ROI': 'V1',
    'Type': 'intrahemispheric',
    'TrgFiles': ['.fmcpr.sm0.self.midgray.00.nb1_rad0-2',
                 '.fmcpr.sm0.self.midgray.00.nb1_rad4-6',
                 '.fmcpr.sm0.self.midgray.00.nb1_rad8-10']
}
hemi = 'rh'

# Define directories
anatFolder = os.path.join('/autofs/space/ardebil_001/users/Shared/good_subjects_anat', f'{subName[:4]}_anat_upsample_B1justsub')
odcFolder = os.path.join('/autofs/space/ardebil_002/users/Shahin/Stereopsis_Upsampled/Subjects_Control_Upsampled/', f'{subName[:4]}_Final_Upsampled_VOL')

#%% Load anatomical surfaces, ODC map and ROI labels

if AnalysisParam['denoised']:
    rsFolder = os.path.join('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Data_Denoised/', f'{subName}', 'bold_Close_Upsampled2')
else:
    rsFolder = os.path.join('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Data/', f'{subName}', 'bold_Close_Upsampled2')

# Load and process data
print(f'Loading {hemi} Data ...')

# Load ODC map
dataFile_odc = os.path.join(odcFolder, f'Stereopsis_TR3_Columnar_Smoothing_0-2.{hemi}', 'R_C', f'{hemi}.ffx.osgm.wls', 'beta.nii')
sig_odc = load_nifti(dataFile_odc)

# Load anatomical surfaces and patch
vtx, fac = read_surf(os.path.join(anatFolder, 'surf', f'{hemi}.inflated'))
if fac is not None:
    fac += 1
occ_patch = read_patch(os.path.join(anatFolder, 'surf', f'{hemi}.occip_patch.flat'))

# Load ROI labels and create V1 patch
roi_label = read_ROIlabel(os.path.join(anatFolder, 'label', 'High_Res_Upsampled', f'{hemi}.V1_Upsampled_Adjusted.label'))

# Initialize an empty V1 patch container
v1_patch = {
    'ind': [],
    'x': [],
    'y': [],
    'z': [],
    'vno': []
}

# Iterate through occ_patch['ind'] and keep only vertices that fall inside the V1 label
for i, ind_value in enumerate(occ_patch['ind']):
    if ind_value in roi_label:
        v1_patch['ind'].append(ind_value)
        v1_patch['x'].append(occ_patch['x'][i])
        v1_patch['y'].append(occ_patch['y'][i])
        v1_patch['z'].append(occ_patch['z'][i])
        v1_patch['vno'].append(occ_patch['vno'][i])

# Convert lists to numpy arrays for consistency
for key in v1_patch:
    v1_patch[key] = np.array(v1_patch[key])
v1_patch['npts'] = len(v1_patch['ind'])

# Create binary masks marking which surface vertices belong to the V1 patch
bin_mask_vtx_v1_patch = np.zeros(vtx.shape[0], dtype=bool)
bin_mask_vtx_v1_patch[v1_patch['ind']] = True

#%% Load and process the resting-state time series

# Collect the session subfolders (resting-state runs) for this subject
rsSessions = [d.name for d in os.scandir(rsFolder) if d.is_dir() and d.name != 'masks']
n_sessions = len(rsSessions)

# Preallocate arrays for all resting-state data
# (sessions x depth-level target files x V1 vertices x timepoints)
all_rs_v1_patch_4V1 = np.empty((n_sessions, len(AnalysisParam['TrgFiles']), v1_patch['npts'], AnalysisParam['TmSeries_Length']))
all_rs_v1_patch_4vtx = np.empty((n_sessions, len(AnalysisParam['TrgFiles']), v1_patch['npts'], AnalysisParam['TmSeries_Length']))

for sessionNum, rsSessionName in enumerate(rsSessions, start=1):
    print(f'Loading Resting-State Data Run {sessionNum}!')

    for trgIdx, trgFile in enumerate(AnalysisParam['TrgFiles']):
        dataFile_rs = os.path.join(rsFolder, rsSessionName, f'{hemi}{trgFile}.mgz')
        rs = nib.load(dataFile_rs).get_fdata().squeeze()
        # Restrict to the requested time window and to V1-patch vertices
        rs = rs[:, AnalysisParam['StrtPnt']:AnalysisParam['StrtPnt'] + AnalysisParam['TmSeries_Length']]
        rs = rs[bin_mask_vtx_v1_patch, :]
        # Calculate percent signal change relative to the per-vertex mean signal
        mean_signal = np.mean(rs,axis=1, keepdims=True)
        zero_mask = (mean_signal == 0)
        mean_signal[zero_mask] = np.nan
        psc = ((rs - mean_signal) / mean_signal) * 100

        # Optional temporal linear detrending of the percent signal change
        if AnalysisParam['detrending']:
            rs_out = detrend(psc.T, type='linear').T

        # Optional temporal high-pass filtering (fs = 1/4 s)
        if AnalysisParam['hpf']:
            rs_out = highpass(rs_out.T, 0.01, 1/4).T

        # Store processed (detrended/filtered) and raw psc data
        all_rs_v1_patch_4vtx[sessionNum-1, trgIdx, :, :] = rs_out
        all_rs_v1_patch_4V1[sessionNum-1, trgIdx, :, :] = psc

#%% Select a single vertex for the per-vertex time-series plot

all_rs_v1_patch_4vtx2 = all_rs_v1_patch_4vtx

#%%

all_rs_v1_patch_4vtx = all_rs_v1_patch_4vtx2

#%%

#vertex_idx = np.where(v1_patch['ind'] == 7998)[0][0] # randomly picked in ODC map

vertex_idx = np.where(v1_patch['ind'] == 8014)[0][0] # randomly picked in ODC map
#vertex_idx = np.where(v1_patch['ind'] == 129674)[0][0] # randomly picked in ODC map
# Keep only the single vertex of interest (sessions x target files x timepoints)
all_rs_v1_patch_4vtx = all_rs_v1_patch_4vtx[:, :, vertex_idx, :]  

#%% Plot time series of a single vertex

# Generate x values (0 to 117)
x = np.arange(all_rs_v1_patch_4vtx.shape[-1])

# Define colors and labels
colors = ['#0072B2', '#E69F00', '#009E73']
labels = ['deep', 'middle', 'superficial']

# Create the original plot for the entire time series, averaged across sessions for the vertex
plt.figure(figsize=(15, 8), facecolor='white')
ax = plt.gca()
ax.set_facecolor('white')

# Plot the mean across vertices for a single vertex
for trgIdx, trgFile in enumerate(AnalysisParam['TrgFiles']):
    plt.plot(x, np.mean(all_rs_v1_patch_4vtx,axis=0)[trgIdx], marker='o', linestyle='-', color=colors[trgIdx], label=labels[trgIdx])

# Add a red transparent rectangle
plt.axvspan(0, 9, color='red', alpha=0.3)

# Adjusting the axis limits
plt.xlim(0, 131)

# Set specific x-ticks and y-ticks
plt.xticks(ticks=[0, 10, 128], labels=['0', '10', '128'])
plt.yticks(ticks=[-2.1, 0, 2.1], labels=['-2.1', '0', '2.1'])

# Labels and Title
plt.xlabel("Time (s)", fontsize=14)
plt.ylabel("Signal Change (%)", fontsize=14)
plt.title("Time Series of Vertex Signal Change", fontsize=16)

# Add legend
plt.legend()

plt.grid(False)

# Save figure
if save_figures:
    file_name = f"Figure_1E.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=600)

# Show the plot
plt.show()

# --- Zoomed-in Plot ---
#%% Plot mean time series across all vertices
plt.figure(figsize=(15, 8), facecolor='white')
ax = plt.gca()
ax.set_facecolor('white')

# Create a similar plot across three target files (mean across vertices)
plt.figure(figsize=(15, 8), facecolor='white')
ax = plt.gca()
ax.set_facecolor('white')

# Calculate mean time series for all vertices
mean_vertex_time_series = np.mean(np.mean(all_rs_v1_patch_4V1, axis=0), axis=1)

# Plot the mean time series across target files
for trgIdx, trgFile in enumerate(AnalysisParam['TrgFiles']):
    plt.plot(x, mean_vertex_time_series[trgIdx], marker='o', linestyle='-', color=colors[trgIdx], label=labels[trgIdx])

# Add a red transparent rectangle
plt.axvspan(0, 9, color='red', alpha=0.3)

# Adjusting the axis limits
plt.xlim(0, 131)

# Set specific x-ticks and y-ticks
plt.xticks(ticks=[0, 10, 128], labels=['0', '10', '128'])
plt.yticks(ticks=[-2.5, 0, 2.5], labels=['-2.5', '0', '2.5'])

# Labels and Title
plt.xlabel("Time (s)", fontsize=14)
plt.ylabel("Signal Change (%)", fontsize=14)
plt.title("Mean Time Series Across Vertices", fontsize=16)

# Add legend
plt.legend()

plt.grid(False)

# Save figure
if save_figures:
    file_name = f"Figure_1D.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=600)

# Show the plot
plt.show()

#%%
#%% Save the plotted time series to CSV for further use

import pandas as pd
import os

save_csv_dir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figure_1_plot_data'
os.makedirs(save_csv_dir, exist_ok=True)

# Save all_rs_v1_patch_4vtx: shape (sessions, targets, timepoints)
n_sessions, n_targets, n_timepoints = all_rs_v1_patch_4vtx.shape
for trgIdx in range(n_targets):
    df = pd.DataFrame(all_rs_v1_patch_4vtx[:, trgIdx, :].reshape(n_sessions, n_timepoints))
    df.to_csv(os.path.join(save_csv_dir, f'all_rs_v1_patch_4vtx_target{trgIdx}.csv'), index=False)

# Save all_rs_v1_patch_4V1: shape (sessions, targets, vertices, timepoints)
# For mean across vertices (used in Figure 1D)
mean_all_rs_v1_patch_4V1 = np.mean(all_rs_v1_patch_4V1, axis=2)  # average across vertices
for trgIdx in range(mean_all_rs_v1_patch_4V1.shape[1]):
    df = pd.DataFrame(mean_all_rs_v1_patch_4V1[:, trgIdx, :])
    df.to_csv(os.path.join(save_csv_dir, f'all_rs_v1_patch_4V1_target{trgIdx}.csv'), index=False)

print(f"CSV data saved to {save_csv_dir}")


