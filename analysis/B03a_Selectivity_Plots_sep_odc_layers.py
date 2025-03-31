#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 09:10:07 2025

@author: ms1454
"""

import os
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.io as sio
import warnings
warnings.filterwarnings("ignore")

#%% specifications

save_figures = 1
if save_figures:
    date_stamp = datetime.now().strftime('%Y-%m-%d')
    save_dir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Control_Analyses/odc_sep_layers/{date_stamp}'
    os.makedirs(save_dir, exist_ok=True)
    
subjectFolders = ['haas',
    'chss',
    'aroo',
    'aman',
    'ylri',
    'rcgr',
    'atib',
    'evad',
    'arak',
    'imyy',
    'auil']
    
ROI = 'V1'

baseDir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results'

layers = [
    f'{ROI}_layers_8-10_intrahemispheric_odc_sep_depths',
    f'{ROI}_layers_4-6_intrahemispheric_odc_sep_depths',
    f'{ROI}_layers_0-2_intrahemispheric_odc_sep_depths'
]

num_layers = len(layers)
num_quantiles = 10 #len(beta_comp)

short_layers = ['8-10', '4-6', '0-2']
short_quantiles = ['10', '9', '8', '7', '6', '5', '4', '3', '2', '1']
# normalized_diff = 1
fisher_transform = 1
eye = 0 # 0 = both, 1 = eye1 (negative beta, right), 2 = eye2 (positive beta, left eye)
exclude_dist_quant = 0 # exclude number of distance quantiles before averaging

vmax_corr = 0.21#0.185 
vmin_corr = 0.14
vmax_diff = 0.025
vmax_diff_norm = 8

avg_values_alike = np.zeros((num_layers, num_quantiles, num_quantiles))
avg_values_unalike = np.zeros((num_layers, num_quantiles, num_quantiles))
avg_values_diff = np.zeros((num_layers, num_quantiles, num_quantiles))

#%% load data

for version_index, (version, short_version) in enumerate(zip(layers, short_layers)):
    dataDir = os.path.join(baseDir, version)
   
    alike_values = np.zeros((len(subjectFolders), num_quantiles, num_quantiles))
    unalike_values = np.zeros((len(subjectFolders), num_quantiles, num_quantiles))
    diff_values = np.zeros((len(subjectFolders), num_quantiles, num_quantiles))

    for sub_index, subName in enumerate(subjectFolders):
        subDir = os.path.join(dataDir, subName)
        data_mat_name = 'CorrelationMtx_Selectivity_subsampled_1.mat'
        data_mat_path = os.path.join(subDir, data_mat_name)

        if os.path.exists(data_mat_path):
            data_mat = sio.loadmat(data_mat_path)
            Data_Combined = np.squeeze(data_mat['Data_Combined'])
            if fisher_transform == 1:
                Data_Combined = Data_Combined[:,:,:,:,:,1]
            else:
                Data_Combined = Data_Combined[:,:,:,:,:,0]
            Data_Combined_test = np.squeeze(data_mat['Data_Combined'])

            Data_Combined2 = Data_Combined
            del Data_Combined
            Data_Combined = Data_Combined2[:,[eye,3],:,:,:]
            
            # Average across hemispheres (axis 2) and distances
            Data_Combined_avg = np.mean(np.mean(Data_Combined[exclude_dist_quant:len(short_quantiles),:,:,:,:], axis=2), axis=0)
            
            # Store values for each quantile
            alike_values[sub_index, :,:] = Data_Combined_avg[0,:,:]
            unalike_values[sub_index, :,:] = Data_Combined_avg[1,:,:]
            diff_values[sub_index, :,:] = alike_values[sub_index, :,:] - unalike_values[sub_index, :,:]
            
    avg_values_alike[version_index, :,:] = np.mean(alike_values, axis=0)
    avg_values_unalike[version_index, :,:] = np.mean(unalike_values, axis=0)
    avg_values_diff[version_index, :,:] = np.mean(diff_values, axis=0)
    
#%% figure with labels

fig, axes = plt.subplots(3, 3, figsize=(14, 16), gridspec_kw={'height_ratios': [1, 1, 1], 'hspace': 0.2, 'wspace': 0.1})

# Titles on the left of each row
fig.text(0.05, 0.77, 'Superficial', fontsize=23, ha='center', va='center', rotation='vertical', weight='bold')
fig.text(0.05, 0.50, 'Middle', fontsize=23, ha='center', va='center', rotation='vertical', weight='bold')
fig.text(0.05, 0.22, 'Deep', fontsize=23, ha='center', va='center', rotation='vertical', weight='bold')

# Titles on top of each column
fig.text(0.25, 0.9, 'Alike', fontsize=23, ha='center', va='center', weight='bold')
fig.text(0.51, 0.9, 'Unalike', fontsize=23, ha='center', va='center', weight='bold')
fig.text(0.78, 0.9, 'Alike-Unalike', fontsize=23, ha='center', va='center', weight='bold')

for version_index, (version, short_version) in enumerate(zip(layers, short_layers)):
    
    # First column: Alike ODC
    cax1 = axes[version_index, 0].imshow(np.transpose(avg_values_alike[version_index, :,::-1]), cmap='gist_heat', vmin=vmin_corr, vmax=vmax_corr)
    axes[version_index, 0].set_xticks(np.arange(0, num_quantiles))
    axes[version_index, 0].set_xticklabels(short_quantiles[::-1], fontsize=25)
    axes[version_index, 0].set_yticks(np.arange(0, num_quantiles))
    axes[version_index, 0].set_yticklabels(short_quantiles, fontsize=25)
    axes[version_index, 0].tick_params(axis='x', labelsize=25)
    axes[version_index, 0].tick_params(axis='y', labelsize=25)
    if version_index == 2:  # Bottom row
        axes[version_index, 0].set_xlabel('Beta quantile 1', fontsize=25)
    if version_index == 0:  # Leftmost column
        axes[version_index, 0].set_ylabel('Beta quantile 2', fontsize=25)

    # Second column: Unalike ODC
    cax2 = axes[version_index, 1].imshow(np.transpose(avg_values_unalike[version_index, :,::-1]), cmap='gist_heat', vmin=vmin_corr, vmax=vmax_corr)
    axes[version_index, 1].set_xticks(np.arange(0, num_quantiles))
    axes[version_index, 1].set_xticklabels(short_quantiles[::-1], fontsize=25)
    axes[version_index, 1].set_yticks(np.arange(0, num_quantiles))
    axes[version_index, 1].set_yticklabels(short_quantiles, fontsize=25)
    axes[version_index, 1].tick_params(axis='x', labelsize=25)
    axes[version_index, 1].tick_params(axis='y', labelsize=25)
    if version_index == 2:  # Bottom row
        axes[version_index, 1].set_xlabel('Beta quantile 1', fontsize=25)
    if version_index == 1:  # Leftmost column
        axes[version_index, 0].set_ylabel('Beta quantile 2', fontsize=25)

    # Third column: Not normalized Difference
    cmap_diff = cm.get_cmap('coolwarm')
    cax3 = axes[version_index, 2].imshow(np.transpose(avg_values_diff[version_index, :,::-1]), cmap=cmap_diff, vmin=-vmax_diff, vmax=vmax_diff)
    axes[version_index, 2].set_xticks(np.arange(0, num_quantiles))
    axes[version_index, 2].set_xticklabels(short_quantiles[::-1], fontsize=25)
    axes[version_index, 2].set_yticks(np.arange(0, num_quantiles))
    axes[version_index, 2].set_yticklabels(short_quantiles, fontsize=25)
    axes[version_index, 2].tick_params(axis='x', labelsize=25)
    axes[version_index, 2].tick_params(axis='y', labelsize=25)
    if version_index == 2:  # Bottom row
        axes[version_index, 2].set_xlabel('Beta quantile 1', fontsize=25)
    if version_index == 2:  # Leftmost column
        axes[version_index, 0].set_ylabel('Beta quantile 2', fontsize=25)

    # Disable gridlines
    for ax in axes[version_index]:
        ax.grid(False)

# Adding shared colorbars at the bottom of each column
cbar_ax1 = fig.add_axes([0.27, 0.01, 0.2, 0.02])  # Bottom of the first column
cbar1 = fig.colorbar(cax1, cax=cbar_ax1, orientation='horizontal')
if fisher_transform == 1:
    cbar1.set_label('z-transformed Pearson r', fontsize=25)
else:
    cbar1.set_label('Pearson r', fontsize=25)
cbar1.set_ticks([vmin_corr, vmin_corr+(vmax_corr-vmin_corr)/2, vmax_corr])
cbar1.ax.tick_params(labelsize=20)  # Increase tick font size

cbar_ax3 = fig.add_axes([0.68, 0.0, 0.2, 0.02])  # Bottom of the third column
cbar3 = fig.colorbar(cax3, cax=cbar_ax3, orientation='horizontal')
cbar3.set_label('Difference', fontsize=25)
cbar3.set_ticks([-vmax_diff, 0, vmax_diff])
cbar3.ax.tick_params(labelsize=20)  # Increase tick font size

plt.tight_layout()

if save_figures:
    file_name = f"Figure_3_{ROI}.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=200)
  
plt.show()

#%% figure without labels

if save_figures:
    
    fig, axes = plt.subplots(3, 3, figsize=(10, 10), gridspec_kw={'height_ratios': [1, 1, 1], 'hspace': 0.1, 'wspace': 0.1})
    
    # Iterate over each version and its corresponding index
    for version_index, (version, short_version) in enumerate(zip(layers, short_layers)):
        
        # First column: Alike ODC       
        cax1 = axes[version_index, 0].imshow(np.transpose(avg_values_alike[version_index, :,::-1]), cmap='gist_heat', vmin=vmin_corr, vmax=vmax_corr)
        axes[version_index, 0].set_yticks(np.arange(0, num_quantiles))
        axes[version_index, 0].set_xticks(np.arange(0, num_quantiles))
    
        # Second column: Unalike ODC
        cax2 = axes[version_index, 1].imshow(np.transpose(avg_values_unalike[version_index, :,::-1]), cmap='gist_heat', vmin=vmin_corr, vmax=vmax_corr)
        axes[version_index, 1].set_xticks(np.arange(0, num_quantiles))
        axes[version_index, 1].set_yticks(np.arange(0, num_quantiles))
    
        # Third column: Not normalized Difference
        cmap_diff = cm.get_cmap('coolwarm')
        cax3 = axes[version_index, 2].imshow(np.transpose(avg_values_diff[version_index, :,::-1]), cmap=cmap_diff, vmin=-vmax_diff, vmax=vmax_diff)
        axes[version_index, 2].set_xticks(np.arange(0, num_quantiles))
        axes[version_index, 2].set_yticks(np.arange(0, num_quantiles))
    
        # Remove labels but keep the ticks
        for ax in axes[version_index]:
            
            ax.set_xticklabels([])  # Remove x-axis labels
            ax.set_yticklabels([])  # Remove y-axis ticks
            ax.grid(False)
            
    plt.tight_layout()

    file_name = f"Figure_3_{ROI}_clean.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=200)



