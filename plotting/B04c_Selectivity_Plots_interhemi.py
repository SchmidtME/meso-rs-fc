#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is plotting correlation matrices for rs-FC and selectivity for each
cortical depth level and the respective V1 subregion between hemispheres. 
It corresponds to Figure 6 & 7 of the manuscript.

@author: Marianna Elisa Schmidt (marianna.schmidt@maxplanckschools.de)
"""

import os
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.io as sio
import warnings
warnings.filterwarnings("ignore")


save_figures = 1
if save_figures:
    date_stamp = datetime.now().strftime('%Y-%m-%d')
    save_dir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Figure_5/{date_stamp}'
    os.makedirs(save_dir, exist_ok=True)
    
ROI = 'V1_Ventral'

saveDir = "/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figure_4-7_plot_data"
os.makedirs(saveDir, exist_ok=True)

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
   'auil',
   'myla']

baseDir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results'

layers = [
    f'{ROI}_layers_8-10_interhemispheric',
    f'{ROI}_layers_4-6_interhemispheric',
    f'{ROI}_layers_0-2_interhemispheric'
]

num_layers = len(layers)
num_quantiles = 10 #len(beta_comp)

short_layers = ['8-10', '4-6', '0-2']
short_quantiles = ['10', '9', '8', '7', '6', '5', '4', '3', '2', '1']
normalized_diff = 0
fisher_transform = 1
eye = 0 # 0 = both, 1 = eye1 (negative beta, left eye), 2 = eye2 (positive beta, right eye)
vmax_corr = 0.22
vmin_corr = 0.12
vmin_diff = -0.01
vmax_diff = 0.01

avg_values_alike = np.zeros((num_layers, num_quantiles, num_quantiles))
avg_values_unalike = np.zeros((num_layers, num_quantiles, num_quantiles))
avg_values_diff = np.zeros((num_layers, num_quantiles, num_quantiles))
avg_values_FC = np.zeros((num_layers, num_quantiles, num_quantiles))

print("Mean values for subjects per version to check for outlier values")
for version_index, (version, short_version) in enumerate(zip(layers, short_layers)):
    dataDir = os.path.join(baseDir, version)
    alike_values = np.zeros((len(subjectFolders), num_quantiles, num_quantiles))
    unalike_values = np.zeros((len(subjectFolders), num_quantiles, num_quantiles))
    diff_values = np.zeros((len(subjectFolders), num_quantiles, num_quantiles))
    
    FC_values = np.zeros((len(subjectFolders), num_quantiles, num_quantiles))
            
    for sub_index, subName in enumerate(subjectFolders):
        
        subDir = os.path.join(dataDir, subName)
        
        # load overall FC
        data_mat_name = 'CorrelationMtx_FC_beta.mat'
        data_mat_path = os.path.join(subDir, data_mat_name)
        
        if not os.path.exists(data_mat_path):
            print(f"File not found: {data_mat_path}")
            continue
        
        data_mat = sio.loadmat(data_mat_path)
        Data_Combined = np.squeeze(data_mat['Data_Combined'])
        
        # Apply Fisher Transform or not
        Data_Combined = Data_Combined[:, :, 1] if fisher_transform else Data_Combined[ :, :, 0]
        
        # Average across hemispheres (axis=1) and distance quantiles (axis=0)
        Data_Combined_avg = Data_Combined
        
        FC_values[sub_index, :, :] = Data_Combined_avg
        
        
        subDir = os.path.join(dataDir, subName)
        data_mat_name = 'CorrelationMtx_Selectivity.mat'
        data_mat_path = os.path.join(subDir, data_mat_name)

        if os.path.exists(data_mat_path):
            data_mat = sio.loadmat(data_mat_path)
            Data_Combined = np.squeeze(data_mat['Data_Combined'])
            if fisher_transform == 1:
                Data_Combined = Data_Combined[:,:,:,1]
            else:
                Data_Combined = Data_Combined[:,:,:,0]
            Data_Combined_test = np.squeeze(data_mat['Data_Combined'])

            Data_Combined2 = Data_Combined
            del Data_Combined
            Data_Combined = Data_Combined2[[eye,3],:,:]

            Data_Combined_avg = Data_Combined
            # Store values for each quantile
            alike_values[sub_index, :,:] = Data_Combined_avg[0,:,:]
            unalike_values[sub_index, :,:] = Data_Combined_avg[1,:,:]
            diff_values[sub_index, :,:] = alike_values[sub_index, :,:] - unalike_values[sub_index, :,:]
            print(np.min(np.min(diff_values[sub_index,:,:], axis=1),axis=0))

    avg_values_alike[version_index, :,:] = np.mean(alike_values, axis=0)
    avg_values_unalike[version_index, :,:] = np.mean(unalike_values, axis=0)
    avg_values_diff[version_index, :,:] = np.mean(diff_values, axis=0)
    avg_values_FC[version_index, :,:] = np.mean(FC_values, axis=0)
    
#%%

fig, axes = plt.subplots(3, 2, figsize=(11, 15), gridspec_kw={'height_ratios': [1, 1, 1], 'hspace': 0.2, 'wspace': 0})

# Titles on the left of each row
fig.text(0.07, 0.77, 'Superficial', fontsize=23, ha='center', va='center', rotation='vertical', weight='bold')
fig.text(0.07, 0.50, 'Middle', fontsize=23, ha='center', va='center', rotation='vertical', weight='bold')
fig.text(0.07, 0.22, 'Deep', fontsize=23, ha='center', va='center', rotation='vertical', weight='bold')

# Titles on top of each column
fig.text(0.32, 0.9, 'rs-FC', fontsize=23, ha='center', va='center', weight='bold')
fig.text(0.7, 0.9, 'Alike - Unalike', fontsize=23, ha='center', va='center', weight='bold')

# Iterate over each version and its corresponding index
for version_index, (version, short_version) in enumerate(zip(layers, short_layers)):
    # First column: Alike ODC
    cax1 = axes[version_index, 0].imshow(np.transpose(avg_values_FC[version_index, :,::-1]), cmap='gist_heat', vmin=vmin_corr, vmax=vmax_corr)
    axes[version_index, 0].set_yticks(np.arange(0, num_quantiles))
    axes[version_index, 0].set_yticklabels(short_quantiles[::-1])
    axes[version_index, 0].set_xticks(np.arange(0, num_quantiles))
    axes[version_index, 0].set_xticklabels(short_quantiles)
    axes[version_index, 0].tick_params(axis='x', labelsize=12)
    axes[version_index, 0].tick_params(axis='y', labelsize=12)
    if version_index == 2:  # Bottom row
        axes[version_index, 0].set_xlabel('Beta quantile 1', fontsize=15)
    if version_index == 0:  # Leftmost column
        axes[version_index, 0].set_ylabel('Beta quantile 2', fontsize=15)
    axes[version_index, 0].grid(False)

    # Third column: Difference
    cmap_diff = cm.get_cmap('coolwarm')
    cax3 = axes[version_index, 1].imshow(np.transpose(avg_values_diff[version_index, :,::-1]), cmap=cmap_diff, vmin=vmin_diff, vmax=vmax_diff)
    axes[version_index, 1].set_yticks(np.arange(0, num_quantiles))
    axes[version_index, 1].set_yticklabels(short_quantiles[::-1])
    axes[version_index, 1].set_xticks(np.arange(0, num_quantiles))
    axes[version_index, 1].set_xticklabels(short_quantiles)
    axes[version_index, 1].tick_params(axis='x', labelsize=12)
    axes[version_index, 1].tick_params(axis='y', labelsize=12)
    if version_index == 1:  # Bottom row
        axes[version_index, 1].set_xlabel('Beta quantile 1', fontsize=15)
    if version_index == 2:  # Leftmost column
        axes[version_index, 0].set_ylabel('Beta quantile 2', fontsize=15)
    axes[version_index, 1].grid(False)

# Adding shared colorbars at the bottom of each column
cbar_ax1 = fig.add_axes([0.17, 0.05, 0.3, 0.02])  # Bottom of the first column
cbar1 = fig.colorbar(cax1, cax=cbar_ax1, orientation='horizontal')
if fisher_transform == 1:
    cbar1.set_label('z-transformed Pearson correlation coefficient', fontsize=15)
else:
    cbar1.set_label('Pearson correlation coefficient', fontsize=15)
cbar1.set_ticks([vmin_corr, vmin_corr+(vmax_corr-vmin_corr)/4, vmin_corr+(vmax_corr-vmin_corr)/2, vmin_corr+(vmax_corr-vmin_corr)/4*3, vmax_corr])
cbar1.ax.tick_params(labelsize=14)  # Increase tick font size

cbar_ax3 = fig.add_axes([0.56, 0.05, 0.3, 0.02])  # Bottom of the third column
cbar3 = fig.colorbar(cax3, cax=cbar_ax3, orientation='horizontal')
cbar3.set_label('Percentage change', fontsize=15)
cbar3.set_ticks([vmin_diff, vmin_diff/2, 0, vmax_diff/2, vmax_diff])
cbar3.ax.tick_params(labelsize=14)  # Increase tick font size

plt.tight_layout()

if save_figures:
    file_name = f"Figure_5_{ROI}.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=200)
    
plt.show()

#%% clean figure

if save_figures:

    fig, axes = plt.subplots(3, 2, figsize=(7, 10), gridspec_kw={'height_ratios': [1, 1, 1], 'hspace': 0.1, 'wspace': 0.1})

    # Iterate over each version and its corresponding index
    for version_index, (version, short_version) in enumerate(zip(layers, short_layers)):
        # First column: Alike ODC
        cax1 = axes[version_index, 0].imshow(np.transpose(avg_values_FC[version_index, :,::-1]), cmap='gist_heat', vmin=vmin_corr, vmax=vmax_corr)
        axes[version_index, 0].set_yticks(np.arange(0, num_quantiles))
        axes[version_index, 0].set_yticklabels(short_quantiles[::-1])
        axes[version_index, 0].set_xticks(np.arange(0, num_quantiles))
        axes[version_index, 0].set_xticklabels(short_quantiles)
        axes[version_index, 0].tick_params(axis='x', labelsize=12)
        axes[version_index, 0].tick_params(axis='y', labelsize=12)
        
        
        # Third column: Difference
        cmap_diff = cm.get_cmap('coolwarm')
        cax3 = axes[version_index, 1].imshow(np.transpose(avg_values_diff[version_index, :,::-1]), cmap=cmap_diff, vmin=vmin_diff, vmax=vmax_diff)
        axes[version_index, 1].set_yticks(np.arange(0, num_quantiles))
        axes[version_index, 1].set_yticklabels(short_quantiles[::-1])
        axes[version_index, 1].set_xticks(np.arange(0, num_quantiles))
        axes[version_index, 1].set_xticklabels(short_quantiles)
        axes[version_index, 1].tick_params(axis='x', labelsize=12)
        axes[version_index, 1].tick_params(axis='y', labelsize=12)
        
        # Remove labels but keep the ticks
        for ax in axes[version_index]:
            ax.set_xticklabels([])  # Remove x-axis labels
            ax.set_yticklabels([])  # Remove y-axis ticks
            ax.grid(False)
            
    plt.tight_layout()

    file_name = f"Figure_5_clean_{ROI}.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=400)
    
# Save data to .mat file
output_path = os.path.join(saveDir, f'Figure_4-7_V1_Peripheral_interhemispheric_data.mat')
sio.savemat(output_path, {
    'avg_values_alike': avg_values_alike,
    'avg_values_unalike': avg_values_unalike,
    'avg_values_diff': avg_values_diff,
    'avg_values_FC': avg_values_FC,
    'layers': short_layers,
    'quantiles': short_quantiles
})

print(f"Saved processed data to: {output_path}")
    