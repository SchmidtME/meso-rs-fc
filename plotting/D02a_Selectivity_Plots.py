#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 15:33:39 2024

This script corresponds to Analysis 2 - the effect of cortical depth, ODI and
ocular polarity on rs-FC and selectivity. It creates Figure 3 oof the 
manuscript. 

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
 
Sbjs = ['haas', 'chss', 'aroo', 'aman', 'ylri', 'rcgr', 'atib', 'evad', 'arak', 'imyy', 'auil']
    
ROI = 'V1'

baseDir = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results'

layers = [
    f'{ROI}_layers_8-10_intrahemispheric',
    f'{ROI}_layers_4-6_intrahemispheric',
    f'{ROI}_layers_0-2_intrahemispheric'
]

num_quantiles = 10 # number of ODI/beta quantiles

save_figures = 1
if save_figures:
    date_stamp = datetime.now().strftime('%Y-%m-%d')
    save_dir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Figure_3/{date_stamp}'
    os.makedirs(save_dir, exist_ok=True)

# ylims for figures
vmax_corr = 0.185 
vmin_corr = 0.14
vmax_diff = 0.01
vmax_diff_norm = 8


#%% load data

avg_values_alike = np.zeros((len(layers), num_quantiles, num_quantiles))
avg_values_unalike = np.zeros((len(layers), num_quantiles, num_quantiles))
avg_values_diff = np.zeros((len(layers), num_quantiles, num_quantiles))

for layer_index, layer in enumerate(layers):
    
    dataDir = os.path.join(baseDir, layer)
   
    alike_values = np.zeros((len(Sbjs), num_quantiles, num_quantiles))
    unalike_values = np.zeros((len(Sbjs), num_quantiles, num_quantiles))
    diff_values = np.zeros((len(Sbjs), num_quantiles, num_quantiles))

    for sub_index, subName in enumerate(Sbjs):
        
        subDir = os.path.join(dataDir, subName)

        data_mat_path = os.path.join(subDir, 'CorrelationMtx_Selectivity_subsampled_1.mat')

        Data_Combined = np.squeeze(sio.loadmat(data_mat_path)['Data_Combined'])
        # shape of data: 10 distance quantiles, 4 (both eyes, left2left, right2right, left2right), 2 hemis, 10, 10 (ODI), 2 (r or z value)

        # take the mean over distance and hemispheres
        Data_Combined_avg = np.mean(Data_Combined[:,[0,3],:,:,:,1], axis=(1,2))
        
        # Store values for each subject
        alike_values[sub_index, :,:] = Data_Combined_avg[0,:,:]
        unalike_values[sub_index, :,:] = Data_Combined_avg[1,:,:]
        diff_values[sub_index, :,:] = alike_values[sub_index, :,:] - unalike_values[sub_index, :,:]
     
    # take the average over subjetcs and store for each layer
    avg_values_alike[layer_index, :,:] = np.mean(alike_values, axis=0)
    avg_values_unalike[layer_index, :,:] = np.mean(unalike_values, axis=0)
    avg_values_diff[layer_index, :,:] = np.mean(diff_values, axis=0)
    
#%% figure with labels

fig, axes = plt.subplots(3, 3, figsize=(14, 13), gridspec_kw={'height_ratios': [1, 1, 1], 'hspace': 0.2, 'wspace': 0.1})

# Titles on the left of each row
fig.text(0.05, 0.77, 'Superficial', fontsize=23, ha='center', va='center', rotation='vertical', weight='bold')
fig.text(0.05, 0.50, 'Middle', fontsize=23, ha='center', va='center', rotation='vertical', weight='bold')
fig.text(0.05, 0.22, 'Deep', fontsize=23, ha='center', va='center', rotation='vertical', weight='bold')

# Titles on top of each column
fig.text(0.25, 0.9, 'Alike', fontsize=23, ha='center', va='center', weight='bold')
fig.text(0.51, 0.9, 'Unalike', fontsize=23, ha='center', va='center', weight='bold')
fig.text(0.78, 0.9, 'Alike-Unalike', fontsize=23, ha='center', va='center', weight='bold')

for layer_index, layer in enumerate(layers):
    
    # First column: Alike ODC
    cax1 = axes[layer_index, 0].imshow(np.transpose(avg_values_alike[layer_index, :,::-1]), cmap='gist_heat', vmin=vmin_corr, vmax=vmax_corr)
    axes[layer_index, 0].set_xticks(np.arange(0, num_quantiles))
    axes[layer_index, 0].set_xticklabels(['10', '9', '8', '7', '6', '5', '4', '3', '2', '1'][::-1], fontsize=25)
    axes[layer_index, 0].set_yticks(np.arange(0, num_quantiles))
    axes[layer_index, 0].set_yticklabels(['10', '9', '8', '7', '6', '5', '4', '3', '2', '1'], fontsize=25)
    axes[layer_index, 0].tick_params(axis='x', labelsize=25)
    axes[layer_index, 0].tick_params(axis='y', labelsize=25)
    if layer_index == 2:  # Bottom row
        axes[layer_index, 0].set_xlabel('Beta quantile 1', fontsize=25)
    if layer_index == 0:  # Leftmost column
        axes[layer_index, 0].set_ylabel('Beta quantile 2', fontsize=25)

    # Second column: Unalike ODC
    cax2 = axes[layer_index, 1].imshow(np.transpose(avg_values_unalike[layer_index, :,::-1]), cmap='gist_heat', vmin=vmin_corr, vmax=vmax_corr)
    axes[layer_index, 1].set_xticks(np.arange(0, num_quantiles))
    axes[layer_index, 1].set_xticklabels(['10', '9', '8', '7', '6', '5', '4', '3', '2', '1'][::-1], fontsize=25)
    axes[layer_index, 1].set_yticks(np.arange(0, num_quantiles))
    axes[layer_index, 1].set_yticklabels(['10', '9', '8', '7', '6', '5', '4', '3', '2', '1'], fontsize=25)
    axes[layer_index, 1].tick_params(axis='x', labelsize=25)
    axes[layer_index, 1].tick_params(axis='y', labelsize=25)
    if layer_index == 2:  # Bottom row
        axes[layer_index, 1].set_xlabel('Beta quantile 1', fontsize=25)
    if layer_index == 1:  # Leftmost column
        axes[layer_index, 0].set_ylabel('Beta quantile 2', fontsize=25)

    # Third column: Not normalized Difference
    cmap_diff = cm.get_cmap('coolwarm')
    cax3 = axes[layer_index, 2].imshow(np.transpose(avg_values_diff[layer_index, :,::-1]), cmap=cmap_diff, vmin=-vmax_diff, vmax=vmax_diff)
    axes[layer_index, 2].set_xticks(np.arange(0, num_quantiles))
    axes[layer_index, 2].set_xticklabels(['10', '9', '8', '7', '6', '5', '4', '3', '2', '1'][::-1], fontsize=25)
    axes[layer_index, 2].set_yticks(np.arange(0, num_quantiles))
    axes[layer_index, 2].set_yticklabels(['10', '9', '8', '7', '6', '5', '4', '3', '2', '1'], fontsize=25)
    axes[layer_index, 2].tick_params(axis='x', labelsize=25)
    axes[layer_index, 2].tick_params(axis='y', labelsize=25)
    if layer_index == 2:  # Bottom row
        axes[layer_index, 2].set_xlabel('Beta quantile 1', fontsize=25)
    if layer_index == 2:  # Leftmost column
        axes[layer_index, 0].set_ylabel('Beta quantile 2', fontsize=25)

    # Disable gridlines
    for ax in axes[layer_index]:
        ax.grid(False)

# Adding shared colorbars at the bottom of each column
cbar_ax1 = fig.add_axes([0.27, 0.01, 0.2, 0.02])  # Bottom of the first column
cbar1 = fig.colorbar(cax1, cax=cbar_ax1, orientation='horizontal')
cbar1.set_label('z-transformed Pearson r', fontsize=25)
cbar1.set_ticks([vmin_corr, vmin_corr+(vmax_corr-vmin_corr)/2, vmax_corr])
cbar1.ax.tick_params(labelsize=20)  # Increase tick font size

cbar_ax3 = fig.add_axes([0.68, 0.01, 0.2, 0.02])  # Bottom of the third column
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

#%% figure without labels to be modified in powerpoint

if save_figures:
    
    fig, axes = plt.subplots(3, 3, figsize=(10, 10), gridspec_kw={'height_ratios': [1, 1, 1], 'hspace': 0.1, 'wspace': 0.1})
    
    # Iterate over each layer and its corresponding index
    for layer_index, layer in enumerate(layers):
        
        # First column: Alike ODC       
        cax1 = axes[layer_index, 0].imshow(np.transpose(avg_values_alike[layer_index, :,::-1]), cmap='gist_heat', vmin=vmin_corr, vmax=vmax_corr)
        axes[layer_index, 0].set_yticks(np.arange(0, num_quantiles))
        axes[layer_index, 0].set_xticks(np.arange(0, num_quantiles))
    
        # Second column: Unalike ODC
        cax2 = axes[layer_index, 1].imshow(np.transpose(avg_values_unalike[layer_index, :,::-1]), cmap='gist_heat', vmin=vmin_corr, vmax=vmax_corr)
        axes[layer_index, 1].set_xticks(np.arange(0, num_quantiles))
        axes[layer_index, 1].set_yticks(np.arange(0, num_quantiles))
    
        # Third column: Not normalized Difference
        cmap_diff = cm.get_cmap('coolwarm')
        cax3 = axes[layer_index, 2].imshow(np.transpose(avg_values_diff[layer_index, :,::-1]), cmap=cmap_diff, vmin=-vmax_diff, vmax=vmax_diff)
        axes[layer_index, 2].set_xticks(np.arange(0, num_quantiles))
        axes[layer_index, 2].set_yticks(np.arange(0, num_quantiles))
    
        # Remove labels but keep the ticks
        for ax in axes[layer_index]:
            
            ax.set_xticklabels([])  # Remove x-axis labels
            ax.set_yticklabels([])  # Remove y-axis ticks
            ax.grid(False)
            
    plt.tight_layout()

    file_name = f"Figure_3_{ROI}_clean.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=200)

