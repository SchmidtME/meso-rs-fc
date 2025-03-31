#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 10:50:46 2025

This script corresponds to Analysis 1 - the effect of distance and ocular polarity
on rs-FC and selectivity. It is generating Figure 2 for the manuscript.

@author: ms1454
"""

import os
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import warnings
warnings.filterwarnings("ignore")

#%% specify input

Sbjs = ['haas', 'chss', 'aroo', 'aman', 'ylri', 'rcgr', 'atib', 'evad', 'arak', 'imyy', 'auil']
ROI = 'V1'
layer = f'{ROI}_layers_0-10_intrahemispheric'

dataDir = os.path.join('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results', layer)

num_quantiles = 10 # number of distance quantiles

save_figures = 1
if save_figures:
    date_stamp = datetime.now().strftime('%Y-%m-%d')
    save_dir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Figure_2/{date_stamp}'
    os.makedirs(save_dir, exist_ok=True)

#%% load data

# Initialize matrices
alike_values, unalike_values, diff_values = (np.zeros((len(Sbjs), num_quantiles)) for _ in range(3))
FC_values = np.zeros((len(Sbjs), num_quantiles))  # Overall FC for normalization
Sel_avg_norm = np.zeros((num_quantiles, 2))

for sub_index, subName in enumerate(Sbjs):
    subDir = os.path.join(dataDir, subName)

    # Load Selectivity Data
    Sel_data = sio.loadmat(os.path.join(subDir, 'CorrelationMtx_Selectivity_subsampled_1.mat'))['Data_Combined']
    # Load Functional Connectivity Data
    FC_data = sio.loadmat(os.path.join(subDir, 'CorrelationMtx_FC_subsampled_1.mat'))['Data_Combined']
    # shape of data: 10 distance quantiles, 4 (both eyes, left2left, right2right, left2right), 2 hemis, 10, 10 (ODI), 2 (r or z value)
    
    # take the mean over ODI quantiles & hemispheres
    FC_avg = np.transpose(np.mean(FC_data[:,[0,3],:,:,:,1], axis=(2,3,4)))
    FC_values[sub_index, :] = FC_avg[:, 0]

    Sel_avg = np.transpose(np.mean(Sel_data[:,[0,3],:,:,:,1], axis=(2,3,4)))
    
    # Normalize and compute difference
    Sel_avg_norm[:, 0] = (Sel_avg[:,0] / FC_values[sub_index, :]) * 200
    Sel_avg_norm[:, 1] = (Sel_avg[:, 1] / FC_values[sub_index, :]) * 200
    
    # Store values for subjects
    alike_values[sub_index, :], unalike_values[sub_index, :] = Sel_avg[:, 0], Sel_avg[:, 1]
    diff_values[sub_index, :] = Sel_avg_norm[:, 0] - Sel_avg_norm[:, 1]

# Compute mean and standard error over subjects
avg_alike, sem_alike = np.mean(alike_values, axis=0), np.std(alike_values, axis=0) / np.sqrt(len(Sbjs))
avg_unalike, sem_unalike = np.mean(unalike_values, axis=0), np.std(unalike_values, axis=0) / np.sqrt(len(Sbjs))
avg_diff, sem_diff = np.mean(diff_values, axis=0), np.std(diff_values, axis=0) / np.sqrt(len(Sbjs))

print("Mean selectivity:", np.mean(avg_diff))
print("SEM selectivity:", np.std(avg_diff) / np.sqrt(len(Sbjs)))

#%% plot rs-FC and selectivity

bar_width = 1
cap_width = 9

x = np.arange(avg_alike.shape[0]) * 2.7  # Adjust spacing for alike/unalike plot
x_diff = np.arange(avg_diff.shape[0]) * 2  # Adjust spacing for selectivity plot

# Create figure with different widths for left and right subplots
fig = plt.figure(figsize=(16, 6), constrained_layout=True)
gs = fig.add_gridspec(1, 2, width_ratios=[1.9, 1.1])  # Left plot wider, right plot narrower

# Alike and Unalike Subplot (left)
ax1 = fig.add_subplot(gs[0])
positions_mid_alike = x - 0.5 * bar_width
positions_mid_unalike = x + 0.5 * bar_width

ax1.errorbar(positions_mid_alike, avg_alike, yerr=sem_alike, fmt='none', 
             color='black', capsize=cap_width, elinewidth=2, zorder=0)
deep_alike_bar = ax1.bar(positions_mid_alike, avg_alike, bar_width, 
                         label='Alike', color='white', edgecolor='black', linewidth=1.5)

ax1.errorbar(positions_mid_unalike, avg_unalike, yerr=sem_unalike, fmt='none', 
             color='black', capsize=cap_width, elinewidth=2, zorder=0)
deep_unalike_bar = ax1.bar(positions_mid_unalike, avg_unalike, bar_width, 
                           label='Unalike', color='black', edgecolor='black', linewidth=1.5)

mid_x_ticks = (positions_mid_alike + positions_mid_unalike) / 2
ax1.set_xticks(mid_x_ticks)
ax1.set_xticklabels(range(1, 11))
ax1.set_xlabel('Distance')
ax1.set_ylim([0.14, 0.40])
ax1.set_xlim(left=positions_mid_alike.min() - 1.2, right=positions_mid_unalike.max() + 1.2)
ax1.set_ylabel('z-transformed Pearson correlation coefficient')
ax1.set_title('Effect of cortical depth and distance on rs-FC strength (Alike vs Unalike)')
ax1.legend(handles=[deep_alike_bar, deep_unalike_bar])

# Set y-ticks for the left plot: Bottom, Middle, and Top
ax1.set_yticks([0.14, 0.27, 0.4])

# Selectivity Subplot (right)
ax2 = fig.add_subplot(gs[1])
x_diff = np.arange(avg_diff.shape[0]) * 1.6
positions_deep_alike = x_diff - 1 * bar_width

ax2.errorbar(positions_deep_alike, avg_diff, yerr=sem_diff, fmt='none', 
             color='black', capsize=cap_width, elinewidth=2, zorder=0)
deep_diff_bar = ax2.bar(positions_deep_alike, avg_diff, bar_width, 
                        label='Selectivity', color='#808080', edgecolor='black', linewidth=1.5)

ax2.set_xticks(positions_deep_alike)
ax2.set_xticklabels(range(1, 11))
ax2.set_xlabel('Distance')
ax2.set_ylim([0, 12])
ylabel_text = ('% change of Fishers z (normalized)')
ax2.set_ylabel(ylabel_text)
ax2.set_title('Effect of cortical depth and distance on selectivity')
ax2.legend(handles=[deep_diff_bar])

# Set y-ticks for the right plot: Bottom, Middle, and Top
ax2.set_yticks([0, 6, 12])

for ax in [ax1, ax2]:
    ax.spines['bottom'].set_linewidth(2)  # Increase x-axis linewidth
    ax.spines['left'].set_linewidth(2)    # Increase y-axis linewidth
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.tick_params(axis='both', labelsize=16, width=2, color='black')

plt.tight_layout()

if save_figures:
    file_name = f"Figure_2.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=200)
    
plt.show()
 
# %% clean version of the plot that is to be modified in powerpoint

# Create figure with different widths for left and right subplots
if save_figures:
    fig = plt.figure(figsize=(16, 6), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.9, 1.1])  # Left plot wider, right plot narrower
    
    # Alike and Unalike Subplot (left)
    ax1 = fig.add_subplot(gs[0])
    positions_mid_alike = x - 0.5 * bar_width
    positions_mid_unalike = x + 0.5 * bar_width
    
    ax1.errorbar(positions_mid_alike, avg_alike, yerr=sem_alike, fmt='none', 
                 color='black', capsize=cap_width, capthick=3, elinewidth=3, zorder=0)
    deep_alike_bar = ax1.bar(positions_mid_alike, avg_alike, bar_width, 
                             label='Alike', color='white', edgecolor='black', linewidth=2)
    
    ax1.errorbar(positions_mid_unalike, avg_unalike, yerr=sem_unalike, fmt='none', 
                 color='black', capsize=cap_width, capthick=3, elinewidth=3, zorder=0)
    deep_unalike_bar = ax1.bar(positions_mid_unalike, avg_unalike, bar_width, 
                               label='Unalike', color='black', edgecolor='black', linewidth=2)
    
    mid_x_ticks = (positions_mid_alike + positions_mid_unalike) / 2
    ax1.set_xticks(mid_x_ticks)
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.set_ylim([0.14, 0.40])
    ax1.set_xlim(left=positions_mid_alike.min() - 1.2, right=positions_mid_unalike.max() + 1.2)
    
    # Set y-ticks for the left plot: Bottom, Middle, and Top
    ax1.set_yticks([0.14, 0.27, 0.4])
    
    # Selectivity Subplot (right)
    ax2 = fig.add_subplot(gs[1])
    x_diff = np.arange(avg_diff.shape[0]) * 1.6
    positions_deep_alike = x_diff - 1 * bar_width
    
    ax2.errorbar(positions_deep_alike, avg_diff, yerr=sem_diff, 
    fmt='none', ecolor='black', elinewidth=3, capsize=cap_width, capthick=3,zorder=0)
    deep_diff_bar = ax2.bar(positions_deep_alike, avg_diff, bar_width, color='#808080', ecolor='black', edgecolor='black', linewidth=2)
    
    ax2.set_xticks(positions_deep_alike)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.set_ylim([0, 12])
    
    # Set y-ticks for the right plot: Bottom, Middle, and Top
    ax2.set_yticks([0, 6, 12])
    
    for ax in [ax1, ax2]:
        ax.spines['bottom'].set_linewidth(3)  # Increase x-axis linewidth
        ax.spines['left'].set_linewidth(3)    # Increase y-axis linewidth
        ax.spines['bottom'].set_color('black')
        ax.spines['left'].set_color('black')
        ax.tick_params(axis='both', labelsize=16, width=2, color='black')
        
    plt.tight_layout()
    file_name = f"Figure_2_clean.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=200)
