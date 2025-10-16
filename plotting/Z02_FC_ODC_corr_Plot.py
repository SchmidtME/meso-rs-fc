#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is plotting correlation correlation coefficient between 
rs-FC rings of varying outer radii and either the corresponding,spatially
shifted or rotated ODC map. It corresponds to Figure 8 of the manuscript.

@author: Marianna Elisa Schmidt (marianna.schmidt@maxplanckschools.de)
"""

import os
from datetime import datetime
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel

#%% Define parameters and directories

# Flag to save generated figures
save_figures = 1
if save_figures:
    # Define directory for saving figures, create if it doesn't exist
    date_stamp = datetime.now().strftime('%Y-%m-%d')
    save_dir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Figure_6/{date_stamp}'
    os.makedirs(save_dir, exist_ok=True)

layer = '0-2'
# Base directory containing functional connectivity data for cortical depths
base_directory = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results/V1_layers_{layer}_intrahemispheric'

# Included subjects for this analysis
sub_names = ['haas',
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

# Initialize lists to store mean absolute R values and their rotated counterparts
mean_abs_r_all = []
mean_abs_r_rotated_all = []

#%% Load data for each subject

for sub in sub_names:
    # Load the mean absolute R values for the subject
    mat_file_path = os.path.join(base_directory, sub, 'meanAbsR_Denoised_excl_3mm.mat')
    if os.path.exists(mat_file_path):
        mat_data = sio.loadmat(mat_file_path)
        mean_abs_r_all.append(mat_data['meanAbsR'])
    
    # Load the rotated mean absolute R values for the subject
    rotated_file_path = os.path.join(base_directory, sub, 'meanAbsR_Denoised_excl_3mm_rotated.mat')
    if os.path.exists(rotated_file_path):
        rotated_data = sio.loadmat(rotated_file_path)
        mean_abs_r_rotated_all.append(rotated_data['meanAbsR'])

# Convert lists to numpy arrays for easier manipulation
mean_abs_r_all = np.stack(mean_abs_r_all, axis=2)
mean_abs_r_rotated_all = np.stack(mean_abs_r_rotated_all, axis=2)

#%% Compute mean and standard error across subjects

# Mean and standard error for the standard data
mean_abs_r = np.mean(np.mean(mean_abs_r_all, axis=2), axis=2)
std_abs_r = np.std(np.mean(mean_abs_r_all, axis=3), axis=2) / np.sqrt(mean_abs_r_all.shape[2])

# Mean and standard error for the rotated data
mean_abs_r_rotated = np.mean(np.mean(mean_abs_r_rotated_all, axis=2), axis=2)
std_abs_r_rotated = np.std(np.mean(mean_abs_r_rotated_all, axis=3), axis=2) / np.sqrt(mean_abs_r_rotated_all.shape[2])

# Number of quantiles (radius levels)
num_quantiles = mean_abs_r.shape[0]
x = np.arange(num_quantiles)  # Positions for quantiles on the x-axis

#%% Plot mean absolute R values and p-values

fig, axs = plt.subplots(1, 2, figsize=(12, 6), sharex=True)

# Plot mean absolute R values
axs[0].plot(x, mean_abs_r[:, 1], label='H1', color='black', linewidth=2, marker='o')
axs[0].fill_between(x, mean_abs_r[:, 1] - std_abs_r[:, 1], mean_abs_r[:, 1] + std_abs_r[:, 1], color='black', alpha=0.2)
axs[0].plot(x, mean_abs_r[:, 0], label='H0 shifted', color='blue', linewidth=2, marker='o')
axs[0].fill_between(x, mean_abs_r[:, 0] - std_abs_r[:, 0], mean_abs_r[:, 0] + std_abs_r[:, 0], color='blue', alpha=0.2)
axs[0].plot(x, mean_abs_r_rotated[:, 0], label='H0 rotated', color='red', linewidth=2, marker='o')
axs[0].fill_between(x, mean_abs_r_rotated[:, 0] - std_abs_r_rotated[:, 0], mean_abs_r_rotated[:, 0] + std_abs_r_rotated[:, 0], color='red', alpha=0.2)
axs[0].set_ylim(0.06, 0.17)
axs[0].set_xticks(x[::2])
axs[0].set_xticklabels(['4', '5', '6', '7', '8', '9', '10'])
axs[0].set_xlabel('Radius (mm)')
axs[0].set_ylabel('Pearson r')
axs[0].set_title('Cortical depth 0-10')
axs[0].legend()

# Compute p-values for paired t-tests
p_values = []
p_values_rot = []
for i in x:
    ar1 = np.mean(mean_abs_r_all[i, 1, :, :], axis=1)
    ar2 = np.mean(mean_abs_r_all[i, 0, :, :], axis=1)
    ar3 = np.mean(mean_abs_r_rotated_all[i, 0, :, :], axis=1)
    t_stat, p_value = ttest_rel(ar1, ar2)
    p_values.append(p_value)
    t_stat_rot, p_value_rot = ttest_rel(ar1, ar3)
    p_values_rot.append(p_value_rot)

# Convert p-values to -log10(p) for plotting
neg_log_p_values = -np.log10(np.array(p_values))
neg_log_p_values_rot = -np.log10(np.array(p_values_rot))

# Plot p-values
axs[1].plot(x, neg_log_p_values, label='H0 shifted', color='blue', linewidth=2, marker='o')
axs[1].plot(x, neg_log_p_values_rot, label='H0 rotated', color='red', linewidth=2, marker='o')
axs[1].axhline(y=-np.log10(0.05), color='black', linestyle='--', linewidth=1.5, label='p = 0.05')
axs[1].set_yticks([1.3, 2, 3, 4, 5, 6])
axs[1].set_xticks(x[::2])
axs[1].set_xticklabels(['4', '5', '6', '7', '8', '9', '10'])
axs[1].set_xlabel('Radius (mm)')
axs[1].set_ylabel('-log10(p-value)')
axs[1].set_title('Cortical depth 0-10')
axs[1].legend()

plt.tight_layout()

# Save the figures

if save_figures:
    file_name = f"Figure_6_{layer}.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=600)
    
plt.show()

#%% Save clean version of figures

if save_figures:
    # Create subplots side by side
    fig, axs = plt.subplots(1, 2, figsize=(12, 6), sharex=True)
    
    # First subplot: Mean Absolute R
    axs[0].plot(x, mean_abs_r[:, 1], label='H1', color='black', linewidth=2, marker='o')
    axs[0].fill_between(x, mean_abs_r[:, 1] - std_abs_r[:, 1], mean_abs_r[:, 1] + std_abs_r[:, 1], color='black', alpha=0.2)
    axs[0].plot(x, mean_abs_r[:, 0], label='H0 shifted', color='blue', linewidth=2, marker='o')
    axs[0].fill_between(x, mean_abs_r[:, 0] - std_abs_r[:, 0], mean_abs_r[:, 0] + std_abs_r[:, 0], color='blue', alpha=0.2)
    axs[0].plot(x, mean_abs_r_rotated[:, 0], label='H0 rotated', color='red', linewidth=2, marker='o')
    axs[0].fill_between(x, mean_abs_r_rotated[:, 0] - std_abs_r_rotated[:, 0], mean_abs_r_rotated[:, 0] + std_abs_r_rotated[:, 0], color='red', alpha=0.2)
    axs[0].set_ylim(0.06, 0.17)
    axs[0].set_xticks(x[::2])
    axs[0].set_xticklabels([])
    axs[0].set_yticklabels([])

    # Plot p-values on -log10 scale
    axs[1].plot(x, neg_log_p_values, color='blue', linewidth=2, marker='o')
    axs[1].plot(x, neg_log_p_values_rot, color='red', linewidth=2, marker='o')
    axs[1].axhline(y=-np.log10(0.05), color='black', linestyle='--', linewidth=1.5)  # Dotted black line at -log10(0.05)

    # Add a corresponding y-tick and label at p = 0.05
    axs[1].set_yticks([1.3, 2, 3, 4, 5, 6])
    axs[1].set_yticklabels([])
    axs[1].set_xticks(x[::2])
    axs[1].set_xticklabels([])

    # Increase linewidth of axes
    for ax in axs:
        ax.spines['bottom'].set_linewidth(2.5)
        ax.spines['left'].set_linewidth(2.5)
        ax.spines['top'].set_linewidth(2.5)
        ax.spines['right'].set_linewidth(2.5)

        # Increase linewidth of ticks
        ax.tick_params(axis='both', which='major', width=2.5, length=6)
        ax.tick_params(axis='both', which='minor', width=2.5, length=3)

    plt.tight_layout()
    
    # save figures
    file_name = f"Figure_6_{layer}_clean.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=400)
