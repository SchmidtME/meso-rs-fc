#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 03:17:34 2025

@author: ms1454
"""

import os
from datetime import datetime
import scipy.io as sio
import h5py
import numpy as np
import matplotlib.pyplot as plt

save_figures = 1
if save_figures:
    date_stamp = datetime.now().strftime('%Y-%m-%d')
    save_dir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Supplementary/Figure_2/{date_stamp}'
    os.makedirs(save_dir, exist_ok=True)


# Define the base path
base_path = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results/V1_layers_0-2_intrahemispheric'

# Get the list of subdirectories (folders) in the base path
subjects = ['haas',
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

# Initialize a list to store averaged results from each subject
mat_data = np.zeros((len(subjects),2,2,10,1))

# Loop through each folder
for s, subject in enumerate(subjects):
    folder_path = os.path.join(base_path, subject)
    mat_file_path = os.path.join(folder_path, 'CorrelationMtx_FC_beta_Params_subsampled.mat')

    # Load the .mat file
    data = h5py.File(mat_file_path, 'r')
    mat_data_eye1_lh = data['AnalysisParam']['betas_eye1_medians'][0,0]
    mat_data[s,0,0,:,:] = np.array(data[mat_data_eye1_lh])
    mat_data_eye1_rh = data['AnalysisParam']['betas_eye1_medians'][1,0]
    mat_data[s,0,1,:,:] = np.array(data[mat_data_eye1_rh])
    mat_data_eye2_lh = data['AnalysisParam']['betas_eye2_medians'][0,0]
    mat_data[s,1,0,:,:] = np.array(data[mat_data_eye2_lh])
    mat_data_eye2_rh = data['AnalysisParam']['betas_eye2_medians'][1,0]
    mat_data[s,1,1,:,:] = np.array(data[mat_data_eye2_rh])
    
mean_data = np.mean(np.mean(np.abs(mat_data), axis=1), axis=1)

mean_mean_data = np.mean(mean_data, axis=0)

#%% positive quantiles

# Create a figure for all subjects
plt.figure(figsize=(10, 6))

# Plot individual subject results
for result in mean_data:
    plt.plot(np.arange(1, 11), result, color='grey', alpha=0.5, linewidth=1)

# Plot mean across subjects
plt.plot(np.arange(1, 11), mean_mean_data[:, 0], color='b', linewidth=3, label='Mean across subjects')

for x, y in enumerate(mean_mean_data[:, 0], start=1):
    plt.scatter(x, y, color='b', s=50)  # s=50 sets the size of the dots
    plt.text(x - 0.08, y + 0.15, f'{y:.2f}', ha='center', va='bottom', fontsize=15)  # Move the text up and increase size

# Set title and labels
plt.title('Beta medians within quantiles across all subjects', fontsize=20)
plt.xlabel('Beta quantiles', fontsize=15)
plt.ylabel('Median beta value', fontsize=15)
plt.xticks(np.arange(1, 11), fontsize=15)
plt.yticks(fontsize=15)

# Increase linewidth of axes
ax = plt.gca()
ax.spines['bottom'].set_linewidth(2.5)
ax.spines['left'].set_linewidth(2.5)
ax.spines['top'].set_linewidth(2.5)
ax.spines['right'].set_linewidth(2.5)

# Increase linewidth of ticks
ax.tick_params(axis='both', which='major', width=2.5, length=6)
ax.tick_params(axis='both', which='minor', width=2.5, length=3)

#plt.grid(False, axis='y', linestyle='--', alpha=0.6)
plt.legend(fontsize=15)

# Save the plot
if save_figures:
    file_name = "Figure_2.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=200)

plt.show()

#%%

if save_figures:
    # Create a figure for all subjects
    plt.figure(figsize=(10, 6))
    
    # Plot individual subject results
    for result in mean_data:
        plt.plot(np.arange(1, 11), result, color='grey', alpha=0.5, linewidth=1)
    
    # Plot mean across subjects
    plt.plot(np.arange(1, 11), mean_mean_data[:, 0], color='b', linewidth=3, label='Mean across subjects')
    
    for x, y in enumerate(mean_mean_data[:, 0], start=1):
        plt.scatter(x, y, color='b', s=50)  # s=50 sets the size of the dots
        #plt.text(x -0.02, y + 0.15, f'{y:.2f}', ha='center', va='bottom', fontsize=15)  # Move the text up and increase size
    
    # Get current axis
    ax = plt.gca()
    
    # Explicitly set all 10 x-ticks
    ax.set_xticks(np.arange(1, 11))
    
    # Remove tick labels but keep tick marks
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    
    # Increase linewidth of axes
    ax.spines['bottom'].set_linewidth(2.5)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['top'].set_linewidth(2.5)
    ax.spines['right'].set_linewidth(2.5)
    
    # Increase linewidth of ticks
    ax.tick_params(axis='both', which='major', width=2.5, length=6)
    ax.tick_params(axis='both', which='minor', width=2.5, length=3)
    
    ax.set_ylim(0, 2.5)
    
    
    # Save the plot
    file_name = "Figure_2_clean.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=400)
    
    plt.show()
