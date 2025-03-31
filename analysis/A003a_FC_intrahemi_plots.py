#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 08:40:52 2025

@author: ms1454
"""


import os
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.io as sio
import warnings
warnings.filterwarnings("ignore")

#%% Input specifications

analysis_title= 'Intrahemispheric rs-FC'
   
ROIs = ['V1', 'V2', 'V3']

layers = ['0-2', '4-6', '8-10']

save_figures = 1
if save_figures:
    date_stamp = datetime.now().strftime('%Y-%m-%d')
    save_dir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Supplementary/Figure_4/{date_stamp}'
    os.makedirs(save_dir, exist_ok=True)

fisher_transform = 1 # r values or z values
    
ylims = [0, 0.25]

#%%

Groups = ['Control']
group_data = {}  # Dictionary to store data for each group

for Group in Groups:
    if Group == 'Control':
        Sbjs = ['haas', 'chss', 'aroo', 'aman', 'ylri', 'rcgr', 'atib', 'evad', 'arak', 'imyy', 'auil']
    
    baseDir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results'
    
    FC_values = np.zeros((len(Sbjs), len(ROIs), len(layers)))
    
    for ROI_index, ROI in enumerate(ROIs):
        
        for version_index, version in enumerate(layers):
            
            targetDir = f'{ROI}_layers_{version}_intrahemispheric'
            
            dataDir = os.path.join(baseDir, targetDir)
            
            for sub_index, subName in enumerate(Sbjs):
                
                subDir = os.path.join(dataDir, subName)

                data_mat = sio.loadmat(os.path.join(subDir, 'CorrelationMtx_justFC.mat'))
                Data_Combined = np.mean(np.squeeze(data_mat['Data_Combined']),axis=0)
                
                if fisher_transform == 1:
                    Data_Combined = Data_Combined[1]
                else:
                    Data_Combined = Data_Combined[0]
    
                FC_values[sub_index, ROI_index, version_index] = Data_Combined
    
        # Save data for the group
        group_data[Group] = {
            'FC_values': FC_values
        }

#%%

# Define grayscale shades for ROIs
roi_colors = ['black', 'gray', 'lightgray']

# Compute mean and SEM
means = np.zeros((len(ROIs), len(layers)))
sems = np.zeros((len(ROIs), len(layers)))
for roi_idx, ROI in enumerate(ROIs):
    FC_values = group_data['Control']['FC_values'][:, roi_idx, :]
    means[roi_idx, :] = np.nanmean(FC_values, axis=0)
    sems[roi_idx, :] = stats.sem(FC_values, axis=0, nan_policy='omit')

# Plot results
fig, ax = plt.subplots(figsize=(8, 5))
x = np.arange(len(layers))  # X positions for layers
width = 0.25  # Bar width

for roi_idx, ROI in enumerate(ROIs):
    ax.bar(x + roi_idx * width, means[roi_idx, :], width=width, yerr=sems[roi_idx, :],
           label=ROI, color=roi_colors[roi_idx], edgecolor='black', capsize=5, alpha=0.8)

ax.set_xticks(x + width)
ax.set_xticklabels(['Deep', 'Middle', 'Superficial'], fontsize=14)
ax.set_xlabel('Cortical Depth', fontsize=15)
ax.set_ylabel('Fisher z transformed r' if fisher_transform else 'Pearson r', fontsize=15)
ax.set_ylim(ylims)
ax.legend(loc='upper right', fontsize=13, frameon=False)
ax.set_title(analysis_title, fontsize=16, fontweight='bold')

plt.tight_layout()
if save_figures:
    file_name = f"Figure.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=100)
plt.show()

