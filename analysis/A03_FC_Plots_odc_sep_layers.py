#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 09:26:04 2025

@author: ms1454
"""

import os
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import matplotlib.cm as cm
import warnings
warnings.filterwarnings("ignore")

#%% specify input

subsample_iter = ''

save_figures = 1
if save_figures:
    date_stamp = datetime.now().strftime('%Y-%m-%d')
    save_dir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Control_Analyses/odc_sep_layers/{date_stamp}'
    os.makedirs(save_dir, exist_ok=True)

subjectFolders = ['haas', 'chss', 'aroo', 'aman', 'ylri', 'rcgr', 'atib', 'evad', 'arak', 'imyy', 'auil']
ROI = 'V1'
baseDir = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results'

layers = [f'{ROI}_layers_0-2_intrahemispheric',
          f'{ROI}_layers_4-6_intrahemispheric',
          f'{ROI}_layers_8-10_intrahemispheric']

num_quantiles = 10 # number of distance quantiles

all_beta = 1 # use avergae over specific beta quantiles or average over all
if all_beta != 1:
    beta_comp = [[[1,2,3],[1,2,3]],[[4,5,6,7],[4,5,6,7]],[[8,9,10],[8,9,10]]]

fisher_transform = 1 # r values or z values

normalization = 'not_weighted' # weighted (alike + unalike), not_weighted (overall FC), none
eye = 0# 0 = both, 1 = eye1 (negative beta, left eye), 2 = eye2 (positive beta, right eye)
ylims = [0.1, 0.3]
ylims_sel = [0, 12]
#%% load data

if all_beta != 1:
    alike_values = np.zeros((len(subjectFolders), len(beta_comp), len(layers), num_quantiles))
    unalike_values = np.zeros((len(subjectFolders), len(beta_comp), len(layers), num_quantiles))
    
    for beta_index, beta in enumerate(beta_comp):
        for version_index, version in enumerate(layers):
            dataDir = os.path.join(baseDir, version)
            
            for sub_index, subName in enumerate(subjectFolders):
                subDir = os.path.join(dataDir, subName)
                data_mat_name = 'CorrelationMtx_Selectivity_subsampled.mat'
                data_mat_path = os.path.join(subDir, data_mat_name)
                data_mat = sio.loadmat(data_mat_path)
                Data_Combined = np.squeeze(data_mat['Data_Combined'])
                
                if fisher_transform == 1:
                    Data_Combined = Data_Combined[:,:,:,:,:,1]
                else:
                    Data_Combined = Data_Combined[:,:,:,:,:,0]
    
                Data_Combined = np.mean(np.mean(Data_Combined[:,:,:, beta[0][0]-1:beta[0][-1],beta[1][0]-1:beta[1][-1]],axis=3),axis=3)
                Data_Combined = Data_Combined[:,[eye,3],:]
                
                Data_Combined_avg = np.mean(Data_Combined, axis=2) # Average across hemispheres (axis 2)

                alike_values[sub_index, beta_index, version_index, :] = Data_Combined_avg[:, 0]
                unalike_values[sub_index, beta_index, version_index, :] = Data_Combined_avg[:, 1]
                
    alike_values = np.mean(alike_values, axis=1) # average over beta comps
    unalike_values = np.mean(unalike_values, axis=1) #average over beta_comps

else: 
    
    alike_values = np.zeros((len(subjectFolders), len(layers), num_quantiles))
    unalike_values = np.zeros((len(subjectFolders), len(layers), num_quantiles))
    diff_values = np.zeros((len(subjectFolders), len(layers), num_quantiles))
    
    FC_values = np.zeros((len(subjectFolders), len(layers), num_quantiles))
    
    for version_index, version in enumerate(layers):
        dataDir = os.path.join(baseDir, version)
        
        
        for sub_index, subName in enumerate(subjectFolders):
            subDir = os.path.join(dataDir, subName)
            
            # load overall FC
            data_mat_name = 'CorrelationMtx_FC_subsampled.mat'
            data_mat_path = os.path.join(subDir, data_mat_name)
            
            if not os.path.exists(data_mat_path):
                print(f"File not found: {data_mat_path}")
                continue
            
            data_mat = sio.loadmat(data_mat_path)
            Data_Combined = np.squeeze(data_mat['Data_Combined'])
            
            # Apply Fisher Transform or not
            Data_Combined = Data_Combined[:, :, :, :, 1] if fisher_transform else Data_Combined[:, :, :, :, 0]
            
            # Average across hemispheres (axis=1) and distance quantiles (axis=0)
            Data_Combined_avg = np.mean(np.mean(np.mean(Data_Combined, axis=4), axis=3), axis=2)
            
            FC_values[sub_index, version_index :] = Data_Combined_avg[:,0]
            
            data_mat_name = 'CorrelationMtx_Selectivity_subsampled.mat'
            data_mat_path = os.path.join(subDir, data_mat_name)
            data_mat = sio.loadmat(data_mat_path)
            Data_Combined = np.squeeze(data_mat['Data_Combined'])
            
            if fisher_transform == 1:
                Data_Combined = Data_Combined[:,:,:,:,:,1]
            else:
                Data_Combined = Data_Combined[:,:,:,:,:,0]

            Data_Combined = np.mean(np.mean(Data_Combined,axis=3),axis=3)
            Data_Combined = Data_Combined[:,[eye,3],:]
            
            Data_Combined = np.mean(Data_Combined, axis=2) # Average across hemispheres (axis 2)

            alike_values[sub_index, version_index, :] = Data_Combined[:, 0]
            unalike_values[sub_index, version_index :] = Data_Combined[:, 1]
            
            if normalization == 'weighted':
                Data_Combined_norm = Data_Combined.copy()

                Data_Combined_alike_unalike = (Data_Combined_norm[:, 0]) + (Data_Combined_norm[:, 1])
                Data_Combined_norm[:, 0] = (Data_Combined_norm[:, 0] / Data_Combined_alike_unalike) * 200
                Data_Combined_norm[:, 1] = (Data_Combined_norm[:, 1] / Data_Combined_alike_unalike) * 200   
                
                diff_values[sub_index, version_index, :] = Data_Combined_norm[:, 0] - Data_Combined_norm[:,1]        

            elif normalization == 'not_weighted':
                
                Data_Combined_norm = Data_Combined.copy()

                Data_Combined_norm[:, 0] = (Data_Combined_norm[:, 0] / FC_values[sub_index,version_index,:]) * 200
                Data_Combined_norm[:, 1] = (Data_Combined_norm[:, 1] / FC_values[sub_index,version_index,:]) * 200   
                
                diff_values[sub_index, version_index, :] = Data_Combined_norm[:, 0] - Data_Combined_norm[:,1] 
                
            else:

                diff_values[sub_index, :,:] = alike_values[sub_index, :,:] - unalike_values[sub_index, :,:]
                    
            
# mean and ste over subjects
avg_values_alike = np.mean(np.mean(alike_values, axis=0),axis=0)
avg_values_unalike = np.mean(np.mean(unalike_values, axis=0),axis=0)
std_values_alike = np.std(np.mean(alike_values,axis=1), axis=0) / np.sqrt(len(subjectFolders))
std_values_unalike = np.std(np.mean(unalike_values,axis=1), axis=0) / np.sqrt(len(subjectFolders))

avg_values_diff = np.mean(np.mean(diff_values, axis=0),axis=0)
std_values_diff = np.std(np.mean(diff_values,axis=1), axis=0) / np.sqrt(len(subjectFolders))

#%%

# Summary statistics
print("Mean of selectivity values:", np.mean(avg_values_diff))
print("SEM of selectivity values:", np.std(avg_values_diff) / np.sqrt(10))

#%%

# Define accessible colors
cud_colors = ['#808080']

bar_width = 1  # Width of the bars
cap_width = 9  # Cap size of error bars
x = np.arange(avg_values_alike.shape[0]) * 2.7  # Adjust spacing for alike/unalike plot
x_diff = np.arange(avg_values_diff.shape[0]) * 2  # Adjust spacing for selectivity plot

# Create figure with different widths for left and right subplots
fig = plt.figure(figsize=(16, 6), constrained_layout=True)
gs = fig.add_gridspec(1, 2, width_ratios=[1.9, 1.1])  # Left plot wider, right plot narrower

# Alike and Unalike Subplot (left)
ax1 = fig.add_subplot(gs[0])
positions_mid_alike = x - 0.5 * bar_width
positions_mid_unalike = x + 0.5 * bar_width

ax1.errorbar(positions_mid_alike, avg_values_alike, yerr=std_values_alike, fmt='none', 
             color='black', capsize=cap_width, elinewidth=2, zorder=0)
deep_alike_bar = ax1.bar(positions_mid_alike, avg_values_alike, bar_width, 
                         label='Alike', color='white', edgecolor='black', linewidth=1.5)

ax1.errorbar(positions_mid_unalike, avg_values_unalike, yerr=std_values_unalike, fmt='none', 
             color='black', capsize=cap_width, elinewidth=2, zorder=0)
deep_unalike_bar = ax1.bar(positions_mid_unalike, avg_values_unalike, bar_width, 
                           label='Unalike', color='black', edgecolor='black', linewidth=1.5)

mid_x_ticks = (positions_mid_alike + positions_mid_unalike) / 2
ax1.set_xticks(mid_x_ticks)
ax1.set_xticklabels(range(1, 11))
ax1.set_xlabel('Distance')
ax1.set_ylim(ylims)
ax1.set_xlim(left=positions_mid_alike.min() - 1.2, right=positions_mid_unalike.max() + 1.2)
ax1.set_ylabel('z-transformed Pearson correlation coefficient' if fisher_transform else 'Pearson correlation coefficient')
ax1.set_title('Effect of cortical depth and distance on rs-FC strength (Alike vs Unalike)')
ax1.legend(handles=[deep_alike_bar, deep_unalike_bar])

# Set y-ticks for the left plot: Bottom, Middle, and Top
ax1.set_yticks([ylims[0], ylims[0]+(ylims[1]-ylims[0])/2, ylims[1]])

# Selectivity Subplot (right)
ax2 = fig.add_subplot(gs[1])
x_diff = np.arange(avg_values_diff.shape[0]) * 1.6
positions_deep_alike = x_diff - 1 * bar_width

ax2.errorbar(positions_deep_alike, avg_values_diff, yerr=std_values_diff, fmt='none', 
             color='black', capsize=cap_width, elinewidth=2, zorder=0)
deep_diff_bar = ax2.bar(positions_deep_alike, avg_values_diff, bar_width, 
                        label='Selectivity', color=cud_colors[0], edgecolor='black', linewidth=1.5)

ax2.set_xticks(positions_deep_alike)
ax2.set_xticklabels(range(1, 11))
ax2.set_xlabel('Distance')
ax2.set_ylim(ylims_sel)
ylabel_text = ('% change of z-transformed Pearson r (weighted normalization)' if fisher_transform and normalization == 'weighted' else
               '% change of z-transformed Pearson r (normalization)' if fisher_transform and normalization == 'not_weighted' else
               '% change of Pearson r (weighted normalization)' if normalization == 'weighted' else
               '% change of Pearson r (normalization)' if normalization == 'not_weighted' else
               'Difference of z-transformed Pearson r' if fisher_transform else
               'Difference of Pearson r')
ax2.set_ylabel(ylabel_text)
ax2.set_title('Effect of cortical depth and distance on selectivity')
ax2.legend(handles=[deep_diff_bar])

# Set y-ticks for the right plot: Bottom, Middle, and Top
ax2.set_yticks([ylims_sel[0], ylims_sel[0]+(ylims_sel[1]-ylims_sel[0])/2, ylims_sel[1]])

for ax in [ax1, ax2]:
    ax.spines['bottom'].set_linewidth(2)  # Increase x-axis linewidth
    ax.spines['left'].set_linewidth(2)    # Increase y-axis linewidth
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.tick_params(axis='both', labelsize=16, width=2, color='black')

plt.tight_layout()

if save_figures:
    file_name = f"Figure.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=200)
    
# Show the final plot
plt.show()

# Create figure with different widths for left and right subplots
if save_figures:
    fig = plt.figure(figsize=(16, 6), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.9, 1.1])  # Left plot wider, right plot narrower
    
    # Alike and Unalike Subplot (left)
    ax1 = fig.add_subplot(gs[0])
    positions_mid_alike = x - 0.5 * bar_width
    positions_mid_unalike = x + 0.5 * bar_width
    
    ax1.errorbar(positions_mid_alike, avg_values_alike, yerr=std_values_alike, fmt='none', 
                 color='black', capsize=cap_width, capthick=3, elinewidth=3, zorder=0)
    deep_alike_bar = ax1.bar(positions_mid_alike, avg_values_alike, bar_width, 
                             label='Alike', color='white', edgecolor='black', linewidth=2)
    
    ax1.errorbar(positions_mid_unalike, avg_values_unalike, yerr=std_values_unalike, fmt='none', 
                 color='black', capsize=cap_width, capthick=3, elinewidth=3, zorder=0)
    deep_unalike_bar = ax1.bar(positions_mid_unalike, avg_values_unalike, bar_width, 
                               label='Unalike', color='black', edgecolor='black', linewidth=2)
    
    mid_x_ticks = (positions_mid_alike + positions_mid_unalike) / 2
    ax1.set_xticks(mid_x_ticks)
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.set_ylim(ylims)
    ax1.set_xlim(left=positions_mid_alike.min() - 1.2, right=positions_mid_unalike.max() + 1.2)
    
    # Set y-ticks for the left plot: Bottom, Middle, and Top
    ax1.set_yticks([ylims[0], ylims[0]+(ylims[1]-ylims[0])/2, ylims[1]])
    
    # Selectivity Subplot (right)
    ax2 = fig.add_subplot(gs[1])
    x_diff = np.arange(avg_values_diff.shape[0]) * 1.6
    positions_deep_alike = x_diff - 1 * bar_width
    
    ax2.errorbar(positions_deep_alike, avg_values_diff, yerr=std_values_diff, 
    fmt='none', ecolor='black', elinewidth=3, capsize=cap_width, capthick=3,zorder=0)
    deep_diff_bar = ax2.bar(positions_deep_alike, avg_values_diff, bar_width, color=cud_colors[0], ecolor='black', edgecolor='black', linewidth=2)
    
    ax2.set_xticks(positions_deep_alike)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.set_ylim([0, 12])
    
    # Set y-ticks for the right plot: Bottom, Middle, and Top
    ax2.set_yticks([ylims_sel[0], ylims_sel[0]+(ylims_sel[1]-ylims_sel[0])/2, ylims_sel[1]])
    
    for ax in [ax1, ax2]:
        ax.spines['bottom'].set_linewidth(3)  # Increase x-axis linewidth
        ax.spines['left'].set_linewidth(3)    # Increase y-axis linewidth
        ax.spines['bottom'].set_color('black')
        ax.spines['left'].set_color('black')
        ax.tick_params(axis='both', labelsize=16, width=2, color='black')
        
    plt.tight_layout()
    file_name = f"Figure_2_clean{subsample_iter}.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=200)

#%%

# Define accessible colors
cud_colors = ['#808080']
bar_width = 1  # Width of the bars
cap_width = 9  # Cap size of error bars
x_diff = np.arange(avg_values_diff.shape[0]) * 2  # Adjust spacing for selectivity plot

# Generate distinct rainbow colors for each subject
num_subjects = len(subjectFolders)
colors = cm.rainbow(np.linspace(0, 1, num_subjects))

# Create standalone figure for selectivity
fig, ax2 = plt.subplots(figsize=(8, 6))
positions_mid_diff = x_diff  # Centered positions for bars

# Plot bars and error bars
ax2.errorbar(positions_mid_diff, avg_values_diff, yerr=std_values_diff, 
             fmt='none', color='black', capsize=cap_width, elinewidth=2, zorder=0)
diff_bar = ax2.bar(positions_mid_diff, avg_values_diff, width=bar_width, color=cud_colors[0], zorder=1)

# Scatter plot for individual subjects
for sub_idx in range(num_subjects):
    ax2.scatter(positions_mid_diff, diff_values[sub_idx, 0, :], 
                color=colors[sub_idx], s=50, label=subjectFolders[sub_idx], edgecolors='black', alpha=0.5, zorder=2)
    ax2.plot(positions_mid_diff, diff_values[sub_idx, 0, :], color=colors[sub_idx], alpha=0.5, linestyle='-', linewidth=2, zorder=1)
# Formatting
ax2.set_xticks(positions_mid_diff)
ax2.set_xticklabels([f'{i+1}' for i in range(len(positions_mid_diff))], fontsize=12)
ax2.set_ylabel('Rs-FC selectivity (%)', fontsize=14)
ax2.set_xlabel('Distance quantile', fontsize=14)
ax2.set_ylim(-ylims_sel[1], ylims_sel[1])
ax2.grid(axis='y', linestyle='--', alpha=0.6)

# Add legend
ax2.legend(loc='upper right', bbox_to_anchor=(1.3, 1), prop={'size': 14})  # Increase size to 14 (adjust as needed)

plt.tight_layout()

if save_figures:
    file_name = f"Figure_2{subsample_iter}_Suppl.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=200)
    
plt.show()