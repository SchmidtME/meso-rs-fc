#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is plotting violin plots to illustrate rs-FC as a function of
distance and ocular polarity and corresponds to Figure 2 of the manuscript.

@author: Marianna Elisa Schmidt (marianna.schmidt@maxplanckschools.de)
"""

import os
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.io as sio
import warnings
warnings.filterwarnings("ignore")

#%% specify input

subsample_iter = ''

save_figures = 1
if save_figures:
    date_stamp = datetime.now().strftime('%Y-%m-%d')
    save_dir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Figure_2/{date_stamp}'
    os.makedirs(save_dir, exist_ok=True)
    
subjectFolders = ['myla', 'haas', 'chss', 'aroo', 'aman', 'ylri', 'rcgr', 'atib', 'evad', 'arak', 'imyy', 'auil']
ROI = 'V1'
baseDir = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results'

layers = [f'{ROI}_layers_0-10_intrahemispheric']

num_quantiles = 10 # number of distance quantiles

fisher_transform = 1 # r values or z values

normalization = 'not_weighted' # weighted (alike + unalike), not_weighted (overall FC), none
eye = 0# 0 = both, 1 = eye1 (negative beta, left eye), 2 = eye2 (positive beta, right eye)
ylims = [0, 0.7] # 0.14
ylims_sel = [-25,60]

ylabel_text = ('% change of z-transformed Pearson r (weighted normalization)' if fisher_transform and normalization == 'weighted' else
               '% change of z-transformed Pearson r (normalization)' if fisher_transform and normalization == 'not_weighted' else
               '% change of Pearson r (weighted normalization)' if normalization == 'weighted' else
               '% change of Pearson r (normalization)' if normalization == 'not_weighted' else
               'Difference of z-transformed Pearson r' if fisher_transform else
               'Difference of Pearson r')

#%% load data


alike_values = np.zeros((len(subjectFolders), len(layers), num_quantiles))
unalike_values = np.zeros((len(subjectFolders), len(layers), num_quantiles))
diff_values = np.zeros((len(subjectFolders), len(layers), num_quantiles))

FC_values = np.zeros((len(subjectFolders), len(layers), num_quantiles))

for version_index, version in enumerate(layers):
    dataDir = os.path.join(baseDir, version)
    
    
    for sub_index, subName in enumerate(subjectFolders):
        subDir = os.path.join(dataDir, subName)
        
        # load overall FC
        data_mat_name = 'CorrelationMtx_FC_subsampled_1.mat'
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
        
        data_mat_name = 'CorrelationMtx_Selectivity_subsampled_1.mat'
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

            Data_Combined_norm[:, 0] = (Data_Combined_norm[:, 0] / FC_values[sub_index,:,:]) * 200
            Data_Combined_norm[:, 1] = (Data_Combined_norm[:, 1] / FC_values[sub_index,:,:]) * 200   
            
            diff_values[sub_index, version_index, :] = Data_Combined_norm[:, 0] - Data_Combined_norm[:,1] 
            
        else:

            diff_values[sub_index, :,:] = alike_values[sub_index, :,:] - unalike_values[sub_index, :,:]
                    
            
# mean and ste over subjects
avg_values_alike = np.mean(alike_values, axis=0)
avg_values_unalike = np.mean(unalike_values, axis=0)
std_values_alike = np.std(alike_values, axis=0) / np.sqrt(len(subjectFolders))
std_values_unalike = np.std(unalike_values, axis=0) / np.sqrt(len(subjectFolders))

avg_values_diff = np.mean(diff_values, axis=0)
std_values_diff = np.std(diff_values, axis=0) / np.sqrt(len(subjectFolders))

#%%

# Summary statistics
print("Mean of selectivity values:", np.mean(avg_values_diff[0, :]))
print("SEM of selectivity values:", np.std(avg_values_diff[0, :]) / np.sqrt(10))

#%%

# Prepare dataframe for alike/unalike
data_list = []
for q in range(num_quantiles):
    for subj_val in alike_values[:, 0, q]:
        data_list.append([subj_val, q+1, 'Alike'])
    for subj_val in unalike_values[:, 0, q]:
        data_list.append([subj_val, q+1, 'Unalike'])

df_alike_unalike = pd.DataFrame(data_list, columns=['Value', 'Quantile', 'Condition'])

# Prepare dataframe for selectivity
data_list_sel = []
for q in range(num_quantiles):
    for subj_val in diff_values[:, 0, q]:
        data_list_sel.append([subj_val, q+1])
df_selectivity = pd.DataFrame(data_list_sel, columns=['Value', 'Quantile'])

#%%

# Brighter, CUD-friendly colors with higher transparency (alpha=0.4)
alike_color   = (1.0, 0.8, 0.3, 0.4)  # yellowish-orange
unalike_color = (0.6, 0.5, 0.9, 0.4)  # bluish-purple

# Blend for difference (average RGB, alpha=0.4)
diff_color = (
    (alike_color[0] + unalike_color[0]) / 2,
    (alike_color[1] + unalike_color[1]) / 2,
    (alike_color[2] + unalike_color[2]) / 2,
    0.4
)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# --- Left: split violin ---
sns.violinplot(
    data=df_alike_unalike,
    x="Quantile", y="Value", hue="Condition", cut=0,
    split=True, inner="quartile", linewidth=1.2, linecolor='black', gap=.2,
    palette={"Alike": alike_color, "Unalike": unalike_color},
    ax=axes[0]
)

# Calculate and plot means for Alike and Unalike
means_alike = df_alike_unalike[df_alike_unalike["Condition"] == "Alike"].groupby("Quantile")["Value"].mean()
means_unalike = df_alike_unalike[df_alike_unalike["Condition"] == "Unalike"].groupby("Quantile")["Value"].mean()

# Overlay means (left half = Alike, right half = Unalike)
for i, q in enumerate(sorted(df_alike_unalike["Quantile"].unique())):
    axes[0].plot(i - 0.15, means_alike[q], 'o', color='black', markersize=6, zorder=3)  # Alike
    axes[0].plot(i + 0.15, means_unalike[q], 'o', color='black', markersize=6, zorder=3)  # Unalike

axes[0].set_xlabel("Distance")
axes[0].set_ylabel("z-transformed Pearson r" if fisher_transform else "Pearson r")
axes[0].set_ylim(ylims)
axes[0].set_title("Alike vs Unalike (Split Violin)")
sns.move_legend(axes[0], "upper right", frameon=False)

# --- Right: selectivity violins ---
sns.violinplot(
    data=df_selectivity,
    x="Quantile", y="Value", cut=0,
    inner="quartile", linewidth=1.2, linecolor='black',
    color=diff_color, ax=axes[1]
)

# Calculate and plot means for selectivity
means_sel = df_selectivity.groupby("Quantile")["Value"].mean()
for i, q in enumerate(sorted(df_selectivity["Quantile"].unique())):
    axes[1].plot(i, means_sel[q], 'o', color='black', markersize=6, zorder=3)

axes[1].set_xlabel("Distance")
axes[1].set_ylabel(ylabel_text)
axes[1].set_ylim(ylims_sel)
axes[1].set_title("Selectivity")

# --- Styling ---
for ax in axes:
    ax.tick_params(axis='both', labelsize=14, width=2, color='black')
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)

plt.tight_layout()

if save_figures:
    file_name = f"Figure_2.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=400)

plt.show()

#%% Clean version with specific y-ticks (no tick labels) and thicker axes/ticks
fig_clean, axes_clean = plt.subplots(1, 2, figsize=(16, 6))

# --- Left: split violin ---
sns.violinplot(
    data=df_alike_unalike,
    x="Quantile", y="Value", hue="Condition", cut=0,
    split=True, inner="quartile", linewidth=1.2, linecolor='black', gap=.2,
    palette={"Alike": alike_color, "Unalike": unalike_color},
    ax=axes_clean[0]
)
for i, q in enumerate(sorted(df_alike_unalike["Quantile"].unique())):
    axes_clean[0].plot(i - 0.15, means_alike[q], 'o', color='black', markersize=6, zorder=3)
    axes_clean[0].plot(i + 0.15, means_unalike[q], 'o', color='black', markersize=6, zorder=3)

# --- Right: selectivity violins ---
sns.violinplot(
    data=df_selectivity,
    x="Quantile", y="Value", cut=0,
    inner="quartile", linewidth=1.2, linecolor='black',
    color=diff_color, ax=axes_clean[1]
)
for i, q in enumerate(sorted(df_selectivity["Quantile"].unique())):
    axes_clean[1].plot(i, means_sel[q], 'o', color='black', markersize=6, zorder=3)

# --- Enforce specific y-ticks (no labels) ---
axes_clean[0].set_yticks([0.00, 0.35, 0.70])    # left panel
axes_clean[1].set_yticks([-30, 15, 60])         # right panel

# --- Styling: thicker axes/ticks, no labels anywhere ---
for ax in axes_clean:
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")
    if ax.get_legend():
        ax.get_legend().remove()

    # Thicker ticks and axes; remove x- and y-tick labels
    ax.tick_params(axis='x', labelbottom=False, width=3.5, length=8, color='black')
    ax.tick_params(axis='y', labelleft=False,  width=3.5, length=8, color='black')

    ax.spines['bottom'].set_linewidth(3.5)
    ax.spines['left'].set_linewidth(3.5)

plt.tight_layout()

if save_figures:
    file_name = f"Figure_2_clean.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=400)
    
plt.show()


#%% Save DataFrames used for plotting
save_data_dir = "/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figure_2_plot_data"
os.makedirs(save_data_dir, exist_ok=True)

df_alike_unalike.to_csv(os.path.join(save_data_dir, "df_alike_unalike.csv"), index=False)
df_selectivity.to_csv(os.path.join(save_data_dir, "df_selectivity.csv"), index=False)
