#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 02:41:00 2025

@author: ms1454
"""

import os
from datetime import datetime
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import sem
import warnings
warnings.filterwarnings("ignore")

#%% specifications

test = '1st_2nd_sess_sep_layers_concat_runs'

save_figures = 1
if save_figures:
    date_stamp = datetime.now().strftime('%Y-%m-%d')
    save_dir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Figure_7/{date_stamp}'
    os.makedirs(save_dir, exist_ok=True)

subjects = ['chss', 'aroo', 'aman', 'auil']
ROI = 'V1'
hemis = ['lh', 'rh']
layers = ['0-2', '4-6', '8-10']
layers_cd = ['Deep', 'Middle', 'Superficial']
classifiers = ['fitrsvm_1', 'fitrsvm_2', 'fitrnet_1', 'fitrnet_2', 'fitrnet_3', 'fitrnet_4',
               'fitrnet_5', 'fitrtree_1', 'fitrtree_2', 'fitrtree_3', 'fitrensemble_1', 'fitrensemble_2']

baseDir = f'/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Prediction/{test}'

#%% load data

all_data = np.zeros((len(subjects), len(layers), len(classifiers), 2, 2))

for layer_index, layer in enumerate(layers):
    
    for sub_index, subject in enumerate(subjects):
        
        data_path = os.path.join(baseDir, f'Regression_test_{subject}_layers_{layer}.mat')
        
        all_data[sub_index, layer_index, :, :, :] = sio.loadmat(data_path)['acc'][6:19,:,:]

#%% plot test data last 12 classifiers

# shape all_data = 5 subjects, 3 layers, 12 classifiers, 2 (test, training acc), 2 (lh, rh hemi)

# Plot settings
sns.set(style="whitegrid")
fig1, axes = plt.subplots(2, 3, figsize=(18, 10), sharey=False)
colors = sns.color_palette("rainbow", 12)  # 18 unique colors for classifiers

for h, hemi in enumerate(hemis):
    for l, layer in enumerate(layers):
        ax = axes[h, l]

        # Extract data for the given hemisphere and layer
        data = all_data[:, l, :, :, h]  # Shape: (subjects, classifiers, acc types)

        # Compute means and SEM
        mean_acc = np.mean(data, axis=0)  # Shape: (18 classifiers, 2 acc types)
        sem_acc = sem(data, axis=0)      # Shape: (18 classifiers, 2 acc types)

        # Classifier indices
        x = np.arange(12)

        # Bar width
        width = 0.8

        for i in range(12):
            # Plot bars with SEM error bars
            ax.errorbar(x[i] - width / 2, mean_acc[i, 0], yerr=sem_acc[i, 0], fmt='none', 
                         color='black', capsize=5, zorder=0)
            ax.bar(x[i] - width / 2, mean_acc[i, 0], width, label="Test", color=colors[i], alpha=1, capsize=5, zorder=1,  edgecolor='black', linewidth=1.5)
            
            # Plot individual subject values as dots
            for subj_idx in range(data.shape[0]):  # Loop through subjects
                ax.scatter(x[i] - width / 2, data[subj_idx, i, 0], color=colors[i], edgecolors='black', s=50, zorder=2)

        ax.grid(True)
        
        # Labels and titles
        ax.set_xticks(x)
        ax.set_xticklabels([f"{classifiers[i]}" for i in x], rotation=45, ha="right", fontsize = 16)
        ax.set_title(f"{layers_cd[l]} cortical depth - {hemi.upper()}", fontsize = 18)
        ax.set_ylim(bottom=0, top=0.6)
        ax.set_ylabel("Correlation coefficient Pearson r", fontsize = 16)
        ax.tick_params(axis='y', labelsize=16)

# Adjust layout
plt.tight_layout()

# Save the  plot
if save_figures:
    file_name = f"Figure_7_{test}_test.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=50)
    
plt.show()

#%% plot train data of last 12 classifiers

# shape all_data = 5 subjects, 3 layers, 18 classifiers, 2 (test, training acc), 2 (lh, rh hemi)

# Plot settings
sns.set(style="whitegrid")
fig2, axes = plt.subplots(2, 3, figsize=(18, 10), sharey=False)
colors = sns.color_palette("rainbow", 12)  # 18 unique colors for classifiers

for h, hemi in enumerate(hemis):
    for l, layer in enumerate(layers):
        ax = axes[h, l]

        # Extract data for the given hemisphere and layer
        data = all_data[:, l, :, :, h]  # Shape: (subjects, classifiers, acc types)

        # Compute means and SEM
        mean_acc = np.mean(data, axis=0)  # Shape: (18 classifiers, 2 acc types)
        sem_acc = sem(data, axis=0)      # Shape: (18 classifiers, 2 acc types)

        # Classifier indices
        x = np.arange(12)

        # Bar width
        width = 0.8

        for i in range(12):
            # Plot bars with SEM error bars
            ax.errorbar(x[i] - width / 2, mean_acc[i, 1], yerr=sem_acc[i, 1], fmt='none', 
                         color='black', capsize=5, zorder=0)
            ax.bar(x[i] - width / 2, mean_acc[i, 1], width, label="Test", color=colors[i], alpha=1, capsize=5, zorder=1,  edgecolor='black', linewidth=1.5)
            
            # Plot individual subject values as dots
            for subj_idx in range(data.shape[0]):  # Loop through subjects
                ax.scatter(x[i] - width / 2, data[subj_idx, i, 1], color=colors[i], edgecolors='black', s=50, zorder=2)

        ax.grid(True)
        
        # Labels and titles
        ax.set_xticks(x)
        ax.set_xticklabels([f"{classifiers[i]}" for i in x], rotation=45, ha="right", fontsize = 16)
        ax.set_title(f"{layers_cd[l]} cortical depth - {hemi.upper()}", fontsize = 18)
        ax.set_ylim(bottom=0, top=1)
        ax.set_ylabel("Training accuracy", fontsize = 16)
        ax.tick_params(axis='y', labelsize=16)

# Adjust layout
plt.tight_layout()

# Save the  plot
if save_figures:
    file_name = f"Figure_7_{test}_train.tiff"
    file_path = os.path.join(save_dir, file_name)
    plt.savefig(file_path, format='tiff', dpi=50)
    
plt.show()

#%% cross-validation

# shape all_data = 5 subjects, 3 layers, 18 classifiers, 2 (test, training acc), 2 (lh, rh hemi)

# Initialize data storage
N = len(subjects)
L = len(layers)
C = len(classifiers)

all_data = 0.5 * np.log((1 + all_data) / (1 - all_data))    
        
# Cross-validation results storage
selected_methods = []
test_results = np.zeros((N, L, 2))  # Store best test results per subject

for test_idx, test_subject in enumerate(subjects):
    methods_wins = np.zeros(C)

    # check which method wins most often for all 6 tests: 3 layers x 2 hemis
    for layer_idx in range(L):
        for hemi_idx in range(2):  # Loop over left and right hemispheres
            # Exclude the test subject
            train_data = np.delete(all_data[:, layer_idx, :, 0, hemi_idx], test_idx, axis=0)  # (N-1, C)

            # Compute mean performance across N-1 subjects
            mean_r = np.mean(train_data, axis=0)  # (C,)
            
            # Select the best method
            best_method_idx = np.argmax(mean_r)
            best_method = classifiers[best_method_idx]
            methods_wins[best_method_idx] += 1

            # Store the test subject's performance with the selected method
            test_results[test_idx, layer_idx, hemi_idx] = all_data[test_idx, layer_idx, best_method_idx, 0, hemi_idx]

    # Store which methods won for this subject across the 6 tests
    selected_methods.append((test_subject, classifiers[np.argmax(methods_wins)]))

# Count how often each method was chosen
method_counts = {classifier: 0 for classifier in classifiers}
for _, method in selected_methods:
    method_counts[method] += 1

# Print summary
print("\nCross-Validation Summary:")
for method, count in sorted(method_counts.items(), key=lambda x: -x[1]):
    print(f"{method}: selected {count}/{N} times")

# Test results for each subject using the best method from LOSO-CV
print("\nTest Subject Performance using Selected Methods:")
for i, (subject, method) in enumerate(selected_methods):
    print(f"{subject}: Best method: {method}, Test Results: {test_results[i, :, :]}")
