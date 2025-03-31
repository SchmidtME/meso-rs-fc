clear all
close all
clc

% Script for Denoising fMRI Data Using NORDIC
%
% This script applies NORDIC denoising to fMRI data for multiple subjects and runs.
% The denoised data is saved in the specified output directory. It handles different
% subject-specific directory structures and performs the denoising with pre-defined 
% parameters in ARG (coded like this to avoid interference with parallel processing
% if enabled. The NORDIC algorithm is used to reduce noise and improve data quality.

% Add NORDIC software path
addpath("/autofs/space/ardebil_001/users/Others/Marianna/Software/NORDIC_Raw")

%% Define paths and variables

% Define paths for original and processed data
dir_orig_data_ups = '/autofs/space/ardebil_003/users/Shahin/HighRes_FC/FC_7T_Coronal/Sbjs_Upsampled/';
dir_orig_data = '/autofs/space/ardebil_003/users/Shahin/HighRes_FC/FC_7T_Coronal/Sbjs/';
dir_data = '/autofs/space/ardebil_001/users/Others/Marianna/FC_7T_Coronal/Data_Bay2/';

% Get list of subject names from the processed data directory
name_subs = arrayfun(@(d0) d0.name, dir([dir_data 'h*']), 'UniformOutput', false);

%% Run denoising

% Loop through each subject
for sub = 1:numel(name_subs)
    % Retrieve the list of runs for the current subject
    runs{sub, :} = arrayfun(@(d0) d0.name, dir(fullfile(dir_data, [name_subs{sub}], 'bold_FLEET', '0*')), 'UniformOutput', false);

    % Process each run for the current subject
    for run = 1:numel(runs{sub})
        % Define the path to the input magnitude image (subject-specific)
        if name_subs{sub} == "chsB"
            fn_mag = fullfile(dir_orig_data, "chss_2", 'bold_Close', runs{sub}{run}, 'f.nii');
        elseif name_subs{sub} == "haaB"
            fn_mag = fullfile('/autofs/space/ardebil_003/users/Shahin/HighRes_FC/FC_7T_Coronal/Sbjs_Adaptation_Upsampled/haas1/S01_P01/bold_Close_Upsampled2', runs{sub}{run}, 'f.nii.gz');
        elseif name_subs{sub} == "rcgr"
            fn_mag = fullfile(dir_orig_data, name_subs{sub}, 'rest', runs{sub}{run}, 'f.nii');
        else
            fn_mag = fullfile(dir_data, name_subs{sub}, 'bold_FLEET', runs{sub}{run}, 'f.nii');
        end

        % Define output filename and directory for denoised data
        fn_out = 'data_denoised10';
        ARG.magnitude_only = 1;  % Process only the magnitude image
        ARG.temporal_phase = 1;  % Enable temporal phase filtering
        ARG.phase_filter_width = 10;  % Width of the temporal phase filter
        ARG.DIROUT = fullfile(dir_data, name_subs{sub}, 'bold_FLEET', runs{sub}{run});

        % Perform NORDIC denoising
        NIFTI_NORDIC(fn_mag, fn_mag, fn_out, ARG);

        % Clear variables except essential ones for the next iteration
        clearvars -except dir_orig_data_ups dir_orig_data dir_data name_subs sub runs
    end
end
