clear all
close all
clc

% Description:
% This script is designed for preprocessing functional MRI (fMRI) data,
% including slice-timing correction, upsampling, FSFAST preprocessing,
% and layer-specific analysis. It generates intermediate and final outputs
% required for analyzing the data in the context of intracortical layers.

%% Set paths, directory, and file names
% Define paths to input/output data and supporting directories
DIR_FS = '_anat_upsample_B1justsub'; % FreeSurfer anatomical data directory
DIR_SESS = 'bold_Close_Upsampled2'; % Session-specific functional data directory
PROJECT_PATH = '/autofs/space/ardebil_001/users/Others/Marianna/FC_7T_Coronal/Data_Denoised'; % Project data path
SUBJECTS_DIR = '/autofs/space/ardebil_001/users/Shared/good_subjects_anat'; % Subjects directory
CODE_PATH = '/autofs/space/ardebil_001/users/Others/Marianna/FC_7T_Coronal/Code/preproc'; % Code directory

% Export SUBJECTS_DIR as an environment variable
system(sprintf('export SUBJECTS_DIR=%s', SUBJECTS_DIR));

%% Analysis settings
% Specify processing steps to be performed
stc = 0;        % Perform slice-timing correction
upsample = 0;   % Perform data upsampling
do_preproc_1 = 0; % Run FreeSurfer preprocessing commands
gen_layers = 0;  % Sample fMRI to 9 intermediate surfaces
vert_smooth = 0; % Smooth across surfaces

funcstem = 'fmcpr'; % Functional data stem name

%% Analysis parameters
% Define subject(s), smoothing parameters, hemispheres, and layer settings
name_subj = {'haas'}; % List of subjects
FWHM = 0;             % Full-width half-maximum for smoothing (mm)
hemispheres = {'rh', 'lh'}; % Hemispheres: right and left
layers = arrayfun(@(x) sprintf('%02d', x), 0:10, 'UniformOutput', false); % Layer indices
layer_avgs = {'0-2', '4-6', '8-10', '0-10'}; % Averaged layer groups
rad_sizes = [3, 3, 3, 10];  % vertical smoothing across how many layers?  Corresponding to layer_avgs

%% Parallel loop over subjects
for i = 1:length(name_subj)
    sub = name_subj{i};
    
    % Set functional directory for the subject
    F_SUBJECTS_PATH = fullfile(PROJECT_PATH, sub);
    
    % Retrieve the names of functional runs for the session
    runs_dir = fullfile(F_SUBJECTS_PATH, DIR_SESS);
    name_runs_struct = dir(runs_dir);
    name_runs_struct = name_runs_struct([name_runs_struct.isdir] & ~startsWith({name_runs_struct.name}, '.') & startsWith({name_runs_struct.name}, '0'));
    name_runs = {name_runs_struct.name};
    
    %% Preprocessing steps
    % Perform slice-timing correction
    if stc
        cd(PROJECT_PATH);
        system(sprintf('stc-sess -s %s -fsd %s -i data_denoised10 -o data_denoised10_stc -so even', sub, DIR_SESS));
    end
    
    % Perform data upsampling
    if upsample
        for j = 1:length(name_runs)
            run = name_runs{j};
            cd(fullfile(F_SUBJECTS_PATH, DIR_SESS, run));
            if stc
                system('mri_convert data_denoised10_stc.nii.gz f.nii.gz --upsample 2');
            else
                system('mri_convert data_denoised10.nii f.nii.gz --upsample 2');
            end
        end
    end
    
    % Run FreeSurfer preprocessing
    if do_preproc_1
        cd(PROJECT_PATH);
        system(sprintf('preproc-sess -s %s -fsd %s -surf self rhlh -fwhm %d -per-run -no-reg', sub, DIR_SESS, FWHM));
    end
    
    %% Layer-specific analysis
    if gen_layers
        % Generate layer-specific data
        cd(PROJECT_PATH);
        for h = 1:length(hemispheres)
            hem = hemispheres{h};
            for j = 1:length(name_runs)
                run = name_runs{j};
                for l = 1:length(layers)
                    layer = layers{l};
                    output_file = fullfile(F_SUBJECTS_PATH, sprintf('%s/%s/%s.%s.sm0.self.midgray.%s.mgz', DIR_SESS, run, hem, funcstem, layer));
                    if ~isfile(output_file)
                        system(sprintf('mri_vol2surf --interp trilin --mov %s/%s/%s/%s.nii.gz --reg %s/%s/%s/register.dof6.lta --surf midgray.B1justsub.%s --hemi %s --o %s', ...
                            sub, DIR_SESS, run, funcstem, sub, DIR_SESS, run, layer, hem, output_file));
                    end
                end
            end
        end
    end

    if vert_smooth
        % Perform intracortical smoothing for averaged layers
        cd(PROJECT_PATH);
        for h = 1:length(hemispheres)
            hem = hemispheres{h};
            for k = 1:length(layer_avgs)
                layer_avg = layer_avgs{k};
                for j = 1:length(name_runs)
                    run = name_runs{j};
                    system(sprintf('mris_smooth_intracortical --surf_dir %s/%s%s/surf/ --surf_name %s.midgray.B1justsub.''??'' --overlay_dir %s/%s/%s/ --overlay_name %s.%s.sm0.self.midgray.''??''.mgz --output_dir %s/%s/%s/ --output_name %s.%s.sm0.self.midgray.00.nb1_rad%s.mgz --tan-size 1 --rad-start %s --rad-size %d', ...
                        SUBJECTS_DIR, sub, DIR_FS, hem, sub, DIR_SESS, run, hem, funcstem, sub, DIR_SESS, run, hem, funcstem, layer_avg, strtok(layer_avg, '-'), rad_sizes(k)));
                end
            end
        end
    end
end
