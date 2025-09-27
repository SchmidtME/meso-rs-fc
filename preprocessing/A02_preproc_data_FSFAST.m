clear all
close all
clc

% Description:
% This script is designed for preprocessing functional MRI (fMRI) data,
% including NORDIC denoising, slice-timing correction, upsampling, 
% FSFAST preprocessing, and layer-specific analysis. It generates intermediate and final outputs
% required for analyzing the data in the context of intracortical layers.

% add paths to necessary software
addpath("/autofs/space/ardebil_001/users/Others/Marianna/Software/NORDIC_Raw")

%% specifications

% Set paths and directory and file names
DIR_FS = '_anat_upsample_B1justsub';
DIR_SESS = 'bold_Close_Upsampled2';
PROJECT_PATH = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Data_Denoised';
SUBJECTS_DIR = '/autofs/space/ardebil_001/users/Shared/good_subjects_anat';

system(sprintf('export SUBJECTS_DIR=/autofs/space/ardebil_001/users/Shared/good_subjects_anat'));
% Analysis steps
denoise = 0; % denoising
stc = 0; % slice timing correction
upsample = 0; % upsampling
do_preproc_1 = 0; % preprco-sess until registration
do_preproc_2 = 0; % updates registration after manual correction
do_preproc_3 = 1; % rest of preproc-sess
gen_layers = 1; % sample fMRI data to 9 intermediate surfaces
vert_smooth = 1; % columnar smoothing smoothing
gen_covariates = 1; % generate wm.dat
config_file = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/CFG_Files/wm.config';
% config file can be created by fcseed-config, here we use one of a previous analysis

funcstem = 'fmcpr';

% Analysis parameters
name_subj = {'oban'}%'myla'}%'main', 'uces'}%'muit', 'alga', 'atic', 'enam', 'etex', 'ntil', 'ndnj', 'nauc'};
amblyopia = 0;
if amblyopia
    dir_orig_data = '/space/zahedan/3/users/Shahin/Marianna/MPM/Subjects/Amblyopes/';
    dir_data = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Amblyopia/Data_Bay2_Denoised/';
else
    dir_orig_data = '/space/ardebil/3/users/Shahin/HighRes_FC/FC_7T_Coronal/Sbjs/';
    dir_data = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Data_Denoised/';
end

FWHM = 0; % smoothing that should be performed during the processing in mm
hemispheres = {'rh', 'lh'};
layers = arrayfun(@(x) sprintf('%02d', x), 0:10, 'UniformOutput', false);
layer_avgs = {'0-10', '0-2', '4-6', '8-10'}; % Averaged layer groups
rad_sizes = [10, 3, 3, 3];  % vertical smoothing across how many layers?  Corresponding to layer_avgs

%% Processing

% Parallel loop over subjects
for i = 1:length(name_subj)
    sub = name_subj{i};
    
    % Set functional directory
    F_SUBJECTS_PATH = fullfile(PROJECT_PATH, sub); 
    
    % Get the name of the runs
    name_runs = arrayfun(@(d0) d0.name, ...
        dir(fullfile(dir_orig_data, sub, 'bold_Close', '0*')), ...
        'UniformOutput', false);

    %% Denoising
    % NORDIC denoising based on magnitude data

    if denoise

        for run = 1:numel(name_runs)
            
            dir_out = fullfile(dir_data, [name_subj{i}], 'bold_Close_Upsampled2', name_runs{run});
            mkdir(dir_out);
    
            % convert input .mgh files to .nii and save at new data structure
            dir_in = fullfile(dir_orig_data, sub, 'bold_Close');
            fn_mag = fullfile(dir_in, name_runs{run}, dir(fullfile(dir_in, name_runs{run}, 'f.nii')).name); 
    
            fn_out = 'data_denoised10';
            
            % parameters for NORDIC denoising of functional magnitude data
            % for parallel processing it seems that it has to be specified within the loop
            ARG = struct();
            ARG.magnitude_only = 1;
            ARG.temporal_phase=1;
            ARG.phase_filter_width=10;
            ARG.DIROUT = [dir_out '/'];
            NIFTI_NORDIC(fn_mag, fn_mag, fn_out, ARG);
    
            % delete input file (conversion of .mgh file)
            delete(fn_mag);

        end

    end

    %% Slice timing correction

    if stc
        cd(PROJECT_PATH);
        system(sprintf('stc-sess -s %s -fsd %s -i data_denoised10 -o data_stc -so odd', sub, DIR_SESS));
    end
    
    %% Upsampling

    if upsample
        for j = 1:length(name_runs)
            run = name_runs{j};
            cd(fullfile(F_SUBJECTS_PATH, DIR_SESS, run));
            %if ~isfile('f.nii.gz');
                system('mri_convert data_stc.nii.gz f.nii.gz --upsample 2');
            %end
        end
    end
    
    %% FSFAST up to registration

    if do_preproc_1
        cd(PROJECT_PATH);

        % create subjectname txt file
        fileID = fopen(fullfile(F_SUBJECTS_PATH, 'subjectname'), 'w');
        fprintf(fileID, '%s\n', [sub, DIR_FS]);
        fclose(fileID);

        % make template
        system(sprintf('mktemplate-sess -s %s -fsd %s -funcstem f', sub, DIR_SESS));
        % run bbregister
        system(sprintf('register-sess -s %s -fsd %s -bold -per-run', sub, DIR_SESS));
    end

    %% Update registration after manual adjustments

    if do_preproc_2
        
        % during manual correction, create register.dof6_2.lta' with updated registration matrix
        cd(fullfile(F_SUBJECTS_PATH, DIR_SESS));
        if isfile(fullfile(F_SUBJECTS_PATH, DIR_SESS, 'register.dof6_2.lta'));
            % rerun registration on subject level
            system(sprintf('bbregister --mov template.nii.gz --reg register.dof6.lta --bold --init-reg %s', fullfile(F_SUBJECTS_PATH, DIR_SESS, 'register.dof6_2.lta')));
            % rerun registration for every run
            for j = 1:length(name_runs)
                run = name_runs{j};
                cd(fullfile(F_SUBJECTS_PATH, DIR_SESS, run));
                system(sprintf('bbregister --mov template.nii.gz --reg register.dof6.lta --bold --init-reg %s', fullfile(F_SUBJECTS_PATH, DIR_SESS, 'register.dof6.lta')));
            end
        end
    end

    %% After second check, run rest of FSFAST

    if do_preproc_3
        cd(PROJECT_PATH);
        % Frun rest of FSFAST pipeline
        system(sprintf('preproc-sess -s %s -fsd %s -surf self rhlh -fwhm 0 -per-run -no-reg', sub, DIR_SESS));
    end
    
    %% Sampling the signal to 9 intermediate surfaces
    % If layer analysis then generate layers
    if gen_layers
        cd(PROJECT_PATH);
        for h = 1:length(hemispheres)
            hem = hemispheres{h};
            for j = 1:length(name_runs)
                run = name_runs{j};
                for l = 1:length(layers)
                    layer = layers{l};
                    %if ~isfile(fullfile(F_SUBJECTS_PATH, sprintf('%s/%s/%s.%s.sm0.self.midgray.%s.mgz', DIR_SESS, run, hem, funcstem, layer)));
                            system(sprintf('mri_vol2surf --interp trilin --mov %s/%s/%s/%s.nii.gz --reg %s/%s/%s/register.dof6.lta --surf midgray.B1justsub.%s --hemi %s --o %s/%s/%s/%s.%s.sm0.self.midgray.%s.mgz', ...
                            sub, DIR_SESS, run, funcstem, sub, DIR_SESS, run, layer, hem, sub, DIR_SESS, run, hem, funcstem, layer));
                    %end
                end
            end
        end
    end

    %% Intracortical smoothing

    if vert_smooth
       
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

    %% Generate covariates (wm signal)

    if gen_covariates
        cd(PROJECT_PATH);
        system(sprintf('/usr/local/freesurfer/stable6_0_0/fsfast/bin/fcseed-sess -s %s -cfg %s', sub, config_file));
    end

end
