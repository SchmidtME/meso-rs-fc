clear all
close all
clc

% Description:
% This script loads denoised fMRI data from V1 (primary visual cortex) for multiple subjects, 
% cortical layers, runs, and hemispheres. It averages the mean fMRI signal across runs and
% hemispheres and performs a repeated-measures ANOVA to assess the effect of 
% cortical depth on the fMRI signal.

%% load fmcpr data

data_dir = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Data_Denoised';
label_dir = '/space/ardebil/1/users/Others/Marianna/good_subjects_anat/retinotopy';

Sbjs = {'aman', 'ylri', 'auil', 'arak', 'aroo', 'atib', 'imyy', 'chss', 'evad', 'haas', 'rcgr'};
layers = {'0-2', '4-6', '8-10'}; 
hemis = {'lh', 'rh'};

for i = 1:length(Sbjs)

    Sbj = Sbjs{i};
    Sbj_dir = fullfile(data_dir, Sbj, 'bold_Close_Upsampled2');
    d = dir(fullfile(Sbj_dir, '0*')); 
    runs = {d.name};
    
    for l = 1:length(layers)

        layer = layers{l};

        for r = 1:length(runs)
            
            run = runs{r};
    
            for h = 1:length(hemis)
    
                hemi = hemis{h};
                                
                % Define file paths
                bold_file = fullfile(Sbj_dir, run, [hemi '.fmcpr.sm0.self.midgray.00.nb1_rad' layer '.mgz']);
                label_file = fullfile(label_dir, Sbj, [hemi '.V1_Upsampled_Rtopy.label']);
                
                % Load the fMRI data
                [bold_data, mghM]= load_mgh(bold_file);
                
                % Load the label file (V1 mask)
                roi_indices = read_ROIlabel(label_file);
                
                % Keep only the bold_data corresponding to the ROI indices
                bold_data_filtered(h,r) = mean(mean(bold_data(roi_indices, :)));
    
            end
            
        end

        data(i,l) = mean(mean(bold_data_filtered));

    end

end

%% rmANOVA to assess effects of cortical depth on fMRI signal

STATS = data;

ts = arrayfun(@(x) sprintf('t%d', x), 1:3, 'UniformOutput', false);

t =array2table(STATS,'VariableNames', ts);

within = table({'A'; 'B'; 'C'}, ...
    'VariableNames', {'Layer'});

rm = fitrm(t,'t1-t3~1','WithinDesign',within)
disp('Comparison between layer in V1 - averaged over hemispheres, runs')
ranovatbl = ranova(rm,'WithinModel','Layer')
