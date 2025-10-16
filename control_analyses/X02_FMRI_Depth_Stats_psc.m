clear all
close all
clc

% Description:
% This script loads denoised fMRI data from V1 (primary visual cortex) for multiple subjects, 
% cortical layers, runs, and hemispheres. It averages the mean fMRI signal across runs and
% hemispheres and performs a repeated-measures ANOVA to assess the effect of 
% cortical depth on the fMRI signal.
% Authors: Marianna Elisa Schmidt (marianna.schmidt@maxplanckschools.de)

%% load fmcpr data

data_dir = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Data_Denoised';
label_dir = '/space/ardebil/1/users/Others/Marianna/good_subjects_anat/retinotopy';

Sbjs = {'myla','aman', 'ylri', 'auil', 'arak', 'aroo', 'atib', 'imyy', 'chss', 'evad', 'haas', 'rcgr'};
layers = {'0-2', '4-6', '8-10'}; 
hemis = {'lh', 'rh'};
ROIs = {'V1_Posterior', 'V1_Anterior', 'V1_Ventral', 'V1_Dorsal', 'V1'};

saveTables = 1;
currentDate = datestr(now, 'yyyy-mm-dd');
savePath = ['/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Figure_1/psc/' currentDate '/'];
if saveTables 
    mkdir(savePath);
end

%%
data=[];

for roi = 5%1:length(ROIs)
    ROI = ROIs{roi};

    for i = 1:length(Sbjs)

        psc = [];
    
        Sbj = Sbjs{i};
        Sbj_dir = fullfile(data_dir, Sbj, 'bold_Close_Upsampled2');
        d = dir(fullfile(Sbj_dir, '0*')); 
        runs = {d.name};
    
        labelFolder = fullfile('/space/ardebil/1/users/Others/Marianna/good_subjects_anat/retinotopy', Sbj);
        
        for l = 1:length(layers)
    
            layer = layers{l};
    
            for r = 1:length(runs)
                
                run = runs{r};
        
                for h = 1:length(hemis)
        
                    hemi = hemis{h};
                                    
                    % Define file paths
                    bold_file = fullfile(Sbj_dir, run, [hemi '.fmcpr.sm0.self.midgray.00.nb1_rad' layer '.mgz']);
                    if strcmp(ROI, 'V1')
                        roi_indices = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Upsampled_Rtopy.label']));
                    elseif strcmp(ROI, 'V1_Posterior')
                        V1_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Upsampled_Rtopy.label']));
                        V1_Posterior_label = read_ROIlabel(fullfile(labelFolder, [hemi '.Posterior_Upsampled_Rtopy.label']));
                        [~, idx] = ismember(V1_label, V1_Posterior_label);
                        roi_indices = V1_label(idx ~= 0, :);
                    elseif strcmp(ROI, 'V1_Anterior')
                        V1_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Upsampled_Rtopy.label']));
                        V1_Dorsal_label = read_ROIlabel(fullfile(labelFolder, [hemi '.Posterior_Upsampled_Rtopy.label']));
                        unique_center_labels = unique(V1_Dorsal_label);
                        rows_to_keep = ~ismember(V1_label, unique_center_labels);
                        roi_indices = V1_label(rows_to_keep, :);
                    elseif strcmp(ROI, 'V1_Dorsal')
                        V1_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Upsampled_Rtopy.label']));
                        V1_Dorsal_label = read_ROIlabel(fullfile(labelFolder, [hemi '.Dorsal_Upsampled_Rtopy.label']));
                        [~, idx] = ismember(V1_label, V1_Dorsal_label);
                        roi_indices = V1_label(idx ~= 0, :);
                    elseif strcmp(ROI, 'V1_Ventral')
                        V1_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Upsampled_Rtopy.label']));
                        V1_Dorsal_label = read_ROIlabel(fullfile(labelFolder, [hemi '.Dorsal_Upsampled_Rtopy.label']));
                        unique_center_labels = unique(V1_Dorsal_label);
                        rows_to_keep = ~ismember(V1_label, unique_center_labels);
                        roi_indices = V1_label(rows_to_keep, :);
                    end
    
                    % Load the fMRI data
                    [bold_data, mghM]= load_mgh(bold_file);

                    bold_data = squeeze(bold_data(:,:,:,10:128));
                    
                    % Keep only the bold_data corresponding to the ROI indices
                    bold_data = bold_data(roi_indices, :);

                    % Compute the mean signal per vertex (across time)
                    mean_signal = mean(bold_data, 2);
                    
                    % Replace zero means with NaN to avoid division by zero
                    mean_signal(mean_signal == 0) = NaN;
                    
                    % Expand mean_signal to match the size of bold_data for element-wise operations
                    mean_signal_expanded = repmat(mean_signal, 1, size(bold_data, 2));
                    
                    % Compute percent signal change
                    psc(h,r) = mean(mean((abs(((bold_data - mean_signal_expanded) ./ mean_signal_expanded) * 100)),2),1);
        
                end
                
            end
    
            % mean over runs
            data(i,l,roi,:) = mean(psc,2);
    
        end
    
    end
end

% average across hemispheres
data = mean(data, 4);

%% rmANOVA to assess effects of cortical depth on fMRI signal - center vs peri

STATS = [data(:,1,1) data(:,2,1) data(:,3,1)...
    data(:,1,2) data(:,2,2) data(:,3,2)];

ts = arrayfun(@(x) sprintf('t%d', x), 1:6, 'UniformOutput', false);

t =array2table(STATS,'VariableNames', ts);

within = table({'A'; 'B'; 'C'; 'A'; 'B'; 'C'}, ... % Layer
    {'AA'; 'AA'; 'AA'; 'BB'; 'BB'; 'CC'}, ... % ROI
    'VariableNames', {'Layer', 'ROI'});

rm = fitrm(t,'t1-t6~1','WithinDesign',within)
disp('Comparison of signal differences in layer between central vs. peripheral V1')
ranovatbl = ranova(rm,'WithinModel','Layer*ROI')

if saveTables
    writetable(ranovatbl, [savePath 'CentvsPeri_anova.xlsx'], 'WriteRowNames', true);
end

% Assume STATS already defined as above
nSubjects = size(STATS, 1);

% Reshape into [nSubjects x nLayers x nROIs]
STATS_reshaped = reshape(STATS, [nSubjects, 3, 2]);

% Compute mean and std over subjects
mean_vals = mean(STATS_reshaped, 1);  % [1 x 3 x 2]
std_vals  = std(STATS_reshaped, 0, 1); % [1 x 3 x 2]

% Squeeze to get [3 x 2]
mean_vals = squeeze(mean_vals);  % [3 layers x 2 ROIs]
std_vals  = squeeze(std_vals);   % [3 layers x 2 ROIs]

% Create grouped bar plot
figure(1);
hold on;
b = bar(mean_vals, 'grouped');  % Each group: a layer; each bar within: an ROI
colors = [0.2 0.6 0.8; 0.8 0.4 0.4]; % Light blue for ROI3, reddish for ROI4

for i = 1:2
    b(i).FaceColor = colors(i,:);
end

% Add error bars
ngroups = size(mean_vals, 1);  % 3 layers
nbars = size(mean_vals, 2);    % 2 ROIs
groupwidth = min(0.8, nbars/(nbars + 1.5));

for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, mean_vals(:,i), std_vals(:,i), 'k', 'linestyle', 'none', 'linewidth', 1);
end

% Aesthetics
xticks(1:3);
xticklabels({'Deep', 'Middle', 'Superficial'});
xlabel('Cortical Layer');
ylabel('Metric');
%ylim([800 1600]);
legend({'Center', 'Periphery'}, 'Location', 'southeast');
title('Mean ± SD across Subjects per Layer and ROI');
set(gca, 'FontSize', 14);
box on;

%% rmANOVA to assess effects of cortical depth on fMRI signal - dorsal vs. ventral

STATS = [data(:,1,3) data(:,2,3) data(:,3,3)...
    data(:,1,4) data(:,2,4) data(:,3,4)];

ts = arrayfun(@(x) sprintf('t%d', x), 1:6, 'UniformOutput', false);

t =array2table(STATS,'VariableNames', ts);

within = table({'A'; 'B'; 'C'; 'A'; 'B'; 'C'}, ... % Layer
    {'AA'; 'AA'; 'AA'; 'BB'; 'BB'; 'CC'}, ... % ROI
    'VariableNames', {'Layer', 'ROI'});

rm = fitrm(t,'t1-t6~1','WithinDesign',within)
disp('Comparison of signal differences in layer between dorsal and ventral V1')
ranovatbl = ranova(rm,'WithinModel','Layer*ROI')

if saveTables
    writetable(ranovatbl, [savePath 'DorsVsVent_anova.xlsx'], 'WriteRowNames', true);
end

% Figure

% Assume STATS already defined as above
nSubjects = size(STATS, 1);

% Reshape into [nSubjects x nLayers x nROIs]
STATS_reshaped = reshape(STATS, [nSubjects, 3, 2]);

% Compute mean and std over subjects
mean_vals = mean(STATS_reshaped, 1);  % [1 x 3 x 2]
std_vals  = std(STATS_reshaped, 0, 1); % [1 x 3 x 2]

% Squeeze to get [3 x 2]
mean_vals = squeeze(mean_vals);  % [3 layers x 2 ROIs]
std_vals  = squeeze(std_vals);   % [3 layers x 2 ROIs]

% Create grouped bar plot
figure(2);
hold on;
b = bar(mean_vals, 'grouped');  % Each group: a layer; each bar within: an ROI
colors = [0.2 0.6 0.8; 0.8 0.4 0.4]; % Light blue for ROI3, reddish for ROI4

for i = 1:2
    b(i).FaceColor = colors(i,:);
end

% Add error bars
ngroups = size(mean_vals, 1);  % 3 layers
nbars = size(mean_vals, 2);    % 2 ROIs
groupwidth = min(0.8, nbars/(nbars + 1.5));

for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, mean_vals(:,i), std_vals(:,i), 'k', 'linestyle', 'none', 'linewidth', 1);
end

% Aesthetics
xticks(1:3);
xticklabels({'Deep', 'Middle', 'Superficial'});
xlabel('Cortical Layer');
ylabel('Metric');
%ylim([800 1600]);
legend({'Ventral', 'Dorsal'}, 'Location', 'southeast');
title('Mean ± SD across Subjects per Layer and ROI');
set(gca, 'FontSize', 14);
box on;

%% rmANOVA to assess effects of cortical depth on fMRI signal - center

STATS = [data(:,1,1) data(:,2,1) data(:,3,1)];

ts = arrayfun(@(x) sprintf('t%d', x), 1:3, 'UniformOutput', false);

t =array2table(STATS,'VariableNames', ts);

within = table({'A'; 'B'; 'C'}, ... 
    'VariableNames', {'Layer'});

rm = fitrm(t,'t1-t3~1','WithinDesign',within)
disp('Comparison of signal differences in layer for central V1')
ranovatbl = ranova(rm,'WithinModel','Layer')

if saveTables
    writetable(ranovatbl, [savePath 'Central_V1_anova.xlsx'], 'WriteRowNames', true);
end

%% rmANOVA to assess effects of cortical depth on fMRI signal - peri

STATS = [data(:,1,2) data(:,2,2) data(:,3,2)];

ts = arrayfun(@(x) sprintf('t%d', x), 1:3, 'UniformOutput', false);

t =array2table(STATS,'VariableNames', ts);

within = table({'A'; 'B'; 'C'}, ... 
    'VariableNames', {'Layer'});

rm = fitrm(t,'t1-t3~1','WithinDesign',within)
disp('Comparison of signal differences in layer for peripheral V1')
ranovatbl = ranova(rm,'WithinModel','Layer')

if saveTables
    writetable(ranovatbl, [savePath 'Peripheral_V1_anova.xlsx'], 'WriteRowNames', true);
end

%% rmANOVA to assess effects of cortical depth on fMRI signal - ventral

STATS = [data(:,1,3) data(:,2,3) data(:,3,3)];

ts = arrayfun(@(x) sprintf('t%d', x), 1:3, 'UniformOutput', false);

t =array2table(STATS,'VariableNames', ts);

within = table({'A'; 'B'; 'C'}, ... 
    'VariableNames', {'Layer'});

rm = fitrm(t,'t1-t3~1','WithinDesign',within)
disp('Comparison of signal differences in layer for ventral V1')
ranovatbl = ranova(rm,'WithinModel','Layer')

if saveTables
    writetable(ranovatbl, [savePath 'Ventral_V1_anova.xlsx'], 'WriteRowNames', true);
end

%% rmANOVA to assess effects of cortical depth on fMRI signal - dorsal

STATS = [data(:,1,4) data(:,2,4) data(:,3,4)];

ts = arrayfun(@(x) sprintf('t%d', x), 1:3, 'UniformOutput', false);

t =array2table(STATS,'VariableNames', ts);

within = table({'A'; 'B'; 'C'}, ... 
    'VariableNames', {'Layer'});

rm = fitrm(t,'t1-t3~1','WithinDesign',within)
disp('Comparison of signal differences in layer for dorsal V1')
ranovatbl = ranova(rm,'WithinModel','Layer')

if saveTables
    writetable(ranovatbl, [savePath 'Dorsal_V1_anova.xlsx'], 'WriteRowNames', true);
end

%% rmANOVA to assess effects of cortical depth on fMRI signal - V1

STATS = [data(:,1,5) data(:,2,5) data(:,3,5)];

ts = arrayfun(@(x) sprintf('t%d', x), 1:3, 'UniformOutput', false);

t =array2table(STATS,'VariableNames', ts);

within = table({'A'; 'B'; 'C'}, ... 
    'VariableNames', {'Layer'});

rm = fitrm(t,'t1-t3~1','WithinDesign',within)
disp('Comparison of signal differences in layer for V1')
ranovatbl = ranova(rm,'WithinModel','Layer')

if saveTables
    writetable(ranovatbl, [savePath 'V1_anova.xlsx'], 'WriteRowNames', true);
end