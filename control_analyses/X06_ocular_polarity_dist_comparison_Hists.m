clc;
close all;
clear;

%% Specifications

saveFigures = 0;

% Define the base directory where the subject folders are stored
baseDirectory = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results';

% Subjects
subNames = {'aman', 'ylri', 'auil', 'arak', 'aroo', 'atib', 'imyy', 'chss', 'evad', 'haas', 'rcgr'};
layers = {'0-2'};
hemis = {'lh', 'rh'};

%% Load data

% Initialize cell arrays to store distances for beta quantiles for all subjects
data_eye1eye2 = cell(length(subNames), length(layers), length(hemis));
data_same_eye = cell(length(subNames), length(layers), length(hemis));

for n = 1:length(subNames)
    for l = 1:length(layers)
        
        subjectDir = fullfile(baseDirectory, ['V1_layers_' layers{l} '_intrahemispheric'], subNames{n});
        matFilePath = fullfile(subjectDir, 'CorrelationMtx_Selectivity_Params_subsampled.mat');
        data = load(matFilePath, 'AnalysisParam');

        for h = 1:length(hemis)

            % Store data in cell arrays
            data_same_eye{n, l, h} = squeeze(mean(mean(data.AnalysisParam.same_eye_median_quants_dist{h},1),2));
            data_eye1eye2{n, l, h} = squeeze(mean(mean(data.AnalysisParam.eye1eye2_median_quants_dist{h},1),2));

        end
    end
end

%% Average over hemispheres

data_eye1eye2_hemis = cell(length(subNames),1);
data_same_eye_hemis = cell(length(subNames),1);

for n = 1:length(subNames)

    % unalike
    max_rows = max(size(data_eye1eye2{n, 1, 1},1), size(data_eye1eye2{n, 1, 2},1));
    max_cols = max(size(data_eye1eye2{n, 1, 1},2), size(data_eye1eye2{n, 1, 2},2));
    padded_hemi_1 = NaN(max_rows(1), max_cols(1));
    padded_hemi_2 = NaN(max_rows(1), max_cols(1));
    padded_hemi_1(1:size(data_eye1eye2{n, 1, 1},1), 1:size(data_eye1eye2{n, 1, 1},2)) = data_eye1eye2{n, 1, 1};
    padded_hemi_2(1:size(data_eye1eye2{n, 1, 2},1), 1:size(data_eye1eye2{n, 1, 2},2)) = data_eye1eye2{n, 1, 2};
    data_eye1eye2_hemis_mean(n) = mean(vertcat(padded_hemi_1, padded_hemi_2), 'all', 'omitnan');
    data_eye1eye2_hemis{n} = vertcat(padded_hemi_1, padded_hemi_2);
    data_eye1eye2_combined(n,:) = (padded_hemi_1 + padded_hemi_2)/2;
    data_eye1eye2_hemis{n}(data_eye1eye2_hemis{n} < 3) = NaN;

    % alike
    max_rows = max(size(data_same_eye{n, 1, 1},1), size(data_same_eye{n, 1, 2},1));
    max_cols = max(size(data_same_eye{n, 1, 1},2), size(data_same_eye{n, 1, 2},2));
    padded_hemi_1 = NaN(max_rows(1), max_cols(1));
    padded_hemi_2 = NaN(max_rows(1), max_cols(1));
    padded_hemi_1(1:size(data_same_eye{n, 1, 1},1), 1:size(data_same_eye{n, 1, 1},2)) = data_same_eye{n, 1, 1};
    padded_hemi_2(1:size(data_same_eye{n, 1, 2},1), 1:size(data_same_eye{n, 1, 2},2)) = data_same_eye{n, 1, 2};
    data_same_eye_hemis_mean(n) = mean(vertcat(padded_hemi_1, padded_hemi_2), 'all', 'omitnan');
    data_same_eye_hemis{n} = vertcat(padded_hemi_1, padded_hemi_2);
    data_same_eye_combined(n,:) = (padded_hemi_1 + padded_hemi_2)/2;
    data_same_eye_hemis{n}(data_same_eye_hemis{n} < 3) = NaN;
    
    max_data(n,1) = max(data_same_eye_hemis{n}, [], 'all');
    max_data(n,2) = max(data_eye1eye2_hemis{n}, [], 'all');
end

%% ttest

[h,p] = ttest(mean(data_same_eye_combined,2), mean(data_eye1eye2_combined,2))

%% plotting

fin_data = (data_same_eye_combined + data_eye1eye2_combined)/2;

sprintf('Median distance within each distance quantile:');
mean(fin_data,1)
