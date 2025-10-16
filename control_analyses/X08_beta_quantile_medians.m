clc;
close all;
clear;

% Description:
% This script plots the median beta for each beta quantile.
% Authors: Marianna Elisa Schmidt (marianna.schmidt@maxplanckschools.de)

%% Specifications

saveFigures = 0;

% Define the base directory where the subject folders are stored
baseDirectory = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results';

% Subjects
subNames = {'aman', 'ylri', 'auil', 'arak', 'aroo', 'atib', 'imyy', 'chss', 'evad', 'haas', 'rcgr'};
Layer = '0-10';
hemis = {'lh', 'rh'};
ROI = 'V1_Ventral';

%% Load data

% Initialize cell arrays to store distances for beta quantiles for all subjects
data_eye1 = zeros(length(subNames), length(hemis), 10);
data_eye2 = zeros(length(subNames), length(hemis), 10);

for n = 1:length(subNames)
        
    subjectDir = fullfile(baseDirectory, [ROI '_layers_' Layer '_intrahemispheric'], subNames{n});
    matFilePath = fullfile(subjectDir, 'CorrelationMtx_Selectivity_Params_subsampled.mat');
    data = load(matFilePath, 'AnalysisParam');

    for h = 1:length(hemis)

        % Store data in cell arrays
        data_eye1(n, h, :) = abs(data.AnalysisParam.betas_eye1_medians{h});
        data_eye2(n, h, :) = abs(data.AnalysisParam.betas_eye2_medians{h});

    end

end

%% Average over hemispheres

mean_data_eye1 = squeeze(mean(data_eye1, 2));
mean_data_eye2 = squeeze(mean(data_eye2, 2));

mean_data = (mean_data_eye1 + mean_data_eye2) / 2;

%% plotting

fin_data = mean(mean_data,1);

sprintf('Median beta within each beta quantile:');
mean(fin_data,1)
