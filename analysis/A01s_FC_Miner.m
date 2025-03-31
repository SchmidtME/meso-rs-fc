close all
clear all
clc

% Description:
% This script corresponds to supplementary analysisof Analysis A -
% The effect of distance (& type) on rs-FC and selectivity..
% This script runs subfunctions to calculate the rs-FC for each subject and at
% different cortical depth in V1 and subregions of V1. 
% It specifies the input parameters to the subfunctions to run additional 
% preprocessing (such as detrending and high-pass filtering), and to save the 
% mean correlation for a number of distances.

%% Specifications

% Define the region of interest (ROI) patch and specific area of interest (e.g., V1, V2, etc.)
AnalysisParam.ROIpatch = "V1_patch.flat"; % Patch file for the ROI
AnalysisParam.ROI = 'V1'; % Choose the ROI: V1, V2, V3, V4, V1_Center, V1_Periphery
AnalysisParam.Type='intrahemispheric_vis_resp'; % Define type of analysis: intrahemispheric (within one hemisphere) or interhemispheric (between hemispheres)

% List of subjects to process
Sbjs = {'arak', 'imyy', 'haas'}%'rcgr'};
% 'ylri', 'aman', 'auil', 'aroo', 'atib', 'chss', 'evad',
% 'arak', 'imyy', 'haas',
% Layers of interest for analysis
layers = {'0-10'};

%% Data Organization and Preprocessing Parameters

% Exclude the first 9 volumes from the time series for each subject
AnalysisParam.StrtPnt = 10; % Start from volume 10
AnalysisParam.TmSeries_Length = 118; % Total time series length (119 volumes)
AnalysisParam.nQuant_dist = 10; % Number of quantiles for grouping distances from the selected vertex
AnalysisParam.nQuant_beta = 10; % Number of quantiles for grouping beta values from the selected vertex
AnalysisParam.distThresh = 3; % Minimum distance from the selected vertex in mm

% Preprocessing steps (detrending and high-pass filtering)
AnalysisParam.detrending = 1; % Enable detrending of the time series
AnalysisParam.hpf = 1; % Enable high-pass filtering

%% Process Data for Each Layer and Subject

% Loop through each layer
for layer = 1:length(layers)
    
    % Define the target file for the current layer
    TrgFile = sprintf(['.fmcpr.sm0.self.midgray.00.nb1_rad','%s'], layers{layer});

    % Set the root directory for results
    Root = sprintf('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results/Supplementary/%s_layers_%s_%s', AnalysisParam.ROI, layers{layer}, AnalysisParam.Type); 

    % Parallel processing for each subject
    for i = 1:length(Sbjs)

        % supplementary analysis A001s: check for differences in rs-FC for vertices with a string vs. weak visual response
        if strcmp(AnalysisParam.Type, 'intrahemispheric_vis_resp')
            for s = 1
                A001s_FC_Proc_Data_subsample_vis_resp(Sbjs{i}, Root, AnalysisParam, TrgFile, s)
            end

        % supplementary analysis A002s: check for differences in rs-FC between 1st and last run
        elseif strcmp(AnalysisParam.Type, 'intrahemispheric_1st_run')
            for s=1
                A002as_FC_Proc_Data_subsample_1st_run(Sbjs{i}, Root, AnalysisParam, TrgFile, s)
            end
        elseif strcmp(AnalysisParam.Type, 'intrahemispheric_last_run')
            for s=1
                A002bs_FC_Proc_Data_subsample_last_run(Sbjs{i}, Root, AnalysisParam, TrgFile, s)
            end

        % supplementary analysis A003s: check for differences in rs-FC if odc map is sampled from the corresponding cortical depth
        elseif strcmp(AnalysisParam.Type, 'intrahemispheric_odc_sep_depths')
            for s=1
                A003s_FC_Proc_Data_subsample_odc_sep_depths(Sbjs{i}, Root, AnalysisParam, TrgFile, s, layers{layer})
            end

        % supplementary analysis A004s:interhemispheric rs-FC instead of vertex2vertex, avg V1 signal 2 avg V1 signal
        elseif strcmp(AnalysisParam.Type, 'interhemispheric_avg')
            A004s_FC_Proc_Data_interhemi_avg(Sbjs{i}, Root, AnalysisParam, TrgFile)
        end
    end
end
