
close all
clear all
clc

% Description:
% Driver script for Analysis Z - the correlation of rs-FC with the
% ODC differential map and to Figure 8 of the manuscript.
% It defines the analysis parameters, loops over layers and subjects and calls
% subfunctions that calculate the correlation of rs-FC 
% with either a spatially shifted or rotated version of the differential ODC map,
% saving the per-subject meanAbsR .mat results.
% Authors: Marianna E. Schmidt (marianna.schmidt@maxplanckschools.de), Iman Aganj, Shahin Nasr

%% specifications
% --- Define the ROI, the null hypothesis type, and subjects/layers to process ---

AnalysisParam.ROIpatch = "V1_patch.flat";
rotated_H0 = 1; % if 1: 180 deg rotated ODC as H0, else spatially shifted ODC as H0
AnalysisParam.ROI='V1' % _Dorsal'; % V1, V2, V3, V4, V1_Center, V1_Periphery
AnalysisParam.Type='intrahemispheric'; % intrahemispheric = within one hemisphere; interhemispheric = between hemispheres
Sbjs = {'myla'}%'aman', 'ylri', 'auil', 'arak', 'aroo', 'atib', 'imyy', 'chss', 'evad', 'haas', 'rcgr'};
layers = {'0-2', '4-6', '8-10', '0-10'}; 

% data organization
% --- Set time-series window, quantile counts, distance threshold and preprocessing flags ---
AnalysisParam.StrtPnt = 10; % excludes volumes 1-9 from timeseries
AnalysisParam.TmSeries_Length = 118; % up to volume 128, so 119 in total
AnalysisParam.nQuant_dist = 10; % number of quantiles for grouping of distances from selected vertex
AnalysisParam.nQuant_beta = 10; % number of quantiles for grouping of betas from selected vertex
AnalysisParam.distThresh = 3; % minimal distance from selected vertex in mm

% additional preprocessing
AnalysisParam.detrending = 1;
AnalysisParam.hpf=1;

%%  Call the subfunction to process data
% --- Outer loop: for each layer, run the per-subject processing function ---
 
for layer = 1:length(layers)
    
    % Define the target file for the current layer (surface overlay name with layer tag)
    TrgFile = sprintf(['.fmcpr.sm0.self.midgray.00.nb1_rad', '%s'], layers{layer});

    % Set the root directory for results (one per ROI / layer / connectivity type)
    Root = sprintf('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results/%s_layers_%s_%s', AnalysisParam.ROI, layers{layer}, AnalysisParam.Type); 

    for i=1:length(Sbjs)

        % Route to the per-subject function based on the chosen null hypothesis
        if ~rotated_H0
            Z001a_FC_ODC_corr_Proc_Data(Sbjs{i}, Root, AnalysisParam, TrgFile) % shifted-ODC null
        else
            Z001b_FC_ODC_corr_Proc_Data_rotated(Sbjs{i}, Root, AnalysisParam, TrgFile) % rotated-ODC null
        end

    end

end

% Shutdown parallel pool
delete(gcp('nocreate'));