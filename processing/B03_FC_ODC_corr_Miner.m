close all
clear all
clc

% Description:
% This script corresponds to Analysis Z - the correlation of rs-FC with the
% ODC differential map to Figure 6 of the manuscript.
% This script calls subfunctions that calculate the correlation of rs-FC 
% with either a spatially shifted or rotated version of the differential ODC map
% and saves the resulting correlation coefficients for each subject.

%% specifications

AnalysisParam.ROIpatch = "V1_patch.flat";
rotated_H0 = 1; % if 1: 180 deg rotated ODC as H0, else spatially shifted ODC as H0
AnalysisParam.ROI='V1' % _Dorsal'; % V1, V2, V3, V4, V1_Center, V1_Periphery
AnalysisParam.Type='intrahemispheric'; % intrahemispheric = within one hemisphere; interhemispheric = between hemispheres
Sbjs = {'aman', 'ylri', 'auil', 'arak', 'aroo', 'atib', 'imyy', 'chss', 'evad', 'haas', 'rcgr'};
layers = {'0-2', '4-6', '8-10', '0-10'}; 

% data organization
AnalysisParam.StrtPnt = 10; % excludes volumes 1-9 from timeseries
AnalysisParam.TmSeries_Length = 118; % up to volume 128, so 119 in total
AnalysisParam.nQuant_dist = 10; % number of quantiles for grouping of distances from selected vertex
AnalysisParam.nQuant_beta = 10; % number of quantiles for grouping of betas from selected vertex
AnalysisParam.distThresh = 3; % minimal distance from selected vertex in mm

% additional preprocessing
AnalysisParam.detrending = 1;
AnalysisParam.hpf=1;

%%  Call the subfunction to process data
 
for layer = 1:length(layers)
    
    TrgFile = sprintf(['.fmcpr.sm0.self.midgray.00.nb1_rad', '%s'], layers{layer});

    Root = sprintf('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results/%s_layers_%s_%s', AnalysisParam.ROI, layers{layer}, AnalysisParam.Type); 

    parfor i=1:length(Sbjs)

        if ~rotated_H0
            Z001a_FC_ODC_corr_Proc_Data(Sbjs{i}, Root, AnalysisParam, TrgFile)
        else
            Z001b_FC_ODC_corr_Proc_Data_rotated(Sbjs{i}, Root, AnalysisParam, TrgFile)
        end

    end

end

% Shutdown parallel pool
delete(gcp('nocreate'));