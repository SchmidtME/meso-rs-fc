close all
clear all
clc

warning('off', 'all');

% Description:
% This script corresponds to Analysis B - The effect of ocular preference strength (beta),
% cortical depth, ROI (V1 subregions) (and type) on rs-FC (and selectivity).
% It corresponds to Figures 3-5 of the manuscript.
% This script runs subfunctions to calculate the rs-FC for each subject and at
% different cortical depth in V1 and within and between subregions of V1. 
% It specifies the input parameters to the subfunctions to run additional 
% preprocessing (such as detrending and high-pass filtering), and to save the 
% mean correlation for a number of distances and beta quantiles.

%% specifications

AnalysisParam.ROIpatch = "V1_patch.flat";
sub2sub = 0; 
AnalysisParam.ROI='V1' % _Dorsal'; % V1, V2, V3, V4, V1_Center, V1_Periphery
AnalysisParam.Type='intrahemispheric_vis_resp'; % intrahemispheric = within one hemisphere; interhemispheric = between hemispheres
Sbjs = {'arak', 'imyy', 'haas'}%{'aman', 'ylri', 'auil', 'aroo', 'atib', 'chss', 'evad', 'rcgr'};
layers = {'0-2', '4-6', '8-10'}; 

% data organization
AnalysisParam.StrtPnt = 10; % excludes volumes 1-9 from timeseries
AnalysisParam.TmSeries_Length = 118; % up to volume 128, so 119 in total
AnalysisParam.nQuant_dist = 10; % number of quantiles for grouping of distances from selected vertex
AnalysisParam.nQuant_beta = 10; % number of quantiles for grouping of betas from selected vertex
AnalysisParam.distThresh = 3; % minimal distance from selected vertex in mm

% additional preprocessing
AnalysisParam.detrending = 1;
AnalysisParam.hpf=1;

%% call function for processing and evtl. enable parallel processing
 
parfor layer = 1:length(layers)
    
    TrgFile = sprintf(['.fmcpr.sm0.self.midgray.00.nb1_rad','%s'], layers{layer});

    Root = sprintf('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results/Supplementary/%s_layers_%s_%s', AnalysisParam.ROI, layers{layer}, AnalysisParam.Type); 

    for i=1:length(Sbjs)

        if ~sub2sub
            
            if strcmp(AnalysisParam.Type, 'intrahemispheric_vis_resp')
                % within V1 or V1 subregions and wthin one hemisphere
                for s=1
                    B011s_Selectivity_Proc_Data_vis_resp(Sbjs{i}, Root, AnalysisParam, TrgFile, s)
                end
            elseif strcmp(AnalysisParam.Type, 'intrahemispheric_1st_run')
                % within V1 or V1 subregions and wthin one hemisphere
                for s=1
                    B012s_Selectivity_Proc_Data_1st_run(Sbjs{i}, Root, AnalysisParam, TrgFile)
                end
            elseif strcmp(AnalysisParam.Type, 'intrahemispheric_last_run')
                % within V1 or V1 subregions and wthin one hemisphere
                for s=1
                    B013s_Selectivity_Proc_Data_last_run(Sbjs{i}, Root, AnalysisParam, TrgFile)
                end
            elseif strcmp(AnalysisParam.Type, 'intrahemispheric_odc_sep_depths')
                % within V1 or V1 subregions and wthin one hemisphere
                for s=1
                    B014s_Selectivity_Proc_Data_odc_sep_depths(Sbjs{i}, Root, AnalysisParam, TrgFile, s, layers{layer})
                end
            end

        elseif sub2sub
            % between V1 subregions
            B015s_Selectivity_Proc_Data_sub2sub(Sbjs{i}, Root, AnalysisParam, TrgFile)
        end
    
    end
end
