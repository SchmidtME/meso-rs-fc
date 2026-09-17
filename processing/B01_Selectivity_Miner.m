close all
clear all
clc

warning('off', 'all');

% Description:
% Driver script for Analysis A & B - The effect of distance, 
% ocular preference strength (beta), cortical depth, ROI (V1 subregions)
% on rs-FC (and selectivity).
% It corresponds to Figures 2-7 of the manuscript.
% It defines the analysis parameters and loops over layers and subjects, routing
% to either the intra-hemispheric or inter-hemispheric per-subject function, and
% runs subfunctions to calculate the rs-FC for each subject, at
% different cortical depth in V1 and within and between subregions of V1 in
% the left and right hemisphere or between the hemispheres. 
% It specifies the input parameters to the subfunctions to run additional 
% preprocessing (such as detrending and high-pass filtering), and to save the 
% mean correlation for a number of distances and beta quantiles as per-subject
% .mat result files.
% Authors: Marianna E. Schmidt (marianna.schmidt@maxplanckschools.de), Iman Aganj, Shahin Nasr

%% specifications
% --- Define the ROI, connectivity type, subjects and layers to process ---

AnalysisParam.ROIpatch = "V1_patch.flat";
sub2sub = 0; 
AnalysisParam.ROI='V1_Posterior' % V1, V1_Posterior (Center), V1_Anterior (Periphery), V1_Dorsal, V1_Ventral
AnalysisParam.Type='interhemispheric'; % intrahemispheric = within one hemisphere; interhemispheric = between hemispheres
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

%% call function for processing and evtl. enable parallel processing
% --- Outer loop: for each layer run the per-subject processing function ---
 
for layer = 1:length(layers)
    
    % Define the target file for the current layer (surface overlay name with layer tag)
    TrgFile = sprintf(['.fmcpr.sm0.self.midgray.00.nb1_rad','%s'], layers{layer});

    % Set the root directory for results (one per ROI / layer / connectivity type)
    Root = sprintf('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results/%s_layers_%s_%s', AnalysisParam.ROI, layers{layer}, AnalysisParam.Type); 

    for i=1:length(Sbjs)
            
	    if strcmp(AnalysisParam.Type, 'intrahemispheric')
		    % within V1 or V1 subregions and wthin one hemisphere
		for s=1 % subsampling iteration
		    B001a_Selectivity_Proc_Data(Sbjs{i}, Root, AnalysisParam, TrgFile,s)
		end
	    
	    elseif strcmp(AnalysisParam.Type, 'interhemispheric')
		    % within V1 or V1 subregions and between hemispheres
		    B001b_Selectivity_Proc_Data_Interhemi(Sbjs{i}, Root, AnalysisParam, TrgFile)
        end
    
    end
end
