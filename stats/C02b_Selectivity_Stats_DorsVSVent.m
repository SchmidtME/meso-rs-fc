close all
clear all
clc

% Description:
% This script corresponds to Analysis B - The effect of beta, type, and layer on rs-FC and selectivity.
% This is the analysis comparing V1 subregions. It corresponds to Figure 4 of the manuscript. 
% The data is averaged over hemispheres and distances. The first part assesses differences between 
% V1_Dorsal and V1_Ventral. The second part assesses differences 
% between V1_Ventral and V1_Dorsal.
% First, an rm ANOVA is computed to assess the effects of beta, layer, ROI and type on rs-FC.
% Then, an rm ANOVA is computed to assess the effects of beta, layer, and ROI on selectivity.

%% Load the data for V1 Dorsal and Ventral

Root = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results';

saveTables = 1;
currentDate = datestr(now, 'yyyy-mm-dd');
savePath = ['/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Figure_4/' currentDate '/'];
if saveTables
    mkdir(savePath);
end

Labels = {'V1_Dorsal', 'V1_Ventral'}; 

Layers = {'0-2', '4-6', '8-10'};
Sbjs = {'aman', 'arak', 'aroo', 'atib', 'auil', 'chss', 'evad', 'haas', 'rcgr', 'ylri', 'imyy'};
Dist=1:10;

for Label = 1:length(Labels)

    for Layer = 1:length(Layers)

        String = [Labels{Label}, '_layers_', Layers{Layer}, '_intrahemispheric/'];

        for Sbj = 1:length(Sbjs)

            load([Root,  '/', String, Sbjs{Sbj},'/CorrelationMtx_Selectivity_subsampled.mat'])
                
            % 10     4     2    10    10     2
            % Distance, Alike_both/alike_Eye1/alike_Eye2/Unalile (1), Hemi, Beta1, Beta2, z(2)

            DATA_ALIKE(Sbj, Dist, Label, Layer, :, :, :) = Data_Combined(Dist, 1, :, :, :, 2);

            DATA_UNALIKE(Sbj, Dist, Label, Layer, :, :, :) = Data_Combined(Dist, 4, :, :, :, 2);

            DATA_DIFF(Sbj, Dist, Label, Layer, :, :, :) = DATA_ALIKE(Sbj, Dist, Label, Layer, :, :, :) - DATA_UNALIKE(Sbj, Dist, Label, Layer, :, :, :);

            clear Data_Combined

        end

    end
end

% size: 11        10         2        3         2              10    10
%       subjecs   distances  labels   layers    hemispheres    betas

% mean over distances, hemispheres for V1
DATA_ALIKE_mean = squeeze(mean(mean(DATA_ALIKE(:,:,:,:,:,:,:), 5), 2));
DATA_UNALIKE_mean = squeeze(mean(mean(DATA_UNALIKE(:,:,:,:,:,:,:), 5), 2));

DATA(:,:,:,:,:,1) = DATA_ALIKE_mean;
DATA(:,:,:,:,:,2) = DATA_UNALIKE_mean;

%% rmANOVA cumulative

STATS = [
    mean(mean(DATA_ALIKE_mean(:,1,1,1:2,1:2),4),5) mean(mean(DATA_ALIKE_mean(:,1,2,1:2,1:2),4),5) mean(mean(DATA_ALIKE_mean(:,1,3,1:2,1:2),4),5)...
    mean(mean(DATA_ALIKE_mean(:,2,1,1:2,1:2),4),5) mean(mean(DATA_ALIKE_mean(:,2,2,1:2,1:2),4),5) mean(mean(DATA_ALIKE_mean(:,2,3,1:2,1:2),4),5)...
    mean(mean(DATA_ALIKE_mean(:,1,1,1:4,1:4),4),5) mean(mean(DATA_ALIKE_mean(:,1,2,1:4,1:4),4),5) mean(mean(DATA_ALIKE_mean(:,1,3,1:4,1:4),4),5)...
    mean(mean(DATA_ALIKE_mean(:,2,1,1:4,1:4),4),5) mean(mean(DATA_ALIKE_mean(:,2,2,1:4,1:4),4),5) mean(mean(DATA_ALIKE_mean(:,2,3,1:4,1:4),4),5)...
    mean(mean(DATA_ALIKE_mean(:,1,1,1:6,1:6),4),5) mean(mean(DATA_ALIKE_mean(:,1,2,1:6,1:6),4),5) mean(mean(DATA_ALIKE_mean(:,1,3,1:6,1:6),4),5)...
    mean(mean(DATA_ALIKE_mean(:,2,1,1:6,1:6),4),5) mean(mean(DATA_ALIKE_mean(:,2,2,1:6,1:6),4),5) mean(mean(DATA_ALIKE_mean(:,2,3,1:6,1:6),4),5)...
    mean(mean(DATA_ALIKE_mean(:,1,1,1:8,1:8),4),5) mean(mean(DATA_ALIKE_mean(:,1,2,1:8,1:8),4),5) mean(mean(DATA_ALIKE_mean(:,1,3,1:8,1:8),4),5)...
    mean(mean(DATA_ALIKE_mean(:,2,1,1:8,1:8),4),5) mean(mean(DATA_ALIKE_mean(:,2,2,1:8,1:8),4),5) mean(mean(DATA_ALIKE_mean(:,2,3,1:8,1:8),4),5)...
    mean(mean(DATA_ALIKE_mean(:,1,1,1:10,1:10),4),5) mean(mean(DATA_ALIKE_mean(:,1,2,1:10,1:10),4),5) mean(mean(DATA_ALIKE_mean(:,1,3,1:10,1:10),4),5)...
    mean(mean(DATA_ALIKE_mean(:,2,1,1:10,1:10),4),5) mean(mean(DATA_ALIKE_mean(:,2,2,1:10,1:10),4),5) mean(mean(DATA_ALIKE_mean(:,2,3,1:10,1:10),4),5)...
    mean(mean(DATA_UNALIKE_mean(:,1,1,1:2,1:2),4),5) mean(mean(DATA_UNALIKE_mean(:,1,2,1:2,1:2),4),5) mean(mean(DATA_UNALIKE_mean(:,1,3,1:2,1:2),4),5)...
    mean(mean(DATA_UNALIKE_mean(:,2,1,1:2,1:2),4),5) mean(mean(DATA_UNALIKE_mean(:,2,2,1:2,1:2),4),5) mean(mean(DATA_UNALIKE_mean(:,2,3,1:2,1:2),4),5)...
    mean(mean(DATA_UNALIKE_mean(:,1,1,1:4,1:4),4),5) mean(mean(DATA_UNALIKE_mean(:,1,2,1:4,1:4),4),5) mean(mean(DATA_UNALIKE_mean(:,1,3,1:4,1:4),4),5)...
    mean(mean(DATA_UNALIKE_mean(:,2,1,1:4,1:4),4),5) mean(mean(DATA_UNALIKE_mean(:,2,2,1:4,1:4),4),5) mean(mean(DATA_UNALIKE_mean(:,2,3,1:4,1:4),4),5)...
    mean(mean(DATA_UNALIKE_mean(:,1,1,1:6,1:6),4),5) mean(mean(DATA_UNALIKE_mean(:,1,2,1:6,1:6),4),5) mean(mean(DATA_UNALIKE_mean(:,1,3,1:6,1:6),4),5)...
    mean(mean(DATA_UNALIKE_mean(:,2,1,1:6,1:6),4),5) mean(mean(DATA_UNALIKE_mean(:,2,2,1:6,1:6),4),5) mean(mean(DATA_UNALIKE_mean(:,2,3,1:6,1:6),4),5)...
    mean(mean(DATA_UNALIKE_mean(:,1,1,1:8,1:8),4),5) mean(mean(DATA_UNALIKE_mean(:,1,2,1:8,1:8),4),5) mean(mean(DATA_UNALIKE_mean(:,1,3,1:8,1:8),4),5)...
    mean(mean(DATA_UNALIKE_mean(:,2,1,1:8,1:8),4),5) mean(mean(DATA_UNALIKE_mean(:,2,2,1:8,1:8),4),5) mean(mean(DATA_UNALIKE_mean(:,2,3,1:8,1:8),4),5)...
    mean(mean(DATA_UNALIKE_mean(:,1,1,1:10,1:10),4),5) mean(mean(DATA_UNALIKE_mean(:,1,2,1:10,1:10),4),5) mean(mean(DATA_UNALIKE_mean(:,1,3,1:10,1:10),4),5)...
    mean(mean(DATA_UNALIKE_mean(:,2,1,1:10,1:10),4),5) mean(mean(DATA_UNALIKE_mean(:,2,2,1:10,1:10),4),5) mean(mean(DATA_UNALIKE_mean(:,2,3,1:10,1:10),4),5)...
];

ts = arrayfun(@(x) sprintf('t%d', x), 1:60, 'UniformOutput', false);

t =array2table(STATS,'VariableNames', ts);

within = table(...
    repmat({'A'; 'B'; 'C'}, 20, 1), ... % layers
    repmat({'AA'; 'AA'; 'AA'; 'BB'; 'BB'; 'BB'}, 10, 1), ... % labels
    [repmat({'AAA'}, 6, 1); repmat({'BBB'}, 6, 1); repmat({'CCC'}, 6, 1); repmat({'DDD'}, 6, 1); repmat({'EEE'}, 6, 1); repmat({'AAA'}, 6, 1); repmat({'BBB'}, 6, 1); repmat({'CCC'}, 6, 1); repmat({'DDD'}, 6, 1); repmat({'EEE'}, 6, 1)], ... % beta
    [repmat({'AAAA'}, 30, 1); repmat({'BBBB'}, 30, 1)], ... % type
    'VariableNames', {'Layer', 'ROI', 'Beta', 'Type'});

rm = fitrm(t,'t1-t60~1','WithinDesign',within)

disp('V1 Dorsal vs Ventral - Effect of type, beta, layer and ROI on rs-FC - averaged over hemispheres, distances - accumulative beta')
ranovatbl = ranova(rm,'WithinModel','Layer+ROI+Beta+Type+Layer*ROI+Layer*Beta+Layer*Type+ROI*Beta+ROI*Type+Beta*Type+Layer*Beta*ROI+Layer*Beta*Type+Beta*ROI*Type+Layer*Beta*ROI*Type')

if saveTables
    writetable(ranovatbl, [savePath 'DorsVSVent_anova_cum.xlsx'], 'WriteRowNames', true);
end

%% LME - Ventral vs. Dorsal, mean over distances

clear DATA
% mean over distances
DATA_ALIKE_mean = squeeze(mean(DATA_ALIKE, 2));
DATA_UNALIKE_mean = squeeze(mean(DATA_UNALIKE, 2));

DATA(:,:,:,:,:,:,1) = DATA_UNALIKE_mean;
DATA(:,:,:,:,:,:,2) = DATA_ALIKE_mean;

sz = size(DATA);
Mask = true(sz) & permute(triu(true(sz(5:6))), [3 4 5 6 1 2 7]);

sizeInd = arrayfun(@(s) 1:s, size(DATA), 'UniformOutput', false);
[Subject, ROI, Layer, Hemi, BQ1, BQ2, TypeODC] = ndgrid(sizeInd{:});
[Subject, ROI, Layer, Hemi, TypeODC] = deal(categorical(Subject), categorical(ROI), categorical(Layer), categorical(Hemi), categorical(TypeODC));

prodBQ = BQ1 .* BQ2;

T = table(DATA(Mask), Subject(Mask), ROI(Mask), Layer(Mask), Hemi(Mask), prodBQ(Mask), TypeODC(Mask), 'VariableNames', {'rsFC', 'Subject', 'ROI', 'Layer', 'Hemi', 'prodBQ', 'TypeODC'});

lme = fitlme(T, 'rsFC ~ prodBQ*Layer*TypeODC*ROI + (1|Subject)');

if saveTables
    [xx, xxx, Coefficients] = fixedEffects(lme, 'Alpha', 0.05);
    resultsTable = table(lme.CoefficientNames', Coefficients.Estimate, Coefficients.SE, Coefficients.tStat, Coefficients.DF, Coefficients.pValue, Coefficients.Upper, Coefficients.Lower, ...
                     'VariableNames', {'Name', 'Estimate', 'SE', 'tStat', 'DF', 'pValue', 'Upper', 'Lower'});
    writetable(resultsTable, [savePath 'V1_DorsVSVent_lme.xlsx']);
end
