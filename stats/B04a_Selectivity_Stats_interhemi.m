close all
clear all
clc

% This script corresponds to Analysis B - The effect of beta, type, and layer 
% on interhemispheric rs-FC - and to Figure 7 of the manuscript. An rm ANOVA 
% is computed to assess the effects of beta, layer, and type on rs-FC.

%% Load the data

Root = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results';

saveTables = 1;

currentDate = datestr(now, 'yyyy-mm-dd');
savePath = ['/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Figure_5/' currentDate '/'];
if saveTables
    mkdir(savePath);
end

Label = 'V1'; 

Layers = {'0-2', '4-6', '8-10'};
Sbjs = {'aman', 'arak', 'aroo', 'atib', 'auil', 'chss', 'evad', 'haas', 'rcgr', 'ylri', 'imyy'};

for Dist=1:10; %1:10%0:0.1:1

    for Layer = 1:length(Layers)

        String = [Label, '_layers_', Layers{Layer}, '_interhemispheric/'];

        for Sbj = 1:length(Sbjs)

            load([Root,  '/', String, Sbjs{Sbj},'/CorrelationMtx_Selectivity.mat']);
                
            % 6   10    10     2
            % alike(both eyes) eye1eye1 eye2eye2 diff eye (concat) eye1eye2 eye2eye1, Beta1, Beta2, z(2)

            DATA_ALIKE(Sbj, Layer, :, :) = Data_Combined(1, :, :, 2);

            DATA_UNALIKE(Sbj, Layer, :, :) = Data_Combined(4, :, :, 2);

            DATA_DIFF(Sbj, Layer, :, :) = DATA_ALIKE(Sbj, Layer, :, :) - DATA_UNALIKE(Sbj, Layer, :, :);

            clear Data_Combined

        end

    end

end

% size: 11        1        1       10    10
%       subjecs   label   layer    betas

% mean over distances, hemispheres not applicable
DATA_ALIKE_mean = squeeze(DATA_ALIKE);
DATA_UNALIKE_mean = squeeze(DATA_UNALIKE);


%% Stats for the effects of type, layers, and beta (accumulative) on rs-FC

STATS = [
    mean(mean(DATA_ALIKE_mean(:,1,1:2,1:2),3),4) mean(mean(DATA_ALIKE_mean(:,2,1:2,1:2),3),4) mean(mean(DATA_ALIKE_mean(:,3,1:2,1:2),3),4)...
    mean(mean(DATA_ALIKE_mean(:,1,1:4,1:4),3),4) mean(mean(DATA_ALIKE_mean(:,2,1:4,1:4),3),4) mean(mean(DATA_ALIKE_mean(:,3,1:4,1:4),3),4)...
    mean(mean(DATA_ALIKE_mean(:,1,1:6,1:6),3),4) mean(mean(DATA_ALIKE_mean(:,2,1:6,1:6),3),4) mean(mean(DATA_ALIKE_mean(:,3,1:6,1:6),3),4)...
    mean(mean(DATA_ALIKE_mean(:,1,1:8,1:8),3),4) mean(mean(DATA_ALIKE_mean(:,2,1:8,1:8),3),4) mean(mean(DATA_ALIKE_mean(:,3,1:8,1:8),3),4)...
    mean(mean(DATA_ALIKE_mean(:,1,1:10,1:10),3),4) mean(mean(DATA_ALIKE_mean(:,2,1:10,1:10),3),4) mean(mean(DATA_ALIKE_mean(:,3,1:10,1:10),3),4)...
    mean(mean(DATA_UNALIKE_mean(:,1,1:2,1:2),3),4) mean(mean(DATA_UNALIKE_mean(:,2,1:2,1:2),3),4) mean(mean(DATA_UNALIKE_mean(:,3,1:2,1:2),3),4)...
    mean(mean(DATA_UNALIKE_mean(:,1,1:4,1:4),3),4) mean(mean(DATA_UNALIKE_mean(:,2,1:4,1:4),3),4) mean(mean(DATA_UNALIKE_mean(:,3,1:4,1:4),3),4)...
    mean(mean(DATA_UNALIKE_mean(:,1,1:6,1:6),3),4) mean(mean(DATA_UNALIKE_mean(:,2,1:6,1:6),3),4) mean(mean(DATA_UNALIKE_mean(:,3,1:6,1:6),3),4)...
    mean(mean(DATA_UNALIKE_mean(:,1,1:8,1:8),3),4) mean(mean(DATA_UNALIKE_mean(:,2,1:8,1:8),3),4) mean(mean(DATA_UNALIKE_mean(:,3,1:8,1:8),3),4)...
    mean(mean(DATA_UNALIKE_mean(:,1,1:10,1:10),3),4) mean(mean(DATA_UNALIKE_mean(:,2,1:10,1:10),3),4) mean(mean(DATA_UNALIKE_mean(:,3,1:10,1:10),3),4)...
];

ts = arrayfun(@(x) sprintf('t%d', x), 1:30, 'UniformOutput', false);

t =array2table(STATS,'VariableNames', ts);

within = table(...
    repmat({'A'; 'B'; 'C'}, 10, 1), ... % effect of layer
    repmat({'AA'; 'AA'; 'AA'; 'BB'; 'BB'; 'BB'; 'CC'; 'CC'; 'CC'; 'DD'; 'DD'; 'DD'; 'EE'; 'EE'; 'EE'}, 2, 1), ... % effect of beta
    [repmat({'AAA'}, 15, 1); repmat({'BBB'}, 15, 1)], ... % effect of type
    'VariableNames', {'Layer', 'Beta', 'Type'});

rm = fitrm(t,'t1-t30~1','WithinDesign',within)
disp('Comparison between beta & layer & type for interhemispheric rs-FC')
ranovatbl = ranova(rm,'WithinModel','Layer+Beta+Type+Layer*Type+Layer*Beta+Beta*Layer+Beta*Type+Beta*Type*Layer')

if saveTables
    writetable(ranovatbl, [savePath 'V1_anova_cum.xlsx'], 'WriteRowNames', true);
end

%% LME, mean over distances, exclude lower triangle

clear DATA

% mean over distances
DATA_ALIKE_mean = squeeze(DATA_ALIKE);
DATA_UNALIKE_mean = squeeze(DATA_UNALIKE);

DATA(:,:,:,:,1) = DATA_UNALIKE_mean;
DATA(:,:,:,:,2) = DATA_ALIKE_mean;

sz = size(DATA);
Mask = true(sz) & permute(triu(true(sz(3:4))), [3 4 1 2 5]);

sizeInd = arrayfun(@(s) 1:s, size(DATA), 'UniformOutput', false);
[Subject, Layer, BQ1, BQ2, TypeODC] = ndgrid(sizeInd{:});
[Subject, Layer, TypeODC] = deal(categorical(Subject), categorical(Layer), categorical(TypeODC));
prodBQ = BQ1 .* BQ2;

T = table(DATA(Mask), Subject(Mask), Layer(Mask), prodBQ(Mask), TypeODC(Mask), 'VariableNames', {'rsFC', 'Subject', 'Layer', 'prodBQ', 'TypeODC'});

lme = fitlme(T, 'rsFC ~ prodBQ*Layer*TypeODC + (1|Subject)');

if saveTables
    [xx, xxx, Coefficients] = fixedEffects(lme, 'Alpha', 0.05);
    resultsTable = table(lme.CoefficientNames', Coefficients.Estimate, Coefficients.SE, Coefficients.tStat, Coefficients.DF, Coefficients.pValue, Coefficients.Upper, Coefficients.Lower, ...
                     'VariableNames', {'Name', 'Estimate', 'SE', 'tStat', 'DF', 'pValue', 'Upper', 'Lower'});
    writetable(resultsTable, [savePath 'V1_lme.xlsx']);
end