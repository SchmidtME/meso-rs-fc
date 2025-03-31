close all
clear all
clc

% Description:
% This script corresponds to Analysis C - the robustness/replicability of rs-FC
% and selectivity and to Figures suppl xxx of the manuscript.
% This script assesses the main effects of session, type and distance as well
% as respective interaction effects on rs-FC and therefore is a test for the 
% robustness or replicability of the rs-FC measure and analysis.

%% Load data

saveTables = 1;
currentDate = datestr(now, 'yyyy-mm-dd');
savePath = ['/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Control_Analyses/odc_sep_layers/' currentDate '/'];
if saveTables 
    mkdir(savePath);
end

Root = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results';

% Specify group labels and subject names for each group
GroupLabels = {'intrahemispheric', 'intrahemispheric_odc_sep_depths'};
Sbjs = {'aman', 'ylri', 'auil', 'aroo', 'atib', 'imyy', 'chss', 'evad', 'haas', 'rcgr'}; % arak just one run, excluded


% Define labels (e.g., brain regions) and layers for analysis
Labels = {'V1'};
Layers = {'0-2', '4-6', '8-10'};

%%

% Loop through groups, distances, labels, and layers to load and organize data
for Group = 1:length(GroupLabels)
    for Dist = 1:10
        for Label = 1:length(Labels)
            LayerCnt = 0;
            for Layer = 1:length(Layers)
                XX = Layers{Layer};
                LayerCnt = LayerCnt + 1;

                % Construct the path to the data based on group and layer
                String = [Labels{Label}, '_layers_', XX, '_' GroupLabels{Group} '/'];

                % Load correlation matrices for each subject
                for Sbj = 1:length(Sbjs)
                    
                    if Group ==1
                        load([Root, '/', String, Sbjs{Sbj}, '/CorrelationMtx_Selectivity_subsampled.mat']);
                    else
                        load([Root, '/', String, Sbjs{Sbj}, '/CorrelationMtx_Selectivity_subsampled_1.mat']);
                    end
                    % Extract data for alike and unalike conditions and organize into arrays
                    DATA_ALIKE{Group}(Sbj, Dist, Label, LayerCnt, :, :, :) = Data_Combined(Dist, 1, :, :, :, 2);
                    DATA_UNALIKE{Group}(Sbj, Dist, Label, LayerCnt, :, :, :) = Data_Combined(Dist, 4, :, :, :, 2);

                    clear Data_Combined
                end
            end
        end
    end
end

% Data dimensions: 
% 11 subjects, 10 distances, 5 labels, 3 layers, 2 hemispheres, 2 betas

%% Effect of distance and type
% Calculate mean values across layers, hemispheres, and betas for V1
DATA_ALIKE_mean_run1 = squeeze(mean(mean(mean(mean(DATA_ALIKE{1}(:,:,1,:,:,:,:), 7), 6), 5), 4));
DATA_UNALIKE_mean_run1 = squeeze(mean(mean(mean(mean(DATA_UNALIKE{1}(:,:,1,:,:,:,:), 7), 6), 5), 4));
STATS_run1 = [DATA_ALIKE_mean_run1 DATA_UNALIKE_mean_run1];

DATA_ALIKE_mean_run2 = squeeze(mean(mean(mean(mean(DATA_ALIKE{2}(:,:,1,:,:,:,:), 7), 6), 5), 4));
DATA_UNALIKE_mean_run2 = squeeze(mean(mean(mean(mean(DATA_UNALIKE{2}(:,:,1,:,:,:,:), 7), 6), 5), 4));
STATS_run2 = [DATA_ALIKE_mean_run2 DATA_UNALIKE_mean_run2];

% Combine statistics across sessions
STATS = [STATS_run1 STATS_run2];

% Define table for repeated measures ANOVA
ts = arrayfun(@(x) sprintf('t%d', x), 1:40, 'UniformOutput', false);
t = array2table(STATS(:, :), 'VariableNames', ts');
within = table(...
    {'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; 'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; 'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; 'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'}, ... % Distance
    {'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'}, ... % Type
    {'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'}, ... % Run
    'VariableNames', {'Distance', 'Type', 'Run'});

% Perform repeated measures ANOVA
rm = fitrm(t, 't1-t40~1', 'WithinDesign', within);
disp('Comparison between distance in V1 - averaged over betas, hemispheres, layers');
ranovatbl = ranova(rm, 'WithinModel', 'Distance+Type+Run+Distance*Type+Distance*Run+Type*Run+Type*Distance*Run')

if saveTables
    writetable(ranovatbl, [savePath 'V1_anova_justDist.xlsx'], 'WriteRowNames', true);
end