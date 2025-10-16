close all
clear all
%clc

% Description:
% This script corresponds to Analysis C - the robustness/replicability of rs-FC
% and selectivity and to Figures S3 of the manuscript.
% This script assesses the main effects of session, type and distance as well
% as respective interaction effects on rs-FC and therefore is a test for the 
% robustness or replicability of the rs-FC measure and analysis.
% Authors: Marianna E. Schmidt (marianna.schmidt@maxplanckschools.de), Iman Aganj, Shahin Nasr

%% Load data
% Define root directories for the new and old control groups
Root_Control_new = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls_2/Results';
Root_Control_old = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results';

% Specify group labels and subject names for each group
GroupLabels = {'Control_old', 'Control_new'};
Sbjs_Control_new = {'aman', 'aroo', 'auil', 'chss','haas'};
Sbjs_Control_old = {'aman', 'aroo', 'auil', 'chss','haas'};

% Define labels (e.g., brain regions) and layers for analysis
Labels = {'V1'};
Layers = {'0-10'};

% Organize groups and roots into cell arrays for iteration
Groups = {Sbjs_Control_old, Sbjs_Control_new};
Roots = {Root_Control_old, Root_Control_new};

subsample_iter = {'_1'};

%%
for s=1:length(subsample_iter)
    % Loop through groups, distances, labels, and layers to load and organize data
    for Group = 1:length(GroupLabels)
        for Dist = 1:10
            for Label = 1:length(Labels)
                LayerCnt = 0;
                for Layer = 1:length(Layers)
                    XX = Layers{Layer};
                    LayerCnt = LayerCnt + 1;
    
                    % Construct the path to the data based on group and layer
                    String = [Labels{Label}, '_layers_', XX, '_intrahemispheric/'];
    
                    % Load correlation matrices for each subject
                    for Sbj = 1:length(Groups{Group})
    
                        load([Roots{Group}, '/', String, Groups{Group}{Sbj}, sprintf('/CorrelationMtx_Selectivity_subsampled%s.mat',subsample_iter{s})]);
    
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
    DATA_ALIKE_mean_sess1 = squeeze(mean(mean(mean(mean(DATA_ALIKE{1}(:,:,1,:,:,:,:), 7), 6), 5), 4));
    DATA_UNALIKE_mean_sess1 = squeeze(mean(mean(mean(mean(DATA_UNALIKE{1}(:,:,1,:,:,:,:), 7), 6), 5), 4));
    STATS_sess1 = [DATA_ALIKE_mean_sess1 DATA_UNALIKE_mean_sess1];
    
    DATA_ALIKE_mean_sess2 = squeeze(mean(mean(mean(mean(DATA_ALIKE{2}(:,:,1,:,:,:,:), 7), 6), 5), 4));
    DATA_UNALIKE_mean_sess2 = squeeze(mean(mean(mean(mean(DATA_UNALIKE{2}(:,:,1,:,:,:,:), 7), 6), 5), 4));
    STATS_sess2 = [DATA_ALIKE_mean_sess2 DATA_UNALIKE_mean_sess2];
    
    % Combine statistics across sessions
    STATS = [STATS_sess1 STATS_sess2];
    
    % Define table for repeated measures ANOVA
    ts = arrayfun(@(x) sprintf('t%d', x), 1:40, 'UniformOutput', false);
    t = array2table(STATS(:, :), 'VariableNames', ts');
    within = table(...
        {'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; 'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; 'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; 'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'}, ... % Distance
        {'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'}, ... % Type
        {'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'AAA'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'; 'BBB'}, ... % Session
        'VariableNames', {'Distance', 'Type', 'Session'});
    
    % Perform repeated measures ANOVA
    rm = fitrm(t, 't1-t40~1', 'WithinDesign', within);
    disp('Comparison between distance in V1 - averaged over betas, hemispheres, layers');
    ranovatbl = ranova(rm, 'WithinModel', 'Distance+Type+Session+Distance*Type+Distance*Session+Type*Session+Type*Distance*Session')
    
    %% Correlation between sessions
    % Calculate correlations between alike/unlike data across sessions
    for i = 1:size(STATS, 1)
        correlation(i) = corr(STATS(i, 1:20)', STATS(i, 21:40)');
        
        % Plot correlations between sessions
        figure(i); subplot(1, 2, 1);
        plot(STATS(i, 1:10)', STATS(i, 21:30)', 'r*'); hold on;
        plot(STATS(i, 11:20)', STATS(i, 31:40)', 'b*'); hold off;
    
        % Plot differences and fit a line
        figure(i); subplot(1, 2, 2);
        correlation2(i) = corr(STATS(i, 1:10)' - STATS(i, 11:20)', STATS(i, 21:30)' - STATS(i, 31:40)');
        plot(STATS(i, 1:10)' - STATS(i, 11:20)', STATS(i, 21:30)' - STATS(i, 31:40)', 'r*'); hold on;
        title(num2str(correlation2(i)));
        lsline;
        pause;
    end
end
