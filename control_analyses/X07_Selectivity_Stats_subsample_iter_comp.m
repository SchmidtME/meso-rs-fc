%% Comparison between 10 iterations of subsampling

close all
clear all
clc

% This script corresponds to Analysis B - The effect of beta, type, and layer 
% on rs-FC - and to Figure 3 of the manuscript. The data is averaged over 
% hemispheres and distances. An rm ANOVA is computed to assess the effects of beta, 
% layer, and type on rs-FC.

%% Load the data

Root = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results';

Labels = {'V1'}; 

Layers = {'0-2', '4-6', '8-10'};
Sbjs = {'aman', 'arak', 'aroo', 'atib', 'auil', 'chss', 'evad', 'haas', 'rcgr', 'ylri', 'imyy'};

subsample_iter = {'','_2', '_3', '_4', '_5', '_6', '_7', '_8', '_9'};

%%

for s = 1:length(subsample_iter)
    for Dist=1:10; %1:10%0:0.1:1
    
        for Label = 1:length(Labels)
            LayerCnt = 0;
            for Layer = 1:length(Layers)
                XX = Layers{Layer};
                LayerCnt = LayerCnt +1;
                String = [Labels{Label}, '_layers_', XX, '_intrahemispheric/'];
    
                for Sbj = 1:length(Sbjs)
    
                    load([Root,  '/', String, Sbjs{Sbj}, sprintf('/CorrelationMtx_Selectivity_subsampled%s.mat',subsample_iter{s})]);
                        
                    % 10     4     2    10    10     2
                    % Distance, Alike_both/alike_Eye1/alike_Eye2/Unalile (1), Hemi, Beta1, Beta2, z(2)
    
                    DATA_ALIKE(Sbj, Dist, Label, LayerCnt, :, :, :) = Data_Combined(Dist, 1, :, :, :, 2);
    
                    DATA_UNALIKE(Sbj, Dist, Label, LayerCnt, :, :, :) = Data_Combined(Dist, 4, :, :, :, 2);
    
                    DATA_DIFF(Sbj, Dist, Label, LayerCnt, :, :, :) = DATA_ALIKE(Sbj, Dist, Label, LayerCnt, :, :, :) - DATA_UNALIKE(Sbj, Dist, Label, LayerCnt, :, :, :);
    
                    clear Data_Combined
    
                end
    
            end
        end
    end
    
    % size: 11        10         1        1         2            10    10
    %       subjecs   distances  label   layer    hemispheres    betas
    
    % mean over distances, hemispheres for V1
    DATA_ALIKE_mean = squeeze(mean(mean(DATA_ALIKE(:,:,1,:,:,:,:), 5), 2));
    DATA_UNALIKE_mean = squeeze(mean(mean(DATA_UNALIKE(:,:,1,:,:,:,:), 5), 2));
    DATA(:,:,:,:,1) = DATA_ALIKE_mean;
    DATA(:,:,:,:,2) = DATA_UNALIKE_mean;
    
    % beta2beta diagonal
    
    STATS = [
        DATA_ALIKE_mean(:,1,1,1) DATA_ALIKE_mean(:,2,1,1) DATA_ALIKE_mean(:,3,1,1)...
        DATA_ALIKE_mean(:,1,2,2) DATA_ALIKE_mean(:,2,2,2) DATA_ALIKE_mean(:,3,2,2)...
        DATA_ALIKE_mean(:,1,3,3) DATA_ALIKE_mean(:,2,3,3) DATA_ALIKE_mean(:,3,3,3)...
        DATA_ALIKE_mean(:,1,4,4) DATA_ALIKE_mean(:,2,4,4) DATA_ALIKE_mean(:,3,4,4)...
        DATA_ALIKE_mean(:,1,5,5) DATA_ALIKE_mean(:,2,5,5) DATA_ALIKE_mean(:,3,5,5)...
        DATA_ALIKE_mean(:,1,6,6) DATA_ALIKE_mean(:,2,6,6) DATA_ALIKE_mean(:,3,6,6)...
        DATA_ALIKE_mean(:,1,7,7) DATA_ALIKE_mean(:,2,7,7) DATA_ALIKE_mean(:,3,7,7)...
        DATA_ALIKE_mean(:,1,8,8) DATA_ALIKE_mean(:,2,8,8) DATA_ALIKE_mean(:,3,8,8)...
        DATA_ALIKE_mean(:,1,9,9) DATA_ALIKE_mean(:,2,9,9) DATA_ALIKE_mean(:,3,9,9)...
        DATA_ALIKE_mean(:,1,10,10) DATA_ALIKE_mean(:,2,10,10) DATA_ALIKE_mean(:,3,10,10)...
        DATA_UNALIKE_mean(:,1,1,1) DATA_UNALIKE_mean(:,2,1,1) DATA_UNALIKE_mean(:,3,1,1)...
        DATA_UNALIKE_mean(:,1,2,2) DATA_UNALIKE_mean(:,2,2,2) DATA_UNALIKE_mean(:,3,2,2)...
        DATA_UNALIKE_mean(:,1,3,3) DATA_UNALIKE_mean(:,2,3,3) DATA_UNALIKE_mean(:,3,3,3)...
        DATA_UNALIKE_mean(:,1,4,4) DATA_UNALIKE_mean(:,2,4,4) DATA_UNALIKE_mean(:,3,4,4)...
        DATA_UNALIKE_mean(:,1,5,5) DATA_UNALIKE_mean(:,2,5,5) DATA_UNALIKE_mean(:,3,5,5)...
        DATA_UNALIKE_mean(:,1,6,6) DATA_UNALIKE_mean(:,2,6,6) DATA_UNALIKE_mean(:,3,6,6)...
        DATA_UNALIKE_mean(:,1,7,7) DATA_UNALIKE_mean(:,2,7,7) DATA_UNALIKE_mean(:,3,7,7)...
        DATA_UNALIKE_mean(:,1,8,8) DATA_UNALIKE_mean(:,2,8,8) DATA_UNALIKE_mean(:,3,8,8)...
        DATA_UNALIKE_mean(:,1,9,9) DATA_UNALIKE_mean(:,2,9,9) DATA_UNALIKE_mean(:,3,9,9)...
        DATA_UNALIKE_mean(:,1,10,10) DATA_UNALIKE_mean(:,2,10,10) DATA_UNALIKE_mean(:,3,10,10)
    ];
    
    ts = arrayfun(@(x) sprintf('t%d', x), 1:60, 'UniformOutput', false);
    
    t =array2table(STATS,'VariableNames', ts);
    
    within = table(...
        repmat({'A'; 'B'; 'C'}, 20, 1), ... % 30 repetitions of 'A', 'B', 'C'
        repmat({'AA'; 'AA'; 'AA'; 'BB'; 'BB'; 'BB'; 'CC'; 'CC'; 'CC'; 'DD'; 'DD'; 'DD'; 'EE'; 'EE'; 'EE'; ...
                'FF'; 'FF'; 'FF'; 'GG'; 'GG'; 'GG'; 'HH'; 'HH'; 'HH'; 'II'; 'II'; 'II'; 'JJ'; 'JJ'; 'JJ'}, 2, 1), ... % 2 repetitions of beta group
        [repmat({'AAA'}, 30, 1); repmat({'BBB'}, 30, 1)], ... % 'AAA' for first 30 rows, 'BBB' for next 30 rows
        'VariableNames', {'Layer', 'Beta', 'Type'});
    
    rm = fitrm(t,'t1-t60~1','WithinDesign',within)
    disp('Comparison between beta & layer & type in V1 - averaged over hemispheres, distances')
    ranovatbl = ranova(rm,'WithinModel','Layer+Beta+Type+Layer*Type+Layer*Beta+Beta*Layer+Type*Beta+Layer*Beta*Type')
    
    
    % LME
    % sizeInd = arrayfun(@(s) 1:s, size(DATA), 'UniformOutput', false);
    % [Subject, Layer, BQ1, BQ2, TypeODC] = ndgrid(sizeInd{:});
    % [Subject, Layer, TypeODC] = deal(categorical(Subject), categorical(Layer), categorical(TypeODC));
    % T = table(DATA(:), Subject(:), Layer(:), BQ1(:), BQ2(:), TypeODC(:), 'VariableNames', {'rsFC', 'Subject', 'Layer', 'BQ1', 'BQ2', 'TypeODC'});
    % %lme = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer + Layer*TypeODC + (1|Subject)'); % With most interactions (except BQ*TypeODC)
    % % lme = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)'); % With all 2nd order interactions
    % lme = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)'); % With 3nd order interactions
    % lme

    pause
end