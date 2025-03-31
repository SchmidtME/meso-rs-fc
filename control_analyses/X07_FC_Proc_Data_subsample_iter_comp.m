close all
clear all
clc

% Descriptin:
% This script corresponds to Analysis A - The effect of distance (& type) on rs-FC and selectivity and to Figure 2 of the manuscript.
% For this version, I used data that was radially smoothed over layers 0-10, assuming that this is the latest version.
% The data is averaged over hemispheres and betas.
% First an rm ANOVA is computed to assess the effects of distance and type on rs-FC.
% Then 3 different rm ANOVAs are computed to assess the effect of distance on not normalized,
% weighted normalized, and unweighted normalized selectivity, respectively.

%% Load the data

% Define the root directory for the data
Root = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results';

% List of labels (regions of interest)
Labels = {'V1'}; 

% List of layers to process
Layers = {'0-10'};

% List of subjects
Sbjs = {'aman', 'arak', 'aroo', 'atib', 'auil', 'chss', 'evad', 'haas', 'rcgr', 'ylri', 'imyy'};

subsample_iter = {'','_2', '_3', '_4', '_5', '_6', '_7', '_8', '_9'};

%%
for s = 1:length(subsample_iter)
    % Loop through distances (1 to 10)
    for Dist=1:10
    
        % Loop through each label (region of interest)
        for Label = 1:length(Labels)
            LayerCnt = 0;
    
            % Loop through each layer
            for Layer = 1:length(Layers)
                % Define the current layer
                XX = Layers{Layer};
                LayerCnt = LayerCnt +1;
                % Construct the directory string for the current layer
                String = [Labels{Label}, '_layers_', XX, '_intrahemispheric/'];
    
                % Loop through each subject
                for Sbj = 1:length(Sbjs)
    
                    % Load the selectivity data for the subject
                    load([Root,  '/', String, Sbjs{Sbj},sprintf('/CorrelationMtx_Selectivity_subsampled%s.mat', subsample_iter{s})]);
                    
                    % Assign the data for the 'Alike' condition (both eyes)
                    DATA_ALIKE(Sbj, Dist, Label, LayerCnt, :, :, :) = Data_Combined(Dist, 1, :, :, :, 2);
    
                    % Assign the data for the 'Unalike' condition (different eyes)
                    DATA_UNALIKE(Sbj, Dist, Label, LayerCnt, :, :, :) = Data_Combined(Dist, 4, :, :, :, 2);
    
                    % Compute the difference between Alike and Unalike conditions
                    DATA_DIFF(Sbj, Dist, Label, LayerCnt, :, :, :) = DATA_ALIKE(Sbj, Dist, Label, LayerCnt, :, :, :) - DATA_UNALIKE(Sbj, Dist, Label, LayerCnt, :, :, :);
    
                    % Apply weighted normalization for the difference data
                    DATA_DIFF_NORM_weighted(Sbj, Dist, Label, LayerCnt, :, :, :) = 200*(DATA_ALIKE(Sbj, Dist, Label, LayerCnt, :, :, :) - DATA_UNALIKE(Sbj, Dist, Label, LayerCnt, :, :, :))./(DATA_ALIKE(Sbj, Dist, Label, LayerCnt, :, :, :) + DATA_UNALIKE(Sbj, Dist, Label, LayerCnt, :, :, :));
    
                    clear Data_Combined
    
                end
    
            end
        end
    end
    
    % Data size explanation:
    %   11        10         1        1         2            10    10
    %   subjects   distances  label   layer    hemispheres    betas
    
    % Calculate the mean over hemispheres and betas for 'Alike' and 'Unalike' conditions
    DATA_ALIKE_mean = squeeze(mean(mean(mean(DATA_ALIKE, 7), 6), 5));
    DATA_UNALIKE_mean = squeeze(mean(mean(mean(DATA_UNALIKE, 7), 6), 5));
    
    DATA(:,:,1) = DATA_ALIKE_mean;
    DATA(:,:,2) = DATA_UNALIKE_mean;
    
    %% Effects of distance and type on rs-FC
    
    % Combine the means of both conditions into one matrix for statistical analysis
    STATS = [DATA_ALIKE_mean DATA_UNALIKE_mean];
    
    % Create a table for the statistical data
    t = array2table(STATS, 'VariableNames', {'t1','t2', 't3', 't4', 't5', 't6', 't7', 't8', 't9', 't10', ...
        't11', 't12', 't13', 't14', 't15', 't16', 't17', 't18', 't19', 't20'});
    
    % Define the within-subject factors for the repeated measures ANOVA
    within = table({'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; 'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'},...
        {'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'AA'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'; 'BB'},...
        'VariableNames', {'Distance', 'Type'});
    
    % Fit the repeated measures model
    rm = fitrm(t, 't1-t20~1', 'WithinDesign', within);
    
    % Display the effects of distance and type on rs-FC
    disp('Effects of distance and type on rs-FC - averaged over betas and hemispheres');
    ranovatbl = ranova(rm, 'WithinModel', 'Distance+Type+Distance*Type')
    
    % LME
    % 
    % sizeInd = arrayfun(@(s) 1:s, size(DATA), 'UniformOutput', false);
    % [Subject, Dist, TypeODC] = ndgrid(sizeInd{:});
    % [Subject, TypeODC] = deal(categorical(Subject), categorical(TypeODC));
    % T = table(DATA(:), Subject(:), Dist(:), TypeODC(:), 'VariableNames', {'rsFC', 'Subject', 'Dist', 'TypeODC'});
    % lme = fitlme(T, 'rsFC ~ Dist*TypeODC + (1|Subject)'); % With most interactions (except BQ*TypeODC)
    % lme
    pause

end