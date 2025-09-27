close all
clear all
clc

% Description:
% This script corresponds to Analysis A - The effect of distance (& type) on 
% rs-FC and selectivity and to Figure 2 of the manuscript.
% For this version, I used data that was radially smoothed over layers 0-10, 
% with a data matrix that was subsampled to match distance distributions of
% vertex pairs with alike and unalike relative ocular polairty.
% The dataa is averaged over hemispheres and betas.
% An rm ANOVA is computed to assess the effects of distance and type on rs-FC.

%% Load the data

% Define the root directory for the data
Root = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results';

% List of labels (regions of interest)
Labels = {'V1'}; 

% List of layers to process
Layers = {'0-10'};

% List of subjects
Sbjs = {'myla','aman', 'arak', 'aroo', 'atib', 'auil', 'chss', 'evad', 'haas', 'rcgr', 'ylri', 'imyy'};

% Loop through distances (1 to 10)
for Dist=1:10

    % Loop through each label (region of interest)
    for Label = 1:length(Labels)

        % Loop through each layer
        for Layer = 1:length(Layers)

            % Construct the directory string for the current layer
            String = [Labels{Label}, '_layers_', Layers{Layer}, '_intrahemispheric/'];

            % Loop through each subject
            for Sbj = 1:length(Sbjs)

                % Load the selectivity data for the subject
                load([Root,  '/', String, Sbjs{Sbj},'/CorrelationMtx_Selectivity_subsampled_1.mat']);
                
                % Assign the data for the 'Alike' condition (both eyes)
                DATA_ALIKE(Sbj, Dist, Label, Layer, :, :, :) = Data_Combined(Dist, 1, :, :, :, 2);

                % Assign the data for the 'Unalike' condition (different eyes)
                DATA_UNALIKE(Sbj, Dist, Label, Layer, :, :, :) = Data_Combined(Dist, 4, :, :, :, 2);

                % Compute the difference between Alike and Unalike conditions
                DATA_DIFF(Sbj, Dist, Label, Layer, :, :, :) = DATA_ALIKE(Sbj, Dist, Label, Layer, :, :, :) - DATA_UNALIKE(Sbj, Dist, Label, Layer, :, :, :);

                clear Data_Combined

            end

        end
    end
end

% Data size explanation:
%   11        10         1        1         2            10    10
%   subjects   distances  label   layer    hemispheres    betas


%% Effects of distance and type on rs-FC

% Calculate the mean over hemispheres and betas for 'Alike' and 'Unalike' conditions
DATA_ALIKE_mean = squeeze(mean(mean(mean(DATA_ALIKE, 7), 6), 5));
DATA_UNALIKE_mean = squeeze(mean(mean(mean(DATA_UNALIKE, 7), 6), 5));

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

%% LME

% Calculate the mean over hemispheres and betas for 'Alike' and 'Unalike' conditions
DATA_ALIKE_mean = squeeze(mean(mean(mean(DATA_ALIKE, 7), 6), 5));
DATA_UNALIKE_mean = squeeze(mean(mean(mean(DATA_UNALIKE, 7), 6), 5));

DATA(:,:,1) = DATA_ALIKE_mean;
DATA(:,:,2) = DATA_UNALIKE_mean;

sizeInd = arrayfun(@(s) 1:s, size(DATA), 'UniformOutput', false);
[Subject, Dist, TypeODC] = ndgrid(sizeInd{:});
[Subject, Dist, TypeODC] = deal(categorical(Subject), categorical(Dist), categorical(TypeODC));
T = table(DATA(:), Subject(:), Dist(:), TypeODC(:), 'VariableNames', {'rsFC', 'Subject', 'Dist', 'TypeODC'});
lme = fitlme(T, 'rsFC ~ Dist*TypeODC + (1|Subject)');
lme
