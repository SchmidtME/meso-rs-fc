clc
close all
clear all

% Description:
% This script evaluates the difference in distance distributions between
% vertex pairs with alike or unalike ocular polarity before or after subsampling.
% Authors: Marianna Elisa Schmidt (marianna.schmidt@maxplanckschools.de)

%% Specifications

% Define the base directory where the subject folders are stored
baseDirectory = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results';

% save figure
saveFolder = fullfile('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Control_Analyses/Distances_ocular_polarity', datestr(now, 'yyyy-mm-dd'));
mkdir(saveFolder);

% subjects
subNames = {'aman', 'ylri', 'auil', 'arak', 'aroo', 'atib', 'imyy', 'chss', 'evad', 'haas', 'rcgr'};
layers = {'0-2'};
hemis = {'lh', 'rh'};

%% Load data

% Initialize cell array to store distances for 10x10 beta quantiles for all subjects
data_eye1eye1 = zeros(length(subNames),length(layers),length(hemis),10,10,10);
data_eye2eye2 = zeros(length(subNames),length(layers),length(hemis),10,10,10);
data_eye1eye2 = zeros(length(subNames),length(layers),length(hemis),10,10,10);
data_same_eye = zeros(length(subNames),length(layers),length(hemis),10,10,10);

for n = 1:length(subNames)
    for l = 1:length(layers)
        for h = 1:length(hemis)

        subjectDir = fullfile(baseDirectory, ['V1_layers_' layers{l} '_intrahemispheric'], subNames{n});
        matFilePath = fullfile(subjectDir, 'CorrelationMtx_Selectivity_Params_subsampled.mat');
        data = load(matFilePath, 'AnalysisParam');
        data_eye1eye1(n,l,h,:,:,:) = data.AnalysisParam.eye1eye1_median_quants_dist{h};
        data_eye2eye2(n,l,h,:,:,:) = data.AnalysisParam.eye2eye2_median_quants_dist{h};
        data_eye1eye2(n,l,h,:,:,:) = data.AnalysisParam.eye1eye2_median_quants_dist{h};
        data_same_eye(n,l,h,:,:,:) = data.AnalysisParam.same_eye_median_quants_dist{h};

        end
    end
end

%% Averaging

data_eye1eye1_mean = squeeze(mean(mean(mean(data_eye1eye1(:,:,:,:,:,1:10),3),4),5));
data_eye2eye2_mean = squeeze(mean(mean(mean(data_eye2eye2(:,:,:,:,:,1:10),3),4),5));
%data_same_eye_mean = (data_eye1eye1_mean+data_eye2eye2_mean)/2;
data_same_eye_mean = squeeze(mean(mean(mean(data_same_eye(:,:,:,:,:,1:10),3),4),5));
data_eye1eye2_mean = squeeze(mean(mean(mean(data_eye1eye2(:,:,:,:,:,1:10),3),4),5));

%% ANOVA for median distances in a quantile

STATS = [data_same_eye_mean(:,1:10) data_eye1eye2_mean(:,1:10)];

ts = arrayfun(@(x) sprintf('t%d', x), 1:20, 'UniformOutput', false);

t =array2table(STATS,'VariableNames', ts);

within = table(...
    [{'A'; 'A'; 'A'; 'A'; 'A'; 'A'; 'A'; 'A'; 'A'; 'A'; ...
    'B'; 'B'; 'B'; 'B'; 'B'; 'B'; 'B'; 'B'; 'B'; 'B'}], ... % effect of type
    [{'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; ...
    'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'}], ... % effect of quantile
    'VariableNames', {'Type', 'Quantile'});

rm = fitrm(t,'t1-t20~1','WithinDesign',within)
disp('Comparison of distances for vertices with alike or unalike ocular polarity')
ranovatbl = ranova(rm,'WithinModel','Type+Quantile+Type*Quantile')

%% plot median distances

% Calculate the mean over subjects for each quantile
mean_same_eye = mean(data_same_eye_mean, 1); % 1x10
mean_eye1eye2 = mean(data_eye1eye2_mean, 1); % 1x10

% Combine data for grouped bar plot
grouped_data = [mean_same_eye; mean_eye1eye2]';

% Create the grouped bar plot
figure;
bar(grouped_data, 'grouped');

% Customize the appearance
colormap([0.3 0.7 0.9; 0.9 0.3 0.3]); % Set custom colors (blue for same eye, red for eye1eye2)
legend({'Alike', 'Unalike'}, 'Location', 'NorthWest', 'FontSize', 12);
xlabel('Distance quantiles', 'FontSize', 14);
ylabel('Median distance (mm)', 'FontSize', 14);
title('Median distance for each distance quantile', 'FontSize', 16);

% Adjust x-axis labels
xticks(1:10); % Set ticks at 1 to 10
xticklabels(arrayfun(@num2str, 1:10, 'UniformOutput', false)); % Label quantiles
grid on;

savePath = fullfile(saveFolder, 'median_distances_quantiles.tiff');
saveas(gcf, savePath);

%% Archive

% %% Load data
% 
% % Initialize cell array to store distances for 10x10 beta quantiles for all subjects
% data_eye1eye1 = zeros(length(subNames),length(layers),length(hemis),10,10,11);
% data_eye2eye2 = zeros(length(subNames),length(layers),length(hemis),10,10,11);
% data_eye1eye2 = zeros(length(subNames),length(layers),length(hemis),10,10,11);
% 
% for n = 1:length(subNames)
%     for l = 1:length(layers)
%         for h = 1:length(hemis)
% 
%         subjectDir = fullfile(baseDirectory, ['V1_layers_' layers{l} '_intrahemispheric'], subNames{n});
%         matFilePath = fullfile(subjectDir, 'CorrelationMtx_Selectivity_Params.mat');
%         data = load(matFilePath, 'AnalysisParam');
%         data_eye1eye1(n,l,h,:,:,:) = data.AnalysisParam.eye1eye1_dist_thresholds{h};
%         data_eye2eye2(n,l,h,:,:,:) = data.AnalysisParam.eye2eye2_dist_thresholds{h};
%         data_eye1eye2(n,l,h,:,:,:) = data.AnalysisParam.eye1eye2_dist_thresholds{h};
% 
%         end
%     end
% end
% 
% %% Averaging
% 
% data_eye1eye1_mean = squeeze(mean(mean(mean(data_eye1eye1(:,:,:,:,:,1:10),3),4),5));
% data_eye2eye2_mean = squeeze(mean(mean(mean(data_eye2eye2(:,:,:,:,:,1:10),3),4),5));
% data_same_eye_mean = (data_eye1eye1_mean+data_eye2eye2_mean)/2;
% data_eye1eye2_mean = squeeze(mean(mean(mean(data_eye1eye2(:,:,:,:,:,1:10),3),4),5));
% 
% %% ANOVA (for shifted version)
% 
% STATS = [data_same_eye_mean(:,1:10) data_eye1eye2_mean(:,1:10)];
% 
% ts = arrayfun(@(x) sprintf('t%d', x), 1:20, 'UniformOutput', false);
% 
% t =array2table(STATS,'VariableNames', ts);
% 
% within = table(...
%     [{'A'; 'A'; 'A'; 'A'; 'A'; 'A'; 'A'; 'A'; 'A'; 'A'; ...
%     'B'; 'B'; 'B'; 'B'; 'B'; 'B'; 'B'; 'B'; 'B'; 'B'}], ... % effect of type
%     [{'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; ...
%     'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'}], ... % effect of quantile
%     'VariableNames', {'Type', 'Quantile'});
% 
% rm = fitrm(t,'t1-t20~1','WithinDesign',within)
% disp('Comparison of distances for vertices with alike or unalike ocular polarity')
% ranovatbl = ranova(rm,'WithinModel','Type+Quantile+Type*Quantile')
% 
% %%
% 
% for q = 1:10
% 
%     [h,p] = ttest(STATS(:,q), STATS(:,10+q))
%     pause
% 
% end
% 
% %%
% disp('Mean of distance between vertex pairs with alike and unalike ocular polarity')
% mean_STATS_same_eye = mean(STATS(:,1:10), 'all')
% mean_STATS_diff_eye = mean(STATS(:,11:20), 'all')
% diff = mean(STATS(:,11:20),1) - mean(STATS(:,1:10), 1)
% 
% % %% Load data
% % 
% % % Initialize cell array to store distances for 10x10 beta quantiles for all subjects
% % data_eye1eye1 = zeros(length(subNames),length(layers),length(hemis),10,10);
% % data_eye2eye2 = zeros(length(subNames),length(layers),length(hemis),10,10);
% % data_eye1eye2 = zeros(length(subNames),length(layers),length(hemis),10,10);
% % 
% % for n = 1:length(subNames)
% %     for l = 1:length(layers)
% %         for h = 1:length(hemis)
% % 
% %         subjectDir = fullfile(baseDirectory, ['V1_layers_' layers{l} '_intrahemispheric'], subNames{n});
% %         matFilePath = fullfile(subjectDir, 'CorrelationMtx_Selectivity_Params.mat');
% %         data = load(matFilePath, 'AnalysisParam');
% %         data_eye1eye1(n,l,h,:,:) = data.AnalysisParam.eye1eye1_mean_dist{h};
% %         data_eye2eye2(n,l,h,:,:) = data.AnalysisParam.eye2eye2_mean_dist{h};
% %         data_eye1eye2(n,l,h,:,:) = data.AnalysisParam.eye1eye2_mean_dist{h};
% % 
% %         end
% %     end
% % end
% % 
% % %% Averaging
% % 
% % data_eye1eye1_mean = squeeze(mean(mean(mean(data_eye1eye1,3),4),5));
% % data_eye2eye2_mean = squeeze(mean(mean(mean(data_eye2eye2,3),4),5));
% % data_same_eye_mean = (data_eye1eye1_mean+data_eye2eye2_mean)/2;
% % data_eye1eye2_mean = squeeze(mean(mean(mean(data_eye1eye2,3),4),5));
% % 
% % %% ANOVA (for shifted version)
% % 
% % STATS = [data_same_eye_mean data_eye1eye2_mean];
% % 
% % ts = arrayfun(@(x) sprintf('t%d', x), 1:2, 'UniformOutput', false);
% % 
% % t =array2table(STATS,'VariableNames', ts);
% % 
% % within = table(...
% %     [{'A'; 'B'}], ... % effect of type
% %     'VariableNames', {'Type'});
% % 
% % rm = fitrm(t,'t1-t2~1','WithinDesign',within)
% % disp('Comparison of distances for vertices with alike or unalike ocular polarity')
% % ranovatbl = ranova(rm,'WithinModel','Type')
% % 
% % %%
% % disp('Mean of distance between vertex pairs with alike and unalike ocular polarity')
% % mean_STATS = mean(STATS)