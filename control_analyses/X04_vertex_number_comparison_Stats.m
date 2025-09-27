clc
close all
clear all

%% Specifications

% Define the base directory where the subject folders are stored
baseDirectory = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results';
% subjects
subNames = {'myla', 'aman', 'ylri', 'auil', 'arak', 'aroo', 'atib', 'imyy', 'chss', 'evad', 'haas', 'rcgr'};
domEye = {'RE','LE', 'RE', 'RE', 'RE', 'LE', 'RE', 'RE', 'RE', 'RE', 'RE', 'RE'};
layers = {'0-2'};
hemis = {'lh', 'rh'};

saveFigures = 1;
if saveFigures
    currentDate = datestr(now, 'yyyy-mm-dd');
    saveDir = sprintf('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Supplementary/Figure_3/%s', currentDate);
    mkdir(saveDir);
end

%% Load data

% Initialize cell array to store distances for 10x10 beta quantiles for all subjects
data_DE = zeros(length(subNames),length(layers),length(hemis),1);
data_NDE = zeros(length(subNames),length(layers),length(hemis),1);

for n = 1:length(subNames)
    for l = 1:length(layers)
        for h = 1:length(hemis)

            subjectDir = fullfile(baseDirectory, ['V1_layers_' layers{l} '_intrahemispheric'], subNames{n});
            matFilePath = fullfile(subjectDir, 'CorrelationMtx_Selectivity_Params.mat');
            data = load(matFilePath, 'AnalysisParam');
            
            if strcmp(domEye{n}, 'RE')
                data_DE(n,l,h,:) = data.AnalysisParam.num_included_vertices_eye1{h};
                data_NDE(n,l,h,:) = data.AnalysisParam.num_included_vertices_eye2{h};
            else
                data_NDE(n,l,h,:) = data.AnalysisParam.num_included_vertices_eye1{h};
                data_DE(n,l,h,:) = data.AnalysisParam.num_included_vertices_eye2{h};
            end
        end
    end
end

%% Averaging

data_DE_mean = squeeze(mean(data_DE,3));
data_NDE_mean = squeeze(mean(data_NDE,3));

%% Stats

vertex_num = data_DE_mean + data_NDE_mean;
vertex_num_1000 = (1000./vertex_num)*100;

vertex_num_1000_mean = mean(vertex_num_1000);
vertex_num_1000_sd = std(vertex_num_1000);

%% ANOVA (for shifted version)

STATS = [data_DE_mean data_NDE_mean];

ts = arrayfun(@(x) sprintf('t%d', x), 1:2, 'UniformOutput', false);

t =array2table(STATS,'VariableNames', ts);

within = table(...
    [{'A'; 'B'}], ... % effect of type
    'VariableNames', {'Type'});

rm = fitrm(t,'t1-t2~1','WithinDesign',within)
disp('Comparison of number of vertices in the masks for the DE and NDE')
ranovatbl = ranova(rm,'WithinModel','Type')

%%

[h,p] = ttest(STATS(:,1),STATS(:,2))

%%

figure('Position', [100, 100, 400, 400]); % Adjust figure size (width 600, height 400)
hold on;
set(gcf, 'Color', 'w'); % Set background to white

% Define an offset to move points slightly to the right
offset = 0.05;

% Adjusted x-positions for better spacing
x_DE = 1 + offset;  
x_NDE = 1.6 + offset;  

% Scatter plots with black outlines
scatter(x_DE * ones(12, 1), data_DE_mean, 100, 'k', 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'LineWidth', 3); % Black outline for DE
scatter(x_NDE * ones(12, 1), data_NDE_mean, 100, [0.5 0.5 0.5], '^', 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 3); % Grey outline for NDE

% Plot lines connecting corresponding data points
for i = 1:length(data_DE_mean)
    plot([x_DE x_NDE], [data_DE_mean(i) data_NDE_mean(i)], 'k-', 'LineWidth', 1);
end

% Adjust the axes and labels
xlim([0.8 2]); % Keep space to the left
ylim([1500 10500]);

% Set y-ticks to just the min and max values
set(gca, 'YTick', [1500 10500]);

% Corrected x-ticks for the new positions
set(gca, 'XTick', [x_DE x_NDE], 'XTickLabel', {'DE', 'NDE'});
set(gca, 'LineWidth', 2); % Increase axis line width
set(gcf,'renderer', 'painters'); 

xlabel('Eye dominance');
ylabel('Vertex count');
%title('Comparison of vertex counts for DE and NDE');

hold off;

if saveFigures
    print([saveDir '/Figure_3'], '-dtiff', '-r400');
end
