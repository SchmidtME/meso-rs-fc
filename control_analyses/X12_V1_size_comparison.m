clc
close all
clear all

%% Specifications

% Define the base directory where the subject folders are stored
baseDirectory = '/space/ardebil/1/users/Others/Marianna/good_subjects_anat/retinotopy';
% subjects
subNames = {'aman', 'ylri', 'auil', 'arak', 'aroo', 'atib', 'imyy', 'chss', 'evad', 'haas', 'rcgr'};

hemis = {'lh', 'rh'};

saveFigures = 0;
if saveFigures
    currentDate = datestr(now, 'yyyy-mm-dd');
    saveDir = sprintf('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Figure_10/%s', currentDate);
    mkdir(saveDir);
end

%% Load data

% Initialize cell array to store distances for 10x10 beta quantiles for all subjects
data_V1 = zeros(length(subNames), length(hemis),1);

for n = 1:length(subNames)
    for h = 1:length(hemis)

        subjectDir = fullfile(baseDirectory, subNames{n});
        roi_label = read_ROIlabel(fullfile(subjectDir, [hemis{h} '.V1_Upsampled_Rtopy.label']));
        data_V1(n,h) = (1000/size(roi_label,1))*100;

    end
end

%%

data = round(mean(data_V1,2));

% Number of bars
numBars = length(data);

% Create a colormap with rainbow colors
cmap = jet(numBars); % 'jet' provides a rainbow-like color scheme

% Create figure
figure;
hold on;

% Plot bars with different colors
for i = 1:numBars
    barHandle = bar(i, data(i), 'FaceColor', cmap(i, :)); % Assign color
end

% Add value labels on top of bars
for i = 1:numBars
    text(i, data(i) + 0.5, sprintf('%.f', data(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
end

% Formatting
xticks(1:numBars);
xticklabels(arrayfun(@num2str, 1:numBars, 'UniformOutput', false)); % Label bars with index numbers
xlabel('Subject');
ylim([0 16])
ylabel('% of V1 label');
title('% of V1 label corresponding to 1000 vertices');
colormap(cmap); % Apply colormap

hold off;

%%

% Calculate Mean
mean_value = mean(data);

% Calculate SEM (Standard Error of the Mean)
STD_value = std(data); %/ sqrt(length(data));

% Display Results
fprintf('Mean: %.3f\n', mean_value);
fprintf('STD: %.3f\n', STD_value);
% if saveFigures
%     print([saveDir '/' 'Figure_10'], '-dtiff', '-r200');
% end
