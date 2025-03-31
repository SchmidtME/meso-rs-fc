close all
clear all
clc

% Description:
% This script corresponds to Analysis B - The effect of beta, type, and layer 
% on rs-FC - and to Figure 3 of the manuscript. The data is averaged over 
% hemispheres and distances. 
% First, an rm ANOVA is computed to assess the effects of beta squares (1-3 to
% 1-3, 4-7 to 4-7, 8-10 to 8-10).
% Then an rmANOVA is computed to assess the effects of beta diagonal (1 to 1, ...)
% Third, an rmAONVA is computed to assess the effects of cumulative beta (1-2 to 1-2,
% 1-4 to 1-4, 1-6 to 1-6, 1-8 to 1-8, 1-10 to 1-10).
% Lastly, an LME is fitted to account for the hierarchical structure of the data.

%% Load the data

Root = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results/Supplementary';

saveTables = 1;
currentDate = datestr(now, 'yyyy-mm-dd');
savePath = ['/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Supplementary/Figure_3/' currentDate '/'];
if saveTables 
    mkdir(savePath);
end

Label = 'V1';
Layers = {'0-2', '4-6', '8-10'};
Sbjs = {'aman', 'arak', 'aroo', 'atib', 'auil', 'chss', 'evad', 'haas', 'rcgr', 'ylri', 'imyy'};

Dist=1:10;

for Layer = 1:length(Layers)

    String = [Label, '_layers_', Layers{Layer}, '_intrahemispheric_vis_resp/'];

    for Sbj = 1:length(Sbjs)

        load([Root,  '/', String, Sbjs{Sbj},'/CorrelationMtx_Selectivity_subsampled_1.mat']);
            
        % 10     4     2    10    10     2
        % Distance, Alike_both/alike_Eye1/alike_Eye2/Unalile (1), Hemi, Beta1, Beta2, z(2)

        DATA_ALIKE(Sbj, Dist, Layer, :, :, :) = Data_Combined(Dist, 1, :, :, :, 2);

        DATA_UNALIKE(Sbj, Dist, Layer, :, :, :) = Data_Combined(Dist, 4, :, :, :, 2);

        DATA_DIFF(Sbj, Dist, Layer, :, :, :) = DATA_ALIKE(Sbj, Dist, Layer, :, :, :) - DATA_UNALIKE(Sbj, Dist, Layer, :, :, :);

        clear Data_Combined

    end

end


% size: 11        10          3         2            10    10
%       subjecs   distances   layers    hemispheres  beta1 beta2

% mean over distances, hemispheres for V1
DATA_ALIKE_mean = squeeze(mean(mean(DATA_ALIKE, 4), 2));
DATA_UNALIKE_mean = squeeze(mean(mean(DATA_UNALIKE, 4), 2));


%% rmANOVA for cumulative beta

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
    repmat({'A'; 'B'; 'C'}, 10, 1), ... % 30 repetitions of 'A', 'B', 'C'
    repmat({'AA'; 'AA'; 'AA'; 'BB'; 'BB'; 'BB'; 'CC'; 'CC'; 'CC'; 'DD'; 'DD'; 'DD'; 'EE'; 'EE'; 'EE'}, 2, 1), ... % 2 repetitions of beta group
    [repmat({'AAA'}, 15, 1); repmat({'BBB'}, 15, 1)], ... % 'AAA' for first 30 rows, 'BBB' for next 30 rows
    'VariableNames', {'Layer', 'Beta', 'Type'});

rm = fitrm(t,'t1-t30~1','WithinDesign',within)
disp('Comparison between beta (cumulative) & layer & type in V1 - averaged over hemispheres, distances')
ranovatbl = ranova(rm,'WithinModel','Layer+Beta+Type+Layer*Type+Layer*Beta+Beta*Layer+Beta*Type+Beta*Type*Layer')

if saveTables
    writetable(ranovatbl, [savePath 'V1_anova_cum.xlsx'], 'WriteRowNames', true);
end

%% LME, mean over distances, exclude lower triangle

% mean over distances
DATA_ALIKE_mean = squeeze(mean(DATA_ALIKE, 2));
DATA_UNALIKE_mean = squeeze(mean(DATA_UNALIKE, 2));

DATA(:,:,:,:,:,1) = DATA_UNALIKE_mean;
DATA(:,:,:,:,:,2) = DATA_ALIKE_mean;

sz = size(DATA);
Mask = true(sz) & permute(triu(true(sz(4:5))), [3 4 5 1 2 6]);

sizeInd = arrayfun(@(s) 1:s, size(DATA), 'UniformOutput', false);
[Subject, Layer, Hemi, BQ1, BQ2, TypeODC] = ndgrid(sizeInd{:});
[Subject, Layer, Hemi, TypeODC] = deal(categorical(Subject), categorical(Layer), categorical(Hemi), categorical(TypeODC));
prodBQ = BQ1 .* BQ2;

T = table(DATA(Mask), Subject(Mask), Layer(Mask), Hemi(Mask), prodBQ(Mask), TypeODC(Mask), 'VariableNames', {'rsFC', 'Subject', 'Layer', 'Hemi', 'prodBQ', 'TypeODC'});

lme = fitlme(T, 'rsFC ~ prodBQ*Layer*TypeODC + (1|Subject)');

if saveTables
    [xx, xxx, Coefficients] = fixedEffects(lme, 'Alpha', 0.05);
    resultsTable = table(lme.CoefficientNames', Coefficients.Estimate, Coefficients.SE, Coefficients.tStat, Coefficients.DF, Coefficients.pValue, Coefficients.Upper, Coefficients.Lower, ...
                     'VariableNames', {'Name', 'Estimate', 'SE', 'tStat', 'DF', 'pValue', 'Upper', 'Lower'});
    writetable(resultsTable, [savePath 'V1_lme.xlsx']);
end

%% Archive
% 
% %% rmANOVA for beta squares
% 
% STATS = [
%     mean(mean(DATA_ALIKE_mean(:,1,1:3,1:3),3),4) mean(mean(DATA_ALIKE_mean(:,2,1:3,1:3),3),4) mean(mean(DATA_ALIKE_mean(:,3,1:3,1:3),3),4)...
%     mean(mean(DATA_ALIKE_mean(:,1,4:7,4:7),3),4) mean(mean(DATA_ALIKE_mean(:,2,4:7,4:7),3),4) mean(mean(DATA_ALIKE_mean(:,3,4:7,4:7),3),4)...
%     mean(mean(DATA_ALIKE_mean(:,1,8:10,8:10),3),4) mean(mean(DATA_ALIKE_mean(:,2,8:10,8:10),3),4) mean(mean(DATA_ALIKE_mean(:,3,8:10,8:10),3),4)...
%     mean(mean(DATA_UNALIKE_mean(:,1,1:3,1:3),3),4) mean(mean(DATA_UNALIKE_mean(:,2,1:3,1:3),3),4) mean(mean(DATA_UNALIKE_mean(:,3,1:3,1:3),3),4)...
%     mean(mean(DATA_UNALIKE_mean(:,1,4:7,4:7),3),4) mean(mean(DATA_UNALIKE_mean(:,2,4:7,4:7),3),4) mean(mean(DATA_UNALIKE_mean(:,3,4:7,4:7),3),4)...
%     mean(mean(DATA_UNALIKE_mean(:,1,8:10,8:10),3),4) mean(mean(DATA_UNALIKE_mean(:,2,8:10,8:10),3),4) mean(mean(DATA_UNALIKE_mean(:,3,8:10,8:10),3),4)...
% ];
% 
% ts = arrayfun(@(x) sprintf('t%d', x), 1:18, 'UniformOutput', false);
% 
% t =array2table(STATS,'VariableNames', ts);
% 
% within = table(...
%     repmat({'A'; 'B'; 'C'}, 6, 1), ... % 30 repetitions of 'A', 'B', 'C'
%     repmat({'AA'; 'AA'; 'AA'; 'BB'; 'BB'; 'BB'; 'CC'; 'CC'; 'CC'}, 2, 1), ... % 2 repetitions of beta group
%     [repmat({'AAA'}, 9, 1); repmat({'BBB'}, 9, 1)], ... % 'AAA' for first 30 rows, 'BBB' for next 30 rows
%     'VariableNames', {'Layer', 'Beta', 'Type'});
% 
% rm = fitrm(t,'t1-t18~1','WithinDesign',within)
% disp('Comparison between beta (squares) & layer & type in V1 - averaged over hemispheres, distances')
% ranovatbl = ranova(rm,'WithinModel','Layer+Beta+Type+Layer*Type+Layer*Beta+Beta*Layer+Beta*Type+Beta*Type*Layer')
% 
% %% rmANOVA for beta diagonal
% 
% STATS = [
%     DATA_ALIKE_mean(:,1,1,1) DATA_ALIKE_mean(:,2,1,1) DATA_ALIKE_mean(:,3,1,1)...
%     DATA_ALIKE_mean(:,1,2,2) DATA_ALIKE_mean(:,2,2,2) DATA_ALIKE_mean(:,3,2,2)...
%     DATA_ALIKE_mean(:,1,3,3) DATA_ALIKE_mean(:,2,3,3) DATA_ALIKE_mean(:,3,3,3)...
%     DATA_ALIKE_mean(:,1,4,4) DATA_ALIKE_mean(:,2,4,4) DATA_ALIKE_mean(:,3,4,4)...
%     DATA_ALIKE_mean(:,1,5,5) DATA_ALIKE_mean(:,2,5,5) DATA_ALIKE_mean(:,3,5,5)...
%     DATA_ALIKE_mean(:,1,6,6) DATA_ALIKE_mean(:,2,6,6) DATA_ALIKE_mean(:,3,6,6)...
%     DATA_ALIKE_mean(:,1,7,7) DATA_ALIKE_mean(:,2,7,7) DATA_ALIKE_mean(:,3,7,7)...
%     DATA_ALIKE_mean(:,1,8,8) DATA_ALIKE_mean(:,2,8,8) DATA_ALIKE_mean(:,3,8,8)...
%     DATA_ALIKE_mean(:,1,9,9) DATA_ALIKE_mean(:,2,9,9) DATA_ALIKE_mean(:,3,9,9)...
%     DATA_ALIKE_mean(:,1,10,10) DATA_ALIKE_mean(:,2,10,10) DATA_ALIKE_mean(:,3,10,10)...
%     DATA_UNALIKE_mean(:,1,1,1) DATA_UNALIKE_mean(:,2,1,1) DATA_UNALIKE_mean(:,3,1,1)...
%     DATA_UNALIKE_mean(:,1,2,2) DATA_UNALIKE_mean(:,2,2,2) DATA_UNALIKE_mean(:,3,2,2)...
%     DATA_UNALIKE_mean(:,1,3,3) DATA_UNALIKE_mean(:,2,3,3) DATA_UNALIKE_mean(:,3,3,3)...
%     DATA_UNALIKE_mean(:,1,4,4) DATA_UNALIKE_mean(:,2,4,4) DATA_UNALIKE_mean(:,3,4,4)...
%     DATA_UNALIKE_mean(:,1,5,5) DATA_UNALIKE_mean(:,2,5,5) DATA_UNALIKE_mean(:,3,5,5)...
%     DATA_UNALIKE_mean(:,1,6,6) DATA_UNALIKE_mean(:,2,6,6) DATA_UNALIKE_mean(:,3,6,6)...
%     DATA_UNALIKE_mean(:,1,7,7) DATA_UNALIKE_mean(:,2,7,7) DATA_UNALIKE_mean(:,3,7,7)...
%     DATA_UNALIKE_mean(:,1,8,8) DATA_UNALIKE_mean(:,2,8,8) DATA_UNALIKE_mean(:,3,8,8)...
%     DATA_UNALIKE_mean(:,1,9,9) DATA_UNALIKE_mean(:,2,9,9) DATA_UNALIKE_mean(:,3,9,9)...
%     DATA_UNALIKE_mean(:,1,10,10) DATA_UNALIKE_mean(:,2,10,10) DATA_UNALIKE_mean(:,3,10,10)
% ];
% 
% ts = arrayfun(@(x) sprintf('t%d', x), 1:60, 'UniformOutput', false);
% 
% t =array2table(STATS,'VariableNames', ts);
% 
% within = table(...
%     repmat({'A'; 'B'; 'C'}, 20, 1), ... % 30 repetitions of 'A', 'B', 'C'
%     repmat({'AA'; 'AA'; 'AA'; 'BB'; 'BB'; 'BB'; 'CC'; 'CC'; 'CC'; 'DD'; 'DD'; 'DD'; 'EE'; 'EE'; 'EE'; ...
%             'FF'; 'FF'; 'FF'; 'GG'; 'GG'; 'GG'; 'HH'; 'HH'; 'HH'; 'II'; 'II'; 'II'; 'JJ'; 'JJ'; 'JJ'}, 2, 1), ... % 2 repetitions of beta group
%     [repmat({'AAA'}, 30, 1); repmat({'BBB'}, 30, 1)], ... % 'AAA' for first 30 rows, 'BBB' for next 30 rows
%     'VariableNames', {'Layer', 'Beta', 'Type'});
% 
% rm = fitrm(t,'t1-t60~1','WithinDesign',within)
% disp('Comparison between beta (diagonal) & layer & type in V1 - averaged over hemispheres, distances')
% ranovatbl = ranova(rm,'WithinModel','Layer+Beta+Type+Layer*Type+Layer*Beta+Beta*Layer+Type*Beta+Layer*Beta*Type')
% 
% 
% %% LME, mean over hemispheres and distances
% clear DATA
% % mean over distances, hemispheres for V1
% DATA_ALIKE_mean = squeeze(mean(mean(DATA_ALIKE, 4), 2));
% DATA_UNALIKE_mean = squeeze(mean(mean(DATA_UNALIKE, 4), 2));
% 
% DATA(:,:,:,:,1) = DATA_ALIKE_mean;
% DATA(:,:,:,:,2) = DATA_UNALIKE_mean;
% 
% sizeInd = arrayfun(@(s) 1:s, size(DATA), 'UniformOutput', false);
% [Subject, Layer, BQ1, BQ2, TypeODC] = ndgrid(sizeInd{:});
% [Subject, Layer, TypeODC] = deal(categorical(Subject), categorical(Layer), categorical(TypeODC));
% T = table(DATA(:), Subject(:), Layer(:), BQ1(:), BQ2(:), TypeODC(:), 'VariableNames', {'rsFC', 'Subject', 'Layer', 'BQ1', 'BQ2', 'TypeODC'});
% %lme = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)'); % With all 2nd order interactions
% %lme = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC - BQ1:BQ2 - BQ1:TypeODC - BQ2:TypeODC - BQ1:Layer - BQ2:Layer - TypeODC:Layer - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)'); % With all 2nd order interactions
% lme = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC -TypeODC -Layer - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)'); % With all 2nd order interactions
% 
% lme
% 
% %% LME, mean over hemispheres and distances, 50% of beta matrix
% 
% % mean over distances, hemispheres for V1
% DATA_ALIKE_mean = squeeze(mean(mean(DATA_ALIKE, 4), 2));
% DATA_UNALIKE_mean = squeeze(mean(mean(DATA_UNALIKE, 4), 2));
% 
% DATA(:,:,:,:,1) = DATA_ALIKE_mean(:,:,[2,4,6,8,10],[2,4,6,8,10]);
% DATA(:,:,:,:,2) = DATA_UNALIKE_mean(:,:,[2,4,6,8,10],[2,4,6,8,10]);
% 
% sizeInd = arrayfun(@(s) 1:s, size(DATA), 'UniformOutput', false);
% [Subject, Layer, BQ1, BQ2, TypeODC] = ndgrid(sizeInd{:});
% [Subject, Layer, TypeODC] = deal(categorical(Subject), categorical(Layer), categorical(TypeODC));
% T = table(DATA(:), Subject(:), Layer(:), BQ1(:), BQ2(:), TypeODC(:), 'VariableNames', {'rsFC', 'Subject', 'Layer', 'BQ1', 'BQ2', 'TypeODC'});
% lme = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)'); % With all 2nd order interactions
% %lme = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC - BQ1:BQ2 - BQ1:TypeODC - BQ2:TypeODC - BQ1:Layer - BQ2:Layer - TypeODC:Layer - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)'); % With all 2nd order interactions
% 
% lme
% 
% %% LME
% 
% clear DATA
% 
% DATA_ALIKE_mean = squeeze(DATA_ALIKE);
% DATA_UNALIKE_mean = squeeze(DATA_UNALIKE);
% 
% DATA(:,:,:,:,:,:,1) = DATA_ALIKE_mean;
% DATA(:,:,:,:,:,:,2) = DATA_UNALIKE_mean;
% 
% sizeInd = arrayfun(@(s) 1:s, size(DATA), 'UniformOutput', false);
% [Subject, Dist, Layer, Hemi, BQ1, BQ2, TypeODC] = ndgrid(sizeInd{:});
% [Subject, Layer, Hemi, TypeODC] = deal(categorical(Subject), categorical(Layer), categorical(Hemi), categorical(TypeODC));
% T = table(DATA(:), Subject(:), Dist(:), Layer(:), Hemi(:), BQ1(:), BQ2(:), TypeODC(:), 'VariableNames', {'rsFC', 'Subject', 'Dist', 'Layer', 'Hemi', 'BQ1', 'BQ2', 'TypeODC'});
% %lme = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject:Hemi)'); % With all 2nd order interactions
% lme = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC - BQ1:BQ2 - BQ1:TypeODC - BQ2:TypeODC - BQ1:Layer - BQ2:Layer - TypeODC:Layer - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject:Hemi)'); % With all 2nd order interactions
% 
% lme
% 
% 
% %% LME, mean over distances
% clear DATA
% % mean over distances
% DATA_ALIKE_mean = squeeze(mean(DATA_ALIKE, 2));
% DATA_UNALIKE_mean = squeeze(mean(DATA_UNALIKE, 2));
% 
% DATA(:,:,:,:,:,1) = DATA_UNALIKE_mean;
% DATA(:,:,:,:,:,2) = DATA_ALIKE_mean;
% 
% 
% sizeInd = arrayfun(@(s) 1:s, size(DATA), 'UniformOutput', false);
% [Subject, Layer, Hemi, BQ1, BQ2, TypeODC] = ndgrid(sizeInd{:});
% [Subject, Layer, Hemi, TypeODC] = deal(categorical(Subject), categorical(Layer), categorical(Hemi), categorical(TypeODC));
% 
% T = table(DATA(:), Subject(:), Layer(:), Hemi(:), BQ1(:), BQ2(:), TypeODC(:), 'VariableNames', {'rsFC', 'Subject', 'Layer', 'Hemi', 'BQ1', 'BQ2', 'TypeODC'});
% %lme = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)'); % With all 2nd order interactions
% % lme_IA = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC -Layer -TypeODC - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)');
% % lme_SN = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC -BQ1 -BQ2 -Layer -TypeODC - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)');
% % lme_SN_2 = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC -BQ1 -BQ2 - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)');
% lme_SN_3 = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC -BQ1 -BQ2 - BQ1:Layer - BQ1:TypeODC - BQ2:Layer - BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC + (1|Subject)');
% % %lme_SN_3 = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC -BQ1 -BQ2 - BQ1:Layer - BQ1:TypeODC - BQ2:Layer - BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)');
% 
% %% LME, mean over distances
% clear DATA
% % mean over distances
% DATA_ALIKE_mean = squeeze(mean(DATA_ALIKE, 2));
% DATA_UNALIKE_mean = squeeze(mean(DATA_UNALIKE, 2));
% 
% DATA(:,:,:,:,:,1) = DATA_UNALIKE_mean;
% DATA(:,:,:,:,:,2) = DATA_ALIKE_mean;
% 
% 
% sizeInd = arrayfun(@(s) 1:s, size(DATA), 'UniformOutput', false);
% [Subject, Layer, Hemi, BQ1, BQ2, TypeODC] = ndgrid(sizeInd{:});
% [Subject, Layer, Hemi, TypeODC] = deal(categorical(Subject), categorical(Layer), categorical(Hemi), categorical(TypeODC));
% prodBQ = BQ1 .* BQ2;
% 
% T = table(DATA(:), Subject(:), Layer(:), Hemi(:), prodBQ(:), TypeODC(:), 'VariableNames', {'rsFC', 'Subject', 'Layer', 'Hemi', 'prodBQ', 'TypeODC'});
% %lme = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)'); % With all 2nd order interactions
% % lme_IA = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC -Layer -TypeODC - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)');
% % lme_SN = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC -BQ1 -BQ2 -Layer -TypeODC - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)');
% % lme_SN_2 = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC -BQ1 -BQ2 - BQ1:BQ2:Layer - BQ1:BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC - BQ1:BQ2:Layer:TypeODC + (1|Subject)');
% %lme_SN_3 = fitlme(T, 'rsFC ~ BQ1*BQ2*Layer*TypeODC -BQ1 -BQ2 - BQ1:Layer - BQ1:TypeODC - BQ2:Layer - BQ2:TypeODC - BQ1:Layer:TypeODC - BQ2:Layer:TypeODC + (1|Subject)');
% lme = fitlme(T, 'rsFC ~ prodBQ*Layer*TypeODC + (1|Subject)');
% 
% if saveTables
%     [xx, xxx, Coefficients] = fixedEffects(lme, 'Alpha', 0.05);
%     resultsTable = table(lme.CoefficientNames', Coefficients.Estimate, Coefficients.SE, Coefficients.tStat, Coefficients.DF, Coefficients.pValue, Coefficients.Upper, Coefficients.Lower, ...
%                      'VariableNames', {'Name', 'Estimate', 'SE', 'tStat', 'DF', 'pValue', 'Upper', 'Lower'});
%     writetable(resultsTable, [savePath 'V1_lme.xlsx']);
% end