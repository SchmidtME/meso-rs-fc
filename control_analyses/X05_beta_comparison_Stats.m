close all
clear all
clc

%%

saveFigures = 1;
if saveFigures
    currentDate = datestr(now, 'yyyy-mm-dd');
    saveDir = sprintf('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Supplementary/Figure_4/%s', currentDate);
    mkdir(saveDir);
end

saveDir_data = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figure_S4_right_plot_data';
mkdir(saveDir_data);

Dest_Folder = ['/autofs/space/ardebil_002/users/Shahin/Stereopsis_Upsampled/Codes/Paper/Stereo_Deficiency/Graphs/', date]; %mkdir(Dest_Folder)
ErrobarPlots = 0; MK_SIZE = 10; LN_SIZE = 4;

addpath('/autofs/space/ardebil_002/users/Shahin/Stereopsis_Upsampled/Codes')
Sbj = {'ylri', 'aroo', 'auil', 'aman', 'arak', 'imyy', 'evad', 'chss', 'haas', 'atib', 'rcgr', ...
       'myla', };                            % SD


Sbj_Grp = {'Control', 'Control', 'Control', 'Control', 'Control', 'Control', 'Control', 'Control', 'Control', 'Control', 'Control',...
    'Control'};      % anisometropic

Excel_Add = '/autofs/space/ardebil_002/users/Shahin/Stereopsis_Upsampled/Demog/Patient_summary_fmri-amblyopia study.xlsx'
Beh_Data = readtable(Excel_Add);

Dom_Eye = {'R', 'L', 'R', 'L', 'R', 'R', 'R', 'R', 'R', 'R', 'R'...
           'R'}                                                   % anisometropic

Type = [1 1 1 1 1 1 1 1 1 1 1 ...calendar
        2 ]        % anisometropic

Anat_Address = '/autofs/space/ardebil_001/users/Shared/good_subjects_anat/'
Func_Address = '/autofs/space/ardebil_002/users/Shahin/Stereopsis_Upsampled/'
Hemis = {'rh', 'lh'}
Labels = {'V1'}
APPLY_ABS =1;

%%

for Lab = 1:length(Labels)
    V1_Data = []

    for Hemi = 1:2

        for i=1:length(Sbj)

            ADD = [Anat_Address, Sbj{i}, '_anat_upsample_B1justsub/label/High_Res_Upsampled/', Hemis{Hemi}, '.', Labels{Lab}, '_Upsampled_Adjusted2.label']
            if ~exist(ADD)
                ADD = [Anat_Address, Sbj{i}, '_anat_upsample_B1justsub/label/High_Res_Upsampled/', Hemis{Hemi}, '.', Labels{Lab}, '_Upsampled.label']
            end

            [num, mtx] = mris_read_label_Full(ADD);
            Vertices = mtx(:, 1) + 1;
    
            if strcmp(Sbj_Grp{i} , 'Control')
                ADD = [Func_Address, 'Subjects_', Sbj_Grp{i}, '_Upsampled/']
            else
                ADD = [Func_Address, 'Subjects_', Sbj_Grp{i}, '_Upsampled/']
            end
            
            for ii=1:10
                ADD2 = [ADD, Sbj{i}, num2str(ii), '/bold_Upsampled2/Stereopsis_TR3_Columnar_Smoothing_0-2.', Hemis{Hemi}, '/2D_R/cespct.nii.gz']
                if ~exist(ADD2)
                    break
                end
                H1 = load_nifti(ADD2);

                ADD2 = [ADD, Sbj{i}, num2str(ii), '/bold_Upsampled2/Stereopsis_TR3_Columnar_Smoothing_0-2.', Hemis{Hemi}, '/2D_C/cespct.nii.gz']
                if ~exist(ADD2)
                    break
                end
                H2 = load_nifti(ADD2);

                if strcmp(Dom_Eye{i}, 'L')
                    TMP(ii, 1:2) = [nanmean((H1.vol(Vertices))) nanmean((H2.vol(Vertices)))];
                else
                    TMP(ii, 1:2) = [nanmean((H2.vol(Vertices))) nanmean((H1.vol(Vertices)))];
                end

            end
    
            if ii==1
                ADD_temp = '/autofs/space/ardebil_002/users/Shahin/ODC/Subjects_Upsampled/';
                Sess = dir([ADD_temp, Sbj{i}, '_red*'])
                for ii=1:length(Sess)
                    ADD2 = [ADD_temp, Sess(ii).name, '/bold_Upsampled2/ODC_Columnar_0-2.', Hemis{Hemi}, '/C1/cespct.nii.gz']                    
                    H1 = load_nifti(ADD2);

                    ADD2 = [ADD_temp, Sess(ii).name, '/bold_Upsampled2/ODC_Columnar_0-2.', Hemis{Hemi}, '/C2/cespct.nii.gz']
                    H2 = load_nifti(ADD2);

                    if strcmp(Dom_Eye{i}, 'L')
                        TMP(ii, 1:2) = [nanmean((H1.vol(Vertices))) nanmean((H2.vol(Vertices)))];
                    else
                        TMP(ii, 1:2) = [nanmean((H2.vol(Vertices))) nanmean((H1.vol(Vertices)))];
                    end
                end
            end
            
            tmp = nanmean(TMP)

            V1_Data(i, 1:2, Hemi) = [tmp(2) tmp(1)];

            % V1_Data
            % % pause

        end
    end

    V1_Data = mean(V1_Data, 3)
end
%%

[h, p, a, stat] = ttest(V1_Data(:, 1), V1_Data(:, 2))

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
scatter(x_DE * ones(12, 1), abs(V1_Data(:,1)), 100, 'k', 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'LineWidth', 3); % Black outline for DE
scatter(x_NDE * ones(12, 1), abs(V1_Data(:,2)), 100, [0.5 0.5 0.5], '^', 'filled', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 3); % Grey outline for NDE

% Plot lines connecting corresponding data points
for i = 1:length(V1_Data(:,1))
    plot([x_DE x_NDE], [abs(V1_Data(i,1)) abs(V1_Data(i,2))], 'k-', 'LineWidth', 1);
end

% Adjust the axes and labels
xlim([0.8 2]); % Keep space to the left
ylim([2 5.5]);

% Set y-ticks to just the min and max values
set(gca, 'YTick', [2 5.5]);

% Corrected x-ticks for the new positions
set(gca, 'XTick', [x_DE x_NDE], 'XTickLabel', {'DE', 'NDE'});
set(gca, 'LineWidth', 2); % Increase axis line width
set(gcf,'renderer', 'painters'); 

xlabel('Eye dominance');
ylabel('Mean percent signal change');
%title('Comparison of vertex counts for DE and NDE');

hold off;

if saveFigures
    print([saveDir '/Figure_4'], '-dtiff', '-r400');
end

%%

writematrix(V1_Data, fullfile(saveDir_data, 'V1_Data.csv'));
disp(['V1_Data saved to ' saveDir_data]);