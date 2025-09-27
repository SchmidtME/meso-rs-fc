clear all
clc
close all
%%


subs = {'aroo'};
hemis = {'rh'};
Layers = {'4-6'};
s = 1; h = 1;

outDir = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Talks/Magdeburg/Figures/'; % Change this to your desired folder

%% get beta thresholds

load('/autofs/space/ardebil_001/users/Others/Marianna/FC_7T_Coronal/Controls/Results/V1_layers_0-10_intrahemispheric/aroo/CorrelationMtx_Selectivity_Params_subsampled_1.mat');

NDE = mean(([AnalysisParam.betas_eye1_medians{1};AnalysisParam.betas_eye1_medians{2}]),1)

DE = mean(([AnalysisParam.betas_eye2_medians{1};AnalysisParam.betas_eye2_medians{2}]),1)

both = [NDE([1,3,5,7,9]) 0 DE([1,3,5,7,9])];

%% 
    
anatFolder = fullfile('/autofs/space/ardebil_001/users/Shared/good_subjects_anat', [subs{s} '_anat_upsample_B1justsub']);
labelFolder = fullfile('/space/ardebil/1/users/Others/Marianna/good_subjects_anat/retinotopy', subs{s});
odcFolder = fullfile('/autofs/space/ardebil_002/users/Shahin/Stereopsis_Upsampled/Subjects_Control_Upsampled/',  [subs{s} '_Final_Upsampled_VOL']);

% load V1 label
[vtx{h}, fac{h}] = read_surf(fullfile(anatFolder, 'surf', [hemis{h} '.inflated'])); 
fac{h} = fac{h} + 1;
occ_patch{h} = read_patch(fullfile(anatFolder, 'surf', sprintf('%s.%s', hemis{h}, 'V1_patch.flat'))); 
occ_patch{h}.ind = occ_patch{h}.ind+1;
bin_mask_vtx_occ_patch{h} = false(size(vtx{h},1),1); 
bin_mask_vtx_occ_patch{h}(occ_patch{h}.ind) = true;    
roi_label = read_ROIlabel(fullfile(labelFolder, [hemis{h} '.V1_Upsampled_Rtopy_4Talk.label']));
v1_patch{h} = struct('npts','ind','x','y','z','vno');
v1_patch{h}.ind = occ_patch{h}.ind(ismember(occ_patch{h}.ind, roi_label));
v1_patch{h}.x = occ_patch{h}.x(ismember(occ_patch{h}.ind, roi_label));
v1_patch{h}.y = occ_patch{h}.y(ismember(occ_patch{h}.ind, roi_label));
v1_patch{h}.z = occ_patch{h}.z(ismember(occ_patch{h}.ind, roi_label));
v1_patch{h}.vno = occ_patch{h}.vno(ismember(occ_patch{h}.ind, roi_label));
v1_patch{h}.npts = length(v1_patch{h}.ind); 
bin_mask_vtx_v1_patch{h} = false(size(vtx{h},1),1); 
bin_mask_vtx_v1_patch{h}(v1_patch{h}.ind) = true;

% just for visualization
% only include fac and vtx within v1 patch
nonbin_mask_vtx_v1_patch = double(bin_mask_vtx_v1_patch{h}); 
nonbin_mask_vtx_v1_patch(bin_mask_vtx_v1_patch{h}) = 1:sum(bin_mask_vtx_v1_patch{h});
% only include fac and vtx within v1 occ
nonbin_mask_vtx_occ_patch = double(bin_mask_vtx_occ_patch{h}); 
nonbin_mask_vtx_occ_patch(bin_mask_vtx_occ_patch{h}) = 1:sum(bin_mask_vtx_occ_patch{h});
fac_occ_patch{h} = nonbin_mask_vtx_occ_patch(fac{h}); 
fac_occ_patch{h} = fac_occ_patch{h}(all(fac_occ_patch{h},2),:);
vtx_occ_patch{h} = nan(size(vtx{h})); 
vtx_occ_patch{h} = [occ_patch{h}.x' occ_patch{h}.y' occ_patch{h}.z'];
fac_v1_patch{h} = nonbin_mask_vtx_v1_patch(fac{h}); 
fac_v1_patch{h} = fac_v1_patch{h}(all(fac_v1_patch{h},2),:);
vtx_v1_patch{h} = nan(size(vtx{h})); 
vtx_v1_patch{h} = [v1_patch{h}.x' v1_patch{h}.y' v1_patch{h}.z'];
% faces mask
fLabelMask = zeros(size(bin_mask_vtx_v1_patch{h})); 
fLabelMask(bin_mask_vtx_v1_patch{h}) = 1:sum(bin_mask_vtx_v1_patch{h});
facesMask{h} = fLabelMask(fac{h}); 
facesMask{h} = facesMask{h}(all(facesMask{h},2),:);

% load odc map
d = dir(fullfile(odcFolder, '*Smoothing_0-2.lh')).name(1:end-3); 
sig_odc{h} = load_nifti(fullfile(odcFolder, [d '.' hemis{h}], 'R_C', [hemis{h} '.ffx.osgm.wls'], 'beta.nii')).vol; 
sigMasked{h} = sig_odc{h}(v1_patch{h}.ind); 

%% load rs map

for l = 1:length(Layers)
    for s = 1:length(subs)
        
        subName = subs{s};
        
        for h = 1:length(hemis)
        
            corrLabel_sessions = [];
            hemi = hemis{h};
            disp(['Hemisphere: ' hemi])
            disp('Loading the data...')
            tic
        
            rsFolder = fullfile('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Data_Denoised', subName, 'bold_Close_Upsampled2');
        
            d = dir(fullfile(rsFolder, '0*')); 
            rsSessions = {d.name};
        
            toc
        
            disp('Resting-state data')
        
           for sessionNum = 1:length(rsSessions)
                
                disp(sprintf('Process run %d', sessionNum))
                tic
        
                wm = []; mcpr = []; rs_v1_patch{h} = [];
        
                % load rs-FC
            
                [rs, mghM{h}] = load_mgh([fullfile(rsFolder, rsSessions{sessionNum}, [hemi sprintf(['.fmcpr.sm0.self.midgray.00.nb1_rad','%s'], Layers{l})]) '.mgz']);
                rs = squeeze(rs);
                rs = rs(:, [10:128]);
                rs = detrend(rs', 2)';
                rs = highpass(rs', 0.01, 1/4)';
                rs_v1_patch{h} = [rs_v1_patch{h} squeeze(rs(bin_mask_vtx_v1_patch{h},:))];
            
                % load covariates
                wm0 = load(fullfile(rsFolder, rsSessions{sessionNum}, 'wm.dat'), '-ascii');
                wm = [wm; wm0([10:128])];
                mcpr0 = load(fullfile(rsFolder, rsSessions{sessionNum}, 'mcprextreg'), '-ascii');
                mcpr = [mcpr; mcpr0([10:128],:)];
                
                corrLabel_sessions(sessionNum, :, :) = partialcorr(rs_v1_patch{h}', [mcpr(:,1:3) wm]);
        
                toc
            
           end

           corrLabel{h} = squeeze(mean(corrLabel_sessions,1));

        end
    end
end

%%

for i = 1:length(both)-1
    binnedData(sigMasked{h} >= both(i) & sigMasked{h} < both(i+1)) = i + 1;
end

% Final bin: values >= last edge
binnedData(sigMasked{h} >= both(end)) = length(both) + 1;

%%
clims = [-2,2];

cmap = linspace(0, 1, 13)';       % Column vector
cmap = repmat(cmap, 1, 3);        % Make it RGB: [gray gray gray]

% create colormap
% r = [0 0.5 1];
% g = [0 0.5 0];
% b = [1 0.5 0];
% x = linspace(-1, 1, 256)';
% redbluecmap = [interp1([-1 0 1], r, x), interp1([-1 0 1], g, x), interp1([-1 0 1], b, x)];

% Number of colors for each range
numColors = 10;

% Define RGB values for the color ranges
negStart = [0, 1, 1];   % Cyan (high negative)
negEnd = [0, 0, 0.5];   % Dark blue (low negative)

posStart = [1, 0, 0];   % Red (low positive)
posEnd = [1, 1, 0];     % Yellow (high positive)

% Generate the colormap for negative and positive values
negColors = [linspace(negStart(1), negEnd(1), numColors)' ...
             linspace(negStart(2), negEnd(2), numColors)' ...
             linspace(negStart(3), negEnd(3), numColors)'];

posColors = [linspace(posStart(1), posEnd(1), numColors)' ...
             linspace(posStart(2), posEnd(2), numColors)' ...
             linspace(posStart(3), posEnd(3), numColors)'];

% Combine negative and positive colormaps (no middle color for zero)
colormapData = [negColors; posColors];
    
%%
% Create an invisible figure
fig = figure('Visible', 'off', 'Color', 'w');

% Tiled layout for clean spacing
t = tiledlayout(fig, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

% Plot
nexttile;
trisurf(facesMask{h}, vtx_v1_patch{h}(:,1), vtx_v1_patch{h}(:,2), vtx_v1_patch{h}(:,3), sigMasked{h}, ...
    'EdgeColor', [0 0 0], 'LineWidth', 0.5);
axis equal off, shading interp, view(2)

clim([-1000, 1000])  % Fixed color limits
set(gca, 'Color', 'w');
colormap(colormapData);  % Your custom colormap

% Save the figure
filename = fullfile(outDir, 'ODC_bin_map.png');  % Change filename as needed

exportgraphics(fig, filename, 'Resolution', 600);  % Save at 300 DPI

% Close the figure to prevent display
close(fig);

%%

% Create an invisible figure
fig = figure('Visible', 'off', 'Color', 'w');

% Tiled layout
t = tiledlayout(fig, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

% Plot
nexttile;
trisurf(facesMask{h}, vtx_v1_patch{h}(:,1), vtx_v1_patch{h}(:,2), vtx_v1_patch{h}(:,3), sigMasked{h}*0, ...
    'EdgeColor', [0 0 0], 'LineWidth', 0.5);
axis equal off, shading interp, view(2)

clim(clims)
set(gca, 'Color', 'w');
colormap(gray);

% Save the figure
filename = fullfile(outDir, 'Distance_map.png'); % Change name as needed

exportgraphics(fig, filename, 'Resolution', 600); % 300 DPI export

% Close the figure to keep it hidden
close(fig);



%%

% Create an invisible figure
fig = figure('Visible', 'off', 'Color', 'w');

% Plotting
t = tiledlayout(fig, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
trisurf(facesMask{h}, vtx_v1_patch{h}(:,1), vtx_v1_patch{h}(:,2), vtx_v1_patch{h}(:,3), sigMasked{h}, ...
    'EdgeColor', [0 0 0], 'LineWidth', 0.5);
axis equal off, shading interp, view(2)
clim(clims)
set(gca, 'Color', 'w');
colormap(colormapData);

% Export the figure
filename = fullfile(outDir, 'ODI_map.png');

% Save with specific resolution (e.g., 300 DPI)
exportgraphics(fig, filename, 'Resolution', 600);

% Close the figure to avoid showing it
close(fig);


%%

for vol = [1:9 11:61]
    % Create an invisible figure
    fig = figure('Visible', 'off', 'Color', 'w');

    % Tiled layout for better spacing
    t = tiledlayout(fig, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    % Plot
    nexttile;
    trisurf(facesMask{h}, ...
            vtx_v1_patch{h}(:,1), ...
            vtx_v1_patch{h}(:,2), ...
            vtx_v1_patch{h}(:,3), ...
            rs_v1_patch{h}(:,vol), ...
            'EdgeColor', [0 0 0], 'LineWidth', 0.5);
    axis equal off;
    shading interp;
    view(2)

    clim([-12, 12]);
    set(gca, 'Color', 'w');
    colormap(redbluecmap);

    % Create filename
    filename = fullfile(outDir, sprintf('RsFC_vol_%02d.png', vol));  % Zero-padded numbering

    % Save figure
    exportgraphics(fig, filename, 'Resolution', 300);

    % Close the figure
    close(fig);
end

%%


%% Archive

% %%
% 
% figure(4);
% set(gcf, 'Color', 'w'); % Set background to white
% 
% t = tiledlayout(1,1, 'Padding', 'compact', 'TileSpacing', 'compact'); % Ensures better spacing
% 
% % First subplot
% nexttile;
% trisurf(facesMask{h}, vtx_v1_patch{h}(:,1), vtx_v1_patch{h}(:,2), vtx_v1_patch{h}(:,3), mean(rs_v1_patch{h},2), ...
%     'EdgeColor', [0 0 0], 'LineWidth', 0.5);axis equal off, shading interp, view(2)
% 
% clim([-2,2])
% set(gca, 'Color', 'w');
% colormap(redbluecmap);
% xlim([-20 -4]);
% ylim([-8 8]);
% 
% %%
% 
% seed_vertex = 131869;
% idx = find(v1_patch{h}.ind == seed_vertex );
% 
% %%
% 
% figure(4);
% set(gcf, 'Color', 'w'); % Set background to white
% 
% t = tiledlayout(1,1, 'Padding', 'compact', 'TileSpacing', 'compact'); % Ensures better spacing
% 
% % First subplot
% nexttile;
% trisurf(facesMask{h}, vtx_v1_patch{h}(:,1), vtx_v1_patch{h}(:,2), vtx_v1_patch{h}(:,3), corrLabel{h}(:,idx), ...
%     'EdgeColor', [0 0 0], 'LineWidth', 0.5);axis equal off, shading interp,view(2)
% 
% clim([-0.5,0.5])
% set(gca, 'Color', 'w');
% colormap(colormapData);
% 
% %%
% 
% figure(5);
% set(gcf, 'Color', 'w'); % Set background to white
% 
% t = tiledlayout(1,1, 'Padding', 'compact', 'TileSpacing', 'compact'); % Ensures better spacing
% 
% % First subplot
% nexttile;
% trisurf(facesMask{h}, vtx_v1_patch{h}(:,1), vtx_v1_patch{h}(:,2), vtx_v1_patch{h}(:,3), corrLabel{h}(:,idx), ...
%     'EdgeColor', [0 0 0], 'LineWidth', 0.5);axis equal off, shading interp, view(2)
% 
% clim([-0.5,0.5])
% set(gca, 'Color', 'w');
% colormap(colormapData);
% xlim([-20 -4]);
% ylim([-8 8]);
% 
% %%
% 
% figure(6);
% set(gcf, 'Color', 'w'); % Set background to white
% 
% t = tiledlayout(1,1, 'Padding', 'compact', 'TileSpacing', 'compact'); % Ensures better spacing
% 
% % First subplot
% nexttile;
% trisurf(facesMask{h}, vtx_v1_patch{h}(:,1), vtx_v1_patch{h}(:,2), vtx_v1_patch{h}(:,3), sign(sigMasked{h}), ...
%     'EdgeColor', [0 0 0], 'LineWidth', 0.5);axis equal off, shading interp, view(2)
% 
% clim([-1,1])
% set(gca, 'Color', 'w');
% colormap(gray);
% 
% xlim([-20 -4]);
% ylim([-8 8]);
% 
% %%
% 
% figure(7);
% set(gcf, 'Color', 'w'); % Set background to white
% 
% t = tiledlayout(1,1, 'Padding', 'compact', 'TileSpacing', 'compact'); % Ensures better spacing
% 
% % First subplot
% nexttile;
% trisurf(facesMask{h}, vtx_v1_patch{h}(:,1), vtx_v1_patch{h}(:,2), vtx_v1_patch{h}(:,3), sigMasked{h}, ...
%     'EdgeColor', [0 0 0], 'LineWidth', 0.5);axis equal off, shading interp, view(2)
% 
% %clim(clims)
% clim([-1000,1000])
% set(gca, 'Color', 'w');
% colormap(colormapData);
% 
% xlim([-20 -4]);
% ylim([-8 8]);
% 
% %%
% 
% figure(8);
% set(gcf, 'Color', 'w'); % Set background to white
% 
% t = tiledlayout(1,1, 'Padding', 'compact', 'TileSpacing', 'compact'); % Ensures better spacing
% 
% % First subplot
% nexttile;
% trisurf(facesMask{h}, vtx_v1_patch{h}(:,1), vtx_v1_patch{h}(:,2), vtx_v1_patch{h}(:,3), 0*sigMasked{h}, ...
%     'EdgeColor', [0 0 0], 'LineWidth', 0.5);axis equal off, shading interp, view(2)
% 
% clim(clims)
% set(gca, 'Color', 'w');
% colormap(gray);
% 
% xlim([-20 -4]);
% ylim([-8 8]);