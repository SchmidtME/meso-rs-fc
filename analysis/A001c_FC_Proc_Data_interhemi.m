function Data_Combined = A001c_FC_Proc_Data_interhemi(subName, Root, AnalysisParam, TrgFile)
%%%
tic
subName

Root = fullfile(Root, subName);
mkdir(Root)
hemis = {'lh', 'rh'};

% TimeSeries Params
StrtPnt = AnalysisParam.StrtPnt;
TmSeries_Length = AnalysisParam.TmSeries_Length;
timeSeriesRange = {StrtPnt:(StrtPnt+TmSeries_Length)};

AnalysisParam.TrgFile = TrgFile;

% additional preprocessing
detrending = AnalysisParam.detrending ;
hpf = AnalysisParam.hpf;

Trg_File = AnalysisParam.TrgFile;

% set directories for the anatomical data, ODC differential maps and resting-state data

anatFolder = fullfile('/autofs/space/ardebil_001/users/Shared/good_subjects_anat', [subName(1:4) '_anat_upsample_B1justsub']);
labelFolder = fullfile('/space/ardebil/1/users/Others/Marianna/good_subjects_anat/retinotopy', [subName(1:4)]);
odcFolder = fullfile('/autofs/space/ardebil_002/users/Shahin/Stereopsis_Upsampled/Subjects_Control_Upsampled/',  [subName(1:4) '_Final_Upsampled_VOL']);
folder = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Data_Denoised';

rsFolder =  fullfile(folder, [subName], 'bold_Close_Upsampled2');

%% Load and show anatomical data and ODC maps

% get folder name where ODC data is located
d = dir(fullfile(odcFolder, '*Smoothing_0-2.lh')).name(1:end-3); %Stereopsis_TR3_Columnar_Smoothing_0-2

% loop over hemispheres
for h = 1:numel(hemis)
    % specify hemisphere
    hemi = hemis{h};
    fprintf('Loading %s Data ... \n', hemi)
    StrtTime = toc;
    % load ODC map
    fprintf('Loading %s ODC Map \n', hemi)
    dataFile_odc = fullfile(odcFolder, [d '.' hemi], 'R_C', [hemi '.ffx.osgm.wls'], 'beta.nii');
    sig_odc{h} = load_nifti(dataFile_odc); 
    sig_odc{h} = sig_odc{h}.vol; 
    fprintf('Done in %s s!\n', num2str(toc- StrtTime))
    
    % load whole brain vertices and faces & adjust index of faces
    fprintf('Loading %s Surface and %s Label\n\n', hemi, AnalysisParam.ROI) 
    StrtTime = toc;
    
    [vtx{h}, fac{h}] = read_surf(fullfile(anatFolder, 'surf', [hemi '.inflated'])); 
    fac{h} = fac{h} + 1;

    % load flattened, 2D occipital patch & adjust index of vertices
    occ_patch{h} = read_patch(fullfile(anatFolder, 'surf', sprintf('%s.%s', hemi,AnalysisParam.ROIpatch))); 
    occ_patch{h}.ind = occ_patch{h}.ind+1;

    % load labels defined on upsampled retinotopy
    if strcmp(AnalysisParam.ROI, 'V1')
        roi_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Upsampled_Rtopy.label']));
    elseif strcmp(AnalysisParam.ROI, 'V2')
        roi_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V2_Upsampled_Rtopy.label']));
    elseif strcmp(AnalysisParam.ROI, 'V3')
        roi_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V3_Upsampled_Rtopy.label']));
    elseif strcmp(AnalysisParam.ROI, 'whole_V1')
        roi_label = read_ROIlabel(fullfile(anatFolder, 'label', [hemi '.V1_exvivo_Upsampled.label']));
    elseif strcmp(AnalysisParam.ROI, 'whole_V2')
        roi_label = read_ROIlabel(fullfile(anatFolder, 'label', [hemi '.V2_exvivo_Upsampled.label']));
    % elseif strcmp(AnalysisParam.ROI, 'V1_Posterior_Posterior')
    %     roi_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Posterior_Posterior_Upsampled_Rtopy.label']));
    % elseif strcmp(AnalysisParam.ROI, 'V1_Posterior')
    %     V1_label = read_ROIlabel(fullfile(anatFolder, 'label', 'High_Res_Upsampled', [hemi '.V1_Upsampled_Adjusted.label']));
    %     V1_Posterior_label = read_ROIlabel(fullfile(anatFolder, 'label', 'High_Res_Upsampled', [hemi '.Posterior_Upsampled.label']));
    %     [~, idx] = ismember(V1_label, V1_Posterior_label);
    %     roi_label = V1_label(idx ~= 0, :);
    elseif strcmp(AnalysisParam.ROI, 'V1_Posterior')
        roi_label_1 = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Posterior_Upsampled_Rtopy.label']));
        roi_label_2 = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Posterior_Posterior_Upsampled_Rtopy.label']));
        roi_label = vertcat(roi_label_1, roi_label_2);
    % elseif strcmp(AnalysisParam.ROI, 'V1_Anterior')
    %     V1_label = read_ROIlabel(fullfile(anatFolder, 'label', 'High_Res_Upsampled', [hemi '.V1_Upsampled_Adjusted.label']));
    %     V1_Dorsal_label = read_ROIlabel(fullfile(anatFolder, 'label', 'High_Res_Upsampled', [hemi '.Posterior_Upsampled.label']));
    %     unique_center_labels = unique(V1_Dorsal_label);
    %     rows_to_keep = ~ismember(V1_label, unique_center_labels);
    %     roi_label = V1_label(rows_to_keep, :);
    elseif strcmp(AnalysisParam.ROI, 'V1_Anterior')
        roi_label_1 = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Anterior_Upsampled_Rtopy.label']));
        roi_label_2 = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Anterior_Anterior_Upsampled_Rtopy.label']));
        roi_label = vertcat(roi_label_1, roi_label_2);
    elseif strcmp(AnalysisParam.ROI, 'V1_Far_Periphery')
        roi_label_1 = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Upsampled_Rtopy.label']));
        roi_label_2 = read_ROIlabel(fullfile(anatFolder, 'label', [hemi '.V1_exvivo_Upsampled.label']));
        unique_center_labels = unique(roi_label_1);
        rows_to_keep = ~ismember(roi_label_2, unique_center_labels);
        roi_label = roi_label_2(rows_to_keep, :);
    end

    % crop occipital patch to V1 patch
    v1_patch{h} = struct('npts','ind','x','y','z','vno');
    v1_patch{h}.ind = occ_patch{h}.ind(ismember(occ_patch{h}.ind, roi_label));
    v1_patch{h}.x = occ_patch{h}.x(ismember(occ_patch{h}.ind, roi_label));
    v1_patch{h}.y = occ_patch{h}.y(ismember(occ_patch{h}.ind, roi_label));
    v1_patch{h}.z = occ_patch{h}.z(ismember(occ_patch{h}.ind, roi_label));
    v1_patch{h}.vno = occ_patch{h}.vno(ismember(occ_patch{h}.ind, roi_label));
    v1_patch{h}.npts = length(v1_patch{h}.ind); 

    % create binary mask of vertices that are included in the roi
    bin_mask_vtx_occ_patch{h} = false(size(vtx{h},1),1); bin_mask_vtx_occ_patch{h}(occ_patch{h}.ind) = true;
    bin_mask_vtx_v1_patch{h} = false(size(vtx{h},1),1); 
    bin_mask_vtx_v1_patch{h}(v1_patch{h}.ind) = true;

end

%% Load and show resting-state functional data
% get rs session names
d = dir(fullfile(rsFolder, '0*')); 
rsSessions = {d.name};
if exist('sessionSubset', 'var') && ~isempty(sessionSubset)
    rsSessions = rsSessions(sessionSubset);
end

fprintf('%s Resting-State Runs were found! \n', num2str(length(rsSessions)))

%%

% load data for every session, do detrend & hpf, compute partialcorr, save it to matrix (sessions, partialcorr), then do averaging
for sessionNum = 1:length(rsSessions)
    fprintf('Loading Restig-Sate Data Run %s! \n', num2str(sessionNum))
    StrtTime = toc;
    clear rs_v1_patch rs_v1_patchBoth

    % initialize wm and mcpr regressors and rs within patch?
    wm = []; mcpr = []; rs_v1_patch{h} = [];

    for h = 1:numel(hemis)
        hemi = hemis{h};
        fprintf('Loading %s! \n', hemi)

        % check which session should be analyzed or if it is all sessions combined
        rsSessionName = rsSessions{sessionNum};
        sessionNumTempAll = sessionNum;
        
        for sessionNumTemp = sessionNumTempAll
            % load the rs data
            fprintf('The target file is : %s \n', [hemi Trg_File])

            dataFile_rs = fullfile(rsFolder, rsSessions{sessionNumTemp}, [hemi Trg_File]);
            [rs, mghM{h}] = load_mgh([dataFile_rs '.mgz']);
            rs = squeeze(rs);

            fprintf('The original resting state data martix is : %d x %d \n', size(rs))
            fprintf('We will use these data points: %d - %d \n', [StrtPnt StrtPnt+TmSeries_Length])
            
            rs = rs(:, [StrtPnt:StrtPnt+TmSeries_Length]);
            fprintf('The used resting state data martix is : %d x %d \n', size(rs))


            % loads volumes and 4x4 vox2ras transform
            if detrending
                rs = detrend(rs', 2)';
            end

            if hpf
                rs = highpass(rs', 0.01, 1/4)';
            end

            % only consider rs within v1 patch
            rs_v1_patch{h} = [rs_v1_patch{h} squeeze(rs(bin_mask_vtx_v1_patch{h},:))];
        end
    end

    % load the covariate data and append for each session
    wm0 = load(fullfile(rsFolder, rsSessions{sessionNumTemp}, 'wm.dat'), '-ascii');
    wm = [wm; wm0(timeSeriesRange{:})];
    mcpr0 = load(fullfile(rsFolder, rsSessions{sessionNumTemp}, 'mcprextreg'), '-ascii');
    mcpr = [mcpr; mcpr0(timeSeriesRange{:},:)];
    % Compute correlations based on the current quantile
    % reorganize the data so that bin_mask_vtx_eye2 vertices are appended below bin_mask_vtx_eye1 vertices
    corr_v1_patch = partialcorr(rs_v1_patch{1}', rs_v1_patch{2}', [mcpr(:,1:3) wm]);
    % Store the correlation matrix for the current quantile
    corr_v1_patch_allSess(sessionNum,:,:) = corr_v1_patch;
    %clear corr_v1_patch
    fprintf('Done in %s s!\n\n', num2str(toc- StrtTime))        
end

%% do analysis for mean across sessions
StrtTime = toc;

Data_Combined = [];

% correlation
corr_v1_patch_mean = squeeze(mean(corr_v1_patch_allSess(:,:,:), 1));
% mean for distance quantiles
mean_quants_corr_v1_patch = mean(corr_v1_patch_mean, 'all');


Data_Combined(1) = mean_quants_corr_v1_patch;
% z-transformation
mean_quants_corr_v1_patch_z = 0.5 * log((1 + mean_quants_corr_v1_patch) ./ (1 - mean_quants_corr_v1_patch));
Data_Combined(2) = mean_quants_corr_v1_patch_z;

fprintf('Done in %s s!\n\n', num2str(toc- StrtTime))

fprintf('Total time elapsed: %s s!\n\n', num2str(toc))

save([Root,'/CorrelationMtx_justFC.mat'], 'Data_Combined', 'AnalysisParam');
