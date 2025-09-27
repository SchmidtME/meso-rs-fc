function Data_Combined = A001b_FC_Proc_Data_subsample_beta(subName, Root, AnalysisParam, TrgFile, s)

warning('off', 'all');

% This script corresponds to Analysis B - The effect of ocular preference strength (beta),
% cortical depth, ROI (V1 subregions) (and type) on rs-FC (and selectivity).
% This function loads anatomical data, ODC (ocular dominance column) maps, 
% and resting-state fMRI data for a given subject. It then processes this data 
% by applying detrending and high-pass filtering (if specified), computes 
% partial correlations within the specified regions of interest (ROI), and 
% calculates the distance matrix between vertices in V1. The function groups 
% distances and beta values into quantiles and computes the correlation between resting-state 
% time series for different quantiles of vertex  and beta. The results are saved 
% for further analysis.

tic
subName

Root = fullfile(Root, subName);
mkdir(Root)
hemis = {'lh', 'rh'};

% indices for subsampling
selected_indices_eye1eye2 = load(fullfile('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results', [AnalysisParam.ROI '_layers_0-10_intrahemispheric'], subName, sprintf('CorrelationMtx_FC_Params_subsampled_%d.mat',s))).AnalysisParam.selected_indices_eye1eye2;
selected_indices_same_eye = load(fullfile('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Results', [AnalysisParam.ROI '_layers_0-10_intrahemispheric'], subName, sprintf('CorrelationMtx_FC_Params_subsampled_%d.mat',s))).AnalysisParam.selected_indices_same_eye;
% TimeSeries Params
StrtPnt = AnalysisParam.StrtPnt;
TmSeries_Length = AnalysisParam.TmSeries_Length;
timeSeriesRange = {StrtPnt:(StrtPnt+TmSeries_Length)};

AnalysisParam.TrgFile = TrgFile;
% number of quantiles for grouping of distances from selected vertex
nQuant = AnalysisParam.nQuant_beta;
quantile_beta_thresholds = 1:nQuant;
% minimal distance from selected vertex in mm
distThresh =  AnalysisParam.distThresh ;

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
    elseif strcmp(AnalysisParam.ROI, 'V1_Posterior')
        V1_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Upsampled_Rtopy.label']));
        V1_Posterior_label = read_ROIlabel(fullfile(labelFolder, [hemi '.Posterior_Upsampled_Rtopy.label']));
        [~, idx] = ismember(V1_label, V1_Posterior_label);
        roi_label = V1_label(idx ~= 0, :);
    elseif strcmp(AnalysisParam.ROI, 'V1_Anterior')
        V1_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Upsampled_Rtopy.label']));
        V1_Dorsal_label = read_ROIlabel(fullfile(labelFolder, [hemi '.Posterior_Upsampled_Rtopy.label']));
        unique_center_labels = unique(V1_Dorsal_label);
        rows_to_keep = ~ismember(V1_label, unique_center_labels);
        roi_label = V1_label(rows_to_keep, :);
    elseif strcmp(AnalysisParam.ROI, 'V1_Dorsal')
        V1_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Upsampled_Rtopy.label']));
        V1_Dorsal_label = read_ROIlabel(fullfile(labelFolder, [hemi '.Dorsal_Upsampled_Rtopy.label']));
        [~, idx] = ismember(V1_label, V1_Dorsal_label);
        roi_label = V1_label(idx ~= 0, :);
    elseif strcmp(AnalysisParam.ROI, 'V1_Ventral')
        V1_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Upsampled_Rtopy.label']));
        V1_Dorsal_label = read_ROIlabel(fullfile(labelFolder, [hemi '.Dorsal_Upsampled_Rtopy.label']));
        unique_center_labels = unique(V1_Dorsal_label);
        rows_to_keep = ~ismember(V1_label, unique_center_labels);
        roi_label = V1_label(rows_to_keep, :);
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

    % create binary masks for the two eyes just within V1 based on quantiles    
    % Extract masked values
    masked_values = sig_odc{h}(bin_mask_vtx_v1_patch{h});   
    % Separate negative and positive values
    negative_values = masked_values(masked_values < 0);
    positive_values = masked_values(masked_values > 0);

    negative_indices = find(masked_values < 0);
    positive_indices = find(masked_values > 0);
    % Compute quantiles for negative and positive values
    neg_quantiles = fliplr([quantile(negative_values, nQuant-1) 0]); % selects as many Quants as for distance
    pos_quantiles = fliplr([0 quantile(positive_values, nQuant-1)]); 
    % Create masks based on quantile threshold
    for quantile_beta_threshold = 1:10
        if quantile_beta_threshold == 10
            selected_neg_mask = negative_values <= neg_quantiles(quantile_beta_threshold);
            selected_pos_mask = positive_values >= pos_quantiles(end - quantile_beta_threshold + 1);
        else
            selected_neg_mask = (neg_quantiles(quantile_beta_threshold+1) <= negative_values) & (negative_values <= neg_quantiles(quantile_beta_threshold));
            selected_pos_mask = (positive_values >= pos_quantiles(end - quantile_beta_threshold + 1)) & (positive_values <= pos_quantiles(end - quantile_beta_threshold));
        end
        
        % Initialize the final masks for each quantile and save them in 3D arrays
        bin_mask_vtx_eye1{h}(:, quantile_beta_threshold) = false(size(roi_label,1), 1);
        bin_mask_vtx_eye2{h}(:, quantile_beta_threshold) = false(size(roi_label,1), 1);
        
        % For the right eye: Select vertices in the negative quantiles
        bin_mask_vtx_eye1{h}(negative_indices(selected_neg_mask), quantile_beta_threshold) = true;
    
        % For the left eye: Select vertices in the positive quantiles
        bin_mask_vtx_eye2{h}(positive_indices(selected_pos_mask), quantile_beta_threshold) = true;
        
        % remove unnecessary empty cells from the mask
        bin_mask_vtx_eye1{h} = bin_mask_vtx_eye1{h}(1:size(roi_label,1),:)
        bin_mask_vtx_eye2{h} = bin_mask_vtx_eye2{h}(1:size(roi_label,1),:)
        
    end
    
    % save number of included vertices for different eyes
    AnalysisParam.num_included_vertices_eye1{h} = sum(bin_mask_vtx_eye1{h},'all');
    AnalysisParam.num_included_vertices_eye2{h} = sum(bin_mask_vtx_eye2{h},'all');

    for q1 = 1:numel(quantile_beta_thresholds) 
        for q2 = 1:numel(quantile_beta_thresholds)
            % sometimes quantiles include one vertex more than the other, so one random vertex is excluded
            if sum(bin_mask_vtx_eye1{h}(:,q1)) > sum(bin_mask_vtx_eye1{h}(:,q2))
                column_data = bin_mask_vtx_eye1{h}(:, q1);
                indices = find(column_data == 1);
                %random_index = indices(randi(numel(indices)));
                random_index = indices(round(numel(indices)/2));
                column_data(random_index) = 0;
                bin_mask_vtx_eye1{h}(:, q1) = column_data;
            elseif sum(bin_mask_vtx_eye1{h}(:,q1)) < sum(bin_mask_vtx_eye1{h}(:,q2))
                column_data = bin_mask_vtx_eye1{h}(:, q2);
                indices = find(column_data == 1);
                %random_index = indices(randi(numel(indices)));
                random_index = indices(round(numel(indices)/2));
                column_data(random_index) = 0;
                bin_mask_vtx_eye1{h}(:, q2) = column_data;
            end
            if sum(bin_mask_vtx_eye2{h}(:,q1)) > sum(bin_mask_vtx_eye2{h}(:,q2))
                column_data = bin_mask_vtx_eye2{h}(:, q1);
                indices = find(column_data == 1);
                %random_index = indices(randi(numel(indices)));
                random_index = indices(round(numel(indices)/2));
                column_data(random_index) = 0;
                bin_mask_vtx_eye2{h}(:, q1) = column_data;
            elseif sum(bin_mask_vtx_eye2{h}(:,q1)) < sum(bin_mask_vtx_eye2{h}(:,q2))
                column_data = bin_mask_vtx_eye2{h}(:, q2);
                indices = find(column_data == 1);
                %random_index = indices(randi(numel(indices)));
                random_index = indices(round(numel(indices)/2));
                column_data(random_index) = 0;
                bin_mask_vtx_eye2{h}(:, q2) = column_data;
            end 
        end
    end

    sig_odc_v1{h} = sig_odc{h}(bin_mask_vtx_v1_patch{h});

    for q = 1:10

        betas_eye1{h}(q) = median(sig_odc_v1{h}(bin_mask_vtx_eye1{h}(:,q)));
        betas_eye2{h}(q) = median(sig_odc_v1{h}(bin_mask_vtx_eye2{h}(:,q)));
 
    end

    AnalysisParam.betas_eye1_medians{h} = betas_eye1{h};
    AnalysisParam.betas_eye2_medians{h} = betas_eye2{h};

    % only include fac and vtx within v1 patch
    nonbin_mask_vtx_v1_patch = double(bin_mask_vtx_v1_patch{h}); 
    nonbin_mask_vtx_v1_patch(bin_mask_vtx_v1_patch{h}) = 1:sum(bin_mask_vtx_v1_patch{h});

    % crop the data
    % only include fac and vtx within v1 occ
    nonbin_mask_vtx_occ_patch = double(bin_mask_vtx_occ_patch{h}); nonbin_mask_vtx_occ_patch(bin_mask_vtx_occ_patch{h}) = 1:sum(bin_mask_vtx_occ_patch{h});
    fac_occ_patch{h} = nonbin_mask_vtx_occ_patch(fac{h}); fac_occ_patch{h} = fac_occ_patch{h}(all(fac_occ_patch{h},2),:);
    vtx_occ_patch{h} = nan(size(vtx{h})); vtx_occ_patch{h} = [occ_patch{h}.x' occ_patch{h}.y' occ_patch{h}.z'];
    
    fac_v1_patch{h} = nonbin_mask_vtx_v1_patch(fac{h}); 
    fac_v1_patch{h} = fac_v1_patch{h}(all(fac_v1_patch{h},2),:);

    vtx_v1_patch{h} = nan(size(vtx{h})); 
    vtx_v1_patch{h} = [v1_patch{h}.x' v1_patch{h}.y' v1_patch{h}.z'];
    
    fprintf('Done in %s s!\n\n', num2str(toc- StrtTime));        
end

%% Create distance map for vertices

fprintf('Creating Distance Matrix! \n')
StrtTime = toc;

for h = 1:numel(hemis)

    % Calculate distance matrix of Euclidean distances between the vertices included in the patch
    dist_v1_patch{h} = sqrt((v1_patch{h}.x - v1_patch{h}.x').^2 + (v1_patch{h}.y - v1_patch{h}.y').^2 + (v1_patch{h}.z - v1_patch{h}.z').^2);    

end

fprintf('Done in %s s!\n\n', num2str(toc- StrtTime))

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

    % Correlation between different eyes, same hemispheres
    %% Initialize Quantile-Based Correlation Storage
    for h = 1:numel(hemis)
        hemi = hemis{h};
        % Compute correlations based on the current quantile
        corr_v1_patch = partialcorr(rs_v1_patch{h}', [mcpr(:,1:3) wm]);
        % Store the correlation matrix for the current quantile
        corr_v1_patch_allSess{h}(sessionNum,:,:) = corr_v1_patch;
        clear corr_v1_patch
    end
    fprintf('Done in %s s!\n\n', num2str(toc- StrtTime))        
end

%% do analysis for mean across sessions
fprintf('Saving the data! \n')
StrtTime = toc;

Data_Combined = zeros(nQuant, 4, 2, nQuant, nQuant,2);
for h = 1:numel(hemis)
     for q1 = 1:numel(quantile_beta_thresholds)  % Assuming quantile_beta_threshold contains thresholds
        for q2 = 1:numel(quantile_beta_thresholds)
            hemi = hemis{h};
        
            corr_v1_patch_mean = squeeze(mean(corr_v1_patch_allSess{h}(:,:,:), 1));
            
            % Find the indices where bin_mask_vtx_eye1 and bin_mask_vtx_eye2 are 1
            idx_eye1_q1 = find(bin_mask_vtx_eye1{h}(:,q1) == 1);
            idx_eye1_q2 = find(bin_mask_vtx_eye1{h}(:,q2) == 1);
            idx_eye2_q1 = find(bin_mask_vtx_eye2{h}(:,q1) == 1);
            idx_eye2_q2 = find(bin_mask_vtx_eye2{h}(:,q2) == 1);
            
            % Filter the correlation matrix based on these rows and columns
            corr_v1_patch_eye1eye1 = corr_v1_patch_mean(idx_eye1_q1, idx_eye1_q2);
            corr_v1_patch_eye2eye2 = corr_v1_patch_mean(idx_eye2_q1, idx_eye2_q2);
            corr_v1_patch_eye1eye2 = [corr_v1_patch_mean(idx_eye1_q1, idx_eye2_q2); corr_v1_patch_mean(idx_eye1_q2, idx_eye2_q1)];

            % filter distance matrix for distances of vertices eye1 with eye1, eye2 with eye2 and eye1 with eye2
            dist_v1_patch_eye1eye1 = dist_v1_patch{h}(idx_eye1_q1, idx_eye1_q2);
            dist_v1_patch_eye2eye2 = dist_v1_patch{h}(idx_eye2_q1, idx_eye2_q2);
            dist_v1_patch_eye1eye2 = [dist_v1_patch{h}(idx_eye1_q1, idx_eye2_q2); dist_v1_patch{h}(idx_eye1_q2, idx_eye2_q1)];

            % concatenate eye1eye1 and eye2eye2
            % pad matrices for concatenation
            max_rows = max(size(corr_v1_patch_eye1eye1), size(corr_v1_patch_eye2eye2));
            max_cols = max(size(corr_v1_patch_eye1eye1), size(corr_v1_patch_eye2eye2));
            padded_eye1 = NaN(max_rows(1), max_cols(1));
            padded_eye2 = NaN(max_rows(1), max_cols(1));
            padded_eye1(1:size(corr_v1_patch_eye1eye1), 1:size(corr_v1_patch_eye1eye1)) = corr_v1_patch_eye1eye1;
            padded_eye2(1:size(corr_v1_patch_eye2eye2), 1:size(corr_v1_patch_eye2eye2)) = corr_v1_patch_eye2eye2;
            corr_v1_patch_same_eye = vertcat(padded_eye1, padded_eye2);
            max_rows = max(size(dist_v1_patch_eye1eye1), size(dist_v1_patch_eye2eye2));
            max_cols = max(size(dist_v1_patch_eye1eye1), size(dist_v1_patch_eye2eye2));
            padded_eye1 = NaN(max_rows(1), max_cols(1));
            padded_eye2 = NaN(max_rows(1), max_cols(1));
            padded_eye1(1:size(dist_v1_patch_eye1eye1), 1:size(dist_v1_patch_eye1eye1)) = dist_v1_patch_eye1eye1;
            padded_eye2(1:size(dist_v1_patch_eye2eye2), 1:size(dist_v1_patch_eye2eye2)) = dist_v1_patch_eye2eye2;
            dist_v1_patch_same_eye = vertcat(padded_eye1, padded_eye2);
            
            % subsampling
            % quantiled for alike and unalike together
            dist_v1_patch_concatenated = vertcat(dist_v1_patch_same_eye(:), dist_v1_patch_eye1eye2(:));
            dist_quantiles_high = [3 quantile(dist_v1_patch_concatenated(dist_v1_patch_concatenated>3), 100-1) inf];
            % number of vertices per quantiles
            quants_number_eye1eye2 = arrayfun(@(i) size(dist_v1_patch_eye1eye2(dist_v1_patch_eye1eye2(:)>dist_quantiles_high(i) & dist_v1_patch_eye1eye2(:)<=dist_quantiles_high(i+1)), 1), 1:100);
            quants_number_same_eye = arrayfun(@(i) size(dist_v1_patch_same_eye(dist_v1_patch_same_eye(:)>dist_quantiles_high(i) & dist_v1_patch_same_eye(:)<=dist_quantiles_high(i+1)), 1), 1:100);
            
            % quantile data
            dist_v1_patch_eye1eye2_quants = cell(100, 1); 
            dist_v1_patch_same_eye_quants = cell(100, 1); 
            corr_v1_patch_eye1eye2_quants = cell(100, 1); 
            corr_v1_patch_same_eye_quants = cell(100, 1); 

            for i = 1:100
                % Filter data within the quantile range
                dist_v1_patch_eye1eye2_quants{i} = ...
                    dist_v1_patch_eye1eye2(dist_v1_patch_eye1eye2 > dist_quantiles_high(i) & ...
                                            dist_v1_patch_eye1eye2 <= dist_quantiles_high(i + 1));
                dist_v1_patch_same_eye_quants{i} = ...
                    dist_v1_patch_same_eye(dist_v1_patch_same_eye > dist_quantiles_high(i) & ...
                                            dist_v1_patch_same_eye<= dist_quantiles_high(i + 1));
                % Filter data within the quantile range
                corr_v1_patch_eye1eye2_quants{i} = ...
                    corr_v1_patch_eye1eye2(dist_v1_patch_eye1eye2 > dist_quantiles_high(i) & ...
                                            dist_v1_patch_eye1eye2 <= dist_quantiles_high(i + 1));
                corr_v1_patch_same_eye_quants{i} = ...
                    corr_v1_patch_same_eye(dist_v1_patch_same_eye > dist_quantiles_high(i) & ...
                                            dist_v1_patch_same_eye <= dist_quantiles_high(i + 1));
            end

            % subsampling
            dist_v1_patch_eye1eye2_subsampled = cell(100, 1);
            dist_v1_patch_same_eye_subsampled = cell(100, 1);
            corr_v1_patch_eye1eye2_subsampled = cell(100, 1);
            corr_v1_patch_same_eye_subsampled = cell(100, 1);
            
            for i = 1:100

                if quants_number_eye1eye2(i) > quants_number_same_eye(i)
                    dist_v1_patch_eye1eye2_subsampled{i} = dist_v1_patch_eye1eye2_quants{i}(selected_indices_eye1eye2{i,h,q1,q2});
                    corr_v1_patch_eye1eye2_subsampled{i} = corr_v1_patch_eye1eye2_quants{i}(selected_indices_eye1eye2{i,h,q1,q2});
                    
                    dist_v1_patch_same_eye_subsampled{i} = dist_v1_patch_same_eye_quants{i};
                    corr_v1_patch_same_eye_subsampled{i} = corr_v1_patch_same_eye_quants{i};
                else
                    dist_v1_patch_eye1eye2_subsampled{i} = dist_v1_patch_eye1eye2_quants{i};
                    corr_v1_patch_eye1eye2_subsampled{i} = corr_v1_patch_eye1eye2_quants{i};
        
                    dist_v1_patch_same_eye_subsampled{i} = dist_v1_patch_same_eye_quants{i}(selected_indices_same_eye{i,h,q1,q2});
                    corr_v1_patch_same_eye_subsampled{i} = corr_v1_patch_same_eye_quants{i}(selected_indices_same_eye{i,h,q1,q2});
                    
                end
        
            end

            dist_v1_patch_eye1eye2_subsampled_concat = vertcat(dist_v1_patch_eye1eye2_subsampled{:});
            dist_v1_patch_same_eye_subsampled_concat = vertcat(dist_v1_patch_same_eye_subsampled{:});

            corr_v1_patch_eye1eye2_subsampled_concat = vertcat(corr_v1_patch_eye1eye2_subsampled{:});
            corr_v1_patch_same_eye_subsampled_concat = vertcat(corr_v1_patch_same_eye_subsampled{:});

            % 10 distance quantiles
            
            dist_concat = vertcat(dist_v1_patch_eye1eye2_subsampled_concat, dist_v1_patch_same_eye_subsampled_concat);
            corr_concat = vertcat(corr_v1_patch_eye1eye2_subsampled_concat, corr_v1_patch_same_eye_subsampled_concat);
            
            dist_quantiles_10_fin = [3 quantile(dist_concat(dist_concat>3), 10-1) inf];

            quants_corr_v1_patch_both{h} = dist_quantiles_10_fin;
            mean_quants_corr_v1_patch_both = arrayfun(@(i) mean(corr_concat(dist_concat>quants_corr_v1_patch_both{h}(i) & dist_concat<=quants_corr_v1_patch_both{h}(i+1)), "omitnan"), 1:nQuant);
            
            % saving distances for distance and beta quantile
            mean_quants_dist_v1_patch_both = arrayfun(@(i) mean(dist_concat(dist_concat>quants_corr_v1_patch_both{h}(i) & dist_concat<=quants_corr_v1_patch_both{h}(i+1)), "omitnan"), 1:nQuant);
            median_quants_dist_v1_patch_both = arrayfun(@(i) median(dist_concat(dist_concat>quants_corr_v1_patch_both{h}(i) & dist_concat<=quants_corr_v1_patch_both{h}(i+1)), "omitnan"), 1:nQuant);
            
            AnalysisParam.mean_quants_dist{h}(q1,q2,:) = mean_quants_dist_v1_patch_both;
            AnalysisParam.median_quants_dist{h}(q1,q2,:) = median_quants_dist_v1_patch_both;
            
            Data_Combined(1:nQuant, 1, h, q1, q2, 1) = mean_quants_corr_v1_patch_both';
            
            % z-transformation
            mean_quants_corr_v1_patch_both_z = 0.5 * log((1 + mean_quants_corr_v1_patch_both) ./ (1 - mean_quants_corr_v1_patch_both));
            Data_Combined(1:nQuant, 1, h, q1, q2, 2) = mean_quants_corr_v1_patch_both_z';

            %clear corr_v1_patch_mean idx_eye1_q1 idx_eye1_q2 idx_eye2_q1 idx_eye2_q2 corr_v1_patch_eye1eye1 corr_v1_patch_eye2eye2 corr_v1_patch_eye1eye2 dist_v1_patch_eye1eye1 dist_v1_patch_eye2eye2 dist_v1_patch_eye1eye2 corr_v1_patch_eye1eye1 corr_v1_patch_eye2eye2 corr_v1_patch_eye1eye2 corr_v1_patch_same_eye corr_v1_patch_same_eye dist_v1_patch_same_eye mean_quants_corr_v1_patch_eye1eye1 mean_quants_corr_v1_patch_eye2eye2 mean_quants_corr_v1_patch_eye1eye2 mean_quants_corr_v1_patch_same_eye mean_quants_corr_v1_patch_eye1eye1_z mean_quants_corr_v1_patch_eye2eye2_z mean_quants_corr_v1_patch_eye1eye2_z mean_quants_corr_v1_patch_same_eye_z
        end
     end
end

if ~exist(sprintf('%s/CorrelationMtx_FC_beta_subsampled.mat', Root))
    save(sprintf('%s/CorrelationMtx_FC_beta_subsampled.mat', Root), 'Data_Combined');
    disp('Correlation matrices saved successfully.');
else
    disp('Correlation matrix file already exists. Skipping save.');
end

if ~exist(sprintf('%s/CorrelationMtx_FC_beta_Params_subsampled.mat', Root))
    save(sprintf('%s/CorrelationMtx_FC_beta_Params_subsampled.mat', Root), 'AnalysisParam', '-v7.3');
    disp('Correlation analysis parameters saved successfully.');
else
    disp('Correlation analysis parameters file already exists. Skipping save.');
end

fprintf('Done in %s s!\n\n', num2str(toc- StrtTime))

fprintf('Total time elapsed: %s s!\n\n', num2str(toc))