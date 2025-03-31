function Data_Combined = B001b_Selectivity_Proc_Data_sub2sub(subName, Root, AnalysisParam, TrgFile)

% This script corresponds to Analysis B - The effect of ocular preference strength (beta),
% cortical depth, ROI (V1 subregions) (and type) on rs-FC (and selectivity).
% This function loads anatomical data, ODC (ocular dominance column) maps, 
% and resting-state fMRI data for a given subject. It then processes this data 
% by applying detrending and high-pass filtering (if specified), computes 
% partial correlations between the specified ROIs (V1 subregions), and 
% calculates the distance matrix between vertices in V1. The function groups 
% distances and beta values into quantiles and computes the correlation between resting-state 
% time series for different quantiles of vertex  and beta. The results are saved 
% for further analysis.

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
% number of quantiles for grouping of distances from selected vertex
nQuant = AnalysisParam.nQuant_beta ;
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
folder = '/autofs/space/ardebil_001/users/Others/Marianna/FC_7T_Coronal/Controls/Data_Denoised'
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
    if strcmp(AnalysisParam.ROI, 'V1_Dorsal_2_Ventral')
        % Dorsal
        roi_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Upsampled_Rtopy.label']));
        V1_Dorsal_label = read_ROIlabel(fullfile(labelFolder, [hemi '.Dorsal_Upsampled_Rtopy.label']));
        [~, idx] = ismember(roi_label, V1_Dorsal_label);
        roi_label_1 = roi_label(idx ~= 0, :);
        % Ventral
        unique_center_labels = unique(V1_Dorsal_label);
        rows_to_keep = ~ismember(roi_label, unique_center_labels);
        roi_label_2 = roi_label(rows_to_keep, :);
    elseif strcmp(AnalysisParam.ROI, 'V1_Posterior_2_Anterior')
        % Posterior
        roi_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Upsampled_Rtopy.label']));
        V1_Posterior_label = read_ROIlabel(fullfile(labelFolder, [hemi '.Posterior_Upsampled_Rtopy.label']));
        [~, idx] = ismember(roi_label, V1_Posterior_label);
        roi_label_1 = roi_label(idx ~= 0, :);
        % Anterior
        unique_center_labels = unique(V1_Posterior_label);
        rows_to_keep = ~ismember(roi_label, unique_center_labels);
        roi_label_2 = roi_label(rows_to_keep, :);
    end

    % crop occipital patch to V1 patch
    v1_patch{h} = struct('npts','ind','x','y','z','vno');
    v1_patch{h}.ind = occ_patch{h}.ind(ismember(occ_patch{h}.ind, roi_label));
    v1_patch{h}.x = occ_patch{h}.x(ismember(occ_patch{h}.ind, roi_label));
    v1_patch{h}.y = occ_patch{h}.y(ismember(occ_patch{h}.ind, roi_label));
    v1_patch{h}.z = occ_patch{h}.z(ismember(occ_patch{h}.ind, roi_label));
    v1_patch_1{h}.x = v1_patch{h}.x(ismember(v1_patch{h}.ind, roi_label_1));
    v1_patch_1{h}.y = v1_patch{h}.y(ismember(v1_patch{h}.ind, roi_label_1));
    v1_patch_1{h}.z = v1_patch{h}.z(ismember(v1_patch{h}.ind, roi_label_1));
    v1_patch_2{h}.x = v1_patch{h}.x(ismember(v1_patch{h}.ind, roi_label_2));
    v1_patch_2{h}.y = v1_patch{h}.y(ismember(v1_patch{h}.ind, roi_label_2));
    v1_patch_2{h}.z = v1_patch{h}.z(ismember(v1_patch{h}.ind, roi_label_2));
    v1_patch{h}.vno = occ_patch{h}.vno(ismember(occ_patch{h}.ind, roi_label));
    v1_patch{h}.npts = length(v1_patch{h}.ind); 

    % create binary mask of vertices that are included in the roi
    bin_mask_vtx_occ_patch{h} = false(size(vtx{h},1),1); bin_mask_vtx_occ_patch{h}(occ_patch{h}.ind) = true;
    bin_mask_vtx_v1_patch{h} = false(size(vtx{h},1),1); 
    bin_mask_vtx_v1_patch_1{h} = false(size(vtx{h},1),1);
    bin_mask_vtx_v1_patch_2{h} = false(size(vtx{h},1),1);  
    bin_mask_vtx_v1_patch{h}(v1_patch{h}.ind) = true;
    bin_mask_vtx_v1_patch_1{h}(roi_label_1) = true;
    bin_mask_vtx_v1_patch_2{h}(roi_label_2) = true;

    % roi1
    % create binary masks for the two eyes just within V1 based on quantiles    
    % Extract masked values
    rows_to_keep = bin_mask_vtx_v1_patch{h} == 1;
    bin_mask_vtx_v1_patch_1{h} = bin_mask_vtx_v1_patch_1{h}(rows_to_keep, :); 
    bin_mask_vtx_v1_patch_2{h} = bin_mask_vtx_v1_patch_2{h}(rows_to_keep, :); 
    
    masked_values = sig_odc{h}(bin_mask_vtx_v1_patch_1{h});   
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
        bin_mask_vtx_eye1_roi1{h}(:, quantile_beta_threshold) = false(size(roi_label,1), 1);
        bin_mask_vtx_eye2_roi1{h}(:, quantile_beta_threshold) = false(size(roi_label,1), 1);
        
        % For the left eye: Select vertices in the most negative quantiles
        bin_mask_vtx_eye1_roi1{h}(negative_indices(selected_neg_mask), quantile_beta_threshold) = true;
    
        % For the right eye: Select vertices in the most positive quantiles
        bin_mask_vtx_eye2_roi1{h}(positive_indices(selected_pos_mask), quantile_beta_threshold) = true;
        
        % remove unnecessary empty cells from the mask
        bin_mask_vtx_eye1_roi1{h} = bin_mask_vtx_eye1_roi1{h}(1:size(roi_label,1),:);
        bin_mask_vtx_eye2_roi1{h} = bin_mask_vtx_eye2_roi1{h}(1:size(roi_label,1),:);
        
        AnalysisParam.num_included_vertices_eye1{h,quantile_beta_threshold} = sum(bin_mask_vtx_eye1_roi1{h});
        AnalysisParam.num_included_vertices_eye2{h,quantile_beta_threshold} = sum(bin_mask_vtx_eye2_roi1{h});
    end

    for q1 = 1:numel(quantile_beta_thresholds) 
        for q2 = 1:numel(quantile_beta_thresholds)
            % sometimes quantiles include one vertex more than the other, so one random vertex is excluded
            if sum(bin_mask_vtx_eye1_roi1{h}(:,q1)) > sum(bin_mask_vtx_eye1_roi1{h}(:,q2))
                column_data = bin_mask_vtx_eye1_roi1{h}(:, q1);
                indices = find(column_data == 1);
                random_index = indices(randi(numel(indices)));
                column_data(random_index) = 0;
                bin_mask_vtx_eye1_roi1{h}(:, q1) = column_data;
            elseif sum(bin_mask_vtx_eye1_roi1{h}(:,q1)) < sum(bin_mask_vtx_eye1_roi1{h}(:,q2))
                column_data = bin_mask_vtx_eye1_roi1{h}(:, q2);
                indices = find(column_data == 1);
                random_index = indices(randi(numel(indices)));
                column_data(random_index) = 0;
                bin_mask_vtx_eye1_roi1{h}(:, q2) = column_data;
            end
            if sum(bin_mask_vtx_eye2_roi1{h}(:,q1)) > sum(bin_mask_vtx_eye2_roi1{h}(:,q2))
                column_data = bin_mask_vtx_eye2_roi1{h}(:, q1);
                indices = find(column_data == 1);
                random_index = indices(randi(numel(indices)));
                column_data(random_index) = 0;
                bin_mask_vtx_eye2_roi1{h}(:, q1) = column_data;
            elseif sum(bin_mask_vtx_eye2_roi1{h}(:,q1)) < sum(bin_mask_vtx_eye2_roi1{h}(:,q2))
                column_data = bin_mask_vtx_eye2_roi1{h}(:, q2);
                indices = find(column_data == 1);
                random_index = indices(randi(numel(indices)));
                column_data(random_index) = 0;
                bin_mask_vtx_eye2_roi1{h}(:, q2) = column_data;
            end 
        end
    end

    % roi2
    % create binary masks for the two eyes just within V1 based on quantiles    
    % Extract masked values
    masked_values = sig_odc{h}(bin_mask_vtx_v1_patch_2{h});   
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
        bin_mask_vtx_eye1_roi2{h}(:, quantile_beta_threshold) = false(size(roi_label,1), 1);
        bin_mask_vtx_eye2_roi2{h}(:, quantile_beta_threshold) = false(size(roi_label,1), 1);
        
        % For the left eye: Select vertices in the most negative quantiles
        bin_mask_vtx_eye1_roi2{h}(negative_indices(selected_neg_mask), quantile_beta_threshold) = true;
    
        % For the right eye: Select vertices in the most positive quantiles
        bin_mask_vtx_eye2_roi2{h}(positive_indices(selected_pos_mask), quantile_beta_threshold) = true;
        
        % remove unnecessary empty cells from the mask
        bin_mask_vtx_eye1_roi2{h} = bin_mask_vtx_eye1_roi2{h}(1:size(roi_label,1),:);
        bin_mask_vtx_eye2_roi2{h} = bin_mask_vtx_eye2_roi2{h}(1:size(roi_label,1),:);
        
        AnalysisParam.num_included_vertices_eye1{h,quantile_beta_threshold} = sum(bin_mask_vtx_eye1_roi2{h});
        AnalysisParam.num_included_vertices_eye2{h,quantile_beta_threshold} = sum(bin_mask_vtx_eye2_roi2{h});
    end

    for q1 = 1:numel(quantile_beta_thresholds) 
        for q2 = 1:numel(quantile_beta_thresholds)
            % sometimes quantiles include one vertex more than the other, so one random vertex is excluded
            if sum(bin_mask_vtx_eye1_roi2{h}(:,q1)) > sum(bin_mask_vtx_eye1_roi2{h}(:,q2))
                column_data = bin_mask_vtx_eye1_roi2{h}(:, q1);
                indices = find(column_data == 1);
                random_index = indices(randi(numel(indices)));
                column_data(random_index) = 0;
                bin_mask_vtx_eye1_roi2{h}(:, q1) = column_data;
            elseif sum(bin_mask_vtx_eye1_roi2{h}(:,q1)) < sum(bin_mask_vtx_eye1_roi2{h}(:,q2))
                column_data = bin_mask_vtx_eye1_roi2{h}(:, q2);
                indices = find(column_data == 1);
                random_index = indices(randi(numel(indices)));
                column_data(random_index) = 0;
                bin_mask_vtx_eye1_roi2{h}(:, q2) = column_data;
            end
            if sum(bin_mask_vtx_eye2_roi2{h}(:,q1)) > sum(bin_mask_vtx_eye2_roi2{h}(:,q2))
                column_data = bin_mask_vtx_eye2_roi2{h}(:, q1);
                indices = find(column_data == 1);
                random_index = indices(randi(numel(indices)));
                column_data(random_index) = 0;
                bin_mask_vtx_eye2_roi2{h}(:, q1) = column_data;
            elseif sum(bin_mask_vtx_eye2_roi2{h}(:,q1)) < sum(bin_mask_vtx_eye2_roi2{h}(:,q2))
                column_data = bin_mask_vtx_eye2_roi2{h}(:, q2);
                indices = find(column_data == 1);
                random_index = indices(randi(numel(indices)));
                column_data(random_index) = 0;
                bin_mask_vtx_eye2_roi2{h}(:, q2) = column_data;
            end 
        end
    end

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
        %dist_v1_patch{h} = sqrt((v1_patch{h}.x - v1_patch{h}.x').^2 + (v1_patch{h}.y - v1_patch{h}.y').^2 + (v1_patch{h}.z - v1_patch{h}.z').^2);    
        dist_v1_patch{h} = sqrt((v1_patch_2{h}.x - v1_patch_1{h}.x').^2 + (v1_patch_2{h}.y - v1_patch_1{h}.y').^2 + (v1_patch_2{h}.z - v1_patch_1{h}.z').^2); 

end
    
    % save([Root,'/DistanceMtx.mat'], 'dist_v1_patch_shuffled', '-v7.3')   
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
            rs_v1_patch_1{h} = [rs_v1_patch{h} squeeze(rs(bin_mask_vtx_v1_patch_1{h},:))];
            rs_v1_patch_2{h} = [rs_v1_patch{h} squeeze(rs(bin_mask_vtx_v1_patch_2{h},:))];
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
        corr_v1_patch = partialcorr(rs_v1_patch_1{h}', rs_v1_patch_2{h}', [mcpr(:,1:3) wm]);
        % Store the correlation matrix for the current quantile
        corr_v1_patch_allSess{h}(sessionNum,:,:) = corr_v1_patch;
        clear corr_v1_patch
    end
    fprintf('Done in %s s!\n\n', num2str(toc- StrtTime))        
end

%% do analysis for mean across sessions

StrtTime = toc;

Data_Combined = zeros(nQuant, 4, 2, nQuant, nQuant,2);
for h = 1:numel(hemis)
     for q1 = 1:numel(quantile_beta_thresholds)  % Assuming quantile_beta_threshold contains thresholds
        for q2 = 1:numel(quantile_beta_thresholds)
            hemi = hemis{h};
        
            corr_v1_patch_mean = squeeze(mean(corr_v1_patch_allSess{h}(:,:,:), 1));

            idx_eye1_q1_roi1 = find(bin_mask_vtx_eye1_roi1{h}(:,q1) == 1);
            idx_eye1_q1_roi2 = find(bin_mask_vtx_eye1_roi2{h}(:,q1) == 1);
            idx_eye1_q2_roi1 = find(bin_mask_vtx_eye1_roi1{h}(:,q2) == 1);
            idx_eye1_q2_roi2 = find(bin_mask_vtx_eye1_roi2{h}(:,q2) == 1);
            idx_eye2_q1_roi1 = find(bin_mask_vtx_eye2_roi1{h}(:,q1) == 1);
            idx_eye2_q1_roi2 = find(bin_mask_vtx_eye2_roi2{h}(:,q1) == 1);
            idx_eye2_q2_roi1 = find(bin_mask_vtx_eye2_roi1{h}(:,q2) == 1);
            idx_eye2_q2_roi2 = find(bin_mask_vtx_eye2_roi2{h}(:,q2) == 1);

            % Filter the correlation matrix based on these rows and columns
            corr_v1_patch_eye1eye1 = [corr_v1_patch_mean(idx_eye1_q1_roi1, idx_eye1_q2_roi2); corr_v1_patch_mean(idx_eye1_q2_roi1, idx_eye1_q1_roi2)];
            corr_v1_patch_eye2eye2 = [corr_v1_patch_mean(idx_eye2_q1_roi1, idx_eye2_q2_roi2); corr_v1_patch_mean(idx_eye2_q2_roi1, idx_eye2_q1_roi2)];
            corr_v1_patch_eye1eye2_1 = [corr_v1_patch_mean(idx_eye1_q1_roi1, idx_eye2_q2_roi2); corr_v1_patch_mean(idx_eye1_q2_roi1, idx_eye2_q1_roi2)]; 
            corr_v1_patch_eye1eye2_2 = [corr_v1_patch_mean(idx_eye1_q1_roi1, idx_eye1_q2_roi2); corr_v1_patch_mean(idx_eye2_q2_roi1, idx_eye1_q1_roi2)];
            
            % Filter the distance matrix based on these rows and columns
            dist_v1_patch_eye1eye1 = [dist_v1_patch{h}(idx_eye1_q1_roi1, idx_eye1_q2_roi2); dist_v1_patch{h}(idx_eye1_q2_roi1, idx_eye1_q1_roi2)];
            dist_v1_patch_eye2eye2 = [dist_v1_patch{h}(idx_eye2_q1_roi1, idx_eye2_q2_roi2); dist_v1_patch{h}(idx_eye2_q2_roi1, idx_eye2_q1_roi2)];
            dist_v1_patch_eye1eye2_1 = [dist_v1_patch{h}(idx_eye1_q1_roi1, idx_eye2_q2_roi2); dist_v1_patch{h}(idx_eye1_q2_roi1, idx_eye2_q1_roi2)]; 
            dist_v1_patch_eye1eye2_2 = [dist_v1_patch{h}(idx_eye1_q1_roi1, idx_eye1_q2_roi2); dist_v1_patch{h}(idx_eye2_q2_roi1, idx_eye1_q1_roi2)];

            % Pad corr_v1_patch_eye1eye2_1 and corr_v1_patch_eye1eye2_2 for concatenation
            [size1_r, size1_c] = size(corr_v1_patch_eye1eye2_1);
            [size2_r, size2_c] = size(corr_v1_patch_eye1eye2_2);
            max_rows = max(size1_r, size2_r);
            max_cols = max(size1_c, size2_c);
            corr_v1_patch_eye1eye2_1_padded = NaN(max_rows, max_cols);
            corr_v1_patch_eye1eye2_1_padded(1:size1_r, 1:size1_c) = corr_v1_patch_eye1eye2_1;
            corr_v1_patch_eye1eye2_2_padded = NaN(max_rows, max_cols);
            corr_v1_patch_eye1eye2_2_padded(1:size2_r, 1:size2_c) = corr_v1_patch_eye1eye2_2;
            corr_v1_patch_eye1eye2 = [corr_v1_patch_eye1eye2_1_padded; corr_v1_patch_eye1eye2_2_padded];
            % same for distances
            [size1_r, size1_c] = size(dist_v1_patch_eye1eye2_1);
            [size2_r, size2_c] = size(dist_v1_patch_eye1eye2_2);
            max_rows = max(size1_r, size2_r);
            max_cols = max(size1_c, size2_c);
            dist_v1_patch_eye1eye2_1_padded = NaN(max_rows, max_cols);
            dist_v1_patch_eye1eye2_1_padded(1:size1_r, 1:size1_c) = dist_v1_patch_eye1eye2_1;
            dist_v1_patch_eye1eye2_2_padded = NaN(max_rows, max_cols);
            dist_v1_patch_eye1eye2_2_padded(1:size2_r, 1:size2_c) = dist_v1_patch_eye1eye2_2;
            dist_v1_patch_eye1eye2 = [dist_v1_patch_eye1eye2_1_padded; dist_v1_patch_eye1eye2_2_padded];

            % concatenate eye1eye1 and eye2eye2
            % pad matrices for concatenation
            [size1_r, size1_c] = size(corr_v1_patch_eye1eye1);
            [size2_r, size2_c] = size(corr_v1_patch_eye2eye2);
            max_rows = max(size1_r, size2_r);
            max_cols = max(size1_c, size2_c);
            corr_v1_patch_eye1eye1_padded = NaN(max_rows, max_cols);
            corr_v1_patch_eye1eye1_padded(1:size1_r, 1:size1_c) = corr_v1_patch_eye1eye1;
            corr_v1_patch_eye2eye2_padded = NaN(max_rows, max_cols);
            corr_v1_patch_eye2eye2_padded(1:size2_r, 1:size2_c) = corr_v1_patch_eye2eye2;
            corr_v1_patch_same_eye = [corr_v1_patch_eye1eye1_padded; corr_v1_patch_eye2eye2_padded];
            % same for distances
            [size1_r, size1_c] = size(dist_v1_patch_eye1eye1);
            [size2_r, size2_c] = size(dist_v1_patch_eye2eye2);
            max_rows = max(size1_r, size2_r);
            max_cols = max(size1_c, size2_c);
            dist_v1_patch_eye1eye1_padded = NaN(max_rows, max_cols);
            dist_v1_patch_eye1eye1_padded(1:size1_r, 1:size1_c) = dist_v1_patch_eye1eye1;
            dist_v1_patch_eye2eye2_padded = NaN(max_rows, max_cols);
            dist_v1_patch_eye2eye2_padded(1:size2_r, 1:size2_c) = dist_v1_patch_eye2eye2;
            dist_v1_patch_same_eye = [dist_v1_patch_eye1eye1_padded; dist_v1_patch_eye2eye2_padded];
            
            % mean for distance quantiles
            quants_corr_v1_patch_eye1eye1{h} = [distThresh quantile(dist_v1_patch_eye1eye1(dist_v1_patch_eye1eye1>distThresh), nQuant-1) inf];
            mean_quants_corr_v1_patch_eye1eye1 = arrayfun(@(i) mean(corr_v1_patch_eye1eye1(dist_v1_patch_eye1eye1>quants_corr_v1_patch_eye1eye1{h}(i) & dist_v1_patch_eye1eye1<=quants_corr_v1_patch_eye1eye1{h}(i+1)), "omitnan"), 1:nQuant);
            
            quants_corr_v1_patch_eye2eye2{h} = [distThresh quantile(dist_v1_patch_eye2eye2(dist_v1_patch_eye2eye2>distThresh), nQuant-1) inf];
            mean_quants_corr_v1_patch_eye2eye2 = arrayfun(@(i) mean(corr_v1_patch_eye2eye2(dist_v1_patch_eye2eye2>quants_corr_v1_patch_eye2eye2{h}(i) & dist_v1_patch_eye2eye2<=quants_corr_v1_patch_eye2eye2{h}(i+1)), "omitnan"), 1:nQuant);
            
            quants_corr_v1_patch_eye1eye2{h} = [distThresh quantile(dist_v1_patch_eye1eye2(dist_v1_patch_eye1eye2>distThresh), nQuant-1) inf];
            mean_quants_corr_v1_patch_eye1eye2 = arrayfun(@(i) mean(corr_v1_patch_eye1eye2(dist_v1_patch_eye1eye2>quants_corr_v1_patch_eye1eye2{h}(i) & dist_v1_patch_eye1eye2<=quants_corr_v1_patch_eye1eye2{h}(i+1)), "omitnan"), 1:nQuant);

            quants_corr_v1_patch_same_eye{h} = [distThresh quantile(dist_v1_patch_same_eye(dist_v1_patch_same_eye>distThresh), nQuant-1) inf];
            mean_quants_corr_v1_patch_same_eye = arrayfun(@(i) mean(corr_v1_patch_same_eye(dist_v1_patch_same_eye>quants_corr_v1_patch_same_eye{h}(i) & dist_v1_patch_same_eye<=quants_corr_v1_patch_same_eye{h}(i+1)), "omitnan"), 1:nQuant);

            Data_Combined(1:nQuant, 1:4, h, q1, q2, 1) = [mean_quants_corr_v1_patch_same_eye' mean_quants_corr_v1_patch_eye1eye1' mean_quants_corr_v1_patch_eye2eye2' mean_quants_corr_v1_patch_eye1eye2'];
            % z-transformation
            mean_quants_corr_v1_patch_same_eye_z = 0.5 * log((1 + mean_quants_corr_v1_patch_same_eye) ./ (1 - mean_quants_corr_v1_patch_same_eye));
            mean_quants_corr_v1_patch_eye1eye1_z = 0.5 * log((1 + mean_quants_corr_v1_patch_eye1eye1) ./ (1 - mean_quants_corr_v1_patch_eye1eye1));
            mean_quants_corr_v1_patch_eye2eye2_z = 0.5 * log((1 + mean_quants_corr_v1_patch_eye2eye2) ./ (1 - mean_quants_corr_v1_patch_eye2eye2));
            mean_quants_corr_v1_patch_eye1eye2_z = 0.5 * log((1 + mean_quants_corr_v1_patch_eye1eye2) ./ (1 - mean_quants_corr_v1_patch_eye1eye2));
            Data_Combined(1:nQuant, 1:4, h, q1, q2, 2) = [mean_quants_corr_v1_patch_same_eye_z' mean_quants_corr_v1_patch_eye1eye1_z' mean_quants_corr_v1_patch_eye2eye2_z' mean_quants_corr_v1_patch_eye1eye2_z'];
            clear corr_v1_patch_mean idx_eye1_q1 idx_eye1_q2 idx_eye2_q1 idx_eye2_q2 corr_v1_patch_eye1eye1 corr_v1_patch_eye2eye2 corr_v1_patch_eye1eye2 dist_v1_patch_eye1eye1 dist_v1_patch_eye2eye2 dist_v1_patch_eye1eye2 corr_v1_patch_eye1eye1 corr_v1_patch_eye2eye2 corr_v1_patch_eye1eye2 corr_v1_patch_same_eye corr_v1_patch_same_eye dist_v1_patch_same_eye mean_quants_corr_v1_patch_eye1eye1 mean_quants_corr_v1_patch_eye2eye2 mean_quants_corr_v1_patch_eye1eye2 mean_quants_corr_v1_patch_same_eye mean_quants_corr_v1_patch_eye1eye1_z mean_quants_corr_v1_patch_eye2eye2_z mean_quants_corr_v1_patch_eye1eye2_z mean_quants_corr_v1_patch_same_eye_z

        end
     end
end
    
save([Root,'/CorrelationMtx_Selectivity.mat'], 'Data_Combined', 'AnalysisParam')
