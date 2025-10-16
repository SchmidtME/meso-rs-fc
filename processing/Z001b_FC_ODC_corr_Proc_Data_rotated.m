function meanAbsR = Z001b_FC_ODC_corr_Proc_Data_rotated(subName, Root, AnalysisParam, TrgFile)

% Description:
% This script corresponds to Analysis Z - the correlation of rs-FC with the
% ODC differential map to Figure 6 of the manuscript.
% This script calculate the rs-FC for 1000 randomly selected vertices with vertices 
% within a ring/donut of variable radii (and excluding the inner vertices with a distance 
% of < 3mm to the seed vertex. The 1000 2D rs-FC maps are individually correlated with 
% a 180 deg rotated version of the differential ODC map. The 1000 correlation
% coefficients are averaged and saved.
% Authors: Marianna E. Schmidt (marianna.schmidt@maxplanckschools.de), Iman Aganj, Shahin Nasr

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
folder = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Data_Denoised'

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
    roi_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Upsampled_Rtopy.label']));

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
            rs_v1_patch{h} = [rs_v1_patch{h} squeeze(rs(bin_mask_vtx_v1_patch{h},:))];
        end
    end

    % load the covariate data and append for each session
    wm0 = load(fullfile(rsFolder, rsSessions{sessionNumTemp}, 'wm.dat'), '-ascii');
    wm = [wm; wm0(timeSeriesRange{:})];
    mcpr0 = load(fullfile(rsFolder, rsSessions{sessionNumTemp}, 'mcprextreg'), '-ascii');
    mcpr = [mcpr; mcpr0(timeSeriesRange{:},:)];

end

%% Correlation

visualizeResults = false
distEdgeRange = 4:0.5:10;
minDist = 3; % mm
numSamplePointsPermm2 = 10; % Number of Halton sample points per mm2.

for h = 1:numel(hemis)
    maskedSig = sig_odc{h}(bin_mask_vtx_v1_patch{h}); % Significance values inside the patch
    for iDistFromEdge = 1:length(distEdgeRange) % mm
        distFromEdge = distEdgeRange(iDistFromEdge);
        N = round(numSamplePointsPermm2 * pi*(distFromEdge^2 - minDist^2));
        clear circR
        ind = find(sum(dist_v1_patch{h} < distFromEdge, 2) > 6*(pi*distFromEdge^2));

        for iHypoth = 1:2
            for k = 1:1000

                randInd(1) = ind(randi(length(ind)));
                randInd(2) = randInd(1);
                
                % Generate Halton points in 2D
                haltonPoints = (net(haltonset(2), round(8*N*(distFromEdge^2)/(pi*(distFromEdge^2 - minDist^2)))) * 2 - 1) * distFromEdge; %(net(haltonset(2), 3*N) * 2 - 1) * distFromEdge; % Generate 3N 2D Halton sample points in [-1 1].
                distSqr = sum(haltonPoints.^2,2);
                haltonPoints = haltonPoints(distSqr <= distFromEdge^2 & distSqr >= minDist^2, :); % Filter Halton points inside the disc.
                haltonPoints = haltonPoints(1:N, :); % Retain the first N points and scale to the radius r.
                
                % Project the rs-fMRI data (circNum=1) and ODC significance map (circNum=2) to Halton points
                for circNum = 1:2
                    discCenter = vtx_v1_patch{h}(randInd(circNum),1:2);
                    haltonPointsMoved{circNum} = haltonPoints + discCenter; % Adding the disc center.
                    nearestIdx = knnsearch(vtx_v1_patch{h}(:, 1:2), haltonPointsMoved{circNum}); % Map Halton points to nearest vertices

                    if circNum == 1
                        projectedHalton{circNum} = rs_v1_patch{h}(nearestIdx,:); % Map rs-fMRI values
                    else
                        if iHypoth == 1 % Null hypothesis, rotation
                            rotatedPoints = -haltonPoints + discCenter(1:2);
                            nearestIdx = knnsearch(vtx_v1_patch{h}(:, 1:2), rotatedPoints); 
                            projectedHalton{circNum} = maskedSig(nearestIdx); % Map rotated ODC significance values
                        else
                            projectedHalton{circNum} = maskedSig(nearestIdx); % Map ODC values
                        end
                    end
                end

                circCorr = partialcorr(projectedHalton{1}', rs_v1_patch{h}(randInd(1),:)', [mcpr(:,1:3) wm]);
                circR(k) = corr(circCorr, projectedHalton{2});

                if visualizeResults
                    figure(iHypoth)
                    subplot(1,3,1)
                    scatter(vtx_v1_patch{h}(:,1), vtx_v1_patch{h}(:,2), [], maskedSig, '.')
                    axis equal tight, colormap turbo; colorbar
                    title('Original ODC Significance Map');

                    subplot(1,3,2)
                    scatter(haltonPointsMoved{2}(:, 1), haltonPointsMoved{2}(:, 2), [], projectedHalton{2}, '.')
                    axis equal tight, colormap turbo; colorbar
                    title('Projection of ODC Significance Map to Halton Points');

                    subplot(1,3,3)
                    scatter(haltonPointsMoved{1}(:, 1), haltonPointsMoved{1}(:, 2), [], circCorr, '.')
                    axis equal tight, colormap turbo; colorbar
                    title('Projection of FC to Halton Points');

                    set(gcf, "Name", ['r = ' num2str(circR(k))])
                end
            end
            meanAbsR(iDistFromEdge, iHypoth, h) = meanabs(circR);
        end
    end
end

save([Root,'/meanAbsR_Denoised_excl_3mm_rotated.mat'], 'meanAbsR', 'distEdgeRange', 'AnalysisParam')

clear meanAbsR