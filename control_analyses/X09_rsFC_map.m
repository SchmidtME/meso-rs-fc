clear all
close all
clc

%% specifications

addpath(genpath('/space/ardebil/1/users/Others/Marianna/Code/MesoVision/MESO-FC-Controls'))

saveDir = '/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Figures/Figure_12';
mkdir(saveDir);

subs = {'chss'}%, 'chss', 'auil', 'aroo'};
Layers = {'4-6'}%, '4-6', '8-10'};
hemis = {'lh'}%, 'rh'};
StrtPnt = 10;
TmSeries_Length = 118;
timeSeriesRange = {StrtPnt:(StrtPnt+TmSeries_Length)};

%% porcessing

for l = 1:length(Layers)
    for s = 1:length(subs)
    
        %% specifications
        
        subName = subs{s};
        
        %% Compute correlations
        
        for h = 1:length(hemis)
        
            corrLabel_sessions = [];
            hemi = hemis{h};
            disp(['Hemisphere: ' hemi])
            disp('Loading the data...')
            tic
        
            anatFolder = fullfile('/autofs/space/ardebil_001/users/Shared/good_subjects_anat', [subName '_anat_upsample_B1justsub']);
            labelFolder = fullfile('/space/ardebil/1/users/Others/Marianna/good_subjects_anat/retinotopy', subName);
            odcFolder = fullfile('/autofs/space/ardebil_002/users/Shahin/Stereopsis_Upsampled/Subjects_Control_Upsampled/',  [subName '_Final_Upsampled_VOL']);
            
            rsFolder = fullfile('/space/ardebil/1/users/Others/Marianna/FC_7T_Coronal/Controls/Data_Denoised', subName, 'bold_Close_Upsampled2');

            % load V1 label
            [vtx{h}, fac{h}] = read_surf(fullfile(anatFolder, 'surf', [hemi '.inflated'])); 
            fac{h} = fac{h} + 1;
            occ_patch{h} = read_patch(fullfile(anatFolder, 'surf', sprintf('%s.%s', hemi, 'V1_patch.flat'))); 
            occ_patch{h}.ind = occ_patch{h}.ind+1;
            bin_mask_vtx_occ_patch{h} = false(size(vtx{h},1),1); 
            bin_mask_vtx_occ_patch{h}(occ_patch{h}.ind) = true;    
            roi_label = read_ROIlabel(fullfile(labelFolder, [hemi '.V1_Upsampled_Rtopy.label']));
            v1_patch{h} = struct('npts','ind','x','y','z','vno');
            v1_patch{h}.ind = occ_patch{h}.ind(ismember(occ_patch{h}.ind, roi_label));
            v1_patch{h}.x = occ_patch{h}.x(ismember(occ_patch{h}.ind, roi_label));
            v1_patch{h}.y = occ_patch{h}.y(ismember(occ_patch{h}.ind, roi_label));
            v1_patch{h}.z = occ_patch{h}.z(ismember(occ_patch{h}.ind, roi_label));
            v1_patch{h}.vno = occ_patch{h}.vno(ismember(occ_patch{h}.ind, roi_label));
            v1_patch{h}.npts = length(v1_patch{h}.ind); 
            bin_mask_vtx_v1_patch{h} = false(size(vtx{h},1),1); 
            bin_mask_vtx_v1_patch{h}(v1_patch{h}.ind) = true;
        
            % load odc map
            d = dir(fullfile(odcFolder, '*Smoothing_0-2.lh')).name(1:end-3); 
            sig_odc{h} = load_nifti(fullfile(odcFolder, [d '.' hemi], 'R_C', [hemi '.ffx.osgm.wls'], 'beta.nii')).vol; 
            sigMasked{h} = sig_odc{h}(v1_patch{h}.ind); 
        
            odcSig = tanh((sigMasked{h}-median(sigMasked{h}))/5); % odcSig = tanh(sig{h}(bin_mask_vtx_occ_patch{h})/5); odcSig = odcSig(bin_mask_vtx_v1_patch{h}(bin_mask_vtx_occ_patch{h}));
        
            medianPchX = median(occ_patch{h}.x(bin_mask_vtx_v1_patch{h}(bin_mask_vtx_occ_patch{h})));
            indLeft = occ_patch{h}.x' <= medianPchX & bin_mask_vtx_v1_patch{h}(bin_mask_vtx_occ_patch{h});
            indRight = occ_patch{h}.x' > medianPchX & bin_mask_vtx_v1_patch{h}(bin_mask_vtx_occ_patch{h});    
        
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
                rs = rs(:, [StrtPnt:StrtPnt+TmSeries_Length]);
                rs = detrend(rs', 2)';
                rs = highpass(rs', 0.01, 1/4)';
                rs_v1_patch{h} = [rs_v1_patch{h} squeeze(rs(bin_mask_vtx_v1_patch{h},:))];
            
                % load covariates
                wm0 = load(fullfile(rsFolder, rsSessions{sessionNum}, 'wm.dat'), '-ascii');
                wm = [wm; wm0(timeSeriesRange{:})];
                mcpr0 = load(fullfile(rsFolder, rsSessions{sessionNum}, 'mcprextreg'), '-ascii');
                mcpr = [mcpr; mcpr0(timeSeriesRange{:},:)];
                
                corrLabel_sessions(sessionNum, :, :) = partialcorr(rs_v1_patch{h}', [mcpr(:,1:3) wm]);
        
                toc
            
           end

           corrLabel{h} = squeeze(mean(corrLabel_sessions,1));

        end
    end
end

%%

seed_vertex = 4236;
idx = find(v1_patch{h}.ind == 4236);

% Save corrLabel in the same format as ODC map
corrLabel_full = zeros(size(vtx{h}, 1), 1); % Initialize with zeros
corrLabel_full(v1_patch{h}.ind) = squeeze(corrLabel{h}(idx,:)); % Assign computed values
disp(['seed vertex: ' num2str(v1_patch{h}.ind(idx))])
% Define output file path
savePath = fullfile(saveDir, sprintf('%s_%s_corrLabel_seed_vtx_%d.mgz', subName, hemi,idx));

% Save as .mgz file
M = mghM{h}; % Use transformation matrix from loaded .mgz file

% Save the full-size correlation map in .mgz format
save_mgh(corrLabel_full, savePath, M);
disp(['Saved: ' savePath])