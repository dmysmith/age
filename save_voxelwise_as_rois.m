%% Load Voxelwise JA data
data_dir = '/space/syn65/1/data/abcd-sync/6.0/imaging_concat/voxelwise';

volinfo = load(fullfile(data_dir,"smri","volinfo.mat"));
JA = load(fullfile(data_dir, 'smri', 'volmat_JA.mat'));

vol = fullvol(JA.volmat, volinfo.vol_mask_aseg); % error

load('/home/dale/tmp/atlas_dspace/atlas_dspace_ABCD3_cor10.mat');

% % debugging...
% atlas = load('/space/amdale/1/tmp/ABCD_cache/ABCD-Atlas/showVolData/Atlas/5.0_ABCD3_cor10/showVolAtlases_5.0_ABCD3_cor10.mat');

%% Load ROI data
% roi_dir = '/space/syn65/1/data/abcd-sync/6.0/imaging_concat/voxelwise/roi';
% 
thalamus = load(fullfile(data_dir, 'roi', "Thalamus_ABCD3_cor10.mat"));
pauli = load(fullfile(data_dir, 'roi', "Pauli_SubcortNuclei_ABCD3_cor10.mat"));
% aseg = load(fullfile(data_dir, 'roi', "volmat_aseg_prob200x200x260.mat"));

vol = subsample_volume(thalamus.vol_labels_reg); % subsample volumes from atlas file

thalamus_prob = zeros(atlas.thalamus.prob{1}, 'single');
thalamus_prob(atlas.thalamus.prob{2}) = atlas.thalamus.prob{3};
thalamus_prob = thalamus_prob(1:2:end, 1:2:end, 1:2:end, :);

% Create output table
roi_sum = struct('sum', [], 'cats', [], 'region_size', []);
roi_vox = struct('vox', [], 'cats', []);

% run function

[roi_sum, roi_vox] = extract_roi_val(vol, atlas.thalamus.roinames, 0.5, atlas.thalamus.roinames, thalamus_prob, roi_sum, roi_vox);


% 
% %% Load basis functions
% tbl_bf = readtable('/space/syn50/1/data/ABCD/d9smith/age/basis.txt'); agevals_tbl = linspace(100,220,121);
% bfmat = table2array(tbl_bf);
% [dummy dbfmat] = gradient(bfmat); dbfmat = dbfmat/(agevals_tbl(2)-agevals_tbl(1));
% % figure; plot(agevals_tbl,bfmat,'LineWidth',2); 
% % figure; plot(agevals_tbl,dbfmat,'LineWidth',2);
% 
% %% Intersect data
% vec1 = strcat(iid,eid);
% vec2 = strcat(tbl_design.src_subject_id,tbl_design.eventname);
% [dummy IA IB] = intersect(vec1,vec2,'stable');
% agevec_tbl = tbl_design.age;
% agevec = agevec_tbl(IB);
% 
% jvec_bf = find(find_cell(regexp(colnames_model,'^bf_','match')));
% jvec_int = find(find_cell(regexp(colnames_model,'^intercept','match')));
% % figure; plot(agevec,X(:,jvec_bf),'*');
% % figure; plot(agevec, sum(X(:,jvec_bf),2), '*');
% 
% intercept = find(find_cell(regexp(colnames_model,'^intercept','match')));
% 
% %% Calculate weighted sum of betas, plus intercept
% vols = NaN([size(mask) length(agevals_tbl)]);
% W = pinv(X);
% cov_beta = W*W'; % "quick and dirty" version - should get cov(beta) from FEMA
% sig2beta = NaN(size(X,2),1);
% for agei = 1:length(agevals_tbl)
%   wvec = bfmat(agei,:)'; % NOTE - we are NOT taking the derivative here
%   valvec = sum(beta_hat([jvec_bf],:).*wvec,1) + beta_hat(intercept,:); % add intercept
%   vols(:,:,:,agei) = fullvol(valvec,mask);
% 
%   % calculate confidence interval
%   % wvec = [bfmat(agei,:) 1]';
%   % sig2beta(agei) = wvec * cov_beta([jvec_bf jvec_int]) * wvec'; % this is a proportion - need to multiply by sig2tvec
% end

%% Set regions of interest
% incl_thalamus_rois = {'Pulvinar', 'Anterior', 'Medio_Dorsal', 'Ventral_Latero_Dorsal','Central_Lateral-Lateral_Posterior-Medial_Pulvinar','Ventral_Anterior','Ventral_Latero_Ventral'};
% incl_pauli_rois = {'Caudate', 'Putamen', 'Nucleus Accumbens', 'Extended Amygdala', 'Globus Pallidus external', 'Globus Pallidus internal', 'Subthalamic Nucleus', 'Red Nucleus', 'Substantia Nigra pars compacta ', 'Substantia Nigra pars reticulata'};
% incl_aseg_rois = {'Cerebellum-White-Matter', 'Cerebellum-Cortex', 'Thalamus-Proper', 'Caudate', 'Putamen', 'Pallidum', 'Hippocampus', 'Amygdala', 'Accumbens-area', 'VentralDC'};

% %% Create stat matrix
% volmat = [];
% for voli = 1:size(vols,4)
%     disp(sprintf('Iteration %d of %d',voli,size(vols,4))); 
%     [vol_roi_sum(voli), vol_roi_vox(voli)] = vol_by_roi(vols(:,:,:,voli),'atlas','ABCD3');
%     volmat(:,voli) = vol_roi_vox(voli).vox;
% end
% groups = vol_roi_vox(1).cats;

% %% Run function
% % aseg is in ABCD3, thalamus and pauli are in ABCD2
% plot_rois(volmat, agevals_tbl/12, incl_aseg_rois, fstem_imaging); 
% % plot_rois(volmat, agevals_tbl/12, incl_thalamus_rois, fstem_imaging);
% % plot_rois(volmat, agevals_tbl/12, incl_pauli_rois, fstem_imaging);


%% Save newly created tabulated data