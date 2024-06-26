%% Load results
designmat_dir = '/space/cluster/1/ABCD/users/d9smith/age/results_2024-01-22/';
designmat_file = 'designMat4_BFsSexIncEducHispPCsScanSoftMotion_bly2y4';
% designMat4_BFsSexIncEducHispPCsScanSoftMotion_bly2y4
% designMat4f_BFsSexIncEducHispPCsScanSoftMotion_bly2y4
% designMat4m_BFsSexIncEducHispPCsScanSoftMotion_bly2y4
% designMat5_BFsSexIncEducHispPCsScanSoftMotionPDS_bly2y4
% designMat5f_BFsSexIncEducHispPCsScanSoftMotionPDS_bly2y4
% designMat5m_BFsSexIncEducHispPCsScanSoftMotionPDS_bly2y4

fname_design = strcat(designmat_dir,designmat_file,'.txt');
tbl_design = readtable(fname_design);
% figure; plot(tbl_design.age,tbl_design.bf_demean_1,'*')

dirname_results = strrep(fname_design,'.txt','');
fstem_imaging = 'RNI';
fname_results = sprintf('%s/FEMA_wrapper_output_voxel_%s.mat',dirname_results, fstem_imaging);
load(fname_results);

%% Load basis functions
tbl_bf = readtable('/space/cluster/1/ABCD/users/d9smith/age/basis.txt'); agevals_tbl = linspace(100,220,121);
bfmat = table2array(tbl_bf);
[dummy dbfmat] = gradient(bfmat); dbfmat = dbfmat/(agevals_tbl(2)-agevals_tbl(1));
% figure; plot(agevals_tbl,bfmat,'LineWidth',2); 
% figure; plot(agevals_tbl,dbfmat,'LineWidth',2);

%% Intersect data
vec1 = strcat(iid,eid);
vec2 = strcat(tbl_design.src_subject_id,tbl_design.eventname);
[dummy IA IB] = intersect(vec1,vec2,'stable');
agevec_tbl = tbl_design.age;
agevec = agevec_tbl(IB);

jvec_bf = find(find_cell(regexp(colnames_model,'^bf_','match')));
jvec_int = find(find_cell(regexp(colnames_model,'^intercept','match')));
% figure; plot(agevec,X(:,jvec_bf),'*');
% figure; plot(agevec, sum(X(:,jvec_bf),2), '*');

intercept = find(find_cell(regexp(colnames_model,'^intercept','match')));

%% Calculate weighted sum of betas, plus intercept
vols = NaN([size(mask) length(agevals_tbl)]);
W = pinv(X);
cov_beta = W*W'; % "quick and dirty" version - should get cov(beta) from FEMA
sig2beta = NaN(size(X,2),1);
for agei = 1:length(agevals_tbl)
  wvec = bfmat(agei,:)'; % NOTE - we are NOT taking the derivative here
  valvec = sum(beta_hat([jvec_bf],:).*wvec,1) + beta_hat(intercept,:); % add intercept
  vols(:,:,:,agei) = fullvol(valvec,mask);

  % calculate confidence interval
  % wvec = [bfmat(agei,:) 1]';
  % sig2beta(agei) = wvec * cov_beta([jvec_bf jvec_int]) * wvec'; % this is a proportion - need to multiply by sig2tvec
end

%% Set regions of interest
incl_thalamus_rois = {'Pulvinar', 'Anterior', 'Medio_Dorsal', 'Ventral_Latero_Dorsal','Central_Lateral-Lateral_Posterior-Medial_Pulvinar','Ventral_Anterior','Ventral_Latero_Ventral'};
incl_pauli_rois = {'Caudate', 'Putamen', 'Nucleus Accumbens', 'Extended Amygdala', 'Globus Pallidus external', 'Globus Pallidus internal', 'Subthalamic Nucleus', 'Red Nucleus', 'Substantia Nigra pars compacta ', 'Substantia Nigra pars reticulata'};
incl_aseg_rois = {'Cerebellum-White-Matter', 'Cerebellum-Cortex', 'Thalamus-Proper', 'Caudate', 'Putamen', 'Pallidum', 'Hippocampus', 'Amygdala', 'Accumbens-area', 'VentralDC'};

%% Create stat matrix
volmat = [];
for voli = 1:size(vols,4)
    disp(sprintf('Iteration %d of %d',voli,size(vols,4))); 
    [vol_roi_sum(voli), vol_roi_vox(voli)] = vol_by_roi(vols(:,:,:,voli),'atlas','ABCD3');
    volmat(:,voli) = vol_roi_vox(voli).vox;
end
groups = vol_roi_vox(1).cats;

%% Run function
% aseg is in ABCD3, thalamus and pauli are in ABCD2
plot_rois(volmat, agevals_tbl/12, incl_aseg_rois, fstem_imaging); 
% plot_rois(volmat, agevals_tbl/12, incl_thalamus_rois, fstem_imaging);
% plot_rois(volmat, agevals_tbl/12, incl_pauli_rois, fstem_imaging);


%% Old code
if 0
    for voli = 1:size(vols,4)
        [thalamus_sum, thalamus_vox] = extract_roi_val(vols(:,:,:,agei), thalamus_names, 0.7, incl_thalamus_rois, prob, roi_sum, roi_vox);
        [pauli_sum, pauli_vox] = extract_roi_val(vols(:,:,:,agei), pauli_names, 0.5,  incl_pauli_rois, pauli_prob, roi_sum, roi_vox);
        
        categories = thalamus_sum.cats;
        thalamus_mean = []; thalamus_ci = []; thalamus_std = [];
        for c = 1:length(categories)
            cat_ind = (ismember(thalamus_vox.cats, categories{c}) & thalamus_vox.vox~=0);
            thalamus_mean = [thalamus_mean, sum(thalamus_vox.vox(cat_ind))/sum(cat_ind)];
            thalamus_std = [thalamus_std, std(thalamus_vox.vox(cat_ind))];
        end
    end
end