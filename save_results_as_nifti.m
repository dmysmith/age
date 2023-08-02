%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save FEMA output in NIfTI format 
%%
%% Diana Smith
%% August 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ADD all ABCD CMIG tools directories to MATLAB path
addpath(genpath('/home/d9smith/github/cmig_tools_internal'));

% Specify where FEMA_wrapper output is saved and load into MATLAB workspace
results_folder='/space/syn50/1/data/ABCD/d9smith/age/results_2023-07-27_completecases';
designmat_folder='designMat2_AgeSexIncEducHispPCsScanSoftMotion_y2y4'; 
% designMat1_AgeSexIncEducHispPCsScanSoftMotion_bly2 
% designMat2_AgeSexIncEducHispPCsScanSoftMotion_y2y4
% designMat3_AgeSexIncEducHispPCsScanSoftMotion_bly2y4

dirname_out = sprintf('%s/%s',results_folder,designmat_folder);

modality = 'dmri';
fstem_imaging='RNI'; %imaging phenotype used for analysis e.g. 'RNI' 'MD' 'JA'

fname_results = sprintf('%s/FEMA_wrapper_output_voxel_%s.mat',dirname_out,fstem_imaging);
% load(fname_results,'vol_beta_hat','vol_z','colnames_model'); % load FEMA output - only need some variables
load(fname_results); % load FEMA output

% NIfTI file will be saved to same location as original results
% fname_nifti = sprintf('%s/FEMA_wrapper_output_voxel_%s.nifti',dirname_out,fstem_imaging);

% Various inputs
dirABCD         = '/space/syn65/1/data/abcd-sync/';
dataRelease     = '5.0';
atlasVersion    = 'ABCD3_cor10';
% dirCode         = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/code/cmig_tools_internal-beta';
% dirOut          = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/FA_voxelWise_FSE';
% dirSupport      = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/support_files/ABCD_rel4.0_unfiltered/';
% dirDesign       = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD';
dirTabulated    = fullfile(dirABCD, dataRelease, 'tabulated', 'released', 'core', 'imaging'); 
dirImaging      = fullfile(dirABCD, dataRelease, 'imaging_concat', 'voxelwise', atlasVersion, modality);
fnameGRM        = fullfile('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/genomics/ABCD_rel4.0_grm.mat'); % still using from 4.0 files
% fnameDesign     = 'GlobalDesignMatrix.csv';

%% Get voxel-wise data
[ymat, iid_concat, eid_concat, ivec_mask, mask, colnames_imaging, pihat] = FEMA_process_data(fstem_imaging, dirTabulated, dirImaging, 'voxel', 'pihat_file', fnameGRM);

%% Saving - code from Pravesh / FEMA_wrapper 
if sum(~mask)>0
    z_tmp        = zeros(size(zmat,1),size(mask,2));
    p_tmp        = zeros(size(logpmat,1),size(mask,2));
    beta_tmp     = zeros(size(beta_hat,1),size(mask,2));
    betase_tmp   = zeros(size(beta_se,1),size(mask,2));
    sig2mat_tmp  = zeros(size(sig2mat,1),size(mask,2));
    sig2tvec_tmp = zeros(size(sig2tvec,1),size(mask,2));
    
    z_tmp(:,ivec_mask)       = zmat;
    p_tmp(:,ivec_mask)       = logpmat;
    beta_tmp(:,ivec_mask)    = beta_hat;
    betase_tmp(:,ivec_mask)  = beta_se;
    sig2mat_tmp(:,ivec_mask) = sig2mat;
    sig2tvec_tmp(:,ivec_mask)= sig2tvec;
    
    zmat     = z_tmp;
    logpmat  = p_tmp;
    beta_hat = beta_tmp;
    beta_se  = betase_tmp;
    sig2mat  = sig2mat_tmp;
    sig2tvec = sig2tvec_tmp;
    
    if 0 % originally: if nperms > 0
        zperm_tmp           = zeros(size(zmat_perm,1),size(mask,2),size(zmat_perm,3));
        betaperm_tmp        = zeros(size(beta_hat_perm,1),size(mask,2),size(beta_hat_perm,3));
        betaseperm_tmp      = zeros(size(beta_se_perm,1),size(mask,2),size(beta_se_perm,3));
        sig2matperm_tmp     = zeros(size(sig2mat_perm,1),size(mask,2),size(sig2mat_perm,3));
        sig2tvecperm_tmp    = zeros(size(sig2tvec_perm,1),size(mask,2),size(sig2tvec_perm,3));
        
        zperm_tmp(:,ivec_mask,:)        = zmat_perm;
        betaperm_tmp(:,ivec_mask,:)     = beta_hat_perm;
        betaseperm_tmp(:,ivec_mask,:)   = beta_se_perm;
        sig2matperm_tmp(:,ivec_mask,:)  = sig2mat_perm;
        sig2tvecperm_tmp(:,ivec_mask,:) = sig2tvec_perm;
        
        zmat_perm       = zperm_tmp;
        beta_hat_perm   = betaperm_tmp;
        beta_se_perm    = betaseperm_tmp;
        sig2mat_perm    = sig2matperm_tmp;
        sig2tvec_perm   = sig2tvecperm_tmp; 
    end
end

%% Additional content for writing NIfTI images
vol_z        = zeros([size(mask) size(zmat,1)]);
vol_logp     = zeros([size(mask) size(zmat,1)]);
vol_beta_hat = zeros([size(mask) size(zmat,1)]);
vol_beta_se  = zeros([size(mask) size(zmat,1)]);

for j = 1:size(zmat,1)
    vol_z(:,:,:,j)          = single(fullvol(zmat(j,:),mask));
    vol_logp(:,:,:,j)       = single(fullvol(logpmat(j,:),mask));
    vol_beta_hat(:,:,:,j)   = single(fullvol(beta_hat(j,:),mask));
    vol_beta_se(:,:,:,j)    = single(fullvol(beta_se(j,:),mask));
end

vol_sig2t            = zeros([size(mask) 1]);
vol_sig2t(ivec_mask) = single(sig2tvec);
vol_sig2             = zeros([size(mask) size(sig2mat,1)]);

for j = 1:size(sig2mat,1)
    vol_sig2(:,:,:,j) = single(fullvol(sig2mat(j,:),mask));
end

% ============================================================================================================================
% == NIFTI Output == FIXME: no longer used for DEAP
results = struct('beta_hat',vol_beta_hat,'beta_se',vol_beta_se,'zmat',vol_z,'logpmat',vol_logp,'sig2tvec',vol_sig2t,'sig2mat',vol_sig2);
writeNIFTI(results, sprintf('%s/%s',dirname_out,'nifti'), fstem_imaging, [], colnames_model);