% Run FEMA_cluster_gather on age project 
% Diana Smith
% Feb 2024
addpath(genpath('/home/d9smith/github/cmig_tools_internal/'));

dirname_out_stem = fullfile('/space/syn50/1/data/ABCD/d9smith/age/results_1000perms_cluster_2024-02-22');
designmat = dir(sprintf('%s/designMat*', dirname_out_stem));
designmat = {designmat.name}';


dirname_out = dirname_out = {fullfile(outDir)}; % filepath to save FEMA output
nperms = 1000;
% modality = {'RNI' 'RNT' 'RND' 'HNT'};
modality = {'RNT'};

for m=1:length(modality)
    fstem_imaging = modality{m};
    clusterjobs_outdir = fullfile(dirname_out_stem, designmat,'dt-voxel_img-RNT_njobs-100_nperms-10');
    outdir = clusterjobs_outdir;
    for d=[2 3] % d=1:length(clusterjobs_outdir)
        [zmat_perm, beta_hat_perm, tfce_perm, colnames_interest, save_params, mask] = FEMA_cluster_gather(outdir{d},clusterjobs_outdir{d},nperms, 'calc_perm_pvals',1)
    end
end

% Run FEMA_cluster_gather
[zmat_perm, beta_hat_perm, tfce_perm, colnames_interest, save_params, mask] = FEMA_cluster_gather_amd(dirname_imaging,dirname_out,nperms, 'calc_perm_pvals',1,'fstem_pheno',fstem_imaging)
[zmat_perm, beta_hat_perm, tfce_perm, colnames_interest, save_params, mask] = FEMA_cluster_gather(dirname_vertexwise_tmp,'/space/amdale/1/tmp/FEMA_cluster_wrapper_out',fstem_pheno,nperms,'calc_perm_pvals',0,'fstem_pheno',fstem_pheno); % Should add whitening option
 