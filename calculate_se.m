% calculate confidence intervals for spline functions
designmat = 'designMat5_red_BFsSexIncEducHispPCsScanSoftMotion_bly2y4';
load(sprintf('/space/syn50/1/data/ABCD/d9smith/age/results_1000perms_cluster_2024-02-22/%s/dt-voxel_img-RNT_njobs-100_nperms-10/nonnullWB_1000perms/FEMA_wrapper_output_voxel_RNT_splinevol_perm_1.mat', designmat));

% calculate standard error
vol_se = std(vol_perm, [], 5) / sqrt(size(vol_perm, 5));

output_file = sprintf('/space/syn50/1/data/ABCD/d9smith/age/results_1000perms_cluster_2024-02-22/%s/dt-voxel_img-RNT_njobs-100_nperms-10/nonnullWB_1000perms/FEMA_wrapper_output_voxel_RNT_splinevol_se_100perms.mat', designmat);
save(output_file, 'vol_se');
