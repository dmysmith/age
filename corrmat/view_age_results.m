% view results from FEMA analysis of age associations with corrmat
% Diana Smith
% March 2024

results_dir = '/space/syn50/1/data/ABCD/d9smith/age/corrmat/results_2024-03-20';

designmat = 'designmat1_AgeSexIncEducPCsScanSoftMotion';
% designmat1_AgeSexIncEducPCsScanSoftMotion
% designmat4_SAgeSexIncEducPCsScanSoftMotion

fstem_imaging = 'rsfMRI_roi_corr';

resultsfile = sprintf('%s/%s/FEMA_wrapper_output_corrmat_%s.mat', results_dir, designmat, fstem_imaging);

load(resultsfile);

index = 1; % index of variable of interest, e.g., colnames_model{index}
dims = size(beta_hat(index,:));
corrmat = reshape(beta_hat(index,:), [dims(1) prod(dims(2:end))]); % BROKEN
figure; imagesc(beta_hat(index,1:582:end));