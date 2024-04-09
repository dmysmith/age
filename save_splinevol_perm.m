%% calculating spline function over 1000 perms
addpath(genpath('/home/d9smith/github/cmig_tools_internal'));

% load results file
% designmat = 'designMat5_red_BFsSexIncEducHispPCsScanSoftMotion_bly2y4';
designmat = 'designMat5_BFsSexIncEducHispPCsScanSoftMotionPDS_bly2y4'; 

results_file = sprintf('/space/syn50/1/data/ABCD/d9smith/age/results_1000perms_cluster_2024-02-22/%s/dt-voxel_img-RNT_njobs-100_nperms-10/nonnullWB_1000perms/FEMA_wrapper_output_voxel_RNT.mat', designmat);
logging('Loading results file from: %s', results_file);
load(results_file, 'colnames_model', 'beta_hat_perm', 'mask');

% specify output file
output_file = sprintf('/space/syn50/1/data/ABCD/d9smith/age/results_1000perms_cluster_2024-02-22/%s/dt-voxel_img-RNT_njobs-100_nperms-10/nonnullWB_1000perms/FEMA_wrapper_output_voxel_RNT_splinevol_perm', designmat);

% load basis functions
tbl_bf = readtable('/space/syn50/1/data/ABCD/d9smith/age/basis.txt');
agevals_tbl = linspace(100,200,101);

%% Combine basis functions to calculate spline function values
jvec_bf = find(find_cell(regexp(colnames_model,'^bf_','match')));

% for perm = 1:size(beta_hat_perm,3)
for perm = 1:101 % shortcut: running for first 100 perms only because otherwise it takes too much memory
    logging(sprintf('Calculating spline function - perm %0d of %0d',perm,size(beta_hat_perm,3)));
    % vol_perm(:,:,:,:,perm) = single(FEMA_convert_splinevols(tbl_bf, agevals_tbl, beta_hat_perm(jvec_bf,:,perm), mask));
    [dummy beta_hat_spline(:,:,perm)] = FEMA_convert_splinevols(tbl_bf, agevals_tbl, beta_hat_perm(jvec_bf,:,perm), mask);
    % every 100 iterations, save variable and clear
    if mod(perm-1, 100) == 0
        % keyboard;
        beta_hat_spline = single(beta_hat_spline);
        filename = sprintf('%s_%s_%0d.mat',output_file, 'beta_hat', (perm-1)/100);
        save(filename, 'beta_hat_spline', '-v7.3');
        logging('Results saved to %s', filename);
        clear beta_hat_spline;
    end
end

