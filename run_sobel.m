% Run Sobel test for puberty as mediator of age effects

FEMA_outfile_reduced='/space/syn50/1/data/ABCD/d9smith/age/results_1000perms_cluster_2024-02-22/designMat5_red_BFsSexIncEducHispPCsScanSoftMotion_bly2y4/dt-voxel_img-RNT_njobs-100_nperms-10/nonnullWB_1000perms/FEMA_wrapper_output_voxel_RNT_splinevol_perm_beta_hat_1.mat';
FEMA_outfile_full='/space/syn50/1/data/ABCD/d9smith/age/results_1000perms_cluster_2024-02-22/designMat5_BFsSexIncEducHispPCsScanSoftMotionPDS_bly2y4/dt-voxel_img-RNT_njobs-100_nperms-10/nonnullWB_1000perms/FEMA_wrapper_output_voxel_RNT_splinevol_perm_beta_hat_1.mat';

model_red=load(FEMA_outfile_reduced);
model_full=load(FEMA_outfile_full);

% Uses save_params from full model output (last model in fname_design)

% For mediation analysis, the user must specify the independent variable of interest that is present in both design matrices as `IVname`.
% Sobel's test is testing for a mediation effect on the association between the specified IV of interest (`IVname`) and the
% imaging phenotype (`fstem_imaging`). The mediator is the variable that is present in the full model, but absent from the reduced model.
IVname='interview_age';
            
alpha=0.05; %generates 95% confidence intervals

tau = permute(model_red.beta_hat_spline,[3,2,1]);
tau_prime = permute(model_full.beta_hat_spline,[3,2,1]);

%Calculate mediation effect and std from bootstrapped distribution
mediation_effect_perms = tau - tau_prime;
mediation_effect = mediation_effect_perms(:,:,1);
mediation_se = std(mediation_effect_perms(:,:,2:end)); % AMD: include permutations only

FEMA_output_file = '/space/syn50/1/data/ABCD/d9smith/age/results_1000perms_cluster_2024-02-22/designMat5_BFsSexIncEducHispPCsScanSoftMotionPDS_bly2y4/dt-voxel_img-RNT_njobs-100_nperms-10/nonnullWB_1000perms/FEMA_wrapper_output_voxel_RNT.mat';
load(FEMA_output_file);
fstem_imaging = 'RNT';


designmat_file = '/space/syn50/1/data/ABCD/d9smith/age/designMat/designMat5_BFsSexIncEducHispPCsScanSoftMotionPDS_bly2y4.txt';
tbl_design = readtable(designmat_file);
agevals_tbl = linspace(100,200,101);

limits = [-0.0008 0.0008];  % Limits for colormap: (1) min (2) max (3) threshold (leave third value empty to have unthresholded maps)
index = [1];                    % Indicates which IVs to plot: colnames_model(index)
showPval = false;               % If showPval==true, creates a map of log10 p-values
cmap = redblackblue;            % Colormap
bg_vol = [];                    % Set background image to overlay statistics. Default: atlas T1
CSFprob = [];                   % Probabilistic threshold for masking CSF e.g. if CSFprob=0.8, only voxels in which 80% of participants labelled that voxel as CSF will be masked (Default: no masking)
interpMethod = [];              % By default uses linear interpolation to go from 2mm to 1mm voxels. Specify 'nearest' to show actual resolution
atlasVersion = 'ABCD3';

iterations = 9:12:size(mediation_effect,1);
for agei = iterations
  statname = sprintf('mediation effect (age: %.1f), %s',agevals_tbl(agei)/12, FEMA_outfile_full);
  vol_stat(:,:,:,agei) = fullvol(mediation_effect(agei,:),mask);
  vols(:,:,agei) = convertFEMAVols(vol_stat(:,:,:,agei),fstem_imaging, statname, colnames_model, index, limits, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion);
end

showVol(vols(:,2,iterations), struct('roiatlas','ABCD3')) % color


% load standard error
se_file = '/space/syn50/1/data/ABCD/d9smith/age/results_1000perms_cluster_2024-02-22/designMat5_BFsSexIncEducHispPCsScanSoftMotionPDS_bly2y4/dt-voxel_img-RNT_njobs-100_nperms-10/nonnullWB_1000perms/FEMA_wrapper_output_voxel_RNT_splinevol_se_100perms.mat';
load(se_file);

%% DIANA START HERE - need to fix for new dimensions of test statistics
%Calculate confidence intervals using:
%1) normal approach
mediation_effect_CInorm(1,:) = mediation_effect-(abs(norminv(alpha/2)).*mediation_se);
mediation_effect_CInorm(2,:) = mediation_effect+(abs(norminv(alpha/2)).*mediation_se);
%2) percentile approach
mediation_effect_CIprctile = prctile(mediation_effect_perms, [alpha/2,1-(alpha/2)],1); 

%Calculate sobel T statistic for every voxel and boostrap
mediation_tstat = mediation_effect_perms./mediation_se; 

mediation_tstat = mediation_tstat - mean(mediation_tstat(2:end,:)); % Remove mean

statvec = colvec(mediation_tstat(2:end,:)).^2;
xvec =  linspace(0,100,10001);
hc = hist(statvec,xvec); chc = cumsum(hc)/sum(hc);  

%pd_fit = makedist('gamma','a',0.5,'b',2); % Hard-code Chi-Square distribution
distname = 'gamma';
%      distname = 'weibull';
statvec=double(statvec);
pd_fit = fitdist(statvec,distname); % Fit to specific distribution 
chc_fit = pd_fit.cdf(xvec);
%figure; plot(xvec,-log10(1-chc),xvec,-log10(1-chc_fit),'LineWidth',2); 

mediation_pvals = pd_fit.cdf(mediation_tstat.^2,'upper');

fields={'beta_hat_perm','zmat_perm'};
model_red = rmfield(model_red, fields);
model_full = rmfield(model_full, fields);

mediation_tstat = mediation_tstat(1,:);
mediation_pvals = mediation_pvals(1,:);