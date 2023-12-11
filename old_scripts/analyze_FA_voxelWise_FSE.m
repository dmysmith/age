%% Analyze dMRI FA measures vertex wise
%% Setup
% Add FEMA code to path
addpath('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/code/cmig_tools_internal-beta/cmig_tools_utils/matlab');
addpath('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/code/cmig_tools_internal-beta/FEMA');

% Specify release number
dataRelease = '4.0';

% Specify atlas version
atlasVersion = 'ABCD2_cor10';

% Various paths
dirABCD         = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/';
dirCode         = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/code/cmig_tools_internal-beta';
dirOut          = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/FA_voxelWise_FSE';
dirSupport      = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/support_files/ABCD_rel4.0_unfiltered/';
dirDesign       = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD';
dirTabulated    = fullfile(dirABCD, dataRelease, 'tabulated', 'released'); 
dirImaging      = fullfile(dirABCD, dataRelease, 'imaging_concat', 'voxelwise', atlasVersion, 'dmri');
fnameGRM        = fullfile(dirABCD, dataRelease, 'genomics', ['ABCD_rel', dataRelease, '_grm.mat']);
fnameDesign     = 'GlobalDesignMatrix.csv';
fnameImaging    = 'FA';

% Prepare output directory
if ~exist(dirOut, 'dir')
    mkdir(dirOut);
end

%% Read various files
% First, read the QC variables
dataQC = readtable(fullfile(dirABCD, dataRelease, 'tabulated', 'img', 'abcd_imgincl01.csv'));

% Next, get some additional info about MRI acquisition
dataInfo = readtable(fullfile(dirABCD, dataRelease, 'tabulated', 'img', 'abcd_mri01.csv'));

% Read various support files - only need allPCs and imgVars
imgVars       = readtable(fullfile(dirSupport, 'ABCD_rel4.0_covars_img_base_2yr.txt'));
allPCs        = readtable(fullfile(dirSupport, 'ABCD_rel4.0_pcs_base_2yr.txt'));

% Unique event names - useful later
allUqEvents = unique(dataQC.eventname);

%% Next, read the design matrix
designMatrix    = readtable(fullfile(dirDesign, fnameDesign));

%% Get voxel-wise data
[ymat, iid_concat, eid_concat, ivec_mask, mask, colnames_imaging, pihat] = FEMA_process_data(fnameImaging, dirTabulated, dirImaging, 'voxel', 'pihat_file', fnameGRM);

%% Within the design matrix, remove subjects who did not pass QC
toDelete = unique(designMatrix.src_subject_id(designMatrix.imgincl_dmri_include == 0));
designMatrix(ismember(designMatrix.src_subject_id, toDelete), :) = [];

%% Attempt to create covariates
% To use: age at recruitment, age delta, sex, scanner, household income, 
% highest parental education, and the first 20 genetic PCs

% Location for PCs
locGenesisPCs = ~cellfun(@isempty, regexpi(designMatrix.Properties.VariableNames, '^genesis'));

intercept       = ones(height(designMatrix), 1);
geneticPCs      = designMatrix{:, locGenesisPCs};

% Compute ageRecruitment - same as baseline age
locBaseline     = strcmpi(designMatrix.eventname, 'baseline_year_1_arm_1');
ageRecruitment  = [designMatrix.interview_age(locBaseline); designMatrix.interview_age(locBaseline)];

% Compute ageDelta
ageDelta                 = zeros(height(designMatrix), 1);
ageDelta(~locBaseline,1) = designMatrix.interview_age(~locBaseline) - designMatrix.interview_age(locBaseline);

% Sex
sex                  = categorical(designMatrix.sex);
namesSex             = categories(sex);
varSex               = dummyvar(sex);

% MRI device number
deviceInfo           = categorical(designMatrix.mri_info_deviceserialnumber);
namesDeviceInfo      = categories(deviceInfo);
varDeviceInfo        = dummyvar(deviceInfo);

% MRI software number
softwareInfo         = categorical(designMatrix.mri_info_softwareversion);
namesSoftwareInfo    = categories(softwareInfo);
varSoftwareInfo      = dummyvar(softwareInfo);

% Household income level
householdIncome      = categorical(designMatrix.household_income);
namesHouseholdIncome = categories(householdIncome);
varHouseholdIncome   = dummyvar(householdIncome);

% Parental educational level
parentalEdu             = categorical(designMatrix.high_educ);
namesParentalEducation  = categories(parentalEdu);
varParentalEducation    = dummyvar(parentalEdu);

% Put covariates together and put covariate names together
% Only keep the first level of sex
% Keep n-2 level of device info
% Keep n-2 level of software info
% Keep n-1 level of household income
% Keep n-1 level of parental education level
% Keeping n-1 levels of device info and software info leads to rank
% deficient matrix - why? Likely because there are too few observations
covariates = [intercept, ageRecruitment, ageDelta, varSex(:,1),         ...
             varDeviceInfo(:, 1:end-2), varSoftwareInfo(:, 1:end-2),    ...
             varHouseholdIncome(:, 1:end-1), varParentalEducation(:, 1:end-1), geneticPCs];
         
%% Ensure no rank deficiency
if rank(covariates) ~= size(covariates,2)
    error('Rank deficient covariates matrix; check covariates');
end
         
% Names of covariates
names_genesis = cellfun(@(x) strrep(x, ' ', ''), strcat({'genensis_PC'}, num2str((1:20)')), 'UniformOutput', false);
covarNames    = [{'Intercept', 'ageBaseline', 'ageDelta'}, namesSex(1)', namesDeviceInfo(1:end-2)', ...
                 namesSoftwareInfo(1:end-2)', namesHouseholdIncome(1:end-1)', namesParentalEducation(1:end-1)', names_genesis'];
             
%% Create a new version of the design matrix with covariates
% First four columns should be: src_subject_id, eventname, rel_family_id, age
toIntersect = cell2table([designMatrix.src_subject_id, designMatrix.eventname, num2cell([designMatrix.rel_family_id, designMatrix.interview_age, covariates])], 'VariableNames', ...
                         [{'src_subject_id', 'eventname', 'rel_family_id', 'age'}, covarNames]);

% Save as T1w design matrix
writetable(toIntersect, fullfile(dirOut, 'DesignMatrix_dMRI_FA.csv'));

%% Intersect the design matrix and the vertex-wise data
[X, iid, eid, fid, agevec, ymatUse, contrasts, colnames_model, pihatmat, PregID, HomeID] = FEMA_intersect_design(fullfile(dirOut, 'DesignMatrix_dMRI_FA.csv'), ...
                                                                                                              ymat, iid_concat, eid_concat, 'pihat', pihat);

%% Standardize ymatUse
ymatUse_std = (ymatUse - mean(ymatUse))./std(ymatUse);

%% Prepare random effects
SubjectEffect = iid;
FamilyEffect  = cellfun(@(x) strrep(x, ' ', ''), strcat({'F'}, num2str(fid)), 'UniformOutput', false);

%% Additional variables for FEMA_fit
nperms          = 0;
RandomEffects   = {'F', 'S', 'E'};
niter           = 1;
nbins           = 20;

%% Save variables
save(fullfile(dirOut, 'vars_analysis.mat'), 'SubjectEffect', 'FamilyEffect', 'X', 'iid', 'fid', 'agevec', 'ymatUse', 'ymatUse_std', 'colnames_model', '-v7.3');

%% Call FEMA_fit
initFEMA = tic;
[beta_hat, beta_se, zmat, logpmat, sig2tvec, sig2mat, Hessmat, logLikvec] =   ...
 FEMA_fit(X, SubjectEffect, eid, FamilyEffect, agevec, ymatUse_std, niter,    ...
          contrasts, nbins, [], 'RandomEffects', RandomEffects, 'nperms', nperms);
elapsedFEMA = toc(initFEMA);

%% Save everything
if not(exist(dirOut, 'dir'))
    mkdir(dirOut);
end
save(fullfile(dirOut, 'Results.mat'), '-v7.3');

%% Additional saving borrwed from FEMA_wrapper
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
    
    if nperms > 0
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

% Save output
base_variables_to_save = {'X','iid','eid','colnames_model','contrasts','zmat','logpmat','beta_hat','beta_se','sig2mat','sig2tvec','mask'};
save(fullfile(dirOut, 'Results_masked.mat'), base_variables_to_save{:}, '-v7.3');

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
writeNIFTI(results, dirOut, fnameImaging, [], colnames_model);

%% Super additional content - maybe useful for showVol
save(fullfile(dirOut, 'Results_maskedVolType.mat'), 'vol_z','vol_beta_hat','logpmat','vol_sig2','vol_sig2t','-v7.3');