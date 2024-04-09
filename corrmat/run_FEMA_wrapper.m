% Run FEMA wrapper for age and rsfmri corrmat analysis using ABCD 5.1 data
% Diana Smith
% Mar 2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify where to store results
dirname_out = fullfile('/space/syn50/1/data/ABCD/d9smith/age/corrmat/results_2024-03-20');

if ~exist(dirname_out, 'dir')
    mkdir(dirname_out)
end

% start diary
diary_path=strcat(dirname_out,'/','diary_', datestr(now, 'yyyy-mm-dd_HHMM'));
diary(diary_path);

% Add cmig_tools to path
addpath(genpath('~/github/cmig_tools_internal'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUTS TO FEMA_wrapper.m

abcd_sync_path=fullfile('/space/syn65/1/data/abcd-sync/5.0');
dirname_tabulated = fullfile(abcd_sync_path,'tabulated/released/core/imaging');

% path to designmat directory
designmat_dir = '/space/syn50/1/data/ABCD/d9smith/age/corrmat/designmat';
designmat_file = dir(sprintf('%s/designmat*.txt', designmat_dir));
designmat_file = {designmat_file.name}';
fname_design = strcat(designmat_dir, '/', designmat_file);

outdir_file = strrep(designmat_file, '.txt', '');
dirname_out = strcat(dirname_out,'/',outdir_file); 

% Optional inputs
ranknorm = 0; % Rank normalizes dependent variables (Y) (default = 0)
contrasts = [];
nperms = 0; % Number of permutations - if wanting to use resampling methods nperms>0
RandomEffects = {'F','S','E'}; % Random effects to include: family, subject, error
mediation = 0; % If wanting to use outputs for a mediation analysis set mediation=1 - ensures same resampling scheme used for each model in fname_design
PermType = 'wildbootstrap-nn'; %Default resampling method is null wild-bootstrap - to run mediation analysis need to use non-null wild-bootstrap ('wildboostrap-nn')
tfce = 0; % If wanting to run threshold free cluster enhancement (TFCE) set tfce=1 (default = 0)
colsinterest=[1]; % Only used if nperms>0. Indicates which IVs (columns of X) the permuted null distribution and TFCE statistics will be saved for (default 1, i.e. column 1)
output = 'mat';

datatype = 'corrmat';
modality = 'rsfmri';

% uses path structure in abcd-sync to automatically find data
dirname_imaging = fullfile(abcd_sync_path, '/imaging_concat/', datatype, modality); % filepath to imaging data

% modality = {'RNI' 'RNT' 'RND' 'RIF' 'RDF' 'HNT' 'HNI' 'HND' 'HIF' 'HDF' 'FNI' 'RD' 'RI' 'RT' 'HD' 'HI' 'HT'};
modality = {'rsfMRI_roi_corr'}; % could also try 'rsfMRI_roi_tseries'

for m=1:length(modality)
    fstem_imaging=modality{m};

    % Run FEMA
    [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm analysis_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
    'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest, 'output', output);
end 
