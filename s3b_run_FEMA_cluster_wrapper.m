%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run FEMA for age analysis using ABCD 5.0 data
%% Create jobs to submit to cluster
%% Diana Smith
%% Feb 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify where to store results
dirname_out_stem = fullfile('/space/syn50/1/data/ABCD/d9smith/age/results_1000perms_cluster_2024-02-22');
batchdir_stem = '/home/d9smith/batchdirs/age';

if ~exist(dirname_out_stem, 'dir')
  mkdir(dirname_out_stem)
end

% start diary
diary_path=strcat(dirname_out_stem,'/','diary_', datestr(now, 'yyyy-mm-dd_HHMM'));
diary(diary_path);

% Add cmig_tools to path
addpath(genpath('~/github/cmig_tools_internal'))

% Get path to local copy of abcd-sync using abcdConfig
abcd_sync_path=fullfile('/space/syn65/1/data/abcd-sync/5.0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SECTION 1 - MEDIATION DESIGN MATRICES 

% Inputs to FEMA_wrapper
dirname_tabulated = fullfile(abcd_sync_path,'tabulated/img'); % directory to tabulated imaging data on abcd-sync 
contrasts=[]; % Contrasts relate to columns in design matrix e.g. [1 -1] will take the difference between cols 1 and 2 in your design matrix (X).  This needs to be padded with zeros at the beginning but not the end.
ranknorm = 1; % Rank normalizes dependent variables (Y) (default = 0)
nperms = 10; % Number of permutations - if wanting to use resampling methods nperms>0
njobs = 100; % Number of permutations - if wanting to use resampling methods nperms>0
RandomEffects = {'F','S','E'}; % Random effects to include: family, subject, error
tfce = 0; % If wanting to run threshold free cluster enhancement (TFCE) set tfce=1 (default = 0)
colsinterest=[1:5]; % Only used if nperms>0. Indicates which IVs (columns of X) the permuted null distribution and TFCE statistics will be saved for (default 1, i.e. column 1)
datatype = 'voxel';
modality='dmri'; % concatenated imaging data stored in directories based on modality (smri, dmri, tfmri, roi)
dirname_imaging = fullfile(abcd_sync_path, '/imaging_concat/voxelwise/', modality); % filepath to imaging data

% Note that FA and MD are not in here
% modality = {'RNI' 'RNT' 'RND' 'RIF' 'RDF' 'HNT' 'HNI' 'HND' 'HIF' 'HDF' 'FNI' 'RD' 'RI' 'RT' 'HD' 'HI' 'HT'};
modality = {'RNI' 'RNT' 'RND' 'HNT'};

%To run multiple design matrices with same imaging data populate each row with path to each design matrix
designmat_dir = '/space/syn50/1/data/ABCD/d9smith/age/designMat';
designmat_file = dir(sprintf('%s/designMat*.txt', designmat_dir));
designmat_file = {designmat_file.name}';

%% SECTION 1 - MEDIATION ANALYSES
mediation = 1; % If wanting to use outputs for a mediation analysis set mediation=1 - ensures same resampling scheme used for each model in fname_design
PermType = 'wildbootstrap-nn';

% select only nested design matrices
idx = find(contains(designmat_file,'designMat5'));
fname_design = strcat(designmat_dir, '/', designmat_file(idx));

outdir_file = strrep(designmat_file(idx), '.txt', '');
dirname_out = strcat(dirname_out_stem,'/',outdir_file); 

% define pairs of design matrices
stems = {'all' 'f' 'm'};
idx_all = find(contains(fname_design, 'designMat5_'));
idx_f = find(contains(fname_design, 'designMat5f_'));
idx_m = find(contains(fname_design, 'designMat5m_'));

designmat_pair_list = {fname_design(idx_all);fname_design(idx_f);fname_design(idx_m)};
outdir_pair_list = {dirname_out(idx_all);dirname_out(idx_f);dirname_out(idx_m)}; 

for m=1:length(modality)
  fstem_imaging=modality{m};

  % One call for each pair of nested design matrices
  for p=1:length(designmat_pair_list)
    batchdir_full = sprintf('%s_%s_%s_%s',batchdir_stem, 'mediation', stems{p}, fstem_imaging); 
    % Run FEMA
    FEMA_cluster_wrapper(fstem_imaging, designmat_pair_list{p}, outdir_pair_list{p}, dirname_tabulated, dirname_imaging, datatype, batchdir_full, njobs,...
      'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest);
  end
end 

%% SECTION 2 - NON MEDIATION
mediation = 0; % If wanting to use outputs for a mediation analysis set mediation=1 - ensures same resampling scheme used for each model in fname_design
PermType = 'wildbootstrap';

% find design matrices that are not part of nested pairs
idx = find(~contains(designmat_file,'designMat5'));
fname_design = strcat(designmat_dir, '/', designmat_file(idx));

outdir_file = strrep(designmat_file(idx), '.txt', '');
dirname_out = strcat(dirname_out_stem,'/',outdir_file); 

for m=1:length(modality)
  fstem_imaging=modality{m};
  batchdir_full = sprintf('%s_%s',batchdir_stem, fstem_imaging);
  % Run FEMA
  FEMA_cluster_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype, batchdir_full, njobs,...
    'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest);
end