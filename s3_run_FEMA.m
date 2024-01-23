%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run FEMA for age analysis using ABCD 5.0 data
%%
%% Diana Smith
%% July 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Note: This analysis is an extension of that found in Palmer 2022, 
%% Microstructural development from 9 to 14 years: Evidence from the ABCD Study
%% doi: https://doi.org/10.1016/j.dcn.2021.101044

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify where to store results
dirname_out = fullfile('/space/syn50/1/data/ABCD/d9smith/age/results_2024-01-22');

if ~exist(dirname_out, 'dir')
      mkdir(dirname_out)
end

% start diary
diary_path=strcat(dirname_out,'/','diary_', datestr(now, 'yyyy-mm-dd_HHMM'));
diary(diary_path);

% Add cmig_tools to path
addpath(genpath('~/github/cmig_tools_internal'))

% Get path to local copy of abcd-sync using abcdConfig
abcd_sync_path=fullfile('/space/syn65/1/data/abcd-sync/5.0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUTS TO FEMA_wrapper.m

dirname_tabulated = fullfile(abcd_sync_path,'tabulated/img'); % directory to tabulated imaging data on abcd-sync 

%To run multiple design matrices with same imaging data populate each row with path to each design matrix
designmat_dir = '/space/syn50/1/data/ABCD/d9smith/age/results_2024-01-22';
designmat_file = dir(sprintf('%s/designMat*.txt', designmat_dir));
designmat_file = {designmat_file.name}';
fname_design = strcat(designmat_dir, '/', designmat_file);

outdir_file = strrep(designmat_file, '.txt', '');
dirname_out = strcat(dirname_out,'/',outdir_file); 

% for running just one design matrix
% fname_design = '/space/syn50/1/data/ABCD/d9smith/age/results_2023-07-27_completecases/designMat4_BFsSexIncEducHispPCsScanSoftMotion_bly2y4.txt'
% dirname_out = '/space/syn50/1/data/ABCD/d9smith/age/results_2023-07-27_completecases/designMat4_BFsSexIncEducHispPCsScanSoftMotion_bly2y4';

% Note: using GRM from 4.0 data release
% fname_pihat = fullfile('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/genomics/ABCD_rel4.0_grm.mat'); 

% Optional inputs for `FEMA_wrapper.m` depending on analysis
contrasts=[]; % Contrasts relate to columns in design matrix e.g. [1 -1] will take the difference between cols 1 and 2 in your design matrix (X).  This needs to be padded with zeros at the beginning but not the end.
ranknorm = 1; % Rank normalizes dependent variables (Y) (default = 0)
nperms = 0; % Number of permutations - if wanting to use resampling methods nperms>0
RandomEffects = {'F','S','E'}; % Random effects to include: family, subject, error
mediation = 0; % If wanting to use outputs for a mediation analysis set mediation=1 - ensures same resampling scheme used for each model in fname_design
PermType = 'wildbootstrap'; %Default resampling method is null wild-bootstrap - to run mediation analysis need to use non-null wild-bootstrap ('wildboostrap-nn')
tfce = 0; % If wanting to run threshold free cluster enhancement (TFCE) set tfce=1 (default = 0)
colsinterest=[1]; % Only used if nperms>0. Indicates which IVs (columns of X) the permuted null distribution and TFCE statistics will be saved for (default 1, i.e. column 1)

% toggle if you just want to do a subset of modalities
do_vertex = 1;
do_smri = 1;
do_dmri = 1;
do_external = 1;

% output = 'nifti'; % toggling output format - default is 'mat'
output = 'mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VERTEXWISE ANALYSES

if do_vertex
  datatype = 'vertex';
  modality='smri'; % concatenated imaging data stored in directories based on modality (smri, dmri, tfmri, roi)

  % uses path structure in abcd-sync to automatically find data
  dirname_imaging = fullfile(abcd_sync_path, '/imaging_concat/vertexwise/', modality); % filepath to imaging data
  modality = {'area_ic5_sm1000' 'thickness_ic5_sm1000' 'sulc_ic5_sm1000'};

  for m=1:length(modality)
    fstem_imaging=modality{m};

    % Run FEMA
    [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm analysis_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
    'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest, 'output', output);
  end 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VOXELWISE ANALYSES

datatype = 'voxel'; % imaging modality selected

%% sMRI
if do_smri
      modality='smri'; % concatenated imaging data stored in directories based on modality (smri, dmri, tfmri, roi)

      % uses path structure in abcd-sync to automatically find data
      dirname_imaging = fullfile(abcd_sync_path, '/imaging_concat/voxelwise/', modality); % filepath to imaging data

      % Note that FA and MD are not in here
      modality = {'FA' 'JA' 'MD' 'T2' 'b0' 'bm' 'di1vol1' 'di2vol1' 'di3vol1' 'gm_new' 'nu' 'wm'};
      modality = {'JA' 'MD' 'FA'}; % just these to start 

      for m=1:length(modality)
            fstem_imaging=modality{m};

            % Run FEMA
            [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm analysis_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
            'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest, 'output', output);
      end 
end

%% dMRI
if do_dmri
      modality='dmri'; % concatenated imaging data stored in directories based on modality (smri, dmri, tfmri, roi)

      % uses path structure in abcd-sync to automatically find data
      dirname_imaging = fullfile(abcd_sync_path, '/imaging_concat/voxelwise/', modality); % filepath to imaging data

      % Note that FA and MD are not in here
      modality = {'RNI' 'RNT' 'RND' 'RIF' 'RDF' 'HNT' 'HNI' 'HND' 'HIF' 'HDF' 'FNI' 'RD' 'RI' 'RT' 'HD' 'HI' 'HT'};

      for m=1:length(modality)
            fstem_imaging=modality{m};

            % Run FEMA
            [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm analysis_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
            'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest, 'output', output);
      end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TABULATED IMAGING DATA ANALYSES
if do_external
  datatype = 'external';
  % modality='smri'; % concatenated imaging data stored in directories based on modality (smri, dmri, tfmri, roi)

  % uses path structure in abcd-sync to automatically find data
  imaging_dir = fullfile(abcd_sync_path, '/tabulated/released/core/imaging'); % filepath to imaging data
  modality = {'dti' 'rsi' 'smr' 'tfmr'};

  % imaging_file = dir(sprintf('%s/mri_y_%s_*.csv', imaging_dir, modality{m}));
  % imaging_file = {imaging_file.name}';
  % fname_design = strcat(designmat_dir, '/', imaging_file);

  for m=1:length(modality)
    dirname_imaging = dir(sprintf('%s/mri_y_%s_*.csv', imaging_dir, modality{m})); 
    imaging_file = {dirname_imaging.name}';
    
    for i=1:length(imaging_file)
      imaging_path = strcat(imaging_dir,'/',imaging_file{i});
      fstem_imaging = strrep(imaging_file{i},'.csv','');
      try
        % Run FEMA
        [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm analysis_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, imaging_path, datatype,...
        'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest, 'output', output);
      catch ME
        fprintf('Error running fema_wrapper using %s: %s', fstem_imaging, ME.message);
        continue; % Jump to next iteration
      end
    end

  end 

end