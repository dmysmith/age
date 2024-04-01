% Run FEMA wrapper for age analysis with volumetric data 
% Diana Smith
% Mar 2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify where to store results
dirname_out = fullfile('/space/cluster/1/ABCD/users/d9smith/age/volume/results_2024-04-01');

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
designmat_dir = '/space/cluster/1/ABCD/users/d9smith/age/volume/designmat';
designmat_file = dir(sprintf('%s/designmat*.txt', designmat_dir));
designmat_file = {designmat_file.name}';
fname_design = strcat(designmat_dir, '/', designmat_file);

outdir_file = strrep(designmat_file, '.txt', '');
dirname_out = strcat(dirname_out,'/',outdir_file); 

% Optional inputs
ranknorm = 1; % Rank normalizes dependent variables (Y) (default = 0)
contrasts = [];
nperms = 0; % Number of permutations - if wanting to use resampling methods nperms>0
RandomEffects = {'F','S','E'}; % Random effects to include: family, subject, error
mediation = 0; % If wanting to use outputs for a mediation analysis set mediation=1 - ensures same resampling scheme used for each model in fname_design
PermType = 'wildbootstrap-nn'; %Default resampling method is null wild-bootstrap - to run mediation analysis need to use non-null wild-bootstrap ('wildboostrap-nn')
tfce = 0; % If wanting to run threshold free cluster enhancement (TFCE) set tfce=1 (default = 0)
colsinterest=[1]; % Only used if nperms>0. Indicates which IVs (columns of X) the permuted null distribution and TFCE statistics will be saved for (default 1, i.e. column 1)
output = 'mat';

doVertex = 1;
doJA = 1;
do_external = 1;
do_manual_phenos = 1;


if doVertex
    datatype = 'vertex';
    modality='smri';

    % uses path structure in abcd-sync to automatically find data
    dirname_imaging = fullfile(abcd_sync_path, '/imaging_concat/vertexwise/', modality); % filepath to imaging data

    modality = {'area_ic5_sm1000' 'thickness_ic5_sm1000'};

    for m=1:length(modality)
        fstem_imaging=modality{m};

        % Run FEMA
        [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
        'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest, 'output', output);
    end 
end

if doJA
    datatype = 'voxel';
    modality='smri';

    % uses path structure in abcd-sync to automatically find data
    dirname_imaging = fullfile(abcd_sync_path, '/imaging_concat/voxelwise/', modality); % filepath to imaging data

    modality = {'JA'};

    for m=1:length(modality)
        fstem_imaging=modality{m};

        % Run FEMA
        [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
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
    modality = {'smr'};
  
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
          [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, imaging_path, datatype,...
          'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest, 'output', output);
        catch ME
          fprintf('Error running fema_wrapper using %s: %s', fstem_imaging, ME.message);
          continue; % Jump to next iteration
        end
      end
  
    end 
  
  end

  if do_manual_phenos
    datatype = 'external';
    imaging_dir = '/space/cluster/1/ABCD/users/d9smith/age/volume/phenofiles';
    dirname_imaging = dir(sprintf('%s/*.csv', imaging_dir));
    imaging_file = {dirname_imaging.name}';
    for i=1:length(imaging_file)
      imaging_path = strcat(imaging_dir,'/',imaging_file{i});
      fstem_imaging = strrep(imaging_file{i},'.csv','');
      try
        % Run FEMA
        [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, imaging_path, datatype,...
        'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest, 'output', output);
      catch ME
        fprintf('Error running fema_wrapper using %s: %s\n', fstem_imaging, ME.message);
        continue; % Jump to next iteration
      end
    end

  end