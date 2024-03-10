% Visualize results of age analysis
% Diana Smith
% March 2024

% path to designmat directory
designmat_dir = '/space/syn50/1/data/ABCD/d9smith/age/mediation/designmat';
designmat_file = dir(sprintf('%s/designmat*.txt', designmat_dir));
designmat_file = {designmat_file.name}';
fname_design = strcat(designmat_dir, '/', designmat_file);

% path to FEMA results directory
dirname_out = fullfile('/space/syn50/1/data/ABCD/d9smith/age/mediation/results_2024-03-07');
outdir_file = strrep(designmat_file, '.txt', '');
dirname_out = strcat(dirname_out,'/',outdir_file); 

% Add cmig_tools to path
addpath(genpath('~/github/cmig_tools_internal'))

modality = {'RNT' 'RNI' 'RND' 'HNT'};