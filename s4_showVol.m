%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run showVol for age analysis using ABCD 5.0 data
%%
%% Diana Smith
%% July 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: This analysis is an extension of that found in Palmer 2022, 
% Microstructural development from 9 to 14 years: Evidence from the ABCD Study
% doi: https://doi.org/10.1016/j.dcn.2021.101044
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify location to save images
plots_folder = '/home/d9smith/projects/age/plots';

% ADD all ABCD CMIG tools directories to MATLAB path
addpath(genpath('/home/d9smith/github/cmig_tools_internal'));

% Specify where FEMA_wrapper output is saved and load into MATLAB workspace

results_folder='/space/syn50/1/data/ABCD/d9smith/age/results_2023-07-23';
designmat_folder='designMat2_AgeSexIncEducHispPCsScanSoftMotion_y2y4'; 
% designMat1_AgeSexIncEducHispPCsScanSoftMotion_bly2 
% designMat2_AgeSexIncEducHispPCsScanSoftMotion_y2y4
% designMat3_AgeSexIncEducHispPCsScanSoftMotion_bly2y4

dirname_out = sprintf('%s/%s',results_folder,designmat_folder);

fstem_imaging='RNI'; %imaging phenotype used for analysis

fname_results = sprintf('%s/FEMA_wrapper_output_voxel_%s.mat',dirname_out,fstem_imaging);
% load(fname_results,'vol_beta_hat','vol_z','colnames_model'); % load FEMA output - only need some variables
load(fname_results); % load FEMA output

% Specify data release and atlas version
dataRelease='5.0';
atlasVersion = 'ABCD2_cor10'; % TODO - update to ABCD3_cor10 ?
  
loadPrerenderedIfNeeded(atlasVersion) % loads atlas specific data
global PRI % saves atlas data as a global variable so does not need to be loaded every time you run showVol

vol_stat = vol_z; %specify output to plot e.g. vol_z, vol_beta_hat etc
stat_name = 'z stat'; %name to label output

nvox=156662;
thresh_bonfcorr=abs(norminv(0.05/nvox));

% 3) Use a helper function to convert FEMA output into volumes for showVol. 

%  for each IV, creates a raw grayscale map and colorized map overlaid on T1 
limits = [-40 40];  % Limits for colormap: (1) min (2) max (3) threshold (leave third value empty to have unthresholded maps)
% limits = [-.01 .01];
index = [1 13];                          % Indicates which IVs to plot: colnames_model(index)
showPval = false;                   % If showPval==true, creates a map of log10 p-values
cmap = redblackblue;                % Colormap
bg_vol = [];                        % Set background image to overlay statistics. Default: atlas T1
CSFprob = [];                       % Probabilistic threshold for masking CSF e.g. if CSFprob=0.8, only voxels in which 80% of participants labelled that voxel as CSF will be masked (Default: no masking)
interpMethod = [];                  % By default uses linear interpolation to go from 2mm to 1mm voxels. Specify 'nearest' to show actual resolution

vols = convertFEMAVols(vol_stat, fstem_imaging, stat_name, colnames_model, index, limits, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion)

% 4) Visualize using showVol
coords = [-7 -13 -4]; % Opens showVol to a specific location e.g. R NAcc
showVol(vols, PRI.ABCD2.T1, PRI.ABCD2.CO, PRI.ABCD2.FOD, struct('roiatlas','ABCD2'), coords)

% create directory (if it doesn't exist) and navigate there
savepath = sprintf('%s/%s',plots_folder, designmat_folder);
if ~exist(savepath,'dir'), mkdir(savepath); end
cd(savepath);

% While running, use these keyboard shortcuts to save images
% -- 'v': cycle to the next volume
% -- 'o': cycle orientations
% -- '!': save screenshot of main axis
