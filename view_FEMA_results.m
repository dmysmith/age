% Visualize results of age analysis
% Diana Smith
% March 2024
desmatnum = '7f'; % which designmat_file you want, e.g. '8m'
modality = 'JA'; % which modality to index from results_file, e.g. 'RNT'

% path to FEMA results directory
dirname_out = fullfile('/space/cluster/1/ABCD/users/d9smith/age/volume/results_2024-04-01');
designmat_file = dir(sprintf('%s/designmat*', dirname_out));
designmat_file = {designmat_file.name}';
dirname_results = strcat(dirname_out, '/', designmat_file);

% path to save figures
savedir = '/home/d9smith/projects/age/volume/plots';

% Add cmig_tools to path
addpath(genpath('~/github/cmig_tools_internal'));
addpath(genpath('/home/d9smith/projects/age'));

% inputs for showVol
loadPrerenderedIfNeeded('ABCD3'); % loads atlas specific data
global PRI;

limits = [-0.15 0.15];  % Limits for colormap
agelim = [-0.01 0.01];
sexlim = [-0.4 0.4];
pdslim = [-0.06 0.06];

showPval = false;               % If showPval==true, creates a map of log10 p-values
cmap = redblackblue;            % Colormap
bg_vol = [];                    % Set background image to overlay statistics. Default: atlas T1
CSFprob = [];                   % Probabilistic threshold for masking CSF e.g. if CSFprob=0.8, only voxels in which 80% of participants labelled that voxel as CSF will be masked (Default: no masking)
interpMethod = [];              % By default uses linear interpolation to go from 2mm to 1mm voxels. Specify 'nearest' to show actual resolution
atlasVersion = 'ABCD3';

% load basis functions for s(PDS)
tbl_bf = readtable('/space/cluster/1/ABCD/users/d9smith/age/basis_pds.txt');
tbl_bf_age = readtable('/space/cluster/1/ABCD/users/d9smith/age/basis.txt');
agevec = linspace(100,220,121);
pdsvec = linspace(0.5,5.5,11);

% loop through designmats
% for des=1:length(designmat_file)
% for des=8
for des = find(contains(designmat_file, sprintf('designmat%s_',desmatnum)))
    if isempty(des) error('Design matrix not found. Check designmat_file to see list of design matrices.'); end
    results_file = dir(sprintf('%s/FEMA_wrapper_output_*.mat', dirname_results{des}));
    results_file = {results_file.name}';
    fname_results = strcat(dirname_results{des}, '/', results_file);

    % loop through modalities
    % for m=1:length(results_file)
    % for m=7
    for m = find(contains(results_file, modality))
        if isempty(m) error('Modality not found. Check results_file to see list of modalities.'); end
        prefix = 'FEMA_wrapper_output_';
        suffix = '_' + wildcardPattern + '.mat';
        datatype = extractBetween(results_file{m}, prefix, suffix); datatype = datatype{:};
        modality = extractBetween(results_file{m}, strcat(datatype,'_'), '.mat'); modality = modality{:};
        load(fname_results{m});

        age_idx = find(strcmp(colnames_model, 'interview_age'));
        sex_idx = find(contains(colnames_model, 'sexM'));
        pds_idx = find(contains(colnames_model, 'pds_y'));
        Sage_idx = find(find_cell(regexp(colnames_model,'^bf_demean','match')));
        Spds_idx = find(find_cell(regexp(colnames_model,'^bf_pds','match')));

        index = [age_idx pds_idx sex_idx]; % Indicates which IVs to plot: colnames_model(index)
        
        if strcmp(datatype,'voxel')
            clear vols;
            statname = sprintf('beta_hat, %s', designmat_file{des});

            % convert vol stats - age
            if ~isempty(age_idx)
                vols_tmp = convertFEMAVols(vol_beta_hat, modality, statname, colnames_model, age_idx, agelim, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion);
                if ~exist('vols') vols = vols_tmp; else vols = cat(2, vols, vols_tmp(1), vols_tmp(2)); end
            end

            % convert vol stats - sex
            if ~isempty(sex_idx)
                vols_tmp = convertFEMAVols(vol_beta_hat, modality, statname, colnames_model, sex_idx, sexlim, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion);
                if ~exist('vols') vols = vols_tmp; else vols = cat(2, vols, vols_tmp(2)); end
            end

            % convert vol stats - pds
            if ~isempty(pds_idx)
                vols_tmp = convertFEMAVols(vol_beta_hat, modality, statname, colnames_model, pds_idx, pdslim, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion);
                if ~exist('vols') vols = vols_tmp; else vols = cat(2, vols, vols_tmp(1), vols_tmp(2)); end
            end

            % calculate spline values
            if any(contains(colnames_model, 'bf_'))
                % pds spline values
                if ~isempty(Spds_idx)
                    [vols_tmp beta_tmp] = FEMA_convert_splinevols(tbl_bf, pdsvec, beta_hat([Spds_idx],:), mask);
                    for slice=1:size(vols_tmp,4)
                        statname = sprintf('beta_hat (PDS = %.1f), %s', pdsvec(slice), designmat_file{des});
                        test = convertFEMAVols(vols_tmp(:,:,:,slice), modality, statname, colnames_model(Spds_idx), 1, limits, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion);
                        if ~exist('vols') vols = test(2); else vols = cat(2, vols, test(2)); end
                    end
                end

                % calculate interaction spline values - beta_hat
                jvec_bf = find(find_cell(regexp(colnames_model,'interview_age_bf_','match')));
                if ~isempty(jvec_bf)
                    [vols_tmp beta_tmp] = FEMA_convert_splinevols(tbl_bf, pdsvec, beta_hat([jvec_bf],:), mask);
                    for slice=1:size(vols_tmp,4)
                        statname = sprintf('beta_hat (age*PDS = %.1f), %s', pdsvec(slice), designmat_file{des});
                        test = convertFEMAVols(vols_tmp(:,:,:,slice), modality, statname, colnames_model(jvec_bf), 1, agelim, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion);
                        if ~exist('vols') vols = test(2); else vols = cat(2, vols, test(2)); end
                    end
                end

                % age spline values
                if ~isempty(Sage_idx)
                    [vols_tmp beta_tmp] = FEMA_convert_splinevols(tbl_bf_age, agevec, beta_hat([Sage_idx],:), mask);
                    for slice=9:12:size(vols_tmp,4)
                        statname = sprintf('beta_hat (age = %.1f), %s', agevec(slice)/12, designmat_file{des});
                        test = convertFEMAVols(vols_tmp(:,:,:,slice), modality, statname, colnames_model(Sage_idx), 1, agelim, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion);
                        if ~exist('vols') vols = test(2); else vols = cat(2, vols, test(2)); end
                    end
                end

                if 0
                    % zmat for interaction - these z scores are generally 2-3 max in designmat 6f, less than 2 in 6m for RNT.
                    % calculate interaction spline values - zmat
                    jvec_bf = find(find_cell(regexp(colnames_model,'interview_age_bf_','match')));
                    [vols_tmp beta_tmp] = FEMA_convert_splinevols(tbl_bf, pdsvec, zmat([jvec_bf],:), mask);
                    for slice=1:size(vols_tmp,4)
                        statname = sprintf('z score (age*PDS = %.1f), %s', pdsvec(slice), designmat_file{des});
                        test = convertFEMAVols(vols_tmp(:,:,:,slice), modality, statname, colnames_model(jvec_bf), 1, [-5 5], showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion);
                        if ~exist('vols') vols = test(1); else vols = cat(2, vols, test(1)); end
                    end
                end
            end

            % run showVol
            coords = [96 87 127];
            showVol(vols, PRI.ABCD3.T1, PRI.ABCD3.CO, struct('roiatlas','ABCD3'), coords);
        end

        if strcmp(datatype,'vertex')


            legendPosition = 'south';
            ico=5; % ico number

            if ~isempty(age_idx)
                vertvals = beta_hat(age_idx,:);
                statname = sprintf('%s ~ %s (beta), %s', modality, colnames_model{age_idx}, designmat_file{des});

                % inputs for plotting
                fmax = min(300,max(abs(beta_tmp))); % max limit for plotting purposes
                fmin = -fmax; % min limit for plotting purposes
                fmid = 0; % middle value for plotting purposes
                fvals = [fmin fmid fmax]; % this will be passed to the SurfView_show_new function
                
                % set colorbar limits - usually [fmin fmax] or [-fmax fmax]
                clim = [fmin fmax];

                FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
                'colormap', blueblackred(), 'legendPosition', legendPosition, 'ico', ico, 'savepath', savepath);

            elseif ~isempty(Sage_idx)                
                [vols_tmp beta_tmp] = FEMA_convert_splinevols(tbl_bf_age, agevec, beta_hat(Sage_idx,:), mask);

                % inputs for plotting
                fmax = min(300,max(beta_tmp(:))); % max limit for plotting purposes
                fmin = -fmax; % min limit for plotting purposes
                fmid = 0; % middle value for plotting purposes
                fvals = [fmin fmid fmax]; % this will be passed to the SurfView_show_new function
                
                % set colorbar limits - usually [fmin fmax] or [-fmax fmax]
                clim = [fmin fmax];

                for slice=9:12:size(beta_tmp,1)
                    vertvals = beta_tmp(slice,:);
                    statname = sprintf('%s ~ age %.1f (beta), %s', modality, agevec(slice)/12, designmat_file{des});
                    savepath = sprintf('%s/%s_%s_age_%.1f.png',savedir, modality, designmat_file{des}, agevec(slice)/12);
                    FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
                    'colormap', blueblackred(), 'legendPosition', legendPosition, 'ico', ico, 'savepath', savepath);
                end

            end
           
        end

    end

end
