fname_design = '/space/syn50/1/data/ABCD/d9smith/age/results_2023-12-11/designMat4_BFsSexIncEducHispPCsScanSoftMotion_bly2y4.txt';
tbl_design = readtable(fname_design);
% figure; plot(tbl_design.age,tbl_design.bf_demean_1,'*')

tbl_bf = readtable('/space/syn50/1/data/ABCD/d9smith/age/basis.txt'); agevals_tbl = linspace(100,200,101);
bfmat = table2array(tbl_bf);
[dummy dbfmat] = gradient(bfmat); dbfmat = dbfmat/(agevals_tbl(2)-agevals_tbl(1));
% figure; plot(agevals_tbl,bfmat,'LineWidth',2); 
% figure; plot(agevals_tbl,dbfmat,'LineWidth',2);

dirname_results = strrep(fname_design,'.txt','');
fstem_imaging = 'RNI';
fname_results = sprintf('%s/FEMA_wrapper_output_voxel_%s.mat',dirname_results, fstem_imaging);
load(fname_results);

% Intersect data -- shouldn't be needed
vec1 = strcat(iid,eid);
vec2 = strcat(tbl_design.src_subject_id,tbl_design.eventname);
[dummy IA IB] = intersect(vec1,vec2,'stable');
agevec_tbl = tbl_design.age;
agevec = agevec_tbl(IB);

jvec_bf = find(find_cell(regexp(colnames_model,'^bf_','match')));
% figure; plot(agevec,X(:,jvec_bf),'*');

% showVol(fullvol(beta_hat(1,:),mask),fullvol(beta_hat(2,:),mask),fullvol(beta_hat(3,:),mask),fullvol(beta_hat(4,:),mask));
if 0
    %valvec0 = sum(beta_hat.*mean(X)');
    vols = NaN([size(mask) length(agevals_tbl)]);
    for agei = 1:length(agevals_tbl)
      wvec = dbfmat(agei,:)';
      valvec = sum(beta_hat([jvec_bf],:).*wvec,1);
      vols(:,:,:,agei) = fullvol(valvec,mask);
    end

    showVol(1e6*vols(:,:,:,1:10:end));
end

%  for each IV, creates a raw grayscale map and colorized map overlaid on T1 
limits = [-0.04 0.04];  % Limits for colormap: (1) min (2) max (3) threshold (leave third value empty to have unthresholded maps)
index = [1];                    % Indicates which IVs to plot: colnames_model(index)
showPval = false;               % If showPval==true, creates a map of log10 p-values
cmap = redblackblue;            % Colormap
bg_vol = [];                    % Set background image to overlay statistics. Default: atlas T1
CSFprob = [];                   % Probabilistic threshold for masking CSF e.g. if CSFprob=0.8, only voxels in which 80% of participants labelled that voxel as CSF will be masked (Default: no masking)
interpMethod = [];              % By default uses linear interpolation to go from 2mm to 1mm voxels. Specify 'nearest' to show actual resolution
atlasVersion = 'ABCD2_cor10';

vol_stat = NaN([size(mask) length(agevals_tbl)]);
for agei = 1:length(agevals_tbl)
  wvec = dbfmat(agei,:)';
  valvec = sum(beta_hat([jvec_bf],:).*wvec,1);
  vol_stat(:,:,:,agei) = fullvol(valvec,mask);
  vols(:,:,agei) = convertFEMAVols(vol_stat(:,:,:,agei),fstem_imaging, 'beta_hat', colnames_model, index, limits, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion);
end

% showVol(vols(:,1,1:10:end), struct('roiatlas','ABCD2')) % black and white

showVol(vols(:,2,1:10:end), struct('roiatlas','ABCD2')) % color

% showVol(vols(:,:,1:10:end), struct('roiatlas','ABCD2')) % all

%% write video
write_video=0; % toggle to write video
if write_video
    zcoord = 100; % specify which coronal slice you want
    cmap = redblackblue; cmap(1,:) = 0; % hack to fix background color to black
    vidfile = VideoWriter('testmovie');
    open(vidfile);
    for ind = 1:size(vols,3)
        im = sc(vols{:,2,ind}.imgs(:,:,zcoord), cmap);
        writeVideo(vidfile, im);
    end
    close(vidfile)
end