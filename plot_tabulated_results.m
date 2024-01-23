% Plot tabulated results from age analysis w/ basis functions
% Diana Smith
% December 2023

% specify imaging outcome of interest
% note: this is only important if you want to visualize just one plot,
% but I was originally using it to steal the iid and eid vecs
img_pheno = 'rsi_rni_aseg';
var_interest = 'dmri_rsirni_scs_vdclh';

fname_basis = '/space/syn50/1/data/ABCD/d9smith/age/basis.txt';
dirname_out = "/home/d9smith/projects/age/plots/tabulated";   

% for running just one design matrix
fname_design = {'/space/syn50/1/data/ABCD/d9smith/age/results_2024-01-06/designMat4_BFsSexIncEducHispPCsScanSoftMotion_bly2y4.txt'};

% for running all design matrices in a directory
designmat_dir = '/space/syn50/1/data/ABCD/d9smith/age/results_2024-01-06'; 

designmat_file = dir(sprintf('%s/designMat*.txt', designmat_dir));
designmat_file = {designmat_file.name}';
fname_design = strcat(designmat_dir, '/', designmat_file);

outdir_file = strrep(designmat_file, '.txt', '');
outpath = strcat(dirname_out,'/',outdir_file); 

for d=1:length(fname_design)
  tbl_design = readtable(fname_design{d});
  % figure; plot(tbl_design.age,tbl_design.bf_demean_1,'*')

  tbl_bf = readtable(fname_basis); agevals_tbl = linspace(100,200,101);
  bfmat = table2array(tbl_bf);
  [dummy dbfmat] = gradient(bfmat); dbfmat = dbfmat/(agevals_tbl(2)-agevals_tbl(1));
  % figure; plot(agevals_tbl,bfmat,'LineWidth',2); 
  % figure; plot(agevals_tbl,dbfmat,'LineWidth',2);

  results_dir = strrep(fname_design{d},'.txt','');
  load(sprintf('%s/FEMA_wrapper_output_external_mri_y_%s.mat',results_dir,img_pheno));
  % char(colnames_imaging)

  % Intersect data
  vec1 = strcat(iid,eid);
  vec2 = strcat(tbl_design.src_subject_id,tbl_design.eventname);
  [dummy IA IB] = intersect(vec1,vec2,'stable');
  agevec_tbl = tbl_design.age;
  agevec = agevec_tbl(IB);

  jvec_bf = find(find_cell(regexp(colnames_model,'^bf_','match')));
  % figure; plot(agevec,X(:,jvec_bf),'*');
  if isempty(jvec_bf)
    disp(sprintf('No basis functions found in design matrix: %s',fname_design{d}));
    continue;
  end

  % compute sum of betas for basis functions
  valvec = NaN([size(beta_hat,2) length(agevals_tbl)]);
  for agei = 1:length(agevals_tbl)
    wvec = dbfmat(agei,:)';
    wvec = bfmat(agei,:)'; % comment me out to get derivative
    valvec(:,agei) = sum(beta_hat([jvec_bf],:).*wvec,1);
  end

  if 1 % to plot just one variable
    col_interest = find(find_cell(regexp(colnames_imaging,var_interest,'match')));

    figure; plot(agevals_tbl/12,valvec(col_interest,:),'LineWidth',2);
    title(var_interest, 'Interpreter', 'none');
    xlabel('Age (years)');
    ylabel('$\sum{\hat{\beta}_{basis}\partial_{basis}}$','Interpreter','latex');
  end

  % distribution of ages in the dataset
  % figure; histogram(tbl_design.age);

  %% create many figures
  results_file_list = dir(sprintf('%s/FEMA_wrapper_output_external_mri_y_*.mat', results_dir));
  results_file = {results_file_list.name}';
  % disp(results_file)

  for filei=1:length(results_file)
    load(sprintf('%s/%s',results_dir,results_file{filei}));
    
    % Intersect data
    vec1 = strcat(iid,eid); % from results file
    vec2 = strcat(tbl_design.src_subject_id,tbl_design.eventname);
    [dummy IA IB] = intersect(vec1,vec2,'stable');
    agevec_tbl = tbl_design.age;
    agevec = agevec_tbl(IB);
    
    jvec_bf = find(find_cell(regexp(colnames_model,'^bf_','match')));
    
    % compute sum of betas for basis functions
    valvec = NaN([size(beta_hat,2) length(agevals_tbl)]);
    for agei = 1:length(agevals_tbl)
      wvec = dbfmat(agei,:)';
      % wvec = bfmat(agei,:)';
      valvec(:,agei) = sum(beta_hat([jvec_bf],:).*wvec,1);
    end
    
    % subdirectory to save figures
    filestr=strrep(strrep(results_file{filei},'FEMA_wrapper_output_external_mri_y_',''),'.mat','');
    savedir = sprintf('%s/%s',outpath{d},filestr);
    if ~exist(savedir, 'dir')
        mkdir(savedir)
    end
    
    % create many figures
    for figi = 1:length(colnames_imaging) 
        figure('visible','off');
        plot(agevals_tbl/12,valvec(figi,:),'LineWidth',2);
        title(colnames_imaging(figi), 'Interpreter', 'none');
        xlabel('Age (years)');
        ylabel('$\sum{\hat{\beta}_{basis}\times\partial_{basis}}$','Interpreter','latex');
        savepath = sprintf('%s/%s.png',savedir,colnames_imaging{figi});
        saveas(gcf,savepath);
        disp(sprintf('Saved figure to %s',savepath));
        clf;
    end
  end
end

% /home/d9smith/projects/age/plots/tabulated/dti_td_fs_gwc_dst/ddtifp_1704.png 