% Plot tabulated results from age analysis w/ basis functions
% Diana Smith
% December 2023

% specify imaging outcome of interest
% note: this is only important if you want to visualize just one plot,
% but I was originally using it to steal the iid and eid vecs
img_pheno = 'smr_vol_aseg';
img_pheno = 'smr_area_dsk';
var_interest = 'smri_vol_scs_tplh';
var_interest = 'smri_area_cdk_ptcatelh';
derivative=0;

fname_basis = '/space/cluster/1/ABCD/users/d9smith/age/basis.txt';
dirname_out = "/home/d9smith/projects/age/volume/plots/tabulated";   

% for running just one design matrix
fname_design = {'/space/cluster/1/ABCD/users/d9smith/age/volume/designmat/designmat1_SAgeSexScanSoft.txt'};

% for running all design matrices in a directory
project_dir = '/space/cluster/1/ABCD/users/d9smith/age/volume';
designmat_dir = sprintf('%s/designmat',project_dir);
results_date = '2024-04-08';

designmat_file = dir(sprintf('%s/designmat*.txt', designmat_dir));
designmat_file = {designmat_file.name}';
fname_design = strcat(designmat_dir, '/', designmat_file);

outdir_file = strrep(designmat_file, '.txt', '');
outpath = strcat(dirname_out,'/',outdir_file); 

for d=1:length(fname_design)
  tbl_design = readtable(fname_design{d});
  % figure; plot(tbl_design.age,tbl_design.bf_demean_1,'*')

  tbl_bf = readtable(fname_basis); agevals_tbl = linspace(100,220,121);
  bfmat = table2array(tbl_bf);
  [dummy dbfmat] = gradient(bfmat); dbfmat = dbfmat/(agevals_tbl(2)-agevals_tbl(1));
  % figure; plot(agevals_tbl,bfmat,'LineWidth',2); 
  % figure; plot(agevals_tbl,dbfmat,'LineWidth',2);

  results_dir = sprintf('%s/results_%s',project_dir,results_date); 
  load(sprintf('%s/%s/FEMA_wrapper_output_external_mri_y_%s.mat', results_dir,strrep(designmat_file{d},'.txt',''),img_pheno));
  % char(colnames_imaging)

  % Intersect data
  vec1 = strcat(iid,eid);
  vec2 = strcat(tbl_design.src_subject_id,tbl_design.eventname);
  [dummy IA IB] = intersect(vec1,vec2,'stable');
  agevec_tbl = tbl_design.age;
  agevec = agevec_tbl(IB);

  jvec_bf = find(find_cell(regexp(colnames_model,'^bf_demean','match')));
  intercept_idx = find(find_cell(regexp(colnames_model,'^intercept','match')));
  % figure; plot(agevec,X(:,jvec_bf),'*');
  if isempty(jvec_bf)
    disp(sprintf('No basis functions found in design matrix: %s',fname_design{d}));
    continue;
  end

  if 1 % to plot just one variable
      col_interest = find(find_cell(regexp(colnames_imaging,var_interest,'match')));

      % calculate spline function and variance
      valvec = NaN([size(beta_hat,2) length(agevals_tbl)]);

      L = zeros([length(agevals_tbl) length(colnames_model)]);
      L(:, jvec_bf) = bfmat; % values of basis functions for each x value, padded with zeros for other variables
      covB = coeffCovar(:,:,col_interest);
 
      spline_sig2 =  NaN([1 length(agevals_tbl)]); % variance of estimate
      spline_se =  NaN([1 length(agevals_tbl)]); % standard error of estimate
      spline_ci =  NaN([2 length(agevals_tbl)]); % 95 % confidence interval of estimate

      for agei = 1:length(agevals_tbl)
        if derivative 
            wvec = dbfmat(agei,:)'; 
            valvec(:,agei) = sum(beta_hat([jvec_bf],:).*wvec,1);
        else 
            wvec = bfmat(agei,:)';
            valvec(:,agei) = sum(beta_hat([jvec_bf intercept_idx],:).*[wvec; 1],1);
            spline_sig2(:,agei) = L(agei,:)*covB*L(agei,:)';
        end
      end
  
      % calculate se and ci from variance estimate
      spline_se = sqrt(spline_sig2);
      spline_ci = [valvec(col_interest,:) - 1.96*spline_se; valvec(col_interest,:) + 1.96*spline_se];
      
      % plot function and 95% CI
      x=agevals_tbl/12; y=valvec(col_interest,:);
      figure; plot(x,y,'LineWidth',2);
      hold on
      patch([x, flip(x)], [y-1.96*spline_se, flip(y+1.96*spline_se)], 'b', 'FaceAlpha',0.25, 'EdgeColor','none');
      hold off
      title(var_interest, 'Interpreter', 'none');
      xlabel('Age (years)');
      if derivative
          ylabel('$\sum{\hat{\beta}_{basis}\partial_{basis}}$','Interpreter','latex'); % if plotting derivatives
      else
          ylabel('Estimate'); % if plotting functions (not derivatives)
      end

  end

  % distribution of ages in the dataset
  % figure; histogram(tbl_design.age);

  %% create many figures
  results_file_list = dir(sprintf('%s/%s/FEMA_wrapper_output_external_mri_y_*.mat', results_dir,strrep(designmat_file{d}, '.txt','')));
  results_file = {results_file_list.name}';
  % disp(results_file)

  for filei=1:length(results_file)
    load(sprintf('%s/%s/%s',results_dir,strrep(designmat_file{d}, '.txt',''),results_file{filei}));
    
    % subdirectory to save figures
    filestr=strrep(strrep(results_file{filei},'FEMA_wrapper_output_external_mri_y_',''),'.mat','');
    savedir = sprintf('%s/%s',outpath{d},filestr);
    if ~exist(savedir, 'dir')
        mkdir(savedir)
    end

    % Intersect data
    vec1 = strcat(iid,eid); % from results file
    vec2 = strcat(tbl_design.src_subject_id,tbl_design.eventname);
    [dummy IA IB] = intersect(vec1,vec2,'stable');
    agevec_tbl = tbl_design.age;
    agevec = agevec_tbl(IB);
    
    jvec_bf = find(find_cell(regexp(colnames_model,'^bf_demean','match')));
    intercept_idx = find(find_cell(regexp(colnames_model,'^intercept','match')));

    % calculate spline function and variance
    valvec = NaN([size(beta_hat,2) length(agevals_tbl)]);

    L = zeros([length(agevals_tbl) length(colnames_model)]);
    L(:, jvec_bf) = bfmat; % values of basis functions for each x value, padded with zeros for other variables
    covB = coeffCovar(:,:,col_interest);
    
    for agei = 1:length(agevals_tbl)
      if derivative 
          wvec = dbfmat(agei,:)'; 
          valvec(:,agei) = sum(beta_hat([jvec_bf],:).*wvec,1);
      else 
          wvec = bfmat(agei,:)';
          valvec(:,agei) = sum(beta_hat([jvec_bf intercept_idx],:).*[wvec; 1],1);
      end
    end
    
    % create many figures
    for figi = 1:length(colnames_imaging)
        
        spline_sig2 =  NaN([1 length(agevals_tbl)]); % variance of estimate
        spline_se =  NaN([1 length(agevals_tbl)]); % standard error of estimate

        for agei=1:length(agevals_tbl)
            spline_sig2(:,agei) = L(agei,:)*covB*L(agei,:)';
        end

        x=agevals_tbl/12; y=valvec(figi,:);
        spline_se = sqrt(spline_sig2);

        % create plot
        figure('visible','off');
        plot(x,y,'LineWidth',2);
        hold on
        patch([x, flip(x)], [y-1.96*spline_se, flip(y+1.96*spline_se)], 'b', 'FaceAlpha',0.25, 'EdgeColor','none');
        hold off
        title(colnames_imaging(figi), 'Interpreter', 'none');
        xlabel('Age (years)');
        if derivative
          ylabel('$\sum{\hat{\beta}_{basis}\times\partial_{basis}}$','Interpreter','latex'); % if plotting derivative
        else
          ylabel('Estimate'); % if plotting functions (not derivatives)
        end
        savepath = sprintf('%s/%s.png',savedir,colnames_imaging{figi});
        saveas(gcf,savepath);
        disp(sprintf('Saved figure to %s',savepath));
        clf;
    end
  end
end