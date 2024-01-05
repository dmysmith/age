fname_design = '/space/syn50/1/data/ABCD/d9smith/age/results_2023-12-11/designMat4_BFsSexIncEducHispPCsScanSoftMotion_bly2y4.txt';
tbl_design = readtable(fname_design);
figure; plot(tbl_design.age,tbl_design.bf_demean_1,'*')

tbl_bf = readtable('/space/syn50/1/data/ABCD/d9smith/age/basis.txt'); agevals_tbl = linspace(100,200,101);
figure; plot(agevals_tbl,table2array(tbl_bf),'LineWidth',2);

dirname_results = strrep(fname_design,'.txt','');
fname_results = sprintf('%s/FEMA_wrapper_output_voxel_RNI.mat',dirname_results);
load(fname_results);


% Intersect data
vec1 = strcat(iid,eid);
vec2 = strcat(tbl_design.src_subject_id,tbl_design.eventname);
[dummy IA IB] = intersect(vec1,vec2,'stable');
agevec_tbl = tbl_design.age;
agevec = agevec_tbl(IB);

jvec_bf = find(find_cell(regexp(colnames_model,'^bf_','match')));

figure; plot(agevec,X(:,jvec_bf),'*'); % Fixed this plot

vol_beta_hat = single(fullvol(beta_hat(1,:),mask),fullvol(beta_hat(2,:),mask),fullvol(beta_hat(3,:),mask),fullvol(beta_hat(4,:),mask));

showVol(fullvol(beta_hat(1,:),mask),fullvol(beta_hat(2,:),mask),fullvol(beta_hat(3,:),mask),fullvol(beta_hat(4,:),mask));

