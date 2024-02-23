addpath(genpath('/home/d9smith/github/cmig_tools_internal'));

fname_design = '/space/syn50/1/data/ABCD/d9smith/age/results_2024-01-22/designMat4_BFsSexIncEducHispPCsScanSoftMotion_bly2y4.txt';
tbl_design = readtable(fname_design);

tbl_bf = readtable('/space/syn50/1/data/ABCD/d9smith/age/basis.txt'); 
agevals_tbl = linspace(100,200,101);

fname_results = '/space/syn50/1/data/ABCD/d9smith/age/results_2024-01-22/designMat4_BFsSexIncEducHispPCsScanSoftMotion_bly2y4/FEMA_wrapper_output_voxel_RNT.mat';
load(fname_results);

%% Intersect data
vec1 = strcat(iid,eid);
vec2 = strcat(tbl_design.src_subject_id,tbl_design.eventname);
[dummy IA IB] = intersect(vec1,vec2,'stable');
agevec_tbl = tbl_design.age;
agevec = agevec_tbl(IB);

jvec_bf = find(find_cell(regexp(colnames_model,'^bf_','match')));

outpath = '/home/d9smith/tmp/save_splinevols_test';

vols = FEMA_convert_splinevols(tbl_bf, agevals_tbl, beta_hat([jvec_bf],:), mask, outpath);

showVol(1e6*vols(:,:,:,1:10:end));



