% Create pheno files for specific tabulated imaging data of interest
% Diana Smith
% March 2024

orig_dir = '/space/syn65/1/data/abcd-sync/5.0/tabulated/released/core/imaging';
savedir = '/space/syn50/1/data/ABCD/d9smith/age/volume/phenofiles';

orig_fname = 'mri_y_smr_vol_aseg.csv';

orig_file = fullfile(orig_dir,orig_fname);
savefile = fullfile(savedir,orig_fname);

roidat = readtable(orig_file); 
iid_concat = roidat.src_subject_id; % Not sure why imaging subject IDs are missing the NDAR_ part
eid_concat = roidat.eventname;
ymat=table2array(roidat(:,3:end));
colnames_imaging=roidat.Properties.VariableNames(3:end);

defvec = isfinite(sum(ymat, 1)) & (sum(ymat,1)~=0);
if any(defvec == 0)
    warning(sprintf('NaNs and/or columns of all zeros detected in ymat (%d columns removed).', length(find(defvec == 0))))
end

if size(ymat, 2) > sum(defvec)
    ymat = ymat(:, defvec); 
    colnames_imaging = colnames_imaging(defvec);
end

y_tbl = array2table(ymat);
y_tbl = renamevars(y_tbl,1:width(y_tbl),colnames_imaging);
y_tbl = [roidat(:,1:2) y_tbl];

writetable(y_tbl, savefile);