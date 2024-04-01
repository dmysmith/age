################################
# Create RDS for age analysis using ABCD 6.0 data
# Diana Smith
# July 2023

# Note: This script was written before all tabulated data was available -- 
# currently only includes imaging related variables.

################################
# The following R packages need to be loaded
for (p in c("tidyverse", "psych", "plyr", "dplyr", "pracma", "PerformanceAnalytics", "splines", "stringr")){
        if(!eval(parse(text=paste("require(",p,")")))) {
                install.packages(p)
                lapply(p,library,character.only=TRUE)
        }
}

################################
# This section defines input and output paths for files and functions called. 

# Define the path to the directory which contains the tabulated ABCD data 
inpath <- '/space/syn65/1/data/abcd-sync/6.0/tabulated/img'

# Define the path to the genetic PCs 
pcfile <- '/space/syn65/1/data/abcd-sync/5.0/genomics/abcd_gen_y_hat.tsv'

# Define the full path to the output RDS file 
outpath <- '/space/cluster/1/ABCD/users/d9smith/age'
fname <- 'nda6.0_bfs'
outmatfile <- paste0(outpath, '/', fname)

# Define the path to tge cmig_utils/r directory, R needs to be able to 
# parse functions from this directory
funcpath <- '/home/d9smith/github/cmig_tools_internal/cmig_tools_utils/r'
# The functionmakeDEAPdemos.R requires the path to the directory which 
# contains the tabulated ABCD data defined explicitly here
# deappath <- '/space/syn65/1/data/abcd-sync/5.0/tabulated/released/core/abcd-general'
deappath <- inpath

# Define the file names for the instruments from which we need to pull
# variables. 
# img_thk_file <- 'imaging/mri_y_smr_thk_dsk.csv'
# img_area_file <- 'imaging/mri_y_smr_area_dsk.csv'
# img_vol_file <- 'imaging/mri_y_smr_vol_aseg.csv'
# MRIinfofile <- 'imaging/mri_y_adm_info.csv'
# imgincfile <- 'imaging/mri_y_qc_incl.csv'
# motionfile <- 'imaging/mri_y_qc_motion.csv'
# pdspfile <- 'physical-health/ph_p_pds.csv'
# pdsyfile <- 'physical-health/ph_y_pds.csv'

MRIinfofile <- 'abcd_mri01.csv'
imgincfile <- 'abcd_imgincl01.csv'

# Define the full paths to these files 
# img_thk_file <- paste0(inpath, '/', img_thk_file)
# img_area_file <- paste0(inpath, '/', img_area_file)
# img_vol_file <- paste0(inpath, '/', img_vol_file)
# imgincfile <- paste0(inpath, '/', imgincfile)
# MRIinfofile <- paste0(inpath, '/', MRIinfofile)
# motionfile <- paste0(inpath, '/', motionfile)
# pdspfile <- paste0(inpath, '/', pdspfile)
# pdsyfile <- paste0(inpath, '/', pdsyfile) 

MRIinfofile <- paste0(inpath, '/', MRIinfofile)
imgincfile <- paste0(inpath, '/', imgincfile)

# pre-release data paths
deap6.0file <- '/space/abcd-sync/1/6.0/data_deap_tabulated.csv'
datadumpfile <- '/space/cluster/1/ABCD/users/d9smith/data/ABCD_6.0_datadump_20240327.csv'

################################

# R needs to parse two functions from cmig_utils/r. One is load.txt
# to load the .txt files and remove unnecessary 2nd row, the other is 
# makeDEAPdemos.R which collates SES variables into the format used by
# DEAP. 
source(paste0(funcpath, '/', 'makeDEAPdemos.R'))

################################
# Create the SES variables as coded by DEAP
# deap <- makeDEAPdemos(deappath)
# deap <- deap[ , -which(names(deap) %in% c("interview_date"))]
# Combine with the previously extracted variables
# outmat <- deap

################################
# Load the MRI info instrument  and extract the device serial number and software version 
# variables which are always needed as covariates when using imaging data  
MRIinfo <- read.csv(MRIinfofile)
MRIinfo <- MRIinfo[,c('src_subject_id','eventname', grep('mri_info', names(MRIinfo), value=TRUE))]
MRIinfo[,'idevent'] <- paste0(MRIinfo$src_subject_id, '_', MRIinfo$eventname)
MRIinfo <- MRIinfo[duplicated(MRIinfo$idevent)==FALSE,]
MRIinfo[which(MRIinfo$mri_info_deviceserialnumber==""),]<-NA
lvl <- unique(MRIinfo$mri_info_deviceserialnumber)
lvl <- lvl[is.na(lvl)==FALSE]
MRIinfo$mri_info_deviceserialnumber<-factor(MRIinfo$mri_info_deviceserialnumber, levels=lvl)
MRIinfo[which(MRIinfo$mri_info_softwareversion==""),]<-NA
lvl <- unique(MRIinfo$mri_info_softwareversion)
lvl <- lvl[is.na(lvl)==FALSE]
MRIinfo$mri_info_softwareversion<-factor(MRIinfo$mri_info_softwareversion, levels=lvl)
MRIinfo <- select(MRIinfo,  -c('idevent'))
# Combine with the previously extracted variables
outmat <- MRIinfo

################################
# Load the genetic PCs file 
pc_mat <- read.delim(pcfile)
# Get just the first 10 PCs and write to a dataframe
pc_names <- paste0('genetic_pc_',c(1:10))
pc <- data.frame(pc_mat[,c('src_subject_id',pc_names)])
names(pc) <- c('src_subject_id',paste0('PC',c(1:10)))
# Combine with the physical health variables. 
outmat <- join(outmat, pc, by='src_subject_id', match = "all")

# ################################
# # Create the SES variables as coded by DEAP
# deap <- makeDEAPdemos(deappath)
# deap <- deap[ , -which(names(deap) %in% c("interview_date"))]
# # Combine with the previously extracted variables
# outmat <- join(outmat, deap, by=c('src_subject_id', 'eventname'))

# ################################
# # Load imaging data files from tabulated data 
# img_thk <- read.csv(img_thk_file)
# # Extract intracranial volume  mean thickness and surface area  
# img_thk_vars <-c('src_subject_id', 'eventname', 'smri_thick_cdk_mean')
# img_thk <- img_thk[, img_thk_vars]
# # Combine with the previously extracted variables
# outmat <- join(outmat, img_thk, by=c('src_subject_id', 'eventname'))

# img_area <- read.csv(img_area_file)
# # Extract total surface area 
# img_area_vars <-c('src_subject_id', 'eventname', 'smri_area_cdk_total')
# img_area <- img_area[, img_area_vars]
# # Combine with the previously extracted variables
# outmat <- join(outmat, img_area, by=c('src_subject_id', 'eventname'))

# img_vol <- read.csv(img_vol_file)
# # Extract total surface area 
# img_vol_vars <-c('src_subject_id', 'eventname', 'smri_vol_scs_intracranialv')
# img_vol <- img_vol[, img_vol_vars]
# # Combine with the previously extracted variables
# outmat <- join(outmat, img_vol, by=c('src_subject_id', 'eventname'))

################################
# Include the MRI QC include/exclude variable 
imginc <- read.csv(imgincfile)

# Reformat "VisitID" variable to match src_subject_id and eventname
visitid = data.frame(str_split_fixed(imginc[,'VisitID'],"_",3))
visitid$src_subject_id = paste0('NDAR_',visitid[,'X2'])
visitid$eventname = case_match(visitid[,'X3'], 
'baseline' ~ 
'baseline_year_1_arm_1', 
'2year' ~ '2_year_follow_up_y_arm_1',
'4year' ~ '4_year_follow_up_y_arm_1',
'6year' ~ '6_year_follow_up_y_arm_1')

imginc = cbind(visitid[,c('src_subject_id', 'eventname')],imginc)

# Exctract the include/exclude variable for all imaging modalities 
imgincvar <- c('src_subject_id', 'eventname', grep('include', names(imginc), value=TRUE))
imginc <- imginc[, imgincvar]

# Combine with the previously extracted variables
outmat <- join(outmat, imginc, by=c('src_subject_id', 'eventname'))

# ################################
# # Include the MRI QC motion variable
# motion <- read.csv(motionfile)
# # Extract intracranial volume  mean thickness and surface area  
# motion_vars <-c('src_subject_id', 'eventname', 'dmri_meanmotion')
# motion <- motion[, motion_vars]
# # Combine with the previously extracted variables
# outmat <- join(outmat, motion, by=c('src_subject_id', 'eventname'))

################################
# hack for 6.0: I don't have family ID for 6.0 so I am using a file 
# from 5.0 (since family ID doesn't change)
famfile = '/space/syn65/1/data/abcd-sync/5.0/support_files/birth_id.txt'
fam = read.delim(famfile, sep = ' ')

# Get just the first 10 PCs and write to a dataframe
famnames <- c('pguid', 'update_family_id')
fam = fam[,famnames]
colnames(fam) = c('src_subject_id', 'rel_family_id')

# Combine with the physical health variables. 
outmat <- join(outmat, fam, by='src_subject_id', match = "all")

# Need to format design matrix so that first three columns are 
# src_subject_id, eventname, rel_family_id
colnames = c('src_subject_id', 'eventname', 'rel_family_id')
outmat = cbind(outmat[,colnames], outmat[,-which(names(outmat) %in% colnames)])

################################
# add variables from DEAP pre-release file
deap = read.csv(deap6.0file)

# rename subject id, eventname, age
deap$src_subject_id = deap$id_redcap
deap$eventname = deap$redcap_event_name
deap$interview_age = deap$age_visit

# sex
sextmp = data.frame(deap[deap$eventname=='baseline_year_1_arm_1',c('src_subject_id','demo_sex_v2b')])
sextmp$sex = recode(as.factor(sextmp$demo_sex_v2b), "1" = "M","2" = "F", "3" = "I")
deap<-join(deap,sextmp[,c('src_subject_id','sex')], by='src_subject_id', match = "all")

# calculate bmi and tmi
weightkg <- deap$anthro_weight_calc*0.453592
heightm <- deap$anthro_height_calc*0.0254
bmi <- weightkg/(heightm^2)
tmi <- weightkg/(heightm^3)
deap$anthro_bmi_calc <- bmi
deap$anthro_tmi_calc <- tmi
# remove biologically implausible values
ulim <- 45
llim <- 11
rm_bmi <- which(deap$anthro_bmi_calc>ulim | deap$anthro_bmi_calc<llim)
anthro_bmi_corr <- deap$anthro_bmi_calc
anthro_tmi_corr <- deap$anthro_tmi_calc
deap$anthro_bmi_corr <- anthro_bmi_corr
deap$anthro_tmi_corr <- anthro_tmi_corr
deap[rm_bmi,'anthro_bmi_corr'] <- NA
deap[rm_bmi,'anthro_tmi_corr'] <- NA

# parental education
deap[,'demo_prnt_ed_p']<-coalesce(deap$demo_prnt_ed_v2b,deap$demo_prnt_ed_v2_l)
deap[,'demo_prnt_ed_p']<-coalesce(deap$demo_prnt_ed_p,deap$demo_prnt_ed_v2_2yr_l)

# commented out because we don't have partner education for 6.0 yet
# deap$demo_prtnr_ed_v2_l = as.integer(deap$demo_prtnr_ed_v2_l)
# deap[,'demo_prtnr_ed_p']<-coalesce(deap$demo_prtnr_ed_v2,deap$demo_prtnr_ed_v2_l)
# deap[,'demo_prtnr_ed_p']<-coalesce(deap$demo_prtnr_ed_p,deap$demo_prtnr_ed_v2_2yr_l)

#highest education: 5 different levels. These levels correspond to the numbers published by the American Community Survey (ACS).
high.educ1 = deap$demo_prnt_ed_p
# high.educ2 = alldems$demo_prtnr_ed_p
high.educ1[which(high.educ1 == "999")] = NA
# high.educ2[which(high.educ2 == "999")] = NA
high.educ1[which(high.educ1 == "777")] = NA
# high.educ2[which(high.educ2 == "777")] = NA
high.educ1[which(high.educ1 == "22" | high.educ1=="23")] = 15 #22 and 23 = some college --> lower level than 18+
# high.educ2[which(high.educ2 == "22" | high.educ2=="23")] = 15
high.educ = pmax(as.numeric(as.character(high.educ1)), na.rm=T) # high.educ = pmax(as.numeric(as.character(high.educ1)), as.numeric(as.character(high.educ2)), na.rm=T)
idx <- which(high.educ %in% 0:12, arr.ind = TRUE)
high.educ[idx] = 1 # "< HS Diploma"
idx <- which(high.educ %in% 13:14, arr.ind = TRUE)
high.educ[idx] = 2 # "HS Diploma/GED"
idx <- which(high.educ %in% c(15:17,22:23), arr.ind = TRUE)
high.educ[idx] = 3 # "Some College"
idx <- which(high.educ == 18, arr.ind = TRUE)
high.educ[idx] = 4 # "Bachelor"
idx <- which(high.educ %in% 19:21, arr.ind = TRUE)
high.educ[idx] = 5 # "Post Graduate Degree"
high.educ[which(high.educ == "999")]=NA
high.educ[which(high.educ == "777")]=NA
deap$high.educ = factor( high.educ, levels= 1:5, labels = c("< HS Diploma","HS Diploma/GED","Some College","Bachelor","Post Graduate Degree") )

# household income
deap[,'demo_comb_income_p']<-coalesce(deap$demo_comb_income_v2b,deap$demo_comb_income_v2_l)

household.income = deap$demo_comb_income_p
household.income[deap$demo_comb_income_p == "1"] = 1 # "[<50K]"
household.income[deap$demo_comb_income_p == "2"] = 1 # "[<50K]"
household.income[deap$demo_comb_income_p == "3"] = 1 # "[<50K]"
household.income[deap$demo_comb_income_p == "4"] = 1 # "[<50K]"
household.income[deap$demo_comb_income_p == "5"] = 1 # "[<50K]"
household.income[deap$demo_comb_income_p == "6"] = 1 # "[<50K]"
household.income[deap$demo_comb_income_p == "7"] = 2 # "[>=50K & <100K]"
household.income[deap$demo_comb_income_p == "8"] = 2 # "[>=50K & <100K]"
household.income[deap$demo_comb_income_p == "9"] = 3 # "[>=100K]"
household.income[deap$demo_comb_income_p == "10"] = 3 # "[>=100K]"
household.income[deap$demo_comb_income_p == "777"] = NA
household.income[deap$demo_comb_income_p == "999"] = NA
household.income[household.income %in% c(NA, "999", "777")] = NA
deap$household.income = factor( household.income, levels= 1:3, labels = c("[<50K]", "[>=50K & <100K]", "[>=100K]") )

### Household income (continuous) - assign value based on middle of category
household.income_cont = deap$demo_comb_income_p
household.income_cont[deap$demo_comb_income_p == "1"] = 2500 # Less than $5,000
household.income_cont[deap$demo_comb_income_p == "2"] = 8500 # $5,000 through $11,999
household.income_cont[deap$demo_comb_income_p == "3"] = 14000 # $12,000 through $15,999
household.income_cont[deap$demo_comb_income_p == "4"] = 20500 # $16,000 through $24,999
household.income_cont[deap$demo_comb_income_p == "5"] = 30000 # $25,000 through $34,999;
household.income_cont[deap$demo_comb_income_p == "6"] = 42500 # $35,000 through $49,999
household.income_cont[deap$demo_comb_income_p == "7"] = 62500 # $50,000 through $74,999
household.income_cont[deap$demo_comb_income_p == "8"] = 87500 # $75,000 through $99,999
household.income_cont[deap$demo_comb_income_p == "9"] = 150000 # $100,000 through $199,999
household.income_cont[deap$demo_comb_income_p == "10"] = 250000 # $200,000 and greater
household.income_cont[deap$demo_comb_income_p == "777"] = NA # Refuse to answer
household.income_cont[deap$demo_comb_income_p == "999"] = NA # Don't know
household.income_cont[household.income_cont %in% c(NA, "999", "777")] = NA
deap$household.income_cont = household.income_cont

# Household income (10 level)
household.income_10level = deap$demo_comb_income_p
household.income_10level[deap$demo_comb_income_p == "777"] = NA # Refuse to answer
household.income_10level[deap$demo_comb_income_p == "999"] = NA # Don't know
household.income_10level[household.income_10level %in% c(NA, "999", "777")] = NA
deap$household.income_10level = household.income_10level

################################
# Pubertal Development, PDS
# average parent and youth report or use whichever report is available if only one informant

#Tanner stage categories
# MALES Prepubertal = 3; early Pubertal = 4 or 5 (no 3-point responses); Midpubertal = 6,7, or 8 (no 4-point responses; Late pubertal = 9-11; Postpubertal = 12. 

## FEMALES Prepubertal = 3; Early Puberty = 3 and no menarche; Midpubertal = 4 and no menarche; Late Puberty = <=7 and menarche; Postpubertal = 8 and menarche

#youth
deap[,'pds_y_ss_category_all']<-coalesce(deap$pds_y_ss_female_category_2, deap$pds_y_ss_male_category_2)
deap$pds_y_ss_category_all<-as.numeric(deap$pds_y_ss_category_all)

#parent
deap[,'pds_p_ss_category_all']<-coalesce(deap$pds_p_ss_female_category_2, deap$pds_p_ss_male_category_2)
deap$pds_p_ss_category_all<-as.numeric(deap$pds_p_ss_category_all)

#take average from parent and youth reports; if one is missing, take the non-missing value
deap$pds_y_p_average <- ifelse(is.na(deap$pds_y_ss_category_all), deap$pds_p_ss_category_all, ifelse(is.na(deap$pds_p_ss_category_all), deap$pds_y_ss_category_all, (deap$pds_y_ss_category_all + deap$pds_p_ss_category_all) / 2))

# merge deapvars with outmat
deap_vars = c('src_subject_id', 'eventname', 'sex', 'interview_age', 'anthro_bmi_corr', 'high.educ', 'household.income', 
'household.income_cont', 'household.income_10level', 'pds_y_p_average', 'pds_y_ss_category_all', 'pds_p_ss_category_all')
deap <- deap[, deap_vars] 

outmat <- join(outmat, deap, by=c('src_subject_id', 'eventname'))

################################
# basis functions for age
agevec = seq(from=100,to=200,length=101) # create age vector (in months)
knots = c(125,150,175) # should be same as default value
# basis <- data.frame(ns(agevec,df=dfs), row.names = agevec) # create basis functions
# colnames(basis) = paste0('bf_',c(1:dfs)) 

# save basis matrix
source('/home/d9smith/github/cmig_tools_internal/cmig_tools_utils/r/createBasis.R')
basis = createBasis(agevec,knots = knots, intercept = TRUE, demean = TRUE)
# library(Matrix)
# rankMatrix(basis) # should be equal to number of columns
write.table(basis, file = paste0(outpath, '/basis.txt'), sep = "\t", row.names = FALSE)

# Apply function to each row of interview_age in outmat
basis_values_df = get_basis_values(outmat,basis,'interview_age')
colnames(basis_values_df) <- colnames(basis)

# add basis functions to outmat
outmat <- cbind(outmat, basis_values_df)
# outmat <- cbind(outmat, bf_demeaned)

################################
# compute basis functions for puberty and save corresponding file
pdsvec = seq(from=0.5,to=5.5,length=11) # create PDS vector

# save basis matrix
source('/home/d9smith/github/cmig_tools_internal/cmig_tools_utils/r/createBasis.R')
basis_pds = createBasis(pdsvec, intercept = TRUE, demean = TRUE)
colnames(basis_pds) <- gsub("demean", "pdsavg", names(basis_pds))

# library(Matrix)
# rankMatrix(basis_pds) # should be equal to number of columns
write.table(basis_pds, file = paste0(outpath, '/basis_pds.txt'), sep = "\t", row.names = FALSE)

# Apply function to each row of interview_age in outmat
basis_pds_df = get_basis_values(outmat,basis_pds,'pds_y_p_average')
colnames(basis_pds_df) <- colnames(basis_pds)

# add basis functions to outmat
outmat <- cbind(outmat, basis_pds_df)

################################
# compute basis functions for BMI and save corresponding file
bmivec = seq(from=11,to=45,length=35) # create PDS vector

# save basis matrix
source('/home/d9smith/github/cmig_tools_internal/cmig_tools_utils/r/createBasis.R')
basis_bmi = createBasis(bmivec, intercept = TRUE, demean = TRUE)
colnames(basis_bmi) <- gsub("demean", "bmi", names(basis_bmi))

# library(Matrix)
# rankMatrix(basis_bmi) # should be equal to number of columns
write.table(basis_bmi, file = paste0(outpath, '/basis_bmi.txt'), sep = "\t", row.names = FALSE)

# Apply function to each row of interview_age in outmat
basis_bmi_df = get_basis_values(outmat,basis_bmi,'anthro_bmi_corr')
colnames(basis_bmi_df) <- colnames(basis_bmi)

# add basis functions to outmat
outmat <- cbind(outmat, basis_bmi_df)

################################
# Save the "outmat"

if ( ! dir.exists(outpath) ) {
        dir.create(outpath, recursive=TRUE)
}

write.table(outmat, file=paste0(outmatfile, '.txt'), sep = "\t", row.names = FALSE)
cat(paste0('File written to ', outmatfile, '.txt\n'))