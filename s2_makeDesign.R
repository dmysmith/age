###############################################################

# Create design matrix for age analysis using ABCD 5.0 data
# Diana Smith
# July 2023

# Note: This analysis is an extension of that found in Palmer 2022, 
# Microstructural development from 9 to 14 years: Evidence from the ABCD Study
# doi: https://doi.org/10.1016/j.dcn.2021.101044

###############################################################
# First need to make sure that R know where to find the makeDesign function 
source('/home/d9smith/github/cmig_tools_internal/cmig_tools_utils/r/makeDesign.R')

# Load the data 
ndafile <- '/space/syn50/1/data/ABCD/d9smith/age/nda5.0_bfs.txt'
# nda <- read.delim(ndafile, header = T, sep = ",")
nda <- read.delim(ndafile)

# only include subjects which pass QC
idx_dmri_inc <- which(nda$imgincl_dmri_include==1|nda$eventname=='4_year_follow_up_y_arm_1')
nda_dmri_inc <- nda[idx_dmri_inc,]

# males and females separately
nda_dmri_inc_male <- nda_dmri_inc[nda_dmri_inc$sex=='M',]
nda_dmri_inc_female <- nda_dmri_inc[nda_dmri_inc$sex=='F',]


# df for subjects with data for all three timepoints
counts = nda_dmri_inc %>% dplyr::count(src_subject_id)
complete_ids = counts[counts$n==3,'src_subject_id']

complete_idx <- which(nda_dmri_inc$src_subject_id %in% complete_ids)
nda_dmri_inc_completecases <- nda_dmri_inc[complete_idx,] 


# Define the path to the directory where you would like to save out your design matrix 
outpath <- '/space/syn50/1/data/ABCD/d9smith/age/results_2024-01-06'

###############################################################
# Design Matrix 1: Replicating model from Clare's paper (baseline and y2 only)
# Note: This is using 5.0 data so there should be more participants than in Palmer 2022

# Define the name of your design matrix file 
fname <- 'designMat1_AgeSexIncEducHispPCsScanSoftMotion_bly2.txt'
# and it's full path to your fave output directory
outfile <- paste0(outpath, '/', fname)

# fixed effects: age, sex, household income, parental education, Hispanic ethnicity 
# top 10 genetic principal components, scanner ID, MRI software version, 
# motion (average frame-wise displacement in mm)

# Continuous variables
contvar <- c('interview_age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10', 'dmri_meanmotion', 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'hisp', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

############################################################
# Design Matrix 2: y2 and y4 only
fname <- 'designMat2_AgeSexIncEducHispPCsScanSoftMotion_y2y4.txt'
outfile <- paste0(outpath, '/', fname) 

# fixed effects: age, sex, household income, parental education, Hispanic ethnicity 
# top 10 genetic principal components, scanner ID, MRI software version, 
# motion (average frame-wise displacement in mm)

# Continuous variables
contvar <- c('interview_age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10', 'dmri_meanmotion', 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'hisp', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

############################################################
# Design Matrix 3: bl, y2 and y4
fname <- 'designMat3_AgeSexIncEducHispPCsScanSoftMotion_bly2y4.txt'
outfile <- paste0(outpath, '/', fname) 

# fixed effects: age, sex, household income, parental education, Hispanic ethnicity 
# top 10 genetic principal components, scanner ID, MRI software version, 
# motion (average frame-wise displacement in mm)

# Continuous variables
contvar <- c('interview_age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10', 'dmri_meanmotion', 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'hisp', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

############################################################
# Design Matrix 4: basis functions, all timepoints
fname <- 'designMat4_BFsSexIncEducHispPCsScanSoftMotion_bly2y4.txt'
outfile <- paste0(outpath, '/', fname) 

# fixed effects: age, sex, household income, parental education, Hispanic ethnicity 
# top 10 genetic principal components, scanner ID, MRI software version, 
# motion (average frame-wise displacement in mm)

# Continuous variables
contvar <- c(paste0('bf_demean_',1:4), 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10', 'dmri_meanmotion', 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'hisp', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

############################################################
# Design Matrix 4f: basis functions, all timepoints, females only
fname <- 'designMat4f_BFsSexIncEducHispPCsScanSoftMotion_bly2y4.txt'
outfile <- paste0(outpath, '/', fname) 

# fixed effects: age, sex, household income, parental education, Hispanic ethnicity 
# top 10 genetic principal components, scanner ID, MRI software version, 
# motion (average frame-wise displacement in mm)

# Continuous variables
contvar <- c(paste0('bf_demean_',1:4), 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10', 'dmri_meanmotion', 'household.income_cont')

# Categorical variables
catvar <- c('high.educ', 'hisp', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_female, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

############################################################
# Design Matrix 4m: basis functions, all timepoints, males only
fname <- 'designMat4m_BFsSexIncEducHispPCsScanSoftMotion_bly2y4.txt'
outfile <- paste0(outpath, '/', fname) 

# fixed effects: age, sex, household income, parental education, Hispanic ethnicity 
# top 10 genetic principal components, scanner ID, MRI software version, 
# motion (average frame-wise displacement in mm)

# Continuous variables
contvar <- c(paste0('bf_demean_',1:4), 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10', 'dmri_meanmotion', 'household.income_cont')

# Categorical variables
catvar <- c('high.educ', 'hisp', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_male, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)


############################################################
# Design Matrix 5: basis functions, all timepoints, plus PDS
fname <- 'designMat5_BFsSexIncEducHispPCsScanSoftMotionPDS_bly2y4.txt'
outfile <- paste0(outpath, '/', fname) 

# fixed effects: age, sex, household income, parental education, Hispanic ethnicity 
# top 10 genetic principal components, scanner ID, MRI software version, 
# motion (average frame-wise displacement in mm)

# Continuous variables
contvar <- c(paste0('bf_demean_',1:4), 'pds_avg', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10', 'dmri_meanmotion', 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'hisp', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

############################################################
# Design Matrix 5f: basis functions, all timepoints, plus PDS, females only
fname <- 'designMat5f_BFsSexIncEducHispPCsScanSoftMotionPDS_bly2y4.txt'
outfile <- paste0(outpath, '/', fname) 

# fixed effects: age, sex, household income, parental education, Hispanic ethnicity 
# top 10 genetic principal components, scanner ID, MRI software version, 
# motion (average frame-wise displacement in mm)

# Continuous variables
contvar <- c(paste0('bf_demean_',1:4), 'pds_avg', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10', 'dmri_meanmotion', 'household.income_cont')

# Categorical variables
catvar <- c('high.educ', 'hisp', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_female, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

############################################################
# Design Matrix 5m: basis functions, all timepoints, plus PDS, males only
fname <- 'designMat5m_BFsSexIncEducHispPCsScanSoftMotionPDS_bly2y4.txt'
outfile <- paste0(outpath, '/', fname) 

# fixed effects: age, sex, household income, parental education, Hispanic ethnicity 
# top 10 genetic principal components, scanner ID, MRI software version, 
# motion (average frame-wise displacement in mm)

# Continuous variables
contvar <- c(paste0('bf_demean_',1:4), 'pds_avg', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10', 'dmri_meanmotion', 'household.income_cont')

# Categorical variables
catvar <- c('high.educ', 'hisp', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_male, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)