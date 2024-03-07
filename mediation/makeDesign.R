###############################################################

# Create design matrix for age and puberty analysis using ABCD 5.0 data
# Diana Smith
# Mar 2024

###############################################################
# install packages
for (p in c("plyr")){
        if(!eval(parse(text=paste("require(",p,")")))) {
                install.packages(p)
                lapply(p,library,character.only=TRUE)
        }
}

# First need to make sure that R know where to find the makeDesign function 
source('/home/d9smith/github/cmig_tools_internal/cmig_tools_utils/r/makeDesign.R')

# Load the data 
ndafile <- '/space/syn50/1/data/ABCD/d9smith/age/nda5.0_bfs.txt'
nda <- read.delim(ndafile)

# only include subjects which pass QC
idx_dmri_inc <- which(nda$imgincl_dmri_include==1)
nda_dmri_inc <- nda[idx_dmri_inc,]

# df for subjects with data for all variables of interest
vars_interest = c('src_subject_id', 'eventname', 'rel_family_id', 'interview_age', 
                paste0('PC',1:10),'household.income_cont','sex', 'high.educ', 
                'mri_info_deviceserialnumber', 'mri_info_softwareversion',
                'pds_y_ss_category_all', names(nda)[grep('bf_',names(nda))])

nda_dmri_inc_completecases = na.omit(nda_dmri_inc[,vars_interest])

# Define the path to the directory where you would like to save out your design matrix 
outpath <- '/space/syn50/1/data/ABCD/d9smith/age/mediation/designmat'

###############################################################
# Design Matrix 1: age only, no PDS
fname <- 'designmat1_AgeSexIncEducPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname) 

# fixed effects: age, sex, household income, parental education, Hispanic ethnicity 
# top 10 genetic principal components, scanner ID, MRI software version, 
# motion (average frame-wise displacement in mm)

# Continuous variables
contvar <- c('interview_age', paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 2: PDS only, no age
fname <- 'designmat2_PDSavgSexIncEducPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname) 

# fixed effects: age, sex, household income, parental education, Hispanic ethnicity 
# top 10 genetic principal components, scanner ID, MRI software version, 
# motion (average frame-wise displacement in mm)

# Continuous variables
contvar <- c('pds_y_ss_category_all', paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 3: age and PDS
fname <- 'designmat3_AgePDSavgSexIncEducPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname) 

# fixed effects: age, sex, household income, parental education, Hispanic ethnicity 
# top 10 genetic principal components, scanner ID, MRI software version, 
# motion (average frame-wise displacement in mm)

# Continuous variables
contvar <- c('interview_age', 'pds_y_ss_category_all', paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 4: s(age) only, no PDS
fname <- 'designmat4_SAgeSexIncEducPCsScanSoft.txt' 
outfile <- paste0(outpath, '/', fname) 

# fixed effects: age, sex, household income, parental education, Hispanic ethnicity 
# top 10 genetic principal components, scanner ID, MRI software version, 
# motion (average frame-wise displacement in mm)

# Continuous variables
contvar <- c(names(nda)[grep('bf_',names(nda))], paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 5: s(age) and PDS
fname <- 'designmat5_SAgePDSavgSexIncEducPCsScanSoft.txt' 
outfile <- paste0(outpath, '/', fname) 

# fixed effects: age, sex, household income, parental education, Hispanic ethnicity 
# top 10 genetic principal components, scanner ID, MRI software version, 
# motion (average frame-wise displacement in mm)

# Continuous variables
contvar <- c('interview_age', 'pds_y_ss_category_all', paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)
