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

nda_dmri_inc_completecases = na.omit(nda_dmri_inc[,vars_interest]) # NB: this includes two intersex subjects
nda_dmri_inc_completemales = nda_dmri_inc_completecases[nda_dmri_inc_completecases$sex=='M',]
nda_dmri_inc_completefemales = nda_dmri_inc_completecases[nda_dmri_inc_completecases$sex=='F',]

# Define the path to the directory where you would like to save out your design matrix 
outpath <- '/space/syn50/1/data/ABCD/d9smith/age/mediation/designmat'

###############################################################
# Design Matrix 1: age only, no PDS
fname <- 'designmat1_AgeSexIncEducPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c('interview_age', paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 1m: age only, no PDS [MALES]
fname <- 'designmat1m_AgeSexIncEducPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname) 
makeDesign(nda_dmri_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 1f: age only, no PDS [FEMALES]
fname <- 'designmat1f_AgeSexIncEducPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname) 
makeDesign(nda_dmri_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 2: PDS only, no age
fname <- 'designmat2_PDSySexIncEducPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c('pds_y_ss_category_all', paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 2m: PDS only, no age [MALES]
fname <- 'designmat2m_PDSySexIncEducPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname)
makeDesign(nda_dmri_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 2f: PDS only, no age [FEMALES]
fname <- 'designmat2f_PDSySexIncEducPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname) 
makeDesign(nda_dmri_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 3: age and PDS
fname <- 'designmat3_AgePDSySexIncEducPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c('interview_age', 'pds_y_ss_category_all', paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 3m: age and PDS [MALES]
fname <- 'designmat3m_AgePDSySexIncEducPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname)
makeDesign(nda_dmri_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 3f: age and PDS [FEMALES]
fname <- 'designmat3f_AgePDSySexIncEducPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname) 
makeDesign(nda_dmri_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 4: s(age) only, no PDS
fname <- 'designmat4_SAgeSexIncEducPCsScanSoft.txt' 
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c(names(nda)[grep('bf_',names(nda))], paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 4m: s(age) only, no PDS [MALES]
fname <- 'designmat4m_SAgeSexIncEducPCsScanSoft.txt' 
outfile <- paste0(outpath, '/', fname)
makeDesign(nda_dmri_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 4f: s(age) only, no PDS [FEMALES]
fname <- 'designmat4f_SAgeSexIncEducPCsScanSoft.txt' 
outfile <- paste0(outpath, '/', fname) 
makeDesign(nda_dmri_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 5: s(age) and PDS
fname <- 'designmat5_SAgePDSySexIncEducPCsScanSoft.txt' 
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c('interview_age', 'pds_y_ss_category_all', paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 5m: s(age) and PDS [MALES]
fname <- 'designmat5m_SAgePDSySexIncEducPCsScanSoft.txt' 
outfile <- paste0(outpath, '/', fname)
makeDesign(nda_dmri_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 5f: s(age) and PDS [FEMALES]
fname <- 'designmat5f_SAgePDSySexIncEducPCsScanSoft.txt' 
outfile <- paste0(outpath, '/', fname) 
makeDesign(nda_dmri_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)