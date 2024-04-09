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

# Define the path to the directory where you would like to save out your design matrix 
outpath <- '/space/syn50/1/data/ABCD/d9smith/age/mediation/designmat'

# Load the data 
ndafile <- '/space/syn50/1/data/ABCD/d9smith/age/nda5.0_bfs.txt'
nda <- read.delim(ndafile)

nda6file <- '/space/syn50/1/data/ABCD/d9smith/age/nda6.0_bfs.txt'
nda6 <- read.delim(nda6file) 

# only include subjects which pass QC
idx_dmri_inc <- which(nda$imgincl_dmri_include==1)
nda_dmri_inc <- nda[idx_dmri_inc,]
idx_dmri_inc <- which(nda6$imgincl_dmri_include==1)
nda6_dmri_inc <- nda6[idx_dmri_inc,]

# df for subjects with data for all variables of interest
vars_interest = c('src_subject_id', 'eventname', 'rel_family_id', 'interview_age', 
                paste0('PC',1:10),'household.income_cont','sex', 'high.educ', 
                'mri_info_deviceserialnumber', 'mri_info_softwareversion',
                'pds_y_p_average', names(nda)[grep('bf_',names(nda))])

nda_dmri_inc_completecases = na.omit(nda_dmri_inc[,vars_interest]) # NB: this includes two intersex subjects
nda_dmri_inc_completemales = nda_dmri_inc_completecases[nda_dmri_inc_completecases$sex=='M',]
nda_dmri_inc_completefemales = nda_dmri_inc_completecases[nda_dmri_inc_completecases$sex=='F',]

nda6_dmri_inc_completecases = na.omit(nda6_dmri_inc[,names(nda6) %in% vars_interest]) 
nda6_dmri_inc_completemales = nda6_dmri_inc_completecases[nda6_dmri_inc_completecases$sex=='M',]
nda6_dmri_inc_completefemales = nda6_dmri_inc_completecases[nda6_dmri_inc_completecases$sex=='F',]

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
fname <- 'designmat1m_AgeIncEducPCsScanSoft_males.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda_dmri_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 1f: age only, no PDS [FEMALES]
fname <- 'designmat1f_AgeIncEducPCsScanSoft_females.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')  
makeDesign(nda_dmri_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 2: PDS only, no age
fname <- 'designmat2_PDSavgSexIncEducPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c('pds_y_p_average', paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 2m: PDS only, no age [MALES]
fname <- 'designmat2m_PDSavgIncEducPCsScanSoft_males.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda_dmri_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 2f: PDS only, no age [FEMALES]
fname <- 'designmat2f_PDSavgIncEducPCsScanSoft_females.txt'
outfile <- paste0(outpath, '/', fname) 
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda_dmri_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 3: age and PDS
fname <- 'designmat3_AgePDSavgSexIncEducPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c('interview_age', 'pds_y_p_average', paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 3m: age and PDS [MALES]
fname <- 'designmat3m_AgePDSavgIncEducPCsScanSoft_males.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda_dmri_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 3f: age and PDS [FEMALES]
fname <- 'designmat3f_AgePDSavgIncEducPCsScanSoft_females.txt'
outfile <- paste0(outpath, '/', fname) 
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda_dmri_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 4: s(PDS) only, no age
fname <- 'designmat4_SPDSavgSexIncEducPCsScanSoft.txt' 
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c(names(nda)[grep('bf_pdsavg',names(nda))], paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 4m: s(PDS) only, no age [MALES]
fname <- 'designmat4m_SPDSavgSexIncEducPCsScanSoft_males.txt' 
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda_dmri_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 4f: s(PDS) only, no age [FEMALES]
fname <- 'designmat4f_SPDSavgSexIncEducPCsScanSoft_females.txt' 
outfile <- paste0(outpath, '/', fname) 
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda_dmri_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 5: age and s(PDS)
fname <- 'designmat5_AgeSPDSavgSexIncEducPCsScanSoft.txt' 
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c('interview_age', names(nda)[grep('bf_pdsavg',names(nda))], paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 5m: age and s(PDS) [MALES]
fname <- 'designmat5m_AgeSPDSavgSexIncEducPCsScanSoft_males.txt' 
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda_dmri_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 5f: age and s(PDS) [FEMALES]
fname <- 'designmat5f_AgeSPDSavgSexIncEducPCsScanSoft_females.txt' 
outfile <- paste0(outpath, '/', fname) 
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda_dmri_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 6: age * s(PDS)
fname <- 'designmat6_AgexSPDSavgSexIncEducPCsScanSoft.txt' 
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c('interview_age', names(nda)[grep('bf_pdsavg',names(nda))], paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

interact <- paste0('interview_age*',names(nda)[grep('bf_pdsavg',names(nda))])

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=interact, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 6m: age and s(PDS) [MALES]
fname <- 'designmat6m_AgexSPDSavgSexIncEducPCsScanSoft_males.txt' 
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda_dmri_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=interact, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 6f: age and s(PDS) [FEMALES]
fname <- 'designmat6f_AgexSPDSavgSexIncEducPCsScanSoft_females.txt' 
outfile <- paste0(outpath, '/', fname) 
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda_dmri_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=interact, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 7: age only, no PDS, in 6.0
fname <- 'designmat7_6.0_AgeSexPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c('interview_age', paste0('PC',1:10))

# Categorical variables
catvar <- c('sex', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1', '6_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 7m: age only, no PDS, in 6.0 [MALES]
fname <- 'designmat7m_6.0_AgeSexPCsScanSoft_males.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda_dmri_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 7f: age only, no PDS, in 6.0 [FEMALES]
fname <- 'designmat7f_6.0_AgeSexPCsScanSoft_females.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')  
makeDesign(nda_dmri_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 8: s(age) only, no PDS, in 6.0
fname <- 'designmat8_6.0_SAgeSexPCsScanSoft.txt'
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c(names(nda6)[grep('bf_demean',names(nda6))], paste0('PC',1:10))

# Categorical variables
catvar <- c('sex', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1', '6_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_dmri_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 8m: s(age) only, no PDS, in 6.0 [MALES]
fname <- 'designmat8m_6.0_SAgeSexPCsScanSoft_males.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda_dmri_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 8f: s(age) only, no PDS, in 6.0 [FEMALES]
fname <- 'designmat8f_6.0_SAgeSexPCsScanSoft_females.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')  
makeDesign(nda_dmri_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)
