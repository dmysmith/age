###############################################################

# Create design matrix for age analysis using volumetric data 
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
outpath <- '/space/cluster/1/ABCD/users/d9smith/age/volume/designmat'

# Load the data 
# ndafile <- '/space/cluster/1/ABCD/users/d9smith/age/nda5.0_bfs.txt'
# nda <- read.delim(ndafile)

nda6file <- '/space/cluster/1/ABCD/users/d9smith/age/nda6.0_bfs.txt'
nda6 <- read.delim(nda6file) 

# only include subjects which pass QC
# idx_t1w_inc <- which(nda$imgincl_t1w_include==1)
# nda_t1w_inc <- nda[idx_t1w_inc,]
idx_t1w_inc <- which(nda6$imgincl_t1w_include==1)
nda6_t1w_inc <- nda6[idx_t1w_inc,]

# df for subjects with data for all variables of interest
vars_interest = c('src_subject_id', 'eventname', 'rel_family_id', 'interview_age', 
                paste0('PC',1:10),'household.income_cont','sex', 'high.educ', 
                'mri_info_deviceserialnumber', 'mri_info_softwareversion',
                'pds_y_p_average', names(nda)[grep('bf_',names(nda))])

# nda_t1w_inc_completecases = na.omit(nda_t1w_inc[,vars_interest]) # NB: this includes two intersex subjects
# nda_t1w_inc_completemales = nda_t1w_inc_completecases[nda_t1w_inc_completecases$sex=='M',]
# nda_t1w_inc_completefemales = nda_t1w_inc_completecases[nda_t1w_inc_completecases$sex=='F',]

# nda6_t1w_inc_completecases = na.omit(nda6_t1w_inc[,names(nda6) %in% vars_interest])
nda6_t1w_inc_completecases = na.omit(nda6_t1w_inc[,vars_interest]) 
nda6_t1w_inc_completemales = nda6_t1w_inc_completecases[nda6_t1w_inc_completecases$sex=='M',]
nda6_t1w_inc_completefemales = nda6_t1w_inc_completecases[nda6_t1w_inc_completecases$sex=='F',]

###############################################################
# Design Matrix 1: brain ~ s(age) + sex + scan + soft
fname <- 'designmat1_SAgeSexScanSoft.txt'
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c(names(nda6)[grep('bf_demean',names(nda6))])

# Categorical variables
catvar <- c('sex', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1', '6_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda6_t1w_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 1m: brain ~ s(age) + scan + soft [MALES]
fname <- 'designmat1m_SAgeSexScanSoft_males.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda6_t1w_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 1f: brain ~ s(age) + scan + soft [FEMALES]
fname <- 'designmat1f_SAgeSexScanSoft_females.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('mri_info_deviceserialnumber', 'mri_info_softwareversion')  
makeDesign(nda6_t1w_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 2: brain ~ s(age) + sex + scan + soft + PCs 
fname <- 'designmat2_SAgeSexScanSoftPCs.txt'
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c(names(nda6)[grep('bf_demean',names(nda6))], paste0('PC',1:10))

# Categorical variables
catvar <- c('sex', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda6_t1w_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 2m: brain ~ s(age) + sex + scan + soft + PCs  [MALES]
fname <- 'designmat2m_SAgeSexScanSoftPCs_males.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda6_t1w_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 2f: brain ~ s(age) + sex + scan + soft + PCs  [FEMALES]
fname <- 'designmat2f_SAgeSexScanSoftPCs_females.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('mri_info_deviceserialnumber', 'mri_info_softwareversion')  
makeDesign(nda6_t1w_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 3: brain ~ s(age) + sex + scan + soft + inc + educ
fname <- 'designmat3_SAgeSexScanSoftIncEduc.txt'
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c(names(nda6)[grep('bf_demean',names(nda6))], 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda6_t1w_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 3m: brain ~ s(age) + sex + scan + soft + inc + educ [MALES]
fname <- 'designmat3m_SAgeSexScanSoftIncEduc_males.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda6_t1w_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 3f: brain ~ s(age) + sex + scan + soft + inc + educ [FEMALES]
fname <- 'designmat3f_SAgeSexScanSoftIncEduc_females.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')  
makeDesign(nda6_t1w_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 4: brain ~ s(age) + sex + scan + soft + inc + educ + PCs
fname <- 'designmat4_SAgeSexScanSoftIncEducPCs.txt'
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c(names(nda6)[grep('bf_demean',names(nda6))], paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda6_t1w_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 4m: brain ~ s(age) + sex + scan + soft + inc + educ + PCs [MALES]
fname <- 'designmat4m_SAgeSexScanSoftIncEducPCs_males.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda6_t1w_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 4f: brain ~ s(age) + sex + scan + soft + inc + educ + PCs [FEMALES]
fname <- 'designmat4f_SAgeSexScanSoftIncEducPCs_females.txt'
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')  
makeDesign(nda6_t1w_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 5: brain ~ s(age) + sex + scan + soft + inc + educ + PCs + s(PDS)
fname <- 'designmat5_SAgeSexScanSoftIncEducPCsSPDSavg.txt' 
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c(names(nda)[grep('bf_demean',names(nda))], names(nda)[grep('bf_pdsavg',names(nda))], paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda6_t1w_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 5m: brain ~ s(age) + sex + scan + soft + inc + educ + PCs + s(PDS) [MALES]
fname <- 'designmat5m_SAgeSexScanSoftIncEducPCsSPDSavg_males.txt' 
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda6_t1w_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 5f: brain ~ s(age) + sex + scan + soft + inc + educ + PCs + s(PDS) [FEMALES]
fname <- 'designmat5f_SAgeSexScanSoftIncEducPCsSPDSavg_females.txt' 
outfile <- paste0(outpath, '/', fname) 
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda6_t1w_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 6: brain ~ s(age) + sex + scan + soft + inc + educ + PCs + s(BMI)
fname <- 'designmat6_SAgeSexScanSoftIncEducPCsSBMI.txt' 
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c(names(nda)[grep('bf_demean',names(nda))], names(nda)[grep('bf_bmi',names(nda))], paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda6_t1w_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 6m: brain ~ s(age) + sex + scan + soft + inc + educ + PCs + s(BMI) [MALES]
fname <- 'designmat6m_SAgeSexScanSoftIncEducPCsSBMI_males.txt' 
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda6_t1w_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 6f: brain ~ s(age) + sex + scan + soft + inc + educ + PCs + s(BMI) [FEMALES]
fname <- 'designmat6f_SAgeSexScanSoftIncEducPCsSBMI_females.txt' 
outfile <- paste0(outpath, '/', fname) 
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda6_t1w_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 7: brain ~ s(age) + sex + scan + soft + inc + educ + PCs + s(PDS) + s(BMI)
fname <- 'designmat7_SAgeSexScanSoftIncEducPCsSPDSavgSBMI.txt' 
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c(names(nda)[grep('bf_demean',names(nda))], names(nda)[grep('bf_pdsavg',names(nda))], names(nda)[grep('bf_bmi',names(nda))], paste0('PC',1:10), 'household.income_cont')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda6_t1w_inc_completecases, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 7m: brain ~ s(age) + sex + scan + soft + inc + educ + PCs + s(PDS) + s(BMI) [MALES]
fname <- 'designmat7m_SAgeSexScanSoftIncEducPCsSBMI_males.txt' 
outfile <- paste0(outpath, '/', fname)
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda6_t1w_inc_completemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 7f: brain ~ s(age) + sex + scan + soft + inc + educ + PCs + s(PDS) + s(BMI) [FEMALES]
fname <- 'designmat7f_SAgeSexScanSoftIncEducPCsSBMI_females.txt' 
outfile <- paste0(outpath, '/', fname) 
catvar <- c('high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion') 
makeDesign(nda6_t1w_inc_completefemales, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)
