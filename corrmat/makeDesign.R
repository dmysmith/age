###############################################################
# Create design matrices for age / corrmat analysis
# Diana Smith
# March 2024

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
outpath <- '/space/syn50/1/data/ABCD/d9smith/age/corrmat/designmat'

# Load the data 
ndafile <- '/space/syn50/1/data/ABCD/d9smith/age/nda5.0_bfs.txt'
nda <- read.delim(ndafile)

# only include subjects which pass QC
idx_rsfmri_inc <- which(nda$imgincl_rsfmri_include==1)
nda_rsfmri_inc <- nda[idx_rsfmri_inc,]

###############################################################
# Design Matrix 1: age only
fname <- 'designmat1_AgeSexIncEducPCsScanSoftMotion.txt'
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c('interview_age', paste0('PC',1:10), 'household.income_cont', 'rsfmri_meanmotion')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_rsfmri_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 4: s(age)
fname <- 'designmat4_SAgeSexIncEducPCsScanSoftMotion.txt' 
outfile <- paste0(outpath, '/', fname) 

# Continuous variables
contvar <- c(names(nda)[grep('bf_',names(nda))], paste0('PC',1:10), 'household.income_cont', 'rsfmri_meanmotion')

# Categorical variables
catvar <- c('sex', 'high.educ', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified; specify in chronoligical order
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1')

# Note that default is set to demean=TRUE (demean continuous variables)
makeDesign(nda_rsfmri_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)