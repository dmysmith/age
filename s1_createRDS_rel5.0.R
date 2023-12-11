################################

# Create RDS for age analysis using ABCD 5.0 data
# Diana Smith
# July 2023

# Note: This analysis is an extension of that found in Palmer 2022, 
# Microstructural development from 9 to 14 years: Evidence from the ABCD Study
# doi: https://doi.org/10.1016/j.dcn.2021.101044

################################
# The following R packages need to be loaded

library(tidyverse)
library(psych)
library(plyr)
library(dplyr)
library(PerformanceAnalytics)
library(pracma)
library(splines)

rm(list=ls())

################################
# This section defines input and output paths for files and functions called. 

# Define the path to the directory which contains the tabulated ABCD data 
inpath <- '/space/syn65/1/data/abcd-sync/5.0/tabulated/released/core'

# Define the path to the genetic PCs 
pcfile <- '/space/syn65/1/data/abcd-sync/5.0/genomics/abcd_gen_y_hat.tsv'

# Define the full path to the output RDS file 
outpath <- '/space/syn50/1/data/ABCD/d9smith/age'
fname <- 'nda5.0'
outmatfile <- paste0(outpath, '/', fname)

# Define the path to tge cmig_utils/r directory, R needs to be able to 
# parse functions from this directory
funcpath <- '/home/d9smith/github/cmig_tools_internal/cmig_tools_utils/r'
# The functionmakeDEAPdemos.R requires the path to the directory which 
# contains the tabulated ABCD data defined explicitly here
deappath <- '/space/syn65/1/data/abcd-sync/5.0/tabulated/released/core/abcd-general'

# Define the file names for the instruments from which we need to pull
# variables. 
img_thk_file <- 'imaging/mri_y_smr_thk_dsk.csv'
img_area_file <- 'imaging/mri_y_smr_area_dsk.csv'
img_vol_file <- 'imaging/mri_y_smr_vol_aseg.csv'
MRIinfofile <- 'imaging/mri_y_adm_info.csv'
imgincfile <- 'imaging/mri_y_qc_incl.csv'
motionfile <- 'imaging/mri_y_qc_motion.csv'
# Define the full paths to these files 
img_thk_file <- paste0(inpath, '/', img_thk_file)
img_area_file <- paste0(inpath, '/', img_area_file)
img_vol_file <- paste0(inpath, '/', img_vol_file)
imgincfile <- paste0(inpath, '/', imgincfile)
MRIinfofile <- paste0(inpath, '/', MRIinfofile)
motionfile <- paste0(inpath, '/', motionfile)

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

################################
# Create the SES variables as coded by DEAP
deap <- makeDEAPdemos(deappath)
deap <- deap[ , -which(names(deap) %in% c("interview_date"))]
# Combine with the previously extracted variables
outmat <- join(outmat, deap, by=c('src_subject_id', 'eventname'))

################################
# Load imaging data files from tabulated data 
img_thk <- read.csv(img_thk_file)
# Extract intracranial volume  mean thickness and surface area  
img_thk_vars <-c('src_subject_id', 'eventname', 'smri_thick_cdk_mean')
img_thk <- img_thk[, img_thk_vars]
# Combine with the previously extracted variables
outmat <- join(outmat, img_thk, by=c('src_subject_id', 'eventname'))

img_area <- read.csv(img_area_file)
# Extract total surface area 
img_area_vars <-c('src_subject_id', 'eventname', 'smri_area_cdk_total')
img_area <- img_area[, img_area_vars]
# Combine with the previously extracted variables
outmat <- join(outmat, img_area, by=c('src_subject_id', 'eventname'))

img_vol <- read.csv(img_vol_file)
# Extract total surface area 
img_vol_vars <-c('src_subject_id', 'eventname', 'smri_vol_scs_intracranialv')
img_vol <- img_vol[, img_vol_vars]
# Combine with the previously extracted variables
outmat <- join(outmat, img_vol, by=c('src_subject_id', 'eventname'))

################################
# Include the MRI QC include/exclude variable 
imginc <- read.csv(imgincfile)
# Exctract the include/exclude variable for all imaging modalities 
imgincvar <- c('src_subject_id', 'eventname', grep('include', names(imginc), value=TRUE))
imginc <- imginc[, imgincvar]

# As of 5.0 release there are 4 src_subject_id/eventname pairs that have multiple observations.
# Two of these have identical data so I can remove them easily:
imginc=distinct(imginc)

# For the other two cases, to be conservative I am assigning a value of 0 if 
# either of the duplicated rows has a value of 0 (i.e. removing the row that 
# has 1s when there is a difference)

# For now I am just hard coding the rows I want to remove.
remove_idx = c(which(imginc$src_subject_id=='NDAR_INV2F729N9A'&imginc$eventname=='baseline_year_1_arm_1'&imginc$imgincl_nback_include==0),
which(imginc$src_subject_id=='NDAR_INVMZZT8FFB'&imginc$eventname=='2_year_follow_up_y_arm_1'&imginc$imgincl_nback_include==0))

imginc = imginc[-remove_idx,]

# Combine with the previously extracted variables
outmat <- join(outmat, imginc, by=c('src_subject_id', 'eventname'))

################################
# Include the MRI QC motion variable
motion <- read.csv(motionfile)
# Extract intracranial volume  mean thickness and surface area  
motion_vars <-c('src_subject_id', 'eventname', 'dmri_meanmotion')
motion <- motion[, motion_vars]
# Combine with the previously extracted variables
outmat <- join(outmat, motion, by=c('src_subject_id', 'eventname'))

################################
# Need to format design matrix so that first three columns are 
# src_subject_id, eventname, rel_family_id
colnames = c('src_subject_id', 'eventname', 'rel_family_id')
outmat = cbind(outmat[,colnames], outmat[,-which(names(outmat) %in% colnames)])

################################
# Save the "outmat" as an RDS 

if ( ! dir.exists(outpath) ) {
        dir.create(outpath, recursive=TRUE)
}

# saveRDS(outmat, file=paste0(outmatfile, '.RDS'))
write.table(outmat, file=paste0(outmatfile, '.txt'), sep = "\t", row.names = FALSE)


################################
# compute basis functions and save corresponding file
agevec = linspace(100,200,101) # create age vector (in months)
dfs=4
basis <- data.frame(ns(agevec,df=dfs), row.names = agevec) # create basis functions
colnames(basis) = paste0('bf_',c(1:dfs)) 

# Function to extract values from basis dataframe based on interview_age values
extract_basis_values <- function(interview_age_value) {
  return(as.data.frame(t(basis[as.character(interview_age_value),])))
}

# Apply the function to each row of interview_age in outmat
basis_values <- lapply(outmat$interview_age, extract_basis_values)

# Combine basis values as columns in outmat dataframe
basis_values_df <- transpose(do.call(cbind, basis_values))
colnames(basis_values_df) <- colnames(basis)

# add basis functions to outmat
outmat <- cbind(outmat, basis_values_df)

# save outmat with bfs
write.table(outmat, file=paste0(outmatfile, '_bfs.txt'), sep = "\t", row.names = FALSE)
