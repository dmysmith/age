# Create table of demographic information for age paper
# Diana Smith
# April 2024

library(tableone)
library(plyr)

# Define the path to the directory where you would like to save out your design matrix 
outfile <- '/space/cluster/1/ABCD/users/d9smith/age/volume/tables/demographic_info.csv'

# Load the data 
nda6file <- '/space/cluster/1/ABCD/users/d9smith/age/nda6.0_bfs.txt'
nda6 <- read.delim(nda6file) 

# only include subjects which pass QC
idx_t1w_inc <- which(nda6$imgincl_t1w_include==1)
nda6_t1w_inc <- nda6[idx_t1w_inc,]

# df for subjects with data for all variables of interest
vars_interest = c('src_subject_id', 'eventname', 'rel_family_id', 'interview_age', 
                paste0('PC',1:10),'household.income_cont','sex', 'high.educ', 
                'mri_info_deviceserialnumber', 'mri_info_softwareversion',
                'pds_y_p_average', 'anthro_bmi_corr')

df = na.omit(nda6_t1w_inc[,vars_interest]) 

table(df$eventname)

# vertexwise data
myvars = c("interview_age", "sex", "household.income_cont","high.educ","pds_y_p_average", "anthro_bmi_corr")

catvars = c("sex","high.educ")

table <- CreateTableOne(data = df, vars = myvars, factorVars = catvars, strata = "eventname")

## Then prepare table for export
table_p <- print(table, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

## Save to a CSV file
write.csv(table_p, file = outfile)

# table1(~ interview_age + sex + household.income_cont + high.educ + pds_y_p_average + anthro_bmi_corr | eventname, data=df)

# table1(~ interview_age, data=df)