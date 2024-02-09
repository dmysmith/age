# Sanity check: run identical models in FEMA and GAMM (using R)
# Diana Smith
# Feb 2024

# Check what version of R we are using
R.version$version.string

# Load packages
for (p in c("plyr", "gamm4", "gratia")){
        if(!eval(parse(text=paste("require(",p,")")))) {
                install.packages(p)
                lapply(p,library,character.only=TRUE)
        }
}

# specify imaging outcome of interest
img_pheno = 'rsi_rni_aseg';
var_interest = 'dmri_rsirni_scs_pllh';

# specify paths
fname_design = '/space/syn50/1/data/ABCD/d9smith/age/results_2024-01-22/designMat4_BFsSexIncEducHispPCsScanSoftMotion_bly2y4.txt';
fname_basis = '/space/syn50/1/data/ABCD/d9smith/age/basis.txt';

outpath = "/home/d9smith/projects/age/plots/r";   

imaging_path = '/space/syn65/1/data/abcd-sync/5.0/tabulated/released/core/imaging'

# load tabulated imaging data and combine with design matrix
imaging_file = paste0(imaging_path, '/mri_y_', img_pheno, '.csv')
imaging = read.delim(imaging_file, sep = ',')

imaging_vars = c('src_subject_id', 'eventname', var_interest)
imaging = imaging[,imaging_vars]

tbl_design = read.delim(fname_design)
tbl_design = join(tbl_design, imaging, by = c('src_subject_id', 'eventname'))

plot_df = na.omit(tbl_design)
attach(plot_df)

# visualize variable of interest and fit smooth spline function
plot(age/12, plot_df[,var_interest])
scatter.smooth(age, plot_df[,var_interest], lpars = list(col = "blue", lwd = 3, lty = 3))
# abline(lm(var_interest~(age), col='red',lwd=3))
# legend('topright', c('Linear','Smoothing'), lty=c(1,2), lwd=c(3,3), col=c('red','blue')) 

# create formula and other inputs for model

# form = paste(var_interest, ' ~ s(age) + ', paste(names(plot_df)[9:76], collapse = ' + '))
# The above results in a failure to converge. Trying a simpler model using s(age)

form = paste(var_interest, ' ~ s(age)')
knots = c(125,150,175)
random='~(1|rel_family_id/src_subject_id)'
dat = plot_df

results = gamm4(formula=formula(form), 
                random=formula(random), 
                data=dat, 
                # knots = knots, 
                family=Gamma(link = "log"))

plot.gam(results$gam)
draw(results$gam)