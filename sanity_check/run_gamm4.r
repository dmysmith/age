# Sanity check: run identical models in FEMA and GAMM (using R)
# Diana Smith
# Feb 2024

# Check what version of R we are using
R.version$version.string

# Load packages
for (p in c("plyr", "gamm4", "ggplot2", "gratia", "R.utils")){
        if(!eval(parse(text=paste("require(",p,")")))) {
                install.packages(p)
                lapply(p,library,character.only=TRUE)
        }
}

# specify imaging outcome of interest
img_pheno = 'smri_vol_aseg';

vars_interest = c('smri_vol_scs_hpuslh', 'smri_vol_scs_hpusrh', 'smri_vol_scs_amygdalalh', 'smri_vol_scs_amygdalarh',
'smri_vol_scs_aal', 'smri_vol_scs_aar', 'smri_vol_scs_vedclh', 'smri_vol_scs_vedcrh',
'smri_vol_scs_ccps', 'smri_vol_scs_ccmidps', 'smri_vol_scs_ccct', 'smri_vol_scs_ccmidat', 'smri_vol_scs_ccat',
'smri_vol_scs_tplh', 'smri_vol_scs_tprh', 'smri_vol_scs_caudatelh', 'smri_vol_scs_caudaterh',
'smri_vol_scs_putamenlh', 'smri_vol_scs_putamenrh', 'smri_vol_scs_pallidumlh', 'smri_vol_scs_pallidumrh');
# vars_interest = c('smri_vol_scs_putamenlh')

# specify paths
fname_design = '/space/cluster/1/ABCD/users/d9smith/age/volume/designmat/designmat1_SAgeSexScanSoft.txt';
fname_basis = '/space/cluster/1/ABCD/users/d9smith/age/basis.txt';

outpath = "/home/d9smith/projects/age/plots/r";   

imaging_path = '/space/syn65/1/data/abcd-sync/6.0/tabulated/img'

# load tabulated imaging data and combine with design matrix
# imaging_file = paste0(imaging_path, '/mri_y_', img_pheno, '.csv')
imaging_file = paste0(imaging_path,'/abcd_smrip102.csv')
imaging = read.delim(imaging_file, sep = ',')

tbl_design_orig = read.delim(fname_design)

for (var in vars_interest){
        cat(paste(Sys.time(), 'Processing', var, '\n', sep = ' '))
        imaging_vars = c('src_subject_id', 'eventname', var)
        imagingvar = imaging[,imaging_vars]

        tbl_design = join(tbl_design_orig, imagingvar, by = c('src_subject_id', 'eventname'))

        plot_df = na.omit(tbl_design)
        attach(plot_df)

        # visualize variable of interest and fit smooth spline function
        # plot(age, plot_df[,var_interest])
        # scatter.smooth(age, plot_df[,var_interest], lpars = list(col = "blue", lwd = 3, lty = 3))
        # abline(lm(var_interest~(age), col='red',lwd=3))
        # legend('topright', c('Linear','Smoothing'), lty=c(1,2), lwd=c(3,3), col=c('red','blue')) 

        # create formula and other inputs for model

        # form = paste(var_interest, ' ~ s(age) + ', paste(names(plot_df)[9:76], collapse = ' + '))
        # The above results in a failure to converge. Trying a simpler model using s(age)

        cov = paste(names(plot_df)[9:(length(names(plot_df))-2)], collapse = ' + ')
        form = paste(var, ' ~ s(age, bs="cs",k=5) + ', cov) # cubic splines, 5 knots
        cat(paste(Sys.time(), 'Formula:', form, '\n', sep = ' '))

        random='~(1|rel_family_id/src_subject_id)'
        dat = plot_df

        cat(paste(Sys.time(), 'Running model on', var, 'using gamm4...\n', sep = ' '))


        tryCatch({
                withTimeout({
                results = gamm4(formula=formula(form), 
                        random=formula(random), 
                        # data=dat, 
                        family=Gamma(link = "log"))


                cat(paste(Sys.time(), 'Solution found. Plotting results...\n', sep = ' '))
                # plot.gam(results$gam)
                # draw(results$gam) # from gratia package

                plottitle = paste(var, ' ~ ', strsplit(strsplit(fname_design,"_")[[1]][2], '.txt'))
                # draw(results$gam) & ggtitle(plottitle)

                png(filename = paste0(outpath, '/', var, '.png'))
                draw(results$gam) & ggtitle(plottitle)
                dev.off()
                cat(paste(Sys.time(), 'File saved to', paste0(outpath, '/', var, '.png\n'), sep = ' '))
                }, timeout = 3600)
        }, TimeoutException = function(ex) {
                message(paste(Sys.time(), 'Timeout. Skipping.', sep = ' '))
        })

        # tryCatch({
        #         results = gamm4(formula=formula(form), 
        #                 random=formula(random), 
        #                 data=dat, 
        #                 # knots = knots, 
        #                 family=Gamma(link = "log"))


        #         cat(paste(Sys.time(), 'Solution found. Plotting results...\n', sep = ' '))
        #         # plot.gam(results$gam)
        #         # draw(results$gam) # from gratia package

        #         plottitle = paste(var, ' ~ ', strsplit(strsplit(fname_design,"_")[[1]][2], '.txt'))
        #         # draw(results$gam) & ggtitle(plottitle)

        #         png(filename = paste0(outpath, '/', var, '.png'))
        #         draw(results$gam) & ggtitle(plottitle)
        #         dev.off()
        #         cat(paste(Sys.time(), 'File saved to', paste0(outpath, '/', var, '.png\n'), sep = ' '))
        # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
       detach(plot_df)
}