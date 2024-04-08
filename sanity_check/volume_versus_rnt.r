# comparison of pallidum volume and RNT
# Diana Smith
# April 2024

require(ggplot2)

# specify path for saving figures
outpath = "/home/d9smith/projects/age/volume";   

fname_design = '/space/cluster/1/ABCD/users/d9smith/age/volume/designmat/designmat1_SAgeSexScanSoft.txt';
df = read.delim(fname_design)

# load tabulated volume data and combine with design matrix
imaging_path = '/space/syn65/1/data/abcd-sync/6.0/tabulated/img'
vol_file = paste0(imaging_path,'/abcd_smrip102.csv')
vol = read.delim(vol_file, sep = ',')

volvars = c('src_subject_id', 'eventname', 'smri_vol_scs_pallidumlh', 'smri_vol_scs_pallidumrh',
            'smri_vol_scs_intracranialv')
vol = vol[,volvars]

df = join(df, vol, by = c('src_subject_id', 'eventname'))

# load tabulated RSI data and combine with design matrix
rsi_file = paste0(imaging_path,'/abcd_drsip301.csv')
rsi = read.delim(rsi_file, sep=',')

rsivars = c('src_subject_id', 'eventname', 'dmri_rsirnt_scs_pllh', 'dmri_rsirnt_scs_plrh')
rsi = rsi[,rsivars]

df = join(df, rsi, by = c('src_subject_id', 'eventname'))

# create global volume-adjusted volume estimate for left and right pallidum
df$pallidum_lh_res <- residuals(lm(smri_vol_scs_pallidumlh ~ smri_vol_scs_intracranialv, data = df))  
df$pallidum_rh_res <- residuals(lm(smri_vol_scs_pallidumrh ~ smri_vol_scs_intracranialv, data = df))

# random sample of 1000 observations
df.sampled = df[sample(1:nrow(df),1000,replace=FALSE),]

# scatterplot of intracranialv-adjusted volumes by rnt
p <- ggplot(df.sampled, aes(dmri_rsirnt_scs_pllh, pallidum_lh_res)) +
  geom_point() +
  geom_smooth(method = "loess") + # Add smooth line
  ggtitle("Left Hemisphere") +
  xlab("RNT") +
  ylab("Adjusted Pallidum Volume")

ggsave(paste0(outpath, "/left_pallidum_vol_vs_rnt.png"), p)

rh <- ggplot(df.sampled, aes(dmri_rsirnt_scs_plrh, pallidum_rh_res)) +
  geom_point() +
  geom_smooth(method = "loess") + # Add smooth line
  ggtitle("Right Hemisphere") +
  xlab("RNT") +
  ylab("Adjusted Pallidum Volume")

ggsave(paste0(outpath, "/right_pallidum_vol_vs_rnt.png"), rh)