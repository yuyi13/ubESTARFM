# this script is about how we calculated metrics for the paper
# no runable code here, just for reference
# the tables indicated in this script are available from the paper
# https://doi.org/10.1016/j.rse.2023.113784

library(stringr)
library(ggsci)

##############################################
# table 4. pre evaluation of modis and landsat
##############################################

library(stringr)
path = '/datasets/work/d61-af-soilmoisture/work/lst_project/metrics/pre_compare'

fl = list.files(path, pattern='.csv$', full.names=TRUE)

sitename = list.files(path, pattern='.csv$')
sitename = str_split(sitename, '_', simplify=TRUE)[,1]

bias_ls_m = c(); bias_ls_l = c(); bias_ls_b = c()
rmse_ls_m = c(); rmse_ls_l = c(); rmse_ls_b = c()
sd_ls_m = c(); sd_ls_l = c(); sd_ls_b = c()

for (k in 1:length(fl)){

    x = read.csv(fl[k])
    xx = na.omit(x)

    # bias
    bias_mod = mean(xx$terra - xx$flux); bias_ls_m = c(bias_ls_m, bias_mod)
    bias_lan = mean(xx$landsat - xx$flux); bias_ls_l = c(bias_ls_l, bias_lan)
    bias_bc = mean(xx$landsat_100m - xx$flux); bias_ls_b = c(bias_ls_b, bias_bc)

    # ubRMSE
    sd_mod = sd(xx$terra - xx$flux); sd_ls_m = c(sd_ls_m, sd_mod)
    sd_lan = sd(xx$landsat - xx$flux); sd_ls_l = c(sd_ls_l, sd_lan)
    sd_bc = sd(xx$landsat_100m - xx$flux); sd_ls_b = c(sd_ls_b, sd_bc)
}

df = data.frame(sitename = sitename, 
                bias_mod = bias_ls_m, bias_lan = bias_ls_l, bias_lan_100m = bias_ls_b,
                ubRMSE_mod = sd_ls_m, ubRMSE_lan = sd_ls_l, ubRMSE_lan_100m = sd_ls_b)

df[,2:ncol(df)] = round(df[,2:ncol(df)] ,2)
write.table(df, '/datasets/work/d61-af-soilmoisture/work/lst_project/metrics/pre_metrics.csv', quote=FALSE, sep=',', row.names=FALSE)

###########################################
# table 5. algorithm performance evaluation
###########################################

library(stringr)

site_df = read.csv('/g/data/os22/users/yu/lst_project/0_code/study_sites.csv')

bias_ls_l = c(); bias_ls_b = c()
sd_ls_l = c(); sd_ls_b = c()
train_number = c()
sample_number = c()

for (k in 1:12){

    terra_lst = paste0('/g/data/os22/users/yu/lst_project/0_terra_data/', site_df$sitename[k], '/')
    landsat_lst = paste0('/g/data/os22/users/yu/lst_project/0_landsat_masked/', site_df$sitename[k], '/')

    # filter out the landsat data for training only
    train_lsat = list.files(landsat_lst, pattern='Landsat_LST_', full.names=TRUE)
    lsat_dates = as.POSIXct(paste0(str_sub(train_lsat, -12, -5), ' 10:30'), format='%Y%m%d %H:%M', tz=site_df$timezone[k])
    attr(lsat_dates, 'tzone') = 'GMT'
    valid_idx = which(file.info(train_lsat)$size > 2.5e6)

    # further check if modis files satisfies size requirements
    mdis_files = paste0(terra_lst, 'MOD11A1_LST_daytime_', str_sub(train_lsat, -12, -5)[valid_idx], '.tif')
    valid_idx = valid_idx[which(file.info(mdis_files)$size > 1e6)]

    train_number = c(train_number, length(valid_idx))

    metr_path = '/g/data/os22/users/yu/lst_project/metrics/estarfm/'

    x = read.csv(paste0(metr_path, site_df$sitename[k], '_validation.csv'))

    t_idx = c()

    for (t in 1:length(lsat_dates[valid_idx])){

        whichT = which(as.POSIXct(x$time_GMT, tz='GMT') < lsat_dates[valid_idx][t] + 7200 & 
                        as.POSIXct(x$time_GMT, tz='GMT') > lsat_dates[valid_idx][t] - 7200)
        t_idx = c(t_idx, whichT)
    }

    # x[t_idx] is the dataframe on the training dates

    xx = na.omit(x[t_idx,-2])
    sample_number = c(sample_number, nrow(xx))

    bias_lan = mean(xx$estarfm_100m - xx$terra); bias_ls_l = c(bias_ls_l, bias_lan)
    bias_bc = mean(xx$estarfm_bc_100m - xx$terra); bias_ls_b = c(bias_ls_b, bias_bc)

    sd_lan = sd(xx$estarfm_100m - xx$terra); sd_ls_l = c(sd_ls_l, sd_lan)
    sd_bc = sd(xx$estarfm_bc_100m - xx$terra); sd_ls_b = c(sd_ls_b, sd_bc)
}


df = data.frame(sitename = site_df$sitename[1:12], bias_estarfm_100m = bias_ls_l, bias_ubestarfm_100m = bias_ls_b,
                ubRMSE_estarfm_100m = sd_ls_l, ubRMSE_ubestarfm_100m = sd_ls_b, train_number = train_number,
                sample_number = sample_number)

df[,2:5] = round(df[,2:5] ,2)
df = df[order(df$sitename),]
write.table(df, '/g/data/os22/users/yu/lst_project/metrics/mod_pefm_metrics.csv', quote=FALSE, sep=',', row.names=FALSE)


########################################
# table 6. overall evaluation of results
########################################

library(stringr)

site_df = read.csv('/datasets/work/d61-af-soilmoisture/work/lst_project/0_code/study_sites.csv')

bias_ls_m = c(); bias_ls_l = c(); bias_ls_b = c()
sd_ls_m = c(); sd_ls_l = c(); sd_ls_b = c()
cor_ls_m = c(); cor_ls_l = c(); cor_ls_b = c()
sample_number = c()

for (k in 1:12){

    terra_lst = paste0('/datasets/work/d61-af-soilmoisture/work/lst_project/0_terra_data/', site_df$sitename[k], '/')
    landsat_lst = paste0('/datasets/work/d61-af-soilmoisture/work/lst_project/0_landsat_masked/', site_df$sitename[k], '/')

    # filter out the landsat data for training
    # in the next step we will exclude these data
    train_lsat = list.files(landsat_lst, pattern='Landsat_LST_', full.names=TRUE)
    lsat_dates = as.POSIXct(paste0(str_sub(train_lsat, -12, -5), ' 10:30'), format='%Y%m%d %H:%M', tz=site_df$timezone[k])
    attr(lsat_dates, 'tzone') = 'GMT'
    valid_idx = which(file.info(train_lsat)$size > 2.5e6)

    # further check if modis files satisfies size requirements
    mdis_files = paste0(terra_lst, 'MOD11A1_LST_daytime_', str_sub(train_lsat, -12, -5)[valid_idx], '.tif')
    valid_idx = valid_idx[which(file.info(mdis_files)$size > 1e6)]

    metr_path = '/datasets/work/d61-af-soilmoisture/work/lst_project/metrics/estarfm/'

    # the validation data, the time step of which is in GMT format
    x = read.csv(paste0(metr_path, site_df$sitename[k], '_validation.csv'))

    t_idx = c()

    for (t in 1:length(lsat_dates[valid_idx])){

        whichT = which(as.POSIXct(x$time_GMT, tz='GMT') < lsat_dates[valid_idx][t] + 7200 & 
                        as.POSIXct(x$time_GMT, tz='GMT') > lsat_dates[valid_idx][t] - 7200)
        t_idx = c(t_idx, whichT)
    }

    # x[t_idx] is the dataframe on the training dates
    # exclude those data used for training 
    xx = na.omit(x[-t_idx,])
    whichT = which(abs(xx$terra - xx$flux) > 10 | abs(xx$estarfm_bc_100m - xx$flux) > 10 ) # exclude outliers

    if (length(whichT) > 0){
      xx = xx[-whichT,]
    }
    
    sample_number = c(sample_number, nrow(xx))

    bias_mod = mean(xx$terra - xx$flux); bias_ls_m = c(bias_ls_m, bias_mod)
    bias_lan = mean(xx$estarfm_100m - xx$flux); bias_ls_l = c(bias_ls_l, bias_lan)
    bias_bc = mean(xx$estarfm_bc_100m - xx$flux); bias_ls_b = c(bias_ls_b, bias_bc)

    sd_mod = sd(xx$terra - xx$flux); sd_ls_m = c(sd_ls_m, sd_mod)
    sd_lan = sd(xx$estarfm_100m - xx$flux); sd_ls_l = c(sd_ls_l, sd_lan)
    sd_bc = sd(xx$estarfm_bc_100m - xx$flux); sd_ls_b = c(sd_ls_b, sd_bc)

    cor_mod = cor(xx$terra, xx$flux); cor_ls_m = c(cor_ls_m, cor_mod)
    cor_lan = cor(xx$estarfm_100m, xx$flux); cor_ls_l = c(cor_ls_l, cor_lan)
    cor_bc = cor(xx$estarfm_bc_100m, xx$flux); cor_ls_b = c(cor_ls_b, cor_bc)

}

df = data.frame(sitename = site_df$sitename[1:12], 
                bias_mod = bias_ls_m, bias_estarfm_100m = bias_ls_l, bias_ubestarfm_100m = bias_ls_b,
                ubRMSE_mod = sd_ls_m, ubRMSE_estarfm_100m = sd_ls_l, ubRMSE_ubestarfm_100m = sd_ls_b, 
                cor_mod = cor_ls_m, cor_estarfm_100m = cor_ls_l, cor_ubestarfm_100m = cor_ls_b,
                sample_number = sample_number)

df[,2:10] = round(df[,2:10] ,2)
df = df[order(df$sitename),]

sum(df$bias_mod * df$sample_number) / sum(df$sample_number)
sum(df$bias_ubestarfm_100m * df$sample_number) / sum(df$sample_number)

sum(df$ubRMSE_mod * df$sample_number) / sum(df$sample_number)
sum(df$ubRMSE_ubestarfm_100m * df$sample_number) / sum(df$sample_number)

write.table(df, '/g/data/os22/users/yu/lst_project/metrics/estarfm_metrics.csv', quote=FALSE, sep=',', row.names=FALSE)
