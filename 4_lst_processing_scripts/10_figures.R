# this script is about how we make plots for the paper
# no runable code here, just for reference
# the figures indicated in this script are available from the paper
# https://doi.org/10.1016/j.rse.2023.113784

######################
# fig.4 pre-evaluation
######################

library(stringr)
library(ggsci)

SCICOL = pal_nejm('default', alpha=0.8)(4)

path = '/g/data/os22/users/yu/lst_project/metrics/pre_compare'
fl  = list.files(path, pattern='.csv', full.names=TRUE)

df = read.csv('/g/data/os22/users/yu/lst_project/metrics/pre_metrics.csv')

png('/g/data/os22/users/yu/lst_project/figures/pre_scatterplot.png', width=1200, height=800)

m = rbind(c(1,2,3,4), c(5,6,7,8), c(9,10,11,12))
layout(m); par(mar = c(2, 2.5, 2, 0.5))

for (k in 1:length(fl)){
    
    x = read.csv(fl[k])
    xx = na.omit(x)

    plot(xx$flux, xx$terra, xlim=c(270,330), ylim=c(270,330), col=SCICOL[1], pch=17, 
         xlab = 'In-situ LST (K)', ylab='Remotely sensed LST (K)', cex=1.5, cex.axis=1.5)
    points(xx$flux, xx$landsat_100m, xlim=c(270,330), ylim=c(270,330), col=SCICOL[2], cex=1.5, pch=19, xlab = '', ylab='')
    abline(a=0, b=1)
    sitename = str_split(fl[k], '/')[[1]][10]; sitename = str_split(sitename, '_')[[1]][1]

    # an exception for legend of Great Western Woodlands
    if (k == 5) {

        legend(x=265, y=333, legend=c(paste0('(', letters[k], ') ', sitename), 
                paste0(' (Sample number: ', nrow(xx), ')')), bty = 'n', cex=1.5)
        legend('bottomright', legend=c('Landsat', paste0('Bias = ', format(df[k,4], nsmall=2)), paste0('ubRMSE = ', format(df[k,10], nsmall=2))), 
                cex=1.5, col = c(SCICOL[2], NA, NA), pch = c(19, NA, NA))
        legend(x=270, y=322, legend=c('MODIS', paste0('Bias = ', format(df[k,2], nsmall=2)), paste0('ubRMSE = ', format(df[k,8], nsmall=2))), 
            cex=1.5, col = c(SCICOL[1], NA, NA), pch = c(17, NA, NA))

    } else {

        legend(x=265, y=333, legend=paste0('(', letters[k], ') ', sitename, ' (Sample number: ', nrow(xx), ')'), bty = 'n', cex=1.5)
        legend('bottomright', legend=c('Landsat', paste0('Bias = ', format(df[k,4], nsmall=2)), paste0('ubRMSE = ', format(df[k,10], nsmall=2))), 
                cex=1.5, col = c(SCICOL[2], NA, NA), pch = c(19, NA, NA))
        legend(x=270, y=326, legend=c('MODIS', paste0('Bias = ', format(df[k,2], nsmall=2)), paste0('ubRMSE = ', format(df[k,8], nsmall=2))), 
                cex=1.5, col = c(SCICOL[1], NA, NA), pch = c(17, NA, NA))
    }
}

dev.off()

############################
# fig.5 training performance
############################

library(ggsci)
library(stringr)

SCICOL = pal_nejm('default', alpha=0.8)(4)

df = read.csv('/datasets/work/d61-af-soilmoisture/work/lst_project/metrics/mod_pefm_metrics.csv')
df = df[order(df$sitename),]

metr_path = '/datasets/work/d61-af-soilmoisture/work/lst_project/metrics/estarfm/'
site_df = read.csv('/datasets/work/d61-af-soilmoisture/work/users/yu/code/study_sites.csv')[1:12,]
site_df = site_df[order(site_df$sitename),]

png('/datasets/work/d61-af-soilmoisture/work/lst_project/figures/model_pefm_scatterplot.png', width=1200, height=800)

m = rbind(c(1,2,3,4), c(5,6,7,8), c(9,10,11,12))
layout(m); par(mar = c(2, 2.5, 2, 0.5))

for (k in 1:12){

    x = read.csv(paste0(metr_path, site_df$sitename[k], '_validation.csv'))
    xx = na.omit(x)

    plot(xx$flux, xx$estarfm_bc_100m, xlim=c(270,330), ylim=c(270,330), col=SCICOL[3], pch=17, 
         xlab = 'In-situ LST (K)', ylab='Remotely sensed LST (K)', cex=1.5, cex.axis=1.5)
    points(xx$flux, xx$estarfm_100m, xlim=c(270,330), ylim=c(270,330), col=SCICOL[4], cex=1.5, pch=19, xlab = '', ylab='')
    abline(a=0, b=1)
    sitename = site_df$sitename[k]

    if (k == 5) {
        
    legend(x=265, y=333, legend=c(paste0('(', letters[k], ') ', sitename), 
            paste0(' (Sample number: ', nrow(xx), ')')), bty = 'n', cex=1.5)

    legend('bottomright', legend=c('ESTARFM', paste0('Bias = ', format(df[k,2], nsmall=2)), paste0('ubRMSE = ', format(df[k,4], nsmall=2))), 
            cex=1.5, col = c(SCICOL[4], NA, NA), pch = c(19, NA, NA))
    legend(x=270, y=322, legend=c('ubESTARFM', paste0('Bias = ', format(df[k,3], nsmall=2)), paste0('ubRMSE = ', format(df[k,5], nsmall=2))), 
           cex=1.5, col = c(SCICOL[3], NA, NA), pch = c(17, NA, NA))

    } else {

    legend(x=265, y=333, legend=paste0('(', letters[k], ') ', sitename, ' (Sample number: ', nrow(xx), ')'), bty = 'n', cex=1.5)
    legend('bottomright', legend=c('ESTARFM', paste0('Bias = ', format(df[k,2], nsmall=2)), paste0('ubRMSE = ', format(df[k,4], nsmall=2))), 
            cex=1.5, col = c(SCICOL[4], NA, NA), pch = c(19, NA, NA))
    legend(x=270, y=326, legend=c('ubESTARFM', paste0('Bias = ', format(df[k,3], nsmall=2)), paste0('ubRMSE = ', format(df[k,5], nsmall=2))), 
           cex=1.5, col = c(SCICOL[3], NA, NA), pch = c(17, NA, NA))
    }
}

dev.off()


###########################################
# fig. 6 spatial examples on training dates
###########################################

library(fields)
source('~/Workspace/RainfallSpectralAnalysis/SpectralAnalysis/function_SetupForGraphics.R')

site_df = read.csv('/g/data/os22/users/yu/lst_project/0_code/study_sites.csv')

png('/g/data/os22/users/yu/lst_project/figures/spatial_example.png', width = 1600, height = 1600)
m = rbind(c(1,2,3,4), c(5,6,7,8), c(9,10,11,12), c(13,14,15,16))
layout(m); par(mar = c(0.5,0.5,0.5,0.5))

for (f in c(1,3)){

    path2esfm = paste0('/scratch/os22/yy4778/lst_project/2_estarfm_output/', site_df$sitename[f], '/')
    path2esfmbc = paste0('/scratch/os22/yy4778/lst_project/2_estarfm_output_bc/', site_df$sitename[f], '/')
    path2mod = paste0('/g/data/os22/users/yu/lst_project/0_terra_data/', site_df$sitename[f], '/')
    path2lan = paste0('/g/data/os22/users/yu/lst_project/0_landsat_masked/', site_df$sitename[f], '/')

    max_val = 325
    min_val = 290

    # Calperum; my initial sequence of sites are not alphabetically ordered
    if (f == 1){

        DOI = as.Date('2014-11-10')
        ex = extent(140.5,140.6,-34.05,-33.95)

        r1 = raster(paste0(path2esfm, format(DOI, 'ESTARFM_LST_daytime_%Y%m%d.tif'))); r1[r1 > max_val] = max_val
        r2 = raster(paste0(path2esfmbc, format(DOI, 'ESTARFM_LST_daytime_%Y%m%d.tif')))
        r3 = raster(paste0(path2mod, format(DOI, 'MOD11A1_LST_daytime_%Y%m%d.tif')))
        r4 = raster(paste0(path2lan, format(DOI, 'Landsat_LST_cloudrm_%Y%m%d.tif'))); r4[r4 > max_val] = max_val

        image(r3, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'MODIS', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(a)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        image(r4, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'Landsat', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(b)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        image(r1, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(c)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        image(r2, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ubESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(d)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        # zoomed area
        r1 = crop(r1, ex); r2 = crop(r2, ex); r3 = crop(r3, ex); r4 = crop(r4, ex)

        image(r3, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'MODIS', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(e)', bty='n', cex=6)
        image(r4, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'Landsat', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(f)', bty='n', cex=6)
        image(r1, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(g)', bty='n', cex=6)
        image(r2, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ubESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(h)', bty='n', cex=6)

    # Samford; my initial sequence of sites are not alphabetically ordered
    } else {
    
        DOI = as.Date('2016-11-29')
        ex = extent(152.8,152.9,-27.4,-27.3)

        r1 = raster(paste0(path2esfm, format(DOI, 'ESTARFM_LST_daytime_%Y%m%d.tif'))); r1[r1 > max_val] = max_val
        r2 = raster(paste0(path2esfmbc, format(DOI, 'ESTARFM_LST_daytime_%Y%m%d.tif')))
        r3 = raster(paste0(path2mod, format(DOI, 'MOD11A1_LST_daytime_%Y%m%d.tif')))
        r4 = raster(paste0(path2lan, format(DOI, 'Landsat_LST_cloudrm_%Y%m%d.tif'))); r4[r4 > max_val] = max_val

        image(r3, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'MODIS', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(i)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        image(r4, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'Landsat', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(j)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        image(r1, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(k)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        image(r2, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ubESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(l)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        # zoomed area
        r1 = crop(r1, ex); r2 = crop(r2, ex); r3 = crop(r3, ex); r4 = crop(r4, ex)

        image(r3, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'MODIS', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(m)', bty='n', cex=6)
        image(r4, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'Landsat', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(n)', bty='n', cex=6)
        image(r1, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(o)', bty='n', cex=6)
        image(r2, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ubESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(p)', bty='n', cex=6)
    }
}

dev.off()

##########################
# fig.7 overall evaluation
##########################

library(stringr)
library(ggsci)

SCICOL = pal_nejm('default', alpha=0.8)(4)

df = read.csv('/g/data/os22/users/yu/lst_project/metrics/estarfm_metrics.csv')
df = df[order(df$sitename),]
df[,2:7] = round(df[,2:7], 2)

site_df = read.csv('/g/data/os22/users/yu/lst_project/0_code/study_sites.csv')[1:12,]
site_df = site_df[order(site_df$sitename),]

png('/g/data/os22/users/yu/lst_project/figures/overall_scatterplot.png', width=1200, height=800)

m = rbind(c(1,2,3,4), c(5,6,7,8), c(9,10,11,12))
layout(m); par(mar = c(2, 2.5, 2, 0.5))

for (k in 1:12){
    
    terra_lst = paste0('/g/data/os22/users/yu/lst_project/0_terra_data/', site_df$sitename[k], '/')
    landsat_lst = paste0('/g/data/os22/users/yu/lst_project/0_landsat_masked/', site_df$sitename[k], '/')

    # landsat data for training, which we will exclude later
    train_lsat = list.files(landsat_lst, pattern='Landsat_LST_', full.names=T)
    lsat_dates = as.POSIXct(paste0(str_sub(train_lsat, -12, -5), ' 10:30'), format='%Y%m%d %H:%M', tz=site_df$timezone[k])
    attr(lsat_dates, 'tzone') = 'GMT'
    valid_idx = which(file.info(train_lsat)$size > 2.5e6)

    # further check if modis files satisfies size requirements
    mdis_files = paste0(terra_lst, 'MOD11A1_LST_daytime_', str_sub(train_lsat, -12, -5)[valid_idx], '.tif')
    valid_idx = valid_idx[which(file.info(mdis_files)$size > 1e6)]

    metr_path = '/g/data/os22/users/yu/lst_project/metrics/estarfm/'

    x = read.csv(paste0(metr_path, site_df$sitename[k], '_validation.csv'))

    t_idx = c()

    for (t in 1:length(lsat_dates[valid_idx])){

        whichT = which(as.POSIXct(x$time_GMT, tz='GMT') < lsat_dates[valid_idx][t] + 7200 & 
                        as.POSIXct(x$time_GMT, tz='GMT') > lsat_dates[valid_idx][t] - 7200)
        t_idx = c(t_idx, whichT)
    }

    # x[t_idx] is the dataframe on the training dates

    xx = na.omit(x[-t_idx,])
    whichT = which(abs(xx$terra - xx$flux) > 10 | abs(xx$estarfm_bc_100m - xx$flux) > 10 )  # find out the abnormal values

    if (length(whichT) > 0){
      xx = xx[-whichT,]
    }

    plot(xx$flux, xx$terra, xlim=c(270,330), ylim=c(270,330), col=SCICOL[1], pch=17, 
         xlab = 'In-situ LST (K)', ylab='Remotely sensed LST (K)', cex.lab=1.5, cex.axis=1.5)
    points(xx$flux, xx$estarfm_100m, xlim=c(270,330), ylim=c(270,330), col=SCICOL[4], pch=19, xlab = '', ylab='')
    points(xx$flux, xx$estarfm_bc_100m, xlim=c(270,330), ylim=c(270,330), col=SCICOL[3], pch=17, xlab = '', ylab='')
    abline(a=0, b=1)
    sitename = site_df$sitename[k]

    if (k == 5) {
        
        legend(x=265, y=333, legend=c(paste0('(', letters[k], ') ', sitename), 
                paste0(' (Sample number: ', nrow(xx), ')')), bty = 'n', cex=1.5)

        legend('bottomright', legend=c('ESTARFM', paste0('Bias = ', format(df[k,3],nsmall=2)), paste0('ubRMSE = ', format(df[k,6],nsmall=2)),
                'ubESTARFM', paste0('Bias = ', format(df[k,4],nsmall=2)), paste0('ubRMSE = ', format(df[k,7],nsmall=2))), 
                cex=1.5, col = c(SCICOL[4], NA, NA, SCICOL[3], NA, NA), pch = c(19, NA, NA, 17, NA, NA))
        legend(x=270, y=322, legend=c('MODIS', paste0('Bias = ', format(df[k,2],nsmall=2)),  paste0('ubRMSE = ', format(df[k,5],nsmall=2))), 
                cex=1.5, col = c(SCICOL[1], NA, NA), pch = c(17, NA, NA))

    } else {

        legend(x=265, y=333, legend=paste0('(', letters[k], ') ', sitename, ' (Sample number: ', nrow(xx), ')'), bty = 'n', cex=1.5)
        legend('bottomright', legend=c('ESTARFM', paste0('Bias = ', format(df[k,3],nsmall=2)), paste0('ubRMSE = ', format(df[k,6],nsmall=2)),
                'ubESTARFM', paste0('Bias = ', format(df[k,4],nsmall=2)), paste0('ubRMSE = ', format(df[k,7],nsmall=2))), 
                cex=1.5, col = c(SCICOL[4], NA, NA, SCICOL[3], NA, NA), pch = c(19, NA, NA, 17, NA, NA))
        legend(x=270, y=326, legend=c('MODIS', paste0('Bias = ', format(df[k,2],nsmall=2)),  paste0('ubRMSE = ', format(df[k,5],nsmall=2))), 
                cex=1.5, col = c(SCICOL[1], NA, NA), pch = c(17, NA, NA))
    }
}

dev.off()

###################
# fig.9 time series
###################

library(stringr)
library(ggsci)

SCICOL = pal_nejm('default', alpha=0.8)(4)

DOI = seq(as.Date('2013-01-01'), as.Date('2021-12-31'), by='day')

df = read.csv('/g/data/os22/users/yu/lst_project/metrics/estarfm_metrics.csv')
df = df[order(df$sitename),]
df[,2:7] = round(df[,2:7], 2)

site_df = read.csv('/g/data/os22/users/yu/lst_project/0_code/study_sites.csv')[1:12,]
site_df = site_df[order(site_df$sitename),]

png('/g/data/os22/users/yu/lst_project/figures/time_series.png', width=1800, height=800)

m = rbind(c(1,2,3,4), c(5,6,7,8), c(9,10,11,12))
layout(m); par(mar = c(2, 2.5, 2, 2))

metr_path = '/g/data/os22/users/yu/lst_project/metrics/estarfm/'

for (k in 1:12){
    
    x = read.csv(paste0(metr_path, site_df$sitename[k], '_validation.csv'))

    whichT = which(abs(x$terra - x$flux) > 10 | abs(x$estarfm_bc_100m - x$flux) > 10 | 
                x$flux > 350)  # find out the abnormal values

    if (length(whichT) > 0){
      x = x[-whichT,]
    }

    df = na.omit(x)

    plot(as.POSIXct(df$time_GMT), df$flux, type='p', xlab = 'Time in GMT', ylab = 'LST (K)', pch=19,
         ylim=c(270,330), cex.axis=2)

    points(as.POSIXct(df$time_GMT), df$terra, col=SCICOL[1], pch=17)
    points(as.POSIXct(df$time_GMT), df$estarfm_100m, col=SCICOL[4], pch=19)
    points(as.POSIXct(df$time_GMT), df$estarfm_bc_100m, col=SCICOL[3], pch=17)

    if (k ==2 | k == 5 | k == 6 | k == 12){
        legend('bottomleft', paste0('(', letters[k], ') ', site_df$sitename[k]), cex=2.5, bty='n')
    } else {
        legend('topleft', paste0('(', letters[k], ') ', site_df$sitename[k]), cex=2.5, bty='n')
    }
    #legend('bottomright', legend = c('Flux', 'MODIS', 'ESTARFM', 'ubESTARFM'), fill= c('black', SCICOL[c(1,4,3)]))

}

dev.off()
xlim=c(1,6), ylim=c(1,3)

# make a legend
png('/g/data/os22/users/yu/lst_project/figures/legend.png', width=800, height=800)

plot(1:10, 1:10, col='white')
legend(1,2, legend = 'Flux', pch=19, col='black', bty='n', cex=2)
legend(2.5,2, legend = 'MODIS', pch=17, col=SCICOL[1], bty='n', cex=2)
legend(4.4,2, legend = 'ESTARFM', pch=19, col=SCICOL[4], bty='n', cex=2)
legend(6.8,2, legend = 'ubESTARFM', pch=17, col=SCICOL[3], bty='n', cex=2)

dev.off()

##############################
# fig.s1 time series sub-plots
##############################

library(ggsci)
SCICOL = pal_nejm('default', alpha=0.8)(4)
DOI = seq(as.Date('2013-01-01'), as.Date('2021-12-31'), by='day')

library(stringr)
df = read.csv('/g/data/os22/users/yu/lst_project/metrics/estarfm_metrics.csv')
df = df[order(df$sitename),]
df[,2:7] = round(df[,2:7], 2)

site_df = read.csv('/g/data/os22/users/yu/lst_project/0_code/study_sites.csv')[1:12,]
site_df = site_df[order(site_df$sitename),]

png('/g/data/os22/users/yu/lst_project/figures/time_series_shortened.png', width=1800, height=800)

m = rbind(c(1,2,3,4), c(5,6,7,8), c(9,10,11,12))
layout(m); par(mar = c(2, 2.5, 2, 2))

metr_path = '/g/data/os22/users/yu/lst_project/metrics/estarfm/'

for (k in 1:12){
    
    x = read.csv(paste0(metr_path, site_df$sitename[k], '_validation.csv'))

    # exclude the abnormal values
    whichT = which(abs(x$terra - x$flux) > 10 | abs(x$estarfm_bc_100m - x$flux) > 10 | x$flux > 350)  
    if (length(whichT) > 0) x = x[-whichT,]

    df = x

    if (k != 7){

        plot(as.POSIXct(df$time_GMT), df$flux, type='l', xlim=c(as.POSIXct('2018-01-01'), as.POSIXct('2018-12-31')),
            xlab = 'Time in GMT', ylab = 'LST (K)', pch=19, lwd=2,
            ylim=c(270,330), cex.axis=2)

    } else { # because Samford (site 7) has no data for 2018

        plot(as.POSIXct(df$time_GMT), df$flux, type='l', xlim=c(as.POSIXct('2017-01-01'), as.POSIXct('2017-12-31')),
            xlab = 'Time in GMT', ylab = 'LST (K)', pch=19, lwd=2,
            ylim=c(270,330), cex.axis=2)
        
    }

    points(as.POSIXct(df$time_GMT), df$flux, pch=17, lwd=2)
    points(as.POSIXct(df$time_GMT), df$terra, col=SCICOL[1], pch=17, lwd=2)
    points(as.POSIXct(df$time_GMT), df$estarfm_100m, col=SCICOL[4], pch=19, lwd=2)
    points(as.POSIXct(df$time_GMT), df$estarfm_bc_100m, col=SCICOL[3], pch=17, lwd=2)

    lines(as.POSIXct(df$time_GMT), df$terra, col=SCICOL[1], pch=17, lwd=2)
    lines(as.POSIXct(df$time_GMT), df$estarfm_100m, col=SCICOL[4], pch=19, lwd=2)
    lines(as.POSIXct(df$time_GMT), df$estarfm_bc_100m, col=SCICOL[3], pch=17, lwd=2)
    
    if (k ==2 | k == 5 | k == 6 | k == 12){
        legend('bottomleft', paste0('(', letters[k], ') ', site_df$sitename[k]), cex=2.5, bty='n')
    } else {
        legend('topleft', paste0('(', letters[k], ') ', site_df$sitename[k]), cex=2.5, bty='n')
    }
}

dev.off()

###############################
# fig. 10 ecostress comparisons
###############################

source('~/Workspace/RainfallSpectralAnalysis/SpectralAnalysis/function_SetupForGraphics.R')

site_df = read.csv('/g/data/os22/users/yu/lst_project/0_code/study_sites.csv')

png('/g/data/os22/users/yu/lst_project/figures/spatial_ecostress.png', width = 1600, height = 1600)
m = rbind(c(1,2,3,4), c(5,6,7,8), c(9,10,11,12), c(13,14,15,16))
layout(m); par(mar = c(0.5,0.5,0.5,0.5))

for (f in c(12,4)){

    path2esfm = paste0('/scratch/os22/yy4778/lst_project/2_estarfm_output/', site_df$sitename[f], '/')
    path2esfmbc = paste0('/scratch/os22/yy4778/lst_project/2_estarfm_output_bc/', site_df$sitename[f], '/')
    path2mod = paste0('/g/data/os22/users/yu/lst_project/0_terra_data/', site_df$sitename[f], '/')
    path2eco = paste0('/g/data/os22/users/yu/lst_project/ecostress_images/', site_df$sitename[f], '/')

    max_val = 305
    min_val = 270

    # Great Western Woodlands; my initial sequence of sites are not alphabetically ordered
    if (f == 12){

        DOI = as.Date('2019-07-11')
        ex = extent(120.6,120.7,-30.25,-30.15)

        r1 = raster(paste0(path2esfm, format(DOI, 'ESTARFM_LST_daytime_%Y%m%d.tif')))
        r2 = raster(paste0(path2esfmbc, format(DOI, 'ESTARFM_LST_daytime_%Y%m%d.tif')))
        r3 = raster(paste0(path2mod, format(DOI, 'MOD11A1_LST_daytime_%Y%m%d.tif')))
        r4 = raster(paste0(path2eco, format(DOI, 'ECOSTRESS_LST_20190711T023821.tif')))

        image(r3, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'MODIS', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(a)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        image(r4, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'Landsat', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(b)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        image(r1, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(c)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        image(r2, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ubESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(d)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        # zoomed area
        r1 = crop(r1, ex); r2 = crop(r2, ex); r3 = crop(r3, ex); r4 = crop(r4, ex)

        image(r3, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'MODIS', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(e)', bty='n', cex=6)
        image(r4, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'Landsat', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(f)', bty='n', cex=6)
        image(r1, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(g)', bty='n', cex=6)
        image(r2, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ubESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(h)', bty='n', cex=6)

    # Tumbarumba; my initial sequence of sites are not alphabetically ordered
    } else {
    
        DOI = as.Date('2020-09-10')
        ex = extent(148.1,148.2,-35.7,-35.6)

        r1 = raster(paste0(path2esfm, format(DOI, 'ESTARFM_LST_daytime_%Y%m%d.tif')))
        r2 = raster(paste0(path2esfmbc, format(DOI, 'ESTARFM_LST_daytime_%Y%m%d.tif')))
        r3 = raster(paste0(path2mod, format(DOI, 'MOD11A1_LST_daytime_%Y%m%d.tif')))
        r4 = raster(paste0(path2eco, format(DOI, 'ECOSTRESS_LST_20200910T001439.tif')))

        image(r3, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'MODIS', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(i)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        image(r4, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'Landsat', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(j)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        image(r1, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(k)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        image(r2, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ubESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        rect(ex[1],ex[3],ex[2],ex[4], lwd = 3); legend(site_df$xmn[f] - 0.2, site_df$ymx[f], '(l)', bty='n', cex=6)
        lines(c(site_df$xmn[f], ex[1]), c(site_df$ymn[f], ex[3]), lty=2); lines(c(site_df$xmx[f], ex[2]), c(site_df$ymn[f], ex[3]), lty=2)

        # zoomed area
        r1 = crop(r1, ex); r2 = crop(r2, ex); r3 = crop(r3, ex); r4 = crop(r4, ex)

        image(r3, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'MODIS', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(m)', bty='n', cex=6)
        image(r4, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'Landsat', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(n)', bty='n', cex=6)
        image(r1, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(o)', bty='n', cex=6)
        image(r2, col=TemperatureRamp(64), zlim=c(min_val, max_val), xaxt='n', yaxt='n', xlab = 'ubESTARFM', ylab = NA, cex.axis=1.5, cex.lab = 1.5)
        legend(ex[1] - 0.02, ex[4], '(p)', bty='n', cex=6)
    }
}

dev.off()
