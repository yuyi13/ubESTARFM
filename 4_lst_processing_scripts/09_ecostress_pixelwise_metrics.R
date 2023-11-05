# this script is to derive pixelwise metrics against ecostress scenes
# no runable code here, just for reference

# packages
source('~/Workspace/RainfallSpectralAnalysis/SpectralAnalysis/function_SetupForGraphics.R')
# more details here: https://github.com/LuigiJR/SpectralAnalysis
library(stringr)

site_df = read.csv('/datasets/work/d61-af-soilmoisture/work/users/yu/code/study_sites.csv')

# create empty lists for all sites
bias_ls_l = c(); bias_ls_b = c()
sd_ls_l = c(); sd_ls_b = c()
cor_ls_l = c(); cor_ls_b = c()
scene_number = c()

for (i in 1:12){

    print(site_df$sitename[i])

    path2terra = paste0('/datasets/work/d61-af-soilmoisture/work/lst_project/0_terra_data/', site_df$sitename[i], '/')
    path2esfm = paste0('/datasets/work/d61-af-soilmoisture/work/lst_project/2_estarfm_output/', site_df$sitename[i], '/')
    path2esfmbc = paste0('/datasets/work/d61-af-soilmoisture/work/lst_project/2_estarfm_output_bc/', site_df$sitename[i], '/')
    path2eco = paste0('/datasets/work/d61-af-soilmoisture/work/ECOSTRESS/flux_subsets/', site_df$sitename[i], '/')

    fl = list.files(path2eco, pattern='.tif$')

    # create empty lists for this site
    bias_this_scene_l = c(); bias_this_scene_b = c()
    sd_this_scene_l = c(); sd_this_scene_b = c()
    cor_this_scene_l = c(); cor_this_scene_b = c()

    for (k in 1:length(fl)){
        
        # read rasters
        eco = raster(paste0(path2eco, fl[k]))
	    date = str_sub(fl[k], 15, 22)
	    esfm = raster(paste0(path2esfm, 'ESTARFM_LST_daytime_', date, '.tif'))
	    bc = raster(paste0(path2esfmbc, 'ESTARFM_LST_daytime_', date, '.tif'))
	    mod = raster(paste0(path2terra, 'MOD11A1_LST_daytime_', date, '.tif'))

        # convert to data frame
	    df = data.frame(x = as.vector(as.matrix(esfm)), y = as.vector(as.matrix(bc)), z = as.vector(as.matrix(eco)))
	    df = na.omit(df)

        # if near clear-sky pixels are more than a half
	    if (nrow(df) > 5e5){
            
            # lan here means estarfm (because it relies on landsat)
            # bc here means ubestarfm
            bias_lan = mean(df$x - df$z); bias_this_scene_l = c(bias_this_scene_l, bias_lan)
            bias_bc = mean(df$y - df$z); bias_this_scene_b = c(bias_this_scene_b, bias_bc)

            sd_lan = sd(df$x - df$z); sd_this_scene_l = c(sd_this_scene_l, sd_lan)
            sd_bc = sd(df$y - df$z); sd_this_scene_b = c(sd_this_scene_b, sd_bc)

            cor_lan = cor(df$x, df$z); cor_this_scene_l = c(cor_this_scene_l, cor_lan)
            cor_bc = cor(df$y, df$z); cor_this_scene_b = c(cor_this_scene_b, cor_bc)
        }
    }

    # count the number of scenes
    scene_number = c(scene_number, length(bias_this_scene_l))
    
    # get average metrics for this site
    bias_ls_l = c(bias_ls_l, mean(bias_this_scene_l)); bias_ls_b = c(bias_ls_b, mean(bias_this_scene_b))
    sd_ls_l = c(sd_ls_l, mean(sd_this_scene_l)); sd_ls_b = c(sd_ls_b, mean(sd_this_scene_b))
    cor_ls_l = c(cor_ls_l, mean(cor_this_scene_l)); cor_ls_b = c(cor_ls_b, mean(cor_this_scene_b))

}

# create the global data frame
df = data.frame(sitename = site_df$sitename[1:12], bias_estarfm_100m = bias_ls_l, bias_ubestarfm_100m = bias_ls_b,
                ubRMSE_estarfm_100m = sd_ls_l, ubRMSE_ubestarfm_100m = sd_ls_b, 
                cor_estarfm_100m = cor_ls_l, cor_ubestarfm_100m = cor_ls_b,
                scene_number = scene_number)

# round the metrics
df[,2:7] = round(df[,2:7] ,2)
df = df[order(df$sitename),]

write.table(df, '/g/data/os22/users/yu/lst_project/metrics/ecostress_metrics.csv', quote=FALSE, sep=',', row.names=FALSE)
