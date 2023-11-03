# this script is to remove cloud of Landsat images
# no runable code here, just for reference

# ga cloud mask
# 0 = null
# 1 = clear
# 2 = cloud
# 3 = cloud shadow
# 4 = snow
# 5 = water

library(stringr)
library(foreach)
library(doParallel)
source('~/Workspace/RainfallSpectralAnalysis/SpectralAnalysis/function_SetupForGraphics.R')

site_df = read.csv('/g/data/os22/users/yu/lst_project/0_code/study_sites.csv')

for (k in 1:12) {
	
	print(site_df$sitename[k])

	roi = raster(extent(site_df$xmn[k], site_df$xmx[k], site_df$ymn[k], site_df$ymx[k]), res=0.001, crs=PROJ_LATLON)

	path2lst = paste0('/g/data/os22/users/yu/lst_project/0_landsat_data/', site_df$tile[k], '/')
	path2ga = paste('/g/data/xu18/ga/ga_ls8c_ard_3', str_replace(site_df$tile[k], '_', '/'), sep='/') # path to geoscience australia landsat dataset
	out_lst = paste0('/g/data/os22/users/yu/lst_project/0_landsat_masked/', site_df$sitename[k], '/')

	if (!dir.exists(out_lst)) dir.create(out_lst)

	# extract landsat time info
	fl = list.files(path2lst, pattern='LC08_L1TP_')
	yyyy = substr(fl, 18, 21)
	mm = substr(fl, 22, 23)
	dd = substr(fl, 24, 25)

	for (i in 1:length(fl)){

		print(fl[i])

		# path2mask actually is the path to ga landsat dataset
		path2mask = paste(path2ga, yyyy[i], mm[i], dd[i], sep='/')
		mask_rst = list.files(path2mask, pattern='_final_fmask.tif')
		mask_rst = raster(paste0(path2mask, '/', mask_rst))

		lst = raster(paste0(path2lst, fl[i]))

		# GreatWesternWoodlands Site has two tiles, so need to do a mosaic
		if (site_df$sitename[k] == 'GreatWesternWoodlands'){
			
			path2mask1 = str_replace(path2mask, '110/080', '110/081')
			mask_rst1 = list.files(path2mask1, pattern='_final_fmask.tif')
			mask_rst1 = raster(paste0(path2mask1, '/', mask_rst1))

			path2lst1 = str_replace(path2lst, '110_080', '110_081')
			onemoretile = list.files(path2lst1, pattern = paste0('110081_',yyyy[i],mm[i],dd[i]))
			
			# get one more tile
			if (length(onemoretile) > 0){
				
				lst1 = raster(paste0(path2lst1,onemoretile[1]))

				mask_rst = mosaic(mask_rst, mask_rst1, fun='mean')
				lst = mosaic(lst, lst1, fun='mean')
			}
		}
		
		# 2 is cloud and 3 is cloud shadow
		masked_lst = mask(lst, mask_rst, maskvalue=c(2,3))
		masked_lst = projectRaster(masked_lst, roi, method='ngb')

		writeRaster(masked_lst, paste0(out_lst, 'Landsat_LST_cloudrm_',
                                       yyyy[i], mm[i], dd[i], '.tif'),overwrite=TRUE)

	}
}

