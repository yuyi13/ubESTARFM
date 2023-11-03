# this script is to extract LST time series for evaluations
# no runable code here, just for reference

# packages
source('~/Workspace/RainfallSpectralAnalysis/SpectralAnalysis/function_SetupForGraphics.R')
# more details here: https://github.com/LuigiJR/SpectralAnalysis
library(foreach)
library(doParallel)

site = read.csv('/g/data/os22/users/yu/lst_project/0_code/study_sites.csv')

cl <- makeCluster(12)
print(cl)
registerDoParallel(cl)

# get the time info
TOIs =  seq(as.Date('2013-01-01'), as.Date('2021-12-31'), by='day')
yyyy = format(TOIs, '%Y'); mm = format(TOIs, '%m'); dd = format(TOIs, '%d')

foreach (i=1:12, .combine=cbind) %dopar% {

	source('~/Workspace/RainfallSpectralAnalysis/SpectralAnalysis/function_SetupForGraphics.R')
	
	# paths to LST products
	path2terra = paste0('/g/data/os22/users/yu/lst_project/0_terra_data/', site$sitename[i], '/') # terra lst
	path2esfm = paste0('/scratch/os22/yy4778/lst_project/2_estarfm_output/', site$sitename[i], '/') # baseline
	path2esfmbc = paste0('/scratch/os22/yy4778/lst_project/2_estarfm_output_bc/', site$sitename[i], '/') # bias corrected
	path2time = '/scratch/os22/yy4778/MODIS_data/LST_daytime_viewtime/' # terra time layer
	flux_lst = read.csv(paste0('/g/data/os22/users/yu/lst_project/flux_data/', site$sitename[i], '_lst.csv'))

	# region of interests
	roi = raster(extent(site$xmn[i], site$xmx[i], site$ymn[i], site$ymx[i]), res=0.01, crs=PROJ_LATLON)

	# time range for evaluation
	local_time = as.POSIXct(flux_lst$local_time, tz=site$timezone[i])

	# specify the output file
	ofile = paste0('/g/data/os22/users/yu/lst_project/metrics/estarfm/', site$sitename[i], '_validation.csv')
	basic = data.frame(matrix(nrow=0,ncol=7))
	colnames(basic) = c('time_GMT','flux','terra','estarfm','estarfm_100m','estarfm_bc','estarfm_bc_100m')
	write.table(basic,ofile,row.names=FALSE,sep=',',quote=FALSE)

	# loop the time of interests
	for (k in 1:length(TOIs)){
	
		print(TOIs[k])
		time_file = paste0(path2time, 'MOD11A1_LST_viewtime_', format(TOIs[k], '%Y%m%d'), '.tif') # specify modis time layer

		# if the time file exists
		if (file.exists(time_file)){

			time_rst = raster(time_file)
			lat = site$lat[i]
			mod_time = time_rst[getCellfromLocation(lat, site$lon[i], time_rst)] # read the modis time value

			# In some cases the modis time is extracted as NA. As in the same swatch the modis overpass time is all the same,
			# we can move the lat northward to find a place that has a time value (max is moving 3 degree)
			while (is.na(mod_time) & lat < site$lat[i] + 3){ 
				lat = lat + 0.05
				mod_time = time_rst[getCellfromLocation(lat, site$lon[i], time_rst)]
			}

			# if finally we get a value of modis time
			if (!is.na(mod_time)){
				
				# convert it to GMT time
				mod_GMT = mod_time - site$lon[i] / 15

				# test if this GMT time is at today or yesterday
				if (mod_GMT >= 0) {

					mod_GMT = ISOdatetime(yyyy[k],mm[k],dd[k],as.integer(mod_GMT),as.integer(mod_GMT%%1 * 60),0,tz='GMT')

				} else {

					mod_GMT = mod_GMT + 24
					mod_GMT = ISOdatetime(yyyy[k],mm[k],dd[k],as.integer(mod_GMT),as.integer(mod_GMT%%1 * 60),0,tz='GMT') - 86400
				}

				# read lst rasters
				lstfile = format(TOIs[k],'ESTARFM_LST_daytime_%Y%m%d.tif')
				lstfile_bc = format(TOIs[k],'ESTARFM_LST_daytime_%Y%m%d.tif')
				terrafile = format(TOIs[k],'MOD11A1_LST_daytime_%Y%m%d.tif')

				# if both fusion results exist
				if (file.exists(paste0(path2esfm, lstfile)) & file.exists(paste0(path2esfmbc, lstfile_bc))) {
					
					# terra 1 km
					terra = raster(paste0(path2terra, terrafile))
					terra_value = terra[getCellfromLocation(site$lat[i],site$lon[i],terra)]

					# estarfm 100 m 
					lst = raster(paste0(path2esfm, lstfile))
					esfm_100m = lst[getCellfromLocation(site$lat[i],site$lon[i],lst)]

					# estarfm 1 km
					lst = projectRaster(lst, roi, method='ngb')
					esfm_value = lst[getCellfromLocation(site$lat[i],site$lon[i],lst)]

					# estarfm bc 100 m 
					bc = raster(paste0(path2esfmbc, lstfile_bc))
					bc_100m = bc[getCellfromLocation(site$lat[i],site$lon[i],bc)]
		
					# estarfm bc 1 km
					bc = projectRaster(bc, roi, method='ngb')
					bc_value = bc[getCellfromLocation(site$lat[i],site$lon[i],bc)]

					# use a +- 30 min window to filter out the flux lst data
					whichT = which(local_time >= (mod_GMT - 1800) & local_time <= (mod_GMT + 1800))
					flux_value = mean(flux_lst$lst[whichT], na.rm=TRUE)

					# write all the values to the output file
					write.table(cbind(format(mod_GMT, '%Y-%m-%d %H:%M'),round(flux_value, 3), 
								round(terra_value, 3), round(esfm_value, 3), round(esfm_100m, 3), 
								round(bc_value, 3), round(bc_100m, 3)),
								ofile,row.names=FALSE,col.names=FALSE,sep=',',quote=FALSE,append=TRUE)
				}
			}
		}
	}
}
