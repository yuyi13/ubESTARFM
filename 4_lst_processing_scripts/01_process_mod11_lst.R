# this script is to process the MODIS LST data downloaded from NASA LPDAAC
# no runable code here, just for reference

# packages
library(raster)
library(foreach)
library(doParallel)

# we have part of the data downloaded in scratch and part of the data archived in u39
# u39 is a project name in NCI (Australian national HPC cluster)
path2modis_u39 = '/g/data/u39/public/data/modis/lpdaac-tiles-c6/MOD11A1.006/'
path2modis_downloaded = '/scratch/os22/yy4778/MODIS_data/LST_daytime_daily/hdf'

# geo-template for australia
PROJ_LATLON = '+proj=longlat +datum=WGS84'
tem = raster(extent(112,154,-45,-10), res=0.01, crs=PROJ_LATLON)

# parallel computing
cl <- makeCluster(20)
print(cl)
registerDoParallel(cl)

# read data from u39
Dates = seq(as.Date('2013-01-01'),as.Date('2019-10-31'),by='day')

foreach (k = 1:length(Dates), .combine=cbind, .packages=c('raster', 'gdalUtils')) %dopar% {

	DOI = paste0(path2modis_u39, format(Dates[k], '%Y.%m.%d/'))
	print(DOI)
	fl = list.files(DOI, pattern='.hdf$', full.names=TRUE) # list files ending with .hdf extension

	# if there are data available
	if (length(fl) > 0){

		lst_list = c()
		time_list = c()

		for (i in 1:length(fl)){
			print(i)
			lst = try(raster(get_subdatasets(fl[i])[1]))
			view_time = try(raster(get_subdatasets(fl[i])[3]))
			lst_list = c(lst_list, lst)
			time_list = c(time_list, view_time)
		}

		# mosaic the data

		half1 = lst_list[1:10]; half1$fun='mean'; aus_lst1 = do.call(mosaic, half1)
		half2 = lst_list[11:length(lst_list)]; half2$fun='mean'; aus_lst2 = do.call(mosaic, half2)

		time1 = time_list[1:10]; time1$fun='mean'; aus_time1 = do.call(mosaic, time1)
		time2 = time_list[11:length(time_list)]; time2$fun='mean'; aus_time2 = do.call(mosaic, time2)

		ll = mosaic(aus_lst1, aus_lst2, fun='mean')
		tt = mosaic(aus_time1, aus_time2, fun='mean')

		# project the mosaiced data to geo-template
		ll = projectRaster(ll, tem, method='ngb',
		 	filename=paste0('/scratch/os22/yy4778/MODIS_data/LST_daytime_daily/MOD11A1_LST_daytime_',format(Dates[k], '%Y%m%d.tif')), overwrite=TRUE)

		tt = projectRaster(tt, tem, method='ngb',
		 	filename=paste0('/scratch/os22/yy4778/MODIS_data/LST_daytime_viewtime/MOD11A1_LST_viewtime_',format(Dates[k], '%Y%m%d.tif')), overwrite=TRUE)
	
	} else {

		# if there is no data available, write a blank raster
		writeRaster(tem, paste0('/scratch/os22/yy4778/MODIS_data/LST_daytime_daily/MOD11A1_LST_daytime_',format(Dates[k], '%Y%m%d.tif')), overwrite=TRUE)
		writeRaster(tem, paste0('/scratch/os22/yy4778/MODIS_data/LST_daytime_viewtime/MOD11A1_LST_viewtime_',format(Dates[k], '%Y%m%d.tif')), overwrite=TRUE)
	}
}

# read downloaded data from scratch 
Dates = seq(as.Date('2019-11-01'),as.Date('2021-12-31'),by='day')

# a similar loop with above
foreach (k = 1:length(Dates), .combine=cbind, .packages=c('raster', 'gdalUtils')) %dopar% {

	DOI = format(Dates[k], 'MOD11A1.A%Y%j')
	print(DOI)
	fl = list.files(path2modis_downloaded, pattern=DOI, full.names=TRUE)

	if (length(fl) > 0){

		lst_list = c()
		time_list = c()

		for (i in 1:length(fl)){
			print(i)
			lst = try(raster(get_subdatasets(fl[i])[1]))
			view_time = try(raster(get_subdatasets(fl[i])[3]))
			lst_list = c(lst_list, lst)
			time_list = c(time_list, view_time)
		}

		part1 = lst_list[1:10]; part1$fun='mean'; aus_lst1 = do.call(mosaic, part1)
		part2 = lst_list[11:20]; part2$fun='mean'; aus_lst2 = do.call(mosaic, part2)
		part3 = lst_list[21:length(lst_list)]; part3$fun='mean'; aus_lst3 = do.call(mosaic, part3)

		time1 = time_list[1:10]; time1$fun='mean'; aus_time1 = do.call(mosaic, time1)
		time2 = time_list[11:20]; time2$fun='mean'; aus_time2 = do.call(mosaic, time2)
		time3 = time_list[21:length(lst_list)]; time3$fun='mean'; aus_time3 = do.call(mosaic, time3)

		ll = mosaic(aus_lst1, aus_lst2, aus_lst3, fun='mean')
		tt = mosaic(aus_time1, aus_time2, aus_time3, fun='mean')

		ll = projectRaster(ll, tem, method='ngb',
		 	filename=paste0('/scratch/os22/yy4778/MODIS_data/LST_daytime_daily/MOD11A1_LST_daytime_',format(Dates[k], '%Y%m%d.tif')), overwrite=TRUE)

		tt = projectRaster(tt, tem, method='ngb',
		 	filename=paste0('/scratch/os22/yy4778/MODIS_data/LST_daytime_viewtime/MOD11A1_LST_viewtime_',format(Dates[k], '%Y%m%d.tif')), overwrite=TRUE)
	
	} else {

		writeRaster(tem, paste0('/scratch/os22/yy4778/MODIS_data/LST_daytime_daily/MOD11A1_LST_daytime_',format(Dates[k], '%Y%m%d.tif')), overwrite=TRUE)
		writeRaster(tem, paste0('/scratch/os22/yy4778/MODIS_data/LST_daytime_viewtime/MOD11A1_LST_viewtime_',format(Dates[k], '%Y%m%d.tif')), overwrite=TRUE)
	}
}
