# this script is to read ECOSTRESS LST data and convert it to tiff files
# no runable code here, just for reference

library(raster)
library(ncdf4)
library(stringr)
library(foreach)
library(doParallel)

path2eco = '/datasets/work/d61-af-soilmoisture/work/ECOSTRESS/'
PROJ_LATLON = '+proj=longlat +datum=WGS84' # proj info

fl = list.files(paste0(path2eco, 'LST'), pattern = '.h5$')

print(paste0('There are how many nodes: ', detectCores()))

cl <- makeCluster(10)
print(cl)
registerDoParallel(cl)

# below are how these ECOSTRESS files look like

# LST and LSE: ECOSTRESS_L2_LSTE_03238_004_20190131T025006_0601_02.h5
# geolocation: ECOSTRESS_L1B_GEO_03238_004_20190131T025006_0601_02.h5
# cloud mask:  ECOSTRESS_L2_CLOUD_03238_004_20190131T025006_0601_02.h5

foreach (k = 1:length(fl), .combine=cbind, .packages=c('raster', 'ncdf4', 'stringr')) %dopar% {

    # find the corresponding geolocation file
    geofile = paste0(path2eco, 'GEOLOCATION/', str_replace(fl[k], 'ECOSTRESS_L2_LSTE', 'ECOSTRESS_L1B_GEO'))

    print(geofile)

    if (file.exists(geofile)){
        
        print(paste0('Doing the ', k, 'th file for ECOSTRESS LST'))

        # open geolocation file
        geo = nc_open(geofile)
        # open lst file
        x = nc_open(paste0(path2eco, 'LST/', fl[k]))

        lon = ncvar_get(geo, 'Geolocation/longitude')
        lat = ncvar_get(geo, 'Geolocation/latitude')
        lst = ncvar_get(x, 'SDS/LST')

        # dim(lst) 5400 5632
        nc_close(geo); nc_close(x)

        # make it a dataframe
        df = data.frame(x=as.vector(lon), y=as.vector(lat), lst=as.vector(lst)) 

        # make the geolocation a rasterlayer 
        rst = raster(extent(df[,1:2]), ncol=5632, nrow=5400, crs=PROJ_LATLON) 
        # involve the lst data
        rst = rasterize(df[, 1:2], rst, df[,3], fun=mean) 

        # write raster with extent being a part of the name
        # so we can use this info to filter out the raster in the next step
        writeRaster(rst, paste0(path2eco, 'tiff/ECOSTRESS_LST_', str_sub(geofile, -26, -12), '_', round(extent(rst)[1], 1), 
                    '_', round(extent(rst)[2], 1), '_', round(extent(rst)[3], 1), '_', round(extent(rst)[4], 1), '.tif'),
                    overwrite=TRUE)
    }
}
