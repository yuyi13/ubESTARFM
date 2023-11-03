# this script is to mask out the ECOSTRESS LST data using the cloud mask data
# no runable code here, just for reference

library(stringr)
library(ncdf4)
library(foreach)
library(doParallel)

# convert integer to binary
number2binary = function(number, noBits) {
       binary_vector = rev(as.numeric(intToBits(number)))
       if(missing(noBits)) {
          return(binary_vector)
       } else {
          binary_vector[-(1:(length(binary_vector) - noBits))]
       }
    }

site = read.csv('/datasets/work/d61-af-soilmoisture/work/lst_project/0_code/study_sites.csv')
path2clod = '/datasets/work/d61-af-soilmoisture/work/ECOSTRESS/CLOUDMASK/'
path2loca = '/datasets/work/d61-af-soilmoisture/work/ECOSTRESS/GEOLOCATION/'
path2tiff = '/datasets/work/d61-af-soilmoisture/work/ECOSTRESS/tiff/'
path2metrics = '/datasets/work/d61-af-soilmoisture/work/ECOSTRESS/estarfm_metrics/' # this contains the MODIS overpass info

PROJ_LATLON = '+proj=longlat +datum=WGS84' # proj info

fl = list.files(path2tiff, pattern='.tif$')

cl <- makeCluster(12)
print(cl)
registerDoParallel(cl)

foreach (f = 1:12, .combine=cbind, .packages=c('stringr', 'ncdf4')) %dopar% {

    source('~/Workspace/RainfallSpectralAnalysis/SpectralAnalysis/function_SetupForGraphics.R')

    # get the region of interests
    roi = raster(extent(site$xmn[f], site$xmx[f], site$ymn[f], site$ymx[f]), res=0.001, crs=PROJ_LATLON)
    # read the metrics file, which has the MODIS time info    
    metrics = read.csv(paste0(path2metrics, site$sitename[f], '_validation.csv'))
    # specify the outpath for each site
    outpath = paste0('/datasets/work/d61-af-soilmoisture/work/ECOSTRESS/flux_subsets/', site$sitename[f])
    print('**************************')
    print(site$sitename[f])
    print('**************************')

    if (!dir.exists(outpath)) dir.create(outpath, recursive=TRUE)

    for (k in 1:length(fl)){
        x = raster(paste0(path2tiff, fl[k]))
        overlap = intersect(extent(x), extent(roi)) # to see if ecostress is overlapping with roi

        # it they are overlapped, then consider next step
        if (!is.null(overlap)){ 
            
            # to see if ecostress overpass time match modis time
            time_str = str_split(fl[k], '_')[[1]][3]
            year = str_sub(time_str, 1, 4); mon = str_sub(time_str, 5, 6); day = str_sub(time_str, 7, 8)
            hour = str_sub(time_str, 10, 11); min = str_sub(time_str, 12, 13); sec = str_sub(time_str, 14, 15)
            Teco = ISOdatetime(year, mon, day, hour, min, sec, tz='GMT') # ECOSTRESS overpass time
            whichT = which(as.POSIXct(metrics$time_GMT, tz='GMT') < Teco + 1800 & as.POSIXct(metrics$time_GMT, tz='GMT') > Teco - 1800)
            
            # test if there is a match between ecostress overpass time and modis time
            # if there is a match, then proceed
            if (length(whichT) > 0){ 
                print(fl[k]); print(metrics$time_GMT[whichT])
                
                # get the ECOSTRESS mask for the date of interest
                mask_file = nc_open(list.files(path2clod, pattern=str_sub(fl[k], 15, 29), full.names=TRUE))
                loca_file = nc_open(list.files(path2loca, pattern=str_sub(fl[k], 15, 29), full.names=TRUE))

                lon = ncvar_get(loca_file, 'Geolocation/longitude')
                lat = ncvar_get(loca_file, 'Geolocation/latitude')
                msk = ncvar_get(mask_file, 'SDS/CloudMask')

                # cloudmask is in integer; convert it to binary
                # The idea is that read bit 1 (the second) and bit 3 (the fourth). If there's any cloud here, 
                # then either the bit 1 or bit 3 will have the value of 1
                # source: https://ecostress.jpl.nasa.gov/faq

                # *** you need to read bit fields from right to left ***
                # *** because hdf files are written in big endian format ***
                # https://www.r-bloggers.com/2012/12/modis-qc-bits/

                # Bit Field: Long Name | Result
                # 0: Cloud Mask Flag | 0=not determined; 1=determined
                # 1: Final Cloud Plus Region-growing | 0=no; 1=yes
                # 2: Final Cloud, either one of bits 2, 3 ,or 4 set | 0=no; 1=yes
                # 3: Band 4 Brightness Threshold Test | 0=no; 1=yes
                # 4: Band 4-5 Thermal Difference Test | 0=no; 1=yes
                # 5: Land/Water Mask | 0=land; 1=water

                # hence the below means we create cloud mask using the bit 4 (thermal diff) and bit 2 (final cloud without region-growing)
                # though this is different with the ecostress tutorial, I believe it is acceptable
                # next time we might just use the best quality pixels as indicated in QC layer (bit 1 = 0 and bit 2 = 0)
                # (then avoid the usage of cloud mask data)

                for (b in 1:63){
                    msk[msk == b] = number2binary(b, 6)[2] + number2binary(b, 6)[4]
                }

                nc_close(mask_file); nc_close(loca_file)

                df = data.frame(x=as.vector(lon), y=as.vector(lat), lst=as.vector(msk)) # make it a dataframe

                rst = raster(extent(df[,1:2]), ncol=5632, nrow=5400, crs=PROJ_LATLON) # make the geolocation a rasterlayer 
                rst = rasterize(df[, 1:2], rst, df[,3], fun=mean) # involve the cloudmask data
                rst[rst > 0] = NA

                # use the generated cloud mask data to mask the ecostress data
                x = mask(x, mask=rst)

                # project to the roi extent
                xx = projectRaster(x, roi, method='ngb',
                     filename = paste0(outpath, '/ECOSTRESS_LST_', time_str, '.tif'), overwrite=TRUE)
            }
        }
    }
}
