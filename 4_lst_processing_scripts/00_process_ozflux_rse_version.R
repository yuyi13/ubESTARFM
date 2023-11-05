# *************************************************************************************************** #
# This is the ozflux processing script used in our RSE paper                                          #
# which used a temporal window to filter in-situ lst data, and considered local daylight saving zones #
# however, a more precise method may be always using the local standard time                          #
# and expilctly claiming the 'second' timestep in the TOI (time of interests)                         #
# hence we have an updated version of this script in this repo                                        #
# *************************************************************************************************** #

# packages
library(foreach)
library(doParallel)
library(itertools)
source('~/Workspace/RainfallSpectralAnalysis/SpectralAnalysis/function_SetupForGraphics.R')
# more details here: https://github.com/LuigiJR/SpectralAnalysis

# site info
ozfluxLL = read.csv('/datasets/work/d61-af-soilmoisture/work/users/yu/git_repo/ubESTARFM/5_ozflux_lst/0_study_sites.csv')
outpath  = '/datasets/work/d61-af-soilmoisture/work/users/yu/git_repo/ubESTARFM/5_ozflux_lst/'

cl <- makeCluster(12)
print(cl)
registerDoParallel(cl)

# do it in parallel
foreach (i=1:nrow(ozfluxLL), .combine=cbind, .packages=c('ncdf4', 'raster')) %dopar% {

    sitename = ozfluxLL$sitename[i]
	print(paste0('start extracting value for site ', sitename))

	ofile = paste0(outpath, sitename, '_lst.csv')
    basic = data.frame(matrix(nrow=0,ncol=5))
	colnames(basic) = c('local_time', 'longwave_u', 'longwave_d', 'emis', 'lst')
    write.table(basic,ofile,row.names=FALSE,sep=',',quote=FALSE)

	# read netcdf files from ozflux server
    x = nc_open(paste0('https://dap.ozflux.org.au/thredds/dodsC/ozflux/sites/',sitename,'/L3/default/',sitename,'_L3.nc'))
    L_u = ncvar_get(x, 'Flu')
    L_d = ncvar_get(x, 'Fld')
    local_time = ncvar_get(x, 'time')
	nc_close(x)

	# *******************************************
	# this is the area where adjustments are made
	# in RSE version, we considered the daylight saving time
	# in the future, we should always use local standard time; i.e., tz=site$standard_timezone[i]
    TZ = ozfluxLL$timezone[i]

	# in the RSE version, we did not explicitly claimed the 'second' timestep in the TOI (time of interests)
    local_time = time_GMT = ISOdatetime(1800,1,1,0,0,0,tz=TZ) + local_time*3600*24
    attr(time_GMT, 'tzone') = 'GMT'
	# *******************************************

	for (k in 1:length(local_time)) {

		# calculate Ts
		month = as.numeric(format(time_GMT[k], '%m'))
		emis_value = as.numeric(ozfluxLL[i, 12 + month])

		flux_lst = ((L_u[k] - (1 - emis_value) * L_d[k]) / (5.670374e-8 * emis_value)) ^ (1/4)
		print(flux_lst)
		write.table(cbind(format(local_time[k], '%Y-%m-%d %H:%M'), L_u[k], L_d[k], emis_value, flux_lst),
					ofile, row.names=FALSE, col.names=FALSE, sep=',', quote=FALSE, append=TRUE)
    }
}
