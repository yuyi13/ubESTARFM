# This is the updated ozflux processing script which always uses the local standard time 
# and expilctly claims the 'second' timestep in the TOI (time of interests)

# packages
library(foreach)
library(doParallel)
library(itertools)
library(lubridate)
source('~/Workspace/RainfallSpectralAnalysis/SpectralAnalysis/function_SetupForGraphics.R')
# more details here: https://github.com/LuigiJR/SpectralAnalysis

# site info
ozfluxLL = read.csv('/datasets/work/d61-af-soilmoisture/work/users/yu/git_repo/ubESTARFM/5_ozflux_lst/0_study_sites.csv')
outpath  = '/datasets/work/d61-af-soilmoisture/work/users/yu/git_repo/ubESTARFM/5_ozflux_lst/'

cl <- makeCluster(12)
print(cl)
registerDoParallel(cl)

# do it in parallel
foreach (i=1:nrow(ozfluxLL), .combine=cbind, .packages=c('ncdf4', 'raster', 'lubridate')) %dopar% {

    sitename = ozfluxLL$sitename[i]
	print(paste0('start extracting value for site ', sitename))

	ofile = paste0(outpath, sitename, '_lst.csv')
    basic = data.frame(matrix(nrow=0,ncol=5))
	colnames(basic) = c('local_time', 'longwave_u', 'longwave_d', 'emis', 'lst')
    write.table(basic,ofile,row.names=FALSE,sep=',',quote=FALSE)

	# read netcdf files from ozflux server
    x = nc_open(paste0('https://dap.ozflux.org.au/thredds/dodsC/ozflux/sites/',sitename,'/L3/default/',sitename,'_L3.nc'))
    L_u = ncvar_get(x, 'Flu'); L_u[L_u == -9999] = NA
    L_d = ncvar_get(x, 'Fld'); L_d[L_d == -9999] = NA
    local_time = ncvar_get(x, 'time')
	nc_close(x)

	# *******************************************
	# this is the area where adjustments are made
	# in the updated version, we only consider the standard time
    TZ = ozfluxLL$standard_timezone[i]

	# in the updated version, we explicitly claim the 'second' timestep in the TOI (time of interests)
    local_time = time_GMT = ISOdatetime(1800,1,1,0,0,0,tz=TZ) + seconds(local_time*3600*24)
    attr(time_GMT, 'tzone') = 'GMT'
	# *******************************************

	for (k in 1:length(local_time)) {

		# calculate Ts
		month = as.numeric(format(time_GMT[k], '%m'))
		emis_value = as.numeric(ozfluxLL[i, 12 + month])

		flux_lst = ((L_u[k] - (1 - emis_value) * L_d[k]) / (5.670374e-8 * emis_value)) ^ (1/4)
		print(flux_lst)
		write.table(cbind(format(local_time[k], '%Y-%m-%d %H:%M'), 
					round(L_u[k], 4), round(L_d[k], 4),
					emis_value, round(flux_lst, 4)),
					ofile, row.names=FALSE, col.names=FALSE, sep=',', quote=FALSE, append=TRUE)
    }
}
