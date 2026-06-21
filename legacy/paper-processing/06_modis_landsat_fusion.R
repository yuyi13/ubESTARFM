# this script is to perform the (ub)ESTARFM to fuse MODIS and Landsat LST data
# no runable code here, just for reference

# load packages
source('0_algorithm/ubESTARFM.R')
library(stringr)
library(foreach)
library(doParallel)
library(ImageFusion) # this is a subtle package for temporal interpolation
# more details here: https://github.com/JohMast/ImageFusion

# site info
site_df = read.csv('/datasets/work/d61-af-soilmoisture/work/users/yu/code/study_sites.csv')

# parallel computing
cl = makeCluster(25)
registerDoParallel(cl)

# get the site index from command line
args = commandArgs(trailingOnly=T)
k = as.numeric(args[1])

# ubESTRAFM
method = 'zero bias'

# if method = 'baseline', that means ESTARFM

PROJ_LATLON = '+proj=longlat +datum=WGS84' # proj info

# modis data paths
terra_lst = paste0('/datasets/work/d61-af-soilmoisture/work/lst_project/0_terra_data/', site_df$sitename[k], '/')
terra_interp = paste0('/datasets/work/d61-af-soilmoisture/work/lst_project/1_terra_interp/', site_df$sitename[k], '/')

# landsat data paths
landsat_lst = paste0('/datasets/work/d61-af-soilmoisture/work/lst_project/0_landsat_masked/', site_df$sitename[k], '/')
landsat_interp = paste0('/datasets/work/d61-af-soilmoisture/work/lst_project/1_landsat_interp/', site_df$sitename[k], '/')

# specify the output and tmp paths
if (method == 'baseline'){

    out_dir = paste0('/datasets/work/d61-af-soilmoisture/work/lst_project/2_estarfm_output/', site_df$sitename[k], '/')
    tmp_path = paste0('/datasets/work/d61-af-soilmoisture/work/tmp/', site_df$sitename[k], '/baseline/')

} else if (method == 'zero bias'){

    out_dir = paste0('/datasets/work/d61-af-soilmoisture/work/lst_project/2_estarfm_output_bc/', site_df$sitename[k], '/')
    tmp_path = paste0('/datasets/work/d61-af-soilmoisture/work/tmp/', site_df$sitename[k], '/bias_correct/')
}

# create directories if not exist
if (!dir.exists(terra_interp)) dir.create(terra_interp); if (!dir.exists(landsat_interp)) dir.create(landsat_interp, recursive=TRUE)
if (!dir.exists(tmp_path)) dir.create(tmp_path, recursive=TRUE); if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# region and date of interests
Roi_1km = raster(extent(site_df$xmn[k], site_df$xmx[k], site_df$ymn[k], site_df$ymx[k]), res=0.01, crs=PROJ_LATLON)
Roi_100m = raster(extent(site_df$xmn[k], site_df$xmx[k], site_df$ymn[k], site_df$ymx[k]), res=0.001, crs=PROJ_LATLON)
DOI = seq(as.Date('2013-01-01'),as.Date('2021-12-31'),by='day')

# do temporal interpolation to fill landsat gaps
# this is just to ensure that the landsat data are more complete, which can benefit the fusion results
# we will remove the interpolated areas later
lsat_fl = list.files(landsat_lst, pattern='Landsat_LST_', full.names=T)
ImageFusion::imginterp_task(filenames = lsat_fl,
                            dates = 1:length(lsat_fl), limit_days = 3,
                            invalid_ranges = '[-inf,250]',
                            out_prefix = landsat_interp)

# do temporal interpolation to fill modis gaps
# this is just to ensure that the modis data are more complete, which can benefit the fusion results
# we will remove the interpolated areas later
terra_fl = list.files(terra_lst, pattern='MOD11A1_LST_daytime', full.names=T)
ImageFusion::imginterp_task(filenames = terra_fl,
                            dates = 1:length(terra_fl), limit_days = 8,
                            invalid_ranges = '[-inf,250]',
                            out_prefix = terra_interp)

# landsat data for training
# we still use the raw landsat lst to select near clear-sky scenes
train_lsat = list.files(landsat_lst, pattern='Landsat_LST_', full.names=T)
lsat_dates = as.Date(str_sub(train_lsat, -12, -5), format='%Y%m%d')
valid_idx = which(file.info(train_lsat)$size > 2.5e6) # select files larger than 2.5MB

# further check if modis files satisfies size requirements
# we still use the raw landsat lst to select near clear-sky scenes
mdis_files = paste0(terra_lst, format(lsat_dates[valid_idx], 'MOD11A1_LST_daytime_%Y%m%d.tif'))
valid_idx = valid_idx[which(file.info(mdis_files)$size > 1e6)] # select files larger than 1MB

# data for training
# use interpolated data
train_lsat = paste0(landsat_interp, format(lsat_dates[valid_idx], 'Landsat_LST_cloudrm_%Y%m%d.tif')) 
train_mdis = paste0(terra_interp, format(lsat_dates[valid_idx], 'MOD11A1_LST_daytime_%Y%m%d.tif'))

# data for prediction
pred_terra_day = list.files(terra_interp, pattern='MOD11A1_LST_daytime')
pred_terra_day_fullname = list.files(terra_interp, pattern='MOD11A1_LST_daytime', full.names=T)

# output names
out_day = paste0(out_dir, str_replace(pred_terra_day, 'MOD11A1_', 'ubESTARFM_'))
pred_time = as.Date(str_sub(pred_terra_day, -12, -5), format='%Y%m%d')

# start ubestarfm
print('start to do the lst fusion')

for (t in 1:(length(train_lsat)-1)){

	# landsat and modis pairs
	l_pair = train_lsat[c(t, t+1)]
    m_pair = train_mdis[c(t, t+1)]

	print(list(l_pair, m_pair))

    # find out dates that are located within the filtered period
    pred_idx = which(pred_time >= lsat_dates[valid_idx[t]] & pred_time < lsat_dates[valid_idx[t+1]])

    for (tt in 1:length(pred_idx)){

    # data fusion using unbiased estarfm
    fused_rst = ubESTARFM(tmp_path=tmp_path, out_path=out_day[pred_idx][tt], method=method,
                rst_fine1 = raster(l_pair[1]), rst_fine2 = raster(l_pair[2]), rst_coarse1 = raster(m_pair[1]), 
                rst_coarse2 = raster(m_pair[2]), rst_coarse0 = raster(pred_terra_day_fullname[pred_idx][tt]))
    }
}

# list all the fused data and apply the modis cloud mask to remove interpolated areas
fl = list.files(out_dir, pattern='.tif$', full.names=TRUE)

foreach (f = 1:length(fl), .combine=cbind, .packages=c('stringr', 'raster')) %dopar% {

    rst = raster(fl[f])

    doi = str_sub(fl[f], -12, -5)
    mod_rst = raster(paste0(terra_lst, 'MOD11A1_LST_daytime_', doi, '.tif'))
    rst = mask(rst, mask=mod_rst)
    writeRaster(rst, fl[f], overwrite=TRUE)
    file.remove(paste0(fl[f], '.aux.xml'))
}
