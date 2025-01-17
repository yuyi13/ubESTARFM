# This is an example to run ubESTARFM in parallel

# load essential packages
requiredPackages = c('raster','foreach','doParallel')
for (k in requiredPackages){
    if (!require(k, character.only = TRUE)) install.packages(k, repos='https://cloud.r-project.org')
    library(k, character.only = TRUE)
}

# load ubESTARFM
source('0_algorithm/ubESTARFM.R')

# set a cluster to allow ubETSRAFM to run in parallel
cl = makeCluster(4)
registerDoParallel(cl)

# set the paths of input, tmp and output
data_path = '1_test_data/'
tmp_path  = '2_tmp_path/'; if (!dir.exists(tmp_path)) dir.create(tmp_path)
out_path  = '3_output/fused_result.tif'; if (!dir.exists('3_output')) dir.create('3_output')

# load the data
rst_fine1   = raster(paste0(data_path, 'Landsat_LST_cloudrm_20160205.tif')) # fine resolution image in pair 1
rst_fine2   = raster(paste0(data_path, 'Landsat_LST_cloudrm_20160308.tif')) # fine resolution image in pair 2
rst_coarse1 = raster(paste0(data_path, 'MOD11A1_LST_cloudrm_20160205.tif')) # coarse resolution image in pair 1
rst_coarse2 = raster(paste0(data_path, 'MOD11A1_LST_cloudrm_20160308.tif')) # coarse resolution image in pair 2
rst_coarse0 = raster(paste0(data_path, 'MOD11A1_LST_cloudrm_20160218.tif')) # coarse resolution image on prediction date

# run the ubESTARFM, and it will write into the out_path
# some parameters are set as default values, you can change them as you see fit
# set if you want to run in baseline mode or unbiased mode
# method = 'zero bias'; or you can set it as 'baseline', which means ESTARFM

ubESTARFM(w = 25, DN_min = 250, DN_max = 350, patch_long = 200, tmp_path, out_path, method='zero bias',
          rst_fine1, rst_fine2, rst_coarse1, rst_coarse2, rst_coarse0)

# stop the cluster
stopCluster(cl)