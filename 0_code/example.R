library(raster)
library(foreach)
library(doParallel)

cl = makeCluster(25)
registerDoParallel(cl)