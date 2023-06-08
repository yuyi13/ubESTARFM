# load required packages
library(raster)
if(!require('fields',character.only = TRUE)) install.packages('fields'); library('fields',character.only = TRUE)

# colour scheme
TemperatureRamp =
colorRampPalette(c("white","darkblue","darkblue","blue","cyan","yellow","orange","red","darkred"))

# read the modis data on prediction date
x = raster('1_test_data/MOD11A1_LST_cloudrm_20160218.tif')

# read the fused result
y = raster('3_output/fused_result.tif')

# generate a png image
png('3_output/visualisation.png', width=900, height=400)

m = rbind(c(1,2))
layout(m); par(mar = c(2,4,2,2))
image.plot(x, col=TemperatureRamp(64), zlim=c(305,315), xlab='', ylab='', main='MODIS LST')
image(y, col=TemperatureRamp(64), zlim=c(305,315), xlab='', ylab='', main='ubESTARFM LST')

dev.off()
