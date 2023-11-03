# this script is to calculate lst from Landsat data using split window algorithm
# no runable code here, just for reference

# module load python3/3.8.5 

import numpy as np
import pandas as pd
import rasterio.plot
import rasterio
from pylandtemp import split_window
import tarfile
import re
import os
import sys

namelist = pd.read_csv('/g/data/os22/users/yu/lst_project/0_code/study_sites.csv')
sitenum = int(sys.argv[1])

for tile_number in namelist['tile'][sitenum:(sitenum+1)]:

    l8_path = '/g/data/da82/AODH/USGS/L1/Landsat/C1/%s/' % tile_number
    extr_path = '/g/data/os22/users/yu/lst_project/0_landsat_data/'

    if not os.path.exists('%s/%s' % (extr_path, tile_number)):
        os.mkdir('%s/%s' % (extr_path, tile_number))

    # extract files from the landsat tars
    fname = [f for f in os.listdir(l8_path) if re.match('LC8', f)]
    for t in range(1, len(fname)):
        print(fname[t])
        for tt in os.listdir('%s%s/' % (l8_path, fname[t])):
            if tt.endswith('.tar'):
                tar_name = '%s%s/%s' % (l8_path, fname[t], tt)
                tar = tarfile.open(tar_name, "r:")
                sub_names = tar.getnames()
                for i in sub_names:
                    # RGB images file names with file extensions
                    if i.endswith('B4.TIF'):
                        tar.extract(i, path=extr_path)
                        redband = i
                    elif i.endswith('B5.TIF'):
                        tar.extract(i, path=extr_path)
                        nirband = i
                    elif i.endswith('B10.TIF'):
                        tar.extract(i, path=extr_path)
                        tempband10 = i
                    elif i.endswith('B11.TIF'):
                        tar.extract(i, path=extr_path)
                        tempband11 = i
                tar.close()

                with rasterio.open(extr_path+redband) as src:
                    redImage = src.read(1).astype('f4')

                with rasterio.open(extr_path+nirband) as src:
                    nirImage = src.read(1).astype('f4')
    
                with rasterio.open(extr_path+tempband10) as src:
                    tempImage10 = src.read(1).astype('f4')

                with rasterio.open(extr_path+tempband11) as src:
                    tempImage11 = src.read(1).astype('f4')

                print('starting calculate lst')

                # calculate lst using split window algorithm
                lst_image_split_window = split_window(
                    tempImage10,
                    tempImage11,
                    redImage,
                    nirImage,
                    lst_method = 'jiminez-munoz',
                    emissivity_method = 'avdan',
                    unit='kelvin'
                )

                with rasterio.open('%s/%s/%s' % (extr_path, tile_number, redband.replace('B4','LST')),'w',driver = 'GTiff',
                        height = lst_image_split_window.shape[0],
                        width = lst_image_split_window.shape[1],
                        count = 1,
                        crs = src.crs,
                        transform = src.transform,
                        dtype = lst_image_split_window.dtype
                ) as dst:
                    dst.write(lst_image_split_window, 1)

                os.remove(extr_path+redband)
                os.remove(extr_path+nirband)
                os.remove(extr_path+tempband10)
                os.remove(extr_path+tempband11)
            