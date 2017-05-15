# -*- coding: utf-8 -*-
"""
Created on Mon May 15 16:45:54 2017

@author: kangwang
"""

import numpy as np
from netCDF4 import Dataset

data_file = '_grib2netcdf-atls17-a562cefde8a29a7288fa0b8b7f9413f7-ts0Hg5.nc';
mask_file = 'mask.nc'

fid  = Dataset(mask_file)
mask = fid.variables['lsm'][0]
mask[np.where(mask==0)] = np.nan
fid.close()

# MAAT:
data = Dataset(data_file)

time = data.variables['time'][:]

t2m  = data.variables['t2m'][:] - 237.15

t2m = np.reshape(t2m, (np.size(t2m,axis=0)/12,12,np.size(t2m,axis=1),np.size(t2m,axis=2)))

ta  = np.mean(t2m, axis=1)