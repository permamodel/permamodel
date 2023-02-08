# -*- coding: utf-8 -*-
"""
Created on Mon May 15 16:45:54 2017

@author: kangwang
"""

import numpy as np
from netCDF4 import Dataset


def write_out_ncfile(output_file, lon, lat, time0, data, longname, units):
    import numpy as np
    from netCDF4 import Dataset

    n_lat = np.size(lat)
    n_lon = np.size(lon)

    #        print np.shape(varname)

    ALT = data + 0.0  # self.mask;
    idx = np.where(np.isnan(ALT))
    ALT[idx] = -999.99

    # print output_file[-1-2]

    # Open a file to save the final result
    w_nc_fid = Dataset(output_file + ".nc", "w", format="NETCDF4")

    # ==== Latitude ====

    w_nc_fid.createDimension("lat", n_lat)  # Create Dimension
    lats = w_nc_fid.createVariable("lat", np.dtype("float32").char, ("lat",))
    lats.units = "degrees_north"
    lats.standard_name = "latitude"
    lats.long_name = "latitude"
    lats.axis = "Y"
    lats[:] = lat

    # ==== Longitude ====

    w_nc_fid.createDimension("lon", n_lon)  # Create Dimension
    lons = w_nc_fid.createVariable("lon", np.dtype("float32").char, ("lon",))
    lons.units = "degrees_east"
    lons.standard_name = "longitude"
    lons.long_name = "longitude"
    lons.axis = "X"
    lons[:] = lon

    # ==== Time ====

    w_nc_fid.createDimension("time", np.size(time0))  # Create Dimension
    time = w_nc_fid.createVariable("time", np.dtype("float32").char, ("time",))
    time.units = "Year"
    #         time.standard_name = 'longitude'
    #         time.long_name = 'longitude'
    time.axis = "Z"
    time[:] = time0

    # ==== Data ====
    temp = w_nc_fid.createVariable(
        "data", np.dtype("float32").char, ("time", "lat", "lon")
    )
    temp.units = units
    temp.missing_value = -999.99
    temp.long_name = longname
    temp[:] = ALT
    #
    w_nc_fid.close()  # close the new file


data_file = "Raw_ERA_Dataset.nc"
mask_file = "mask.nc"


fid = Dataset(mask_file)
mask = fid.variables["lsm"][0]
mask[np.where(mask == 0)] = np.nan
fid.close()

# MAAT:
data = Dataset(data_file)

time = data.variables["time"][:]
lat = data.variables["latitude"][:]
lon = data.variables["longitude"][:]

time0 = np.arange(2014, 2017)

t2m = data.variables["t2m"][:]

t2m = t2m - 273.15

t2m = np.reshape(
    t2m, (np.size(t2m, axis=0) / 12, 12, np.size(t2m, axis=1), np.size(t2m, axis=2))
)

ta = np.mean(t2m, axis=1)

aa = 0.5 * (np.max(t2m, axis=1) - np.min(t2m, axis=1))

sd = data.variables["sd"][:]

sd = np.reshape(
    sd, (np.size(sd, axis=0) / 12, 12, np.size(sd, axis=1), np.size(sd, axis=2))
)
sd = np.mean(sd, axis=1)

sden = data.variables["rsn"][:]
sden = np.reshape(
    sden, (np.size(sden, axis=0) / 12, 12, np.size(sden, axis=1), np.size(sden, axis=2))
)
sden = np.mean(sden, axis=1)

vwc = data.variables["swvl1"][:]
vwc = np.reshape(
    vwc, (np.size(vwc, axis=0) / 12, 12, np.size(vwc, axis=1), np.size(vwc, axis=2))
)
vwc = np.mean(vwc, axis=1)

# mask = sd; mask[np.where(sd>0.)] = 1

ta = ta * mask
aa = aa * mask
sd = sd * mask
sden = sden * mask
vwc = vwc * mask

write_out_ncfile("../ta", lon, lat, time0, ta, "MAAT", "deg C")
write_out_ncfile("../aa", lon, lat, time0, aa, "Ampt", "deg C")
write_out_ncfile("../snd", lon, lat, time0, sd, "SnowDepth", "m")
write_out_ncfile("../rsn", lon, lat, time0, sden, "SnowDensity", "kg m-3")
write_out_ncfile("../vwc", lon, lat, time0, vwc, "Vol.Water.Content", "m3 m-3")
