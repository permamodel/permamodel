"""
couple_cru_KUGeo.py

Example that couples cruAKtemp to KUGeo!

"""

from __future__ import print_function

import numpy as np

from cruAKtemp.bmi_cruAKtemp import BmiCruAKtempMethod
from permamodel.components.bmi_Ku_component import BmiKuMethod


# Initial both components with local config files
cru_config_file = "./cruAKtemp_component.cfg"
kug_config_file = "./KuGeo_component.cfg"

cru = BmiCruAKtempMethod()
kug = BmiKuMethod()

cru.initialize(cru_config_file)
kug.initialize(kug_config_file)

lat = cru.get_value('latitude')
lon = cru.get_value('longitude')

nx = lat.shape[0]
ny = lon.shape[1]

onset = cru._model.first_date.year
final = cru._model.last_date.year

ntime = final-onset+1

kug.set_value('latitude', lat)
kug.set_value('longitude', lon)
kug.set_value('datetime__start',onset)
kug.set_value('datetime__end',final)

# The following two variables are required only for writing netcdf outputs.

kug.output_alt = np.zeros((ntime,nx,ny))* 0 -999.9
kug.output_tps = np.zeros((ntime,nx,ny))* 0 -999.9

T_air = np.zeros((nx,ny))

for i in np.arange(3*12):
    
    cru.update()
    
    T_air = cru.get_value("atmosphere_bottom_air__temperature_year") + T_air
    
    if np.mod(i+1,12) == 0: # When get the 12th monthly data, Ku component will be implemented.
        
        kug.set_value('atmosphere_bottom_air__temperature', T_air/12.)
        kug.update()
        T_air = np.zeros((nx,ny))
#
kug.finalize()
#
print(kug.get_value("soil__active_layer_thickness"))