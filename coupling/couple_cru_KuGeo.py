"""
couple_cru_FNGeo.py

Example that couples cruAKtemp to FNGeo!

"""

from __future__ import print_function

import numpy as np

from cruAKtemp.bmi_cruAKtemp import BmiCruAKtempMethod
from permamodel.components.bmi_frost_number_Geo import BmiFrostnumberGeoMethod
from permamodel.components.bmi_Ku_component import BmiKuMethod


# Initial both components with local config files
cru_config_file = "./cruAKtemp_component.cfg"
kug_config_file = "./KuGeo_component.cfg"

cru = BmiCruAKtempMethod()
kug = BmiKuMethod()

cru.initialize(cru_config_file)
kug.initialize(kug_config_file)

T_air = cru.get_value("atmosphere_bottom_air__temperature_year")

kug.set_value('atmosphere_bottom_air__temperature', T_air)
kug.set_value('latitude', np.linspace(68,71,20))
kug.set_value('longitude', np.linspace(-165,-150,30))
#kug.set_value('soil__active_layer_thickness', T_air)

print(kug.get_value('latitude'))
print(kug._model.lat)

kug._model.lat = np.linspace(68,71,20)
kug._model.lon = np.linspace(-165,-150,30)
kug._model.Zal = T_air * 0 -999.9
kug._model.T_air = T_air
kug.output_alt = np.zeros((kug._model.end_year-kug._model.start_year+1,20,30))* 0 -999.9
kug.output_tps = np.zeros((kug._model.end_year-kug._model.start_year+1,20,30))* 0 -999.9

print(kug.get_value('latitude'))
print(kug._model.lat)

#kug.initialize(kug_config_file,wmt_option='yes')

#kug.set_value(kug.output_alt = 

#kug.set_value('atmosphere_bottom_air__temperature', T_air)

#print(kug.get_value('atmosphere_bottom_air__temperature'))

# Set the initial values of the Jan & Jul air temperature fields
#fng.set_value('atmosphere_bottom_air__temperature_mean_jan',
#    cru.get_value('atmosphere_bottom_air__temperature_mean_jan'))

#fng.set_value('atmosphere_bottom_air__temperature_mean_jul',
#    cru.get_value('atmosphere_bottom_air__temperature_mean_jul'))

# Update the 'output' of FNGeo, but don't increment the timestep
kug.update()
#kug.update()
#kug.update()

#print(kug._model.Zal)
#print(kug.get_value('soil__active_layer_thickness'))

kug.finalize()

#print(kug.get_value("soil__active_layer_thickness"))


#
## Write out the "answer" for this timestep as a raw binary array
#fn_values = fng.get_value('frostnumber__air')
#fn_values.tofile('fn_air%s.dat' % str(fng._model._date_current.year))
#
#print("Initial values are for year:    %s (%s)" %
#      (str(fng.get_current_time()),  str(fng._model._date_current.year)))
#
## Loop through the whole model run!
#while fng.get_current_time() < fng.get_end_time():
#
#    cru.update()
#
#    fng.set_value('atmosphere_bottom_air__temperature_mean_jan',
#        cru.get_value('atmosphere_bottom_air__temperature_mean_jan'))
#    fng.set_value('atmosphere_bottom_air__temperature_mean_jul',
#        cru.get_value('atmosphere_bottom_air__temperature_mean_jul'))
#
#    fng.update()
#
#    fn_values = fng.get_value('frostnumber__air')
#    fn_values.tofile('fn_air%s.dat' % str(fng._model._date_current.year))
#
#    print("Calculated values for timestep: %s (%s)" %
#          (str(fng.get_current_time()), str(fng._model._date_current.year)))
