"""
couple_cru_FNGeo.py

Example that couples cruAKtemp to FNGeo!

"""

from __future__ import print_function

from cruAKtemp.bmi_cruAKtemp import BmiCruAKtempMethod
from permamodel.components.bmi_frost_number_Geo import BmiFrostnumberGeoMethod

# Initial both components with local config files
cru_config_file = "./cruAKtemp_component.cfg"
fng_config_file = "./FNGeo_component.cfg"

cru = BmiCruAKtempMethod()
fng = BmiFrostnumberGeoMethod()

cru.initialize(cru_config_file)
fng.initialize(fng_config_file)

# Set the initial values of the Jan & Jul air temperature fields
fng.set_value('atmosphere_bottom_air__temperature_mean_jan',
    cru.get_value('atmosphere_bottom_air__temperature_mean_jan'))

fng.set_value('atmosphere_bottom_air__temperature_mean_jul',
    cru.get_value('atmosphere_bottom_air__temperature_mean_jul'))

# Update the 'output' of FNGeo, but don't increment the timestep
fng.update_frac(0)

# Write out the "answer" for this timestep as a raw binary array
fn_values = fng.get_value('frostnumber__air')
fn_values.tofile('fn_air%s.dat' % str(fng._model._date_current.year))

print("Initial values are for year:    %s (%s)" %
      (str(fng.get_current_time()),  str(fng._model._date_current.year)))

# Loop through the whole model run!
while fng.get_current_time() < fng.get_end_time():

    cru.update()

    fng.set_value('atmosphere_bottom_air__temperature_mean_jan',
        cru.get_value('atmosphere_bottom_air__temperature_mean_jan'))
    fng.set_value('atmosphere_bottom_air__temperature_mean_jul',
        cru.get_value('atmosphere_bottom_air__temperature_mean_jul'))

    fng.update()

    fn_values = fng.get_value('frostnumber__air')
    fn_values.tofile('fn_air%s.dat' % str(fng._model._date_current.year))

    print("Calculated values for timestep: %s (%s)" %
          (str(fng.get_current_time()), str(fng._model._date_current.year)))
