"""
couple_cru_FNGeo.py

Example that couples cruAKtemp to FNGeo!

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

kug._model.lat = np.linspace(68,71,20)
kug._model.lon = np.linspace(-165,-150,30)
kug._model.Zal = np.zeros((20,30)) -999.9

kug.output_alt = np.zeros((kug._model.end_year-kug._model.start_year+1,20,30))* 0 -999.9
kug.output_tps = np.zeros((kug._model.end_year-kug._model.start_year+1,20,30))* 0 -999.9

print(kug.get_value('latitude'))
print(kug._model.lat)

for i in np.arange(3):
    cru.update()
    T_air = cru.get_value("atmosphere_bottom_air__temperature_year")
    kug.set_value('atmosphere_bottom_air__temperature', T_air)
    kug.update()

kug.finalize()

print(kug.get_value("soil__active_layer_thickness"))