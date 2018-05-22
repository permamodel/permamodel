#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
A simple BMI KuModel example. Run with:

  $ python -m permamodel.tests.run_test_BMI_ku_model

Created on Tue Jan 10 10:56:16 2017
@author: kangwang
"""

import os
from permamodel.components import bmi_Ku_component
from permamodel import examples_directory
import numpy as np


cfg_file = os.path.join(examples_directory, 'Ku_method.cfg')
x = bmi_Ku_component.BmiKuMethod()

x.initialize(cfg_file)

print(x.get_value('atmosphere_bottom_air__temperature'), x._model.T_air)

x.set_value('atmosphere_bottom_air__temperature',-5.)

print(x.get_value('atmosphere_bottom_air__temperature'), x._model.T_air)

x.update()
x.finalize()

print x.get_value('soil__active_layer_thickness')



#print x._var_name_map['atmosphere_bottom_air__temperature']

#setattr(x, x._var_name_map['atmosphere_bottom_air__temperature'],-15.)

#print x._var_name_map['atmosphere_bottom_air__temperature']
