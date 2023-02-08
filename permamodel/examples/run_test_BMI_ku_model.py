#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
A simple BMI KuModel example. Run with:

  $ python -m permamodel.tests.run_test_BMI_ku_model

Created on Tue Jan 10 10:56:16 2017
@author: kangwang
"""
from __future__ import print_function

import os

from permamodel import examples_directory
from permamodel.components import bmi_Ku_component

cfg_file = os.path.join(examples_directory, "Ku_method.cfg")
x = bmi_Ku_component.BmiKuMethod()

x.initialize(cfg_file)
x.update()
x.finalize()

# print x._values["ALT"][:]
print(x.get_value("soil__active_layer_thickness"))
