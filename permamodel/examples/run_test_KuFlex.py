#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
A simple BMI KuModel example. Run with:

  $ python -m permamodel.tests.run_test_BMI_ku_model

Created on Tue Jan 10 10:56:16 2017
@author: kangwang
"""
from __future__ import print_function

import numpy as np

import os
#from permamodel.components import bmi_Ku_component
from permamodel import examples_directory
from permamodel.components import KuFlex_method

ku = KuFlex_method.KuFlex_method()

ku.initialize('KuFlex_method.cfg')

#ku.T_air = 2
#ku.A_air = 12

for i in np.arange(4):

    ku.update()
    
    #for i in dir(ku): 
    #    print (i)
    print('Tair=',ku.T_air)
    print('Aps=',ku.Aps)
    print('Tps=',ku.Tps)
    print('Tps_numerator=', ku.Tps_numerator)
    print('Zal=',ku.Zal)