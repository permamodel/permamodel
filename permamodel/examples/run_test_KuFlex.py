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
import matplotlib.pyplot as plt
from permamodel.components import bmi_KuFlex_component

KuFlex_method_on = False

kuflex = bmi_KuFlex_component.BmiKuFlexMethod()

kuflex.initialize('KuFlex_method.cfg')

#print(kuflex._model)

#for lst in kuflex.get_input_var_names():
#    
#    print(lst,':---')
#    print(kuflex.get_value(lst))
#    print(kuflex.get_var_grid(lst))
#    print(kuflex.get_grid_type(kuflex.get_var_grid(lst)))
#    print('rank=',kuflex.get_grid_rank(kuflex.get_var_grid(lst)))
#    print(kuflex.get_grid_shape(kuflex.get_var_grid(lst)))
#    print(kuflex.get_grid_size(kuflex.get_var_grid(lst)))
#    print(kuflex.get_grid_spacing(kuflex.get_var_grid(lst)))
#    print(kuflex.get_grid_origin(kuflex.get_var_grid(lst)))

for i in np.arange(3):
    
#    kuflex.set_value_at_indices('atmosphere_bottom_air__temperature', -3, 0)
#    kuflex.set_value('atmosphere_bottom_air__temperature', -3)
    print('t=',kuflex.get_current_time())
    print('MAAT=',kuflex.get_value('atmosphere_bottom_air__temperature')) 
#    print(kuflex.get_var_itemsize('atmosphere_bottom_air__temperature'))
        
    kuflex.update()
    
    for lst in kuflex.get_output_var_names():
        print(lst,':---')
        print(np.nanmin(kuflex.get_value(lst)))
        
kuflex.finalize()

#kuflex.initialize('KuFlex_method_ts.cfg')
#
#kuflex.update_until(7)
##
#kuflex.finalize()

#
if KuFlex_method_on == True:
    
    ku = KuFlex_method.KuFlex_method()
    
    ku.initialize('KuFlex_method.cfg')
        
    for i in np.arange(3):
    
        ku.update()
        
        #for i in dir(ku): 
        #    print (i)
        print('Tair=',np.nanmean(ku.T_air))
        print('Aair=',np.nanmean(ku.A_air))
        print('Aps=',np.nanmean(ku.Aps))
    #    print('Tps=',ku.Tps)
    #    print('Tps_numerator=', ku.Tps_numerator)
    #    print('Zal=',ku.Zal)
        
        plt.figure(figsize=[7,6])
        
        plt.subplot(2,1,1)
        plt.imshow(ku.Tps)
        plt.title(str(i))
        plt.colorbar(orientation='horizontal')
        cs = plt.contour(ku.Tps, [0])
        plt.clabel(cs, inline=1, fontsize=10, fmt = '%0d')
        
        plt.subplot(2,1,2)
        plt.imshow(ku.Zal)
        plt.title(str(i))
        plt.colorbar(orientation='horizontal')
        cs2 = plt.contour(ku.Zal, [0.5, 1.0, 2.0], colors = 'k')
        plt.clabel(cs2, inline=1, fontsize=10, fmt = '%0.2f')
        
        plt.show()