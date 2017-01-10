#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 10:56:16 2017

@author: kangwang
"""

import sys
sys.path.append('permamodel/')

from permamodel.components import bmi_frost_number
x=bmi_frost_number.BmiFrostnumberMethod()

x.initialize('permamodel/examples/Frostnumber_example_singlesite_singleyear.cfg')
print x.initialize.im_self.status

x.update()

x.finalize()

print x.status