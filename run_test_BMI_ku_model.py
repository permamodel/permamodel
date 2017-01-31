#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 10:56:16 2017

@author: kangwang
"""

import os
import sys
from permamodel.components import bmi_Ku_component
from permamodel.tests import examples_directory


cfg_file = os.path.join(examples_directory, 'Ku_method.cfg')
x = bmi_Ku_component.BmiKuMethod()

x.initialize(cfg_file)
x.update()
x.finalize()

print x._values["ALT"][:]
