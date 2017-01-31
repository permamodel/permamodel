#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 10:56:16 2017

@author: kangwang
"""

import sys


from permamodel.components import bmi_Ku_component
x=bmi_Ku_component.BmiKuMethod()

x.initialize()
x.update()
x.finalize()

print x._values["ALT"][:]
