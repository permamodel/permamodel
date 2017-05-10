# -*- coding: utf-8 -*-
"""  Frost Number by Nelson and Outcalt 1983.
DOI: 10.2307/1551363. http://www.jstor.org/stable/1551363

*The MIT License (MIT)*
Copyright (c) 2016 permamodel
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*

"""

import os
import numpy as np
from permamodel.utils import model_input
from permamodel.components import perma_base
from .. import examples_directory
from nose.tools import assert_greater_equal


class FrostnumberMethod(perma_base.PermafrostComponent):
    """ Provides 1D Frostnumber component """
    def __init__(self):
        # These are set in open_input_files
        self.air_frost_number = -99.0
        self.surface_frost_number = -99.0
        self.stefan_frost_number = -99.0
        self.T_air_min = -99.0
        self.T_air_min_file = ""
        self.T_air_min_unit = self.dummy_file()
        self.T_air_max = -99.0
        self.T_air_max_file = ""
        self.T_air_max_unit = self.dummy_file()
        self.year = -1
        self.start_year = -1
        self.end_year = -1
        self.dt = 1
        self.output = {}
        self.T_average = -99.0
        self.T_amplitude = -99.0
        self.ddt = []
        self.ddf = []
        self.h_snow = -99.0
        self.c_snow = -99.0
        self.Fplus = -99.0
        self.Twplus = -99.0
        self.Zfplus = -99.0
        self.Z_tot = -99.0
        self.stefan_number = -99.0
        self.fn_out_filename = ""

        # Not sure why these aren't set elsewhere
        self.DEBUG = True
        self.SILENT = True
        self.cfg_prefix = ""
        self.T_air = -1.0
        self.A_air = -1.0
        self.lat = -999
        self.lon = -999
        self.rho_snow = 0.0
        self.vwc_H2O = 0.0
        self.Hvgf = 0.0
        self.Hvgt = 0.0
        self.Dvf = 0.0
        self.Dvt = 0.0

    def dummy_file(self, instring=""):
        """ dummy file class, so can declare empty variable in __init__()
            but not receive "no close() member error by pylint
        """
        self.description = "this is a dummy class"

        def close():
            """ This is the method we need for a dummy_file to pass pylint
            """
            pass

        if instring == "":
            close()

    def open_input_files(self):
        """ Overrides open_input_files() from perma_base """
        self.T_air_min_file = os.path.join(examples_directory,
                                           'fn_t_air_min.dat')
        self.T_air_min_unit = open(self.T_air_min_file, "r")

        self.T_air_max_file = os.path.join(examples_directory,
                                           'fn_t_air_max.dat')
        self.T_air_max_unit = open(self.T_air_max_file, "r")

    def read_input_files(self):
        """ Overrides read_input_files() from perma_base """
        T_air_min = model_input.read_next_modified(self.T_air_min_unit,
                                                   self.T_air_min_type)
        if T_air_min is not None:
            self.T_air_min = T_air_min

        T_air_max = model_input.read_next_modified(self.T_air_max_unit,
                                                   self.T_air_max_type)
        if T_air_max is not None:
            self.T_air_max = T_air_max

    def initialize(self, cfg_file=None):
        """ Set starting values for frost number """
        self.initialize_from_config_file(cfg_file=cfg_file)
        self.initialize_frostnumber_component()

    def initialize_frostnumber_component(self):
        """ Set the starting values for the frostnumber component """
        self.air_frost_number = np.float32(-1.0)
        self.surface_frost_number = np.float32(-1.0)
        self.stefan_frost_number = np.float32(-1.0)

        self.year = self.start_year

        try:
            assert_greater_equal(self.end_year, self.start_year)
        except AssertionError:
            self.end_year = self.start_year

        # Create a dictionary to hold the output values
        # (this is unique to frost_number()
        self.output = {}

        # Here, we should calculate the initial values of all the frost numbers
        self.calculate_frost_numbers()

    def calculate_frost_numbers(self):
        """ Calculate frost numbers at the current timestep """
        # Calculate all the frost numbers using the current data
        self.calculate_air_frost_number()
        self.calculate_surface_frost_number()
        self.calculate_stefan_frost_number()

        # Add these frost numbers to the output dictionary
        self.output[self.year] = ("%5.3f" % self.air_frost_number,
                                  "%5.3f" % self.surface_frost_number,
                                  "%5.3f" % self.stefan_frost_number)

    def print_frost_numbers(self, year=-1):
        """ Print output to screen """
        # if year is -1, then use the current year of self
        # otherwise, use the specified year
        if year > 0:
            print("Year: %d  F_air=%5.3f  F_surface=%5.3f  F_stefan=%5.3f" %
                  (self.year,
                   self.air_frost_number, self.surface_frost_number,
                   self.stefan_frost_number))
        else:
            for year in sorted(self.output.keys()):
                print("Year: %d  output=%s" % (year, self.output[year]))

    def calculate_air_frost_number(self):
        """ Air frost number requires degree days before calculation """
        self.compute_degree_days()
        self.compute_air_frost_number()

    def calculate_surface_frost_number(self):
        """ Dummy value for surface frost number """
        # For now, a dummy value
        self.surface_frost_number = np.float32(-1.0)

    def calculate_stefan_frost_number(self):
        """ Dummy value for Stefan frost number """
        self.stefan_frost_number = np.float32(-1.0)

    def compute_degree_days(self):
        """
        From the min and max temperatures, compute degree freezing
        and thawing days as per Outcalt paper

        Input: T_hot (avg temp of warmest month)
               T_cold (avg temp of coldest month)

        Output: ddf (degree freezing days)
                ddt (degree thawing days)
        """

        # In the first test case, we used T_air_max and T_air_min
        T_cold = self.T_air_min
        T_hot = self.T_air_max

        assert_greater_equal(T_hot, T_cold)
        T_avg = (T_hot + T_cold) / 2.0

        # Note that these conditions should cover T_hot == T_cold
        if T_hot <= 0:
            # Always freezing
            # Negative sign because ddf is + and T_avg (here) is -
            ddf = -365.0 * T_avg
            ddt = 0
            L_winter = 365.0
            L_summer = 0.0
            T_average = (T_hot + T_cold) / 2.0
            T_amplitude = (T_hot - T_cold) / 2.0
        elif T_cold > 0:
            # Never freezing
            ddf = 0
            ddt = 365.0 * T_avg
            L_winter = 0.0
            L_summer = 365.0
            T_average = (T_hot + T_cold) / 2.0
            T_amplitude = (T_hot - T_cold) / 2.0

            """ this section shows how to read a series of temperatures
        elif (self.T_air_type != 'Scalar'):
            #wk = np.loadtxt('examples/temp_copy.txt', skiprows=1,unpack=False)
            temperature_filename = self.permafrost_dir +\
                "permamodel/examples/temp_copy.txt"
            wk = np.loadtxt(temperature_filename, skiprows=1,unpack=False)
            t_month = wk[:,0]
            T_month = wk[:,1]
            Th=max(T_month)
            Tc=min(T_month)
            T=(Th+Tc)/2                                      #(eqn. 2.1)
            A=(Th-Tc)/2                                      #(eqn. 2.2)
            Beta=np.arccos(-T/A)                             #(eqn. 2.3)
            Ts=T+A*np.sin(Beta/Beta)                         #(eqn. 2.4)
            Tw=T-A*np.sin(Beta/(np.pi-Beta))                 #(eqn. 2.5)
            L_summer=365*(Beta/np.pi)                        #(eqn. 2.6)
            L_winter = 365-L_summer                          #(eqn. 2.7)
            print('winter length:',Lw,'summer length:',Ls)
            ddt = Ts*L_summer                                #(eqn. 2.8)
            ddf = -Tw*L_winter                               #(eqn. 2.9)
            print Th,Tc
        """
        else:
            # Assume cosine fit for temp series
            T_average = (T_hot + T_cold) / 2.0
            T_amplitude = (T_hot - T_cold) / 2.0
            Beta = np.arccos(-T_average / T_amplitude)
            T_summer = T_average + T_amplitude * np.sin(Beta) / Beta
            T_winter = T_average - T_amplitude * np.sin(Beta) / (np.pi - Beta)
            L_summer = 365.0 * Beta / np.pi
            L_winter = 365.0 - L_summer
            ddt = T_summer * L_summer
            ddf = -T_winter * L_winter
        self.T_average = T_average
        self.T_amplitude = T_amplitude
        self.ddt = ddt
        self.ddf = ddf

    def compute_air_frost_number(self):
        """
        Calculating Reduced Air Frost Number (pages 280-281).
        The reduced frost number is close 0 for long summers
        and close to 1 for long winters.
        """
        self.air_frost_number = \
            np.sqrt(self.ddf) / \
            (np.sqrt(self.ddf) + np.sqrt(self.ddt))

    def close_input_files(self):
        """ As per topoflow, close the input files to finalize() """
        if self.T_air_min_type != 'Scalar':
            self.T_air_min_unit.close()
        if self.T_air_max_type != 'Scalar':
            self.T_air_max_unit.close()

    def write_output_to_file(self):
        """ Part of finalize, write the output to file(s) """
        # Write the output to a file
        with open(self.fn_out_filename, 'w') as f_out:
            for year in sorted(self.output.keys()):
                f_out.write("Year: %d  output=%s\n" %
                            (year, self.output[year]))
