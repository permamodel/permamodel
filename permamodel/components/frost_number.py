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
from __future__ import print_function

import os

import numpy as np

from permamodel import examples_directory
from permamodel.components import perma_base
from permamodel.utils import model_input


class FrostnumberMethod(perma_base.PermafrostComponent):
    """Provides 1D Frostnumber component"""

    def __init__(self):
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
        """dummy file class, so can declare empty variable in __init__()
        but not receive "no close() member error by pylint
        """
        self.description = "this is a dummy class"

        def close():
            """This is the method we need for a dummy_file to pass pylint"""
            pass

        if instring == "":
            close()

    def initialize(self, cfg_file=None):
        """Set starting values for frost number"""
        if cfg_file is not None:
            self.initialize_from_config_file(cfg_file)
        else:
            # Default configuration (one value)
            self.start_year = 2000
            self.end_year = 2000
            self.T_air_min_type = "scalar"
            self.T_air_max_type = "scalar"
            self.T_air_min = np.array([-20.0], dtype=np.float32)
            self.T_air_max = np.array([10.0], dtype=np.float32)
            self.fn_out_filename = "default_fn_config_outfile.dat"

        self.initialize_frostnumber_component()

    def initialize_from_config_file(self, cfg_file):
        """Use oldstyle configuration file format"""
        self._configuration = self.get_config_from_oldstyle_file(cfg_file)

        # Easy configuration values, simply passed
        self.start_year = self._configuration["start_year"]
        self.end_year = self._configuration["end_year"]
        self.fn_out_filename = self._configuration["fn_out_filename"]

        T_air_min_type = self._configuration["T_air_min_type"]
        if T_air_min_type.lower() == "scalar":
            self.T_air_min = np.array(
                [self._configuration["T_air_min"]], dtype=np.float32
            )
        else:
            in_dir = os.path.realpath(self._configuration["in_directory"])
            self.T_air_min = np.array(
                np.loadtxt(
                    os.path.join(in_dir, self._configuration["T_air_min"]),
                    skiprows=0,
                    unpack=False,
                ),
                dtype=np.float32,
            )

        T_air_max_type = self._configuration["T_air_max_type"]
        if T_air_max_type.lower() == "scalar":
            self.T_air_max = np.array(
                [self._configuration["T_air_max"]], dtype=np.float32
            )
        else:
            in_dir = os.path.realpath(self._configuration["in_directory"])
            self.T_air_max = np.array(
                np.loadtxt(
                    os.path.join(in_dir, self._configuration["T_air_max"]),
                    skiprows=0,
                    unpack=False,
                ),
                dtype=np.float32,
            )

    def initialize_frostnumber_component(self):
        """Set the starting values for the frostnumber component"""
        self.air_frost_number = np.float32(-1.0)
        self.surface_frost_number = np.float32(-1.0)
        self.stefan_frost_number = np.float32(-1.0)

        self.year = self.start_year

        if self.start_year > self.end_year:
            self.end_year = self.start_year

        # Create a dictionary to hold the output values
        # (this is unique to frost_number()
        self.output = {}

        # Here, we should calculate the initial values of all the frost numbers
        self.calculate_frost_numbers()

    def calculate_frost_numbers(self):
        """Calculate frost numbers at the current timestep"""
        # Calculate all the frost numbers using the current data
        self.calculate_air_frost_number()
        self.calculate_surface_frost_number()
        self.calculate_stefan_frost_number()

        # Add these frost numbers to the output dictionary
        self.output[self.year] = (
            "%5.3f" % self.air_frost_number,
            "%5.3f" % self.surface_frost_number,
            "%5.3f" % self.stefan_frost_number,
        )

    def print_frost_numbers(self, year=-1):
        """Print output to screen"""
        # if year is -1, then use the current year of self
        # otherwise, use the specified year
        if year > 0:
            print(
                (
                    "Year: %d  F_air=%5.3f  F_surface=%5.3f  F_stefan=%5.3f"
                    % (
                        self.year,
                        self.air_frost_number,
                        self.surface_frost_number,
                        self.stefan_frost_number,
                    )
                )
            )
        else:
            for year in sorted(self.output.keys()):
                print(("Year: %d  output=%s" % (year, self.output[year])))

    def calculate_air_frost_number(self):
        """Air frost number requires degree days before calculation"""
        self.compute_degree_days()
        self.compute_air_frost_number()

    def calculate_surface_frost_number(self):
        """Dummy value for surface frost number"""
        # For now, a dummy value
        self.surface_frost_number = np.float32(-1.0)

    def calculate_stefan_frost_number(self):
        """Dummy value for Stefan frost number"""
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
        if len(self.T_air_min) == 1:
            T_cold = self.T_air_min[0]
        else:
            T_cold = self.T_air_min[int(self.year - self.start_year)]
        if len(self.T_air_max) == 1:
            T_hot = self.T_air_max[0]
        else:
            T_hot = self.T_air_max[int(self.year - self.start_year)]

        assert T_hot >= T_cold
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
        self.air_frost_number = np.sqrt(self.ddf) / (
            np.sqrt(self.ddf) + np.sqrt(self.ddt)
        )

    def close_input_files(self):
        """As per topoflow, close the input files to finalize()"""
        if self.T_air_min_type != "Scalar":
            self.T_air_min_unit.close()
        if self.T_air_max_type != "Scalar":
            self.T_air_max_unit.close()

    def write_output_to_file(self):
        """Part of finalize, write the output to file(s)"""
        # Write the output to a file
        # If file is written, return value is True
        # If permission is denied, return value is False
        try:
            with open(self.fn_out_filename, "w") as f_out:
                for year in sorted(self.output.keys()):
                    f_out.write("Year: %d  output=%s\n" % (year, self.output[year]))
        except IOError:
            print(
                ("WARNING: Unable to write output to {}".format(self.fn_out_filename))
            )
            return False

        return True

    def get_config_from_oldstyle_file(self, cfg_filename):
        """Modified from that in _Geo code"""
        cfg_struct = {}
        try:
            with open(cfg_filename, "r") as cfg_file:
                # this was originally modeled after read_config_file()
                # in BMI_base.py as modified for cruAKtemp.py
                while True:
                    # Read lines from config file until no more remain
                    line = cfg_file.readline()
                    if line == "":
                        break

                    # Comments start with '#'
                    COMMENT = line[0] == "#"

                    words = line.split("|")
                    if (len(words) == 4) and (not COMMENT):
                        var_name = words[0].strip()
                        value = words[1].strip()
                        var_type = words[2].strip()

                        # Process the variables based on variable name
                        if var_name[-4:] == "date":
                            # date variables end with "_date"
                            cfg_struct[var_name] = datetime.datetime.strptime(
                                value, "%Y-%m-%d"
                            ).date()
                            # datetime.datetime.strptime(value, "%Y-%m-%d")
                        elif var_type == "int":
                            # Convert integers to int
                            cfg_struct[var_name] = int(value)
                        elif var_type == "long":
                            # Convert longs to int
                            cfg_struct[var_name] = int(value)
                        elif var_type == "float":
                            # Convert integers to float
                            cfg_struct[var_name] = float(value)
                        else:
                            # Everything else is just passed as a string
                            assert var_type == "string"
                            cfg_struct[var_name] = value

        except:
            print(
                "\nError opening configuration file in\
                  get_config_from_oldstyle_file()"
            )
            raise

        return cfg_struct

    def update(self, frac=None):
        """Move to the next timestep and update calculations"""
        if frac is not None:
            # print(
            #     "Fractional times not yet permitted, rounding to nearest int")
            time_change = self.dt * int(frac + 0.5)
        else:
            time_change = self.dt

        self.year += time_change
        self.calculate_frost_numbers()


if __name__ == "__main__":
    # Run the code
    fn = FrostnumberMethod()
    fn.initialize(cfg_file="../examples/Frostnumber_example_timeseries.cfg")
    fn.update()
    fn.update()
    fn.update()
    fn.write_output_to_file()
