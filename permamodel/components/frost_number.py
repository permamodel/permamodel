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


class FrostnumberMethod(perma_base.PermafrostComponent):
    """ Provides 1D Frostnumber component """
    def __init__(self):
        # These are set in open_input_files
        self._model = ""
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
        self.r_snow = -99.0
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
        self.SILENT = False
        self.cfg_prefix = ""
        self.T_air = -1.0
        self.A_air = -1.0
        self.lat = -999
        self.lon = -999
        self.h_snow = 0.0
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

        # lat and lon not implemented yet

    def read_input_files(self):
        """ Overrides read_input_files() from perma_base """
        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------
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
        SILENT = True

        # Note: Initialized from initialize() in perma_base.py
        if not SILENT:
            print("Initializing for FrostnumberMethod")

        self._model = 'FrostNumber'

        # Here, initialize the variables which are unique to the
        # frost_number component

        # Initialize the output variables (internal names)
        self.air_frost_number = np.float32(-1.0)
        self.surface_frost_number = np.float32(-1.0)
        self.stefan_frost_number = np.float32(-1.0)

        # Initialize the year to the start year
        #  or to zero if it doesn't exist
        try:
            self.year = self.start_year
        except AttributeError:
            self.year = 0
            self.start_year = 0
            self.end_year = 0

        # Ensure that the end_year is not before the start_year
        # If no end_year is given,
        #   it is assumed that this will run for one year
        #   so the end_year is the same as the start_year
        try:
            assert self.end_year >= self.start_year
        except AttributeError:
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
        """ From the min and max temperatures, compute degree freezing
            and thawing days as per Outcalt paper
        """

        # Input: T_hot (avg temp of warmest month)
        #        T_cold (avg temp of coldest month)

        # Output: ddf (degree freezing days)
        #         ddt (degree thawing days)

        # In the first test case, we used T_air_max and T_air_min
        T_cold = self.T_air_min
        T_hot = self.T_air_max

        assert T_hot > T_cold
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

    def update_snow_prop(self):
        """ Snow properties are needed for Stefan frost number """
        # find indexes for which temp > 0 and make precip = 0
        if self.T_air_type != 'Scalar':  # if not should stop
            precipitation_filename = self.permafrost_dir +\
                "permamodel/examples/prec.txt"
            wk = np.loadtxt(precipitation_filename, skiprows=1, unpack=False)
            # t_month = wk[:, 0]
            prec_month = wk[:, 1]

        pos_temp_ind = np.array(np.where(self.ta_month > 0))
        prec_month[pos_temp_ind] = 0
        neg_temp_ind = np.array(np.where(self.ta_month <= 0))

        if not pos_temp_ind.any():
            # monthly temp is always below zero
            # i.e. it constantly snows over whole year
            # the point associated with glaciaer and needs to excluded
            print('snows constatly: stop!')

        m = np.size(neg_temp_ind)
        # assume only 50% of precip change to at the beg
        # and end of the snow season
        pp = 0.5

        # this is portions of the code assumes a perfect winter season
        # needs to be used with care when there is a warm month
        # during snow season
        if m == 1:
            s_idx = neg_temp_ind[:, 0]
            e_idx = neg_temp_ind[:, m - 1]
            prec_month[s_idx] = prec_month[s_idx] * pp
        else:
            s_idx = neg_temp_ind[:, 0]
            e_idx = neg_temp_ind[:, m - 1]
            prec_month[s_idx] = prec_month[s_idx] * pp
            prec_month[e_idx] = prec_month[e_idx] * pp

        # sum up precip to get SWE
        j = 0
        s = 0
        swe = np.zeros(m)
        for i in range(s_idx, e_idx + 1):
            s = s + prec_month[i]
            swe[j] = s
            j = j + 1

        #calculating snow density, depth and thermal counductivity
        r_snow = np.zeros(m)  # snow density in kg/m3
        h_snow = np.zeros(m)  # snow depth in m
        c_snow = np.zeros(m)  # snow depth in W/mK
        # allowed min and max snow density
        rho_sn_min = 200
        rho_sn_max = 300
        # e-folding value (see Verseghy, 1991)
        tauf = 0.24

        s = rho_sn_min
        s = ((s - rho_sn_max) * np.exp(-tauf)) + rho_sn_max
        r_snow[0] = s
        for i in range(1, m):
            # starting from month 2 tauf should be multpled by the 30 days
            # otherwise snow thermal conductivity can be low and
            # insulate ground well enough over the season
            # usually we assume constant max snow thermal conductivity
            # over snow season
            s = ((s - rho_sn_max) * np.exp(-tauf)) + rho_sn_max
            r_snow[i] = s

        h_snow = (swe / (r_snow * 0.001))
        # snow thermal conductivity according to M. Sturm, 1997.
        c_snow = (0.138 - 1.01 * r_snow + 3.233 * (r_snow**2)) * 1e-6

        self.r_snow = r_snow
        self.h_snow = h_snow
        self.c_snow = c_snow

    def update_surface_frost_number(self):
        """ re-compute surface frost number values """
        # phi [scalar]: sites latitude
        # Zs [scalar]: an average winter snow thickness
        # Zss [scalar]: a damping depth in snow
        # P [scalar]: length of an annual temperature cycle
        # k [scalar]: number of winter months
        # rho_s [scalar]: density of snow [kg m-3]
        # lambda_s [scalar]: snow thermal conductivity [W m-1 C-1]
        # c_s [scalar]: snow specific heat capacity [J kg-1 C-1]
        # alpha_s [scalar]: thermal diffusivity of snow
        # Uw [scalar]: mean winter wind speed [m s-1]
        # Aplus [scalar]: temperature amplitude at the surface with snow
        # Twplus [scalar]: the mean winter surface temperature
        # DDFplus [scalar]: freezing index at the surface
        # Tplus [scalar]: mean annual tempratures at the surface
        # Fplus [scalar] : surface frost number

        rho_s = np.mean(self.r_snow)
        lambda_s = np.mean(self.c_snow)
        Zs = np.mean(self.h_snow)
        # i am not sure what they mean by length of the annual temprature cycle
        # Something worthwhile discussing
        P = 2 * np.pi / 365

        #(eqn. 7)
        c_s = 7.79 * self.Tw + 2115
        #(eqn. 8)
        alpha_s = lambda_s / (c_s * rho_s)
        #(eqn. 10)
        Zss = np.sqrt(alpha_s * P / np.pi)
        #(eqn. 9)
        Aplus = self.A_air * np.exp(-Zs / Zss)
        #(eqn. 11)
        Twplus = self.T_air - Aplus * np.sin(self.beta / (np.pi - self.beta))
        # Twplus is a mean winter surface temprature, I think,
        # should be warmer than air temperature?
        # Here is another problem. DDFplus degree days winter
        # should be positive.
        # The way it is written in the paper is wrong. I added a minus sign
        # to fix it (see eqn. 2.9)
        #(eqn. 12)
        DDFplus = -Twplus * self.Lw
        #(eqn. 13)
        # Tplus = (self.ddt - DDFplus) / 365
        #Nevertheless the surface frost number is smaller than air
        # which looks resonable to me.
        #(eqn. 14)
        self.Fplus = np.sqrt(DDFplus) / (np.sqrt(self.ddt) + np.sqrt(DDFplus))
        self.Twplus = Twplus

    #   update_surface_frost_number()
    #-------------------------------------------------------------------
    def update_stefan_frost_number(self):
        """ Bring Stefan frost number calculation up to date """
        # Zfplus [scalar]: the depth [m] to which forst extends
        # lambda_f [scalar]: frozen soil thermal conductivity [W m-1 C-1]
        # S [scalar]: is a const scalar factor [s d-1]
        # rho_d [scalar]: dry density of soil [kg m-3]
        # wf [scalar] : soil water content (proportion of dry weight)
        # L [scalar] : is a latent heat of fusion of water [J kg-1]

        lambda_f = 1.67  # some dummy thermal conductivity
        # https://shop.bgs.ac.uk/GeoReports/examples/modules/C012.pdf
        sec_per_day = 86400
        rho_d = 2.798  # dry density of silt
        wf = 0.4       # tipical for silty soils
        denominator = rho_d * wf * self.Lf
        #(eqn. 15)
        self.Zfplus = np.sqrt(2 * lambda_f * sec_per_day *
                              np.abs(self.Twplus) * self.Lw / denominator)
        print('Zfplus = %f' % self.Zfplus)

        # assuming 3 soil layers with thickness 0.25, 0.5 and 1.75 m
        # and thermal conductivities 0.08, 1.5, and 2.1
        soil_thick = np.array([0.25, 0.5, 1.75])
        lambda_b = np.array([0.08, 1.5, 2.1])
        # resistivity R
        R = soil_thick / lambda_b
        # volumetric latent heat, 1000 s a density of water
        QL = self.Lf / 1000
        #partial freezing thawing index DD
        DD = np.zeros(3)
        Z = np.zeros(3)
        DD[0] = 0.5 * QL * soil_thick[0] * R[0] / sec_per_day
        S = 0
        for i in range(1, 3):
            S = R[i] + S
            DD[i] = (S + 0.5 * R[i - 1]) * QL * soil_thick[i] / sec_per_day
            #The depth of the frost thaw penetration
        S = 0
        Z_tot = 0
        for i in range(0, 3):
            #The depth of the frost thaw penetration
            Z[i] = np.sqrt(2 * lambda_b[i] * sec_per_day * DD[i] / QL +
                           lambda_b[i]**2 * S ** 2) - lambda_b[i] * S
            S = R[i] + S
            Z_tot = Z_tot + Z[i]

        self.Z_tot = Z_tot
        self.stefan_number = np.sqrt(self.Fplus) / \
            (np.sqrt(self.Fplus) + np.sqrt(self.Z_tot))

    def close_input_files(self):
        """ As per topoflow, close the input files to finalize() """
        if self.T_air_min_type != 'Scalar':
            self.T_air_min_unit.close()
        if self.T_air_max_type != 'Scalar':
            self.T_air_max_unit.close()

    def write_output_to_file(self, SILENT=True):
        """ Part of finalize, write the output to file(s) """
        # Write the output to the screen unless we're silent
        if not SILENT:
            self.print_frost_numbers(self.year)

        # Write the output to a file
        with open(self.fn_out_filename, 'w') as f_out:
            for year in sorted(self.output.keys()):
                f_out.write("Year: %d  output=%s\n" %
                            (year, self.output[year]))
