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
from permamodel import examples_directory
from nose.tools import assert_greater_equal, assert_true, assert_equal
from .. import data_directory

class FrostnumberMethod(perma_base.PermafrostComponent):
    """ Provides 1D Frostnumber component """
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
        self.Csn    = -99.0
        self.Fplus = -99.0
        self.Twplus = -99.0
        self.T_winter_plus = -99.0
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
        self.rho_snow = -99.0
        self.vwc_H2O = 0.0
        self.Hvgf = 0.0
        self.Hvgt = 0.0
        self.Dvf = 0.0
        self.Dvt = 0.0
        
        self.surface_frost_number_on = False
        self.stefan_frost_number_on  = False
        self.sec_per_year            = 365.0*24.0*3600.

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

    def initialize(self, cfg_file=None):
        """ Set starting values for frost number """
        if cfg_file is not None:
            self.initialize_from_config_file(cfg_file)
        else:
            # Default configuration (one value)
            self.start_year = 2000
            self.end_year = 2000
            self.T_air_min_type = 'scalar'
            self.T_air_max_type = 'scalar'
            self.T_air_min = np.array([-20.0], dtype=np.float32)
            self.T_air_max = np.array([10.0], dtype=np.float32)
            self.fn_out_filename = "default_fn_config_outfile.dat"

        self.initialize_frostnumber_component()

    def initialize_from_config_file(self, cfg_file):
        """ Use oldstyle configuration file format """
        self._configuration = self.get_config_from_oldstyle_file(cfg_file)

        # Easy configuration values, simply passed
        self.start_year = self._configuration['start_year']
        self.end_year = self._configuration['end_year']
        self.fn_out_filename = self._configuration['fn_out_filename']
        
        # Snow cover:
        try:
            self.h_snow   = self._configuration['h_snow']
            self.rho_snow = self._configuration['rho_snow']
        except:
            self.h_snow    = -99
            self.rho_snow  = -99
            
        # Site location:
        try:
            self.lat      = self._configuration['lat']
            self.lon      = self._configuration['lon']
        except:
            self.lat    = -999
            self.lon    = -999 
        
        # Soil water content:
        try:
            self.vwc_H2O      = self._configuration['vwc_H2O']
        except:
            self.vwc_H2O      = 0.0

        # These don't need to be used after this routine
        T_air_min_type = self._configuration['T_air_min_type']
        T_air_max_type = self._configuration['T_air_max_type']

        if self.start_year == self.end_year:
            # Only one year specified, so inputs should be scalar
            assert_equal(T_air_min_type.lower(), 'scalar')
            assert_equal(T_air_max_type.lower(), 'scalar')
            self.T_air_min = np.array([self._configuration['T_air_min']],
                                 dtype=np.float32)
            self.T_air_max = np.array([self._configuration['T_air_max']],
                                 dtype=np.float32)
        else:
            # Several years specified, should be timesteps
            in_dir = os.path.realpath(self._configuration['in_directory'])

            fname = os.path.join(in_dir,
                                 self._configuration['T_air_min'])
            assert_true(os.path.isfile(fname))
            Tvalues = np.loadtxt(fname, skiprows=0, unpack=False)
            self.T_air_min = np.array(Tvalues, dtype=np.float32)

            fname = os.path.join(in_dir,
                                 self._configuration['T_air_max'])
            assert_true(os.path.isfile(fname))
            Tvalues = np.loadtxt(fname, skiprows=0, unpack=False)
            self.T_air_max = np.array(Tvalues, dtype=np.float32)
            
        # Launch surface FN when there are inputs of snow:    
        if ((self.h_snow > 0) & (self.rho_snow > 0)):
            self.surface_frost_number_on = True
        
        # Launch stefan FN when there are inputs of lat-lon:
        if ((self.lat > -999) & (self.lon > -999) & (self.vwc_H2O > 0)):
            
            # if there are inputs of snow,
            # launch both. 
            if ((self.h_snow > 0) & (self.rho_snow > 0)):
                self.surface_frost_number_on = True    
                self.stefan_frost_number_on  = True
            else:
            # if there are no inputs of snow,
            # launch both, but set snow depth to zero.
                print ('No snow inputs, set snow to zero')
                self.h_snow   = 0.
                self.rho_snow = 240.
                self.surface_frost_number_on = True    
                self.stefan_frost_number_on  = True
        else:
            
            print ('Warning: The model needs soil water content, latitude, and longitude for Stefan frost number')
            print ('Otherwise, Stefan FN will be skipped!')

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
        
    def estimate_snow_damping(self):

        #--------------------------------------------------
        
        rho_sn=self.rho_snow

        self.Ksn = rho_sn**3 * 2.2E-9 + rho_sn * 4.2E-4 + 2.1E-2; # Unit: (W m-1 C-1)

        self.Csn = 2115 + 7.79 * self.T_winter;
        
        self.alpha_sn = self.Ksn / (self.Csn * rho_sn)
        
        self.Z_sn_star = np.sqrt(self.alpha_sn * self.sec_per_year/np.pi ) 
        
        self.A_plus = self.T_amplitude * np.exp(-1.*self.h_snow / self.Z_sn_star)
        
        self.T_winter_plus = self.T_average - self.A_plus *  \
                             np.sin(self.Beta) / (np.pi - self.Beta)
        
        self.ddf_plus = -self.T_winter_plus * self.L_winter
        
    def obtain_soil_parameters(self):
        
        self.read_whole_soil_texture_from_GSD() 
        self.adjusting_soil_texture()
        self.calculate_soil_bulk_density()
        self.calculate_soil_thermal_conductivity()
        

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
        
        if (self.surface_frost_number_on == True):
            self.estimate_snow_damping()
            self.surface_frost_number = np.sqrt(self.ddf_plus) /\
                          (np.sqrt(self.ddf_plus) + np.sqrt(self.ddt))

    def calculate_stefan_frost_number(self):
        
        """ Dummy value for Stefan frost number """
        
        self.stefan_frost_number = np.float32(-1.0)
        
        self.obtain_soil_parameters()
        
        self.Lf = 333E3
        
        self.Z_Frost_Plus = np.sqrt(2.0 * self.Kf * self.sec_per_year * \
                                    np.abs(self.T_winter_plus) * self.L_winter \
                                    /((0.05) * 1000.* self.Bulk_Density * self.Lf)) 
        
        self.Z_Thaw_Plus = np.sqrt(2.0 * self.Kt * self.sec_per_year * \
                                    np.abs(self.T_summer) * self.L_summer \
                                    /(self.vwc_H2O * 1000. * self.Bulk_Density* self.Lf)) 
        
        self.stefan_frost_number = self.Z_Frost_Plus/ (self.Z_Frost_Plus + self.Z_Thaw_Plus)
        
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
        T_cold = self.T_air_min[int(self.year - self.start_year)]
        T_hot  = self.T_air_max[int(self.year - self.start_year)]

        assert_greater_equal(T_hot, T_cold)
        T_avg = (T_hot + T_cold) / 2.0
        T_winter = -99.
        T_summer = -99.
        Beta     = -99.

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
        
        self.T_winter = T_winter
        self.Beta     = Beta
        self.L_winter = L_winter
        
        self.T_summer = T_summer
        self.L_summer = L_summer

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
        # If file is written, return value is True
        # If permission is denied, return value is False
        try:
            with open(self.fn_out_filename, 'w') as f_out:
                for year in sorted(self.output.keys()):
                    f_out.write("Year: %d  output=%s\n" %
                                (year, self.output[year]))
        except IOError:
            print('WARNING: Unable to write output to {}'.format(
                self.fn_out_filename))
            return False

        return True

    def get_config_from_oldstyle_file(self, cfg_filename):
        """ Modified from that in _Geo code """
        cfg_struct = {}
        try:
            with open(cfg_filename, 'r') as cfg_file:
                # this was originally modeled after read_config_file()
                # in BMI_base.py as modified for cruAKtemp.py
                while True:
                    # Read lines from config file until no more remain
                    line = cfg_file.readline()
                    if line == "":
                        break

                    # Comments start with '#'
                    COMMENT = (line[0] == '#')

                    words = line.split('|')
                    if (len(words) ==4) and (not COMMENT):
                        var_name = words[0].strip()
                        value = words[1].strip()
                        var_type = words[2].strip()

                        # Process the variables based on variable name
                        if var_name[-4:] == 'date':
                            # date variables end with "_date"
                            cfg_struct[var_name] = \
                                datetime.datetime.strptime(
                                    value, "%Y-%m-%d").date()
                                #datetime.datetime.strptime(value, "%Y-%m-%d")
                        elif var_type == 'int':
                            # Convert integers to int
                            cfg_struct[var_name] = int(value)
                        elif var_type == 'long':
                            # Convert longs to int
                            cfg_struct[var_name] = int(value)
                        elif var_type == 'float':
                            # Convert integers to float
                            cfg_struct[var_name] = float(value)
                        else:
                            # Everything else is just passed as a string
                            assert_equal(var_type, 'string')
                            cfg_struct[var_name] = value

        except:
            print("\nError opening configuration file in\
                  get_config_from_oldstyle_file()")
            raise

        return cfg_struct

    def update(self, frac=None):
        """ Move to the next timestep and update calculations """
        if frac is not None:
            # print(
            #     "Fractional times not yet permitted, rounding to nearest int")
            time_change = self.dt * int(frac + 0.5)
        else:
            time_change = self.dt

        if (self.year + time_change) > self.end_year:
            print("Update would have incremented past last year")
            print("So setting to end year")
            self.year = self.end_year
        else:
            self.year += time_change

        self.calculate_frost_numbers()
        
# =============================================================================
# ===============functions for soils===========================================
#==============================================================================

    def adjusting_soil_texture(self):
        
        self.Extract_Soil_Texture_Loops()
        
        # Adjusting percent of sand, silt, clay and peat ==
        self.tot_percent = self.p_sand+self.p_clay+self.p_silt+self.p_peat

        self.percent_sand = self.p_sand / self.tot_percent
        self.percent_clay = self.p_clay / self.tot_percent
        self.percent_silt = self.p_silt / self.tot_percent
        self.percent_peat = self.p_peat / self.tot_percent
        
    def Extract_Soil_Texture_Loops(self):
        
        n_lat = np.size(self.lat)
        n_lon = np.size(self.lon)
        
        n_grid = n_lat*n_lon     
        
        if n_grid > 1:
            
            p_clay_list = np.zeros((n_lat,n_lon));
            p_sand_list = np.zeros((n_lat,n_lon));
            p_silt_list = np.zeros((n_lat,n_lon));
            p_peat_list = np.zeros((n_lat,n_lon));
            
#            lon = np.reshape(self.lon, (n_grid,1))
#            lat = np.reshape(self.lat, (n_grid,1))
                    
            for i in range(n_lon):
                for j in range(n_lat):
                
                    input_lat   = self.lat[j]
                    input_lon   = self.lon[i]
                    
                    [p_clay0, p_sand0, p_silt0, p_peat0] = self.Extract_Soil_Texture(input_lat, input_lon);
                
                    p_clay_list[j,i] = p_clay0            
                    p_sand_list[j,i] = p_sand0        
                    p_silt_list[j,i] = p_silt0        
                    p_peat_list[j,i] = p_peat0
        else:
            
            input_lat   = self.lat
            input_lon   = self.lon
                
            [p_clay0, p_sand0, p_silt0, p_peat0] = self.Extract_Soil_Texture(input_lat, input_lon);
            
            p_clay_list = p_clay0            
            p_sand_list = p_sand0        
            p_silt_list = p_silt0        
            p_peat_list = p_peat0*0.
                         
                    
        self.p_clay = p_clay_list
        self.p_sand = p_sand_list
        self.p_silt = p_silt_list
        self.p_peat = p_peat_list
        
    def Extract_Soil_Texture(self, input_lat, input_lon): 
    
        """ 
        The function is to extract the grid value from matrix,
        according to input of latitude and longitude;
        
        INPUTs:
                input_lat: Latitude;
                input_lon: Longitude;
                lon_grid : Array of longitude
                lat_grid : Array of latitude
                p_data   : Matrix of data (from NetCDF file)
                
        OUTPUTs:
                q_data: grid value (SINGLE)   
                        
        DEPENDENTs:
                None 
        """
            
        import numpy as np
        
        lon_grid_scale = 0.05;
        lat_grid_scale = 0.05;
        
        lon_grid_top = self.lon_soil_grid + lon_grid_scale / 2.0;
        lat_grid_top = self.lat_soil_grid + lat_grid_scale / 2.0;
        
        lon_grid_bot = self.lon_soil_grid - lon_grid_scale / 2.0;
        lat_grid_bot = self.lat_soil_grid - lat_grid_scale / 2.0;
        
        # Get the index of input location acccording to lat and lon inputed
        
        idx_lon = np.where((input_lon <= lon_grid_top) & (input_lon >= lon_grid_bot))          
        idx_lat = np.where((input_lat <= lat_grid_top) & (input_lat >= lat_grid_bot))
        
        idx_lon = np.array(idx_lon)
        idx_lat = np.array(idx_lat)
        
        if np.size(idx_lon) >= 1 and np.size(idx_lat) >= 1:
            clay_perc  = self.Clay_percent[idx_lat[0,0], idx_lon[0,0]]
            sand_perc  = self.Sand_percent[idx_lat[0,0], idx_lon[0,0]]
            silt_perc  = self.Silt_percent[idx_lat[0,0], idx_lon[0,0]]
            peat_perc  = self.Peat_percent[idx_lat[0,0], idx_lon[0,0]]
        else:
            clay_perc  = np.nan;
            sand_perc  = np.nan;
            silt_perc  = np.nan;
            peat_perc  = np.nan;            
    
        return clay_perc, sand_perc, silt_perc, peat_perc        
        
    def import_ncfile(self, input_file, lonname,  latname,  varname): 
                                           
        from netCDF4 import Dataset
        
        # Read the nc file 
        
        fh = Dataset(input_file, mode='r')
        
        # Get the lat and lon
        
        lon_grid = fh.variables[lonname][:]; 
        lat_grid = fh.variables[latname][:];
        
        p_data  = fh.variables[varname][:];
        
        return lat_grid,lon_grid,p_data   
             
    def read_whole_soil_texture_from_GSD(self):
        
        Clay_file = self.get_param_nc4_filename("T_CLAY")
        Sand_file = self.get_param_nc4_filename("T_SAND")
        Silt_file = self.get_param_nc4_filename("T_SILT")
        Peat_file = self.get_param_nc4_filename("T_OC")

        lonname    = 'lon';
        latname    = 'lat';
        
        varname    = 'T_CLAY';                                        
        [lat_grid, lon_grid, Clay_percent] = self.import_ncfile(Clay_file, 
                                                lonname, latname, varname)
        varname    = 'T_SAND';                                        
        [lat_grid, lon_grid, Sand_percent] = self.import_ncfile(Sand_file, 
                                                lonname, latname, varname)
         
        varname    = 'T_SILT';                                        
        [lat_grid, lon_grid, Silt_percent] = self.import_ncfile(Silt_file, 
                                                lonname, latname, varname)
                                                
        varname    = 'T_OC';                                        
        [lat_grid, lon_grid, Peat_percent] = self.import_ncfile(Peat_file, 
                                                lonname, latname, varname)
        
        self.Clay_percent      = Clay_percent;
        self.Sand_percent      = Sand_percent;
        self.Silt_percent      = Silt_percent;
        self.Peat_percent      = Peat_percent;
        self.lon_soil_grid     = lon_grid;
        self.lat_soil_grid     = lat_grid;
        
        self.thermal_parameters_file = os.path.join(data_directory,
                                                    'Typical_Thermal_Parameters.csv')        
        self.thermal_data = np.genfromtxt(self.thermal_parameters_file,
                                          names = True,
                                          delimiter=',',
                                          dtype=None)

    def calculate_soil_bulk_density(self):
        
        Bulk_Density_Texture = self.thermal_data['Bulk_Density']
        
        self.Bulk_Density  =  Bulk_Density_Texture[2]*self.percent_clay + \
                         Bulk_Density_Texture[1]*self.percent_sand + \
                         Bulk_Density_Texture[0]*self.percent_silt + \
                         Bulk_Density_Texture[3]*self.percent_peat        # Unit: kg m-3      
        
    def calculate_soil_thermal_conductivity(self):

        #---------------------------------------------------------
        # Notes: we need a better documentation of this subroutine here
        #
        #
        #---------------------------------------------------------
        # Note: need to update frozen and thawed (kf,kt)
        #       thermal conductivities yearly
        #       this methods overriddes the method in the perma_base
        #
        #
        #--------------------------------------------------
        #input_file = 'Parameters/Typical_Thermal_Parameters.csv'

        vwc=self.vwc_H2O
                
        KT_DRY = self.thermal_data['KT_DRY'] # DRY soil thermal conductivity in THAWED states
        KT_WET = self.thermal_data['KT_WET'] # WET soil thermal conductivity in THAWED states
        KF_DRY = self.thermal_data['KF_DRY'] # DRY soil thermal conductivity in FROZEN states 
        KF_WET = self.thermal_data['KF_WET'] # WET soil thermal conductivity in FROZEN states
        
        KT_DRY = KT_DRY * 1;
        KT_WET = KT_WET * 1;
        KF_DRY = KF_DRY * 1;
        KF_WET = KF_WET * 1;
        
        kt_dry_silt = KT_DRY[0]        
        kt_dry_sand = KT_DRY[1]
        kt_dry_clay = KT_DRY[2]
        kt_dry_peat = KT_DRY[3]
        
        #===
        
        kf_dry_silt = KF_DRY[0]        
        kf_dry_sand = KF_DRY[1]
        kf_dry_clay = KF_DRY[2]
        kf_dry_peat = KF_DRY[3]

        #=== Estimate soil thermal conductivity according to water content:           
        # Estimate thermal conductivity for composed soil
                     
        Kt_Soil_dry = kt_dry_silt**self.percent_silt * \
                   kt_dry_clay**self.percent_clay * \
                   kt_dry_sand**self.percent_sand * \
                   kt_dry_peat**self.percent_peat 
                   
        uwc = 0.05;

        Kt_Soil = Kt_Soil_dry**(1.0-vwc)*0.54**vwc;

        Kf_Soil_dry = kf_dry_silt**self.percent_silt * \
                   kf_dry_clay**self.percent_clay * \
                   kf_dry_sand**self.percent_sand * \
                   kf_dry_peat**self.percent_peat

        Kf_Soil = Kf_Soil_dry**(1.0-vwc)*2.35**(vwc-uwc)*0.54**(uwc);
            

#            Kf_Soil = Kf_Soil*0.+1.38
#            Kt_Soil = Kf_Soil*0.+0.85
            
        # Consider the effect of water content on thermal conductivity        
        
        self.Kt = Kt_Soil;
        self.Kf = Kf_Soil;
                 
if __name__ == "__main__":
    # Run the code
    fn = FrostnumberMethod()
    fn.initialize(cfg_file='../examples/Frostnumber_example_timeseries.cfg')
    fn.update()
    fn.update()
    fn.update()
    fn.write_output_to_file()
