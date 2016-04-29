#
#  Copyright (c) 2001-2014, Scott D. Peckham
#
#
#-----------------------------------------------------------------------
#  NOTES:  This file defines a "base class" for permafrost
#          components as well as functions used by most or
#          all permafrost methods.  The methods of this class
#          should be over-ridden as necessary for different
#          methods of modeling snowmelt.

#          update_snow_vars() in precip.py sets values here.
#-----------------------------------------------------------------------
#
#  class snow_component    (inherits from BMI_base.py)
#
#      set_constants()
#      -----------------------
#      initialize()
#      update()
#      finalize()
#      set_computed_input_vars()
#      --------------------------
#      check_input_types()
#      initialize_computed_vars()
#      ----------------------------
#      update_meltrate()
#      enforce_max_meltrate()
#      update_SM_integral()
#      update_swe()
#      update_depth()
#      -----------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#      -----------------------
#      update_outfile_names()
#      open_output_files()
#      write_output_files()
#      close_output_files()
#      save_grids()
#      save_pixel_values()
#
#-----------------------------------------------------------------------

import numpy as np
import os

from permamodel.utils import BMI_base
from permamodel.utils import model_input
#from permamodel.utils import model_output



#-----------------------------------------------------------------------
class permafrost_component( BMI_base.BMI_component ):


    #------------------------------------------------------------
    # Notes: rho_H2O, Cp_snow, rho_air and Cp_air are currently
    #        hardwired (not adjustable with GUI).
    #------------------------------------------------------------
    def set_constants(self):

        #-----------------------------------
        # Constants not changeable by user
        #---------------------------------------------------------
        # Cp_snow = mass-specific isobaric heat capacity of snow
        #           value from:  NCAR CSM Flux Coupler web page
        #---------------------------------------------------------
        # Lf = latent heat of fusion for water [J kg -1]
        #---------------------------------------------------------
        self.Cp_snow  = np.float64( 2090.0 )  # [J kg-1 K-1]
        self.Lf       = np.float64( 334000 )  # [J kg-1]

        #--------------------------------------
        # Not a constant; read from CFG file.
        #--------------------------------------
        ## self.rho_snow = np.float64(300)
        ## self.rho_H2O  = np.float64(1000)  # (See initialize() method.)

    #   set_constants()
    #-------------------------------------------------------------------

    def open_input_files(self):

        #------------------------------------------------------
        # Each component that inherits from snow_base.py must
        # implement its own versions of these.
        #------------------------------------------------------
        print 'ERROR: open_input_files() for permafrost component'
        print '       has not been implemented.'
        print ' Ignoring this message for now.'

    #   open_input_files()
    #-------------------------------------------------------------------

    def extract_grid_value_from_GSD(self,input_file,
                            lonname, lon_grid_scale,
                            latname, lat_grid_scale,
                            varname):
        """
        The function is to extract the grid value from NetCDF file,
        according to input of latitude and longitude;

        INPUTs:
            input_lat: Latitude;
            input_lon: Longitude;
            input_file: grid data file (NetCDF file)
            lonname: name of longitude in "input_file"
            latname: name of latitude in "input_file"
            lon_grid_scale: grid size of longitude
            lat_grid_scale: grid size of latitude
            varname: name of variable should be extracted.

        OUTPUTs:
            p_data: grid value
        """
        print input_file

        from netCDF4 import Dataset
        #import numpy as np

        # Read the nc file

        fh = Dataset(input_file, mode='r')

        # Get the lat and lon
        #   Set the grid size for lat. and lon. (here is 0.5 degree)

        lon_grid = fh.variables[lonname][:];
        lat_grid = fh.variables[latname][:];

        # Get boundary of each grid
        #    Including top and bottom of latitude
        #              top and bottom of longitude

        lon_grid_top = lon_grid + lon_grid_scale / 2.0;
        lat_grid_top = lat_grid + lat_grid_scale / 2.0;

        lon_grid_bot = lon_grid - lon_grid_scale / 2.0;
        lat_grid_bot = lat_grid - lat_grid_scale / 2.0;

        # Get the index of input location acccording to lat and lon inputed
        idx_lon = np.where((self.lon <= lon_grid_top) & (self.lon > lon_grid_bot))
        idx_lat = np.where((self.lat <= lat_grid_top) & (self.lat > lat_grid_bot))

        idx_lon = np.array(idx_lon)
        idx_lat = np.array(idx_lat)

        p_data  = fh.variables[varname][idx_lat[0,0], idx_lon[0,0]]

        #fh.close()

        return p_data
    #   extract_grid_value_from_GSD()
    #-------------------------------------------------------------------

    def initialize_soil_texture_from_GSD(self):

        # NOTE: this part is a hardcoded
        # Maybe there is abetter way of organizing it

        Clay_file = './Parameters/T_CLAY.nc4'
        Sand_file = './Parameters/T_SAND.nc4'
        Silt_file = './Parameters/T_SILT.nc4'
        Peat_file = './Parameters/T_OC.nc4'

        # Kang please add a file check method here
        # to check that all files do exist
        # if not it should warn user that files are not found

        lonname     = 'lon'; lon_grid_scale = 0.05;
        latname     = 'lat'; lat_grid_scale = 0.05;

        self.p_clay = self.extract_grid_value_from_GSD(Clay_file,
                         lonname, lon_grid_scale,
                         latname, lat_grid_scale,
                         'T_CLAY')
        self.p_sand = self.extract_grid_value_from_GSD(Sand_file,
                         lonname, lon_grid_scale,
                         latname, lat_grid_scale,
                         'T_SAND')
        self.p_peat = self.extract_grid_value_from_GSD(Peat_file,
                         lonname, lon_grid_scale,
                         latname, lat_grid_scale,
                         'T_OC')
        self.p_silt = self.extract_grid_value_from_GSD(Silt_file,
                         lonname, lon_grid_scale,
                         latname, lat_grid_scale,
                         'T_SILT')
        self.p_peat=0
    #   initialize_soil_texture_from_GSD()
    #-------------------------------------------------------------------

    def read_input_files(self):

        print 'ERROR: read_input_files() for permafrost component'
        print '       has not been implemented.'
        print ' Ignoring this message for now. Can not find the rti class'

    #   read_input_files()
    #-------------------------------------------------------------------
    def check_input_types(self):

        #----------------------------------------------------
        # Notes: Usually this will be overridden by a given
        #        method of computing snow meltrate.
        #----------------------------------------------------
        # Notes: rho_H2O, Cp_snow, rho_air and Cp_air are
        #        currently always scalars.
        #----------------------------------------------------
        are_scalars = np.array([
                          self.is_scalar('lat'),
                          self.is_scalar('lon'),
                          self.is_scalar('T_air'),
                          self.is_scalar('A_air'),
                          self.is_scalar('h_snow'),
                          self.is_scalar('rho_snow'),
                          #----------------------------------
                          self.is_scalar('vwc_H2O'),
                          self.is_scalar('Hvgf'),
                          self.is_scalar('Hvgt'),
                          self.is_scalar('Dvf'),
                          self.is_scalar('Dvt') ])

        self.ALL_SCALARS = np.all(are_scalars)

    #   check_input_types()
    #-------------------------------------------------------------------

    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        #---------------------------------------------------------
        # Notes:  Need to make sure than h_swe matches h_snow ?
        #         User may have entered incompatible values.
        #---------------------------------------------------------
        # (3/14/07) If the Energy Balance method is used for ET,
        # then we must initialize and track snow depth even if
        # there is no snowmelt method because the snow depth
        # affects the ET rate.  Otherwise, return to caller.
        #---------------------------------------------------------
        if not(SILENT):
            print ' '
            print 'Permafrost component: Initializing...'

        self.status     = 'initializing'  # (OpenMI 2.0 convention)
        self.mode       = mode
        self.cfg_file   = cfg_file
        #print mode, cfg_file

        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        self.set_constants()
        self.initialize_config_vars()
        # At this stage we are going to ignore read_grid_info b/c
        # we do not have rti file associated with our model
        # we also skipping the basin_vars which calls the outlets
        #self.read_grid_info()
        #self.initialize_basin_vars()

        #---------------------------------------------
        # Extract soil texture from Grid Soil Database (Netcdf files)
        # according to locations
        #---------------------------------------------
        self.initialize_soil_texture_from_GSD()

        self.initialize_time_vars()

        if (self.comp_status == 'Disabled'):
            #########################################
            #  DOUBLE CHECK THIS; SEE NOTES ABOVE
            #########################################
               ####### and (ep.method != 2):  ??????
            if not(SILENT):
                print 'Permafrost component: Disabled.'
            self.lat    = self.initialize_scalar(0, dtype='float64')
            self.lon    = self.initialize_scalar(0, dtype='float64')
            self.T_air  = self.initialize_scalar(0, dtype='float64')
            self.h_snow = self.initialize_scalar(0, dtype='float64')
            self.vwc_H2O= self.initialize_scalar(0, dtype='float64')
            self.Hvgf   = self.initialize_scalar(0, dtype='float64')
            self.Hvgt   = self.initialize_scalar(0, dtype='float64')
            self.Dvf    = self.initialize_scalar(0, dtype='float64')
            self.Dvt    = self.initialize_scalar(0, dtype='float64')
            self.DONE   = True
            self.status = 'initialized'
            return

        #---------------------------------------------
        # Open input files needed to initialize vars
        #---------------------------------------------
        self.open_input_files()
        self.read_input_files()

        #---------------------------
        # Initialize computed vars
        #---------------------------
        self.check_input_types()  # (maybe not used yet)

        self.status = 'initialized'

    #   initialize()
    #-------------------------------------------------------------------
    def update_ALT(self):

        #---------------------------------------------------------
        # Note: We don't need to update any variables if
        #       the soil_thermal_conductivity method is None.  But we need
        #       to make sure that thermal_conductivity = 0.0.
        #       This "method" will be over-ridden by a
        #       particular snowmelt method.
        #--------------------------------------------------
        print 'ERROR: update_ALT() method for Permafrost component'
        print '       has not been implemented.'

    #   update_ALT()
    #-------------------------------------------------------------------
    def update_ground_temperatures(self):

       #---------------------------------------------------------
        # Note: We don't need to update any variables if
        #       the soil_thermal_conductivity method is None.  But we need
        #       to make sure that thermal_conductivity = 0.0.
        #       This "method" will be over-ridden by a
        #       particular snowmelt method.
        #--------------------------------------------------
        print 'ERROR: update_ground_temperatures method for Permafrost component'
        print '       has not been implemented.'

    #   update_ground_temperatures()
    #-------------------------------------------------------------------

    ## def update(self, dt=-1.0, time_seconds=None):
    def update(self, dt=-1.0):

        #----------------------------------------------------------
        # Note: The read_input_files() method is first called by
        #       the initialize() method.  Then, the update()
        #       method is called one or more times, and it calls
        #       other update_*() methods to compute additional
        #       variables using input data that was last read.
        #       Based on this pattern, read_input_files() should
        #       be called at end of update() method as done here.
        #       If the input files don't contain any additional
        #       data, the last data read persists by default.
        #----------------------------------------------------------

        #-------------------------------------------------
        # Note: self.SM already set to 0 by initialize()
        #-------------------------------------------------
        if (self.comp_status == 'Disabled'): return
        self.status = 'updating'  # (OpenMI)

        #-------------------------
        # Update computed values
        #-------------------------
        self.update_ground_temperatures()
        self.update_ALT()

        #-----------------------------------------
        # Read next perm vars from input files ? NOTE: does not work see the read_input_files()
        #-------------------------------------------
        # Note that read_input_files() is called
        # by initialize() and these values must be
        # used for "update" calls before reading
        # new ones.
        #-------------------------------------------
        if (self.time_index > 0):
            self.read_input_files()

        #----------------------------------------------
        # Write user-specified data to output files ?
        #----------------------------------------------
        # Components use own self.time_sec by default.
        #-----------------------------------------------
        self.write_output_files()

        #-----------------------------
        # Update internal clock
        # after write_output_files()
        #-----------------------------
        self.update_time( dt )
        self.status = 'updated'  # (OpenMI)

    #   update()
    #-------------------------------------------------------------------

    def write_output_files(self, time_seconds=None):

        #-----------------------------------------
        # Allows time to be passed from a caller
        #-----------------------------------------
        if (time_seconds is None):
            time_seconds = self.time_sec
        model_time = int(time_seconds)

        #----------------------------------------
        # Save computed values at sampled times
        #----------------------------------------
        if (model_time % int(self.save_grid_dt) == 0):
            self.save_grids()
        if (model_time % int(self.save_pixels_dt) == 0):
            self.save_pixel_values()

        #----------------------------------------
        # Save computed values at sampled times
        #----------------------------------------
##        if ((self.time_index % self.grid_save_step) == 0):
##             self.save_grids()
##        if ((self.time_index % self.pixel_save_step) == 0):
##             self.save_pixel_values()

    #   write_output_files()
    #-------------------------------------------------------------------
    def save_grids(self):
        # Saves the grid values based on the prescribed ones in cfg file

        #if (self.SAVE_MR_GRIDS):
        #    model_output.add_grid( self, self.T_air, 'T_air', self.time_min )

        if (self.SAVE_HS_GRIDS):
            model_output.add_grid( self, self.h_snow, 'h_snow', self.time_min )

        #if (self.SAVE_SW_GRIDS):
        #    model_output.add_grid( self, self.Tps, 'Tps', self.time_min )

        #if (self.SAVE_CC_GRIDS):
        #    model_output.add_grid( self, self.Zal, 'Zal', self.time_min )

    #   save_grids()
    #-------------------------------------------------------------------
    def save_pixel_values(self):

        #IDs  = self.outlet_IDs  : not sure what is this
        time = self.time_min   ###

        #if (self.SAVE_MR_PIXELS):
        #    model_output.add_values_at_IDs( self, time, self.SM, 'mr', IDs )

        if (self.SAVE_HS_PIXELS):
            model_output.add_values_at_IDs( self, time, self.h_snow, 'h_snow', IDs )

        #if (self.SAVE_SW_PIXELS):
        #    model_output.add_values_at_IDs( self, time, self.h_swe, 'sw', IDs )

        #if (self.SAVE_CC_PIXELS):
        #    model_output.add_values_at_IDs( self, time, self.Ecc, 'cc', IDs )

    #   save_pixel_values()
    #-------------------------------------------------------------------
    def finalize(self):

        self.status = 'finalizing'  # (OpenMI)
        self.close_input_files()   ##  TopoFlow input "data streams"
        self.close_output_files()
        self.status = 'finalized'  # (OpenMI)

        self.print_final_report(comp_name='Permafrost component')

        #---------------------------
        # Release all of the ports
        #----------------------------------------
        # Make this call in "finalize()" method
        # of the component's CCA Imple file
        #----------------------------------------
        # self.release_cca_ports( port_names, d_services )

    #   finalize()
    #-------------------------------------------------------------------
    def close_output_files(self):

        #if (self.SAVE_MR_GRIDS): model_output.close_gs_file( self, 'mr')
        if (self.SAVE_HS_GRIDS): model_output.close_gs_file( self, 'hs')
        #if (self.SAVE_SW_GRIDS): model_output.close_gs_file( self, 'sw')
        #if (self.SAVE_CC_GRIDS): model_output.close_gs_file( self, 'cc')
        #-----------------------------------------------------------------
        #if (self.SAVE_MR_PIXELS): model_output.close_ts_file( self, 'mr')
        #if (self.SAVE_HS_PIXELS): model_output.close_ts_file( self, 'hs')
        #if (self.SAVE_SW_PIXELS): model_output.close_ts_file( self, 'sw')
        #if (self.SAVE_CC_PIXELS): model_output.close_ts_file( self, 'cc')
    #-------------------------------------------------------------------
    def close_input_files(self):

        print 'ERROR: close_input_files() for Snow component'
        print '       has not been implemented.'

    #   close_input_files()
    #-------------------------------------------------------------------
