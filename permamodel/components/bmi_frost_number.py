# -*- coding: utf-8 -*-
"""  Frost Number by Nelson and Outcalt 1983. DOI: 10.2307/1551363. http://www.jstor.org/stable/1551363
"""

import numpy as np
from permamodel.utils import model_input
from permamodel.components import perma_base
from permamodel.tests import examples_directory
import os
#import gdal
#from gdalconst import *  # Import standard constants, such as GA_ReadOnly
#import osr
#from pyproj import Proj, transform

"""
class FrostnumberMethod( frost_number.BmiFrostnumberMethod ):
    _thisname = 'this name'
"""

class BmiFrostnumberMethod( perma_base.PermafrostComponent ):

    """ Implement the Nelson-Outcalt Frost numbers """

    # Set up the name of this permafrost module
    _name = 'Frost number module'

    # Indicate the CSDMS standard names of input and output variables
    _input_var_names = ('land_surface_air__temperature',
                        'land_surface__latitude',
                        'land_surface__longitude',
                       )
                       # other standard names that might apply?
                       # land_surface__temperature
                       # model__initial_time_step
                       # model__max_allowed_time_step
                       # model__min_allowed_time_step
                       # model__run_time
                       # model__spinup_time
                       # model__start_time
                       # model__stop_time
                       # model__time
                       # model__time_step
                       # model__time_step_count
    _output_var_names = ('frost_number_air',
                         'frost_number_surface',
                         'frost_number_stefan',
                        )
                        # other standard names that might apply?
                        # soil_permafrost__thickness
                        # soil_permafrost_top__depth
                        # soil_permafrost_bottom__depth

    #-------------------------------------------------------------------
    _att_map = {
    # NOTE: this will change in the future
        'model_name':         'PermaModel_frostnumber_method',
        'version':            '0.1',
        'author_name':        'J. Scott Stewart and Elchin Jafarov',
        'grid_type':          'none',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'frostnumber',
        'model_family':       'PermaModel',
        'cfg_extension':      '_frostnumber_model.cfg',
        'cmt_var_prefix':     '/input/',
        'gui_yaml_file':      '/input/frostnumber_model.yaml',
        'time_units':         'years' }

    _input_var_names = [
        'latitude',
        'longitude',
        'atmosphere_bottom_air__temperature_min',
        'atmosphere_bottom_air__temperature_max',
        'datetime__start',
        'datetime__end']

    _output_var_names = [
        'frostnumber__air',            # Air Frost number
        'frostnumber__surface',        # Surface Frost number
        'frostnumber__stefan' ]        # Stefan Frost number

    _var_name_map = {
    # NOTE: we need to look up for the corresponding standard names
        'latitude':                                  'lat',
        'longitude':                                 'lon',
        'atmosphere_bottom_air__temperature_min':    'T_air_min',
        'atmosphere_bottom_air__temperature_max':    'T_air_max',
        'datetime__start':                           'start_year',
        'datetime__end':                             'end_year',
        'frostnumber__air':                          'frostnumber_air',
        'frostnumber__surface':                      'frostnumber_surface',
        'frostnumber__stefan':                       'frostnumber_stefan'}

    _var_units_map = {
        'latitude':                                           'deg',
        'longitude':                                          'deg',
        'atmosphere_bottom_air__temperature_min':             'deg_C',
        'atmosphere_bottom_air__temperature_max':             'deg_C',
        'datetime__start':                                    'year',
        'datetime__end':                                      'year',
        'frostnumber__air':                                   '',
        'frostnumber__surface':                               '',
        'frostnumber__stefan':                                '' }

    #-------------------------------------------------------------------

    def get_attribute(self, att_name):

        try:
            return self._att_map[ att_name.lower() ]
        except:
            print '###################################################'
            print ' ERROR: Could not find attribute: ' + att_name
            print '###################################################'
            print ' '

    #   get_attribute()
    #-------------------------------------------------------------------
    def get_input_var_names(self):

        #--------------------------------------------------------
        # Note: These are currently variables needed from other
        #       components vs. those read from files or GUI.
        #--------------------------------------------------------
        return self._input_var_names

    #   get_input_var_names()
    #-------------------------------------------------------------------
    def get_output_var_names(self):

        return self._output_var_names

    #   get_output_var_names()
    #-------------------------------------------------------------------
    def get_var_name(self, long_var_name):

        return self._var_name_map[ long_var_name ]

    #   get_var_name()
    #-------------------------------------------------------------------
    def get_var_units(self, long_var_name):

        return self._var_units_map[ long_var_name ]

    #   get_var_units()
    #-------------------------------------------------------------------
    def check_input_types(self):
        """
        this functionality is not used for frostnumber_method
        """
        return

        #--------------------------------------------------
        # Notes: rho_H2O, Cp_snow, rho_air and Cp_air are
        #        currently always scalars.
        #--------------------------------------------------
        """
        are_scalars = np.array([
                          self.is_scalar('lat'),
                          self.is_scalar('lon'),
                          self.is_scalar('T_air_min'),
                          self.is_scalar('T_air_max'),
                          self.is_scalar('start_year'),
                          self.is_scalar('end_year') ])

        self.ALL_SCALARS = np.all(are_scalars)
        """

    #   check_input_types()
    #-------------------------------------------------------------------

    def get_current_time(self):
        # For the frostnumber component, the time is simply the year
        return self.year

    def update(self):
        # Ensure that we've already initialized the run
        assert(self.status == 'initialized')

        # Update the time
        self.year += self.dt

        # Get new input values
        self.read_input_files()

        # Calculate the new frost number values
        self.calculate_frost_numbers()

    def update_until(self, stop_year):
        # Ensure that stop_year is at least the current year
        if stop_year < self.year:
            print("Warning: update_until year is less than current year")
            print("  no update run")
            return

        if stop_year > self.end_year:
            print("Warning: update_until year was greater than end_year")
            print("  setting stop_year to end_year")
            stop_year = self.end_year

        # Implement the loop to update until stop_year
        year = self.year
        while year < stop_year:
            self.update()
            year = self.year

    def finalize(self, SILENT=True):
        # Finish with the run
        self.status = 'finalizing'  # (OpenMI)

        # Close the input files
        self.close_input_files()   # Close any input files

        # Write output last output
        self.write_output_to_file(SILENT=True)

        # Close the output files
        self.close_output_files()

        # Done finalizing  
        self.status = 'finalized'  # (OpenMI)

        # Print final report, as desired
        if not SILENT:
            self.print_final_report(\
                    comp_name='Permamodel FrostNumber component')

    def get_end_time(self):
        return self.end_year

    # ----------------------------------
    # Functions added to pass bmi-tester
    # ----------------------------------
    def get_grid_type(arg1, arg2):
        print("arg1: %s" % arg1)
        print("arg2: %s" % arg2)
        return 'scalar'

