# -*- coding: utf-8 -*-
"""  Frost Number by Nelson and Outcalt 1983.
     DOI: 10.2307/1551363. http://www.jstor.org/stable/1551363

     This version is 2D
"""

import numpy as np
from permamodel.utils import model_input
from permamodel.components import perma_base
from permamodel.components import frost_number_2D
from permamodel.components.perma_base import *
from permamodel.tests import examples_directory
import os

class BmiFrostnumber2DMethod( perma_base.PermafrostComponent ):

    """ Implement the Nelson-Outcalt Frost numbers in 2D"""

    # Set up the name of this permafrost module
    _name = 'Frost number module, 2D version'

    #-------------------------------------------------------------------
    _att_map = {
    # NOTE: this will change in the future
        'model_name':         'PermaModel_frostnumber_2D_method',
        'version':            '0.1',
        'author_name':        'J. Scott Stewart',
        'grid_type':          'none',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'frostnumber2D',
        'model_family':       'PermaModel',
        'cfg_extension':      '_frostnumber2D_model.cfg',
        'cmt_var_prefix':     '/input/',
        'gui_yaml_file':      '/input/frostnumber2D_model.yaml',
        'time_units':         'years' }

    _input_var_names = (
        'atmosphere_bottom_air__temperature',
        )

    _output_var_names = (
        'frostnumber__air',            # Air Frost number
        'frostnumber__surface',        # Surface Frost number
        'frostnumber__stefan' )        # Stefan Frost number

    _var_name_map = {
        # These are the corresponding CSDMS standard names
        # NOTE: we need to look up for the corresponding standard names
        'atmosphere_bottom_air__temperature':        'T_air',
        'datetime__start':                           'start_year',
        'datetime__end':                             'end_year',
        'frostnumber__air':                          'air_frost_number_2D',
        'frostnumber__surface':                      'surface_frost_number_2D',
        'frostnumber__stefan':                       'stefan_frost_number_2D'}


    _var_units_map = {
        # These are the links to the model's variables' units
        'atmosphere_bottom_air__temperature':                 'deg_C',
        'frostnumber__air':                                   '1',
        'frostnumber__surface':                               '1',
        'frostnumber__stefan':                                '1' }

    #-------------------------------------------------------------------
    def __init__(self):
        self._model = None
        self._values = {}
        self._var_units = {}
        self._grids = {}
        self._grid_type = {}

    def initialize(self, cfg_file=None):
        self._model = frost_number_2D.Frostnumber2DMethod()

        self._model.initialize_from_config_file(cfg_file=cfg_file)
        self._model.initialize_frostnumber_component()

        # Set the name of this component
        self._name = "Permamodel Frostnumber2D Component"

        # Verify that all input and output variable names are in the
        # variable name and the units map
        for varname in self._input_var_names:
            assert(varname in self._var_name_map)
            assert(varname in self._var_units_map)
            #print("Input var %s is in the name map and the units map"\
            #      % varname)
        for varname in self._output_var_names:
            assert(varname in self._var_name_map)
            assert(varname in self._var_units_map)
            #print("Output var %s is in the name map and the units map"\
            #      % varname)

        # Set the Frost Number grids, based on input and output variables
        #print("Number of input variables: %d" % len(self._input_var_names))
        #print("Number of output variables: %d" % len(self._output_var_names))

        # Set the names and types of the grids
        # Note: A single value is a uniform rectilinear grid of shape (1)
        #       and size 1
        gridnumber = 0
        for varname in self._input_var_names:
            self._grids[gridnumber] = varname
            self._grid_type[gridnumber] = 'uniform_rectilinear'
            #self._grid_type[gridnumber] = 'scalar'
            gridnumber += 1
        for varname in self._output_var_names:
            self._grids[gridnumber] = varname
            self._grid_type[gridnumber] = 'uniform_rectilinear'
            #self._grid_type[gridnumber] = 'scalar'
            gridnumber += 1

        # Set the internal (frost number) variables that correspond
        # to the input and output variable names
        # Note: since we used Topoflow's _var_name_map for this, it is that
        self._values = _values = {
        # These are the links to the model's variables and
        # should be consistent with _var_name_map
            'atmosphere_bottom_air__temperature':    self._model.T_air,
            'datetime__start':          self._model.start_year,
            'datetime__end':            self._model.end_year,
            'frostnumber__air':         self._model.air_frost_number_2D,
            'frostnumber__surface':     self._model.surface_frost_number_2D,
            'frostnumber__stefan':      self._model.stefan_frost_number_2D}

        # initialize() tasks complete.  Update status.
        self.status = 'initialized'

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

    def update(self):
        # Ensure that we've already initialized the run
        assert(self._model.status == 'initialized')

        # Calculate the new frost number values
        self._model.calculate_frost_numbers_2D()
        self._values['frostnumber__air'] = self._model.air_frost_number_2D

        # Update the time
        self._model.year += self._model.dt

        # Get new input values
        self._model.read_input_files()

    def update_frac(self, time_fraction):
        # Only increment the time by a partial time step
        # Ensure that we've already initialized the run
        assert(self._model.status == 'initialized')

        # Determine which year the model is currently in
        current_model_year = int(self._model.year)

        # Update the time with a partial time step
        self._model.year += time_fraction * self._model.dt

        # Determine if the model year is now different
        new_model_year = int(self._model.year)

        # If the year has changed, change the values
        if new_model_year > current_model_year:
            # Get new input values
            self._model.read_input_files()

            # Calculate the new frost number values
            self._model.calculate_frost_numbers_2D()
            self._values['frostnumber__air'] = self._model.air_frost_number_2D

    def update_until(self, stop_year):
        # Ensure that stop_year is at least the current year
        if stop_year < self._model.year:
            print("Warning: update_until year is less than current year")
            print("  no update run")
            return

        if stop_year > self._model.end_year:
            print("Warning: update_until year was greater than end_year")
            print("  setting stop_year to end_year")
            stop_year = self.end_year

        # Implement the loop to update until stop_year
        year = self._model.year
        while year < stop_year:
            self.update()
            year = self._model.year

    def finalize(self):
        SILENT = True

        # Finish with the run
        self._model.status = 'finalizing'  # (OpenMI)

        # Close the input files
        self._model.close_input_files()   # Close any input files

        # Write output last output
        self._model.write_output_to_file(SILENT=True)

        # Close the output files
        self._model.close_output_files()

        # Done finalizing  
        self._model.status = 'finalized'  # (OpenMI)

        # Print final report, as desired
        if not SILENT:
            self._model.print_final_report(\
                    comp_name='Permamodel FrostNumber2D component')

    def get_start_time(self):
        return 0.0

    def get_current_time(self):
        return self._model.year - self._model.start_year

    def get_end_time(self):
        return self._model.end_year - self._model.start_year + 1.0

    # ----------------------------------
    # Functions added to pass bmi-tester
    # ----------------------------------
    def get_grid_type(self, grid_number):
        return self._grid_type[grid_number]

    def get_time_step(self):
        return self._model.dt

    # Note: get_value_ref() copied from bmi_heat.py
    def get_value_ref(self, var_name):
        """Reference to values.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.

        Returns
        -------
        array_like
            Value array.
        """
        return self._values[var_name]

    def set_value(self, var_name, new_var_values):
        self._values[var_name] = new_var_values

    def set_value_at_indices(self, var_name, new_var_values, indices):
        self.get_value_ref(var_name).flat[indices] = new_var_values

    def get_var_itemsize(self, var_name):
        return np.asarray(self.get_value_ref(var_name)).flatten()[0].nbytes

    def get_value_at_indices(self, var_name, indices):
        return self.get_value_ref(var_name).take(indices)

    def get_var_nbytes(self, var_name):
        return np.asarray(self.get_value_ref(var_name)).nbytes

    def get_value(self, var_name):
        """Copy of values.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.

        Returns
        -------
        array_like
            Copy of values.
        """
        # Original version: from bmi_heat.py
        #return self.get_value_ref(var_name).copy()

        # Version to convert to numpy array for bmi-tester compliance
        # Note: converting to np arrays on the fly here
        # Note: float values don't have a copy() function
        #try:
        #    return np.array(self.get_value_ref(var_name).copy())
        #except AttributeError:
        #    return np.array(self.get_value_ref(var_name))

        # This version is simpler than above, but it may break when
        #   using scalars because the resulting variable doesn't
        #   have a shape
        return np.asarray(self.get_value_ref(var_name))


    def get_var_type(self, var_name):
        """Data type of variable.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.

        Returns
        -------
        str
            Data type.
        """
        return str(self.get_value_ref(var_name).dtype)

    def get_component_name(self):
        return self._name

    # Copied from bmi_heat.py
    def get_var_grid(self, var_name):
        """Grid id for a variable.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.

        Returns
        -------
        int
            Grid id.
        """
        for grid_id, var_name_list in self._grids.items():
            if var_name in var_name_list:
                return grid_id

    def get_grid_shape(self, grid_id):
        """Number of rows and columns of uniform rectilinear grid."""
        var_name = self._grids[grid_id]
        value = np.array(self.get_value_ref(var_name)).shape
        return value

    def get_grid_size(self, grid_id):
        """Size of grid.

        Parameters
        ----------
        grid_id : int
            Identifier of a grid.

        Returns
        -------
        int
            Size of grid.

        """
        grid_size = self.get_grid_shape(grid_id)
        if grid_size == ():
            return 1
        else:
            return int(np.prod(grid_size))

    def get_grid_rank(self, var_id):
        """Rank of grid.

        Parameters
        ----------
        grid_id : int
            Identifier of a grid.

        Returns
        -------
        int
            Rank of grid.
        """
        return len(self.get_grid_shape(self.get_var_grid(var_id)))
