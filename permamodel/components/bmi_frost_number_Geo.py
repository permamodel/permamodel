# -*- coding: utf-8 -*-
"""  Frost Number by Nelson and Outcalt 1983.
     DOI: 10.2307/1551363. http://www.jstor.org/stable/1551363

     This is the Geo version
"""

import numpy as np
from permamodel.components import perma_base
from permamodel.components import frost_number_Geo
from .. import examples_directory
from nose.tools import (assert_is_instance, assert_greater_equal,
                        assert_less_equal, assert_almost_equal,
                        assert_greater, assert_in, assert_true,
                        assert_false, assert_equal, assert_raises)
import datetime

import os

class BmiFrostnumberGeoMethod( perma_base.PermafrostComponent ):

    """ Implement the Nelson-Outcalt Frost numbers
        for a geographic region"""

    # Note: The definitions of _att_map, _var_name_map, etc.
    # have been moved inside the __init__(self) function because
    # otherwise changing thme once changes them for future
    # instantiations, .e.g during testing!


    #-------------------------------------------------------------------
    def __init__(self):
        self._model = None
        self._values = {}
        self._var_units = {}
        self._grids = {}
        self._grid_type = {}
        # Set up the name of this permafrost module
        self._name = 'Frost number module, Geo version'

        #-------------------------------------------------------------------
        self._att_map = {
        # NOTE: this will change in the future
            'model_name':         'PermaModel_frostnumber_Geo_method',
            'version':            '0.1',
            'author_name':        'J. Scott Stewart',
            'grid_type':          'none',
            'time_step_type':     'fixed',
            'step_method':        'explicit',
            #-------------------------------------------------------------
            'comp_name':          'frostnumberGeo',
            'model_family':       'PermaModel',
            'time_units':         'days' }

        self._input_var_names = (
            'atmosphere_bottom_air__temperature',
            )

        self._output_var_names = (
            'frostnumber__air',            # Air Frost number
            'frostnumber__surface',        # Surface Frost number
            'frostnumber__stefan' )        # Stefan Frost number

        self._var_name_map = {
            # These are the corresponding CSDMS standard names
            # NOTE: we need to look up for the corresponding standard names
            'atmosphere_bottom_air__temperature': '_temperature_current',
            'datetime__start':                    '_start_date',
            'datetime__end':                      '_end_date',
            'frostnumber__air':                   'air_frost_number_Geo',
            'frostnumber__surface':               'surface_frost_number_Geo',
            'frostnumber__stefan':                'stefan_frost_number_Geo'}


        self._var_units_map = {
            # These are the links to the model's variables' units
            'atmosphere_bottom_air__temperature':                 'deg_C',
            'frostnumber__air':                                   '1',
            'frostnumber__surface':                               '1',
            'frostnumber__stefan':                                '1' }

    def initialize(self, cfg_file=None):
        self._model = frost_number_Geo.FrostnumberGeoMethod()

        self._model.status = 'initializing'
        self._model.initialize_frostnumberGeo_component()

        # Set the name of this component
        self._name = "Permamodel FrostnumberGeo Component"

        # Set the internal (frost number) variables that correspond
        # to the input and output variable names
        # Note: since we used Topoflow's _var_name_map for this, it is that
        self._values = _values = {
        # These are the links to the model's variables and
        # should be consistent with _var_name_map
            'atmosphere_bottom_air__temperature':
                                        self._model._temperature_current,
            'datetime__start':          self._model._start_date,
            'datetime__end':            self._model._end_date,
            'frostnumber__air':         self._model.air_frost_number_Geo,
            'frostnumber__surface':     self._model.surface_frost_number_Geo,
            'frostnumber__stefan':      self._model.stefan_frost_number_Geo}

        # Surface and Stefan numbers are not necessarily calculated
        if not self._model._calc_surface_fn and \
                'frostnumber__surface' in self._values.keys():
            del self._var_name_map['frostnumber__surface']
            del self._var_units_map['frostnumber__surface']
            del self._values['frostnumber__surface']
            # _output_var_names is an (immutable)
            self._output_var_names = tuple(v for v in self._output_var_names
                                            if v != 'frostnumber__surface')

        if not self._model._calc_stefan_fn and \
                'frostnumber__stefan' in self._values.keys():
            del self._var_name_map['frostnumber__stefan']
            del self._var_units_map['frostnumber__stefan']
            del self._values['frostnumber__stefan']
            self._output_var_names = tuple(v for v in self._output_var_names
                                            if v != 'frostnumber__stefan')

        # Verify that all input and output variable names are in the
        # variable name and the units map
        for varname in self._input_var_names:
            assert_in(varname, self._var_name_map)
            assert_in(varname, self._var_units_map)
        for varname in self._output_var_names:
            assert_in(varname, self._var_name_map)
            assert_in(varname, self._var_units_map)

        # Set the names and types of the grids
        # Note: A single value is a uniform rectilinear grid of shape (1)
        #       and size 1
        gridnumber = 0
        for varname in self._input_var_names:
            self._grids[gridnumber] = varname
            self._grid_type[gridnumber] = 'uniform_rectilinear'
            gridnumber += 1
        for varname in self._output_var_names:
            self._grids[gridnumber] = varname
            self._grid_type[gridnumber] = 'uniform_rectilinear'
            gridnumber += 1

        # initialize() tasks complete.  Update status.
        self._model.status = 'initialized'

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
        assert self._model.status == 'initialized'

        self._model.update()

    def update_frac(self, time_fraction):
        # Only increment the time by a partial time step
        # Ensure that we've already initialized the run
        assert(self._model.status == 'initialized')

        self._model.update_frac(frac=time_fraction)

    def update_until(self, stop_date):
        # Ensure that stop_year is at least the current year
        if stop_date < self._model._date_current:
            print("Warning: update_until year is less than current year")
            print("  no update run")
            return

        if stop_date > self._model._end_date:
            print("Warning: update_until year was greater than end_year")
            print("  setting stop_date to end_date")
            stop_date = self._end_date

        # Implement the loop to update until stop_year
        while self._model._date_current < stop_date:
            self.update()

    def finalize(self):
        SILENT = True

        # Finish with the run
        self._model.status = 'finalizing'  # (OpenMI)

        # frost_number_Geo has a finalize() method
        self._model.finalize()   # Close any input files

        # Done finalizing
        self._model.status = 'finalized'  # (OpenMI)

        # Print final report, as desired
        if not SILENT:
            self._model.print_final_report(\
                    comp_name='Permamodel FrostNumberGeo component')

    def get_start_time(self):
        return 0.0

    def get_current_time(self):
        """ Am assuming current time is number of timesteps since start """
        return((self._model._date_current - \
               self._model._start_date).total_seconds()/ \
               self._model._timestep_duration.total_seconds())

    def get_end_time(self):
        return((self._model._end_date - \
               self._model._start_date).total_seconds()/ \
               self._model._timestep_duration.total_seconds())

    # ----------------------------------
    # Functions added to pass bmi-tester
    # ----------------------------------
    def get_grid_type(self, grid_number):
        return self._grid_type[grid_number]

    def get_time_step(self):
        return self._model._timestep_duration.total_seconds() / \
               datetime.timedelta(days=1).total_seconds()

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

    #def set_value_at_indices(self, var_name, new_var_values, indices):
    def set_value_at_indices(self, var_name, indices, new_var_values):
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
