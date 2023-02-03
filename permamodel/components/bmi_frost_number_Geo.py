# -*- coding: utf-8 -*-
"""  Frost Number by Nelson and Outcalt 1983.
     DOI: 10.2307/1551363. http://www.jstor.org/stable/1551363

     This is the Geo version
"""
from __future__ import print_function

import numpy as np

from permamodel.components import frost_number_Geo, perma_base


class BmiFrostnumberGeoMethod(perma_base.PermafrostComponent):
    """Implement the Nelson-Outcalt Frost numbers
    for a geographic region"""

    # Note: The definitions of _att_map, _var_name_map, etc.
    # have been moved inside the __init__(self) function because
    # otherwise changing thme once changes them for future
    # instantiations, .e.g during testing!

    # -------------------------------------------------------------------
    def __init__(self):
        self._model = None
        self._values = {}
        self._var_units = {}
        self._grids = {}
        self._grid_type = {}
        self.ngrids = 0
        # Set up the name of this permafrost module
        self._name = "Frost number module, Geo version"

        # -------------------------------------------------------------------
        self._att_map = {
            "model_name": "PermaModel_frostnumber_Geo_method",
            "version": "0.1",
            "author_name": "J. Scott Stewart",
            "grid_type": "rectilinear",
            "time_step_type": "fixed",
            "step_method": "explicit",
            # -------------------------------------------------------------
            "comp_name": "frostnumberGeo",
            "model_family": "PermaModel",
            "time_units": "years",
        }

        self._input_var_names = (
            "atmosphere_bottom_air__temperature",
            "atmosphere_bottom_air__temperature_mean_jan",
            "atmosphere_bottom_air__temperature_mean_jul",
        )

        self._output_var_names = (
            "frostnumber__air",  # Air Frost number
            "frostnumber__surface",  # Surface Frost number
            "frostnumber__stefan",
        )  # Stefan Frost number

        self._var_name_map = {
            # These are the corresponding CSDMS standard names
            # NOTE: we need to look up for the corresponding standard names
            "atmosphere_bottom_air__temperature": "_temperature_current",
            "atmosphere_bottom_air__temperature_mean_jan": "_temperature_jan",
            "atmosphere_bottom_air__temperature_mean_jul": "_temperature_jul",
            "datetime__start": "_start_date",
            "datetime__end": "_end_date",
            "frostnumber__air": "air_frost_number_Geo",
            "frostnumber__surface": "surface_frost_number_Geo",
            "frostnumber__stefan": "stefan_frost_number_Geo",
        }

        self._var_units_map = {
            # These are the links to the model's variables' units
            "atmosphere_bottom_air__temperature": "deg_C",
            "atmosphere_bottom_air__temperature_mean_jan": "deg_C",
            "atmosphere_bottom_air__temperature_mean_jul": "deg_C",
            "frostnumber__air": "1",
            "frostnumber__surface": "1",
            "frostnumber__stefan": "1",
        }

    def initialize(self, cfg_file=None):
        self._model = frost_number_Geo.FrostnumberGeoMethod(cfg_file)

        self._model.initialize_frostnumberGeo_component()

        # Set the name of this component
        self._name = "Permamodel FrostnumberGeo Component"

        # Set the internal (frost number) variables that correspond
        # to the input and output variable names
        # Note: since we used Topoflow's _var_name_map for this, it is that
        self._values = {
            # These are the links to the model's variables and
            # should be consistent with _var_name_map
            "atmosphere_bottom_air__temperature": self._model._temperature_current,
            "atmosphere_bottom_air__temperature_mean_jan": self._model._temperature_jan,
            "atmosphere_bottom_air__temperature_mean_jul": self._model._temperature_jul,
            "datetime__start": self._model._start_date,
            "datetime__end": self._model._end_date,
            "frostnumber__air": self._model.air_frost_number_Geo,
            "frostnumber__surface": self._model.surface_frost_number_Geo,
            "frostnumber__stefan": self._model.stefan_frost_number_Geo,
        }

        # Surface and Stefan numbers are not necessarily calculated
        if (
            not self._model._calc_surface_fn
            and "frostnumber__surface" in self._values.keys()
        ):
            del self._var_name_map["frostnumber__surface"]
            del self._var_units_map["frostnumber__surface"]
            del self._values["frostnumber__surface"]
            # _output_var_names is an (immutable)
            self._output_var_names = tuple(
                v for v in self._output_var_names if v != "frostnumber__surface"
            )

        if (
            not self._model._calc_stefan_fn
            and "frostnumber__stefan" in self._values.keys()
        ):
            del self._var_name_map["frostnumber__stefan"]
            del self._var_units_map["frostnumber__stefan"]
            del self._values["frostnumber__stefan"]
            self._output_var_names = tuple(
                v for v in self._output_var_names if v != "frostnumber__stefan"
            )

        # Verify that all input and output variable names are in the
        # variable name and the units map
        for varname in self._input_var_names:
            assert varname in self._var_name_map
            assert varname in self._var_units_map
        for varname in self._output_var_names:
            assert varname in self._var_name_map
            assert varname in self._var_units_map

        # Set the names and types of the grids
        # Note: A single value is a uniform rectilinear grid of shape (1)
        #       and size 1
        gridnumber = 0
        for varname in self._input_var_names:
            self._grids[gridnumber] = varname
            self._grid_type[gridnumber] = "uniform_rectilinear"
            gridnumber += 1
        for varname in self._output_var_names:
            self._grids[gridnumber] = varname
            self._grid_type[gridnumber] = "uniform_rectilinear"
            gridnumber += 1
        self.ngrids = gridnumber

    def get_attribute(self, att_name):
        try:
            return self._att_map[att_name.lower()]
        except KeyError:
            raise KeyError("No attribute %s" % str(att_name))

    def get_input_var_names(self):
        return self._input_var_names

    def get_output_var_names(self):
        return self._output_var_names

    def get_var_name(self, long_var_name):
        return self._var_name_map[long_var_name]

    def get_var_units(self, long_var_name):
        return self._var_units_map[long_var_name]

    def get_var_location(self, long_var_name):
        return "node"

    def update(self):
        self._model.update()
        self._values["frostnumber__air"] = self._model.air_frost_number_Geo

    def update_frac(self, time_fraction):
        # Only increment the time by a partial time step
        # Currently, only integer times are permitted
        self._model.update(frac=time_fraction)
        self._values["frostnumber__air"] = self._model.air_frost_number_Geo

    def update_until(self, time):
        """Advance model state until the given time.

        Parameters
        ----------
        time : float
          A model time value.

        """
        stop_year = time + self._model._start_date.year

        # Ensure that stop_year is at least the current year
        if stop_year < self._model._date_current.year:
            print("Warning: update_until year is less than current year")
            print("  no update run")
            return

        if stop_year > self._model._end_date.year:
            print("Warning: update_until year was greater than end_year")
            print("  setting stop_year to end_date")
            stop_year = self._model._end_date.year

        # Implement the loop to update until stop_year
        while self._model._date_current.year < stop_year:
            self.update()

    def finalize(self):
        # frost_number_Geo has a finalize() method
        self._model.finalize()  # Close any input files

    def get_start_time(self):
        return 0.0

    def get_current_time(self):
        """Am assuming current time is number of timesteps since start"""
        return self._model._date_current.year - self._model._start_date.year

    def get_end_time(self):
        return self._model._end_date.year - self._model._start_date.year

    # ----------------------------------
    # Functions added to pass bmi-tester
    # ----------------------------------
    def get_grid_type(self, grid_number):
        return self._grid_type[grid_number]

    def get_time_step(self):
        return float(self._model._timestep_duration)

    def get_value_ref(self, var_name):
        return self._values[var_name]

    def set_value(self, var_name, new_var_values):
        val = self.get_value_ref(var_name)
        val.flat = new_var_values

    def set_value_at_indices(self, var_name, indices, new_var_values):
        self.get_value_ref(var_name).flat[indices] = new_var_values

    def get_var_itemsize(self, var_name):
        return np.asarray(self.get_value_ref(var_name)).flatten()[0].nbytes

    def get_value_at_indices(self, var_name, indices):
        return self.get_value_ref(var_name).take(indices)

    def get_var_nbytes(self, var_name):
        return np.asarray(self.get_value_ref(var_name)).nbytes

    def get_value(self, var_name, out=None):
        value = self.get_value_ref(var_name).reshape((-1,))
        if out is None:
            out = value.copy()
        else:
            out[...] = value
        return out

    def get_var_type(self, var_name):
        return str(self.get_value_ref(var_name).dtype)

    def get_component_name(self):
        return self._name

    def get_var_grid(self, var_name):
        for grid_id, var_name_list in self._grids.items():
            if var_name in var_name_list:
                return grid_id

    def get_grid_shape(self, grid_id, out=None):
        var_name = self._grids[grid_id]
        grid_shape = np.array(self.get_value_ref(var_name)).shape
        if out is None:
            out = grid_shape
        else:
            out[:] = grid_shape
        return out

    def get_grid_size(self, grid_id):
        grid_size = self.get_grid_shape(grid_id)
        if grid_size == ():
            return 1
        else:
            return int(np.prod(grid_size))

    def get_grid_spacing(self, grid_id, out=None):
        assert grid_id < self.ngrids
        var_name = self._grids[grid_id]
        if out is None:
            ndim = np.array(self.get_value_ref(var_name)).ndim
            out = np.full(ndim, 1.0)
        else:
            out[:] = 1.0
        return out

    # Todo: Revise once we can work with georeferenced data in the CMF.
    def get_grid_origin(self, grid_id, out=None):
        var_name = self._grids[grid_id]
        if out is None:
            ndim = np.array(self.get_value_ref(var_name)).ndim
            out = np.full(ndim, 0.0)
        else:
            out[:] = 0.0
        return out

    def get_grid_rank(self, var_id):
        return len(self.get_grid_shape(var_id))
