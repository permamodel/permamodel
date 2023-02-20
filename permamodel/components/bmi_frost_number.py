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

import numpy as np

from permamodel.components import frost_number, perma_base


class BmiFrostnumberMethod(perma_base.PermafrostComponent):
    """Implement BMI interface to the Nelson-Outcalt Frost number code"""

    def __init__(self):
        """This overwrites __init__() method of PermafrostComponent"""
        self._name = "Permamodel Frostnumber Component"
        self._model = None
        self.model = None
        self._values = {}
        self._var_units = {}
        self._grids = {}
        self._grid_type = {}

        self._att_map = {
            "model_name": "PermaModel_frostnumber_method",
            "version": "0.1",
            "author_name": "J. Scott Stewart and Elchin Jafarov",
            "grid_type": "none",
            "time_step_type": "fixed",
            "step_method": "explicit",
            # -------------------------------------------------------------
            "comp_name": "frostnumber",
            "model_family": "PermaModel",
            "cfg_extension": "_frostnumber_model.cfg",
            "cmt_var_prefix": "/input/",
            "gui_yaml_file": "/input/frostnumber_model.yaml",
            "time_units": "years",
        }

        self._input_var_names = (
            "atmosphere_bottom_air__time_min_of_temperature",
            "atmosphere_bottom_air__time_max_of_temperature",
        )

        self._output_var_names = (
            "frostnumber__air",  # Air Frost number
            "frostnumber__surface",  # Surface Frost number
            "frostnumber__stefan",
        )  # Stefan Frost number

        self._var_name_map = {
            "atmosphere_bottom_air__time_min_of_temperature": "T_air_min",
            "atmosphere_bottom_air__time_max_of_temperature": "T_air_max",
            "frostnumber__air": "air_frost_number",
            "frostnumber__surface": "surface_frost_number",
            "frostnumber__stefan": "stefan_frost_number",
        }

        self._var_units_map = {
            "atmosphere_bottom_air__time_min_of_temperature": "deg_C",
            "atmosphere_bottom_air__time_max_of_temperature": "deg_C",
            "frostnumber__air": "1",
            "frostnumber__surface": "1",
            "frostnumber__stefan": "1",
        }

    def initialize(self, cfg_file=None):
        """this overwrites initialize() in PermafrostComponent"""
        self._model = frost_number.FrostnumberMethod()
        # This allows testing to not use protected access to _model
        self.model = self._model

        self._model.initialize_from_config_file(cfg_file=cfg_file)
        self._model.initialize_frostnumber_component()

        # Verify that all input and output variable names are in the
        # variable name and the units map
        for varname in self._input_var_names:
            assert varname in self._var_name_map
            assert varname in self._var_units_map
        for varname in self._output_var_names:
            assert varname in self._var_name_map
            assert varname in self._var_units_map

        # Set the Frost Number grids, based on input and output variables
        # Set the names and types of the grids
        gridnumber = 0
        for varname in self._input_var_names:
            self._grids[gridnumber] = varname
            self._grid_type[gridnumber] = "scalar"
            gridnumber += 1
        for varname in self._output_var_names:
            self._grids[gridnumber] = varname
            self._grid_type[gridnumber] = "scalar"
            gridnumber += 1

        # Set the internal (frost number) variables that correspond
        # to the input and output variable names
        # Note: since we used Topoflow's _var_name_map for this, it is that
        self._values = {
            # These are the links to the model's variables and
            # should be consistent with _var_name_map
            "atmosphere_bottom_air__time_min_of_temperature": self._model.T_air_min,
            "atmosphere_bottom_air__time_max_of_temperature": self._model.T_air_max,
            "frostnumber__air": self._model.air_frost_number,
            "frostnumber__surface": self._model.surface_frost_number,
            "frostnumber__stefan": self._model.stefan_frost_number,
        }

    def get_attribute(self, att_name):
        try:
            return self._att_map[att_name.lower()]
        except KeyError:
            raise KeyError("No attribute %s" % str(att_name))

    def get_input_var_names(self):
        # --------------------------------------------------------
        # Note: These are currently variables needed from other
        #       components vs. those read from files or GUI.
        # --------------------------------------------------------
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
        """This overwrites update() in Permafrost component
        and has different number of arguments"""
        # Calculate the new frost number values
        self._model.update()
        self._values["frostnumber__air"] = self._model.air_frost_number

    def update_frac(self, time_fraction):
        """This is for BMI compliance"""
        # Only increment the time by a partial time step
        # Currently, only non-fractions are allowed, but this could be
        #  0, 1, 2, ...
        self._model.update(frac=time_fraction)
        self._values["frostnumber__air"] = self._model.air_frost_number

    def get_year_from_timestep(self, this_timestep):
        """given the timestep, return the year"""
        return int(self.get_time_step() * this_timestep + self._model.start_year)

    def update_until(self, stop_time):
        """BMI-required, run until a specified time"""
        stop_year = self.get_year_from_timestep(stop_time)

        # Ensure that stop_year is at least the current year
        if stop_year < self._model.year:
            print("Warning: update_until year is less than current year")
            print("  no update run")
            return

        if stop_year > self._model.end_year:
            print("Warning: update_until year was greater than end_year")
            print("  setting stop_year to end_year")
            stop_year = self._model.end_year

        # Implement the loop to update until stop_year
        # Note: running for the first timestep is a special case,
        #       for which the timestep should not be updated
        if stop_year == self.get_year_from_timestep(self.get_start_time()):
            # print("This year is end_year, so update_frac(0)")
            self.update_frac(0)
        else:
            year = self._model.year
            while year < stop_year:
                self.update()
                year = self._model.year

    def finalize(self):
        """BMI-required, wrap up all things including writing output"""
        # Write output last output
        self._model.write_output_to_file()

    def get_start_time(self):
        """BMI-required although all WMT models start at time=0"""
        return 0.0

    def get_current_time(self):
        """Number of timesteps (years) since start of model run"""
        return float(self._model.year - self._model.start_year)

    def get_end_time(self):
        """Number of timestesp (years) so that last year is included"""
        return float(self._model.end_year - self._model.start_year)

    def get_grid_type(self, grid_number):
        """BMI: inquire about the type of this grid"""
        return self._grid_type[grid_number]

    def get_time_step(self):
        """BMI: return the time step value (years)"""
        return float(self._model.dt)

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
        """BMI: allow external set-access to model variable"""
        array = getattr(self._model, self._var_name_map[var_name])
        if len(array) == 1:
            array[0] = new_var_values
        else:
            array[int(self._model.year - self._model.start_year) + 1] = new_var_values

        # setattr(self._model, self._var_name_map[var_name], new_var_values)
        # self._values[var_name] = new_var_values

    def set_value_at_indices(self, var_name, new_var_values, indices):
        """BMI: allow external set-access to model array variable"""
        self.get_value_ref(var_name).flat[indices] = new_var_values

    def get_var_itemsize(self, var_name):
        """BMI: determine how many bytes a variable uses"""
        return np.asarray(self.get_value_ref(var_name)).flatten()[0].nbytes

    def get_value_at_indices(self, var_name, indices):
        """BMI: allow external get-access to model array variable"""
        return self.get_value_ref(var_name).take(indices)

    def get_var_nbytes(self, var_name):
        """BMI: number of bytes of a variable"""
        return np.asarray(self.get_value_ref(var_name)).nbytes

    def get_value(self, var_name, out=None):
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
        # return self.get_value_ref(var_name).copy()

        # Version to convert to numpy array for bmi-tester compliance
        # Note: converting to np arrays on the fly here
        # Note: float values don't have a copy() function
        # try:
        #    return np.array(self.get_value_ref(var_name).copy())
        # except AttributeError:
        #    return np.array(self.get_value_ref(var_name))
        if out is None:
            out = self.get_value_ref(var_name).copy()
        else:
            out[...] = self.get_value_ref(var_name)
        return out

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
        """BMI: provide this component's name"""
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
        return 1

    def get_grid_rank(self, grid_id):
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
        return 0
