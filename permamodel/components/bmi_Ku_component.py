# -*- coding: utf-8 -*-
"""  Kudryavtsev model by Anisimov et al. (1997). 

Anisimov, O. A., Shiklomanov, N. I., & Nelson, F. E. (1997).
Global warming and active-layer thickness: results from transient general circulation models. 
Global and Planetary Change, 15(3-4), 61-77. DOI:10.1016/S0921-8181(97)00009-X

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

import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

# from permamodel.components.perma_base import *
# from permamodel.tests import examples_directory
import os

import numpy as np

from permamodel.components import Ku_method, perma_base
from permamodel.utils import model_input

"""
class BmiKuMethod( perma_base.PermafrostComponent ):
    _thisname = 'this name'
"""


class BmiKuMethod(perma_base.PermafrostComponent):

    """Implement the Ku model"""

    # Set up the name of this permafrost module
    _name = "Ku module"

    """ Note: all of these are defined below instead!
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
    """

    # -------------------------------------------------------------------
    _att_map = {
        # NOTE: this will change in the future
        "model_name": "PermaModel_Kudryavtsev_method",
        "version": "0.1",
        "author_name": "Kang Wang and Elchin Jafarov",
        "grid_type": "none",
        "time_step_type": "fixed",
        "step_method": "explicit",
        # -------------------------------------------------------------
        "comp_name": "Ku_model",
        "model_family": "PermaModel",
        "cfg_extension": "_ku_model.cfg",
        "cmt_var_prefix": "/input/",
        "gui_yaml_file": "/input/ku_model.yaml",
        "time_units": "years",
    }

    # This used to be [...] instead of (...)
    _input_var_names = (
        "latitude",
        "longitude",
        "datetime__start",
        "datetime__end",
        "atmosphere_bottom_air__temperature",
        "atmosphere_bottom_air__temperature_amplitude",
        "snowpack__depth",
        "snowpack__density",
        "water-liquid__volumetric-water-content-soil",
        "vegetation__Hvgf",
        "vegetation__Hvgt",
        "vegetation__Dvf",
        "vegetation__Dvt",
    )

    _output_var_names = (
        "soil__temperature",  # Tps
        "soil__active_layer_thickness",
    )  # Zal

    _var_name_map = {
        "latitude": "lat",
        "longitude": "lon",
        "datetime__start": "start_year",
        "datetime__end": "end_year",
        "atmosphere_bottom_air__temperature": "T_air",
        "atmosphere_bottom_air__temperature_amplitude": "A_air",
        "snowpack__depth": "h_snow",
        "snowpack__density": "rho_snow",
        "water-liquid__volumetric-water-content-soil": "vwc_H2O",
        "vegetation__Hvgf": "Hvgf",
        "vegetation__Hvgt": "Hvgt",
        "vegetation__Dvf": "Dvf",
        "vegetation__Dvt": "Dvt",
        "soil__temperature": "Tps",
        "soil__active_layer_thickness": "Zal",
    }

    _var_units_map = {
        # These are the links to the model's variables' units
        "latitude": "degree_north",
        "longitude": "degree_east",
        "datetime__start": "year",
        "datetime__end": "year",
        "atmosphere_bottom_air__temperature": "deg_C",
        "atmosphere_bottom_air__temperature_amplitude": "deg_C",
        "snowpack__depth": "m",
        "snowpack__density": "kg m-3",
        "water-liquid__volumetric-water-content-soil": "m3 m-3",
        "vegetation__Hvgf": "m",
        "vegetation__Hvgt": "m",
        "vegetation__Dvf": "m2 s-1",
        "vegetation__Dvt": "m2 s-1",
        "soil__temperature": "deg_C",
        "soil__active_layer_thickness": "m",
    }

    # -------------------------------------------------------------------
    def __init__(self):
        self._model = None
        self._values = {}
        self._var_units = {}
        self._grids = {}
        self._grid_type = {}

    def initialize(self, cfg_file=None):
        self._model = Ku_method.Ku_method()

        self._name = "Permamodel Ku Component"
        self._model.initialize(cfg_file=cfg_file)

        # make 2 vars to store each results and used for write out.
        n_lat = np.size(self._model.lat)
        n_lon = np.size(self._model.lon)
        n_time = self._model.end_year - self._model.start_year + 1

        self.output_alt = np.zeros((n_time, n_lat, n_lon)) * np.nan
        self.output_tps = np.zeros((n_time, n_lat, n_lon)) * np.nan

        # Verify that all input and output variable names are in the
        # variable name and the units map
        for varname in self._input_var_names:
            assert varname in self._var_name_map
            assert varname in self._var_units_map
            # print("Input var %s is in the name map and the units map"\
            #      % varname)
        for varname in self._output_var_names:
            assert varname in self._var_name_map
            assert varname in self._var_units_map

        self._model.cont = -1

        gridnumber = 0
        for varname in self._input_var_names:
            self._grids[gridnumber] = varname
            # self._grid_type[gridnumber] = 'uniform_rectilinear'
            self._grid_type[gridnumber] = "scalar"
            gridnumber += 1
        for varname in self._output_var_names:
            self._grids[gridnumber] = varname
            # self._grid_type[gridnumber] = 'uniform_rectilinear'
            self._grid_type[gridnumber] = "scalar"
            gridnumber += 1

        self._values = {
            # These are the links to the model's variables and
            # should be consistent with _var_name_map
            "latitude": self._model.lat,
            "longitude": self._model.lon,
            "datetime__start": self._model.start_year,
            "datetime__end": self._model.end_year,
            # 'atmosphere_bottom_air__temperature': "T_air",
            "atmosphere_bottom_air__temperature": self._model.T_air,
            "atmosphere_bottom_air__temperature_amplitude": self._model.A_air,
            "snowpack__depth": self._model.h_snow,
            "snowpack__density": self._model.rho_snow,
            "water-liquid__volumetric-water-content-soil": self._model.vwc_H2O,
            "vegetation__Hvgf": self._model.Hvgf,
            "vegetation__Hvgt": self._model.Hvgt,
            "vegetation__Dvf": self._model.Dvf,
            "vegetation__Dvt": self._model.Dvt,
            "soil__temperature": self._model.Tps,
            "soil__active_layer_thickness": self._model.Zal,
        }

        # Set the cfg file if it exists, otherwise, a default

    #        if cfg_file==None:
    #
    #        print self.cfg_file

    def get_attribute(self, att_name):
        try:
            return self._att_map[att_name.lower()]
        except:
            print("###################################################")
            print(" ERROR: Could not find attribute: " + att_name)
            print("###################################################")
            print(" ")

    #   get_attribute()
    # -------------------------------------------------------------------
    def get_input_var_names(self):
        # --------------------------------------------------------
        # Note: These are currently variables needed from other
        #       components vs. those read from files or GUI.
        # --------------------------------------------------------
        return self._input_var_names

    #   get_input_var_names()
    # -------------------------------------------------------------------
    def get_output_var_names(self):
        return self._output_var_names

    #   get_output_var_names()
    # -------------------------------------------------------------------
    def get_var_name(self, long_var_name):
        return self._var_name_map[long_var_name]

    #   get_var_name()
    # -------------------------------------------------------------------

    def get_var_units(self, long_var_name):
        return self._var_units_map[long_var_name]

    def get_var_location(self, long_var_name):
        return "node"

    def update(self):
        #        self._model.update(self._model.dt)
        # Ensure that we've already initialized the run
        assert self._model.status == "initialized"

        # Calculate the new frost number values
        self._model.update_ground_temperatures()
        self._model.update_ALT()

        self._values["soil__active_layer_thickness"] = self._model.Zal
        self._values["soil__temperature"] = self._model.Tps

        # Update the time
        self._model.year += self._model.dt

        self._model.cont = self._model.cont + 1

        #        self.output_alt = np.append(self.output_alt, self._model.Zal)
        #        self.output_tps = np.append(self.output_tps, self._model.Tps)
        self.output_alt[self._model.cont, :, :] = self._model.Zal
        self.output_tps[self._model.cont, :, :] = self._model.Tps

        # Get new input values
        self._model.read_input_files()

    def update_frac(self, time_fraction):
        time_step = self.get_time_step()
        self._model.dt = time_fraction * time_step
        self.update()
        self._model.dt = time_step

    def update_until(self, then):
        n_steps = (then - self.get_current_time()) / self.get_time_step()
        for _ in range(int(n_steps)):
            self.update()
        self.update_frac(n_steps - int(n_steps))

    def finalize(self):
        SILENT = True

        # Finish with the run
        self._model.status = "finalizing"  # (OpenMI)

        # Close the input files
        self._model.close_input_files()  # Close any input files

        # Write output last output
        self.save_grids()

        # Done finalizing
        self._model.status = "finalized"  # (OpenMI)

    def get_start_time(self):
        return 0.0

    def get_current_time(self):
        return float(self._model.year - self._model.start_year)

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
        setattr(self._model, self._var_name_map[var_name], new_var_values)
        # self._values[var_name] = new_var_values

    def set_value_at_indices(self, var_name, new_var_values, indices):
        self.get_value_ref(var_name).flat[indices] = new_var_values

    def get_var_itemsize(self, var_name):
        return np.asarray(self.get_value_ref(var_name)).flatten()[0].nbytes

    def get_value_at_indices(self, var_name, indices):
        return self.get_value_ref(var_name).take(indices)

    def get_var_nbytes(self, var_name):
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
        return 0

    def save_grids(self):
        # Saves the grid values based on the prescribed ones in cfg file

        # if (self.SAVE_MR_GRIDS):
        #    model_output.add_grid( self, self.T_air, 'T_air', self.time_min )
        #        self.ALT_file  = self.out_directory + self.ALT_file
        try:
            assert self._model.SAVE_ALT_GRIDS
        except:
            print("NO OUTPUT of ALT")
        try:
            assert self._model.SAVE_TPS_GRIDS
        except:
            print("NO OUTPUT of TPS")

        if self._model.SAVE_ALT_GRIDS:
            self._model.write_out_ncfile(self._model.ALT_file, self.output_alt)

        #        self.TPS_file  = self.out_directory + self.TPS_file

        if self._model.SAVE_TPS_GRIDS:
            self._model.write_out_ncfile(self._model.TPS_file, self.output_tps)

        print("***")
        print("Writing output finished!")
        print(
            "Please look at"
            + self._model.ALT_file
            + ".nc and "
            + self._model.TPS_file
            + ".nc"
        )
