"""Kudryavtsev permafrost model, adapted from Anisimov et al. (1997).

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

import numpy as np

from permamodel.components.Ku import Ku_model


class BmiKuModel:
    """Basic model interface for the Kudryavstev permafrost model."""

    _name = "Kudryavtsev Permafrost Model"

    _input_var_names = [
        "air_temperature",
        "temperature_amplitude",
        "snow_thickness",
        "snow_density",
        "soil_water_content",
        "frozen_vegetation_height",
        "thawed_vegetation_height",
        "frozen_vegetation_diffusivity",
        "thawed_vegetation_diffusivity",
    ]

    _output_var_names = ["permafrost_temperature", "active_layer_thickness"]

    _var_units_map = {
        "air_temperature": "degrees C",
        "temperature_amplitude": "degrees C",
        "snow_thickness": "m",
        "snow_density": "kg m-3",
        "soil_water_content": "m3 m-3",  # water / soil
        "frozen_vegetation_height": "m",
        "thawed_vegetation_height": "m",
        "frozen_vegetation_diffusivity": "m2 s-1",
        "thawed_vegetation_diffusivity": "m2 s-1",
        "permafrost_temperature": "degrees C",
        "active_layer_thickness": "m",
    }

    def __init__(self):
        """Initialize the Basic Model Interface."""
        self._model = None
        self._values = {}
        self._var_units = {}
        self._var_loc = {}
        self._grids = {}
        self._grid_type = {}

        self._start_time = 0.0
        self._end_time = None
        self._current_time = 0.0

        # Using 'years' as a time unit is generally not preferred
        # However, the Ku model does not support ANY time step other than 1 year
        self._time_units = "years"

    def initialize(self, config_file: str):
        """Initialize the Kudryavstev permafrost model.

        Args:
            config_file: str
                Path to the configuration file.
        """
        self._model = Ku_model()

        # Initialization routines
        self._model.initialize(config_file)

        self._values = {
            var: getattr(self._model, var).values for var in self._input_var_names
        }
        for var in self._output_var_names:
            self._values[var] = np.empty(
                (
                    self._model.number_of_years,
                    self._model.grid_shape[0],
                    self._model.grid_shape[1],
                )
            )

        self._var_units = self._var_units_map.copy()

        self._var_loc = {
            var: "node" for var in self._input_var_names + self._output_var_names
        }

        self._grids = {
            i: list(self._var_units_map.keys())[i]
            for i in range(len(self._var_units_map.keys()))
        }

        self._grid_type = {
            i: "uniform_rectilinear" for i in range(len(self._var_units_map.keys()))
        }

        self._start_time = 0.0
        self._end_time = self._model.number_of_years
        self._grid_shape = (
            self._model.number_of_years,
            self._model.grid_shape[0],
            self._model.grid_shape[1],
        )

    def update(self):
        """Run the model for the current time step."""
        self._model.run_one_step(self._current_time)
        self._current_time += 1.0

    def update_until(self, end_time: int):
        """Update the model until a certain year (inclusive)."""
        years = self._current_time + end_time

        for t in np.arange(self._current_time, end_time + 1, 1):
            self._model.run_one_step(t)

        self._current_time = end_time

    def finalize(self, path_to_output=None):
        """If specified, write output to a netcdf file."""

        if path_to_output is not None:
            self._model.write_output(
                path_to_output, vars_to_write=self._output_var_names
            )

    def get_component_name(self) -> str:
        """Return the name of the component."""
        return self._name

    def get_input_item_count(self) -> int:
        """Return the number of input variables."""
        return len(self._input_var_names)

    def get_output_item_count(self) -> int:
        """Return the number of output variables."""
        return len(self._output_var_names)

    def get_input_var_names(self) -> list:
        """Return a list of input variables."""
        return self._input_var_names

    def get_output_var_names(self) -> list:
        """Return a list of output variables."""
        return self._output_var_names

    def get_var_grid(self, var: str) -> int:
        """Return the grid ID for a variable."""
        for grid, variables in self._grids.items():
            if var in variables:
                return grid

    def get_var_type(self, var: str) -> str:
        """Return the data type of a variable."""
        return str(self.get_value_ptr(var).dtype)

    def get_var_units(self, var: str) -> str:
        """Return the units of a variable."""
        return self._var_units[var]

    def get_var_nbytes(self, var: str) -> int:
        """Return the size of a variable in bytes."""
        return self.get_value_ptr(var).nbytes

    def get_var_itemsize(self, var: str) -> int:
        """Return the size of one element of a variable in bytes."""
        return np.dtype(self.get_var_type(var)).itemsize

    def get_var_location(self, var: str) -> str:
        """Returns the location of a variable on the grid: either 'node', 'edge', or 'face'."""
        return "node"

    def get_current_time(self) -> int:
        """Return the current time."""
        return self._current_time

    def get_start_time(self) -> int:
        """Return the start time."""
        return self._start_time

    def get_end_time(self) -> int:
        """Return the end time."""
        return self._end_time

    def get_time_units(self) -> str:
        """Return the time units of the model."""
        return self._time_units

    def get_time_step(self) -> str:
        """Return the model's time step."""
        return 1.0

    def get_value_ptr(self, var: str) -> np.ndarray:
        """Returns a reference to the values of a variable."""
        return getattr(self._model, var).values[:]

    def get_value(self, var: str, dest: np.ndarray) -> np.ndarray:
        """Returns an array with a copy of the values of a variable."""
        dest[:] = np.ravel(self.get_value_ptr(var), order="C")
        return dest

    def get_value_at_indices(
        self, var: str, dest: np.ndarray, indices: np.ndarray
    ) -> np.ndarray:
        """Returns an array with a copy of the values of a variable at the given indices."""
        dest[:] = np.ravel(self.get_value_ptr(var), order="C").take(indices)
        return dest

    def set_value(self, var: str, values: np.ndarray):
        """Set the value of a variable."""
        vals = np.reshape(values, self._grid_shape)
        getattr(self._model, var).values[:] = vals

    def set_value_at_indices(self, var: str, indices: np.ndarray, values: np.ndarray):
        """Set the value of a variable at the given indices."""
        current_values = getattr(self._model, var).values[:]
        flat_values = np.ravel(current_values, order="C")
        flat_values[indices] = values

    def get_grid_type(self, grid: int) -> str:
        """Given a grid ID, return the type of grid."""
        return self._grid_type[grid]

    def get_grid_rank(self, grid: int) -> int:
        """Given a grid ID, return the number of dimensions of the grid."""
        return len(self._grid_shape)

    def get_grid_size(self, grid: int) -> int:
        """Given a grid ID, return the number of elements of the grid."""
        return int(np.prod(self._grid_shape))

    # TODO assign all following fns to arrays
    def get_grid_shape(self, grid: int, shape: np.ndarray) -> np.ndarray:
        """Given a grid ID, return the shape of the grid."""
        shape[:] = self._grid_shape
        return shape

    def get_grid_spacing(self, grid: int, spacing: np.ndarray) -> np.ndarray:
        """Given a grid ID, return the distance between grid elements in each direction."""
        var = self._grids[grid]

        dims = getattr(self._model, var).dims
        coords = [getattr(self._model, var)[dim] for dim in dims]
        diffs = [np.diff(array)[0] for array in coords]

        spacing[:] = diffs
        return spacing

    def get_grid_origin(self, grid: int, origin: np.ndarray) -> np.ndarray:
        """Given a grid ID, return the [y, x] coordinates of the lower-left corner."""
        var = self._grids[grid]

        ydim = getattr(self._model, var).dims[1]
        xdim = getattr(self._model, var).dims[2]

        y0 = getattr(self._model, var)[ydim][0, 0]
        x0 = getattr(self._model, var)[xdim][0, 0]

        origin[:] = [y0, x0]
        return origin

    def get_grid_x(self, grid: int, xs: np.ndarray) -> np.ndarray:
        """Given a grid ID, return an array of x-coordinates."""
        var = self._grids[grid]
        xs[:] = np.ravel(getattr(self._model, var).dims[2])
        return xs

    def get_grid_y(self, grid: int, ys: np.ndarray) -> np.ndarray:
        """Given a grid ID, return an array of y-coordinates."""
        var = self._grids[grid]
        ys[:] = np.ravel(getattr(self._model, var).dims[1])
        return ys

    def get_grid_z(self, grid: int, zs: np.ndarray):
        raise NotImplementedError("get_grid_z() not implemented.")

    def get_grid_node_count(self, grid: int):
        """Given a grid ID, return the number of nodes."""
        return self.get_grid_size(grid)

    def get_grid_face_count(self, grid: int):
        raise NotImplementedError("get_grid_face_count() not implemented.")

    def get_grid_edge_nodes(self, grid: int, edge_nodes: np.ndarray):
        raise NotImplementedError("get_grid_edge_nodes() not implemented.")

    def get_grid_face_edges(self, grid: int, face_edges: np.ndarray):
        raise NotImplementedError("get_grid_face_edges() not implemented.")

    def get_grid_face_nodes(self, grid: int, face_nodes: np.ndarray):
        raise NotImplementedError("get_grid_face_nodes() not implemented.")

    def get_grid_nodes_per_face(self, grid: int, nodes_per_face: np.ndarray):
        raise NotImplementedError("get_grid_nodes_per_face() not implemented.")
