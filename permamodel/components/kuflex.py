# -*- coding: utf-8 -*-
import numpy as np
from bmipy import Bmi

from .KuFlex_method import KuFlex_method


class KuFlex(Bmi):
    _name = "KuFlex module"

    _input_var_names = (
        "datetime__start",
        "datetime__end",
        "atmosphere_bottom_air__temperature",
        "atmosphere_bottom_air__temperature_amplitude",
        "snowpack__depth",
        "snowpack__density",
        "snow__thermal_conductivity",
        "snow__volume-specific_isochoric_heat_capacity",
        "water__volume_latent_fusion_heat",
        "soil-frozen__volume-specific_isochoric_heat_capacity",
        "soil-thaw__volume-specific_isochoric_heat_capacity",
        "soil-frozen__thermal_conductivity",
        "soil-thaw__thermal_conductivity",
        "vegetation__Hvgf",
        "vegetation__Hvgt",
        "vegetation__Dvf",
        "vegetation__Dvt",
    )

    _output_var_names = (
        "soil_surface__temperature",
        "soil_surface__temperature_amplitude",
        "soil__temperature",  # Tps
        "soil__active_layer_thickness",  # Zal
    )

    _var_name = {
        "datetime__start": "start_year",
        "datetime__end": "end_year",
        "atmosphere_bottom_air__temperature": "T_air",
        "atmosphere_bottom_air__temperature_amplitude": "A_air",
        "snowpack__depth": "h_snow",
        "snowpack__density": "rho_snow",
        "snow__thermal_conductivity": "k_snow",
        "snow__volume-specific_isochoric_heat_capacity": "c_snow",
        "water__volume_latent_fusion_heat": "L",
        "soil-frozen__volume-specific_isochoric_heat_capacity": "Cf",
        "soil-thaw__volume-specific_isochoric_heat_capacity": "Ct",
        "soil-frozen__thermal_conductivity": "Kf",
        "soil-thaw__thermal_conductivity": "Kt",
        "vegetation__Hvgf": "Hvgf",
        "vegetation__Hvgt": "Hvgt",
        "vegetation__Dvf": "Dvf",
        "vegetation__Dvt": "Dvt",
        "soil_surface__temperature": "Tgs",
        "soil_surface__temperature_amplitude": "Ags",
        "soil__temperature": "Tps",
        "soil__active_layer_thickness": "Zal",
    }

    _var_units = {
        # These are the links to the model's variables' units
        "datetime__start": "year",
        "datetime__end": "year",
        "atmosphere_bottom_air__temperature": "deg_C",
        "atmosphere_bottom_air__temperature_amplitude": "deg_C",
        "snowpack__depth": "m",
        "snowpack__density": "kg m-3",
        "snow__thermal_conductivity": "W m-1 K-1",
        "snow__volume-specific_isochoric_heat_capacity": "J m-3 K-1",
        "water__volume_latent_fusion_heat": "J m-3",
        "soil-frozen__volume-specific_isochoric_heat_capacity": "J m-3 K-1",
        "soil-thaw__volume-specific_isochoric_heat_capacity": "J m-3 K-1",
        "soil-frozen__thermal_conductivity": "W m-1 K-1",
        "soil-thaw__thermal_conductivity": "W m-1 K-1",
        "vegetation__Hvgf": "m",
        "vegetation__Hvgt": "m",
        "vegetation__Dvf": "m2 s-1",
        "vegetation__Dvt": "m2 s-1",
        "soil_surface__temperature": "deg_C",
        "soil_surface__temperature_amplitude": "deg_C",
        "soil__temperature": "deg_C",
        "soil__active_layer_thickness": "m",
    }

    def finalize(self):
        """Perform tear-down tasks for the model.

        Perform all tasks that take place after exiting the model's time
        loop. This typically includes deallocating memory, closing files and
        printing reports.
        """
        self._model.status = "finalizing"
        self._model.close_input_files()
        self._model.close_output_files()
        self._model.status = "finalized"

    def get_component_name(self):
        """Name of the component.

        Returns
        -------
        str
            The name of the component.
        """
        return self._name

    def get_current_time(self):
        """Current time of the model.

        Returns
        -------
        float
            The current model time.
        """
        return float(self._model.year - self._model.start_year)

    def get_end_time(self):
        """End time of the model.

        Returns
        -------
        float
            The maximum model time.
        """
        return self._model.end_year - self._model.start_year + 1.0

    def get_grid_edge_count(self, grid):
        """Get the number of edges in the grid.

        Parameters
        ----------
        grid : int
            A grid identifier.

        Returns
        -------
        int
            The total number of grid edges.
        """
        raise NotImplementedError("get_grid_edge_count")

    def get_grid_edge_nodes(self, grid, edge_nodes):
        """Get the edge-node connectivity.

        Parameters
        ----------
        grid : int
            A grid identifier.
        edge_nodes : ndarray of int, shape *(2 x nnodes,)*
            A numpy array to place the edge-node connectivity. For each edge,
            connectivity is given as node at edge tail, followed by node at
            edge head.

        Returns
        -------
        ndarray of int
            The input numpy array that holds the edge-node connectivity.
        """
        raise NotImplementedError("get_grid_edge_nodes")

    def get_grid_face_count(self, grid):
        """Get the number of faces in the grid.

        Parameters
        ----------
        grid : int
            A grid identifier.

        Returns
        -------
        int
            The total number of grid faces.
        """
        raise NotImplementedError("get_grid_face_count")

    def get_grid_face_nodes(self, grid, face_nodes):
        """Get the face-node connectivity.

        Parameters
        ----------
        grid : int
            A grid identifier.
        face_nodes : ndarray of int
            A numpy array to place the face-node connectivity. For each face,
            the nodes (listed in a counter-clockwise direction) that form the
            boundary of the face.

        Returns
        -------
        ndarray of int
            The input numpy array that holds the face-node connectivity.
        """
        raise NotImplementedError("get_grid_face_nodes")

    def get_grid_node_count(self, grid):
        """Get the number of nodes in the grid.

        Parameters
        ----------
        grid : int
            A grid identifier.

        Returns
        -------
        int
            The total number of grid nodes.
        """
        raise NotImplementedError("get_grid_node_count")

    def get_grid_nodes_per_face(self, grid, nodes_per_face):
        """Get the number of nodes for each face.

        Parameters
        ----------
        grid : int
            A grid identifier.
        nodes_per_face : ndarray of int, shape *(nfaces,)*
            A numpy array to place the number of edges per face.

        Returns
        -------
        ndarray of int
            The input numpy array that holds the number of nodes per edge.
        """
        raise NotImplementedError("get_grid_nodes_per_face")

    def get_grid_origin(self, grid, origin):
        """Get coordinates for the lower-left corner of the computational grid.

        Parameters
        ----------
        grid : int
            A grid identifier.
        origin : ndarray of float, shape *(ndim,)*
            A numpy array to hold the coordinates of the lower-left corner of
            the grid.

        Returns
        -------
        ndarray of float
            The input numpy array that holds the coordinates of the grid's
            lower-left corner.
        """
        origin[:] = 0.0
        return origin

    def get_grid_rank(self, grid):
        """Get number of dimensions of the computational grid.

        Parameters
        ----------
        grid : int
            A grid identifier.

        Returns
        -------
        int
            Rank of the grid.
        """
        if grid == 0:
            return 0
        else:
            return 2

    def get_grid_shape(self, grid, shape):
        """Get dimensions of the computational grid.

        Parameters
        ----------
        grid : int
            A grid identifier.
        shape : ndarray of int, shape *(ndim,)*
            A numpy array into which to place the shape of the grid.

        Returns
        -------
        ndarray of int
            The input numpy array that holds the grid's shape.
        """
        shape[:] = self._model.grid_shape
        return shape

    def get_grid_size(self, grid):
        """Get the total number of elements in the computational grid.

        Parameters
        ----------
        grid : int
            A grid identifier.

        Returns
        -------
        int
            Size of the grid.
        """
        return np.prod(self._model.grid_shape)

    def get_grid_spacing(self, grid, spacing):
        """Get distance between nodes of the computational grid.

        Parameters
        ----------
        grid : int
            A grid identifier.
        spacing : ndarray of float, shape *(ndim,)*
            A numpy array to hold the spacing between grid rows and columns.

        Returns
        -------
        ndarray of float
            The input numpy array that holds the grid's spacing.
        """
        spacing[:] = 1.0
        return spacing

    def get_grid_type(self, grid):
        """Get the grid type as a string.

        Parameters
        ----------
        grid : int
            A grid identifier.

        Returns
        -------
        str
            Type of grid as a string.
        """
        if grid == 0:
            return "none"
        else:
            return "uniform_rectilinear"

    def get_grid_x(self, grid, x):
        """Get coordinates of grid nodes in the x direction.

        Parameters
        ----------
        grid : int
            A grid identifier.
        x : ndarray of float, shape *(nrows,)*
            A numpy array to hold the x-coordinates of the grid node columns.

        Returns
        -------
        ndarray of float
            The input numpy array that holds the grid's column x-coordinates.
        """
        raise NotImplementedError("get_grid_x")

    def get_grid_y(self, grid, y):
        """Get coordinates of grid nodes in the y direction.

        Parameters
        ----------
        grid : int
            A grid identifier.
        y : ndarray of float, shape *(ncols,)*
            A numpy array to hold the y-coordinates of the grid node rows.

        Returns
        -------
        ndarray of float
            The input numpy array that holds the grid's row y-coordinates.
        """
        raise NotImplementedError("get_grid_y")

    def get_grid_z(self, grid, z):
        """Get coordinates of grid nodes in the z direction.

        Parameters
        ----------
        grid : int
            A grid identifier.
        z : ndarray of float, shape *(nlayers,)*
            A numpy array to hold the z-coordinates of the grid nodes layers.

        Returns
        -------
        ndarray of float
            The input numpy array that holds the grid's layer z-coordinates.
        """
        raise NotImplementedError("get_grid_z")

    def get_input_var_names(self):
        """List of a model's input variables.

        Input variable names must be CSDMS Standard Names, also known
        as *long variable names*.

        Returns
        -------
        list of str
            The input variables for the model.

        Notes
        -----
        Standard Names enable the CSDMS framework to determine whether
        an input variable in one model is equivalent to, or compatible
        with, an output variable in another model. This allows the
        framework to automatically connect components.

        Standard Names do not have to be used within the model.
        """
        return self._input_var_names

    def get_output_var_names(self):
        """List of a model's output variables.

        Output variable names must be CSDMS Standard Names, also known
        as *long variable names*.

        Returns
        -------
        list of str
            The output variables for the model.
        """
        return self._output_var_names

    def get_start_time(self):
        """Start time of the model.

        Model times should be of type float.

        Returns
        -------
        float
            The model start time.
        """
        return 0.0

    def get_time_step(self):
        """Current time step of the model.

        The model time step should be of type float.

        Returns
        -------
        float
            The time step used in model.
        """
        return self._model.dt

    def get_time_units(self):
        """Time units of the model.

        Returns
        -------
        float
            The model time unit; e.g., `days` or `s`.

        Notes
        -----
        CSDMS uses the UDUNITS standard from Unidata.
        """
        return "years"

    def get_value(self, name, dest):
        """Get a copy of values of the given variable.

        This is a getter for the model, used to access the model's
        current state. It returns a *copy* of a model variable, with
        the return type, size and rank dependent on the variable.

        Parameters
        ----------
        name : str
            An input or output variable name, a CSDMS Standard Name.
        dest : ndarray
            A numpy array into which to place the values.

        Returns
        -------
        ndarray
            The same numpy array that was passed as an input buffer.
        """
        array = self.get_value_ptr(name)
        dest[:] = array
        return dest

    def get_value_at_indices(self, name, dest, inds):
        """Get values at particular indices.

        Parameters
        ----------
        name : str
            An input or output variable name, a CSDMS Standard Name.
        dest : ndarray
            A numpy array into which to place the values.
        indices : array_like
            The indices into the variable array.

        Returns
        -------
        array_like
            Value of the model variable at the given location.
        """
        array = self.get_value_ptr(name)
        dest[:] = array[inds]
        return dest

    def get_value_ptr(self, name):
        """Get a reference to values of the given variable.

        This is a getter for the model, used to access the model's
        current state. It returns a reference to a model variable,
        with the return type, size and rank dependent on the variable.

        Parameters
        ----------
        name : str
            An input or output variable name, a CSDMS Standard Name.

        Returns
        -------
        array_like
            A reference to a model variable.
        """
        return getattr(self._model, self._var_name[name]).reshape(-1)

    def get_var_grid(self, name):
        """Get grid identifier for the given variable.

        Parameters
        ----------
        name : str
            An input or output variable name, a CSDMS Standard Name.

        Returns
        -------
        int
          The grid identifier.
        """
        try:
            array = self.get_value_ptr(name)
        except AttributeError:
            raise RuntimeError("no grid for {name}".format(name=name))
        else:
            if array.size > 1:
                grid = 1
            else:
                grid = 0
        return grid

    def get_var_itemsize(self, name):
        """Get memory use for each array element in bytes.

        Parameters
        ----------
        name : str
            An input or output variable name, a CSDMS Standard Name.

        Returns
        -------
        int
            Item size in bytes.
        """
        return self.get_value_ref(name).itemsize

    def get_var_location(self, name):
        """Get the grid element type that the a given variable is defined on.

        The grid topology can be composed of *nodes*, *edges*, and *faces*.

        *node*
            A point that has a coordinate pair or triplet: the most
            basic element of the topology.

        *edge*
            A line or curve bounded by two *nodes*.

        *face*
            A plane or surface enclosed by a set of edges. In a 2D
            horizontal application one may consider the word "polygon",
            but in the hierarchy of elements the word "face" is most common.

        Parameters
        ----------
        name : str
            An input or output variable name, a CSDMS Standard Name.

        Returns
        -------
        str
            The grid location on which the variable is defined. Must be one of
            `"node"`, `"edge"`, or `"face"`.

        Notes
        -----
        CSDMS uses the `ugrid conventions`_ to define unstructured grids.

        .. _ugrid conventions: http://ugrid-conventions.github.io/ugrid-conventions
        """
        return "node"

    def get_var_nbytes(self, name):
        """Get size, in bytes, of the given variable.

        Parameters
        ----------
        name : str
            An input or output variable name, a CSDMS Standard Name.

        Returns
        -------
        int
            The size of the variable, counted in bytes.
        """
        return self.get_value_ref(name).nbytes

    def get_var_type(self, name):
        """Get data type of the given variable.

        Parameters
        ----------
        name : str
            An input or output variable name, a CSDMS Standard Name.

        Returns
        -------
        str
            The Python variable type; e.g., ``str``, ``int``, ``float``.
        """
        return str(self.get_value_ref(name).dtype)

    def get_var_units(self, name):
        """Get units of the given variable.

        Standard unit names, in lower case, should be used, such as
        ``meters`` or ``seconds``. Standard abbreviations, like ``m`` for
        meters, are also supported. For variables with compound units,
        each unit name is separated by a single space, with exponents
        other than 1 placed immediately after the name, as in ``m s-1``
        for velocity, ``W m-2`` for an energy flux, or ``km2`` for an
        area.

        Parameters
        ----------
        name : str
            An input or output variable name, a CSDMS Standard Name.

        Returns
        -------
        str
            The variable units.

        Notes
        -----
        CSDMS uses the `UDUNITS`_ standard from Unidata.

        .. _UDUNITS: http://www.unidata.ucar.edu/software/udunits
        """
        return self._var_units[name]

    def initialize(self, config_file):
        """Perform startup tasks for the model.

        Perform all tasks that take place before entering the model's time
        loop, including opening files and initializing the model state. Model
        inputs are read from a text-based configuration file, specified by
        `filename`.

        Parameters
        ----------
        config_file : str, optional
            The path to the model configuration file.

        Notes
        -----
        Models should be refactored, if necessary, to use a
        configuration file. CSDMS does not impose any constraint on
        how configuration files are formatted, although YAML is
        recommended. A template of a model's configuration file
        with placeholder values is used by the BMI.
        """
        self._model = KuFlex_method()

        self._name = "Permamodel Ku Component"
        self._model.initialize(cfg_file=config_file)

        # make 2 vars to store each results and used for write out.
        n_lat = self._model.grid_shape[0]
        n_lon = self._model.grid_shape[1]
        n_time = self._model.end_year - self._model.start_year + 1

        self.output_alt = np.full((n_time, n_lat, n_lon), np.nan)
        self.output_tps = np.full((n_time, n_lat, n_lon), np.nan)
        self.output_tgs = np.full((n_time, n_lat, n_lon), np.nan)
        self.output_ags = np.full((n_time, n_lat, n_lon), np.nan)

        # Verify that all input and output variable names are in the
        # variable name and the units map
        var_names = set(self._input_var_names + self._output_var_names)
        missing = var_names - set(self._var_name)
        if missing:
            raise ValueError(
                "_var_name is missing {missing}".format(", ".join(missing=missing))
            )

        missing = var_names - set(self._var_units)
        if missing:
            raise ValueError(
                "_var_units is missing {missing}".format(", ".join(missing=missing))
            )

    def set_value(self, name, values):
        """Specify a new value for a model variable.

        This is the setter for the model, used to change the model's
        current state. It accepts, through *src*, a new value for a
        model variable, with the type, size and rank of *src*
        dependent on the variable.

        Parameters
        ----------
        var_name : str
            An input or output variable name, a CSDMS Standard Name.
        src : array_like
            The new value for the specified variable.
        """
        self.get_value_ptr(name)[:] = values

    def set_value_at_indices(self, name, inds, src):
        """Specify a new value for a model variable at particular indices.

        Parameters
        ----------
        var_name : str
            An input or output variable name, a CSDMS Standard Name.
        indices : array_like
            The indices into the variable array.
        src : array_like
            The new value for the specified variable.
        """
        self.get_value_ptr(name)[inds] = src

    def update(self):
        """Advance model state by one time step.

        Perform all tasks that take place within one pass through the model's
        time loop. This typically includes incrementing all of the model's
        state variables. If the model's state variables don't change in time,
        then they can be computed by the :func:`initialize` method and this
        method can return with no action.
        """
        self._model.update()

    def update_frac(self, time_fraction):
        time_step = self.get_time_step()
        self._model.dt = time_fraction * time_step
        self.update()
        self._model.dt = time_step

    def update_until(self, then):
        if then > self.get_end_time():
            raise RuntimeError("{then} exceeds end time".format(then=then))
        n_steps = (then - self.get_current_time()) / self.get_time_step()
        for _ in range(int(n_steps)):
            self.update()
