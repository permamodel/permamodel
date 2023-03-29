"""Kudryavtsev permafrost model, adapted from Anisimov et al. (1997).

Input variables:
    air_temperature: Mean annual air temperature (degrees Celsius)
    temperature_amplitude: Amplitude of annual air temperature (degrees Celsius)
    snow_thickness: Mean thickness of winter snow cover (meters)
    snow_density: Mean density of winter snowpack (kilograms per cubic meter)
    soil_water_content: Mean volumetric water content (cubic meters of water per cubic meter of soil)
    frozen_vegetation_height: Mean height of vegetation during the frozen season(s) (meters)
    thawed_vegetation_height: Mean height of vegetation during the thawed season(s) (meters)
    frozen_vegetation_diffusivity: Mean thermal diffusivity of vegetation during the frozen season(s) (square meters per second)
    thawed_vegetation_diffusivity: Mean thermal diffusivity of vegetation during the thawed season(s) (square meters per second)

Output variables:
    permafrost_temperature: Mean annual temperature at the top of the permafrost layer (degrees Celsius)
    active_layer_thickness: Mean annual active layer thickness (meters)

Authors: Kang Wang, Elchin Jafarov, Ethan Pierce, Irina Overeem

References:
Anisimov, O. A., Shiklomanov, N. I., & Nelson, F. E. (1997).
    Global warming and active-layer thickness: results from transient general circulation models.
    Global and Planetary Change, 15(3), 61-77.
Romanovsky, V. E., & Osterkamp, T. E. (1997).
    Thawing of the active layer on the coastal plain of the Alaskan Arctic.
    Permafrost and Periglacial processes, 8(1), 1-22.
Sazonova, T. S., & Romanovsky, V. E. (2003).
    A model for regional‐scale estimation of temporal and spatial variability of active layer thickness and mean annual ground temperatures.
    Permafrost and Periglacial Processes, 14(2), 125-139.
Sturm, M., Holmgren, J., König, M., & Morris, K. (1997).
    The thermal conductivity of seasonal snow. Journal of Glaciology, 43(143), 26-41.
Ling, F., & Zhang, T. (2004).
    A numerical model for surface energy balance and thermal regime of the active layer and permafrost containing unfrozen water.
    Cold Regions Science and Technology, 38(1), 1-15.

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
import toml
import xarray as xr


class Ku_model:
    """The Kudryavtsev permafrost model.

    Typical usage example:
        Ku = Ku_model()
        Ku.read_config(path_to_config_file)
        Ku.read_input_files()
        Ku.construct_results(vars_to_save = ['active_layer_thickness'])
        Ku.run_all_steps()
        Ku.write_output(path_to_output_file, ['active_layer_thickness'])

    """

    def __init__(self):
        """Initialize the model with default (blank) values.

        Args:
            None
        """

        # Used by read_config()
        self.experiment = ""
        self.inputs_dir = ""
        self.outputs_dir = ""
        self.number_of_years = 0
        self.grid_shape = [0, 0]
        self.input_files = {}
        self.constants = {}
        self.results = {}

        # Used by read_input_files()
        self.air_temperature = None
        self.temperature_amplitude = None
        self.snow_thickness = None
        self.snow_density = None
        self.soil_water_content = None
        self.frozen_vegetation_height = None
        self.thawed_vegetation_height = None
        self.frozen_vegetation_diffusivity = None
        self.thawed_vegetation_diffusivity = None

        # Used by update_soil_heat_capacity()
        self.bulk_thawed_heat_capacity = None
        self.bulk_frozen_heat_capacity = None

        # Used by update_soil_thermal_conductivity()
        self.bulk_thawed_conductivity = None
        self.bulk_frozen_conductivity = None

        # Used by update_snow_thermal_properties()
        self.snow_thermal_conductivity = None
        self.snow_thermal_diffusivity = None

        # Used by update_season_durations()
        self.length_of_cold_season = None
        self.length_of_warm_season = None

        # Used by update_ground_surface_temperature()
        self.snow_insulation = None
        self.snow_damping = None
        self.temperature_at_vegetation = None
        self.amplitude_at_vegetation = None
        self.winter_vegetation_effect = None
        self.summer_vegetation_effect = None
        self.vegetation_insulation = None
        self.vegetation_damping = None
        self.ground_surface_temperature = None
        self.ground_surface_amplitude = None

        # Used by update_permafrost_temperature()
        self.permafrost_temperature = None
        self.soil_conductivity = None
        self.soil_heat_capacity = None

        # Used by update_active_layer()
        self.active_layer_thickness = None
        self.critical_depth = None
        self.permafrost_amplitude = None

    ##############
    # Initialize #
    ##############

    def initialize(self, config_file: str):
        """Initialize the model, read input files, and prepare an object to save results."""
        self.read_config(config_file)
        self.read_input_files()
        self.construct_results()

    def read_config(self, config_file: str):
        """Read the configuration file and populate the indicated attributes.

        Args:
            config_file: str
                Path to the configuration file.
        """

        # Open the configuration file with a TOML parser
        with open(config_file, "r") as file:
            config = toml.load(file)

        # Variables to help when organizing multiple model runs
        self.experiment = config["experiment"]
        self.inputs_dir = config["directories"]["inputs_dir"]
        self.outputs_dir = config["directories"]["outputs_dir"]

        # The domain for this model instance
        self.number_of_years = config["domain"]["number_of_years"]
        self.grid_shape = config["domain"]["grid_shape"]

        # Paths to the input data files
        self.input_files = {var: ncfile for var, ncfile in config["files"].items()}

        # Soil properties
        self.soils = {soil: props for soil, props in config["soils"].items()}

        # Physical constants
        self.constants = {var: val for var, val in config["constants"].items()}

    def read_input_files(self):
        """Read input data files and store fields as instance variables."""

        # Raise an error if there are no input files to read.
        if len(self.input_files) == 0:
            raise ValueError(
                "No input files to read: did you call read_config() first?"
            )

        # Read input data from each file.
        for key, file in self.input_files.items():
            var = key.replace("_file", "")
            data = xr.open_dataarray(self.inputs_dir + file)

            data = self.broadcast(data)

            setattr(self, var, data)

        # Read soil properties, if not defined directly in the config file.
        for soil, props in self.soils.items():
            if len(props["nc_file"]) > 0:
                data = xr.open_dataarray(prop["nc_file"])

                data = self.broadcast(data)

            else:
                data = xr.full_like(self.air_temperature, props["scalar_fraction"])

            self.soils[soil]["fraction"] = data

    def broadcast(self, data: xr.DataArray) -> xr.DataArray:
        """Broadcast an xarray DataArray to a shape that matches the model's domain.

        Returns:
            data:
                The input DataArray, with new dimensions of (x, y, time).
        """
        if data.shape == (self.number_of_years, self.grid_shape[0], self.grid_shape[1]):
            pass

        elif data.shape == (self.grid_shape[0], self.grid_shape[1]):
            data = data.expand_dims({"time": self.number_of_years}, axis=0)

        elif data.shape == (self.number_of_years,):
            data = data.expand_dims(
                {"x": self.grid_shape[1], "y": self.grid_shape[0]}, axis=[1, 2]
            )

        elif data.shape == (1,):
            dim_name = str(data.dims[0])
            data = data.expand_dims(
                {
                    "time": self.number_of_years,
                    "x": self.grid_shape[1],
                    "y": self.grid_shape[0],
                },
                axis=[0, 1, 2],
            )
            data = data.squeeze(dim_name)

        else:
            raise ValueError(
                var
                + " data cannot be broadcast to shape "
                + str((self.number_of_years, self.grid_shape[0], self.grid_shape[1]))
            )

        return data

    def construct_results(
        self,
        vars_to_save: list = ["permafrost_temperature", "active_layer_thickness"],
        template: str = "air_temperature",
    ):
        """Construct a dictionary of empty arrays to store the model results.

        Args:
            vars_to_save: list
                A list of variable names to write out at every time step.
            template: str
                The name of a variable to use as the template for empty arrays.
        """

        self.results = {
            var: np.empty(
                [self.number_of_years, self.grid_shape[0], self.grid_shape[1]]
            )
            for var in vars_to_save
        }

        for var in vars_to_save:
            data = xr.full_like(getattr(self, template), 0.0)
            data = self.broadcast(data)
            setattr(self, var, data)

    ##########
    # Update #
    ##########

    def update_soil_heat_capacity(self, t: int):
        """Calculate the heat capacity of the soil in each column."""

        total_soil_fraction = np.add.reduce(
            [props["fraction"][t, :, :] for soil, props in self.soils.items()]
        )

        weighted_heat_capacity = np.add.reduce(
            [
                props["heat_capacity"]
                * props["fraction"][t, :, :]
                / total_soil_fraction
                for soil, props in self.soils.items()
            ]
        )

        weighted_bulk_density = np.add.reduce(
            [
                props["bulk_density"] * props["fraction"][t, :, :] / total_soil_fraction
                for soil, props in self.soils.items()
            ]
        )

        # Anisimov et al. (1997)
        self.bulk_thawed_heat_capacity = (
            weighted_heat_capacity * weighted_bulk_density
            + 4190.0 * self.soil_water_content[t, :, :]
        )

        self.bulk_frozen_heat_capacity = (
            weighted_heat_capacity * weighted_bulk_density
            + 2025.0 * self.soil_water_content[t, :, :]
        )

    def update_soil_thermal_conductivity(self, t: int):
        """Calculate the thermal conductivity of the soil in each column."""

        total_soil_fraction = np.add.reduce(
            [props["fraction"][t, :, :] for soil, props in self.soils.items()]
        )

        dry_thawed_conductivity = np.multiply.reduce(
            [
                props["conductivity_thawed_dry"]
                ** (props["fraction"][t, :, :] / total_soil_fraction)
                for soil, props in self.soils.items()
            ]
        )

        self.bulk_thawed_conductivity = dry_thawed_conductivity ** (
            1 - self.soil_water_content[t, :, :]
        ) * 0.54 ** (self.soil_water_content[t, :, :])

        dry_frozen_conductivity = np.multiply.reduce(
            [
                props["conductivity_frozen_dry"]
                ** (props["fraction"][t, :, :] / total_soil_fraction)
                for soil, props in self.soils.items()
            ]
        )

        self.bulk_frozen_conductivity = (
            dry_frozen_conductivity ** (1 - self.soil_water_content[t, :, :])
            * 2.35 ** (self.soil_water_content[t, :, :] - self.constants["uwc"])
            * 0.54 ** (self.constants["uwc"])
        )

    def update_snow_thermal_properties(self, t: int):
        """Calculate the thermal conductivity and diffusivity of the snow layer."""

        # Sturm et al. (1997), eq. (4)
        self.snow_thermal_conductivity = (
            0.138
            - 1.01 * (self.snow_density[t, :, :] / 1000)
            + 3.233 * (self.snow_density[t, :, :] / 1000) ** 2
        )

        self.snow_thermal_diffusivity = self.snow_thermal_conductivity / (
            self.snow_density[t, :, :] * self.constants["snow_heat_capacity"]
        )

    def update_season_durations(self, t: int):
        """Estimate the duration of the frozen and thawed seasons."""

        self.length_of_cold_season = self.constants["sec_per_a"] * (
            0.5
            - (1.0 / np.pi)
            * np.arcsin(
                self.air_temperature[t, :, :] / self.temperature_amplitude[t, :, :]
            )
        )
        self.length_of_warm_season = (
            self.constants["sec_per_a"] - self.length_of_cold_season
        )

    def update_ground_surface_temperature(self, t: int):
        """Calculate the temperature at the base of the snow and vegetation layers."""

        # Anisimov et al. (1997), eq. (7)
        inner_eq7 = np.exp(
            -1.0
            * self.snow_thickness[t, :, :]
            * np.sqrt(
                np.pi / (self.constants["sec_per_a"] * self.snow_thermal_diffusivity)
            )
        )
        self.snow_insulation = self.temperature_amplitude[t, :, :] * (1 - inner_eq7)
        self.snow_damping = self.snow_insulation * 2.0 / np.pi
        self.temperature_at_vegetation = (
            self.air_temperature[t, :, :] + self.snow_insulation
        )
        self.amplitude_at_vegetation = (
            self.temperature_amplitude[t, :, :] - self.snow_damping
        )

        # Anisimov et al. (1997), eq. (10)
        inner_eq10 = 1.0 - np.exp(
            -1.0
            * self.frozen_vegetation_height[t, :, :]
            * np.sqrt(
                np.pi
                / (
                    2
                    * self.frozen_vegetation_diffusivity[t, :, :]
                    * self.length_of_cold_season[:, :]
                )
            )
        )
        self.winter_vegetation_effect = (
            self.amplitude_at_vegetation - self.temperature_at_vegetation
        ) * inner_eq10

        # Anisimov et al. (1997), eq. (11)
        inner_eq11 = 1.0 - np.exp(
            -1.0
            * self.thawed_vegetation_height[t, :, :]
            * np.sqrt(
                np.pi
                / (
                    2
                    * self.thawed_vegetation_diffusivity[t, :, :]
                    * self.length_of_warm_season[:, :]
                )
            )
        )
        self.summer_vegetation_effect = (
            self.amplitude_at_vegetation + self.temperature_at_vegetation
        ) * inner_eq11

        self.vegetation_insulation = (
            (
                self.winter_vegetation_effect * self.length_of_cold_season
                + self.summer_vegetation_effect * self.length_of_warm_season
            )
            / self.constants["sec_per_a"]
            * (2.0 / np.pi)
        )

        self.vegetation_damping = (
            self.winter_vegetation_effect * self.length_of_cold_season
            + self.summer_vegetation_effect * self.length_of_warm_season
        ) / self.constants["sec_per_a"]

        self.ground_surface_temperature = (
            self.temperature_at_vegetation + self.vegetation_insulation
        )
        self.ground_surface_amplitude = (
            self.amplitude_at_vegetation - self.vegetation_damping
        )

    def update_permafrost_temperature(self, t: int):
        """Calculate the temperature at the top of the permafrost layer."""

        # Anisimov et al. (1997), eq. (14)
        first_term_eq14 = (
            0.5
            * self.ground_surface_temperature[:, :]
            * (self.bulk_frozen_conductivity + self.bulk_thawed_conductivity)
        )
        second_term_eq14 = (
            self.ground_surface_amplitude
            * (self.bulk_thawed_conductivity - self.bulk_frozen_conductivity)
            / np.pi
        )
        inner_eq14 = (
            self.ground_surface_temperature
            / self.ground_surface_amplitude
            * np.arcsin(self.ground_surface_temperature / self.ground_surface_amplitude)
            + np.sqrt(1.0 - (np.pi**2 / self.ground_surface_amplitude**2))
        )
        numerator_eq14 = first_term_eq14 + second_term_eq14 * inner_eq14

        self.soil_conductivity = np.where(
            numerator_eq14 > 0.0,
            self.bulk_thawed_conductivity,
            self.bulk_frozen_conductivity,
        )

        self.soil_heat_capacity = np.where(
            numerator_eq14 > 0.0,
            self.bulk_thawed_heat_capacity,
            self.bulk_frozen_heat_capacity,
        )

        self.permafrost_temperature = numerator_eq14 / self.soil_conductivity

    def update_active_layer(self, t: int):
        """Calculate the active layer thickness."""

        self.latent_heat = (
            self.constants["latent_heat"] * 1e3 * self.soil_water_content[t, :, :]
        )

        # Romanovsky et al. (1997), eq. (4)
        self.permafrost_amplitude = (
            self.ground_surface_amplitude - np.abs(self.permafrost_temperature)
        ) / np.log(
            (
                self.ground_surface_amplitude
                + self.latent_heat / (2 * self.soil_heat_capacity)
            )
            / (
                np.abs(self.permafrost_temperature)
                + self.latent_heat / (2 * self.soil_heat_capacity)
            )
        ) - self.latent_heat / (
            2 * self.soil_heat_capacity
        )

        # Romanovsky et al. (1997), eq. (5)
        self.critical_depth = (
            2
            * (self.ground_surface_amplitude - np.abs(self.permafrost_temperature))
            * np.sqrt(
                (
                    self.soil_conductivity
                    * self.soil_heat_capacity
                    * self.constants["sec_per_a"]
                )
                / np.pi
            )
            / (
                2 * self.permafrost_amplitude * self.soil_heat_capacity
                + self.latent_heat
            )
        )

        # Romanovsky et al. (1997), eq. (3)
        self.active_layer_thickness = (
            2
            * (self.ground_surface_amplitude - np.abs(self.permafrost_temperature))
            * np.sqrt(
                self.soil_conductivity
                * self.soil_heat_capacity
                * self.constants["sec_per_a"]
                / np.pi
            )
            + (
                2
                * self.permafrost_amplitude
                * self.soil_heat_capacity
                * self.critical_depth
                + self.latent_heat * self.critical_depth
            )
            * self.latent_heat
            * np.sqrt(
                (self.soil_conductivity * self.constants["sec_per_a"])
                / (self.soil_heat_capacity * np.pi)
            )
            / (
                2
                * self.permafrost_amplitude
                * self.soil_heat_capacity
                * self.critical_depth
                + self.latent_heat * self.critical_depth
                + (
                    2 * self.permafrost_amplitude * self.soil_heat_capacity
                    + self.latent_heat
                )
                * np.sqrt(
                    (self.soil_conductivity * self.constants["sec_per_a"])
                    / (self.soil_heat_capacity * np.pi)
                )
            )
        ) / (2 * self.permafrost_amplitude * self.soil_heat_capacity + self.latent_heat)

    def run_one_step(self, t: int):
        """Run one step of the model.

        Args:
            t: int
                The year to run, indicated by its index. All datasets have shape (years, x, y).
        """

        self.update_soil_heat_capacity(t)
        self.update_soil_thermal_conductivity(t)
        self.update_snow_thermal_properties(t)
        self.update_season_durations(t)
        self.update_ground_surface_temperature(t)
        self.update_permafrost_temperature(t)
        self.update_active_layer(t)

        # Save the output fields to a results table.
        for var in self.results.keys():
            self.results[var][t] = getattr(self, var).values[:, :]

    def run_all_steps(self):
        """Run all years of the simulation."""

        for t in range(self.number_of_years):
            self.run_one_step(t)

    ############
    # Finalize #
    ############

    def write_output(self, path: str, vars_to_write: list):
        """Write out a netCDF file with all indicated output variables.

        Args:
            path: str
                The file path where the output file will be written.
            vars_to_write: list
                A list of output variables that should be included in the netCDF file.
        """

        arrays = {}

        for var in vars_to_write:
            data = self.results[var]
            dataarray = xr.DataArray(data=data, dims=["time", "x", "y"])
            arrays[var] = dataarray

        dataset = xr.Dataset(arrays)

        dataset.to_netcdf(path, mode="w", format="NETCDF4")
