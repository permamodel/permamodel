# -*- coding: utf-8 -*-
""" Frost Number by Nelson and Outcalt 1983.
    DOI: 10.2307/1551363. http://www.jstor.org/stable/1551363

    Geo version
"""
from __future__ import print_function

import datetime
import os

import numpy as np
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset

from permamodel import data_directory, examples_directory
from permamodel.components import perma_base

default_frostnumberGeo_config_filename = "FrostnumberGeo_Default.cfg"


class FrostnumberGeoMethod(perma_base.PermafrostComponent):
    # Note: the current version interprets timesteps in years
    def __init__(self, cfgfile=None):
        """Initial definitions and assignments"""
        self._name = "FrostNumberGeo"
        self._configuration = {}

        self._config_description = ""
        self._run_description = ""

        self._calc_surface_fn = False
        self._calc_stefan_fn = False

        self._using_WMT = False
        self._using_Files = False
        self._using_Default = False
        self._using_ConfigVals = False

        self._grid_description_filename = ""
        self._run_duration_filename = ""
        self._temperature_config_filename = ""
        self._temperature_source_filename = ""
        self._soil_properties_config_filename = ""
        self._precipitation_config_filename = ""
        self.output_filename = ""
        self.output = {}

        # For now, we are assuming year-timesteps with nominal date of
        #   12/15/<year>
        self.month = 12
        self.day = 15

        self._output_fid = -1
        self._nc_time = datetime.date(1900, self.month, self.day)

        self._temperature_first_date = datetime.date(1900, self.month, self.day)
        self._temperature_last_date = datetime.date(1900, self.month, self.day)

        self._grid_region = ""
        self._grid_resolution = ""
        self._grid_type = ""
        self._grid_shape = ()
        self._grid_i0 = -1
        self._grid_j0 = -1
        self._grid_i1 = -1
        self._grid_j1 = -1
        self._grid_iskip = -1
        self._grid_jskip = -1
        self._grid_xdim = -1
        self._grid_ydim = -1

        self._reference_date = datetime.date(1900, self.month, self.day)
        self._start_date = datetime.date(1900, self.month, self.day)
        self._end_date = datetime.date(1900, self.month, self.day)
        self._timestep_duration = -1

        self._temperature_dates = []
        self._temperature_datacube = []
        self._precipitation_dates = []
        self._precipitation_datacube = []
        self._soilproperties_dates = []
        self._soilproperties_datacube = []

        self.T_air = []
        self.T_air_min = np.zeros([1])
        self.T_air_max = np.zeros([1])

        self._temperature_current = np.zeros([1])
        self.Precip = np.zeros([1])
        self.SoilProperties = np.zeros([1])

        self._dd_method = self.__init__
        self.ddf = np.zeros([1])
        self.ddt = np.zeros([1])
        self.air_frost_number_Geo = np.zeros([1])
        self.surface_frost_number_Geo = np.zeros([1])
        self.stefan_frost_number_Geo = np.zeros([1])

        self._nc_time = datetime.date(1900, self.month, self.day)
        self._nc_reference_time = datetime.date(1900, self.month, self.day)
        self._nc_last_time_index = 0
        self._nc_x = np.zeros([1])
        self._nc_y = np.zeros([1])

        self._nc_afn = np.zeros([1])
        self._nc_sfn = np.zeros([1])
        self._nc_stfn = np.zeros([1])
        self.surface_frost_number = np.float32(-1.0)
        self.stefan_frost_number = np.float32(-1.0)

        self._timestep_first = -1
        self._timestep_last = -1
        self._timestep_current = -1
        self._date_current = datetime.date(1900, self.month, self.day)

        self._temperature_dataset = np.zeros([1])

        # There is a default configuration file in the examples directory
        if cfgfile is None:
            self._config_filename = os.path.join(
                examples_directory, default_frostnumberGeo_config_filename
            )
        else:
            self._config_filename = cfgfile

    def initialize_frostnumberGeo_component(self):
        # Read in the overall configuration from the configuration file
        assert os.path.isfile(self._config_filename)
        self._configuration = self.get_config_from_oldstyle_file(self._config_filename)
        # self.get_config_from_yaml_file(self._config_filename)
        # Ensure that this config file is for this type of Method
        assert (
            self._configuration["config_for_method"]
            == self.__class__.__name__.split(".")[-1]
        )

        # Set some model description strings
        self._config_description = self._configuration["config_description"]
        self._run_description = self._configuration["run_description"]

        # Determine whether surface and stefan numbers will be generated
        self._calc_surface_fn = (
            self._configuration["calc_surface_frostnumber"] == "True"
        )
        self._calc_stefan_fn = self._configuration["calc_stefan_frostnumber"] == "True"

        # This model can be run such that input variables are either
        #   read from Files
        #   provided by WMT or
        #   completely described in the config file (default)
        #
        #   This is set by the configuration value of 'input_var_source'
        #   as "Files", "WMT", or "Default".

        # If this is "Files", then which of the frost numbers will be
        # calculated is determined by which filenames are given.
        # Determine which variables will be calculable
        if self._configuration["input_var_source"] == "Files":
            self._using_WMT = False
            self._using_Files = True
            self._using_ConfigVals = False

            # Read information about the grid
            # Either from a file or directly
            if "grid_description_filename" in self._configuration.keys():
                self._grid_description_filename = os.path.join(
                    examples_directory, self._configuration["grid_description_filename"]
                )
                grid_config = self.get_config_from_yaml_file(
                    self._grid_description_filename
                )
                self._grid_region = grid_config["grid_region"]
                self._grid_resolution = grid_config["grid_resolution"]
                self._grid_type = grid_config["grid_type"]
                self._grid_shape = grid_config["grid_shape"]
                self._grid_i0 = grid_config["i_ul"]
                self._grid_j0 = grid_config["j_ul"]
                self._grid_iskip = grid_config["i_skip"]
                self._grid_jskip = grid_config["j_skip"]
            else:
                self._grid_region = self._configuration["grid_region"]
                self._grid_resolution = self._configuration["grid_resolution"]
                self._grid_type = self._configuration["grid_type"]
                self._grid_shape = self._configuration["grid_shape"]
                self._grid_i0 = self._configuration["i_ul"]
                self._grid_j0 = self._configuration["j_ul"]
                self._grid_iskip = self._configuration["i_skip"]
                self._grid_jskip = self._configuration["j_skip"]

            # Read information about the model run
            # Either from a file or directly
            if "run_duration_filename" in self._configuration.keys():
                self._run_duration_filename = os.path.join(
                    examples_directory, self._configuration["run_duration_filename"]
                )
                run_config = self.get_config_from_yaml_file(self._run_duration_filename)
                self._reference_date = run_config["model_reference_date"]
                self._start_date = run_config["model_start_date"]
                self._end_date = run_config["model_end_date"]
                self._timestep_duration = run_config["model_timestep"]
            else:
                self._reference_date = self._configuration["model_reference_date"]
                self._start_date = self._configuration["model_start_date"]
                self._end_date = self._configuration["model_end_date"]
                self._timestep_duration = self._configuration["model_timestep"]

            # Set up the grid subset
            if len(self._grid_shape) == 2:
                (self._grid_ydim, self._grid_xdim) = self._grid_shape
            else:
                raise ValueError(
                    "cannot handle grid of shape %s" % str(self._grid_shape)
                )
            self._grid_i1 = self._grid_i0 + self._grid_iskip * self._grid_xdim
            self._grid_j1 = self._grid_j0 + self._grid_jskip * self._grid_ydim

            # At a minimum, there must be temperature information
            if "temperature_config_filename" in self._configuration.keys():
                self._temperature_config_filename = os.path.join(
                    examples_directory,
                    self._configuration["temperature_config_filename"],
                )
                assert os.path.isfile(self._temperature_config_filename)
                temperature_configuration = self.get_config_from_yaml_file(
                    self._temperature_config_filename
                )
                self._temperature_source_filename = os.path.join(
                    data_directory,
                    temperature_configuration["temperature_source_filename"],
                )
                self._temperature_first_date = temperature_configuration[
                    "dataset_first_date"
                ]
                self._temperature_last_date = temperature_configuration[
                    "dataset_last_date"
                ]
            else:
                self._temperature_source_filename = os.path.join(
                    data_directory, self._configuration["temperature_source_filename"]
                )
                self._temperature_first_date = self._configuration["dataset_first_date"]
                self._temperature_last_date = self._configuration["dataset_last_date"]

            if self._configuration["precipitation_config_filename"] is not None:
                self._precipitation_config_filename = os.path.join(
                    examples_directory,
                    self._configuration["precipitation_config_filename"],
                )
                assert os.path.isfile(self._precipitation_config_filename)
                self._calc_surface_fn = True
            else:
                self._calc_surface_fn = False
                self._calc_stefan_fn = False

            if (
                self._calc_surface_fn
                and self._configuration["soil_properties_config_filename"] is not None
            ):
                self._soil_properties_config_filename = os.path.join(
                    examples_directory,
                    self._configuration["soil_properties_config_filename"],
                )
                assert os.path.isfile(self._soil_properties_config_filename)
                self._calc_stefan_fn = True
            else:
                self._calc_stefan_fn = False

            self.initialize_input_vars_from_files()

        # If this is "WMT", then a flag must be set in the config file
        # to determine which
        elif self._configuration["input_var_source"] == "WMT":
            self._using_WMT = True
            self._using_Files = False
            self._using_ConfigVals = False

            # Run duration
            self._reference_date = self._configuration["model_reference_date"]
            self._start_date = self._configuration["model_start_date"]
            self._end_date = self._configuration["model_end_date"]
            self._timestep_duration = self._configuration["model_timestep"]

            # Grid shape
            self._grid_type = self._configuration["grid_type"]
            self._grid_shape = self._configuration["grid_shape"]
            if len(self._grid_shape) == 2:
                (self._grid_ydim, self._grid_xdim) = self._grid_shape
            else:
                raise ValueError(
                    "cannot handle grid of shape %s" % str(self._grid_shape)
                )

            # WMT provides an array called self._temperature_current[]
            # which will be set elsewhere and used in get_
            # ... this is set after this mode check so the code can
            #     pass the bmitester ...

        # If initialized completely from the config file, this is 'Default'
        elif self._configuration["input_var_source"] == "Default":
            self._using_WMT = False
            self._using_Files = False
            self._using_ConfigVals = True

            # 'Default' has the same values as 'WMT',
            # but it also specifies a 3D (x, y, time) data array
            # from which the grids will derive their values

            # Run duration
            self._reference_date = self._configuration["model_reference_date"]
            self._start_date = self._configuration["model_start_date"]
            self._end_date = self._configuration["model_end_date"]
            self._timestep_duration = self._configuration["model_timestep"]

            # Grid shape
            self._grid_type = self._configuration["grid_type"]
            self._grid_shape = self._configuration["grid_shape"]
            if len(self._grid_shape) == 2:
                (self._grid_ydim, self._grid_xdim) = self._grid_shape
            else:
                raise ValueError(
                    "cannot handle grid of shape %s" % str(self._grid_shape)
                )

            # Set the valid dates for which temperature is available
            self._temperature_first_date = self.datefrom(
                self._configuration["temperature_grid_date_0"]
            )
            last_date_varname = "temperature_grid_date_%d" % (
                int(self._configuration["n_temperature_grid_fields"]) - 1
            )
            self._temperature_last_date = eval(
                "self.datefrom(self._configuration['%s'])" % last_date_varname
            )

            # Parse the (x, y, time) data cube from which the values
            # will be drawn at the appropriate time
            (
                self._temperature_dates,
                self._temperature_datacube,
            ) = self.initialize_datacube("temperature", self._configuration)
            if self._configuration["n_precipitation_grid_fields"] > 0:
                (
                    self._precipitation_dates,
                    self._precipitation_datacube,
                ) = self.initialize_datacube("precipitation", self._configuration)
            if self._configuration["n_soilproperties_grid_fields"] > 0:
                (
                    self._soilproperties_dates,
                    self._soilproperties_datacube,
                ) = self.initialize_datacube("soilproperties", self._configuration)

        # Initialize the temperature min/max arrays
        self.T_air_min = np.zeros(self._grid_shape, dtype=np.float32)
        self.T_air_min.fill(np.nan)
        self.T_air_max = np.zeros(self._grid_shape, dtype=np.float32)
        self.T_air_max.fill(np.nan)

        # Initialize the current temperature array
        # Even though this is only used in WMT-mode, it needs to be set
        # generally to pass the bmi-tester
        self._temperature_current = np.zeros(self._grid_shape, dtype=np.float32)
        self._temperature_current.fill(np.nan)

        # Initialize the Jan and Jul arrays
        self._temperature_jan = np.zeros(self._grid_shape, dtype=np.float32)
        self._temperature_jan.fill(np.nan)
        self._temperature_jul = np.zeros(self._grid_shape, dtype=np.float32)
        self._temperature_jul.fill(np.nan)

        # Initialize the Precip and SoilProperites arrays if needed
        if self._calc_surface_fn:
            self.Precip = np.zeros(self._grid_shape, dtype=np.float32)
            self.Precip.fill(np.nan)
        else:
            self.Precip = None

        if self._calc_stefan_fn:
            self.SoilProperties = np.zeros(self._grid_shape, dtype=np.float32)
            self.SoilProperties.fill(np.nan)
        else:
            self.SoilProperties = None

        # There are different ways of computing degree days.  Ensure that
        # the method specified has been coded
        self._dd_method = self._configuration["degree_days_method"]
        self.ddf = np.zeros(self._grid_shape, dtype=np.float32)
        self.ddf.fill(np.nan)
        self.ddt = np.zeros(self._grid_shape, dtype=np.float32)
        self.ddt.fill(np.nan)

        # Initialize the output grids
        self.air_frost_number_Geo = np.zeros(self._grid_shape, dtype=np.float32)
        if self._calc_surface_fn:
            self.surface_frost_number_Geo = np.zeros(self._grid_shape, dtype=np.float32)
        else:
            self.surface_frost_number_Geo = None
        if self._calc_stefan_fn:
            self.stefan_frost_number_Geo = np.zeros(self._grid_shape, dtype=np.float32)
        else:
            self.stefan_frost_number_Geo = None

        # Initialize the time information from a configuration file
        self.initialize_model_time()

        # Initialize the output filename
        self.initialize_output(
            self._configuration["output_directory"],
            self._configuration["output_filename"],
        )

        """
        # Note: initialization only allocates for values; it doesn't
        #       calculate the first timestep's values

        # This is to allow another routine, e.g. from WMT, to set the
        # necessary values
        """

    def datefrom(self, datestring):
        if isinstance(datestring, str):
            return datetime.datetime.strptime(datestring, "%Y-%m-%d").date()
        else:
            return datestring

    def initialize_datacube(self, gridname, config):
        # Determine the number of lines for this grid
        ngridlines = 0
        thisdate = datetime.date(
            1900, self.month, self.day
        )  # This is a dummy init value
        ngridlines = config["n_{0}_grid_fields".format(gridname)]

        # Create the datelist for this grid
        datelist = []
        for n in range(ngridlines):
            # This version converts a string to a datetime.date object
            # exec("thisdate = datetime.datetime.strptime(
            # config['%s_grid_date_%d' % (gridname, n)], '%Y-%m-%d').date()")
            thisdate = config["{0}_grid_date_{1}".format(gridname, n)]
            if isinstance(thisdate, str):
                thisdate = self.datefrom(thisdate)
            datelist.append(thisdate)

        # Create the datacube for this grid
        slicelist = []
        for n in range(ngridlines):
            config_arg = eval(config["%s_grid_data_%d" % (gridname, n)])
            slicelist.append(np.array(config_arg))
        datacube = np.stack(slicelist)

        return datelist, datacube

    def get_datacube_slice(self, thisdate, datacube, cubedates):
        # Data are invalid until the first date in the datacube
        if thisdate < cubedates[0]:
            raise ValueError(
                "Date %s is before first valid date (%s) in datacube"
                % (thisdate, cubedates[0])
            )

        if thisdate > cubedates[-1]:
            raise ValueError(
                "Date %s is after last valid date (%s) in datacube"
                % (thisdate, cubedates[-1])
            )

        use_this_date = 0
        for date_okay, checkdate in enumerate(cubedates):
            if thisdate > checkdate:
                use_this_date = date_okay
            else:
                break

        return datacube[use_this_date, :]

    def initialize_output(self, outdirname, outfilename):
        """Initialize the output file

        There are some substitutions that can be made in the
        dirname and filename

        dirname:
            EXAMPLES:  use the .../examples/ subdirectory
            DATA    :  use the .../data/ subdirectory
        filename:
            STARTDATE:  replace this with the start date
            ENDDATE:    replace this with the end date
            TIMESTEP:   replace this with the timestep (in days)
        """
        # Set the output directory name
        if outdirname is None:
            outdirname = "."
        elif outdirname == "DATA":
            outdirname = data_directory
        elif outdirname == "EXAMPLES":
            outdirname = examples_directory
        assert os.path.isdir(outdirname)

        # Replace substitutions in filename
        outfilename = outfilename.replace("STARTDATE", str(self._start_date))
        outfilename = outfilename.replace("ENDDATE", str(self._end_date))
        outfilename = outfilename.replace("TIMESTEP", str(self._timestep_duration))

        self.output_filename = os.path.join(outdirname, outfilename)

        # Open the output netcdf file
        # Note: Be sure to close this, e.g. in self.finalize(), because
        # otherwise multiple openings, e.g. when testing, may fail
        try:
            self._output_fid = Dataset(self.output_filename, "w", format="NETCDF4")
        except IOError:
            # Note: fid will have type integer if it fails to open
            print(
                "WARNING: Could not open file for output: {}".format(
                    self.output_filename
                )
            )

        if self._output_fid != -1:
            # Add attributes
            #  including comments from config files
            setattr(
                self._output_fid,
                "company",
                "Community Surface Dynamics Modeling System",
            )
            setattr(self._output_fid, "Permafrost Component", "FrostnumberGeo")

            ### Init dimensions
            # Time dimension
            tdim = 0
            self._output_fid.createDimension("time", tdim)
            self._nc_time = self._output_fid.createVariable(
                "time", "i", ("time",), zlib=True
            )
            setattr(self._nc_time, "time_long_name", "time")
            setattr(self._nc_time, "time_standard_name", "time")
            setattr(self._nc_time, "time_units", "months since 1900-01-01 00:00:00")
            self._nc_reference_time = datetime.date(1900, 1, 15)
            setattr(self._nc_time, "time_format", "modified julian day (MJD)")
            setattr(self._nc_time, "time_time_zone", "UTC")
            setattr(self._nc_time, "time__FillValue", "-9999")
            self._nc_last_time_index = 0

            ### For now, the X and Y dimensions are just the indexes
            # X dimension
            self._output_fid.createDimension("x", self._grid_xdim)
            self._nc_x = self._output_fid.createVariable("x", "f", ("x",), zlib=True)
            setattr(self._nc_x, "x_long_name", "projected x direction")
            setattr(self._nc_x, "x_standard_name", "x")
            setattr(self._nc_x, "x_units", "meters")
            setattr(self._nc_x, "x__FillValue", "NaN")
            # fill the x- values
            for x in range(self._grid_xdim):
                self._nc_x[x] = x

            # Y dimension
            self._output_fid.createDimension("y", self._grid_ydim)
            self._nc_y = self._output_fid.createVariable("y", "f", ("y",), zlib=True)
            setattr(self._nc_y, "y_long_name", "projected y direction")
            setattr(self._nc_y, "y_standard_name", "y")
            setattr(self._nc_y, "y_units", "meters")
            setattr(self._nc_y, "y__FillValue", "NaN")
            # fill the y- values
            for y in range(self._grid_ydim):
                self._nc_y[y] = y

            ### Init grids with sizes
            # Allocate air frost number field
            self._nc_afn = self._output_fid.createVariable(
                "air_fn", "f", ("time", "y", "x"), zlib=True
            )
            setattr(self._nc_afn, "afn_long_name", "Air Frost Number")
            setattr(self._nc_afn, "afn_standard_name", "Frostnumber_air")
            setattr(self._nc_afn, "afn_units", "none")
            setattr(self._nc_afn, "afn__FillValue", "-99")

            # Allocate surface frost number field, if computing it
            if self._calc_surface_fn:
                self._nc_sfn = self._output_fid.createVariable(
                    "surface_fn", "f", ("time", "y", "x"), zlib=True
                )
                setattr(self._nc_sfn, "sfn_long_name", "Surface Frost Number")
                setattr(self._nc_sfn, "sfn_standard_name", "Frostnumber_surface")
                setattr(self._nc_sfn, "sfn_units", "none")
                setattr(self._nc_sfn, "sfn__FillValue", "-99")

            # Allocate Stefan frost number field, if computing it
            if self._calc_stefan_fn:
                self._nc_stfn = self._output_fid.createVariable(
                    "stefan_fn", "f", ("time", "y", "x"), zlib=True
                )
                setattr(self._nc_stfn, "stfn_long_name", "Stefan Frost Number")
                setattr(self._nc_stfn, "stfn_standard_name", "Frostnumber_stefan")
                setattr(self._nc_stfn, "stfn_units", "none")
                setattr(self._nc_stfn, "stfn__FillValue", "-99")

    def finalize(self):
        # Define this so we don't call the permamodel base class version
        self.finalize_frostnumber_Geo()

    def finalize_frostnumber_Geo(self):
        if self._output_fid != -1:
            self._output_fid.close()

    def check_whether_output_timestep(self, this_timestep):
        # Only output on the 15th of each month
        this_date = self.get_date_from_timestep(this_timestep)
        return this_date.day == 15

    def add_to_output(self):
        do_add = self.check_whether_output_timestep(self._timestep_current)
        if do_add and self._output_fid != -1:
            # I think this assumes that the file has the same time reference as
            # the model run...
            this_date = self.get_date_from_timestep(self._timestep_current)
            years = this_date.year - self._nc_reference_time.year
            months = this_date.month - self._nc_reference_time.month
            time_index = 12 * years + months
            # print("Adding output at index %d for date %s" % \
            #       (time_index, str(this_date)))
            self._nc_afn[time_index, :] = self.air_frost_number_Geo[:]
            if self._calc_surface_fn:
                self._nc_sfn[time_index, :] = self.surface_frost_number_Geo[:]
            if self._calc_stefan_fn:
                self._nc_stfn[time_index, :] = self.stefan_frost_number_Geo[:]

    def initial_update(self):
        # Increment the model for the first time step
        # This is a separate function because the input_vars may be
        # set by WMT and this can't be called until that is done
        self.get_input_vars()
        self.compute_degree_days()
        self.calculate_frost_numbers_Geo()
        self.add_to_output()

    def update(self, frac=None):
        """Compute the values for the current time, then update the time"""
        years_change = 0
        if frac is not None:
            # print("Fractional times not yet permitted, rounding to nearest int")
            years_change = self._timestep_duration * int(frac + 0.5)
        else:
            years_change = self._timestep_duration

        self._date_current += relativedelta(years=years_change)

        self._timestep_current = self.get_timestep_from_date(self._date_current)

        self.get_input_vars()
        self.compute_degree_days()
        self.calculate_frost_numbers_Geo()
        self.add_to_output()

    def update_until_timestep(self, stop_timestep):
        while self._timestep_current < stop_timestep:
            self.update()

    def initialize_model_time(self):
        # The model run duration configuration file has information
        # about the reference time, the start and end times, and
        # the timestep
        # Determine the first and last timesteps, ensure ordering!
        self._timestep_first = self.get_timestep_from_date(self._start_date)
        self._timestep_last = self.get_timestep_from_date(self._end_date)
        assert np.all(self._timestep_last >= self._timestep_first)

        # Set the current timestep to the first timestep
        self.set_current_date_and_timestep_with_timestep(self._timestep_first)

    def set_current_date_and_timestep_with_timestep(self, this_timestep):
        self._timestep_current = this_timestep
        self._date_current = self.get_date_from_timestep(self._timestep_current)

    def get_date_from_timestep(self, timestep):
        return self._reference_date + relativedelta(
            years=timestep * self._timestep_duration
        )

    def get_timestep_from_date(self, this_date):
        # If the timestep is an integer number of days, this would be fine:
        # return (this_date-self._reference_date).days /\
        #    self._timestep_duration.days
        # Compute from seconds:
        # return int((this_date-self._reference_date).total_seconds() / \
        #           (self._timestep_duration.total_seconds()) + 0.5)
        if this_date.month != self.month:
            raise RuntimeWarning(
                "Current month (%d) is not model month(%d)"
                % (this_date.month, self.month)
            )
        if this_date.day != self.day:
            raise RuntimeWarning(
                "Current day (%d) is not model day(%d)" % (this_date.day, self.day)
            )
        return this_date.year - self._reference_date.year

    def initialize_input_vars_from_files(self):
        # If the model does not have its input variables set by an
        # external runner, e.g. WMT, then the values will have to
        # be read from data files here.  This routine sets up the
        # files needed to do this.
        # Currently, this is done by reading values from a netcdf
        # file with the appropriate monthly values for
        self._temperature_dataset = Dataset(
            self._temperature_source_filename, "r", mmap=True
        )
        assert self._temperature_dataset is not None

    def get_input_vars(self):
        if self._using_WMT:
            # With WMT, the input variables will be set externally via BMI
            # Here, we assume the values of the preceding year will be in
            # the self._temperature_[jan|jul] arrays.  And we will interpret
            # January (index 0) as coldest, and July (index 6) as warmest
            if self._dd_method == "MinJanMaxJul":
                self.T_air_min = self._temperature_jan
                self.T_air_max = self._temperature_jul
            else:
                raise ValueError(
                    "Degree days method %s not recognized" % self._dd_method
                )
            return
        elif self._using_Files or self._using_ConfigVals:
            # In standalone mode, variables must be set locally
            # All frost number types need temperature data
            if self._dd_method == "MinJanMaxJul":
                # For the MinJanMaxJul method, need min and max temp fields
                (mindate, maxdate) = self.get_min_and_max_dates(self._date_current)
                self.T_air_min = self.get_temperature_field(mindate)
                self.T_air_max = self.get_temperature_field(maxdate)
            else:
                raise ValueError(
                    "Degree days method %s not recognized" % self._dd_method
                )
        else:
            raise ValueError(
                "Frostnumber must use either Files, ConfigVals \
                              or WMT to get input variables"
            )

    def get_temperature_field(self, t_date=None):
        # By default, return the temperature field at the date of the current
        # timestep
        if t_date is None:
            t_date = self.get_date_from_timestep(self._timestep_current)

        if (
            t_date >= self._temperature_first_date
            and t_date <= self._temperature_last_date
        ):
            t_index = self.get_temperature_month_index(t_date)

            if self._using_Files:
                # Files uses netcdf input
                temperature_subregion = self._temperature_dataset.variables["temp"][
                    t_index
                ][
                    self._grid_j0 : self._grid_j1 : self._grid_jskip,
                    self._grid_i0 : self._grid_i1 : self._grid_iskip,
                ]
            elif self._using_ConfigVals:
                # ConfigVals is the Default, where the grids are in cfg file
                temperature_subregion = self.get_datacube_slice(
                    t_date, self._temperature_datacube, self._temperature_dates
                ).astype(np.float32)
            elif self._using_WMT:
                # Note: I don't think this functionality is currently
                # used because get_input_vars doesn't call this routine
                # for WMT
                # The standalone versions read from a data file
                # or a datacube.
                #
                # Using WMT, this value will be set elsewhere -- via
                # BMI functions -- which will modify the array
                # self._temperature_current[]
                temperature_subregion = self._temperature_current

            assert np.all(temperature_subregion.shape == self._grid_shape)

            # Fill the field with Nan if value is 'missing value'
            nan_locations = temperature_subregion < -90
            temperature_subregion[nan_locations] = np.nan
        else:
            # print("Date is outside valid date range of temperature data")
            temperature_subregion = np.zeros(self._grid_shape, dtype=np.float32)
            temperature_subregion.fill(np.nan)

        return temperature_subregion

    def get_min_and_max_dates(self, this_date):
        """Return the 15th of the month
        of the nearest preceding January and July"""

        if this_date.month >= 7:
            maxyear = this_date.year
        else:
            maxyear = this_date.year - 1

        return (datetime.date(this_date.year, 1, 15), datetime.date(maxyear, 7, 15))

    def get_temperature_month_index(self, t_date):
        """Returns the time index to the temperature dataset for this date"""
        year_offset = t_date.year - self._temperature_first_date.year
        month_offset = t_date.month - self._temperature_first_date.month
        t_index = 12 * year_offset + month_offset
        return t_index

    def get_config_from_yaml_file(self, cfg_filename):
        print("Ignoring %s" % cfg_filename)
        raise RuntimeError("config from YAML not currently supported")
        """
        cfg_struct = None
        try:
            with open(cfg_filename, 'r') as cfg_file:
                cfg_struct = yaml.load(cfg_file)
        except:
            print("\nError opening configuration file in\
                  get_config_from_yaml_file()")
            raise

        return cfg_struct
        """

    def get_config_from_oldstyle_file(self, cfg_filename):
        cfg_struct = {}
        grid_struct = {}
        try:
            with open(cfg_filename, "r") as cfg_file:
                # this was originally modeled after read_config_file()
                # in BMI_base.py as modified for cruAKtemp.py
                while True:
                    # Read lines from config file until no more remain
                    line = cfg_file.readline()
                    if line == "":
                        break

                    # Comments start with '#'
                    COMMENT = line[0] == "#"

                    words = line.split("|")
                    if (len(words) == 4) and (not COMMENT):
                        var_name = words[0].strip()
                        value = words[1].strip()
                        var_type = words[2].strip()

                        # Process the variables based on variable name
                        if var_name[-4:] == "date":
                            # date variables end with "_date"
                            # Note: these should be years
                            assert int(value) <= 2100
                            assert int(value) >= 1800
                            cfg_struct[var_name] = datetime.date(
                                int(value), self.month, self.day
                            )
                        elif var_name[0:4] == "grid":
                            # grid variables are processed after cfg file read
                            grid_struct[var_name] = value
                        elif var_name == "timestep" or var_name == "model_timestep":
                            # timestep is a number of years
                            cfg_struct[var_name] = int(value)
                        elif var_type == "int":
                            # Convert integers to int
                            cfg_struct[var_name] = int(value)
                        else:
                            # Everything else is just passed as a string
                            assert var_type == "string"
                            cfg_struct[var_name] = value

        except:
            print(
                "\nError opening configuration file in\
                  get_config_from_yaml_file()"
            )
            raise

        # Process the grid information
        # I think I had rows and columns switched in cruAKtemp!
        # cfg_struct['grid_shape'] = (int(grid_struct['grid_columns']),
        #                            int(grid_struct['grid_rows']))
        cfg_struct["grid_shape"] = (
            int(grid_struct["grid_rows"]),
            int(grid_struct["grid_columns"]),
        )
        cfg_struct["grid_type"] = grid_struct["grid_type"]

        # for keyname in cfg_struct.keys():
        #    print(keyname)
        cfg_struct["grids"] = {"temperature": "np.float"}
        if cfg_struct["n_precipitation_grid_fields"] > 0:
            cfg_struct["grids"] = {"precipitation": "np.float"}
            self._calc_surface_fn = True
        else:
            self._calc_surface_fn = False
        if cfg_struct["n_soilproperties_grid_fields"] > 0:
            cfg_struct["grids"] = {"soilproperties": "np.float"}
            self._calc_stefan_fn = True
        else:
            self._calc_stefan_fn = False

        return cfg_struct

    def calculate_frost_numbers_Geo(self):
        # Calculate all the frost numbers using the current data
        self.calculate_air_frost_number_Geo()
        self.calculate_surface_frost_number_Geo()
        self.calculate_stefan_frost_number_Geo()

        # Add these frost numbers to the output dictionary
        # self.update_output()

    def print_frost_numbers(self, year=-1):
        # if year is -1, then use the current year of self
        # otherwise, use the specified year
        if year > 0:
            print(
                "Year: %d  F_air=%5.3f  F_surface=%5.3f  F_stefan=%5.3f"
                % (
                    year,
                    self.air_frost_number_Geo,
                    self.surface_frost_number_Geo,
                    self.stefan_frost_number_Geo,
                )
            )
        else:
            for year in sorted(self.output.keys()):
                print("Year: %d  output=%s" % (year, self.output[year]))

    def calculate_air_frost_number_Geo(self):
        self.compute_degree_days()
        self.compute_air_frost_number_Geo()

    def calculate_surface_frost_number_Geo(self):
        # For now, a dummy value
        self.surface_frost_number = np.float32(-1.0)

    def calculate_stefan_frost_number_Geo(self):
        self.stefan_frost_number = np.float32(-1.0)

    def compute_degree_days(self):
        if self._dd_method == "MinJanMaxJul":
            self.compute_degree_days_MinJanMaxJul()

    def compute_degree_days_MinJanMaxJul(self):
        # Input: T_hot (avg temp of warmest month)
        #        T_cold (avg temp of coldest month)

        # Output: ddf (degree freezing days)
        #         ddt (degree thawing days)

        # Set to Nan where appropriate
        for (j, i), maxvalue in np.ndenumerate(self.T_air_max):
            minvalue = self.T_air_min[j, i]

            if minvalue == np.nan or maxvalue == np.nan:
                # Any non-numbers invalidate the dd calcs
                self.ddf[j, i] = np.nan
                self.ddt[j, i] = np.nan
            elif maxvalue < minvalue:
                # Can't have min temp > max temp!
                self.ddf[j, i] = np.nan
                self.ddt[j, i] = np.nan
            elif minvalue >= 0.0:
                # Never freezes
                self.ddf[j, i] = 0.0
                self.ddt[j, i] = 365.0 * (minvalue + maxvalue) / 2.0
            elif maxvalue <= 0.0:
                # Never thaws
                self.ddf[j, i] = -365.0 * (minvalue + maxvalue) / 2.0
                self.ddt[j, i] = 0.0
            else:
                # Freezes in winter, thaws in summer
                T_average = (minvalue + maxvalue) / 2.0
                T_amplitude = (maxvalue - minvalue) / 2.0
                Beta = np.arccos(-T_average / T_amplitude)
                T_summer = T_average + T_amplitude * np.sin(Beta) / Beta
                T_winter = T_average - T_amplitude * np.sin(Beta) / (np.pi - Beta)
                L_summer = 365.0 * Beta / np.pi
                L_winter = 365.0 - L_summer
                self.ddt[j, i] = T_summer * L_summer
                self.ddf[j, i] = -T_winter * L_winter

            if (self.ddt[j, i] == 0.0) and (self.ddf[j, i] == 0.0):
                # This shouldn't happen with real values
                self.ddf[j, i] = np.nan
                self.ddt[j, i] = np.nan

    def compute_air_frost_number_Geo(self):
        # Calculating Reduced Air Frost Number (pages 280-281).
        # The reduced frost number is close 0 for long summers
        #   and close to 1 for long winters.
        # self.air_frost_number_Geo = np.sqrt(self.ddf) / \
        #    (np.sqrt(self.ddf) + np.sqrt(self.ddt))
        where_nan = np.isnan(self.ddf + self.ddt)
        where_notnan = np.logical_not(np.isnan(self.ddf + self.ddt))

        self.air_frost_number_Geo[where_nan] = np.nan

        self.air_frost_number_Geo[where_notnan] = np.sqrt(self.ddf[where_notnan]) / (
            np.sqrt(self.ddf[where_notnan]) + np.sqrt(self.ddt[where_notnan])
        )


if __name__ == "__main__":
    # Run the FrostnumberGeo model
    # Currently, this just runs the defaults
    fn_geo = FrostnumberGeoMethod()
    fn_geo.initialize_frostnumberGeo_component()
    fn_geo.initial_update()
    fn_geo.update_until_timestep(fn_geo._timestep_last)
    fn_geo.finalize()
