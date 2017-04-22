# -*- coding: utf-8 -*-
""" Frost Number by Nelson and Outcalt 1983.
    DOI: 10.2307/1551363. http://www.jstor.org/stable/1551363

    Geo version
"""

import numpy as np
import ast
from permamodel.utils import model_input
from permamodel.components import perma_base
from permamodel import examples_directory, data_directory
import os
import yaml
from nose.tools import (assert_is_instance, assert_greater_equal,
                        assert_less_equal, assert_almost_equal,
                        assert_greater, assert_in, assert_true,
                        assert_equal)
import datetime
from netCDF4 import Dataset

default_frostnumberGeo_config_filename = "FrostnumberGeo_Default.cfg"

class FrostnumberGeoMethod( perma_base.PermafrostComponent ):
    def __init__(self, cfgfile=None):
        # There is a default configuration file in the examples directory
        if cfgfile is None:
            self._config_filename = \
                    os.path.join(examples_directory,
                    default_frostnumberGeo_config_filename)
        else:
            self._config_filename = cfgfile

    def initialize_frostnumberGeo_component(self):
        SILENT = True

        # Note: Initialized from initialize() in perma_base.py
        if not SILENT:
            print("Initializing for FrostnumberGeoMethod")

        self._model = 'FrostNumberGeo'

        # Read in the overall configuration from the configuration file
        assert_true(os.path.isfile(self._config_filename))
        self._configuration = \
                self.get_config_from_oldstyle_file(self._config_filename)
                # self.get_config_from_yaml_file(self._config_filename)
        # Ensure that this config file is for this type of Method
        assert_equal(self._configuration['config_for_method'],
                     str(self.__class__).split('.')[-1])

        # Set some model description strings
        self._config_description = self._configuration['config_description']
        self._run_description = self._configuration['run_description']

        # Determine whether surface and stefan numbers will be generated
        self._calc_surface_fn = \
            self._configuration['calc_surface_frostnumber'] == 'True'
        self._calc_stefan_fn = \
            self._configuration['calc_stefan_frostnumber'] == 'True'

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
        if self._configuration['input_var_source'] == 'Files':
            self._using_WMT = False
            self._using_Files = True
            self._using_ConfigVals = False

            # Read information about the grid
            # Either from a file or directly
            if 'grid_description_filename' in self._configuration.keys():
                self._grid_description_filename = \
                    os.path.join(examples_directory,
                    self._configuration['grid_description_filename'])
                grid_config = \
                    self.get_config_from_yaml_file(
                        self._grid_description_filename)
                self._grid_region = grid_config['grid_region']
                self._grid_resolution = grid_config['grid_resolution']
                self._grid_type = grid_config['grid_type']
                self._grid_shape = grid_config['grid_shape']
                self._grid_i0 = grid_config['i_ul']
                self._grid_j0 = grid_config['j_ul']
                self._grid_iskip = grid_config['i_skip']
                self._grid_jskip = grid_config['j_skip']
            else:
                self._grid_region = self._configuration['grid_region']
                self._grid_resolution = self._configuration['grid_resolution']
                self._grid_type = self._configuration['grid_type']
                self._grid_shape = self._configuration['grid_shape']
                self._grid_i0 = self._configuration['i_ul']
                self._grid_j0 = self._configuration['j_ul']
                self._grid_iskip = self._configuration['i_skip']
                self._grid_jskip = self._configuration['j_skip']

            # Read information about the model run
            # Either from a file or directly
            if 'run_duration_filename' in self._configuration.keys():
                self._run_duration_filename = \
                    os.path.join(examples_directory,
                    self._configuration['run_duration_filename'])
                run_config = \
                    self.get_config_from_yaml_file(self._run_duration_filename)
                self._reference_date = run_config['model_reference_date']
                self._start_date = run_config['model_start_date']
                self._end_date = run_config['model_end_date']
                self._timestep_duration = run_config['model_timestep']
            else:
                self._reference_date = \
                        self._configuration['model_reference_date']
                self._start_date = self._configuration['model_start_date']
                self._end_date = self._configuration['model_end_date']
                self._timestep_duration = self._configuration['model_timestep']

            # Set up the grid subset
            if len(self._grid_shape) == 2:
                (self._grid_ydim, self._grid_xdim) = self._grid_shape
            else:
                raise ValueError("cannot handle grid of shape %s" %
                                str(self._grid_shape))
            self._grid_i1 = self._grid_i0 + self._grid_iskip * self._grid_xdim
            self._grid_j1 = self._grid_j0 + self._grid_jskip * self._grid_ydim

            # At a minimum, there must be temperature information
            if 'temperature_config_filename' in self._configuration.keys():
                self._temperature_config_filename = \
                    os.path.join(examples_directory,
                    self._configuration['temperature_config_filename'])
                assert_true(os.path.isfile(self._temperature_config_filename))
                temperature_configuration = \
                        self.get_config_from_yaml_file(self._temperature_config_filename)
                self._temperature_source_filename = \
                        os.path.join(data_directory,
                        temperature_configuration['temperature_source_filename'])
                self._temperature_first_date= \
                        temperature_configuration['dataset_first_date']
                self._temperature_last_date = \
                        temperature_configuration['dataset_last_date']
            else:
                self._temperature_source_filename = \
                        os.path.join(data_directory,
                        self._configuration['temperature_source_filename'])
                self._temperature_first_date= \
                        self._configuration['dataset_first_date']
                self._temperature_last_date = \
                        self._configuration['dataset_last_date']

            if self._configuration['precipitation_config_filename'] is not None:
                self._precipitation_config_filename = \
                        os.path.join(examples_directory,
                        self._configuration['precipitation_config_filename'])
                assert_true(os.path.isfile(\
                    self._precipitation_config_filename))
                self._calc_surface_fn = True
            else:
                self._calc_surface_fn = False
                self._calc_stefan_fn = False

            if self._calc_surface_fn and \
               self._configuration['soil_properties_config_filename'] \
                    is not None:
                self._soil_properties_config_filename = \
                        os.path.join(examples_directory,
                        self._configuration['soil_properties_config_filename'])
                assert_true(os.path.isfile(\
                    self._soil_properties_config_filename))
                self._calc_stefan_fn = True
            else:
                self._calc_stefan_fn = False

            self.initialize_input_vars_from_files()

        # If this is "WMT", then a flag must be set in the config file
        # to determine which
        elif self._configuration['input_var_source'] == 'WMT':
            self._using_WMT = True
            self._using_Files = False
            self._using_ConfigVals = False

            # Run duration
            self._reference_date = self._configuration['model_reference_date']
            self._start_date = self._configuration['model_start_date']
            self._end_date = self._configuration['model_end_date']
            self._timestep_duration = self._configuration['model_timestep']

            # Grid shape
            self._grid_type = self._configuration['grid_type']
            self._grid_shape = self._configuration['grid_shape']

        # If initialized completely from the config file, this is 'Default'
        elif self._configuration['input_var_source'] == 'Default':
            self._using_WMT = False
            self._using_Files = False
            self._using_ConfigVals = True

            # 'Default' has the same values as 'WMT',
            # but it also specifies a 3D (x, y, time) data array
            # from which the grids will derive their values

            # Run duration
            self._reference_date = self._configuration['model_reference_date']
            self._start_date = self._configuration['model_start_date']
            self._end_date = self._configuration['model_end_date']
            self._timestep_duration = self._configuration['model_timestep']

            # Grid shape
            self._grid_type = self._configuration['grid_type']
            self._grid_shape = self._configuration['grid_shape']
            if len(self._grid_shape) == 2:
                (self._grid_ydim, self._grid_xdim) = self._grid_shape
            else:
                raise ValueError("cannot handle grid of shape %s" %
                                str(self._grid_shape))

            # Set the valid dates for which temperature is available
            self._temperature_first_date = \
                self.datefrom(self._configuration['temperature_grid_date_0'])
            last_date_varname = 'temperature_grid_date_%d' % \
                (int(self._configuration['n_temperature_grid_fields']) - 1)
            self._temperature_last_date = \
                eval("self.datefrom(self._configuration['%s'])" % last_date_varname)

            # Parse the (x, y, time) data cube from which the values
            # will be drawn at the appropriate time
            self._temperature_dates, self._temperature_datacube = \
                self.initialize_datacube('temperature',
                                         self._configuration)
            if self._configuration['n_precipitation_grid_fields'] > 0:
                self._precipitation_dates, self._precipitation_datacube = \
                    self.initialize_datacube('precipitation',
                                             self._configuration)
            if self._configuration['n_soilproperties_grid_fields'] > 0:
                self._soilproperties_dates, self._soilproperties_datacube = \
                    self.initialize_datacube('soilproperties',
                                             self._configuration)

        # There are different ways of computing degree days.  Ensure that
        # the method specified has been coded
        self._dd_method = self._configuration['degree_days_method']
        if self._dd_method == 'MinJanMaxJul':
            self.initialize_MinJanMaxJul()
        self.ddf = np.zeros(self._grid_shape, dtype=np.float32)
        self.ddf.fill(np.nan)
        self.ddt = np.zeros(self._grid_shape, dtype=np.float32)
        self.ddt.fill(np.nan)

        # Initialize the output grids
        self.air_frost_number_Geo = \
                np.zeros(self._grid_shape, dtype=np.float32)
        if self._calc_surface_fn:
            self.surface_frost_number_Geo = \
                    np.zeros(self._grid_shape, dtype=np.float32)
        else:
            self.surface_frost_number_Geo = None
        if self._calc_stefan_fn:
            self.stefan_frost_number_Geo = \
                    np.zeros(self._grid_shape, dtype=np.float32)
        else:
            self.stefan_frost_number_Geo = None

        # Initialize the time information from a configuration file
        self.initialize_model_time()

        # Initialize the output filename
        self.initialize_output(self._configuration['output_directory'],
                               self._configuration['output_filename'])

        """
        # Note: initialization only allocates for values; it doesn't
        #       calculate the first timestep's values
        """

    def datefrom(self, datestring):
        if isinstance(datestring, str):
            return datetime.datetime.strptime(datestring, '%Y-%m-%d').date()
        else:
            return datestring


    def initialize_datacube(self, gridname, config):
        # Determine the number of lines for this grid
        exec("ngridlines = config['n_%s_grid_fields']" % gridname)

        # Create the datelist for this grid
        datelist = []
        for n in range(ngridlines):
            # This version converts a string to a datetime.date object
            # exec("thisdate = datetime.datetime.strptime(config['%s_grid_date_%d' % (gridname, n)], '%Y-%m-%d').date()")
            exec("thisdate = config['%s_grid_date_%d' % (gridname, n)]")
            if isinstance(thisdate, str):
                thisdate = self.datefrom(thisdate)
            datelist.append(thisdate)

        # Create the datacube for this grid
        slicelist = []
        for n in range(ngridlines):
            config_arg = eval(config['%s_grid_data_%d' % (gridname, n)])
            slicelist.append(np.array(config_arg))
        # In NumPy 1.10 later, we can use np.stack()
        # But for v1.09 -- which pip seems to install -- we need np.dstack()
        # datacube = np.stack(slicelist)
        datacube = self.np_stack(slicelist)

        return datelist, datacube

    def np_stack(self, arrays, axis=0):
        """
        In order to provide this function from numpy v1.10, this is copied
        from:
           https://github.com/numpy/numpy/blob/
           f4cc58c80df5202a743bddd514a3485d5e4ec5a4/numpy/core/shape_base.py
        Join a sequence of arrays along a new axis.
        The `axis` parameter specifies the index of the new axis in the dimensions
        of the result. For example, if ``axis=0`` it will be the first dimension
        and if ``axis=-1`` it will be the last dimension.
        .. versionadded:: 1.10.0
        Parameters
        ----------
        arrays : sequence of array_like
            Each array must have the same shape.
        axis : int, optional
            The axis in the result array along which the input arrays are stacked.
        Returns
        -------
        stacked : ndarray
            The stacked array has one more dimension than the input arrays.
        See Also
        --------
        concatenate : Join a sequence of arrays along an existing axis.
        split : Split array into a list of multiple sub-arrays of equal size.
        Examples
        Note: Replace '>>>' with '>>' so it doesn't doctest
        --------
        >> arrays = [np.random.randn(3, 4) for _ in range(10)]
        >> np.stack(arrays, axis=0).shape
        (10, 3, 4)
        >> np.stack(arrays, axis=1).shape
        (3, 10, 4)
        >> np.stack(arrays, axis=2).shape
        (3, 4, 10)
        >> a = np.array([1, 2, 3])
        >> b = np.array([2, 3, 4])
        >> np.stack((a, b))
        array([[1, 2, 3],
               [2, 3, 4]])
        >> np.stack((a, b), axis=-1)
        array([[1, 2],
               [2, 3],
               [3, 4]])
        """
        arrays = [np.asanyarray(arr) for arr in arrays]
        if not arrays:
            raise ValueError('need at least one array to stack')

        shapes = set(arr.shape for arr in arrays)
        if len(shapes) != 1:
            raise ValueError('all input arrays must have the same shape')

        result_ndim = arrays[0].ndim + 1
        if not -result_ndim <= axis < result_ndim:
            msg = 'axis {0} out of bounds [-{1}, {1})'.format(axis, result_ndim)
            raise IndexError(msg)
        if axis < 0:
            axis += result_ndim

        sl = (slice(None),) * axis + (np.newaxis,)
        expanded_arrays = [arr[sl] for arr in arrays]

        # Changed name from that of copied code
        # return _nx.concatenate(expanded_arrays, axis=axis)
        return np.concatenate(expanded_arrays, axis=axis)

    def get_datacube_slice(self, thisdate, datacube, cubedates):
        # Data are invalid until the first date in the datacube
        if thisdate < cubedates[0]:
            raise ValueError(
                "Date %s is before first valid date (%s) in datacube" %
                (thisdate, cubedates[0]))

        if thisdate > cubedates[-1]:
            raise ValueError(
                "Date %s is after last valid date (%s) in datacube" %
                (thisdate, cubedates[-1]))

        for date_okay, checkdate in enumerate(cubedates):
            if thisdate > checkdate:
                pass
            else:
                break

        date_okay -= 1
        return datacube[date_okay, :]


    def initialize_output(self, outdirname, outfilename):
        """ Initialize the output file

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
            outdirname = '.'
        elif outdirname == 'DATA':
            outdirname = data_directory
        elif outdirname == 'EXAMPLES':
            outdirname = examples_directory
        assert_true(os.path.isdir(outdirname))

        # Replace substitutions in filename
        outfilename = outfilename.replace('STARTDATE', str(self._start_date))
        outfilename = outfilename.replace('ENDDATE', str(self._end_date))
        outfilename = outfilename.replace('TIMESTEP',
                                          str(self._timestep_duration.days))

        self.output_filename = os.path.join(outdirname, outfilename)

        # Open the output netcdf file
        # Note: Be sure to close this, e.g. in self.finalize(), because
        # otherwise multiple openings, e.g. when testing, may fail
        self._output_fid = Dataset(self.output_filename, 'w', format='NETCDF4')

        # Add attributes
        #  including comments from config files
        setattr(self._output_fid, 'company',
                'Community Surface Dynamics Modeling System')
        setattr(self._output_fid, 'Permafrost Component', 'FrostnumberGeo')

        ### Init dimensions
        # Time dimension
        tdim = 0
        self._output_fid.createDimension("time", tdim)
        self._nc_time = self._output_fid.createVariable('time', 'i',
                                                        ('time',),zlib=True)
        setattr(self._nc_time, 'time_long_name', 'time')
        setattr(self._nc_time, 'time_standard_name', 'time')
        setattr(self._nc_time, 'time_units', 'months since 1900-01-01 00:00:00')
        self._nc_reference_time = datetime.date(1900, 1, 15)
        setattr(self._nc_time, 'time_format', 'modified julian day (MJD)')
        setattr(self._nc_time, 'time_time_zone', 'UTC')
        setattr(self._nc_time, 'time__FillValue', '-9999')
        self._nc_last_time_index = 0

        ### For now, the X and Y dimensions are just the indexes
        # X dimension
        self._output_fid.createDimension("x", self._grid_xdim)
        self._nc_x = self._output_fid.createVariable('x', 'f', ('x',),zlib=True)
        setattr(self._nc_x, 'x_long_name', 'projected x direction')
        setattr(self._nc_x, 'x_standard_name', 'x')
        setattr(self._nc_x, 'x_units', 'meters')
        setattr(self._nc_x, 'x__FillValue', 'NaN')
        # fill the x- values
        for x in range(self._grid_xdim):
            self._nc_x[x] = x

        # Y dimension
        self._output_fid.createDimension("y", self._grid_ydim)
        self._nc_y = self._output_fid.createVariable('y', 'f', ('y',),zlib=True)
        setattr(self._nc_y, 'y_long_name', 'projected y direction')
        setattr(self._nc_y, 'y_standard_name', 'y')
        setattr(self._nc_y, 'y_units', 'meters')
        setattr(self._nc_y, 'y__FillValue', 'NaN')
        # fill the y- values
        for y in range(self._grid_ydim):
            self._nc_y[y] = y

        ### Init grids with sizes
        # Allocate air frost number field
        self._nc_afn= \
                self._output_fid.createVariable('air_fn', 'f', ('time', 'y',
                                                                'x'),zlib=True)
        setattr(self._nc_afn, 'afn_long_name', 'Air Frost Number')
        setattr(self._nc_afn, 'afn_standard_name', 'Frostnumber_air')
        setattr(self._nc_afn, 'afn_units', 'none')
        setattr(self._nc_afn, 'afn__FillValue', '-99')

        # Allocate surface frost number field, if computing it
        if self._calc_surface_fn:
            self._nc_sfn= \
                    self._output_fid.createVariable('surface_fn', 'f', ('time',
                                                                        'y',
                                                                        'x'),zlib=True)
            setattr(self._nc_sfn, 'sfn_long_name', 'Surface Frost Number')
            setattr(self._nc_sfn, 'sfn_standard_name', 'Frostnumber_surface')
            setattr(self._nc_sfn, 'sfn_units', 'none')
            setattr(self._nc_sfn, 'sfn__FillValue', '-99')

        # Allocate Stefan frost number field, if computing it
        if self._calc_stefan_fn:
            self._nc_stfn= \
                    self._output_fid.createVariable(
                        'stefan_fn', 'f', ('time', 'y', 'x'),zlib=True)
            setattr(self._nc_stfn, 'stfn_long_name', 'Stefan Frost Number')
            setattr(self._nc_stfn, 'stfn_standard_name', 'Frostnumber_stefan')
            setattr(self._nc_stfn, 'stfn_units', 'none')
            setattr(self._nc_stfn, 'stfn__FillValue', '-99')


    def finalize(self):
        # Define this so we don't call the permamodel base class version
        self.finalize_frostnumber_Geo()

    def finalize_frostnumber_Geo(self):
        self._output_fid.close()

    def check_whether_output_timestep(self, this_timestep):
        # Only output on the 15th of each month
        this_date = self.get_date_from_timestep(this_timestep)
        if this_date.day == 15:
            return True
        else:
            return False

    def add_to_output(self):
        do_add = self.check_whether_output_timestep(self._timestep_current)
        if do_add:
            # I think this assumes that the file has the same time reference as
            # the model run...
            this_date = self.get_date_from_timestep(self._timestep_current)
            years = this_date.year - self._nc_reference_time.year
            months = this_date.month - self._nc_reference_time.month
            time_index = 12*years + months
            print("Adding output at index %d for date %s" % \
                  (time_index, str(this_date)))
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
        self.calculate_frost_numbers()
        self.add_to_output()

    def update(self):
        # Increment the model one time step
        self._date_current += self._timestep_duration
        self._timestep_current = \
                self.get_timestep_from_date(self._date_current)

        self.get_input_vars()
        self.compute_degree_days()
        self.calculate_frost_numbers()
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
        assert_greater_equal(self._timestep_last, self._timestep_first)

        # Set the current timestep to the first timestep
        self.set_current_date_and_timestep_with_timestep(self._timestep_first)

    def set_current_date_and_timestep_with_timestep(self, this_timestep):
        self._timestep_current = this_timestep
        self._date_current = \
                self.get_date_from_timestep(self._timestep_current)

    def get_date_from_timestep(self, timestep):
        return self._reference_date + timestep*self._timestep_duration

    def get_timestep_from_date(self, this_date):
        # If the timestep is an integer number of days, this would be fine:
        #return (this_date-self._reference_date).days /\
        #    self._timestep_duration.days
        return int(  (this_date-self._reference_date).total_seconds() / \
                     (self._timestep_duration.total_seconds()       ) +0.5 )

    def initialize_input_vars_from_files(self):
        # If the model does not have its input variables set by an
        # external runner, e.g. WMT, then the values will have to
        # be read from data files here.  This routine sets up the
        # files needed to do this.
        # Currently, this is done by reading values from a netcdf
        # file with the appropriate monthly values for
        self._temperature_dataset = \
                Dataset(self._temperature_source_filename, 'r', mmap=True)
        assert_true(self._temperature_dataset is not None)

    def get_input_vars(self):
        if self._using_WMT:
            # With WMT, the input variables will be set externally via BMI
            return
        elif self._using_Files or self._using_ConfigVals:
            # In standalone mode, variables must be set locally
            # All frost number types need temperature data
            if self._dd_method == 'MinJanMaxJul':
                # For the MinJanMaxJul method, need min and max temp fields
                (mindate, maxdate) = \
                        self.get_min_and_max_dates(self._date_current)
                self.T_air_min = self.get_temperature_field(mindate)
                self.T_air_max = self.get_temperature_field(maxdate)
        else:
            raise ValueError("Frostnumber must use either Files, ConfigVals \
                              or WMT to get input variables")

    def get_temperature_field(self, t_date = None):
        # By default, return the temperature field at the date of the current
        # timestep
        if t_date is None:
            t_date = self.get_date_from_timestep(self._timestep_current)

        if t_date >= self._temperature_first_date and \
           t_date <= self._temperature_last_date:

            #print("Getting temperature field for: %s" % str(t_date))

            t_index = self.get_temperature_month_index(t_date)

            if self._using_Files:
                # Files uses netcdf input
                temperature_subregion = \
                    self._temperature_dataset.variables['temp']\
                        [t_index]\
                        [self._grid_j0:self._grid_j1:self._grid_jskip,\
                        self._grid_i0:self._grid_i1:self._grid_iskip]
            elif self._using_ConfigVals:
                # ConfigVals is the Default, where the grids are in cfg file
                temperature_subregion = \
                        self.get_datacube_slice(t_date,
                                                self._temperature_datacube,
                                                self._temperature_dates).astype(np.float32)
            elif self._using_WMT:
                # Using WMT, this value will be set elsewhere
                pass

            assert_equal(temperature_subregion.shape, self._grid_shape)

            # Fill the field with Nan if value is 'missing value'
            nan_locations = temperature_subregion < -90
            temperature_subregion[nan_locations] = np.nan
        else:
            #print("Date is outside valid date range of temperature data")
            temperature_subregion = np.zeros(self._grid_shape, dtype=np.float32)
            temperature_subregion.fill(np.nan)

        return temperature_subregion

    def get_min_and_max_dates(self, this_date):
        """ Return the 15th of the month
        of the nearest preceding January and July """

        if this_date.month >= 7:
            maxyear = this_date.year
        else:
            maxyear = this_date.year - 1

        return (datetime.date(this_date.year, 1, 15),
                datetime.date(maxyear, 7, 15))

    def get_temperature_month_index(self, t_date):
        """ Returns the time index to the temperature dataset for this date """
        year_offset = t_date.year - self._temperature_first_date.year
        month_offset = t_date.month - self._temperature_first_date.month
        t_index = 12*year_offset + month_offset
        return t_index

    def initialize_MinJanMaxJul(self):
        # Initialize the temperature (and as needed, the precip and soil) grids
        self.T_air_min = np.zeros(self._grid_shape, dtype=np.float32)
        self.T_air_min.fill(np.nan)
        self.T_air_max = np.zeros(self._grid_shape, dtype=np.float32)
        self.T_air_max.fill(np.nan)

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

    def get_config_from_yaml_file(self, cfg_filename):
        raise RuntimeError("config from YAML not currently supported")
        cfg_struct = None
        try:
            with open(cfg_filename, 'r') as cfg_file:
                cfg_struct = yaml.load(cfg_file)
        except:
            print("\nError opening configuration file in\
                  get_config_from_yaml_file()")
            raise

        return cfg_struct

    def get_config_from_oldstyle_file(self, cfg_filename):
        cfg_struct = {}
        grid_struct = {}
        try:
            with open(cfg_filename, 'r') as cfg_file:
                # this was originally modeled after read_config_file()
                # in BMI_base.py as modified for cruAKtemp.py
                while True:
                    # Read lines from config file until no more remain
                    line = cfg_file.readline()
                    if line == "":
                        break

                    # Comments start with '#'
                    COMMENT = (line[0] == '#')

                    words = line.split('|')
                    if (len(words) ==4) and (not COMMENT):
                        var_name = words[0].strip()
                        value = words[1].strip()
                        var_type = words[2].strip()

                        # Process the variables based on variable name
                        if var_name[-4:] == 'date':
                            # date variables end with "_date"
                            cfg_struct[var_name] = \
                                datetime.datetime.strptime(
                                    value, "%Y-%m-%d").date()
                                #datetime.datetime.strptime(value, "%Y-%m-%d")
                        elif var_name[0:4] == 'grid':
                            # grid variables are processed after cfg file read
                            grid_struct[var_name] = value
                        elif var_name == 'timestep' \
                                or var_name == 'model_timestep':
                            # timestep is a timedelta object
                            cfg_struct[var_name] = \
                                datetime.timedelta(days=int(value))
                        elif var_type == 'int':
                            # Convert integers to int
                            cfg_struct[var_name] = int(value)
                        else:
                            # Everything else is just passed as a string
                            assert_equal(var_type, 'string')
                            cfg_struct[var_name] = value

        except:
            print("\nError opening configuration file in\
                  get_config_from_yaml_file()")
            raise

        # Process the grid information
        # I think I had rows and columns switched in cruAKtemp!
        #cfg_struct['grid_shape'] = (int(grid_struct['grid_columns']),
        #                            int(grid_struct['grid_rows']))
        cfg_struct['grid_shape'] = (int(grid_struct['grid_rows']),
                                    int(grid_struct['grid_columns']))
        cfg_struct['grid_type'] = grid_struct['grid_type']

        #for keyname in cfg_struct.keys():
        #    print(keyname)
        cfg_struct['grids'] = {'temperature': 'np.float'}
        if cfg_struct['n_precipitation_grid_fields'] > 0:
            cfg_struct['grids'] = {'precipitation': 'np.float'}
            self._calc_surface_fn = True
        else:
            self._calc_surface_fn = False
        if cfg_struct['n_soilproperties_grid_fields'] > 0:
            cfg_struct['grids'] = {'soilproperties': 'np.float'}
            self._calc_stefan_fn = True
        else:
            self._calc_stefan_fn = False

        return cfg_struct

    def calculate_frost_numbers(self):
        # Calculate all the frost numbers using the current data
        self.calculate_air_frost_number_Geo()
        self.calculate_surface_frost_number_Geo()
        self.calculate_stefan_frost_number_Geo()

        # Add these frost numbers to the output dictionary
        #self.update_output()

    def print_frost_numbers(self, year=-1):
        # if year is -1, then use the current year of self
        # otherwise, use the specified year
        if year > 0:
            print("Year: %d  F_air=%5.3f  F_surface=%5.3f  F_stefan=%5.3f" %
              (self.year, self.air_frost_number_Geo,
               self.surface_frost_number_Geo,
               self.stefan_frost_number_Geo))
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
        if self._dd_method == 'MinJanMaxJul':
            self.compute_degree_days_MinJanMaxJul()

    def compute_degree_days_MinJanMaxJul(self):

        # Input: T_hot (avg temp of warmest month)
        #        T_cold (avg temp of coldest month)

        # Output: ddf (degree freezing days)
        #         ddt (degree thawing days)

        ### Note: this doesn't work for more general cases
        ## Set values where NaN to NaN
        #nan_locations = np.isnan(self.T_air_min) | np.isnan(self.T_air_max)
        #self.ddf[nan_locations] = np.nan
        #self.ddt[nan_locations] = np.nan


        # Try looping...
        for (j, i), maxvalue in np.ndenumerate(self.T_air_max):
            minvalue = self.T_air_min[j, i]

            if minvalue==np.nan or maxvalue==np.nan:
                # Any non-numbers invalidate the dd calcs
                self.ddf[j, i] = np.nan
                self.ddt[j, i] = np.nan
            elif maxvalue < minvalue:
                # Can't have min temp > max temp!
                self.ddf[j, i] = np.nan
                self.ddt[j, i] = np.nan
            elif minvalue > 0.0:
                # Never freezes
                self.ddf[j, i] = 0.0
                self.ddt[j, i] = 365.0 * (minvalue+maxvalue)/2.0
            elif maxvalue <= 0.0:
                # Never thaws
                self.ddf[j, i] = -365.0 * (minvalue+maxvalue)/2.0
                self.ddt[j, i] = 0.0
            else:
                # Freezes in winter, thaws in summer
                T_average = (minvalue+maxvalue) / 2.0
                T_amplitude = (maxvalue-minvalue) / 2.0
                Beta = np.arccos(-T_average / T_amplitude)
                T_summer = T_average + T_amplitude * np.sin(Beta) / Beta
                T_winter = T_average - T_amplitude * np.sin(Beta) / (np.pi - Beta)
                L_summer = 365.0 * Beta / np.pi
                L_winter = 365.0 - L_summer
                self.ddt[j, i] = T_summer * L_summer
                self.ddf[j, i] = -T_winter * L_winter
    #   compute_degree_days()
    #-------------------------------------------------------------------
    def compute_air_frost_number_Geo(self):
        # Calculating Reduced Air Frost Number (pages 280-281).
        # The reduced frost number is close 0 for long summers and close to 1 for long winters.
        self.air_frost_number_Geo = np.sqrt(self.ddf) / ( np.sqrt( self.ddf) + np.sqrt( self.ddt) )

    #   update_air_frost_number_Geo()
    #-------------------------------------------------------------------

    def update_snow_prop(self):
        # find indexes for which temp > 0 and make precip = 0
        if (self.T_air_type != 'Scalar'): # if not should stop
            #wk = np.loadtxt('examples/prec.txt', skiprows=1,unpack=False)
            precipitation_filename = self.permafrost_dir +\
                "permamodel/examples/prec.txt"
            wk = np.loadtxt(precipitation_filename, skiprows=1,unpack=False)
            t_month = wk[:,0]
            prec_month = wk[:,1]

        pos_temp_ind=np.array(np.where(self.ta_month>0))
        prec_month[pos_temp_ind]=0
        neg_temp_ind=np.array(np.where(self.ta_month<=0))

        if not pos_temp_ind.any():
        # monthly temp is always below zero
        # i.e. it constantly snows over whole year
        # the point associated with glaciaer and needs to excluded
            print 'snows constatly: stop!'

        m=np.size(neg_temp_ind)
        pp=0.5; # assume only 50% of precip change to at the beg and end of the snow season

        # this is portions of the code assumes a perfect winter season
        # needs to be used with care when there is a warm month during snow season
        if (m==1):
            s_idx=neg_temp_ind[:,0]
            e_idx=neg_temp_ind[:,m-1]
            prec_month[s_idx]=prec_month[s_idx]*pp
        else:
            s_idx=neg_temp_ind[:,0]
            e_idx=neg_temp_ind[:,m-1]
            prec_month[s_idx]=prec_month[s_idx]*pp
            prec_month[e_idx]=prec_month[e_idx]*pp

        # sum up precip to get SWE
        j=0; s=0; swe=np.zeros(m);
        for i in range(s_idx,e_idx+1):
            s=s+prec_month[i]
            swe[j]=s
            j=j+1

        #calculating snow density, depth and thermal counductivity
        r_snow=np.zeros(m); # snow density in kg/m3
        h_snow=np.zeros(m); # snow depth in m
        c_snow=np.zeros(m); # snow depth in W/mK

        rho_sn_min=200; rho_sn_max=300 # allowed min and max snow density
        tauf=0.24 # e-folding value (see Verseghy, 1991)

        s=rho_sn_min
        s=((s - rho_sn_max)*np.exp(-tauf)) + rho_sn_max
        r_snow[0] = s
        for i in range(1,m):
        # starting from month 2 tauf should be multpled by the 30 days
        # otherwise snow thermal conductivity can be low and insulate ground well enough over the season
        # usually we assume constant max snow thermal conductivity over snow season
            s=((s - rho_sn_max)*np.exp(-tauf)) + rho_sn_max
            r_snow[i] = s

        h_snow  = (swe/(r_snow*0.001))
        # snow thermal conductivity according to M. Sturm, 1997.
        c_snow = (0.138-1.01*r_snow + 3.233*(r_snow**2))*1e-6

        self.r_snow=r_snow
        self.h_snow=h_snow
        self.c_snow=c_snow

    #   update_snow_prop()
    #-------------------------------------------------------------------
    def update_surface_frost_number_Geo(self):
        # phi [scalar]: sites latitude
        # Zs [scalar]: an average winter snow thickness
        # Zss [scalar]: a damping depth in snow
        # P [scalar]: length of an annual temperature cycle
        # k [scalar]: number of winter months
        # rho_s [scalar]: density of snow [kg m-3]
        # lambda_s [scalar]: snow thermal conductivity [W m-1 C-1]
        # c_s [scalar]: snow specific heat capacity [J kg-1 C-1]
        # alpha_s [scalar]: thermal diffusivity of snow
        # Uw [scalar]: mean winter wind speed [m s-1]
        # Aplus [scalar]: temperature amplitude at the surface with snow
        # Twplus [scalar]: the mean winter surface temperature
        # DDFplus [scalar]: freezing index at the surface
        # Tplus [scalar]: mean annual tempratures at the surface
        # Fplus [scalar] : surface frost number

        rho_s=np.mean(self.r_snow)
        lambda_s=np.mean(self.c_snow)
        Zs=np.mean(self.h_snow)
        # i am not sure what they mean by length of the annual temprature cycle
        # Something worthwhile discussing
        P=2*np.pi/365;

        c_s=7.79*self.Tw+2115                 #(eqn. 7)
        alpha_s=lambda_s/(c_s*rho_s)          #(eqn. 8)
        Zss=np.sqrt(alpha_s*P/np.pi)          #(eqn. 10)
        Aplus=self.A_air*np.exp(-Zs/Zss)      #(eqn. 9)
        Twplus=self.T_air-Aplus*np.sin(self.beta/(np.pi-self.beta)) #(eqn. 11)
        # Twplus is a mean winter surface temprature, I think, should be warmer than air temperature?
        # Here is another problem. DDFplus degree days winter should be positive.
        # The way it is written in the paper is wrong. I added a minus sign to fix it (see eqn. 2.9)
        DDFplus=-Twplus*self.Lw                                             #(eqn. 12)
        Tplus=(self.ddt-DDFplus)/365                                        #(eqn. 13)
        #Nevertheless the surface frost number is smaller than air which looks resonable to me.
        self.Fplus=np.sqrt(DDFplus)/(np.sqrt(self.ddt)+np.sqrt(DDFplus))    #(eqn. 14)
        self.Twplus=Twplus

    #   update_surface_frost_number()
    #-------------------------------------------------------------------
    def update_stefan_frost_number_Geo(self):
        # Zfplus [scalar]: the depth [m] to which forst extends
        # lambda_f [scalar]: frozen soil thermal conductivity [W m-1 C-1]
        # S [scalar]: is a const scalar factor [s d-1]
        # rho_d [scalar]: dry density of soil [kg m-3]
        # wf [scalar] : soil water content (proportion of dry weight)
        # L [scalar] : is a latent heat of fusion of water [J kg-1]

        lambda_f=1.67 # some dummy thermal conductivity
        # https://shop.bgs.ac.uk/GeoReports/examples/modules/C012.pdf
        sec_per_day=86400
        rho_d=2.798  # dry density of silt
        wf=0.4       # tipical for silty soils
        denominator=rho_d*wf*self.Lf
        self.Zfplus=np.sqrt(2*lambda_f*sec_per_day*np.abs(self.Twplus)*self.Lw/denominator)               #(eqn. 15)
        print 'Zfplus=',self.Zfplus

        # assuming 3 soil layers with thickness 0.25, 0.5 and 1.75 m
        # and thermal conductivities 0.08, 1.5, and 2.1
        soil_thick=np.array([0.25, 0.5 , 1.75])
        lambda_b=np.array([0.08, 1.5, 2.1])
        # resistivity R
        R=soil_thick/lambda_b
        QL=self.Lf/1000 # volumetric latent heat, 1000 s a density of water
        #partial freezing thawing index DD
        DD=np.zeros(3)
        Z=np.zeros(3)
        DD[0]=0.5*QL*soil_thick[0]*R[0]/sec_per_day
        S=0;
        for i in range(1,3):
            S=R[i]+S
            DD[i]=(S+0.5*R[i-1])*QL*soil_thick[i]/sec_per_day
            #The depth of the frost thaw penetration
        S=0; Z_tot=0
        for i in range(0,3):
            #The depth of the frost thaw penetration
            Z[i]=np.sqrt(2*lambda_b[i]*sec_per_day*DD[i]/QL + lambda_b[i]**2*S**2) \
                - lambda_b[i]*S
            S=R[i]+S
            Z_tot=Z_tot + Z[i]

        self.Z_tot=Z_tot
        self.stefan_number = np.sqrt(self.Fplus) / ( np.sqrt( self.Fplus) + np.sqrt( self.Z_tot) )

    #   update_stefan_frost_number()
    #-------------------------------------------------------------------

if __name__ == "__main__":
    # Run the FrostnumberGeo model
    # Currently, this just runs the defaults
    fn_geo = FrostnumberGeoMethod()
    fn_geo.initialize_frostnumberGeo_component()
    fn_geo.initial_update()
    fn_geo.update_until_timestep(fn_geo._timestep_last)
    fn_geo.finalize()

