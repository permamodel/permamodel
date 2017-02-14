"""
test_frost_number_2D.py
  tests of the 2D version of the frost_number component of permamodel
"""

from permamodel.components import frost_number_2D
import os
import numpy as np
import datetime
from .. import permamodel_directory, examples_directory
from nose.tools import (assert_is_instance, assert_greater_equal,
                        assert_less_equal, assert_almost_equal,
                        assert_greater, assert_less, assert_in,
                        assert_false, assert_true, assert_equal)

# List of files to be removed after testing is complete
# use files_to_remove.append(<filename>) to add to it
files_to_remove = []

def setup_module():
    """ Standard fixture called before any tests in this file are performed """
    pass

def teardown_module():
    """ Standard fixture called after all tests in this file are performed """
    for f in files_to_remove:
        if os.path.exists(f):
            os.remove(f)

# ---------------------------------------------------
# Tests that ensure we are reaching this testing file
# ---------------------------------------------------
def test_testing_2D():
    # This should pass as long as this routine is getting called
    assert(True)

def test_can_initialize_2D_frostnumber_module():
    fn2 = frost_number_2D.Frostnumber2DMethod()

def test_2D_frostnumber_has_default_config_file():
    fn2 = frost_number_2D.Frostnumber2DMethod()
    assert_true(fn2._config_filename is not None)

def test_2D_frostnumber_can_be_passed_config_filename():
    fn2 = frost_number_2D.Frostnumber2DMethod(cfgfile="a file")
    assert_true(fn2._config_filename == "a file")

def test_2D_frostnumber_initializes_from_default_config_file():
    fn2 = frost_number_2D.Frostnumber2DMethod()
    assert_true(os.path.isfile(fn2._config_filename))
    fn2.initialize_frostnumber2D_component()
    assert_true(fn2._grid_type == 'uniform rectilinear')
    assert_true(fn2._calc_surface_fn is not None)
    assert_true(fn2._calc_stefan_fn is not None)
    assert_in(fn2._dd_method, ('ObservedMinMax', 'MinJanMaxJul',
            'MonthlyAverages', 'DailyValues'))
    if fn2._dd_method == 'MinJanMaxJul':
        assert_equal(fn2.T_air_min.shape, fn2._grid_shape)
        assert_equal(fn2.T_air_max.shape, fn2._grid_shape)
        if fn2._using_Files:
            assert_true(fn2._temperature_dataset is not None)

    fn2.finalize_frostnumber_2D()
    files_to_remove.append(fn2.output_filename)

def test_2D_frostnumber_can_compute_real_date_from_timestep():
    fn2 = frost_number_2D.Frostnumber2DMethod()
    fn2.initialize_frostnumber2D_component()

    # Timestep should be one day
    assert_equal(fn2._timestep_duration, datetime.timedelta(days=1))

    # Model reference date should have timestep of zero
    assert_equal(fn2.get_timestep_from_date(fn2._reference_date),
                 datetime.timedelta(days=0).days)

    # Timestep of first date should be first timestep
    assert_equal(fn2._timestep_first,
                 fn2.get_timestep_from_date(fn2._start_date))

    # Last date should be date of last timestep
    assert_equal(fn2._end_date,
                 fn2.get_date_from_timestep(fn2._timestep_last))

    fn2.finalize_frostnumber_2D()

def test_2D_frostnumber_can_return_temperature_field():
    fn2 = frost_number_2D.Frostnumber2DMethod()
    fn2.initialize_frostnumber2D_component()

    current_temperature_field = fn2.get_temperature_field()
    assert_equal(current_temperature_field.shape, fn2._grid_shape)

    bad_date_field = fn2.get_temperature_field(datetime.date(100, 1, 1))
    assert_true(np.all(np.isnan(bad_date_field)))

    fn2.finalize_frostnumber_2D()

def test_2D_frostnumber_get_latest_min_max_months():
    fn2 = frost_number_2D.Frostnumber2DMethod()
    fn2.initialize_frostnumber2D_component()

    this_date = datetime.date(1905, 1, 1)
    (mindate, maxdate) = fn2.get_min_and_max_dates(this_date)
    assert_equal(mindate, datetime.date(1905, 1, 15))
    assert_equal(maxdate, datetime.date(1904, 7, 15))

    this_date = datetime.date(2004, 8, 1)
    (mindate, maxdate) = fn2.get_min_and_max_dates(this_date)
    assert_equal(mindate, datetime.date(2004, 1, 15))
    assert_equal(maxdate, datetime.date(2004, 7, 15))

    fn2.finalize_frostnumber_2D()

def test_2D_frostnumber_compute_array_of_degree_days():
    fn2 = frost_number_2D.Frostnumber2DMethod()
    fn2.initialize_frostnumber2D_component()

    if fn2._using_WMT:
        print("Can't test WMT in standalone mode")
    elif not fn2._using_Files:
        raise ValueError("Frostnumber must use either Files \
                         or WMT to get input variables")

    # After this, we assume we are using Files for input data
    assert_true(fn2._using_Files)

    # Test that we get NaN-filled arrays when no input vars available
    fn2.set_current_date_and_timestep_with_timestep(-100)
    fn2.get_input_vars()
    assert_less(fn2._date_current, fn2._temperature_first_date)
    if fn2._dd_method == 'MinJanMaxJul':
        # A timestep out of bounds should yield nans for temp and dd
        assert_true(np.all(np.isnan(fn2.T_air_min)))
        assert_true(np.all(np.isnan(fn2.T_air_max)))

        # Calculating degree days on all NaNs yields all NaNs
        fn2.compute_degree_days()
        assert_true(np.all(np.isnan(fn2.ddt)))
        assert_true(np.all(np.isnan(fn2.ddf)))

        # Calculating air frost number on NaNs yields NaNs
        fn2.compute_air_frost_number_2D()
        assert_true(np.all(np.isnan(fn2.air_frost_number_2D)))

    # Test that we get real values
    fn2.set_current_date_and_timestep_with_timestep(1000)
    fn2.get_input_vars()
    assert_greater_equal(fn2._date_current, fn2._temperature_first_date)
    if fn2._dd_method == 'MinJanMaxJul':
        # The default should have no missing temperature values
        assert_false(np.any(np.isnan(fn2.T_air_min)))
        assert_false(np.any(np.isnan(fn2.T_air_max)))
        #fn2.T_air_min.tofile("Tairmin.dat")
        #fn2.T_air_max.tofile("Tairmax.dat")

        # Calculating degree days on all NaNs yields all NaNs
        fn2.compute_degree_days()
        assert_false(np.any(np.isnan(fn2.ddt)))
        assert_false(np.any(np.isnan(fn2.ddf)))
        #fn2.ddt.tofile("ddt.dat")
        #fn2.ddf.tofile("ddf.dat")

        # Calculating air frost number on NaNs yields NaNs
        fn2.compute_air_frost_number_2D()
        assert_false(np.any(np.isnan(fn2.air_frost_number_2D)))
        #fn2.air_frost_number_2D.tofile("afn2d.dat")
    fn2.finalize_frostnumber_2D()

def test_2D_frostnumber_output_a_netcdf_file():
    fn2 = frost_number_2D.Frostnumber2DMethod()
    fn2.initialize_frostnumber2D_component()
    assert_true(fn2.output_filename is not None)
    assert_true(fn2.output_filename[-3:] == '.nc')
    fn2.initial_update()
    for t in range(0, 10):
        fn2.get_input_vars()
        fn2.compute_degree_days()
        fn2.calculate_frost_numbers()
        fn2.add_to_output()
        fn2.update()
    fn2.finalize()

# This test is long, so doesn't need to run too often
# It is the default operation of calling the routine without arguments
#def test_2D_frostnumber_runs_the_test_config_file():
#    fn2 = frost_number_2D.Frostnumber2DMethod()
#    fn2.initialize_frostnumber2D_component()
#    fn2.initial_update()
#    fn2.update_until_timestep(fn2._timestep_last)
#    fn2.finalize()

