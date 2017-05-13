"""
test_frost_number_bmi.py
  tests of the frost_number component of permamodel using bmi API
"""

from permamodel.components import frost_number_Geo
from permamodel.components import bmi_frost_number_Geo
from permamodel.components import perma_base
from dateutil.relativedelta import relativedelta
import os
import numpy as np
from .. import examples_directory
from nose.tools import (assert_is_instance, assert_greater_equal,
                        assert_less_equal, assert_almost_equal,
                        assert_greater, assert_in, assert_true,
                        assert_false, assert_equal, assert_raises)
from numpy.testing import assert_array_equal
import datetime


# Set the file names for the example cfg files
config_filename = \
        os.path.join(examples_directory,
                     'FrostnumberGeo_Default.cfg')

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
# Tests that ensure we have bmi functionality
# Note: the netcdf functions seem to lead to errors if the netcdf
# files are not closed.  This is done in the finalize() routine
# which should therefore be called after every test
# ---------------------------------------------------
def test_frost_number_Geo_has_initialize():
    # Can we call an initialize function?
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize(cfg_file=config_filename)
    fng.finalize()  # Must have this or get IOError later

def test_frost_number_initialize_sets_year():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize(cfg_file=config_filename)

    # Assert the values from the cfg file
    assert_equal(fng._model._date_current.year, 1901)
    fng.finalize()  # Must have this or get IOError later

def test_frost_number_initialize_sets_air_min_and_max():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize(cfg_file=config_filename)

    # The temperature arrays are initialized to all NaNs
    nan_array = np.zeros((3, 2), dtype=np.float32)
    nan_array.fill(np.nan)
    assert_array_equal(fng._model.T_air_min, nan_array)
    assert_array_equal(fng._model.T_air_max, nan_array)
    fng.finalize()  # Must have this or get IOError later

def test_frost_number_update_increments_time():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()
    fng.update()

    assert_equal(fng._model._date_current,
                 fng._model._start_date + \
                 relativedelta(years=fng._model._timestep_duration))
    fng.finalize()  # Must have this or get IOError later

def test_frost_number_update_changes_air_frost_number():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()

    fng.update()
    afn0 = fng._model.air_frost_number_Geo
    fng.update_until(fng._model._date_current + datetime.timedelta(days=370))
    afn1 = fng._model.air_frost_number_Geo
    try:
        # If this is none, then the first array was all NaNs
        assert_true(assert_array_equal(afn0, afn1) is None)
    except AssertionError:
        # If this is raised, then the arrays are different, which is nood
        pass
    fng.finalize()  # Must have this or get IOError later

def test_frost_number_get_current_time_returns_scalar():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()
    current_time = fng.get_current_time()
    assert_true(isinstance(current_time, float) \
                or isinstance(current_time, int))
    fng.finalize()  # Must have this or get IOError later

def test_frost_number_get_end_time_returns_scalar():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()
    end_time = fng.get_end_time()
    assert(isinstance(end_time, float) \
           or isinstance(end_time, int))
    fng.finalize()  # Must have this or get IOError later

def test_frost_number_implements_update_until():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()

    fng.update_until(fng._model._end_date)
    assert_true(fng._model._date_current, fng._model._end_date)

    fng.finalize()  # Must have this or get IOError later

def test_FNGeo_computes_default_values():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()
    for n in range(5):
        if (fng._model.T_air_min[0, 0] + fng._model.T_air_max[0, 0]) == 0:
            assert_equal(fng._values['frostnumber__air'][0, 0], 0.5)
        if (fng._model.T_air_min[0, 0] <= 0.0) and \
           (fng._model.T_air_max[0, 0] <= 0.0):
            assert_equal(fng._values['frostnumber__air'][0, 0], 1.0)
        if (fng._model.T_air_min[0, 0] > 0.0) and \
           (fng._model.T_air_max[0, 0] > 0.0):
            assert_equal(fng._values['frostnumber__air'][0, 0], 0.0)
        #print("FNGeo frostnumber__air: %d" % n)
        #print(fng._model.T_air_min)
        #print(fng._model.T_air_max)
        #print(fng._values['frostnumber__air'])
        fng.update()
    fng.finalize()  # Must have this or get IOError later

def test_FNGeo_get_attribute():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()
    # Check an attribute that exists
    this_att = fng.get_attribute('time_units')
    assert_equal(this_att, 'years')

    # Check an attribute that doesn't exist
    assert_raises(KeyError, fng.get_attribute, 'not_an_attribute')

    fng.finalize()  # Must have this or get IOError later

def test_FNGeo_get_input_var_names():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()
    assert_in('atmosphere_bottom_air__temperature', fng.get_input_var_names())

    fng.finalize()  # Must have this or get IOError later

def test_FNGeo_get_output_var_names():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()
    assert_in('frostnumber__air', fng.get_output_var_names())

    fng.finalize()  # Must have this or get IOError later

def test_FNGeo_get_set_value():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()

    airtempref = fng.get_value_ref('atmosphere_bottom_air__temperature')
    airtempval = fng.get_value('atmosphere_bottom_air__temperature')
    airtempnew = 123 * np.ones_like(airtempval)
    fng.set_value('atmosphere_bottom_air__temperature', airtempnew)
    assert_raises(AssertionError, assert_array_equal, airtempval, airtempnew)

    fng.finalize()  # Must have this or get IOError later

def test_FNGeo_get_grid_shape():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()

    airtemp_gridid = fng.get_var_grid('atmosphere_bottom_air__temperature')
    assert_equal((3, 2), fng.get_grid_shape(airtemp_gridid))

    fng.finalize()  # Must have this or get IOError later

def test_FNGeo_get_grid_size():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()

    airtemp_gridid = fng.get_var_grid('atmosphere_bottom_air__temperature')
    assert_equal(6, fng.get_grid_size(airtemp_gridid))

    fng.finalize()  # Must have this or get IOError later

def test_FNGeo_get_grid_spacing():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()

    airtemp_gridid = fng.get_var_grid('atmosphere_bottom_air__temperature')
    assert_array_equal(
        np.array([1.0, 1.0]),
        fng.get_grid_spacing(airtemp_gridid))

    fng.finalize()  # Must have this or get IOError later

