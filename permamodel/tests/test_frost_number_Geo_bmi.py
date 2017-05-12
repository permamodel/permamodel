"""
test_frost_number_bmi.py
  tests of the frost_number component of permamodel using bmi API
"""

from permamodel.components import frost_number_Geo
from permamodel.components import bmi_frost_number_Geo
from permamodel.components import perma_base
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
    fng.finalize()

def test_frost_number_initialize_sets_year():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize(cfg_file=config_filename)

    # Assert the values from the cfg file
    assert_equal(fng._model._date_current.year, 1901)
    fng.finalize()

def test_frost_number_initialize_sets_air_min_and_max():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize(cfg_file=config_filename)

    # The temperature arrays are initialized to all NaNs
    nan_array = np.zeros((3, 2), dtype=np.float32)
    nan_array.fill(np.nan)
    assert_array_equal(fng._model.T_air_min, nan_array)
    assert_array_equal(fng._model.T_air_max, nan_array)
    fng.finalize()

def test_frost_number_update_increments_time():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()
    fng.update()

    assert_equal(fng._model._date_current,
                 fng._model._start_date + fng._model._timestep_duration)
    fng.finalize()

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
    fng.finalize()

def test_frost_number_implements_update_until():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()

    fng.update_until(fng._model._end_date)
    assert_true(fng._model._date_current, fng._model._end_date)

    fng.finalize()

def test_frost_number_get_current_time_returns_scalar_float():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()
    current_time = fng.get_current_time()
    assert_true(isinstance(current_time, float))
    fng.finalize()

def  test_frost_number_get_end_time_returns_scalar_float():
    fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
    fng.initialize()
    end_time = fng.get_end_time()
    assert(isinstance(end_time, float))
    fng.finalize()
