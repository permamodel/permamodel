"""
test_frost_number.py
  tests of the frost_number component of permamodel
"""

from permamodel.components import frost_number
import os
import numpy as np

# ---------------------------------------------------
# Tests that ensure we are reaching this testing file
# ---------------------------------------------------
def test_testing():
    # This should pass as long as this routine is getting called
    assert(True)

# ---------------------------------------------------
# Tests that the frost_number module is importing
# ---------------------------------------------------
def test_can_initialize_frost_number_module():
    fn = frost_number.frostnumber_method
    assert(True)

def test_have_output_var_names():
    fn = frost_number.frostnumber_method
    assert(fn._output_var_names != None)


# ---------------------------------------------------
# Tests that input data is being read correctly
# ---------------------------------------------------
def test_can_initialize_frostnumber_method_from_file():
    fn = frost_number.frostnumber_method()
    fn.initialize(cfg_file="/home/scotts/permamodel/permamodel/examples/Fairbanks_frostnumber_method.cfg", SILENT=True)

def test_frostnumber_method_has_date_and_location():
    fn = frost_number.frostnumber_method()
    fn.initialize(cfg_file="/home/scotts/permamodel/permamodel/examples/Fairbanks_frostnumber_method.cfg", SILENT=True)
    assert(fn.year >= 0)
    assert(fn.year == fn.start_year)

def test_can_get_temperature_filename_from_date_and_location():
    fn = frost_number.frostnumber_method()
    fn.initialize(cfg_file="/home/scotts/permamodel/permamodel/examples/Fairbanks_frostnumber_method.cfg",
                 SILENT=True)
    # Given year and month, calculate the tiff_filename
    fname = fn.get_temperature_tiff_filename(fn.year, 6)
    assert(fname != None)

def test_get_temperature_from_cru_indexes():
    fn = frost_number.frostnumber_method()
    fn.initialize(cfg_file="/home/scotts/permamodel/permamodel/examples/Fairbanks_frostnumber_method.cfg",
                 SILENT=True)

    # Test using data file for June of 2009
    month = 6
    year = 2009

    # We know that the temperature at x=2640 y=924 is 15.4 
    i = 2640
    j = 924
    temp = fn.get_temperature_from_cru_indexes(i, j, month, year)
    np.testing.assert_almost_equal(temp, 15.4, decimal=3)

    # We know that the temperature at x=2640 y=924 is 12.8
    i = 2612
    j = 1328
    temp = fn.get_temperature_from_cru_indexes(i, j, month, year)
    np.testing.assert_almost_equal(temp, 12.8, decimal=3)

def test_get_cru_indexes_from_lon_lat():
    fn = frost_number.frostnumber_method()
    fn.initialize(cfg_file="/home/scotts/permamodel/permamodel/examples/Fairbanks_frostnumber_method.cfg",SILENT=True)

    # Test using CRU tiff file for June 2009
    month = 6
    year = 2009

    # Test midpoint
    # Difficult to test because even number of x points
    lon = -160.55371
    lat = 62.382381
    (i, j) = fn.get_cru_indexes_from_lon_lat(lon, lat, month, year)
    assert(j==1278)

    # Test upper left corner
    lon = 157.44516
    lat = 63.90865
    (i, j) = fn.get_cru_indexes_from_lon_lat(lon, lat, month, year)
    assert(i==0)
    assert(j==0)

    # Test lower right corner
    lon = -132.17774
    lat = 51.460417
    (i, j) = fn.get_cru_indexes_from_lon_lat(lon, lat, month, year)
    assert(i==4761)
    assert(j==2556)

def test_get_lon_lat_from_cru_indexes():
    fn = frost_number.frostnumber_method()
    fn.initialize(cfg_file="/home/scotts/permamodel/permamodel/examples/Fairbanks_frostnumber_method.cfg",SILENT=True)

    # Try to get the upper left corner indexes
    month = 6
    year = 2009
    i = -0.5
    j = -0.5
    (lon, lat) = fn.get_lon_lat_from_cru_indexes(i, j, month, year)
    assert(abs(157.44516 - lon) < 1e-3)
    assert(abs(63.90865 - lat) < 1e-3)

    # Try to get the lower right corner indexes
    month = 6
    year = 2009
    i = 4761.5
    j = 2556.5
    (lon, lat) = fn.get_lon_lat_from_cru_indexes(i, j, month, year)
    assert(abs(-132.17774 - lon) < 1e-3)
    assert(abs(51.460417 - lat) < 1e-3)


def test_get_temperature_from_cru():
    fn = frost_number.frostnumber_method()
    fn.initialize(cfg_file="/home/scotts/permamodel/permamodel/examples/Fairbanks_frostnumber_method.cfg",SILENT=True)
    month = 6
    year = 2009
    lon = -160.55371
    lat = 62.382381
    temperature = fn.get_temperature_from_cru(lon, lat, month, year)
    print("Temperature: %f" % temperature)

