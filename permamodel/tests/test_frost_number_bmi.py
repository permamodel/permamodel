"""
test_frost_number_bmi.py
  tests of the frost_number component of permamodel using bmi API
"""

from permamodel.components import frost_number
from permamodel.components import perma_base
import os
import numpy as np

# ---------------------------------------------------
# Tests that ensure we have bmi functionality
# ---------------------------------------------------
def test_frost_number_has_initialize():
    # Can we call an initialize function?
    fn = frost_number.frostnumber_method()
    # With hard-coded cfg filename
    #fn.initialize(cfg_file='/home/scotts/permamodel/permamodel/examples/Fairbanks_frostnumber_method.cfg', SILENT=True)
    # With relative cfg filename
    fn.initialize(cfg_file='./permamodel/examples/Frostnumber_example_singlesite_singleyear.cfg', SILENT=True)

def test_frost_number_initialize_sets_year():
    fn = frost_number.frostnumber_method()
    fn.initialize(cfg_file='./permamodel/examples/Frostnumber_example_singlesite_singleyear.cfg', SILENT=True)

    # Assert the values from the cfg file
    assert(fn.year == 2000)

def test_frost_number_initialize_sets_air_min_and_max():
    fn = frost_number.frostnumber_method()
    fn.initialize(cfg_file='./permamodel/examples/Frostnumber_example_singlesite_singleyear.cfg', SILENT=True)

    # Assert the values from the cfg file
    assert(fn.T_air_min == -20.0)
    assert(fn.T_air_max == 10.0)

def test_frost_number_update_increments_year():
    fn = frost_number.frostnumber_method()
    fn.initialize(cfg_file='./permamodel/examples/Frostnumber_example_singlesite_multiyear.cfg', SILENT=True)

    fn.update(dt=fn.dt)
    assert(fn.year == fn.start_year + fn.dt)
    assert(fn.year != fn.start_year)

def test_frost_number_update_changes_air_frost_number():
    fn = frost_number.frostnumber_method()
    fn.initialize(cfg_file='./permamodel/examples/Frostnumber_example_singlesite_multiyear.cfg', SILENT=True)

    afn0 = fn.air_frost_number
    fn.update(dt=fn.dt)
    afn1 = fn.air_frost_number
    assert(afn0 != afn1)

def test_frost_number_runs_several_years():
    fn = frost_number.frostnumber_method()
    fn.initialize(cfg_file='./permamodel/examples/Frostnumber_example_singlesite_multiyear.cfg', SILENT=True)

    while fn.year < fn.end_year:
        fn.update(dt=fn.dt)

    assert(fn.output is not None)

    # Ensure that each year exists in the output dictionary
    year = fn.start_year
    while year < fn.end_year:
        assert(year in fn.output.keys())
        year += 1

    # print the output to check it
    #fn.print_frost_numbers()



""" ------------DELETED CODE---------------------------
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
# Test that environment variables have been set
# ---------------------------------------------------
def test_environment_variables_set():
    env_var_to_test = "PERMAMODEL_EXAMPLEDIR"
    if not os.environ.get(env_var_to_test):
        raise ValueError('Environment variable %s not set', env_var_to_test)

    env_var_to_test = "PERMAMODEL_DATADIR"
    if not os.environ.get(env_var_to_test):
        raise ValueError('Environment variable %s not set', env_var_to_test)

# ---------------------------------------------------
# Tests that input data is being read correctly
# ---------------------------------------------------
def test_can_initialize_frostnumber_method_from_file():
    fn = frost_number.frostnumber_method()
    cfg_file = os.path.join(os.environ.get('PERMAMODEL_EXAMPLEDIR',\
                                '/examples/'),
                 'Fairbanks_frostnumber_method.cfg')
    fn.initialize(cfg_file=cfg_file, SILENT=True)

def test_frostnumber_method_has_date_and_location():
    fn = frost_number.frostnumber_method()
    cfg_file = os.path.join(os.environ.get('PERMAMODEL_EXAMPLEDIR',\
                                '/examples/'),
                 'Fairbanks_frostnumber_method.cfg')
    fn.initialize(cfg_file=cfg_file, SILENT=True)
    assert(fn.year >= 0)
    assert(fn.year == fn.start_year)

def test_can_get_temperature_filename_from_date_and_location():
    fn = frost_number.frostnumber_method()
    cfg_file = os.path.join(os.environ.get('PERMAMODEL_EXAMPLEDIR',\
                                '/examples/'),
                 'Fairbanks_frostnumber_method.cfg')
    fn.initialize(cfg_file=cfg_file, SILENT=True)
    cfg_file = os.path.join(os.environ.get('PERMAMODEL_EXAMPLEDIR',\
                                '/examples/'),
                 'Fairbanks_frostnumber_method.cfg')
    fn.initialize(cfg_file=cfg_file, SILENT=True)
    # Given year and month, calculate the tiff_filename
    fname = fn.get_temperature_tiff_filename(fn.year, 6)
    assert(fname != None)

def test_get_temperature_from_cru_indexes():
    fn = frost_number.frostnumber_method()
    cfg_file = os.path.join(os.environ.get('PERMAMODEL_EXAMPLEDIR',\
                                '/examples/'),
                 'Fairbanks_frostnumber_method.cfg')
    fn.initialize(cfg_file=cfg_file, SILENT=True)

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
    cfg_file = os.path.join(os.environ.get('PERMAMODEL_EXAMPLEDIR',\
                                '/examples/'),
                 'Fairbanks_frostnumber_method.cfg')
    fn.initialize(cfg_file=cfg_file, SILENT=True)

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

    # Test lower left corner
    lon = 175.49893
    lat = 49.106936
    (i, j) = fn.get_cru_indexes_from_lon_lat(lon, lat, month, year)
    assert(i==0)
    assert(j==2556)


def test_get_lon_lat_from_cru_indexes():
    fn = frost_number.frostnumber_method()
    cfg_file = os.path.join(os.environ.get('PERMAMODEL_EXAMPLEDIR',\
                                '/examples/'),
                 'Fairbanks_frostnumber_method.cfg')
    fn.initialize(cfg_file=cfg_file, SILENT=True)

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
    cfg_file = os.path.join(os.environ.get('PERMAMODEL_EXAMPLEDIR',\
                                '/examples/'),
                 'Fairbanks_frostnumber_method.cfg')
    fn.initialize(cfg_file=cfg_file, SILENT=True)
    month = 6
    year = 2009

    # I can see the rightmost value on the geotiff at (4761, 1973)
    i = 4761
    j = 1973
    (lon, lat) = fn.get_lon_lat_from_cru_indexes(i, j, month, year)
    temperature = fn.get_temperature_from_cru(lon, lat, month, year)
    assert(abs(temperature - 8.6) < 1e-3)

"""
