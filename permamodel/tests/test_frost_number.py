"""
test_frost_number.py
  tests of the frost_number component of permamodel
"""
from __future__ import print_function

import os

from pytest import approx

from permamodel.components import frost_number

from .. import examples_directory

# List of files to be removed after testing is complete
# use files_to_remove.append(<filename>) to add to it
files_to_remove = []


def setup_module():
    """Standard fixture called before any tests in this file are performed"""
    pass


def teardown_module():
    """Standard fixture called after all tests in this file are performed"""
    for f in files_to_remove:
        if os.path.exists(f):
            os.remove(f)


# ---------------------------------------------------
# Tests that input data is being read correctly
# ---------------------------------------------------
def test_initialize_frostnumber_method():
    """Test fn initialization without config file"""
    fn = frost_number.FrostnumberMethod()
    fn.initialize()


def test_end_year_before_start_year_error():
    """Test end_year before start_year"""
    fn = frost_number.FrostnumberMethod()
    fn.initialize()
    fn.start_year = 2000
    fn.end_year = 1990
    fn.initialize_frostnumber_component()
    assert fn.start_year == fn.end_year


def test_can_initialize_frostnumber_method_from_scalar_file():
    """Test fn initialization from config file (scalar)"""
    fn = frost_number.FrostnumberMethod()
    cfg_file = os.path.join(examples_directory, "Frostnumber_example_scalar.cfg")
    fn.initialize(cfg_file=cfg_file)
    files_to_remove.append(fn.fn_out_filename)


def test_frostnumber_method_has_date_info():
    """Test that fn has time values"""
    fn = frost_number.FrostnumberMethod()
    cfg_file = os.path.join(examples_directory, "Frostnumber_example_scalar.cfg")
    fn.initialize(cfg_file=cfg_file)
    assert fn.year >= 0
    assert fn.year == fn.start_year


def test_can_initialize_frostnumber_method_from_timeseries_file():
    """Test fn initialization from config file (time series)"""
    fn = frost_number.FrostnumberMethod()
    cfg_file = os.path.join(examples_directory, "Frostnumber_example_timeseries.cfg")
    fn.initialize(cfg_file=cfg_file)
    files_to_remove.append(fn.fn_out_filename)


def test_frostnumber_method_calculates_fn():
    """Test fn gets calculated"""
    fn = frost_number.FrostnumberMethod()
    fn.initialize()
    assert fn.air_frost_number == approx(0.6326749410343562)


def test_frostnumber_method_calculates_exact_fn():
    """Test fn gets calculated"""
    fn = frost_number.FrostnumberMethod()
    fn.initialize()
    fn.T_air_min = [5.0]
    fn.T_air_max = [15.0]
    fn.calculate_air_frost_number()
    assert fn.air_frost_number == approx(0.0)
    fn.T_air_min = [-25.0]
    fn.T_air_max = [-5.0]
    fn.calculate_air_frost_number()
    assert fn.air_frost_number == approx(1.0)


def test_frostnumber_method_updates():
    """Test fn update()"""
    fn = frost_number.FrostnumberMethod()
    cfg_file = os.path.join(examples_directory, "Frostnumber_example_timeseries.cfg")
    fn.initialize(cfg_file=cfg_file)
    fn.update()
    assert fn.year == 2001
    assert fn.air_frost_number == approx(0.5)


def test_frostnumber_generates_output():
    """Test fn generates output file"""
    fn = frost_number.FrostnumberMethod()
    cfg_file = os.path.join(examples_directory, "Frostnumber_example_timeseries.cfg")
    fn.initialize(cfg_file=cfg_file)
    fn.update()
    output_written = fn.write_output_to_file()
    if output_written:
        assert os.path.isfile(fn.fn_out_filename)
        files_to_remove.append(fn.fn_out_filename)
    else:
        print("Unable to test output to: {}".format(fn.fn_out_filename))
