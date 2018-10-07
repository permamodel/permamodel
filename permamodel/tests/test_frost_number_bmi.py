"""
test_frost_number_bmi.py
  tests of the frost_number component of permamodel using bmi API
"""
import os

import pytest

from permamodel import examples_directory
from permamodel.components import bmi_frost_number

# Set the file names for the example cfg files
onesite_oneyear_filename = \
        os.path.join(examples_directory,
                     'Frostnumber_example_scalar.cfg')
onesite_multiyear_filename = \
        os.path.join(examples_directory,
                     'Frostnumber_example_timeseries.cfg')

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
# ---------------------------------------------------
def test_can_initialize_bmi_frost_number_module():
    """ Test that the BMI frost number module can be initialized """
    fn = bmi_frost_number.BmiFrostnumberMethod
    assert fn is not None

def test_frost_number_has_initialize():
    """ Test that the BMI frost number module has initialize() """
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_oneyear_filename)

def test_frost_number_initialize_sets_year():
    """ help ensure proper initialization by checking for first year value """
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_oneyear_filename)

    # Assert the values from the cfg file
    assert fn.model.year == 2000

def test_frost_number_initialize_sets_air_min_and_max():
    """ Verify initialied values """
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_oneyear_filename)

    # Assert the values from the cfg file
    assert fn.model.T_air_min == -20.0
    assert fn.model.T_air_max == 10.0

def test_frost_number_update_increments_year():
    """ Test update increments time """
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    fn.update()
    assert fn.model.year == fn.model.start_year + fn.model.dt
    assert fn.model.year != fn.model.start_year

def test_frost_number_update_changes_air_frost_number():
    """ Test that value changes with update """
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    fn.update()
    afn0 = fn.model.air_frost_number
    fn.update()
    afn1 = fn.model.air_frost_number
    assert afn0 != afn1

def test_frost_number_runs_several_years():
    """ Test that frostnumber advances over several years """
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    while fn.model.year < fn.model.end_year:
        fn.update()

    assert fn.model.output is not None

    # Ensure that each year exists in the output dictionary
    year = fn.model.start_year
    while year < fn.model.end_year:
        assert year in fn.model.output.keys()
        year += 1

def test_frost_number_implements_update_until():
    """ Test that can update to arbitrary time """
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    fn.update_until(fn.model.end_year)

def test_frost_number_implements_finalize():
    """ Test BMI-required finalize() routine """
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    fn.update_until(fn.model.end_year)
    fn.finalize()
    files_to_remove.append(fn.model.fn_out_filename)

def test_frost_number_get_current_time_returns_scalar_float():
    """ Test that time is floating point """
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    current_time = fn.get_current_time()
    assert isinstance(current_time, float)

def test_frost_number_get_end_time_returns_scalar_float():
    """ Test that end time is floating point """
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    end_time = fn.get_end_time()
    assert isinstance(end_time, float)

def test_frostnumber_get_attribute():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    # Check an attribute that exists
    this_att = fn.get_attribute('time_units')
    assert this_att == 'years'

    # Check an attribute that doesn't exist
    with pytest.raises(KeyError):
        fn.get_attribute('not_an_attribute')

def test_frostnumber_update():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    fn.update()

def test_frostnumber_update_frac():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    fn.update_frac(2)

def test_bmi_fn_get_input_var_names():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    assert 'atmosphere_bottom_air__time_min_of_temperature' in fn.get_input_var_names()

def test_bmi_fn_get_output_var_names():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    assert 'frostnumber__air' in fn.get_output_var_names()

def test_bmi_fn_get_var_name():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    assert 'air_frost_number' in fn.get_var_name('frostnumber__air')

def test_bmi_fn_get_var_units():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    assert 'deg' in fn.get_var_units('atmosphere_bottom_air__time_min_of_temperature')
