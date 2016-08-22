"""
test_frost_number_bmi.py
  tests of the frost_number component of permamodel using bmi API
"""

from permamodel.components import frost_number
from permamodel.components import perma_base
import os
import numpy as np
from ..tests import examples_directory

# Set the file names for the example cfg files
onesite_oneyear_filename = \
        os.path.join(examples_directory,
                     'Frostnumber_example_singlesite_singleyear.cfg')
        #'./permamodel/examples/Frostnumber_example_singlesite_singleyear.cfg'
onesite_multiyear_filename = \
        os.path.join(examples_directory,
                     'Frostnumber_example_singlesite_multiyear.cfg')
        #'./permamodel/examples/Frostnumber_example_singlesite_multiyear.cfg'

# ---------------------------------------------------
# Tests that ensure we have bmi functionality
# ---------------------------------------------------
def test_frost_number_has_initialize():
    # Can we call an initialize function?
    fn = frost_number.FrostnumberMethod()
    # With hard-coded cfg filename
    #fn.initialize(cfg_file='/home/scotts/permamodel/permamodel/examples/Fairbanks_frostnumber_method.cfg', SILENT=True)
    # With relative cfg filename
    fn.initialize(cfg_file=onesite_oneyear_filename)

def test_frost_number_initialize_sets_year():
    fn = frost_number.FrostnumberMethod()
    fn.initialize(cfg_file=onesite_oneyear_filename)

    # Assert the values from the cfg file
    assert(fn.year == 2000)

def test_frost_number_initialize_sets_air_min_and_max():
    fn = frost_number.FrostnumberMethod()
    fn.initialize(cfg_file=onesite_oneyear_filename)

    # Assert the values from the cfg file
    assert(fn.T_air_min == -20.0)
    assert(fn.T_air_max == 10.0)

def test_frost_number_update_increments_year():
    fn = frost_number.FrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    fn.update()
    assert(fn.year == fn.start_year + fn.dt)
    assert(fn.year != fn.start_year)

def test_frost_number_update_changes_air_frost_number():
    fn = frost_number.FrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    afn0 = fn.air_frost_number
    fn.update()
    afn1 = fn.air_frost_number
    assert(afn0 != afn1)

def test_frost_number_runs_several_years():
    fn = frost_number.FrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    while fn.year < fn.end_year:
        fn.update()

    assert(fn.output is not None)

    # Ensure that each year exists in the output dictionary
    year = fn.start_year
    while year < fn.end_year:
        assert(year in fn.output.keys())
        year += 1

    # print the output to check it
    #fn.print_frost_numbers()

def test_frost_number_implements_update_until():
    fn = frost_number.FrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    fn.update_until(fn.end_year)

    # Check the output, if desired
    #fn.print_frost_numbers()

def test_frost_number_implements_finalize():
    fn = frost_number.FrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    fn.update_until(fn.end_year)
    fn.finalize()

def test_frost_number_get_current_time_returns_scalar_float():
    fn = frost_number.FrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    current_time = fn.get_current_time()
    assert(isinstance(current_time, float))

def  test_frost_number_get_end_time_returns_scalar_float():
    fn = frost_number.FrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    end_time = fn.get_end_time()
    assert(isinstance(end_time, float))




