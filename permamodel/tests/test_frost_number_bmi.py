"""
test_frost_number_bmi.py
  tests of the frost_number component of permamodel using bmi API
"""

from permamodel.components import frost_number
from permamodel.components import bmi_frost_number
from permamodel.components import perma_base
import os
import numpy as np
from .. import examples_directory

# Set the file names for the example cfg files
onesite_oneyear_filename = \
        os.path.join(examples_directory,
                     'Frostnumber_example_singlesite_singleyear.cfg')
        #'./permamodel/examples/Frostnumber_example_singlesite_singleyear.cfg'
onesite_multiyear_filename = \
        os.path.join(examples_directory,
                     'Frostnumber_example_singlesite_multiyear.cfg')
        #'./permamodel/examples/Frostnumber_example_singlesite_multiyear.cfg'

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
def test_frost_number_has_initialize():
    # Can we call an initialize function?
    fn = bmi_frost_number.BmiFrostnumberMethod()
    # With hard-coded cfg filename
    #fn.initialize(cfg_file='/home/scotts/permamodel/permamodel/examples/Fairbanks_frostnumber_method.cfg', SILENT=True)
    # With relative cfg filename
    fn.initialize(cfg_file=onesite_oneyear_filename)

def test_frost_number_initialize_sets_year():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_oneyear_filename)

    # Assert the values from the cfg file
    assert(fn._model.year == 2000)

def test_frost_number_initialize_sets_air_min_and_max():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_oneyear_filename)

    # Assert the values from the cfg file
    assert(fn._model.T_air_min == -20.0)
    assert(fn._model.T_air_max == 10.0)

def test_frost_number_update_increments_year():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    fn.update()
    assert(fn._model.year == fn._model.start_year + fn._model.dt)
    assert(fn._model.year != fn._model.start_year)

def test_frost_number_update_changes_air_frost_number():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    fn.update()
    afn0 = fn._model.air_frost_number
    fn.update()
    afn1 = fn._model.air_frost_number
    assert(afn0 != afn1)

def test_frost_number_runs_several_years():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    while fn._model.year < fn._model.end_year:
        fn.update()

    assert(fn._model.output is not None)

    # Ensure that each year exists in the output dictionary
    year = fn._model.start_year
    while year < fn._model.end_year:
        assert(year in fn._model.output.keys())
        year += 1

    # print the output to check it
    #fn.print_frost_numbers()

def test_frost_number_implements_update_until():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    fn.update_until(fn._model.end_year)

    # Check the output, if desired
    #fn.print_frost_numbers()

def test_frost_number_implements_finalize():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    fn.update_until(fn._model.end_year)
    fn.finalize()

def test_frost_number_get_current_time_returns_scalar_float():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    current_time = fn.get_current_time()
    assert(isinstance(current_time, float))

def  test_frost_number_get_end_time_returns_scalar_float():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    end_time = fn.get_end_time()
    assert(isinstance(end_time, float))




