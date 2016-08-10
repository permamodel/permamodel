"""
test_frost_number.py
  tests of the frost_number component of permamodel
"""

from permamodel.components import frost_number
import os
import numpy as np
import pprint

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
    fn.initialize(cfg_file=cfg_file)

def test_frostnumber_method_has_date_and_location():
    fn = frost_number.frostnumber_method()
    cfg_file = os.path.join(os.environ.get('PERMAMODEL_EXAMPLEDIR',\
                                '/examples/'),
                 'Fairbanks_frostnumber_method.cfg')
    fn.initialize(cfg_file=cfg_file)
    assert(fn.year >= 0)
    assert(fn.year == fn.start_year)

