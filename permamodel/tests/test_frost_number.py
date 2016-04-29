"""
test_frost_number.py
  tests of the frost_number component of permamodel
"""

from permamodel.components import frost_number

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
def test_initialized_model_has_time_and_location():
    fn = frost_number.frostnumber_method()
    fn.initialize(cfg_file="./examples/Fairbanks_Ku_method.cfg")
    assert(False)


