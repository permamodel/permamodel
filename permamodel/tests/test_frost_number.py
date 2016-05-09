"""
test_frost_number.py
  tests of the frost_number component of permamodel
"""

from permamodel.components import frost_number
import os

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
    fn.initialize(cfg_file="/home/scotts/permamodel/permamodel/examples/Fairbanks_frostnumber_method.cfg")

def test_frostnumber_method_has_date_and_location():
    fn = frost_number.frostnumber_method()
    fn.initialize(cfg_file="/home/scotts/permamodel/permamodel/examples/Fairbanks_frostnumber_method.cfg")
    #print("start_year: %d" % fn.start_year)
    #print("year: %d" % fn.year)
    #print("lon: %f" % fn.lon)
    #print("lat: %f" % fn.lat)
    assert(fn.year >= 0)
    assert(fn.year == fn.start_year)

def test_can_get_temperature_filename_from_date_and_location():
    fn = frost_number.frostnumber_method()
    fn.initialize(cfg_file="/home/scotts/permamodel/permamodel/examples/Fairbanks_frostnumber_method.cfg")
    # Given year and month, calculate the tiff_filename
    fname = fn.get_temperature_tiff_filename(fn.year, 6)
    print("fname: %s" % fname)



