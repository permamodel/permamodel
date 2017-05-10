"""
test_frost_number.py
  tests of the frost_number component of permamodel
"""

from permamodel.components import frost_number
from permamodel.components import bmi_frost_number
import os
import numpy as np
from .. import permamodel_directory, examples_directory

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
# Tests that input data is being read correctly
# ---------------------------------------------------
def test_can_initialize_frostnumber_method_from_file():
    """ Test fn initialization from config file """
    fn = bmi_frost_number.BmiFrostnumberMethod()
    cfg_file = os.path.join(examples_directory,
                            'Fairbanks_frostnumber_method.cfg')
    fn.initialize(cfg_file=cfg_file)
    files_to_remove.append(fn._model.fn_out_filename)

def test_frostnumber_method_has_date_info():
    """ Test that fn has time values """
    fn = bmi_frost_number.BmiFrostnumberMethod()
    cfg_file = os.path.join(examples_directory,
                            'Fairbanks_frostnumber_method.cfg')
    fn.initialize(cfg_file=cfg_file)
    assert(fn._model.year >= 0)
    assert(fn._model.year == fn._model.start_year)
