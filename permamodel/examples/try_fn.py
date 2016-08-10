"""
try_fn.py

A test file to see if permamodel is installed correctly

Usage:
    python try_fn.py

"""

from permamodel.components import frost_number

# Set the file names for the example cfg files
onesite_oneyear_filename = \
        './permamodel/examples/Frostnumber_example_singlesite_singleyear.cfg'
onesite_multiyear_filename = \
        './permamodel/examples/Frostnumber_example_singlesite_multiyear.cfg'


fn = frost_number.frostnumber_method()
fn.initialize(cfg_file=onesite_multiyear_filename)
fn.update_until(fn.get_end_time())
fn.finalize()

print('Successfully finished running the Permamodel FrostNumber component!')
