"""
try_fn.py

A test file to see if permamodel is installed correctly

Usage:
    python try_fn.py

"""

import os

from permamodel import examples_directory
from permamodel.components import bmi_frost_number, frost_number

# Set the file names for the example cfg files
onesite_singleyear_filename = os.path.join(
    examples_directory, "Frostnumber_example_singlesite_singleyear.cfg"
)
onesite_multiyear_filename = os.path.join(
    examples_directory, "Frostnumber_example_singlesite_multiyear.cfg"
)


fn = bmi_frost_number.BmiFrostnumberMethod()
fn.initialize(cfg_file=onesite_singleyear_filename)
fn.update_until(fn.get_end_time())
fn.finalize()
print("Successfully finished running singleyear FrostNumber component")

fn = bmi_frost_number.BmiFrostnumberMethod()
fn.initialize(cfg_file=onesite_multiyear_filename)
fn.update_until(fn.get_end_time())
fn.finalize()
print("Successfully finished running multiyear FrostNumber component")
