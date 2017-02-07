"""
test_perma_base.py
  tests of the perma_base component of permamodel
"""

from permamodel.components import frost_number
import os
import numpy as np
from .. import permamodel_directory, data_directory, examples_directory

def test_directory_names_are_set():
    assert(permamodel_directory is not None)

