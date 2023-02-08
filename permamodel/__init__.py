""" permamodel package

Provides:
    permamodel_directory: the 'root' directory of the permamodel module
    data_directory: the directory where installed data files are found
    examples_directory: the directory where installed data files are found
    tests_directory: location of unit tests
"""

import os

from ._version import __version__

permamodel_directory = os.path.dirname(__file__)
data_directory = os.path.join(permamodel_directory, "data")
examples_directory = os.path.join(permamodel_directory, "examples")
tests_directory = os.path.join(permamodel_directory, "tests")

__all__ = ["__version__"]
