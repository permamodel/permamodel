""" tests subpackage for permamodel

Provides:
    permamodel_directory: the 'root' directory of the permamodel module
    data_directory: the directory where installed data files are found
    examples_directory: the directory where installed data files are found
"""

import os

tests_directory = os.path.dirname(__file__)
permamodel_directory = os.path.join(tests_directory, '..')
data_directory = os.path.join(permamodel_directory, 'data')
examples_directory = os.path.join(permamodel_directory, 'examples')
