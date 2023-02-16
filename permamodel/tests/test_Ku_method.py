"""Test Ku model"""
import os

import pytest

from permamodel import examples_directory
from permamodel.components.Ku_method import Ku_method

CFG_FILENAME = "Ku_method.cfg"
CFG_FILE = os.path.join(examples_directory, CFG_FILENAME)


def test_instantiate_Ku():
    m = Ku_method()
    assert type(m) is Ku_method


def test_not_initialized_Ku():
    with pytest.raises(AttributeError):
        m = Ku_method()
        m.time


def test_initialize_Ku():
    m = Ku_method()
    m.initialize(CFG_FILE)
    assert m.time == 0.0


# See #78
# def test_update_Ku():
#     m = Ku_method()
#     m.initialize(CFG_FILE)
#     m.update()
#     assert m.time == 1.0


def test_finalize_Ku():
    m = Ku_method()
    m.initialize(CFG_FILE)
    m.finalize()
