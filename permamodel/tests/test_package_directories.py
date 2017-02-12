"""Tests directories set in the permamodel package definition file."""

import os
from nose.tools import assert_true
from .. import (permamodel_directory, data_directory,
                examples_directory, tests_directory)


def test_permamodel_directory_is_set():
    assert(permamodel_directory is not None)


def test_data_directory_is_set():
    assert(data_directory is not None)


def test_examples_directory_is_set():
    assert(examples_directory is not None)


def test_tests_directory_is_set():
    assert(tests_directory is not None)


def test_permamodel_directory_exists():
    assert_true(os.path.isdir(permamodel_directory))


def test_data_directory_exists():
    assert_true(os.path.isdir(data_directory))


def test_examples_directory_exists():
    assert_true(os.path.isdir(examples_directory))


def test_tests_directory_exists():
    assert_true(os.path.isdir(tests_directory))
