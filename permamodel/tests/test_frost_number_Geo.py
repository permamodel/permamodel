"""
test_frost_number_Geo.py
  tests of the Geo version of the frost_number component of permamodel
"""
import datetime
import os

import numpy as np
import pytest

from permamodel.components import frost_number_Geo

from .. import examples_directory, permamodel_directory

# List of files to be removed after testing is complete
# use files_to_remove.append(<filename>) to add to it
files_to_remove = []
files_cfg_file = os.path.join(examples_directory, "FrostnumberGeo_Files.cfg")


def setup_module():
    """Standard fixture called before any tests in this file are performed"""
    pass


def teardown_module():
    """Standard fixture called after all tests in this file are performed"""
    for f in files_to_remove:
        if os.path.exists(f):
            os.remove(f)


def test_can_initialize_Geo_frostnumber_module():
    fn_geo = frost_number_Geo.FrostnumberGeoMethod()


def test_Geo_frostnumber_has_default_config_file():
    fn_geo = frost_number_Geo.FrostnumberGeoMethod()
    assert fn_geo._config_filename is not None


def test_Geo_frostnumber_can_be_passed_config_filename():
    fn_geo = frost_number_Geo.FrostnumberGeoMethod(cfgfile="a file")
    assert fn_geo._config_filename == "a file"


# Configuration from "Files" not supported for WMT version
# def test_Geo_frostnumber_initializes_from_files_config_file():
#    fn_geo = frost_number_Geo.FrostnumberGeoMethod(cfgfile=files_cfg_file)
#    assert os.path.isfile(fn_geo._config_filename)
#    fn_geo.initialize_frostnumberGeo_component()
#    assert fn_geo._grid_type == 'uniform rectilinear'
#    assert fn_geo._calc_surface_fn is not None
#    assert fn_geo._calc_stefan_fn is not None
#    assert fn_geo._dd_method in ('ObservedMinMax', 'MinJanMaxJul', 'MonthlyAverages', 'DailyValues')
#    if fn_geo._dd_method == 'MinJanMaxJul':
#        assert fn_geo.T_air_min.shape == fn_geo._grid_shape
#        assert fn_geo.T_air_max.shape == fn_geo._grid_shape
#        if fn_geo._using_Files:
#            assert fn_geo._temperature_dataset is not None
#
#    fn_geo.finalize_frostnumber_Geo()
#    files_to_remove.append(fn_geo.output_filename)


def test_Geo_frostnumber_initialize_datacube():
    fn_geo = frost_number_Geo.FrostnumberGeoMethod()
    config_dict = {
        "n_temperature_grid_fields": 4,
        "temperature_grid_date_0": "1901-01-01",
        "temperature_grid_data_0": "((-10, -5), (-20, -15), (0, 5))",
        "temperature_grid_date_1": "1901-07-01",
        "temperature_grid_data_1": "((10, 15), (0, 5), (20, 15))",
        "temperature_grid_date_2": "1902-01-01",
        "temperature_grid_data_2": "((-7, -2), (-17, -12), (3, 8))",
        "temperature_grid_date_3": "1902-07-01",
        "temperature_grid_data_3": "((7, 2), (17, 12), (23, 28))",
    }
    dates, cube = fn_geo.initialize_datacube("temperature", config_dict)

    assert datetime.date(1901, 7, 1) in dates  # second date

    assert cube[0, 0, 0] == -10  # very first value
    assert cube[3, 2, 1] == 28  # very last value
    assert cube[2, 1, 0] == -17  # 3rd date, 2nd set, 1st value


def test_Geo_frostnumber_get_datacube_slice():
    fn_geo = frost_number_Geo.FrostnumberGeoMethod()
    config_dict = {
        "n_temperature_grid_fields": 4,
        "temperature_grid_date_0": datetime.date(1901, 1, 1),
        "temperature_grid_data_0": "((-10, -5), (-20, -15), (0, 5))",
        "temperature_grid_date_1": datetime.date(1901, 7, 1),
        "temperature_grid_data_1": "((10, 15), (0, 5), (20, 15))",
        "temperature_grid_date_2": datetime.date(1902, 1, 1),
        "temperature_grid_data_2": "((-7, -2), (-17, -12), (3, 8))",
        "temperature_grid_date_3": datetime.date(1902, 7, 1),
        "temperature_grid_data_3": "((7, 2), (17, 12), (23, 28))",
    }
    dates, cube = fn_geo.initialize_datacube("temperature", config_dict)

    test_date = datetime.date(1901, 2, 27)
    test_field = fn_geo.get_datacube_slice(test_date, cube, dates)
    assert test_field[0, 0] == -10
    assert test_field[1, 1] == -15

    test_date = datetime.date(1902, 3, 10)
    test_field = fn_geo.get_datacube_slice(test_date, cube, dates)
    assert test_field[1, 0] == -17
    assert test_field[2, 1] == 8

    test_date = datetime.date(1900, 3, 10)
    with pytest.raises(ValueError):
        fn_geo.get_datacube_slice(test_date, cube, dates)

    test_date = datetime.date(1980, 3, 10)
    with pytest.raises(ValueError):
        fn_geo.get_datacube_slice(test_date, cube, dates)


def test_Geo_frostnumber_initializes_from_default_config_file():
    fn_geo = frost_number_Geo.FrostnumberGeoMethod()
    assert os.path.isfile(fn_geo._config_filename)
    fn_geo.initialize_frostnumberGeo_component()
    assert fn_geo._grid_type == "uniform rectilinear"
    assert fn_geo._calc_surface_fn is not None
    assert fn_geo._calc_stefan_fn is not None
    assert fn_geo._dd_method in (
        "ObservedMinMax",
        "MinJanMaxJul",
        "MonthlyAverages",
        "DailyValues",
    )
    if fn_geo._dd_method == "MinJanMaxJul":
        assert fn_geo.T_air_min.shape == fn_geo._grid_shape
        assert fn_geo.T_air_max.shape == fn_geo._grid_shape
        if fn_geo._using_Files:
            assert fn_geo._temperature_dataset is not None

    fn_geo.finalize_frostnumber_Geo()
    files_to_remove.append(fn_geo.output_filename)


def test_Geo_frostnumber_can_compute_real_date_from_timestep():
    fn_geo = frost_number_Geo.FrostnumberGeoMethod()
    fn_geo.initialize_frostnumberGeo_component()

    # Timestep should be one year
    assert fn_geo._timestep_duration == 1

    # Model reference date should have timestep of zero
    assert fn_geo.get_timestep_from_date(fn_geo._reference_date) == 0

    # Timestep of first date should be first timestep
    assert fn_geo._timestep_first == fn_geo.get_timestep_from_date(fn_geo._start_date)

    # Last date should be date of last timestep
    assert fn_geo._end_date == fn_geo.get_date_from_timestep(fn_geo._timestep_last)

    fn_geo.finalize_frostnumber_Geo()


def test_Geo_frostnumber_can_return_temperature_field():
    fn_geo = frost_number_Geo.FrostnumberGeoMethod()
    fn_geo.initialize_frostnumberGeo_component()

    current_temperature_field = fn_geo.get_temperature_field()
    assert current_temperature_field.shape == fn_geo._grid_shape

    bad_date_field = fn_geo.get_temperature_field(datetime.date(100, 1, 1))
    assert np.all(np.isnan(bad_date_field))

    fn_geo.finalize_frostnumber_Geo()


def test_Geo_frostnumber_get_latest_min_max_months():
    fn_geo = frost_number_Geo.FrostnumberGeoMethod()
    fn_geo.initialize_frostnumberGeo_component()

    this_date = datetime.date(1905, 1, 1)
    (mindate, maxdate) = fn_geo.get_min_and_max_dates(this_date)
    assert mindate == datetime.date(1905, 1, 15)
    assert maxdate == datetime.date(1904, 7, 15)

    this_date = datetime.date(2004, 8, 1)
    (mindate, maxdate) = fn_geo.get_min_and_max_dates(this_date)
    assert mindate == datetime.date(2004, 1, 15)
    assert maxdate == datetime.date(2004, 7, 15)

    fn_geo.finalize_frostnumber_Geo()


def test_Geo_frostnumber_compute_array_of_degree_days():
    fn_geo = frost_number_Geo.FrostnumberGeoMethod()
    fn_geo.initialize_frostnumberGeo_component()

    if fn_geo._using_WMT:
        print("Can't test WMT in standalone mode")
        return
    elif not (fn_geo._using_Files or fn_geo._using_ConfigVals):
        raise ValueError(
            "Frostnumber must use either Files, ConfigVals, \
                         or WMT to get input variables"
        )

    # Test that we get NaN-filled arrays when no input vars available
    fn_geo.set_current_date_and_timestep_with_timestep(-100)
    fn_geo.get_input_vars()
    assert fn_geo._date_current < fn_geo._temperature_first_date
    if fn_geo._dd_method == "MinJanMaxJul":
        # A timestep out of bounds should yield nans for temp and dd
        assert np.all(np.isnan(fn_geo.T_air_min))
        assert np.all(np.isnan(fn_geo.T_air_max))

        # Calculating degree days on all NaNs yields all NaNs
        fn_geo.compute_degree_days()
        assert np.all(np.isnan(fn_geo.ddt))
        assert np.all(np.isnan(fn_geo.ddf))

        # Calculating air frost number on NaNs yields NaNs
        fn_geo.compute_air_frost_number_Geo()
        assert np.all(np.isnan(fn_geo.air_frost_number_Geo))

    # Test that we get real values
    fn_geo.set_current_date_and_timestep_with_timestep(2)
    fn_geo.get_input_vars()
    assert fn_geo._date_current >= fn_geo._temperature_first_date
    if fn_geo._dd_method == "MinJanMaxJul":
        # The default should have no missing temperature values
        assert not np.any(np.isnan(fn_geo.T_air_min))
        assert not np.any(np.isnan(fn_geo.T_air_max))
        # fn_geo.T_air_min.tofile("Tairmin.dat")
        # fn_geo.T_air_max.tofile("Tairmax.dat")

        # Calculating degree days on all NaNs yields all NaNs
        fn_geo.compute_degree_days()
        assert not np.any(np.isnan(fn_geo.ddt))
        assert not np.any(np.isnan(fn_geo.ddf))
        # fn_geo.ddt.tofile("ddt.dat")
        # fn_geo.ddf.tofile("ddf.dat")

        # Calculating air frost number on NaNs yields NaNs
        fn_geo.compute_air_frost_number_Geo()
        assert not np.any(np.isnan(fn_geo.air_frost_number_Geo))
        # fn_geo.air_frost_number_Geo.tofile("afn_geod.dat")
    fn_geo.finalize_frostnumber_Geo()


def test_Geo_frostnumber_output_a_netcdf_file():
    fn_geo = frost_number_Geo.FrostnumberGeoMethod()
    fn_geo.initialize_frostnumberGeo_component()
    assert fn_geo.output_filename is not None
    assert fn_geo.output_filename[-3:] == ".nc"
    fn_geo.initial_update()
    for t in range(0, 10):
        fn_geo.get_input_vars()
        fn_geo.compute_degree_days()
        fn_geo.calculate_frost_numbers_Geo()
        fn_geo.add_to_output()
        fn_geo.update()
    fn_geo.finalize()


def test_Geo_frostnumber_update_until_timestep():
    fn_geo = frost_number_Geo.FrostnumberGeoMethod()
    fn_geo.initialize_frostnumberGeo_component()
    fn_geo.initial_update()
    fn_geo.update_until_timestep(fn_geo._timestep_last)
    fn_geo.finalize()
