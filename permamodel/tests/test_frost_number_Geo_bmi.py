"""
test_frost_number_bmi.py
  tests of the frost_number component of permamodel using bmi API
"""
import datetime
import os

import numpy as np
import pytest
from dateutil.relativedelta import relativedelta
from pytest import approx

from permamodel.components import bmi_frost_number_Geo, frost_number_Geo, perma_base

from .. import examples_directory

# Set the file names for the example cfg files
config_filename = os.path.join(examples_directory, "FrostnumberGeo_Default.cfg")

# List of files to be removed after testing is complete
# use files_to_remove.append(<filename>) to add to it
files_to_remove = []


def setup_module():
    """Standard fixture called before any tests in this file are performed"""
    pass


def teardown_module():
    """Standard fixture called after all tests in this file are performed"""
    # If need to remove files that are created:
    # for f in files_to_remove:
    #    if os.path.exists(f):
    #        os.remove(f)
    pass


# ---------------------------------------------------
# Tests that ensure we have bmi functionality
# Note: the netcdf functions seem to lead to errors if the netcdf
# files are not closed.  This is done in the finalize() routine
# which should therefore be called after every test
# ---------------------------------------------------
def test_frost_number_Geo_has_initialize(tmpdir):
    # Can we call an initialize function?
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize(cfg_file=config_filename)
        fng.finalize()  # Must have this or get IOError later


def test_frost_number_initialize_sets_year(tmpdir):
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize(cfg_file=config_filename)

        # Assert the values from the cfg file
        assert fng._model._date_current.year == 1901
        fng.finalize()  # Must have this or get IOError later


def test_frost_number_initialize_sets_air_min_and_max(tmpdir):
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize(cfg_file=config_filename)

        # The temperature arrays are initialized to all NaNs
        nan_array = np.zeros((3, 2), dtype=np.float32)
        nan_array.fill(np.nan)
        assert np.all(np.isnan(fng._model.T_air_min))  # == nan_array)
        fng.finalize()  # Must have this or get IOError later


def test_frost_number_update_increments_time(tmpdir):
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()
        fng.update()

        assert fng._model._date_current == fng._model._start_date + relativedelta(
            years=fng._model._timestep_duration
        )
        fng.finalize()  # Must have this or get IOError later


def test_frost_number_update_changes_air_frost_number(tmpdir):
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()

        afn0 = fng._model.air_frost_number_Geo.copy()
        fng.update()
        afn1 = fng._model.air_frost_number_Geo.copy()
        assert np.any(afn0 != afn1)
        fng.update_until(1.0)
        afn2 = fng._model.air_frost_number_Geo.copy()
        assert np.all(afn1 == afn2)

        fng.finalize()  # Must have this or get IOError later


def test_frost_number_get_current_time_returns_scalar(tmpdir):
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()
        current_time = fng.get_current_time()
        assert isinstance(current_time, (float, int))
        fng.finalize()  # Must have this or get IOError later


def test_frost_number_get_end_time_returns_scalar(tmpdir):
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()
        end_time = fng.get_end_time()
        assert isinstance(end_time, float) or isinstance(end_time, int)
        fng.finalize()  # Must have this or get IOError later


def test_frost_number_implements_update_until(tmpdir):
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()

        fng.update_until(fng.get_end_time())
        assert fng._model._date_current == fng._model._end_date

        fng.finalize()  # Must have this or get IOError later


def test_FNGeo_computes_default_values(tmpdir):
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()
        for n in range(5):
            if (fng._model.T_air_min[0, 0] + fng._model.T_air_max[0, 0]) == 0:
                assert fng._values["frostnumber__air"][0, 0] == 0.5
            if (fng._model.T_air_min[0, 0] <= 0.0) and (
                fng._model.T_air_max[0, 0] <= 0.0
            ):
                assert fng._values["frostnumber__air"][0, 0] == 1.0
            if (fng._model.T_air_min[0, 0] > 0.0) and (
                fng._model.T_air_max[0, 0] > 0.0
            ):
                assert fng._values["frostnumber__air"][0, 0] == 0.0
            # print("FNGeo frostnumber__air: %d" % n)
            # print(fng._model.T_air_min)
            # print(fng._model.T_air_max)
            # print(fng._values['frostnumber__air'])
            fng.update()
        fng.finalize()  # Must have this or get IOError later


def test_FNGeo_get_attribute(tmpdir):
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()
        # Check an attribute that exists
        this_att = fng.get_attribute("time_units")
        assert this_att == "years"

        # Check an attribute that doesn't exist
        with pytest.raises(KeyError):
            fng.get_attribute("not_an_attribute")

        fng.finalize()  # Must have this or get IOError later


def test_FNGeo_get_input_var_names(tmpdir):
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()
        assert "atmosphere_bottom_air__temperature" in fng.get_input_var_names()

        fng.finalize()  # Must have this or get IOError later


def test_FNGeo_get_output_var_names(tmpdir):
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()
        assert "frostnumber__air" in fng.get_output_var_names()

        fng.finalize()  # Must have this or get IOError later


def test_FNGeo_get_set_value(tmpdir):
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()

        airtempval = fng.get_value("atmosphere_bottom_air__temperature")
        airtempnew = 123 * np.ones_like(airtempval)
        fng.set_value("atmosphere_bottom_air__temperature", airtempnew)
        assert not np.all(airtempval == approx(airtempnew))

        fng.finalize()  # Must have this or get IOError later


def test_FNGeo_get_grid_shape(tmpdir):
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()

        airtemp_gridid = fng.get_var_grid("atmosphere_bottom_air__temperature")
        assert (3, 2) == fng.get_grid_shape(airtemp_gridid)

        fng.finalize()  # Must have this or get IOError later


def test_FNGeo_get_grid_size(tmpdir):
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()

        airtemp_gridid = fng.get_var_grid("atmosphere_bottom_air__temperature")
        assert fng.get_grid_size(airtemp_gridid) == 6

        fng.finalize()  # Must have this or get IOError later


def test_FNGeo_get_grid_spacing(tmpdir):
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()

        airtemp_gridid = fng.get_var_grid("atmosphere_bottom_air__temperature")
        assert np.all(np.array([1.0, 1.0]) == fng.get_grid_spacing(airtemp_gridid))

        fng.finalize()  # Must have this or get IOError later


def test_FNGeo_jan_and_jul_temperatures_are_grids(tmpdir):
    """Test that FNGeo BMI has input variables for jan and jul data"""
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()

        airtemp_gridid = fng.get_var_grid("atmosphere_bottom_air__temperature")
        jan_airtemp_gridid = fng.get_var_grid(
            "atmosphere_bottom_air__temperature_mean_jan"
        )
        assert jan_airtemp_gridid is not None
        jul_airtemp_gridid = fng.get_var_grid(
            "atmosphere_bottom_air__temperature_mean_jul"
        )
        assert jul_airtemp_gridid is not None

        fng.finalize()  # Must have this or get IOError later


def test_FNGeo_can_set_current_and_jan_temperatures(tmpdir):
    """Test that FNGeo BMI can set jan temperature field"""
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()

        airtemp_values = fng.get_value("atmosphere_bottom_air__temperature")
        airtemp_values = np.ones_like(airtemp_values)
        fng.set_value("atmosphere_bottom_air__temperature", airtemp_values)

        jan_airtemp_values = fng.get_value(
            "atmosphere_bottom_air__temperature_mean_jan"
        )
        jan_airtemp_values = np.ones_like(jan_airtemp_values)
        fng.set_value("atmosphere_bottom_air__temperature_mean_jan", jan_airtemp_values)
        assert np.all(
            fng.get_value("atmosphere_bottom_air__temperature")
            == fng.get_value("atmosphere_bottom_air__temperature_mean_jan")
        )

        fng.finalize()  # Must have this or get IOError later


def test_FNGeo_update_zero_fraction_does_not_change_time(tmpdir):
    """Test that running update_frac(0) does not change the time"""
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        fng.initialize()
        init_time = fng.get_current_time()
        fng.update_frac(0)
        plus_zero_time = fng.get_current_time()
        assert init_time == plus_zero_time
        fng.update()
        one_update_time = fng.get_current_time()
        assert not np.all(init_time == one_update_time)

        fng.finalize()  # Must have this or get IOError later


def test_FNGeo_simulate_WMT_run(tmpdir):
    """Test that we can set values as if running in WMT"""
    with tmpdir.as_cwd():
        fng = bmi_frost_number_Geo.BmiFrostnumberGeoMethod()
        wmt_cfg_file = os.path.join(examples_directory, "FNGeo_WMT_testing.cfg")
        fng.initialize(wmt_cfg_file)
        assert fng._name == "Permamodel FrostnumberGeo Component"

        # Until set, e.g. by WMT with cru values, all vals are NaN
        # Note: these are actually setting references to the underlying
        #       model arrays!
        airtemp_values = fng.get_value("atmosphere_bottom_air__temperature")
        jan_airtemp_values = fng.get_value(
            "atmosphere_bottom_air__temperature_mean_jan"
        )
        jul_airtemp_values = fng.get_value(
            "atmosphere_bottom_air__temperature_mean_jul"
        )

        # In WMT mode, must set monthly temperature values, then run update_frac(0)
        #   to get valid values in frost number array
        # use January as 'coldest'
        # use July as 'warmest'
        airtemps_of_one = np.ones_like(airtemp_values)
        jan_airtemp_values[:] = -10.0 * airtemps_of_one
        jul_airtemp_values[:] = 10.0 * airtemps_of_one
        fng.set_value("atmosphere_bottom_air__temperature_mean_jan", jan_airtemp_values)
        fng.set_value("atmosphere_bottom_air__temperature_mean_jul", jul_airtemp_values)
        fng.update_frac(0)
        airfn_values = fng.get_value("frostnumber__air")

        assert np.all(0.5 * airtemps_of_one == airfn_values)

        fng.finalize()  # Must have this or get IOError later
