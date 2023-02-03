"""
test_frost_number_bmi.py
  tests of the frost_number component of permamodel using bmi API
"""
import os

import pytest
from pytest import approx

from permamodel import examples_directory
from permamodel.components import bmi_frost_number

# Set the file names for the example cfg files
onesite_oneyear_filename = os.path.join(
    examples_directory, "Frostnumber_example_scalar.cfg"
)
onesite_multiyear_filename = os.path.join(
    examples_directory, "Frostnumber_example_timeseries.cfg"
)

# List of files to be removed after testing is complete
# use files_to_remove.append(<filename>) to add to it
files_to_remove = []


def setup_module():
    """Standard fixture called before any tests in this file are performed"""
    pass


def teardown_module():
    """Standard fixture called after all tests in this file are performed"""
    for f in files_to_remove:
        if os.path.exists(f):
            os.remove(f)


_CONFIG_FILE_TEMPLATE = """
start_year      | {start_year} | int      | begining of the simulation time [year]
end_year        | {end_year}   | int      | begining of the simulation time [year]
dt              | 1            | int      | timestep for permafrost process [year]
T_air_min_type  | Scalar       | string   | allowed input types
T_air_min       | {t_air_min}  | float    | Mean annual air temperature [C]
T_air_max_type  | Scalar       | string   | allowed input types
T_air_max       | {t_air_max}  | float    | Mean annual air temperature [C]
#===============================================================================
# Output
fn_out_filename | fn_test_output.dat | string   | Name of the file to output the frostnumber data to
"""


def write_config_file(
    config_file, start_year=2000, end_year=2010, t_air_min=-20, t_air_max=20
):
    """Write a default configuration file for FrostNumberMethod."""
    values = {
        "start_year": start_year,
        "end_year": end_year,
        "t_air_min": t_air_min,
        "t_air_max": t_air_max,
    }
    with open(config_file, "w") as fp:
        fp.write(_CONFIG_FILE_TEMPLATE.format(**values))

    return values


# ---------------------------------------------------
# Tests that ensure we have bmi functionality
# ---------------------------------------------------
def test_can_initialize_bmi_frost_number_module():
    """Test that the BMI frost number module can be initialized"""
    fn = bmi_frost_number.BmiFrostnumberMethod
    assert fn is not None


def test_frost_number_has_initialize():
    """Test that the BMI frost number module has initialize()"""
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_oneyear_filename)


def test_frost_number_initialize_sets_year():
    """help ensure proper initialization by checking for first year value"""
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_oneyear_filename)

    # Assert the values from the cfg file
    assert fn.model.year == 2000


def test_frost_number_initialize_sets_air_min_and_max():
    """Verify initialied values"""
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_oneyear_filename)

    # Assert the values from the cfg file
    assert fn.model.T_air_min == -20.0
    assert fn.model.T_air_max == 10.0


def test_frost_number_update_increments_year():
    """Test update increments time"""
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    fn.update()
    assert fn.model.year == fn.model.start_year + fn.model.dt
    assert fn.model.year != fn.model.start_year


def test_frost_number_update_changes_air_frost_number():
    """Test that value changes with update"""
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    fn.update()
    afn0 = fn.model.air_frost_number
    fn.update()
    afn1 = fn.model.air_frost_number
    assert afn0 != afn1


def test_frost_number_runs_several_years():
    """Test that frostnumber advances over several years"""
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    while fn.model.year < fn.model.end_year:
        fn.update()

    assert fn.model.output is not None

    # Ensure that each year exists in the output dictionary
    year = fn.model.start_year
    while year < fn.model.end_year:
        assert year in fn.model.output.keys()
        year += 1


def test_frost_number_implements_update_until():
    """Test that can update to arbitrary time"""
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)

    fn.update_until(fn.model.end_year)


def test_frost_number_implements_finalize():
    """Test BMI-required finalize() routine"""
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    fn.update_until(fn.model.end_year)
    fn.finalize()
    files_to_remove.append(fn.model.fn_out_filename)


def test_frost_number_get_current_time_returns_scalar_float():
    """Test that time is floating point"""
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    current_time = fn.get_current_time()
    assert isinstance(current_time, float)


def test_frost_number_get_end_time_returns_scalar_float():
    """Test that end time is floating point"""
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    end_time = fn.get_end_time()
    assert isinstance(end_time, float)


def test_frostnumber_get_attribute():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    # Check an attribute that exists
    this_att = fn.get_attribute("time_units")
    assert this_att == "years"

    # Check an attribute that doesn't exist
    with pytest.raises(KeyError):
        fn.get_attribute("not_an_attribute")


def test_frostnumber_update():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    fn.update()


def test_frostnumber_update_frac():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    fn.update_frac(2)


def test_bmi_fn_get_input_var_names():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    assert "atmosphere_bottom_air__time_min_of_temperature" in fn.get_input_var_names()


def test_bmi_fn_get_output_var_names():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    assert "frostnumber__air" in fn.get_output_var_names()


def test_bmi_fn_get_var_name():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    assert "air_frost_number" in fn.get_var_name("frostnumber__air")


def test_bmi_fn_get_var_units():
    fn = bmi_frost_number.BmiFrostnumberMethod()
    fn.initialize(cfg_file=onesite_multiyear_filename)
    assert "deg" in fn.get_var_units("atmosphere_bottom_air__time_min_of_temperature")


def test_frostnumber_updating_with_scalars(tmpdir):
    """Test updating past one year with scalars."""
    start_year, end_year = 2000, 2010
    with tmpdir.as_cwd():
        write_config_file("frost_number.cfg", start_year=start_year, end_year=end_year)

        fn = bmi_frost_number.BmiFrostnumberMethod()
        fn.initialize(cfg_file="frost_number.cfg")
        for year in range(start_year, end_year):
            assert fn.get_value("frostnumber__air") == approx(0.5)
            assert fn.get_current_time() == (year - start_year)
            fn.update()


def test_frostnumber_set_value_with_scalars(tmpdir):
    """Test set_values changes the frost number to the correct value."""
    with tmpdir.as_cwd():
        write_config_file("frost_number.cfg", t_air_max=10.0)

        fn = bmi_frost_number.BmiFrostnumberMethod()
        fn.initialize(cfg_file="frost_number.cfg")
        fn.update()
        expected_value = fn.get_value("frostnumber__air")

        write_config_file("frost_number.cfg", t_air_max=999.0)
        fn = bmi_frost_number.BmiFrostnumberMethod()
        fn.initialize(cfg_file="frost_number.cfg")
        fn.update()

        fn.set_value("atmosphere_bottom_air__time_max_of_temperature", 10)
        fn.update()
        assert fn.get_value("frostnumber__air") == approx(expected_value)
