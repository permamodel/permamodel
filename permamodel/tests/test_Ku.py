import os

import pytest
from netCDF4 import Dataset
from numpy.testing import assert_approx_equal, assert_array_equal

from permamodel.components.Ku import Ku_model


def test_always_passes():
    assert True


@pytest.fixture
def Ku():
    return Ku_model()


config = "permamodel/examples/Ku_example_config.toml"


def test_read_config(Ku):
    Ku.read_config(config)

    assert Ku.experiment == "test"
    assert Ku.inputs_dir == "permamodel/data/test_directory/inputs/"
    assert Ku.outputs_dir == "permamodel/data/test_directory/outputs/"
    assert Ku.number_of_years == 100
    assert Ku.grid_shape == [100, 100]
    assert len(Ku.input_files) > 0
    assert Ku.soils["sand"]["heat_capacity"] == 1500
    assert len(Ku.constants) > 0


class TestReadInputs:
    def test_throw_error_if_no_files_list(self, Ku):
        with pytest.raises(ValueError):
            Ku.read_input_files()

    def test_read_inputs(self, Ku):
        Ku.read_config(config)
        Ku.read_input_files()

        assert Ku.snow_thickness.shape == (100, 100, 100)
        assert Ku.soils["sand"]["fraction"].mean() == 0.25

    def test_initialize(self, Ku):
        Ku.initialize(config)

        assert Ku.snow_thickness.shape == (100, 100, 100)
        assert Ku.soils["sand"]["fraction"].mean() == 0.25
        assert "permafrost_temperature" in Ku.results.keys()
        assert "active_layer_thickness" in Ku.results.keys()


@pytest.fixture
def Kutest():
    K = Ku_model()
    K.read_config(config)
    K.read_input_files()
    return K


class TestUpdatePhysicalProperties:
    def test_bulk_heat_capacity(self, Kutest):
        Kutest.update_soil_heat_capacity(0)

        assert_approx_equal(
            Kutest.bulk_thawed_heat_capacity[0, 0], 1.415e6, significant=4
        )
        assert_approx_equal(
            Kutest.bulk_frozen_heat_capacity[0, 0], 1.414e6, significant=4
        )

    def test_bulk_thermal_conductivity(self, Kutest):
        Kutest.update_soil_thermal_conductivity(0)

        assert_approx_equal(
            Kutest.bulk_thawed_conductivity[0, 0], 0.5880, significant=4
        )
        assert_approx_equal(
            Kutest.bulk_frozen_conductivity[0, 0], 1.0970, significant=4
        )

    def test_update_snow_thermal_properties(self, Kutest):
        Kutest.update_snow_thermal_properties(0)

        assert_approx_equal(
            Kutest.snow_thermal_conductivity[0, 0], 0.08182, significant=4
        )


class TestUpdateSurfaceTemperature:
    def test_update_season_durations(self, Kutest):
        Kutest.update_season_durations(0)

        assert_approx_equal(Kutest.length_of_cold_season[0, 0], 1.951e7, significant=4)
        assert_approx_equal(Kutest.length_of_warm_season[0, 0], 1.205e7, significant=4)

    def test_update_snow_and_veg_insulation(self, Kutest):
        Kutest.update_snow_thermal_properties(0)
        Kutest.update_season_durations(0)
        Kutest.update_ground_surface_temperature(0)

        assert_approx_equal(Kutest.snow_insulation[0, 0], 5.610, significant=4)
        assert_approx_equal(Kutest.snow_damping[0, 0], 3.571, significant=4)
        assert_approx_equal(
            Kutest.temperature_at_vegetation[0, 0], -1.297, significant=4
        )
        assert_approx_equal(Kutest.amplitude_at_vegetation[0, 0], 15.48, significant=4)

    def test_update_vegetation_effects(self, Kutest):
        Kutest.update_snow_thermal_properties(0)
        Kutest.update_season_durations(0)
        Kutest.update_ground_surface_temperature(0)

        assert_approx_equal(
            Kutest.winter_vegetation_effect[0, 0], 0.3051, significant=4
        )
        assert_approx_equal(
            Kutest.summer_vegetation_effect[0, 0], 0.4896, significant=4
        )
        assert_approx_equal(Kutest.vegetation_damping[0, 0], 0.3756, significant=4)
        assert_approx_equal(Kutest.vegetation_insulation[0, 0], 0.2391, significant=4)
        assert_approx_equal(
            Kutest.ground_surface_temperature[0, 0], -1.058, significant=4
        )
        assert_approx_equal(Kutest.ground_surface_amplitude[0, 0], 15.10, significant=4)


class TestUpdatePermafrostTemperature:
    def test_update_permafrost_temperature(self, Kutest):
        Kutest.update_soil_heat_capacity(0)
        Kutest.update_soil_thermal_conductivity(0)
        Kutest.update_snow_thermal_properties(0)
        Kutest.update_season_durations(0)
        Kutest.update_ground_surface_temperature(0)
        Kutest.update_permafrost_temperature(0)

        assert_approx_equal(Kutest.permafrost_temperature[0, 0], -3.006, significant=4)


class TestUpdateActiveLayer:
    def test_update_active_layer(self, Kutest):
        Kutest.update_soil_heat_capacity(0)
        Kutest.update_soil_thermal_conductivity(0)
        Kutest.update_snow_thermal_properties(0)
        Kutest.update_season_durations(0)
        Kutest.update_ground_surface_temperature(0)
        Kutest.update_permafrost_temperature(0)
        Kutest.update_active_layer(0)

        assert_approx_equal(Kutest.permafrost_amplitude[0, 0], 8.783, significant=4)
        assert_approx_equal(Kutest.critical_depth[0, 0], 0.7504, significant=4)
        assert_approx_equal(Kutest.active_layer_thickness[0, 0], 1.226, significant=4)


class TestRun:
    def test_run_one_step(self, Kutest):
        Kutest.run_one_step(0)

        assert_approx_equal(Kutest.permafrost_amplitude[0, 0], 8.783, significant=4)
        assert_approx_equal(Kutest.critical_depth[0, 0], 0.7504, significant=4)
        assert_approx_equal(Kutest.active_layer_thickness[0, 0], 1.226, significant=4)

    def test_run_all_steps(self, Kutest):
        Kutest.number_of_years = 10
        Kutest.run_all_steps()

        assert_approx_equal(Kutest.permafrost_amplitude[0, 0], 8.995, significant=4)
        assert_approx_equal(Kutest.critical_depth[0, 0], 0.7525, significant=4)
        assert_approx_equal(Kutest.active_layer_thickness[0, 0], 1.224, significant=4)


class TestWriteOutput:
    def test_construct_results(self, Kutest):
        Kutest.construct_results()
        assert "permafrost_temperature" in Kutest.results
        assert "active_layer_thickness" in Kutest.results

        Kutest.run_one_step(0)

        assert_array_equal(
            Kutest.permafrost_temperature, Kutest.results["permafrost_temperature"][0]
        )
        assert_array_equal(
            Kutest.active_layer_thickness, Kutest.results["active_layer_thickness"][0]
        )

    @pytest.fixture
    def output_file(self, tmp_path):
        target_output = os.path.join(tmp_path, "output.nc")
        return target_output

    def test_write_output(self, Kutest, output_file):
        Kutest.construct_results(vars_to_save=["permafrost_temperature"])
        Kutest.run_one_step(0)
        Kutest.write_output(output_file, ["permafrost_temperature"])

        data = Dataset(output_file)

        assert "permafrost_temperature" in data.variables.keys()
        assert data["permafrost_temperature"].shape == (100, 100, 100)
