"""Regression tests for invalid cycle inputs and the computed COP."""
import math
import pytest

import cycle as cycle


@pytest.mark.parametrize(
    "invalid_cycle_inputs, expected_error_message",
    [
        ({"evaporator_cooling_kw": 0, "condenser_temperature_k": 313.15}, "evaporator_cooling_kw must be greater than 0"),
        ({"evaporator_cooling_kw": -10, "condenser_temperature_k": 313.15}, "evaporator_cooling_kw must be greater than 0"),
        ({"evaporator_cooling_kw": 5000, "condenser_temperature_k": 270.0}, "must be higher than evaporator_temperature_k"),
        ({"evaporator_cooling_kw": 5000, "condenser_temperature_k": 313.15, "evaporator_temperature_k": 313.15}, "must be higher than evaporator_temperature_k"),
        ({"evaporator_cooling_kw": 5000, "condenser_temperature_k": 313.15, "concentration_difference": 0}, "concentration_difference must be between 0 and 1"),
        ({"evaporator_cooling_kw": 5000, "condenser_temperature_k": 313.15, "concentration_difference": 1.2}, "concentration_difference must be between 0 and 1"),
    ],
)
def test_invalid_inputs_fail_with_clear_message(invalid_cycle_inputs, expected_error_message):
    with pytest.raises(ValueError, match=expected_error_message):
        cycle.solve_cycle(**invalid_cycle_inputs)


def test_cop_is_calculated_from_cycle_energy_inputs():
    cycle_result = cycle.solve_cycle(5000.0, 313.15)
    expected_coefficient_of_performance = cycle_result["Qevap"] / (cycle_result["Qgen"] + cycle_result["Wp"])
    assert math.isclose(cycle_result["COP"], expected_coefficient_of_performance, rel_tol=1e-12)
    assert cycle_result["Qgen"] > 0
    assert cycle_result["Wp"] >= 0
