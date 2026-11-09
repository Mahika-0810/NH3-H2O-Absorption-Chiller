"""Regression tests for cycle outputs, component balances, and generator/rectifier closure."""
import pytest
import cycle
import cycle_study as project_method
import solution_props

CASES = [
    (5000.0, 313.15, 0.16, 275.15),
    (5000.0, 318.15, 0.16, 275.15),
    (1000.0, 303.15, 0.12, 275.15),
    (2500.0, 308.15, 0.20, 278.15),
]

@pytest.mark.parametrize("cooling_kw, condenser_k, dx, evaporator_k", CASES)
def test_cycle_outputs_finite_and_positive(cooling_kw, condenser_k, dx, evaporator_k):
    result = cycle.solve_cycle(cooling_kw, condenser_k, concentration_difference=dx, evaporator_temperature_k=evaporator_k)
    for key in ("COP", "Qgen", "Qcond", "Qabs", "Qrect", "Wp", "mfg"):
        assert result[key] > 0 and result[key] == pytest.approx(result[key])
    assert 0 < result["yg"] < 1
    assert result["bal_m"] == pytest.approx(0, abs=1e-8)
    assert result["bal_NH3"] == pytest.approx(0, abs=1e-8)
    assert result["shx_E"] == pytest.approx(0, abs=1e-6)


def test_generator_vapour_equilibrium_uses_strong_solution_basis():
    result = cycle.solve_cycle(5000.0, 313.15, concentration_difference=0.16, evaporator_temperature_k=275.15)
    generator_pressure_kpa = result["Phigh"]
    strong_mole_fraction = solution_props.mass_fraction_to_mole_fraction(result["x_strong"])
    expected_vapour_mass_fraction = solution_props.mole_fraction_to_mass_fraction(solution_props.y_eq(generator_pressure_kpa, strong_mole_fraction))
    assert result["yg"] == pytest.approx(expected_vapour_mass_fraction, rel=1e-10)
    assert result["Ttop"] == pytest.approx(solution_props.T_bub(generator_pressure_kpa, strong_mole_fraction), rel=1e-10)


def test_state_point_table_has_all_ten_states_and_correct_pressure_levels():
    result = cycle.solve_cycle(5000.0, 313.15)
    rows = project_method._state_point_rows(result)
    assert [row["state"] for row in rows] == list(range(1, 11))
    assert rows[0]["pressure_kPa"] == pytest.approx(result["Phigh"])
    assert rows[2]["pressure_kPa"] == pytest.approx(result["Plow"])
    assert rows[9]["pressure_kPa"] == pytest.approx(result["Plow"])


def test_study_component_balances_close():
    result = cycle.solve_cycle(5000.0, 313.15)
    rows = project_method._component_balance_rows(result)
    assert len(rows) >= 8
    assert all(row["result"] == "PASS" for row in rows)
