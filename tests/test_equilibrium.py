"""Regression tests for the direct Patek-Klomfar Eq. (8) correlation."""
import pytest
import equilibrium


@pytest.mark.parametrize("pressure_mpa", [0.1, 0.5, 1.0, 1.5, 1.9])
@pytest.mark.parametrize("liquid_ammonia_mole_fraction", [0.1, 0.2, 0.33, 0.5, 0.8])
def test_eq8_remains_physical_and_agrees_with_eq6_eq7(pressure_mpa, liquid_ammonia_mole_fraction):
    y_eq8 = float(equilibrium.y_of_px_eq8(pressure_mpa, liquid_ammonia_mole_fraction))
    y_eq67 = float(equilibrium.y_of_px(pressure_mpa, liquid_ammonia_mole_fraction))
    assert 0.0 < y_eq8 < 1.0
    assert abs(y_eq8 - y_eq67) / y_eq67 <= 0.055


# Property anchors, round-trip checks, and physical trends.
def close(actual_value, expected_value, allowed_difference):
    assert abs(actual_value - expected_value) <= allowed_difference, f"{actual_value} vs {expected_value} (tol {allowed_difference})"

def test_pure_water_boiling():                                      
    close(equilibrium.T_bub(0.1, 0.0), 372.76, 2.5)
    close(equilibrium.T_bub(1.0, 0.0), 453.03, 2.5)
    close(equilibrium.T_bub(2.0, 0.0), 485.53, 2.5)

def test_pure_ammonia_boiling():                                               
    close(equilibrium.T_bub(0.1013, 1.0), 239.8, 2.0)
    close(equilibrium.T_bub(0.4294, 1.0), 273.15, 2.0)                        
    close(equilibrium.T_bub(1.5549, 1.0), 313.15, 2.0)                          

def test_dew_pure_limits():                            
    close(equilibrium.T_dew(0.1, 0.0), 372.76, 2.5)
    close(equilibrium.T_dew(0.4294, 1.0), 273.15, 2.0)
    close(equilibrium.T_dew(1.5549, 1.0), 313.15, 2.0)

def test_dew_above_bubble():                                   
    for pressure_mpa in (0.5, 1.0, 1.5):
        for liquid_ammonia_mole_fraction in (0.2, 0.4, 0.6):
            vapour_ammonia_mole_fraction = equilibrium.y_of_px(pressure_mpa, liquid_ammonia_mole_fraction)
            assert vapour_ammonia_mole_fraction > liquid_ammonia_mole_fraction
            close(equilibrium.T_dew(pressure_mpa, vapour_ammonia_mole_fraction), equilibrium.T_bub(pressure_mpa, liquid_ammonia_mole_fraction), 1e-6)

def test_liquid_enthalpy_anchors():
    close(equilibrium.h_liq(373.15, 0.0), 419.1, 6)                     
    close(equilibrium.h_liq(273.16, 1.0), 0.0, 1.0)                                 
    close(equilibrium.h_liq(273.16, 0.0), 0.0, 1.0)

def test_gas_enthalpy_pure_nh3():
                                                 
    close(equilibrium.h_gas(273.15, 1.0) - equilibrium.h_liq(273.15, 1.0), 1262.0, 25.0)

def test_roundtrips():
    pressure_mpa = 1.2; liquid_ammonia_mole_fraction = 0.35
    bubble_temperature_k = equilibrium.T_bub(pressure_mpa, liquid_ammonia_mole_fraction)
    close(equilibrium.x_bub(bubble_temperature_k, pressure_mpa), liquid_ammonia_mole_fraction, 1e-9)
    close(equilibrium.p_bub(bubble_temperature_k, liquid_ammonia_mole_fraction), pressure_mpa, 1e-6)
    vapour_mole_fraction = equilibrium.y_of_px(pressure_mpa, liquid_ammonia_mole_fraction)
    close(equilibrium.y_dew(equilibrium.T_dew(pressure_mpa, vapour_mole_fraction), pressure_mpa), vapour_mole_fraction, 1e-9)
    ammonia_mass_fraction = 0.4
    close(equilibrium.mole_to_mass(equilibrium.mass_to_mole(ammonia_mass_fraction)), ammonia_mass_fraction, 1e-12)

def test_physical_trends():
                                                                     
    assert equilibrium.T_bub(1.5, 0.3) > equilibrium.T_bub(0.5, 0.3)
    assert equilibrium.T_bub(1.0, 0.5) < equilibrium.T_bub(1.0, 0.2)
                                              
    assert equilibrium.h_liq(350, 0.3) > equilibrium.h_liq(300, 0.3)

