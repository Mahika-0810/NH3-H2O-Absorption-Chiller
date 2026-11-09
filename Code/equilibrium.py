"""
equilibrium.py
================
Own implementation of the Patek & Klomfar (1995) NH3-H2O equilibrium
correlations (Int. J. Refrig. 18(4) 228-234), Eqs. (6)-(10):

    T(p,x)   bubble temperature         Eq. (6)
    T(p,y)   dew temperature            Eq. (7)
    y(p,x)   vapour mole fraction       Eq. (8)
    h_l(T,x) liquid enthalpy            Eq. (9)
    h_g(T,y) vapour enthalpy            Eq. (10)

plus inverse helpers needed by the cycle solver (all root-finding based on the
explicit equations above, no external property library).

Units (as in the paper)
-----------------------
    p  : MPa          T : K          h : kJ/kg
    x  : NH3 mole fraction, liquid   y : NH3 mole fraction, vapour
    xi : NH3 MASS fraction (helpers ending in *_mass)

Reference state of the raw correlations
---------------------------------------
    Liquid NH3 : h = 0 at 273.16 K  (verified: Eq. 9 at x=1,T=273.16 gives ~0)
    Liquid H2O : h = 0 at the triple point
Conversion to the project reference (ideal gas, 273.15 K, 101.325 kPa, h=0)
is done in `to_project_reference()` with per-species shifts - see
`REF_SHIFT_*` and the docstring there.

!!! COEFFICIENT NOTE (handout typos, verified numerically) !!!
The tables printed in the assignment handout contain three transcription
errors with respect to the original paper. Each was identified by checking
the correlation against known saturation data (see tests/test_equilibrium.py):
    Eq.(6) a13 : handout 4.895960E-02  ->  paper 4.89596 (x10^1 in print).
                 With E-02, pure water at 2 MPa boils at 1.6 K (nonsense);
                 with 4.89596 it boils at 486.3 K (true 485.5 K).
    Eq.(7) a3  : handout -4.356240E-02 ->  paper +4.35624E-02.
                 Pure NH3 at 0.1 MPa: +sign 238.6 K (true 239.7 K),
                 -sign 160 K.
    Eq.(7) a10 : handout -2.107800E+01 ->  paper -2.01780E+01.
                 Pure water vapour at 0.1 MPa: 373.0 K (true 372.8 K);
                 handout value gives 283 K.
The corrected values are used below.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq

                                                                            
           
                                                                            
MOLECULAR_WEIGHT_AMMONIA_KG_KMOL = 17.03026                        
MOLECULAR_WEIGHT_WATER_KG_KMOL = 18.01528                        
CORRELATION_REFERENCE_PRESSURE_MPA = 2.0                                 

                                                                            
                                     
                                                                            
                    
BUBBLE_TEMPERATURE_COEFFICIENTS = np.array([
    (0, 0, 3.223020E+00), (0, 1, -3.842060E-01), (0, 2, 4.609650E-02),
    (0, 3, -3.789450E-03), (0, 4, 1.356100E-04), (1, 0, 4.877550E-01),
    (1, 1, -1.201080E-01), (1, 2, 1.061540E-02), (2, 3, -5.335890E-04),
    (4, 0, 7.850410E+00), (5, 0, -1.159410E+01), (5, 1, -5.231500E-02),
    (6, 0, 4.895960E+00),                                                 
    (13, 1, 4.210590E-02)])
                    
DEW_TEMPERATURE_COEFFICIENTS = np.array([
    (0, 0, 3.240040E+00), (0, 1, -3.959200E-01),
    (0, 2, 4.356240E-02),                                      
    (0, 3, -2.189430E-03), (1, 0, -1.435260E+00), (1, 1, 1.052560E+00),
    (1, 2, -7.192810E-02), (2, 0, 1.223620E+01), (2, 1, -2.243680E+00),
    (3, 0, -2.017800E+01),                                                  
    (3, 1, 1.108340E+00), (4, 0, 1.453990E+01), (4, 2, 6.443120E-01),
    (5, 0, -2.212640E+00), (5, 2, -7.562660E-01), (6, 0, -1.355290E+00),
    (7, 2, 1.835410E-01)])
        
VAPOUR_COMPOSITION_COEFFICIENTS = np.array([
    (0, 0, 1.98022017E+01), (0, 1, -1.18092669E+01), (0, 6, 2.77479980E+01),
    (0, 7, -2.88634277E+01), (1, 0, -5.91616608E+01), (2, 1, 5.78091305E+02),
    (2, 2, -6.21736743E+00), (3, 2, -3.42198402E+03), (4, 3, 1.19403127E+04),
    (5, 4, -2.45413777E+04), (6, 5, 2.91591865E+04), (7, 6, -1.84782290E+04),
    (7, 7, 2.34819434E+01),                                           
    (8, 7, 4.80310617E+03)])
                                       
LIQUID_ENTHALPY_COEFFICIENTS = np.array([
    (0, 1, -7.61080E+00), (0, 4, 2.56905E+01), (0, 8, -2.47092E+02),
    (0, 9, 3.25952E+02), (0, 12, -1.58854E+02), (0, 14, 6.19084E+01),
    (1, 0, 1.14314E+01), (1, 1, 1.18157E+00), (2, 1, 2.84179E+00),
    (3, 3, 7.41609E+00), (5, 3, 8.91844E+02), (5, 4, -1.61309E+03),
    (5, 5, 6.22106E+02), (6, 2, -2.07588E+02), (6, 4, -6.87393E+00),
    (8, 0, 3.50716E+00)])
                                      
VAPOUR_ENTHALPY_COEFFICIENTS = np.array([
    (0, 0, 1.28827E+00), (1, 0, 1.25247E-01), (2, 0, -2.08748E+00),
    (3, 0, 2.17696E+00), (0, 2, 2.35687E+00), (1, 2, -8.86987E+00),
    (2, 2, 1.02635E+01), (3, 2, -2.37440E+00), (0, 3, -6.70515E+00),
    (1, 3, 1.64508E+01), (2, 3, -9.36849E+00), (0, 4, 8.42254E+00),
    (1, 4, -8.58807E+00), (0, 5, -2.77049E+00), (4, 6, -9.61248E-01),
    (2, 7, 9.88009E-01), (1, 10, 3.08482E-01)])


def _cols(coefficient_table):
    return coefficient_table[:, 0], coefficient_table[:, 1], coefficient_table[:, 2]


                                                                            
                         
                                                                            
def mole_to_mass(ammonia_mole_fraction):
    ammonia_mole_fraction = np.asarray(ammonia_mole_fraction, float)
    return ammonia_mole_fraction * MOLECULAR_WEIGHT_AMMONIA_KG_KMOL / (ammonia_mole_fraction * MOLECULAR_WEIGHT_AMMONIA_KG_KMOL + (1.0 - ammonia_mole_fraction) * MOLECULAR_WEIGHT_WATER_KG_KMOL)


def mass_to_mole(ammonia_mass_fraction):
    ammonia_mass_fraction = np.asarray(ammonia_mass_fraction, float)
    return (ammonia_mass_fraction / MOLECULAR_WEIGHT_AMMONIA_KG_KMOL) / (ammonia_mass_fraction / MOLECULAR_WEIGHT_AMMONIA_KG_KMOL + (1.0 - ammonia_mass_fraction) / MOLECULAR_WEIGHT_WATER_KG_KMOL)


                                                                            
                                
                                                                            
def T_bub(pressure_mpa, liquid_ammonia_mole_fraction):
    pressure_exponents, composition_exponents, correlation_coefficients = _cols(BUBBLE_TEMPERATURE_COEFFICIENTS)
    log_pressure_ratio = np.log(CORRELATION_REFERENCE_PRESSURE_MPA / np.asarray(pressure_mpa, float))
    liquid_ammonia_mole_fraction = np.asarray(liquid_ammonia_mole_fraction, float)
    return 100.0 * sum(coefficient * (1.0 - liquid_ammonia_mole_fraction) ** pressure_exponent * log_pressure_ratio ** composition_exponent
                       for pressure_exponent, composition_exponent, coefficient in zip(pressure_exponents, composition_exponents, correlation_coefficients))


def T_dew(pressure_mpa, vapour_ammonia_mole_fraction):
    pressure_exponents, composition_exponents, correlation_coefficients = _cols(DEW_TEMPERATURE_COEFFICIENTS)
    log_pressure_ratio = np.log(CORRELATION_REFERENCE_PRESSURE_MPA / np.asarray(pressure_mpa, float))
    vapour_ammonia_mole_fraction = np.asarray(vapour_ammonia_mole_fraction, float)
    return 100.0 * sum(coefficient * (1.0 - vapour_ammonia_mole_fraction) ** (pressure_exponent / 4.0) * log_pressure_ratio ** composition_exponent
                       for pressure_exponent, composition_exponent, coefficient in zip(pressure_exponents, composition_exponents, correlation_coefficients))


def y_of_px_eq8(pressure_mpa, liquid_ammonia_mole_fraction):
    pressure_exponents, composition_exponents, correlation_coefficients = _cols(VAPOUR_COMPOSITION_COEFFICIENTS)
    reduced_pressure = np.asarray(pressure_mpa, float) / CORRELATION_REFERENCE_PRESSURE_MPA
    liquid_ammonia_mole_fraction = np.asarray(liquid_ammonia_mole_fraction, float)
    composition_exponent_sum = sum(coefficient * reduced_pressure ** pressure_exponent * liquid_ammonia_mole_fraction ** (composition_exponent / 3.0) for pressure_exponent, composition_exponent, coefficient in zip(pressure_exponents, composition_exponents, correlation_coefficients))
    return 1.0 - np.exp(np.log(1.0 - liquid_ammonia_mole_fraction) * composition_exponent_sum)


def y_of_px(pressure_mpa, liquid_ammonia_mole_fraction):
    bubble_temperature_k = float(T_bub(pressure_mpa, liquid_ammonia_mole_fraction))
    return y_dew(bubble_temperature_k, float(pressure_mpa))


def h_liq(temperature_k, liquid_ammonia_mole_fraction):
    temperature_exponents, composition_exponents, correlation_coefficients = _cols(LIQUID_ENTHALPY_COEFFICIENTS)
    reduced_temperature = np.asarray(temperature_k, float) / 273.16 - 1.0
    liquid_ammonia_mole_fraction = np.asarray(liquid_ammonia_mole_fraction, float)
    return 100.0 * sum(coefficient * reduced_temperature ** temperature_exponent * liquid_ammonia_mole_fraction ** composition_exponent for temperature_exponent, composition_exponent, coefficient in zip(temperature_exponents, composition_exponents, correlation_coefficients))


def h_gas(temperature_k, vapour_ammonia_mole_fraction):
    temperature_exponents, composition_exponents, correlation_coefficients = _cols(VAPOUR_ENTHALPY_COEFFICIENTS)
    reduced_temperature = 1.0 - np.asarray(temperature_k, float) / 324.0
    vapour_ammonia_mole_fraction = np.asarray(vapour_ammonia_mole_fraction, float)
    return 1000.0 * sum(coefficient * reduced_temperature ** temperature_exponent * (1.0 - vapour_ammonia_mole_fraction) ** (composition_exponent / 4.0)
                        for temperature_exponent, composition_exponent, coefficient in zip(temperature_exponents, composition_exponents, correlation_coefficients))


                                                                            
                                                  
                                                                            
def x_bub(temperature_k, pressure_mpa, root_tolerance=1e-12):
    bubble_temperature_residual = lambda x: T_bub(pressure_mpa, x) - temperature_k
    lower_mole_fraction_bound, upper_mole_fraction_bound = 0.0, 1.0
    lower_bound_residual, upper_bound_residual = bubble_temperature_residual(lower_mole_fraction_bound), bubble_temperature_residual(upper_mole_fraction_bound)
    if lower_bound_residual * upper_bound_residual > 0:
        raise ValueError(f"x_bub: T={temperature_k:.2f} K outside [{T_bub(pressure_mpa,0):.1f}, "
                         f"{T_bub(pressure_mpa,1):.1f}] K at p={pressure_mpa:.4f} MPa")
    return brentq(bubble_temperature_residual, lower_mole_fraction_bound, upper_mole_fraction_bound, xtol=root_tolerance)


def p_bub(temperature_k, liquid_ammonia_mole_fraction, minimum_pressure_mpa=0.005, maximum_pressure_mpa=3.0):
    temperature_residual = lambda lp: T_bub(np.exp(lp), liquid_ammonia_mole_fraction) - temperature_k
    log_pressure_mpa = brentq(temperature_residual, np.log(minimum_pressure_mpa), np.log(maximum_pressure_mpa), xtol=1e-13)
    return float(np.exp(log_pressure_mpa))


def y_dew(temperature_k, pressure_mpa, root_tolerance=1e-12):
    dew_temperature_residual = lambda y: T_dew(pressure_mpa, y) - temperature_k
    return brentq(dew_temperature_residual, 0.0, 1.0, xtol=root_tolerance)


def mass_balance_split(strong_solution_ammonia_mass_fraction, weak_solution_ammonia_mass_fraction, refrigerant_ammonia_mass_fraction):
    solution_circulation_ratio = (refrigerant_ammonia_mass_fraction - weak_solution_ammonia_mass_fraction) / (strong_solution_ammonia_mass_fraction - weak_solution_ammonia_mass_fraction)
    return solution_circulation_ratio, solution_circulation_ratio - 1.0                                          


                                                                            
                                                             
                                                                            
def h_liq_mass(temperature_k, ammonia_mass_fraction):
    return h_liq(temperature_k, mass_to_mole(ammonia_mass_fraction))


def h_gas_mass(temperature_k, vapour_ammonia_mass_fraction):
    return h_gas(temperature_k, mass_to_mole(vapour_ammonia_mass_fraction))


def sat_liquid_xi(temperature_k, pressure_mpa):
    return float(mole_to_mass(x_bub(temperature_k, pressure_mpa)))


def vapour_xi_over_liquid(temperature_k, pressure_mpa):
    return float(mole_to_mass(y_dew(temperature_k, pressure_mpa)))


                                                                            
                                                         
                                                                            
                                                          
                                                                              
                                                                           
                                                                             
                                   
                                                                                 
                                                                                 
                                                                               
                                                                               
                                             
REF_SHIFT_NH3 = -1281.5
REF_SHIFT_H2O = -2500.9


def set_ref_shift_nh3(enthalpy_reference_shift_kj_kg):
    global REF_SHIFT_NH3
    REF_SHIFT_NH3 = float(enthalpy_reference_shift_kj_kg)


def to_project_reference(raw_correlation_enthalpy_kj_kg, ammonia_mass_fraction):
    return raw_correlation_enthalpy_kj_kg + ammonia_mass_fraction * REF_SHIFT_NH3 + (1.0 - ammonia_mass_fraction) * REF_SHIFT_H2O
