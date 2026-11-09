"""nh3_props.py - pure NH3: Peng-Robinson EOS + Lee-Kesler/Pitzer vapour pressure + ideal-gas Cp0 (B&S Table A.6).
Reference: T_ref = 273.15 K, P_ref = 101.325 kPa, ideal gas, h_ref = s_ref = 0.  Units: K, kPa, kJ/kg, kJ/kg-K."""
import numpy as np

CRITICAL_TEMPERATURE_K, CRITICAL_PRESSURE_KPA, ACENTRIC_FACTOR, MOLECULAR_WEIGHT_KG_KMOL = 405.56, 11333.0, 0.256, 17.031
SPECIFIC_GAS_CONSTANT_KJ_KG_K = 8.314462 / MOLECULAR_WEIGHT_KG_KMOL                                
REFERENCE_TEMPERATURE_K, REFERENCE_PRESSURE_KPA = 273.15, 101.325
IDEAL_GAS_CP_COEFFICIENTS = (4.238, -4.215, 2.041, -0.2126)                                                                       

def cp0(temperature_k):
    reduced_temperature = temperature_k / 1000.0
    return IDEAL_GAS_CP_COEFFICIENTS[0] + IDEAL_GAS_CP_COEFFICIENTS[1] * reduced_temperature + IDEAL_GAS_CP_COEFFICIENTS[2] * reduced_temperature**2 + IDEAL_GAS_CP_COEFFICIENTS[3] * reduced_temperature**3

def _h0(temperature_k):                                  
    enthalpy_primitive = lambda T: 1000 * (IDEAL_GAS_CP_COEFFICIENTS[0] * (T/1000) + IDEAL_GAS_CP_COEFFICIENTS[1] * (T/1000)**2 / 2 + IDEAL_GAS_CP_COEFFICIENTS[2] * (T/1000)**3 / 3 + IDEAL_GAS_CP_COEFFICIENTS[3] * (T/1000)**4 / 4)
    return enthalpy_primitive(temperature_k) - enthalpy_primitive(REFERENCE_TEMPERATURE_K)

def _s0(temperature_k):                                    
    entropy_primitive = lambda T: IDEAL_GAS_CP_COEFFICIENTS[0] * np.log(T) + IDEAL_GAS_CP_COEFFICIENTS[1] * (T/1000) + IDEAL_GAS_CP_COEFFICIENTS[2] * (T/1000)**2 / 2 + IDEAL_GAS_CP_COEFFICIENTS[3] * (T/1000)**3 / 3
    return entropy_primitive(temperature_k) - entropy_primitive(REFERENCE_TEMPERATURE_K)

def psat(temperature_k):
    reduced_temperature = temperature_k / CRITICAL_TEMPERATURE_K
    simple_fluid_pressure_term = 5.92714 - 6.09648 / reduced_temperature - 1.28862 * np.log(reduced_temperature) + 0.169347 * reduced_temperature**6
    acentric_factor_pressure_term = 15.2518 - 15.6875 / reduced_temperature - 13.4721 * np.log(reduced_temperature) + 0.43577 * reduced_temperature**6
    return CRITICAL_PRESSURE_KPA * np.exp(simple_fluid_pressure_term + ACENTRIC_FACTOR * acentric_factor_pressure_term)

                         
_kap = 0.37464 + 1.54226 * ACENTRIC_FACTOR - 0.26992 * ACENTRIC_FACTOR**2
_ac = 0.45724 * SPECIFIC_GAS_CONSTANT_KJ_KG_K**2 * CRITICAL_TEMPERATURE_K**2 / CRITICAL_PRESSURE_KPA
_b = 0.07780 * SPECIFIC_GAS_CONSTANT_KJ_KG_K * CRITICAL_TEMPERATURE_K / CRITICAL_PRESSURE_KPA
def _a(temperature_k):  return _ac * (1 + _kap * (1 - np.sqrt(temperature_k / CRITICAL_TEMPERATURE_K)))**2
def _da(temperature_k): return -_ac * _kap * (1 + _kap * (1 - np.sqrt(temperature_k / CRITICAL_TEMPERATURE_K))) / np.sqrt(temperature_k * CRITICAL_TEMPERATURE_K)

def Z_roots(temperature_k, pressure_kpa):
    eos_parameter_a = _a(temperature_k) * pressure_kpa / (SPECIFIC_GAS_CONSTANT_KJ_KG_K * temperature_k)**2; eos_parameter_b = _b * pressure_kpa / (SPECIFIC_GAS_CONSTANT_KJ_KG_K * temperature_k)
    cubic_coefficients = [1, -(1 - eos_parameter_b), eos_parameter_a - 3 * eos_parameter_b**2 - 2 * eos_parameter_b, -(eos_parameter_a * eos_parameter_b - eos_parameter_b**2 - eos_parameter_b**3)]
    compressibility_roots = np.roots(cubic_coefficients); compressibility_roots = np.sort(compressibility_roots[abs(compressibility_roots.imag) < 1e-9].real); compressibility_roots = compressibility_roots[compressibility_roots > eos_parameter_b]
    return compressibility_roots, eos_parameter_a, eos_parameter_b

def _dep(temperature_k, pressure_kpa, compressibility_factor, eos_parameter_b):
    logarithmic_eos_term = np.log((compressibility_factor + (1 + np.sqrt(2)) * eos_parameter_b) / (compressibility_factor + (1 - np.sqrt(2)) * eos_parameter_b))
    enthalpy_departure_kj_kg = SPECIFIC_GAS_CONSTANT_KJ_KG_K * temperature_k * (compressibility_factor - 1) + (temperature_k * _da(temperature_k) - _a(temperature_k)) / (2 * np.sqrt(2) * _b) * logarithmic_eos_term
    entropy_departure_kj_kg_k = SPECIFIC_GAS_CONSTANT_KJ_KG_K * np.log(compressibility_factor - eos_parameter_b) + _da(temperature_k) / (2 * np.sqrt(2) * _b) * logarithmic_eos_term
    return enthalpy_departure_kj_kg, entropy_departure_kj_kg_k

def props(temperature_k, pressure_kpa, phase_name='vap'):
    compressibility_roots, eos_parameter_a, eos_parameter_b = Z_roots(temperature_k, pressure_kpa)
    compressibility_factor = compressibility_roots[-1] if phase_name == 'vap' else compressibility_roots[0]
    enthalpy_departure_kj_kg, entropy_departure_kj_kg_k = _dep(temperature_k, pressure_kpa, compressibility_factor, eos_parameter_b)
    specific_enthalpy_kj_kg = _h0(temperature_k) + enthalpy_departure_kj_kg
    specific_entropy_kj_kg_k = _s0(temperature_k) - SPECIFIC_GAS_CONSTANT_KJ_KG_K * np.log(pressure_kpa / REFERENCE_PRESSURE_KPA) + entropy_departure_kj_kg_k
    return dict(h=specific_enthalpy_kj_kg, s=specific_entropy_kj_kg_k, v=compressibility_factor * SPECIFIC_GAS_CONSTANT_KJ_KG_K * temperature_k / pressure_kpa, Z=compressibility_factor)

def sat(temperature_k):
    saturation_pressure_kpa = psat(temperature_k); saturated_liquid_properties = props(temperature_k, saturation_pressure_kpa, 'liq'); saturated_vapour_properties = props(temperature_k, saturation_pressure_kpa, 'vap')
    return dict(P=saturation_pressure_kpa, hf=saturated_liquid_properties['h'], hg=saturated_vapour_properties['h'], hfg=saturated_vapour_properties['h'] - saturated_liquid_properties['h'], sf=saturated_liquid_properties['s'], sg=saturated_vapour_properties['s'], vf=saturated_liquid_properties['v'], vg=saturated_vapour_properties['v'], Zg=saturated_vapour_properties['Z'], Zf=saturated_liquid_properties['Z'])

def h_vap(temperature_k, pressure_kpa): return props(temperature_k, pressure_kpa, 'vap')['h']
def h_liq_sat(temperature_k): return sat(temperature_k)['hf']
def T_sat(pressure_kpa):
    from scipy.optimize import brentq
    return brentq(lambda T: psat(T) - pressure_kpa, 200.0, CRITICAL_TEMPERATURE_K - 1)

if __name__ == '__main__':
    tab = {0: (429.6, 1261.8), 40: (1554.9, 1098.8)}
    print('T(C)  Psat   tab    %err  | hfg    tab    %err  | Zg')
    for Tc_ in (0, 40):
        s = sat(Tc_ + 273.15); P0, h0 = tab[Tc_]
        print('%3d %8.1f %7.1f %5.2f | %7.1f %7.1f %5.2f | %.3f' % (Tc_, s['P'], P0, 100 * (s['P'] / P0 - 1), s['hfg'], h0, 100 * (s['hfg'] / h0 - 1), s['Zg']))
    g0, f40 = sat(273.15), sat(313.15)
    qe = g0['hg'] - f40['hf']; print('q_e = %.1f kJ/kg ; m_r for 5000 kW = %.3f kg/s' % (qe, 5000 / qe))
