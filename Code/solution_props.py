"""solution_props.py - adapter exposing the Patek & Klomfar (1995) module with the interface cycle.py expects
(same function names/units as the old nh3h2o_fallback): P in kPa, T in K, *_mole = NH3 mole fraction, *_m / x_mass = NH3 MASS fraction.
All solution and pure-refrigerant enthalpies use the same Patek–Klomfar (P&K) correlation family and are shifted by constant per-species offsets. The ammonia offset is anchored to the Peng–Robinson saturated-vapour enthalpy at 273.15 K; this is a practical cross-reference anchor, not a claim that the saturated vapour enthalpy is zero. The water offset uses the stated liquid-water reference. Constant species offsets cancel in closed cycle balances when the species balances close."""
import numpy as np
from scipy.optimize import brentq
import equilibrium as equilibrium_model
import nh3_props as ammonia_properties

PROPERTY_MODEL_DESCRIPTION = 'Patek & Klomfar (1995) Eqs. 6, 7, 9, 10 (own implementation)'
MOLECULAR_WEIGHT_AMMONIA_KG_KMOL, MOLECULAR_WEIGHT_WATER_KG_KMOL = equilibrium_model.MOLECULAR_WEIGHT_AMMONIA_KG_KMOL, equilibrium_model.MOLECULAR_WEIGHT_WATER_KG_KMOL
                                                                              
                                                                            
                                                                     
                                                                                 
                                                                               
                                                        
equilibrium_model.set_ref_shift_nh3(ammonia_properties.sat(273.15)['hg'] - float(equilibrium_model.h_gas(273.15, 1.0)))

mole_fraction_to_mass_fraction = lambda x: float(equilibrium_model.mole_to_mass(x))
mass_fraction_to_mole_fraction = lambda xm: float(equilibrium_model.mass_to_mole(xm))
def MW_mix(ammonia_mole_fraction): return ammonia_mole_fraction * MOLECULAR_WEIGHT_AMMONIA_KG_KMOL + (1 - ammonia_mole_fraction) * MOLECULAR_WEIGHT_WATER_KG_KMOL

def T_bub(pressure_kpa, liquid_ammonia_mole_fraction):  return float(equilibrium_model.T_bub(pressure_kpa / 1e3, liquid_ammonia_mole_fraction))
def T_dew(pressure_kpa, vapour_ammonia_mole_fraction):  return float(equilibrium_model.T_dew(pressure_kpa / 1e3, vapour_ammonia_mole_fraction))
def x_bub(pressure_kpa, temperature_k):  return float(equilibrium_model.x_bub(temperature_k, pressure_kpa / 1e3))
def y_dew(pressure_kpa, temperature_k):  return float(equilibrium_model.y_dew(temperature_k, pressure_kpa / 1e3))
def y_eq(pressure_kpa, liquid_ammonia_mole_fraction):   return float(equilibrium_model.y_of_px(pressure_kpa / 1e3, liquid_ammonia_mole_fraction))

def h_NH3_liq(temperature_k):  return h_liquid(0.0, 1.0, temperature_k)                                                                     
def h_NH3_vap(temperature_k):  return h_vapour(0.0, 1.0, temperature_k)                                                                     

def h_liquid(pressure_kpa, liquid_ammonia_mass_fraction, temperature_k):
    return float(equilibrium_model.to_project_reference(equilibrium_model.h_liq_mass(temperature_k, liquid_ammonia_mass_fraction), liquid_ammonia_mass_fraction))
def h_vapour(pressure_kpa, vapour_ammonia_mass_fraction, temperature_k):
    return float(equilibrium_model.to_project_reference(equilibrium_model.h_gas_mass(temperature_k, vapour_ammonia_mass_fraction), vapour_ammonia_mass_fraction))
def T_from_hL(pressure_kpa, liquid_ammonia_mass_fraction, target_enthalpy_kj_kg, minimum_temperature_k=255.0, maximum_temperature_k=445.0):
    return brentq(lambda T: h_liquid(pressure_kpa, liquid_ammonia_mass_fraction, T) - target_enthalpy_kj_kg, minimum_temperature_k, maximum_temperature_k)
def rho_liquid(temperature_k, liquid_ammonia_mass_fraction):
    ammonia_density_kg_m3 = 639.0 - 1.5 * (temperature_k - 273.15); water_density_kg_m3 = 1000.0 - 0.25 * max(temperature_k - 293.15, 0)
    return 1.0 / (liquid_ammonia_mass_fraction / ammonia_density_kg_m3 + (1 - liquid_ammonia_mass_fraction) / water_density_kg_m3)
