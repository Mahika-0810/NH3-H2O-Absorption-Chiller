"""air_handler.py - ideal-gas moist-air relations + one-zone AHU model (Project 2, Tasks 1/3/4 air side).
Reference: 0 C, ideal gas, h = 1.006 T + w (2501 + 1.86 T) [kJ/kg da, T in C].
CoolProp is used ONLY for P_sat,H2O (and h_g water); if CoolProp is absent an own Wagner-Pruss-type P_sat correlation is used (<0.02 % error)."""
import numpy as np
from scipy.optimize import brentq

ATMOSPHERIC_PRESSURE_KPA = 101.325       
try:
    from CoolProp.CoolProp import PropsSI
    def psat(dry_bulb_temperature_c):               
        return PropsSI('P', 'T', dry_bulb_temperature_c + 273.15, 'Q', 0, 'Water') / 1e3
    def hg_w(dry_bulb_temperature_c):                                                            
        return 2501.0 + 1.86 * dry_bulb_temperature_c
    PSAT_SOURCE = 'CoolProp (water only, as permitted)'
except Exception:                                                                                                                
    def psat(dry_bulb_temperature_c):
        temperature_k = dry_bulb_temperature_c + 273.15; reduced_temperature = 1.0 - temperature_k / 647.096
        saturation_coefficients = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)
        saturation_exponents = (1.0, 1.5, 3.0, 3.5, 4.0, 7.5)
        return 22064.0 * np.exp(647.096 / temperature_k * sum(coefficient * reduced_temperature ** exponent for coefficient, exponent in zip(saturation_coefficients, saturation_exponents)))        
    PSAT_SOURCE = 'own Wagner-Pruss-type correlation (CoolProp not installed)'
    def hg_w(dry_bulb_temperature_c):
        return 2501.0 + 1.86 * dry_bulb_temperature_c

def w_from_pv(water_vapour_partial_pressure_kpa, atmospheric_pressure_kpa=ATMOSPHERIC_PRESSURE_KPA):  return 0.622 * water_vapour_partial_pressure_kpa / (atmospheric_pressure_kpa - water_vapour_partial_pressure_kpa)
def pv_from_w(humidity_ratio_kg_per_kg_dry_air, atmospheric_pressure_kpa=ATMOSPHERIC_PRESSURE_KPA):   return humidity_ratio_kg_per_kg_dry_air * atmospheric_pressure_kpa / (0.622 + humidity_ratio_kg_per_kg_dry_air)
def w_from_RH(dry_bulb_temperature_c, relative_humidity_fraction, atmospheric_pressure_kpa=ATMOSPHERIC_PRESSURE_KPA): return w_from_pv(relative_humidity_fraction * psat(dry_bulb_temperature_c), atmospheric_pressure_kpa)
def RH_from_w(dry_bulb_temperature_c, humidity_ratio_kg_per_kg_dry_air, atmospheric_pressure_kpa=ATMOSPHERIC_PRESSURE_KPA):  return pv_from_w(humidity_ratio_kg_per_kg_dry_air, atmospheric_pressure_kpa) / psat(dry_bulb_temperature_c)
def wsat(dry_bulb_temperature_c, atmospheric_pressure_kpa=ATMOSPHERIC_PRESSURE_KPA):          return w_from_pv(psat(dry_bulb_temperature_c), atmospheric_pressure_kpa)
def h_moist(dry_bulb_temperature_c, humidity_ratio_kg_per_kg_dry_air):             return 1.006 * dry_bulb_temperature_c + humidity_ratio_kg_per_kg_dry_air * (2501.0 + 1.86 * dry_bulb_temperature_c)
def T_from_hw(moist_air_enthalpy_kj_kg_dry_air, humidity_ratio_kg_per_kg_dry_air):           return (moist_air_enthalpy_kj_kg_dry_air - 2501.0 * humidity_ratio_kg_per_kg_dry_air) / (1.006 + 1.86 * humidity_ratio_kg_per_kg_dry_air)
def hf_w(condensate_temperature_c):                   return 4.186 * condensate_temperature_c                                               

def mix(stream_1_temperature_c, stream_1_humidity_ratio, stream_2_temperature_c, stream_2_humidity_ratio, stream_2_dry_air_mass_fraction):
    mixed_air_humidity_ratio = (1 - stream_2_dry_air_mass_fraction) * stream_1_humidity_ratio + stream_2_dry_air_mass_fraction * stream_2_humidity_ratio
    mixed_air_enthalpy_kj_kg_dry_air = (1 - stream_2_dry_air_mass_fraction) * h_moist(stream_1_temperature_c, stream_1_humidity_ratio) + stream_2_dry_air_mass_fraction * h_moist(stream_2_temperature_c, stream_2_humidity_ratio)
    return T_from_hw(mixed_air_enthalpy_kj_kg_dry_air, mixed_air_humidity_ratio), mixed_air_humidity_ratio, mixed_air_enthalpy_kj_kg_dry_air

def Twb(dry_bulb_temperature_c, humidity_ratio_kg_per_kg_dry_air, atmospheric_pressure_kpa=ATMOSPHERIC_PRESSURE_KPA):
    inlet_air_enthalpy_kj_kg_dry_air = h_moist(dry_bulb_temperature_c, humidity_ratio_kg_per_kg_dry_air)
    saturation_energy_residual = lambda wet_bulb_temperature_c: h_moist(wet_bulb_temperature_c, wsat(wet_bulb_temperature_c, atmospheric_pressure_kpa)) - (wsat(wet_bulb_temperature_c, atmospheric_pressure_kpa) - humidity_ratio_kg_per_kg_dry_air) * hf_w(wet_bulb_temperature_c) - inlet_air_enthalpy_kj_kg_dry_air
    return brentq(saturation_energy_residual, -10, dry_bulb_temperature_c)

                                              
ZONE_MODEL_SETTINGS = dict(envelope_heat_transfer_w_per_k=80.0, solar_heat_gain_w=1200.0, occupant_count=4, sensible_heat_per_occupant_w=75.0, latent_heat_per_occupant_w=55.0, lighting_heat_w_per_m2=20.0, floor_area_m2=40.0,
         dry_air_mass_flow_kg_s=0.40, coil_bypass_factor=0.15, room_temperature_c=22.0, maximum_room_relative_humidity=0.75, zone_count=400, distribution_heat_gain_fraction=0.05, pump_heat_gain_fraction=0.02)
WEATHER_CONDITIONS = {                          
    'Summer day': (40, .50, 1), 'Summer night': (30, .50, 0),
    'Monsoon day': (35, .85, 1), 'Monsoon night': (30, .85, 0),
    'Winter day': (28, .50, 1), 'Winter night': (16, .50, 0)}

def zone(condition_name, outdoor_air_fraction=0.30, forced_apparatus_dew_point_c=None, target_reheat_temperature_c=None, humidity_ratio_tolerance=1e-10):
    outdoor_temperature_c, outdoor_relative_humidity, solar_gain_factor = WEATHER_CONDITIONS[condition_name]; outdoor_humidity_ratio = w_from_RH(outdoor_temperature_c, outdoor_relative_humidity); room_temperature_c = ZONE_MODEL_SETTINGS['room_temperature_c']
    sensible_room_load_w = (ZONE_MODEL_SETTINGS['envelope_heat_transfer_w_per_k'] * (outdoor_temperature_c - room_temperature_c) + solar_gain_factor * ZONE_MODEL_SETTINGS['solar_heat_gain_w'] + ZONE_MODEL_SETTINGS['occupant_count'] * ZONE_MODEL_SETTINGS['sensible_heat_per_occupant_w'] + ZONE_MODEL_SETTINGS['lighting_heat_w_per_m2'] * ZONE_MODEL_SETTINGS['floor_area_m2'])              
    latent_room_load_w = ZONE_MODEL_SETTINGS['occupant_count'] * ZONE_MODEL_SETTINGS['latent_heat_per_occupant_w']                                                                            
    room_moisture_generation_kg_s = latent_room_load_w / 1e3 / 2501.0                                                                                         
    dry_air_mass_flow_kg_s = ZONE_MODEL_SETTINGS['dry_air_mass_flow_kg_s']; coil_bypass_factor = ZONE_MODEL_SETTINGS['coil_bypass_factor']
    room_humidity_ratio = w_from_RH(room_temperature_c, 0.5); humidity_iteration_count = 0
    while True:                                                            
        humidity_iteration_count += 1
        mixed_air_temperature_c, mixed_air_humidity_ratio, mixed_air_enthalpy_kj_kg_dry_air = mix(room_temperature_c, room_humidity_ratio, outdoor_temperature_c, outdoor_humidity_ratio, outdoor_air_fraction)
                                        
        initial_supply_humidity_ratio = mixed_air_humidity_ratio
        supply_temperature_c = room_temperature_c - sensible_room_load_w / 1e3 / (dry_air_mass_flow_kg_s * (1.006 + 1.86 * initial_supply_humidity_ratio))                                     
        apparatus_dew_point_c = (supply_temperature_c - coil_bypass_factor * mixed_air_temperature_c) / (1 - coil_bypass_factor)
        humidity_ratio_at_apparatus_dew_point = wsat(apparatus_dew_point_c)
        coil_is_wet = humidity_ratio_at_apparatus_dew_point < mixed_air_humidity_ratio
        if forced_apparatus_dew_point_c is not None and coil_is_wet:
            apparatus_dew_point_c = forced_apparatus_dew_point_c; humidity_ratio_at_apparatus_dew_point = wsat(apparatus_dew_point_c)
            coil_leaving_air_temperature_c = apparatus_dew_point_c + coil_bypass_factor * (mixed_air_temperature_c - apparatus_dew_point_c)                                                          
        elif coil_is_wet:
            coil_leaving_air_temperature_c = supply_temperature_c
        else:
            coil_leaving_air_temperature_c = supply_temperature_c
        supply_humidity_ratio = humidity_ratio_at_apparatus_dew_point + coil_bypass_factor * (mixed_air_humidity_ratio - humidity_ratio_at_apparatus_dew_point) if coil_is_wet else mixed_air_humidity_ratio
        supply_temperature_c = room_temperature_c - sensible_room_load_w / 1e3 / (dry_air_mass_flow_kg_s * (1.006 + 1.86 * supply_humidity_ratio))
        if coil_is_wet and forced_apparatus_dew_point_c is None:
            apparatus_dew_point_c = (supply_temperature_c - coil_bypass_factor * mixed_air_temperature_c) / (1 - coil_bypass_factor); humidity_ratio_at_apparatus_dew_point = wsat(apparatus_dew_point_c); supply_humidity_ratio = humidity_ratio_at_apparatus_dew_point + coil_bypass_factor * (mixed_air_humidity_ratio - humidity_ratio_at_apparatus_dew_point); coil_leaving_air_temperature_c = supply_temperature_c
        updated_room_humidity_ratio = supply_humidity_ratio + room_moisture_generation_kg_s / dry_air_mass_flow_kg_s
        if abs(updated_room_humidity_ratio - room_humidity_ratio) < humidity_ratio_tolerance or humidity_iteration_count > 500: room_humidity_ratio = updated_room_humidity_ratio; break
        room_humidity_ratio = updated_room_humidity_ratio
    room_relative_humidity = RH_from_w(room_temperature_c, room_humidity_ratio)
    coil_leaving_air_enthalpy_kj_kg_dry_air = h_moist(coil_leaving_air_temperature_c, supply_humidity_ratio)
    condensate_mass_flow_kg_s = dry_air_mass_flow_kg_s * (mixed_air_humidity_ratio - supply_humidity_ratio) if coil_is_wet else 0.0                                                        
    condensate_temperature_c = apparatus_dew_point_c if coil_is_wet else np.nan
    coil_load_kw = dry_air_mass_flow_kg_s * (mixed_air_enthalpy_kj_kg_dry_air - coil_leaving_air_enthalpy_kj_kg_dry_air) - condensate_mass_flow_kg_s * hf_w(condensate_temperature_c if coil_is_wet else 0.0)        
    reheat_power_kw = dry_air_mass_flow_kg_s * (1.006 + 1.86 * supply_humidity_ratio) * (supply_temperature_c - coil_leaving_air_temperature_c) if supply_temperature_c > coil_leaving_air_temperature_c + 1e-9 else 0.0       
    supply_air_enthalpy_kj_kg_dry_air = h_moist(supply_temperature_c, supply_humidity_ratio)
    air_handler_energy_balance_residual_kw = dry_air_mass_flow_kg_s * (mixed_air_enthalpy_kj_kg_dry_air - supply_air_enthalpy_kj_kg_dry_air) - (coil_load_kw + condensate_mass_flow_kg_s * hf_w(condensate_temperature_c if coil_is_wet else 0.0) - reheat_power_kw)                                                                               
    return dict(name=condition_name, To=outdoor_temperature_c, RHo=outdoor_relative_humidity, wo=outdoor_humidity_ratio, Tm=mixed_air_temperature_c, wm=mixed_air_humidity_ratio, hm=mixed_air_enthalpy_kj_kg_dry_air, ADP=apparatus_dew_point_c if coil_is_wet else np.nan,
                wet=coil_is_wet, Ts=supply_temperature_c, Tcoil=coil_leaving_air_temperature_c, ws=supply_humidity_ratio, wr=room_humidity_ratio, RH=room_relative_humidity, viol=room_relative_humidity > ZONE_MODEL_SETTINGS['maximum_room_relative_humidity'] + 1e-9,
                Qs=sensible_room_load_w, Ql=latent_room_load_w, Qcoil=coil_load_kw, ahu_bal=air_handler_energy_balance_residual_kw, Qheat_equiv=max(-coil_load_kw, 0.0), Qcool_demand=max(coil_load_kw, 0.0), mcond=condensate_mass_flow_kg_s * 3600, Qreheat=reheat_power_kw, iters=humidity_iteration_count,
                Qevap=max(coil_load_kw, 0) * ZONE_MODEL_SETTINGS['zone_count'] * (1 + ZONE_MODEL_SETTINGS['distribution_heat_gain_fraction'] + ZONE_MODEL_SETTINGS['pump_heat_gain_fraction']),                               
                cond_m3_day=condensate_mass_flow_kg_s * 3600 * 24 / 1000)

def fix_ceiling(condition_name, outdoor_air_fraction=0.30):
    uncontrolled_zone_result = zone(condition_name, outdoor_air_fraction)
    if not uncontrolled_zone_result['viol']: return uncontrolled_zone_result
    humidity_limit_residual = lambda adp: zone(condition_name, outdoor_air_fraction, forced_apparatus_dew_point_c=adp)['RH'] - ZONE_MODEL_SETTINGS['maximum_room_relative_humidity'] + 1e-4
    apparatus_dew_point_c = brentq(humidity_limit_residual, -1.0, uncontrolled_zone_result['ADP'])
    return zone(condition_name, outdoor_air_fraction, forced_apparatus_dew_point_c=apparatus_dew_point_c)

if __name__ == '__main__':
             
    Pr, Po = psat(22), psat(40)
    wr, wo = w_from_RH(22, .75), w_from_RH(40, .5)
    Tm, wm, hm = mix(22, wr, 40, wo, .3)
    print('Psat22 %.3f Psat40 %.3f | wr %.5f hr %.2f | wo %.5f ho %.2f | Tm %.2f wm %.5f hm %.2f'
          % (Pr, Po, wr, h_moist(22, wr), wo, h_moist(40, wo), Tm, wm, hm))
    print('Twb summer day %.2f monsoon day %.2f' % (Twb(40, wo), Twb(35, w_from_RH(35, .85))))
    for n in WEATHER_CONDITIONS:
        zone_result = zone(n); print('%-14s Tm %.2f ADP %.2f ws %.5f RH %.1f%% Qcoil %.3f kW cond %.2f kg/h Qevap %.0f kW %s'
            % (n, zone_result['Tm'], zone_result['ADP'], zone_result['ws'], zone_result['RH'] * 100, zone_result['Qcoil'], zone_result['mcond'], zone_result['Qevap'], 'VIOL' if zone_result['viol'] else ''))


def zone_float(condition_name, outdoor_air_fraction=0.30):
    outdoor_temperature_c, outdoor_relative_humidity, solar_gain_factor = WEATHER_CONDITIONS[condition_name]; outdoor_humidity_ratio = w_from_RH(outdoor_temperature_c, outdoor_relative_humidity)
    sensible_room_load_w = ZONE_MODEL_SETTINGS['envelope_heat_transfer_w_per_k'] * (outdoor_temperature_c - ZONE_MODEL_SETTINGS['room_temperature_c']) + solar_gain_factor * ZONE_MODEL_SETTINGS['solar_heat_gain_w'] + ZONE_MODEL_SETTINGS['occupant_count'] * ZONE_MODEL_SETTINGS['sensible_heat_per_occupant_w'] + ZONE_MODEL_SETTINGS['lighting_heat_w_per_m2'] * ZONE_MODEL_SETTINGS['floor_area_m2']                                                                      
    room_moisture_generation_kg_s = ZONE_MODEL_SETTINGS['occupant_count'] * ZONE_MODEL_SETTINGS['latent_heat_per_occupant_w'] / 1e3 / 2501.0; dry_air_mass_flow_kg_s = ZONE_MODEL_SETTINGS['dry_air_mass_flow_kg_s']
    def f(trial_room_temperature_c):
        room_humidity_ratio = w_from_RH(trial_room_temperature_c, 0.5)
        for _ in range(200):
            mixed_air_temperature_c, mixed_air_humidity_ratio, mixed_air_enthalpy_kj_kg_dry_air = mix(trial_room_temperature_c, room_humidity_ratio, outdoor_temperature_c, outdoor_humidity_ratio, outdoor_air_fraction); room_humidity_ratio = mixed_air_humidity_ratio + room_moisture_generation_kg_s / dry_air_mass_flow_kg_s
        sensible_room_load_at_trial_w = ZONE_MODEL_SETTINGS['envelope_heat_transfer_w_per_k'] * (outdoor_temperature_c - trial_room_temperature_c) + solar_gain_factor * ZONE_MODEL_SETTINGS['solar_heat_gain_w'] + ZONE_MODEL_SETTINGS['occupant_count'] * ZONE_MODEL_SETTINGS['sensible_heat_per_occupant_w'] + ZONE_MODEL_SETTINGS['lighting_heat_w_per_m2'] * ZONE_MODEL_SETTINGS['floor_area_m2']
        return trial_room_temperature_c - (mixed_air_temperature_c + sensible_room_load_at_trial_w / 1e3 / (dry_air_mass_flow_kg_s * (1.006 + 1.86 * mixed_air_humidity_ratio))), room_humidity_ratio
    room_temperature_c = brentq(lambda dry_bulb_temperature_c: f(dry_bulb_temperature_c)[0], 10.0, 40.0); return dict(T_room_float=room_temperature_c, RH_float=RH_from_w(room_temperature_c, f(room_temperature_c)[1]))
