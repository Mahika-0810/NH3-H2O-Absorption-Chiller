"""condenser.py - Task 2: air-cooled finned-tube vs evaporative condenser, design point, off-design (UA fixed) and feasibility checks."""
import numpy as np
from scipy.optimize import brentq
import air_handler as air_handler_model
import cycle as absorption_cycle_model

OVERALL_HEAT_TRANSFER_COEFFICIENT_W_M2_K = {'air': 30.0, 'evap': 600.0}                      
DESIGN_TEMPERATURE_APPROACH_K = {'air': 12.0, 'evap': 8.0}                                                    
DESIGN_SINK_TEMPERATURE_RISE_K = {'air': 8.0, 'evap': 4.0}                                                                                             
EVAPORATIVE_CONCENTRATION_CYCLES = 4.0
SINK_SIDE_UA_FLOW_EXPONENT = 0.8                                                                                                                         

def sink_T(condition_name, condenser_type):
    outdoor_dry_bulb_temperature_c, outdoor_relative_humidity, _ = air_handler_model.WEATHER_CONDITIONS[condition_name]
    return outdoor_dry_bulb_temperature_c if condenser_type == 'air' else air_handler_model.Twb(outdoor_dry_bulb_temperature_c, air_handler_model.w_from_RH(outdoor_dry_bulb_temperature_c, outdoor_relative_humidity))

def governing(condenser_type):
    return max(air_handler_model.WEATHER_CONDITIONS, key=lambda n: sink_T(n, condenser_type))

def design(condenser_type, cooling_load_by_condition):
    governing_condition = governing(condenser_type); heat_sink_temperature_c = sink_T(governing_condition, condenser_type); design_condenser_temperature_k = heat_sink_temperature_c + DESIGN_TEMPERATURE_APPROACH_K[condenser_type] + 273.15
    cycle_result = absorption_cycle_model.solve_cycle(cooling_load_by_condition[governing_condition], design_condenser_temperature_k); condenser_heat_rejection_kw = cycle_result['Qcond']                      
    hot_end_temperature_difference_k, cold_end_temperature_difference_k = DESIGN_TEMPERATURE_APPROACH_K[condenser_type], DESIGN_TEMPERATURE_APPROACH_K[condenser_type] - DESIGN_SINK_TEMPERATURE_RISE_K[condenser_type]
    log_mean_temperature_difference_k = (hot_end_temperature_difference_k - cold_end_temperature_difference_k) / np.log(hot_end_temperature_difference_k / cold_end_temperature_difference_k); overall_conductance_w_k = condenser_heat_rejection_kw * 1e3 / log_mean_temperature_difference_k                            
    sink_heat_capacity_rate_kw_k = condenser_heat_rejection_kw / DESIGN_SINK_TEMPERATURE_RISE_K[condenser_type]                                                                                
    design_result = dict(type=condenser_type, governing=governing_condition, T_sink=heat_sink_temperature_c, Tcond=design_condenser_temperature_k - 273.15, Phigh=cycle_result['Phigh'], Qcond=condenser_heat_rejection_kw, LMTD=log_mean_temperature_difference_k, UA=overall_conductance_w_k, area=overall_conductance_w_k / OVERALL_HEAT_TRANSFER_COEFFICIENT_W_M2_K[condenser_type], C=sink_heat_capacity_rate_kw_k, cycle=cycle_result)
    if condenser_type == 'air':
        design_result['air_kg_s'] = sink_heat_capacity_rate_kw_k / 1.006
    else:
        outdoor_dry_bulb_temperature_c, outdoor_relative_humidity, _ = air_handler_model.WEATHER_CONDITIONS[governing_condition]; inlet_air_humidity_ratio = air_handler_model.w_from_RH(outdoor_dry_bulb_temperature_c, outdoor_relative_humidity); inlet_air_enthalpy_kj_kg_dry_air = air_handler_model.h_moist(outdoor_dry_bulb_temperature_c, inlet_air_humidity_ratio)
        outlet_air_temperature_c = heat_sink_temperature_c + DESIGN_SINK_TEMPERATURE_RISE_K['evap']; outlet_air_enthalpy_kj_kg_dry_air = air_handler_model.h_moist(outlet_air_temperature_c, air_handler_model.wsat(outlet_air_temperature_c))                                      
        design_result['air_kg_s'] = condenser_heat_rejection_kw / (outlet_air_enthalpy_kj_kg_dry_air - inlet_air_enthalpy_kj_kg_dry_air)
        evaporation_rate_kg_s = condenser_heat_rejection_kw / (2501.0 - 2.36 * heat_sink_temperature_c)                                                                                   
        design_result['evap_kg_s'] = evaporation_rate_kg_s; design_result['makeup_kg_s'] = evaporation_rate_kg_s * EVAPORATIVE_CONCENTRATION_CYCLES / (EVAPORATIVE_CONCENTRATION_CYCLES - 1); design_result['blowdown_kg_s'] = evaporation_rate_kg_s / (EVAPORATIVE_CONCENTRATION_CYCLES - 1)
        design_result['makeup_m3_day'] = design_result['makeup_kg_s'] * 86.4
    design_result['feasible_Pmax'] = cycle_result['chk']['Phigh_ok']
    return design_result

def _eps_C(condenser_design, size_factor, sink_flow_fraction):
    sink_heat_capacity_rate_kw_k = condenser_design['C'] * size_factor * sink_flow_fraction; overall_conductance_w_k = condenser_design['UA'] * size_factor * sink_flow_fraction ** SINK_SIDE_UA_FLOW_EXPONENT
    return 1 - np.exp(-overall_conductance_w_k / 1e3 / sink_heat_capacity_rate_kw_k), sink_heat_capacity_rate_kw_k

def _pack(condenser_design, condition_name, evaporator_cooling_kw, condenser_temperature_k, heat_sink_temperature_k, size_factor, sink_flow_fraction):
    cycle_result = absorption_cycle_model.solve_cycle(evaporator_cooling_kw, condenser_temperature_k)
    return dict(cond=condition_name, status='ok', Qevap=evaporator_cooling_kw, Tsink=heat_sink_temperature_k - 273.15, Tcond=condenser_temperature_k - 273.15, Phigh_bar=cycle_result['Phigh'] / 100, Tgen=cycle_result['Tgen'] - 273.15,
                COP=cycle_result['COP'], Qgen=cycle_result['Qgen'], Qcond=cycle_result['Qcond'], mfg=cycle_result['mfg'], Tfg_out=cycle_result['Tfg_out'] - 273.15, sink_flow_fraction=sink_flow_fraction, size_factor=size_factor,
                UA_kW_K=condenser_design['UA'] / 1e3 * size_factor * sink_flow_fraction ** SINK_SIDE_UA_FLOW_EXPONENT, C_kW_K=condenser_design['C'] * size_factor * sink_flow_fraction,
                a_Pmax=cycle_result['chk']['Phigh_ok'], b_gen_ok=cycle_result['chk']['Tgen_deliverable'], flue_limit_ok=cycle_result['chk']['mfg_ok'], overall_thermal_feasible=(cycle_result['chk']['Phigh_ok'] and cycle_result['chk']['Tgen_deliverable'] and cycle_result['chk']['mfg_ok'] and cycle_result['chk']['Tcond_ok']), fg_floor_active=cycle_result['fg_floor_active'], eps_gen_eff=cycle_result['eps_gen_eff'], c_Tcond=cycle_result['chk']['Tcond_ok'], cycle=cycle_result)

def offdesign(condenser_design, condition_name, evaporator_cooling_kw, size_factor=1.0, sink_flow_fraction=1.0):
    condenser_type = condenser_design['type']; heat_sink_temperature_k = sink_T(condition_name, condenser_type) + 273.15; heat_exchanger_effectiveness, sink_heat_capacity_rate_kw_k = _eps_C(condenser_design, size_factor, sink_flow_fraction); evaporator_cooling_kw = evaporator_cooling_kw
    if evaporator_cooling_kw <= 1e-6: return dict(cond=condition_name, status='no cooling load')
    heat_balance_residual = lambda trial_condenser_temperature_k: absorption_cycle_model.solve_cycle(evaporator_cooling_kw, trial_condenser_temperature_k)['Qcond'] - heat_exchanger_effectiveness * sink_heat_capacity_rate_kw_k * (trial_condenser_temperature_k - heat_sink_temperature_k)
    condenser_temperature_k = brentq(heat_balance_residual, heat_sink_temperature_k + 0.3, heat_sink_temperature_k + 45.0, xtol=1e-9)
    return _pack(condenser_design, condition_name, evaporator_cooling_kw, condenser_temperature_k, heat_sink_temperature_k, size_factor, sink_flow_fraction)

def offdesign_at_Tcond(condenser_design, condition_name, evaporator_cooling_kw, target_condenser_temperature_c, size_factor=1.0, minimum_sink_flow_fraction=0.10):
    heat_sink_temperature_k = sink_T(condition_name, condenser_design['type']) + 273.15; target_condenser_temperature_k = target_condenser_temperature_c + 273.15; required_condenser_heat_rejection_kw = absorption_cycle_model.solve_cycle(evaporator_cooling_kw, target_condenser_temperature_k)['Qcond']
    heat_rejection_capacity_at_flow = lambda trial_sink_flow_fraction: (lambda exchanger_effectiveness, sink_capacity_rate_kw_k: exchanger_effectiveness * sink_capacity_rate_kw_k * (target_condenser_temperature_k - heat_sink_temperature_k))(*_eps_C(condenser_design, size_factor, trial_sink_flow_fraction))
    if heat_rejection_capacity_at_flow(minimum_sink_flow_fraction) > required_condenser_heat_rejection_kw: return None
    sink_flow_fraction = brentq(lambda trial_sink_flow_fraction: heat_rejection_capacity_at_flow(trial_sink_flow_fraction) - required_condenser_heat_rejection_kw, minimum_sink_flow_fraction, 1.0, xtol=1e-12)
    return _pack(condenser_design, condition_name, evaporator_cooling_kw, target_condenser_temperature_k, heat_sink_temperature_k, size_factor, sink_flow_fraction)

def controlled(condenser_design, condition_name, evaporator_cooling_kw, minimum_condenser_temperature_c=25.0, size_factor=1.0):
    uncontrolled_result = offdesign(condenser_design, condition_name, evaporator_cooling_kw, size_factor=size_factor)
    if uncontrolled_result['status'] != 'ok' or uncontrolled_result['Tcond'] >= minimum_condenser_temperature_c - 1e-6:
        uncontrolled_result['action'] = 'none (floating T_cond >= %.0f C)' % minimum_condenser_temperature_c if uncontrolled_result['status'] == 'ok' else 'chiller off'; return uncontrolled_result
    controlled_result = offdesign_at_Tcond(condenser_design, condition_name, evaporator_cooling_kw, minimum_condenser_temperature_c, size_factor=size_factor)
    if controlled_result is None:
        controlled_result = offdesign(condenser_design, condition_name, evaporator_cooling_kw, size_factor=size_factor, sink_flow_fraction=0.10); controlled_result['action'] = 'fans at minimum, T_cond still < %.0f C (needs water-side bypass)' % minimum_condenser_temperature_c; return controlled_result
    controlled_result['action'] = 'sink flow throttled to %.0f %% (fan speed / spray-pump staging) to hold T_cond = %.2f C' % (100 * controlled_result['sink_flow_fraction'], minimum_condenser_temperature_c); return controlled_result

def qevap_at_fg_limit(condenser_design, condition_name, maximum_flue_gas_mass_flow_kg_s=80.0, size_factor=1.0, maximum_trial_cooling_kw=20000.0):
    residual = lambda trial_cooling_duty_kw: offdesign(condenser_design, condition_name, trial_cooling_duty_kw, size_factor=size_factor)['mfg'] - maximum_flue_gas_mass_flow_kg_s
    trial_grid = np.linspace(200.0, maximum_trial_cooling_kw, 30)
    residual_values = [residual(float(trial_duty)) for trial_duty in trial_grid]
    if residual_values[0] > 0:
        return 0.0
    for grid_index in range(1, len(trial_grid)):
        if residual_values[grid_index] >= 0:
            return brentq(residual, float(trial_grid[grid_index - 1]), float(trial_grid[grid_index]), xtol=1e-3)
    return float(maximum_trial_cooling_kw)

def size_for_fg_limit(condenser_design, condition_name, evaporator_cooling_kw, maximum_flue_gas_mass_flow_kg_s=80.0):
    residual = lambda trial_size_factor: offdesign(condenser_design, condition_name, evaporator_cooling_kw, size_factor=trial_size_factor)['mfg'] - maximum_flue_gas_mass_flow_kg_s
    if residual(1.0) <= 0:
        return 1.0
    trial_grid = np.linspace(1.0, 6.0, 26)
    previous_factor = float(trial_grid[0]); previous_residual = residual(previous_factor)
    for trial_factor in trial_grid[1:]:
        trial_factor = float(trial_factor); trial_residual = residual(trial_factor)
        if trial_residual <= 0:
            return brentq(residual, previous_factor, trial_factor, xtol=1e-4)
        previous_factor, previous_residual = trial_factor, trial_residual
    return 6.0
