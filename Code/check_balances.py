"""check_engineering_balances.py - automated validation (pytest) + Validation Table.   Run:  pytest -q check_engineering_balances.py     or     python check_engineering_balances.py
Hand values are from the project PDF (Section 2.4). Tolerances: moist air 0.5 % (w,h) / 0.1 K (T); ammonia bands from the PDF (Psat +-5 %, hfg +-8 %).
Check-3 balances and COP bands are tested on the cycle built on the own Patek-Klomfar module (solution_props); any band miss is reported, not hidden."""
import numpy as np, pandas as pd, os
try:
    import pytest
except ImportError:                                                                                   
    pytest = None
import air_handler as air_handler_model, nh3_props as ammonia_properties, cycle as absorption_cycle, solution_props as solution_properties, equilibrium as equilibrium_model

def rel(calculated_value, reference_value): return 100 * (calculated_value - reference_value) / reference_value
def _check1():
    room_humidity_ratio, outdoor_humidity_ratio = air_handler_model.w_from_RH(22, .75), air_handler_model.w_from_RH(40, .5); mixed_air_temperature_c, mixed_air_humidity_ratio, mixed_air_enthalpy_kj_kg_dry_air = air_handler_model.mix(22, room_humidity_ratio, 40, outdoor_humidity_ratio, .3)
    monsoon_outdoor_humidity_ratio = air_handler_model.w_from_RH(35, .85)
    return [('Pv room 22C/75% [kPa]', 1.984, .75 * air_handler_model.psat(22), .5), ('w room', .01242, room_humidity_ratio, .5), ('h room [kJ/kg]', 53.71, air_handler_model.h_moist(22, room_humidity_ratio), .5),
            ('Pv out 40C/50% [kPa]', 3.692, .5 * air_handler_model.psat(40), .5), ('w outdoor', .02352, outdoor_humidity_ratio, .5), ('h outdoor', 100.82, air_handler_model.h_moist(40, outdoor_humidity_ratio), .5),
            ('w mixed', .01575, mixed_air_humidity_ratio, .5), ('h mixed', 67.84, mixed_air_enthalpy_kj_kg_dry_air, .5), ('T mixed [C] (abs 0.1 K)', 27.48, mixed_air_temperature_c, 'K0.1'),
            ('Twb summer day [C] (abs 0.15 K)', 30.3, air_handler_model.Twb(40, outdoor_humidity_ratio), 'K0.15'), ('Twb monsoon day [C] (abs 0.15 K)', 32.6, air_handler_model.Twb(35, monsoon_outdoor_humidity_ratio), 'K0.15')]
def _check2():
    ammonia_saturation_at_0c, ammonia_saturation_at_40c = ammonia_properties.sat(273.15), ammonia_properties.sat(313.15); evaporator_specific_cooling_kj_kg = ammonia_saturation_at_0c['hg'] - ammonia_saturation_at_40c['hf']
    return [('NH3 Psat 0C [kPa] (band +-5%)', 429.6, ammonia_saturation_at_0c['P'], 5), ('NH3 Psat 40C [kPa] (band +-5%)', 1554.9, ammonia_saturation_at_40c['P'], 5),
            ('NH3 hfg 0C [kJ/kg] (band +-8%)', 1261.8, ammonia_saturation_at_0c['hfg'], 8), ('NH3 hfg 40C [kJ/kg] (band +-8%)', 1098.8, ammonia_saturation_at_40c['hfg'], 8),
            ('q_e = hg(0C)-hf(40C) in 1050-1100', 1075.0, evaporator_specific_cooling_kj_kg, 'rng:1050:1100'), ('m_r @5000 kW in 4.5-4.8 kg/s', 4.65, 5000 / evaporator_specific_cooling_kj_kg, 'rng:4.5:4.8')]
def _check3():
    cycle_result = absorption_cycle.solve_cycle(5000., 313.15); checks = []
    checks += [('Energy bal: (Qgen+Qevap+W)-(Qcond+Qabs+Qrect) /Qgen', 0.0, cycle_result['bal_E_rel'], 'abs1e-6'), ('Mass bal: ms-(mr+mw) [kg/s]', 0.0, cycle_result['bal_m'], 'abs1e-9'),
            ('NH3 bal: ms*xs-(mw*xw+mr*0.995)', 0.0, cycle_result['bal_NH3'], 'abs1e-9'), ('SHX energy bal [kW]', 0.0, cycle_result['shx_E'], 'abs1e-6')]
    generator_temperature_k, condenser_temperature_k, evaporator_temperature_k = 403.15, 318.15, 275.15; checks.append(('COP_rev(Te=2C,Tc=45C,Tg=130C) ~1.35', 1.35, (evaporator_temperature_k / (condenser_temperature_k - evaporator_temperature_k)) * ((generator_temperature_k - condenser_temperature_k) / generator_temperature_k), 1.0))
    cycle_at_45c = absorption_cycle.solve_cycle(5000., condenser_temperature_k); checks.append(('COP < COP_rev (Tc=45C)', 1.0, float(cycle_at_45c['COP'] < cycle_at_45c['COP_rev']), 'rng:0.5:1.5'))
    cop_values = [absorption_cycle.solve_cycle(5000., T + 273.15)['COP'] for T in (30, 35, 40, 45)]; checks.append(('COP falls as Tcond rises (30..45C)', 1.0, float(all(np.diff(cop_values) < 0)), 'rng:0.5:1.5'))
    checks.append(('COP screening band 0.45-0.65 at Tcond=40C (review only; model-dependent)', .55, absorption_cycle.solve_cycle(5000., 313.15)['COP'], 'rng:0.45:0.65'))
    checks.append(('COP screening band 0.45-0.65 at Tcond=45C (review only; model-dependent)', .55, cycle_at_45c['COP'], 'rng:0.45:0.65'))
    return checks
def _props():
    ammonia_mole_fraction_values = np.linspace(0.02, 0.99, 25); liquid_enthalpy_function = lambda T, x: solution_properties.h_liquid(0, x, T)
    return [('P&K water boils at 101.325 kPa [K] (abs 2.5 K, paper band)', 373.12, solution_properties.T_bub(101.325, 0.0), 'K2.5'),
            ('Cross-model: P&K h_fg NH3 0C vs EOS (band +-5 %)', ammonia_properties.sat(273.15)['hfg'], solution_properties.h_NH3_vap(273.15) - solution_properties.h_NH3_liq(273.15), 5),
            ('Cross-model: P&K h_fg NH3 40C vs EOS (band +-5 %)', ammonia_properties.sat(313.15)['hfg'], solution_properties.h_NH3_vap(313.15) - solution_properties.h_NH3_liq(313.15), 5),
            ('Cycle q_e (P&K scale, 0C vap - 40C liq) in 1050-1100', 1075.0, solution_properties.h_NH3_vap(273.15) - solution_properties.h_NH3_liq(313.15), 'rng:1050:1100'),
            ('Cycle states 2-4 share the P&K pure-ammonia enthalpy basis (abs 1e-9)', absorption_cycle.solve_cycle(5000., 313.15)['h']['h4'] - absorption_cycle.solve_cycle(5000., 313.15)['h']['h2'], solution_properties.h_NH3_vap(275.15) - solution_properties.h_NH3_liq(313.15), 'abs1e-9'),
            ('P&K pure NH3 end = EOS Tsat @1 MPa [K] (abs 2 K)', ammonia_properties.T_sat(1000.), solution_properties.T_bub(1000., 1.0), 'K2'),
            ('P&K T_bub monotonic in x (1.7 MPa)', 1.0, float(all(np.diff([solution_properties.T_bub(1700., x) for x in ammonia_mole_fraction_values]) < 0)), 'rng:0.5:1.5'),
            ('P&K bubble/dew round trip dT [K]', 0.0, solution_properties.T_dew(1500., solution_properties.y_eq(1500., .3)) - solution_properties.T_bub(1500., .3), 'abs1e-6'),
            ('P&K strong x at 40C/Plow(2C) in 0.35-0.50 mass', .42, solution_properties.mole_fraction_to_mass_fraction(solution_properties.x_bub(ammonia_properties.psat(275.15), 313.15)), 'rng:0.35:0.50'),
            ('P&K h_liq water 100C vs h_f 419.1 [kJ/kg] (abs 6, PK scale)', 419.1, equilibrium_model.h_liq(373.15, 0.0), 'K6'),
            ('P&K h_liq NH3 at 273.16 K = 0 (abs 1, PK scale)', 0.0, equilibrium_model.h_liq(273.16, 1.0), 'K1'),
            ('P&K h_fg NH3 0C vs 1261.8 [kJ/kg] (abs 25)', 1261.8, equilibrium_model.h_gas(273.15, 1.0) - equilibrium_model.h_liq(273.15, 1.0), 'K25'),
            ('Ref state: P&K sat. vapour NH3 at 0C = EOS on project reference (abs 1)', ammonia_properties.sat(273.15)['hg'], solution_properties.h_NH3_vap(273.15), 'K1')]
def _ops():
    import condenser as cd
    zones_30_percent_fresh_air = {condition_name: air_handler_model.fix_ceiling(condition_name, 0.30) for condition_name in air_handler_model.WEATHER_CONDITIONS}; zones_40_percent_fresh_air = {condition_name: air_handler_model.fix_ceiling(condition_name, 0.40) for condition_name in air_handler_model.WEATHER_CONDITIONS}; cooling_loads_30_percent = {condition_name: zones_30_percent_fresh_air[condition_name]['Qevap'] for condition_name in zones_30_percent_fresh_air}
    condenser_design = cd.design('evap', cooling_loads_30_percent)
                                                                              
                                                                          
    condenser_size_factor = cd.size_for_fg_limit(
        condenser_design, 'Monsoon day',
        zones_40_percent_fresh_air['Monsoon day']['Qevap'],
        maximum_flue_gas_mass_flow_kg_s=76.0,
    )
    conditions_with_cooling_load = [condition_name for condition_name in zones_30_percent_fresh_air if zones_30_percent_fresh_air[condition_name]['Qevap'] > 1e-6]
    controlled_results_30_percent = {
        condition_name: cd.controlled(
            condenser_design, condition_name,
            zones_30_percent_fresh_air[condition_name]['Qevap'],
            size_factor=condenser_size_factor,
        )
        for condition_name in conditions_with_cooling_load
    }
    conditions_with_cooling_load_40_percent = [condition_name for condition_name in zones_40_percent_fresh_air if zones_40_percent_fresh_air[condition_name]['Qevap'] > 1e-6]
    controlled_results_40_percent = {condition_name: cd.controlled(condenser_design, condition_name, zones_40_percent_fresh_air[condition_name]['Qevap'], size_factor=condenser_size_factor) for condition_name in conditions_with_cooling_load_40_percent}
    monsoon_day_result = controlled_results_40_percent['Monsoon day']
    throttled_condenser_temperatures = [cycle_result['Tcond'] for cycle_result in controlled_results_30_percent.values() if 'throttled' in cycle_result['action']]
    all_controlled_zone_results = [air_handler_model.fix_ceiling(condition_name, fresh_air_fraction) for condition_name in air_handler_model.WEATHER_CONDITIONS for fresh_air_fraction in (0.30, 0.40)]
    checks = [('Final resized design, 30 % fresh: max flue gas <= 80 kg/s', 70.0, max(cycle_result['mfg'] for cycle_result in controlled_results_30_percent.values()), 'rng:0:80'),
            ('Final design: max P_high <= 20 bar (controlled)', 15.0, max(cycle_result['Phigh_bar'] for cycle_result in controlled_results_30_percent.values()), 'rng:0:20'),
            ('Final design: min controlled T_cond >= 25 C (tol 1e-6 K)', 25.0, min(cycle_result['Tcond'] for cycle_result in controlled_results_30_percent.values()), 'rng:24.999999:60'),
            ('Throttled states: controlled T_cond = 25.00 C setpoint (abs 1e-6)', 25.0, float(np.mean(throttled_condenser_temperatures)), 'abs1e-6'),
            ('Winter night is a zero-load (chiller OFF) state', 0.0, zones_30_percent_fresh_air['Winter night']['Qevap'], 'abs1e-6'),
            ('Winter night: raw coil balance < 0 is NOT booked as cooling demand (Qcool = 0)', 0.0, zones_30_percent_fresh_air['Winter night']['Qcool_demand'], 'abs1e-9'),
            ('Final design: check (b) T_gen deliverable (flue gas leaves >= 140 C and >= 5 K above T_gen), all 30 % states', 1.0, float(all(cycle_result['b_gen_ok'] for cycle_result in controlled_results_30_percent.values())), 'rng:0.5:1.5'),
            ('Monsoon day 40 % fresh, resized condenser: flue gas <= 80 kg/s', 76.0, monsoon_day_result['mfg'], 'rng:0:80'),
            ('Monsoon day 40 % fresh, resized condenser: P_high <= 20 bar', 16.0, monsoon_day_result['Phigh_bar'], 'rng:0:20'),
            ('40 % fresh, ONE resized condenser, all conditions: max flue gas <= 80 kg/s', 76.0, max(cycle_result['mfg'] for cycle_result in controlled_results_40_percent.values()), 'rng:0:80'),
            ('40 % fresh, ONE resized condenser, all conditions: check (b) T_gen deliverable', 1.0, float(all(cycle_result['b_gen_ok'] for cycle_result in controlled_results_40_percent.values())), 'rng:0.5:1.5'),
            ('40 % fresh, ONE resized condenser, all conditions: max P_high <= 20 bar', 16.0, max(cycle_result['Phigh_bar'] for cycle_result in controlled_results_40_percent.values()), 'rng:0:20'),
            ('Room RH <= 75 % in all 12 zone cases (30 / 40 % fresh, controlled)', 75.0, 100 * max(z['RH'] for z in all_controlled_zone_results), 'rng:0:75'),
            ('AHU enthalpy balance incl. reheat + condensate, max |residual| [kW]', 0.0, max(abs(z['ahu_bal']) for z in all_controlled_zone_results), 'abs1e-9')]
    import os, pandas as pd
    if os.path.exists('results/exhaust_gas_limit_solutions.csv') and os.path.exists('results/resized_condenser_feasibility.csv'):                                                      
        remedy_csv_row = pd.read_csv('results/exhaust_gas_limit_solutions.csv').set_index('Condition').loc['Monsoon day']; resized_condenser_csv_row = pd.read_csv('results/resized_condenser_feasibility.csv').set_index('Condition').loc['Monsoon day']
        checks += [('Mitigation CSVs agree (same rerun): flue gas, remedy vs design table [kg/s]', remedy_csv_row.optA_mfg_kg_s, resized_condenser_csv_row.MIT_flue_kg_s, 'abs1e-6'),
                ('Mitigation CSVs agree: T_cond [C]', remedy_csv_row.optA_Tcond_C, resized_condenser_csv_row.MIT_Tcond_C, 'abs1e-6'), ('Mitigation CSVs agree: P_high [bar]', remedy_csv_row.optA_Phigh_bar, resized_condenser_csv_row.MIT_Phigh_bar, 'abs1e-6'),
                ('Mitigation CSVs agree: UA [kW/K]', remedy_csv_row.optA_UA_kW_K, resized_condenser_csv_row.MIT_UA_kW_K, 'abs1e-6'), ('Mitigation CSVs agree: area [m2]', remedy_csv_row.optA_area_m2, resized_condenser_csv_row.condenser_area_m2, 'abs1e-6')]
    return checks
ENGINEERING_CHECKS = _check1() + _check2() + _check3() + _props() + _ops()

def _ok(reference_value, calculated_value, tolerance):
    if isinstance(tolerance, (int, float)): return abs(rel(calculated_value, reference_value)) <= tolerance, rel(calculated_value, reference_value)
    if tolerance.startswith('K'): absolute_tolerance = float(tolerance[1:]); return abs(calculated_value - reference_value) <= absolute_tolerance, calculated_value - reference_value
    if tolerance.startswith('abs'): return abs(calculated_value - reference_value) <= float(tolerance[3:]), calculated_value - reference_value
    lower_limit, upper_limit = map(float, tolerance.split(':')[1:]); return lower_limit <= calculated_value <= upper_limit, rel(calculated_value, reference_value)

                                                                                                                   
                                                                                                                           
                                                                                                                       
CHECK_CATEGORIES = {'A': 'A Independent reference', 'B': 'B Plausibility band', 'C': 'C Cross-model', 'D': 'D Internal consistency', 'E': 'E Design-limit / assumption-dependent'}
INTERNAL_CONSISTENCY_LABELS = ('COP_rev', 'AHU enthalpy balance', 'Mitigation CSVs agree', 'Energy bal', 'Mass bal', 'NH3 bal', 'SHX energy bal', 'COP < COP_rev', 'COP falls', 'Cycle states 2-4', 'P&K T_bub monotonic', 'P&K bubble/dew round trip', 'P&K h_liq NH3 at 273.16')
CROSS_MODEL_LABELS = ('Cross-model', 'Cycle q_e (P&K', 'P&K pure NH3 end', 'Ref state:')
PLAUSIBILITY_LABELS = ('COP band', 'COP screening band', 'P&K strong x')
DESIGN_LIMIT_LABELS = ('Final design', 'Throttled states', 'Winter night', 'Monsoon day 40 %', '40 % fresh, ONE', 'Room RH')
def category(check_name):
    for category_code, category_prefixes in (('D', INTERNAL_CONSISTENCY_LABELS), ('C', CROSS_MODEL_LABELS), ('B', PLAUSIBILITY_LABELS), ('E', DESIGN_LIMIT_LABELS)):
        if check_name.startswith(category_prefixes): return CHECK_CATEGORIES[category_code]
    return CHECK_CATEGORIES['A']

def table():
    table_rows = []
    for check_name, reference_value, calculated_value, tolerance in ENGINEERING_CHECKS:
        passed, difference = _ok(reference_value, calculated_value, tolerance); display_reference, display_calculated, display_difference = reference_value, calculated_value, difference
        if isinstance(tolerance, str) and tolerance.startswith('rng'):                                                                       
            lower_limit, upper_limit = map(float, tolerance.split(':')[1:])
            if (lower_limit, upper_limit) == (0.5, 1.5): display_reference, display_calculated, display_difference = 'True', ('True' if passed else 'False'), 'n/a (logical check)'
            else:
                display_reference = ('<= %g' % upper_limit) if lower_limit == 0 else ('>= 25' if (upper_limit == 60 and abs(lower_limit - 25) < 0.01) else '%g - %g' % (lower_limit, upper_limit)); display_difference = 'inside band' if passed else 'OUTSIDE band'
        format_value = lambda v: float('%.6g' % v) if isinstance(v, (float, np.floating)) else v
        result_status = 'Pass' if passed else ('Review' if category(check_name).startswith(('B ', 'E ')) else 'FAIL')
        table_rows.append(dict(Category=category(check_name), Quantity=check_name, Hand_value=format_value(display_reference), Code_output=format_value(display_calculated), pct_or_abs_diff=format_value(display_difference), Result=result_status))
    results_table = pd.DataFrame(table_rows); results_table['_o'] = results_table.Category.str[0]
    return results_table.sort_values('_o', kind='stable').drop(columns='_o').reset_index(drop=True)

if pytest is not None:
    @pytest.mark.parametrize('name,hand,code,tol', ENGINEERING_CHECKS, ids=[c[0] for c in ENGINEERING_CHECKS])
    def test_check(check_name, reference_value, calculated_value, tolerance):
        assert _ok(reference_value, calculated_value, tolerance)[0], f'{check_name}: hand {reference_value} code {calculated_value}'


if __name__ == '__main__':
    os.makedirs('results', exist_ok=True); t = table(); t.to_csv('results/model_check_summary.csv', index=False); print(t.to_string())
