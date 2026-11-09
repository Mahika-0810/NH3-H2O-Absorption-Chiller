"""run_cycle.py - single reproducible driver for Project 2: python run_cycle.py   (writes ./results/*.csv, *.png, results.json)
PROPERTY MODEL = own Patek-Klomfar (1995) Eqs. 6,7,9,10 (corrected coefficients, see equilibrium.py)."""
import os, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

# Consistent, readable plot formatting across the project. These settings change
# presentation only; they do not change any calculated values.
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.labelsize": 9,
    "axes.edgecolor": "#334155",
    "axes.linewidth": 0.8,
    "axes.facecolor": "#f7f9fc",
    "figure.facecolor": "white",
    "grid.color": "#94a3b8",
    "grid.linestyle": "--",
    "grid.alpha": 0.28,
    "legend.frameon": False,
    "lines.linewidth": 2.0,
    "savefig.facecolor": "white",
})
plt.rcParams["axes.prop_cycle"] = matplotlib.cycler(
    color=["#2457A7", "#0F8B8D", "#D97706", "#C2414B", "#7656A8", "#475569"]
)
import air_handler as air_handler_model, cycle as absorption_cycle_model, condenser as condenser_model, nh3_props as ammonia_properties, solution_props as solution_properties
RESULTS_FOLDER = 'results'; os.makedirs(RESULTS_FOLDER, exist_ok=True)
WEATHER_CONDITIONS = list(air_handler_model.WEATHER_CONDITIONS); ZONE_COUNT = air_handler_model.ZONE_MODEL_SETTINGS['zone_count']; HOURS_PER_SEASON = 12 * 120                                                       
REPORT_DATA = {}                                                                                                                    

def build_zone_results(fresh): return {condition_name: air_handler_model.fix_ceiling(condition_name, fresh) for condition_name in WEATHER_CONDITIONS}                                                              
ZONES_30_PERCENT_FRESH_AIR, ZONES_40_PERCENT_FRESH_AIR = build_zone_results(0.30), build_zone_results(0.40)
RAW_ZONES_30_PERCENT_FRESH_AIR = {condition_name: air_handler_model.zone(condition_name, 0.30) for condition_name in WEATHER_CONDITIONS}
COOLING_LOADS_30_PERCENT_FRESH_AIR = {condition_name: ZONES_30_PERCENT_FRESH_AIR[condition_name]['Qevap'] for condition_name in WEATHER_CONDITIONS}; DESIGN_COOLING_LOAD_KW = max(COOLING_LOADS_30_PERCENT_FRESH_AIR.values()); GOVERNING_CONDITION = max(COOLING_LOADS_30_PERCENT_FRESH_AIR, key=COOLING_LOADS_30_PERCENT_FRESH_AIR.get)

                                                                 
air_handler_rows = []
for condition_name in WEATHER_CONDITIONS:
    raw_air_handler_result, controlled_cycle_result = RAW_ZONES_30_PERCENT_FRESH_AIR[condition_name], ZONES_30_PERCENT_FRESH_AIR[condition_name]
    air_handler_rows.append(dict(Condition=condition_name, Tdb=raw_air_handler_result['To'], RH_out=raw_air_handler_result['RHo'], T_mix=raw_air_handler_result['Tm'], ADP_uncontrolled=raw_air_handler_result['ADP'], w_supply=raw_air_handler_result['ws'], T_supply=raw_air_handler_result['Ts'], w_room=raw_air_handler_result['wr'],
                     RH_room_uncontrolled=100 * raw_air_handler_result['RH'], ceiling_violated=bool(raw_air_handler_result['viol']), ADP_controlled=controlled_cycle_result['ADP'], RH_room_controlled=100 * controlled_cycle_result['RH'],
                     reheat_kW_zone=controlled_cycle_result['Qreheat'], reheat_kW_district=controlled_cycle_result['Qreheat'] * ZONE_COUNT, T_supply_after_reheat=controlled_cycle_result['Ts'], T_coil_leaving=controlled_cycle_result['Tcoil'], w_supply_controlled=controlled_cycle_result['ws'], w_room_controlled=controlled_cycle_result['wr'],
                     Qcoil_kW_zone=controlled_cycle_result['Qcoil'], Qcoil_raw_note=('RAW coil balance m*(h_mix-h_supply) < 0: mixed air is colder than the supply T required for 22 C; this is a small HEATING need, NOT a cooling load' if controlled_cycle_result['Qcoil'] < 0 else ''),
                     Qcool_demand_kW_zone=controlled_cycle_result['Qcool_demand'], Qheat_equiv_kW_zone=controlled_cycle_result['Qheat_equiv'], coil_mode=('OFF / bypassed (no cooling demand)' if controlled_cycle_result['Qcoil'] <= 0 else ('wet coil + reheat' if controlled_cycle_result['Qreheat'] > 1e-9 else 'wet coil')),
                     T_room_if_no_coil_C=(air_handler_model.zone_float(condition_name, 0.30)['T_room_float'] if controlled_cycle_result['Qcoil'] <= 0 else np.nan), AHU_energy_residual_kW=controlled_cycle_result['ahu_bal'],
                     Qevap_kW=controlled_cycle_result['Qevap'], condensate_kg_h_zone=controlled_cycle_result['mcond']))
air_handler_results = pd.DataFrame(air_handler_rows); air_handler_results.to_csv(f'{RESULTS_FOLDER}/air_handler_calculation_results.csv', index=False); REPORT_DATA['ahu'] = air_handler_results.round(4).to_dict('records'); REPORT_DATA['Gov'] = GOVERNING_CONDITION; REPORT_DATA['Qdes'] = DESIGN_COOLING_LOAD_KW

                                                                                 
def build_cycle_summary_row(condition_name, evaporator_cooling_kw, condenser_temperature_k):
    if evaporator_cooling_kw <= 1e-6: return dict(Condition=condition_name, note='no cooling load (dry coil / below turndown)')
    cycle_result = absorption_cycle_model.solve_cycle(evaporator_cooling_kw, condenser_temperature_k)
    return dict(Condition=condition_name, Tcond_C=condenser_temperature_k - 273.15, Phigh_bar=cycle_result['Phigh'] / 100, Qevap_kW=evaporator_cooling_kw, Qgen_kW=cycle_result['Qgen'], Qrect_kW=cycle_result['Qrect'], Qabs_kW=cycle_result['Qabs'], Qcond_kW=cycle_result['Qcond'], Wp_kW=cycle_result['Wp'],
                COP=cycle_result['COP'], COP_rev=cycle_result['COP_rev'], mr_kg_s=cycle_result['mr'], x_strong=cycle_result['x_strong'], x_weak=cycle_result['x_weak'], m_strong=cycle_result['ms'], m_weak=cycle_result['mw'], f_circ=cycle_result['f'], Tgen_C=cycle_result['Tgen'] - 273.15,
                pump_m3_h=cycle_result['pump_m3h'], pump_head_kPa=cycle_result['pump_head_kPa'], evap_vapour_m3_h=cycle_result['evap_vap_m3h'], mfg_kg_s=cycle_result['mfg'], Tfg_out_C=cycle_result['Tfg_out'] - 273.15,
                turndown=evaporator_cooling_kw / DESIGN_COOLING_LOAD_KW, below_25pct=evaporator_cooling_kw / DESIGN_COOLING_LOAD_KW < 0.25, P_over_20bar=not cycle_result['chk']['Phigh_ok'], fg_over_80=not cycle_result['chk']['mfg_ok'], balE_rel=cycle_result['bal_E_rel'], balNH3=cycle_result['bal_NH3'])
first_pass_cycle_results = pd.DataFrame([build_cycle_summary_row(condition_name, COOLING_LOADS_30_PERCENT_FRESH_AIR[condition_name], air_handler_model.WEATHER_CONDITIONS[condition_name][0] + 273.15 + 12.0) for condition_name in WEATHER_CONDITIONS]); first_pass_cycle_results.to_csv(f'{RESULTS_FOLDER}/initial_cycle_estimates.csv', index=False); REPORT_DATA['cycle_first'] = first_pass_cycle_results.round(4).to_dict('records')

                                                     
condenser_designs = {condenser_type: condenser_model.design(condenser_type, COOLING_LOADS_30_PERCENT_FRESH_AIR) for condenser_type in ('air', 'evap')}
WET_BULB_BY_CONDITION = {condition_name: condenser_model.sink_T(condition_name, 'evap') for condition_name in WEATHER_CONDITIONS}
condenser_comparison_rows = []
for condenser_type, condenser_design in condenser_designs.items():
    condenser_comparison_rows.append(dict(type=condenser_type, governing=condenser_design['governing'], T_sink_C=condenser_design['T_sink'], Tcond_design_C=condenser_design['Tcond'], Phigh_bar=condenser_design['Phigh'] / 100, within_20bar=condenser_design['feasible_Pmax'], verdict=('ACCEPTED - meets 20 bar vessel limit' if condenser_design['feasible_Pmax'] else 'REJECTED - exceeds 20 bar vessel limit at design point'), Qcond_kW=condenser_design['Qcond'],
                         LMTD_K=condenser_design['LMTD'], UA_kW_K=condenser_design['UA'] / 1e3, area_m2=condenser_design['area'], air_kg_s=condenser_design['air_kg_s'], makeup_m3_day=condenser_design.get('makeup_m3_day', np.nan)))
condenser_comparison_results = pd.DataFrame(condenser_comparison_rows); condenser_comparison_results.to_csv(f'{RESULTS_FOLDER}/condenser_option_comparison.csv', index=False)
CHOSEN_CONDENSER_TYPE = 'evap' if condenser_designs['evap']['feasible_Pmax'] and not condenser_designs['air']['feasible_Pmax'] else ('air' if condenser_designs['air']['feasible_Pmax'] else 'evap')
SELECTED_CONDENSER_DESIGN = condenser_designs[CHOSEN_CONDENSER_TYPE]; REPORT_DATA['cond_types'] = condenser_comparison_results.round(3).to_dict('records'); REPORT_DATA['chosen'] = CHOSEN_CONDENSER_TYPE; REPORT_DATA['twb'] = {result_key: round(result_value, 2) for result_key, result_value in WET_BULB_BY_CONDITION.items()}
def build_operating_summary_row(condenser_design, condition_name, evaporator_cooling_kw):
    uncontrolled_result = condenser_model.offdesign(condenser_design, condition_name, evaporator_cooling_kw)
    if uncontrolled_result['status'] == 'no cooling load': return dict(Condition=condition_name, status='CHILLER OFF - zero cooling load')
    controlled_result = condenser_model.controlled(condenser_design, condition_name, evaporator_cooling_kw)
    return dict(Condition=condition_name, status='operating', Qevap_kW=uncontrolled_result['Qevap'], Tcond_C=uncontrolled_result['Tcond'], Phigh_bar=uncontrolled_result['Phigh_bar'], Tgen_C=uncontrolled_result['Tgen'], COP=uncontrolled_result['COP'], Qgen_kW=uncontrolled_result['Qgen'], mfg_kg_s=uncontrolled_result['mfg'],
                a_P_le_20bar=uncontrolled_result['a_Pmax'], b_Tgen_deliverable=bool(uncontrolled_result['b_gen_ok']), b_fg_floor_active=bool(uncontrolled_result['fg_floor_active']), b_eps_gen_eff=uncontrolled_result['eps_gen_eff'], b_Tfg_out_C=uncontrolled_result['Tfg_out'], c_Tcond_ge_25C=uncontrolled_result['c_Tcond'], control_action=controlled_result['action'], sink_flow_frac_ctl=controlled_result.get('sink_flow_fraction', 1.0), UA_ctl_kW_K=controlled_result['UA_kW_K'], sink_C_ctl_kW_K=controlled_result['C_kW_K'],
                Tcond_ctl_C=controlled_result['Tcond'], Phigh_ctl_bar=controlled_result['Phigh_bar'], Tgen_ctl_C=controlled_result['Tgen'], COP_ctl=controlled_result['COP'], COP_ctl_in_band_045_065=bool(0.45 <= controlled_result['COP'] <= 0.65), b_ctl_Tgen_deliverable=bool(controlled_result['b_gen_ok']), b_ctl_fg_floor_active=bool(controlled_result['fg_floor_active']), b_ctl_eps_gen_eff=controlled_result['eps_gen_eff'], Qgen_ctl_kW=controlled_result['Qgen'], mfg_ctl_kg_s=controlled_result['mfg'])
off_design_results = pd.DataFrame([build_operating_summary_row(SELECTED_CONDENSER_DESIGN, condition_name, COOLING_LOADS_30_PERCENT_FRESH_AIR[condition_name]) for condition_name in WEATHER_CONDITIONS]); off_design_results.to_csv(f'{RESULTS_FOLDER}/seasonal_operating_point_results.csv', index=False); REPORT_DATA['offdesign'] = off_design_results.round(3).to_dict('records')
closed_loop_rows = []
for condition_name in WEATHER_CONDITIONS:
    if COOLING_LOADS_30_PERCENT_FRESH_AIR[condition_name] <= 1e-6: closed_loop_rows.append(dict(Condition=condition_name, status='CHILLER OFF - zero cooling load')); continue
    controlled_cycle_result = condenser_model.controlled(SELECTED_CONDENSER_DESIGN, condition_name, COOLING_LOADS_30_PERCENT_FRESH_AIR[condition_name])['cycle']
    closed_loop_rows.append(dict(Condition=condition_name, status='operating', Tcond_C=controlled_cycle_result['Tcond'] - 273.15, Phigh_bar=controlled_cycle_result['Phigh'] / 100, Plow_bar=controlled_cycle_result['Plow'] / 100, Qevap_kW=controlled_cycle_result['Qevap'], Qgen_kW=controlled_cycle_result['Qgen'], COP=controlled_cycle_result['COP'], mr_kg_s=controlled_cycle_result['mr'],
                        x_strong=controlled_cycle_result['x_strong'], x_weak=controlled_cycle_result['x_weak'], m_strong_kg_s=controlled_cycle_result['ms'], m_weak_kg_s=controlled_cycle_result['mw'], f_circ=controlled_cycle_result['f'], Tgen_C=controlled_cycle_result['Tgen'] - 273.15, pump_m3_h=controlled_cycle_result['pump_m3h'], pump_head_kPa=controlled_cycle_result['pump_head_kPa'],
                        pump_W_kW=controlled_cycle_result['Wp'], evap_vapour_m3_h=controlled_cycle_result['evap_vap_m3h'], mfg_kg_s=controlled_cycle_result['mfg']))
pd.DataFrame(closed_loop_rows).to_csv(f'{RESULTS_FOLDER}/coupled_cycle_results.csv', index=False); REPORT_DATA['closedloop'] = pd.DataFrame(closed_loop_rows).round(4).to_dict('records')
OFF_DESIGN_BY_CONDITION = {operating_result_row['Condition']: operating_result_row for operating_result_row in REPORT_DATA['offdesign']}

                                                                                                    
seasonal_energy_rows = []
for condition_name in WEATHER_CONDITIONS:
    operating_condition_result = OFF_DESIGN_BY_CONDITION[condition_name]
    if 'Qgen_ctl_kW' in operating_condition_result and not np.isnan(operating_condition_result.get('Qgen_ctl_kW', np.nan)) and operating_condition_result['Qevap_kW'] > 1e-6 and 'turndown' not in operating_condition_result['status']:
        seasonal_energy_rows.append(dict(Condition=condition_name, Qgen_MWh_season=operating_condition_result['Qgen_ctl_kW'] * HOURS_PER_SEASON / 1e3, Qevap_MWh_season=operating_condition_result['Qevap_kW'] * HOURS_PER_SEASON / 1e3))
REPORT_DATA['annual'] = seasonal_energy_rows
REPORT_DATA['annual_total'] = dict(Qgen_GWh=sum(seasonal_condition_energy['Qgen_MWh_season'] for seasonal_condition_energy in seasonal_energy_rows) / 1e3 * 1.0, Qevap_GWh=sum(seasonal_condition_energy['Qevap_MWh_season'] for seasonal_condition_energy in seasonal_energy_rows) / 1e3)                                            

                                                    
condensate_results = pd.DataFrame([dict(Condition=condition_name, kg_h_zone=ZONES_30_PERCENT_FRESH_AIR[condition_name]['mcond'], kg_h_district=ZONES_30_PERCENT_FRESH_AIR[condition_name]['mcond'] * ZONE_COUNT, m3_day_district=ZONES_30_PERCENT_FRESH_AIR[condition_name]['cond_m3_day'] * ZONE_COUNT, kg_h_riser=ZONES_30_PERCENT_FRESH_AIR[condition_name]['mcond'] * 20) for condition_name in WEATHER_CONDITIONS])
condensate_results.to_csv(f'{RESULTS_FOLDER}/condensate_recovery_estimates.csv', index=False); REPORT_DATA['cond'] = condensate_results.round(3).to_dict('records'); REPORT_DATA['gov_decanter'] = condensate_results.loc[condensate_results.kg_h_zone.idxmax(), 'Condition']

                                                          
def district_performance_for_condition(zone_results, condition_name):
    if zone_results[condition_name]['Qevap'] <= 1e-6: return dict(Qcoil=zone_results[condition_name]['Qcoil'] * ZONE_COUNT, Qevap=0, Qgen=0, mr=0, mfg=0, cond=zone_results[condition_name]['mcond'] * ZONE_COUNT, P=np.nan, Tgen=np.nan)
    controlled_operating_result = condenser_model.controlled(SELECTED_CONDENSER_DESIGN, condition_name, zone_results[condition_name]['Qevap'])
    return dict(Qcoil=zone_results[condition_name]['Qcoil'] * ZONE_COUNT, Qevap=zone_results[condition_name]['Qevap'], Qgen=controlled_operating_result['Qgen'], mr=controlled_operating_result['cycle']['mr'], mfg=controlled_operating_result['mfg'], cond=zone_results[condition_name]['mcond'] * ZONE_COUNT, P=controlled_operating_result['Phigh_bar'], Tgen=controlled_operating_result['Tgen'])
fresh_air_comparison_results = []
for condition_name in WEATHER_CONDITIONS:
    case_30_percent, case_40_percent = district_performance_for_condition(ZONES_30_PERCENT_FRESH_AIR, condition_name), district_performance_for_condition(ZONES_40_PERCENT_FRESH_AIR, condition_name); percent_change = lambda result_key: (100 * (case_40_percent[result_key] / case_30_percent[result_key] - 1) if case_30_percent[result_key] > 1e-9 else np.nan)
    raw_cooling_coil_load_30_kw, raw_cooling_coil_load_40_kw = ZONES_30_PERCENT_FRESH_AIR[condition_name]['Qcoil'] * ZONE_COUNT, ZONES_40_PERCENT_FRESH_AIR[condition_name]['Qcoil'] * ZONE_COUNT
    fresh_air_comparison_results.append(dict(Condition=condition_name, Qcoil30_kW_raw=raw_cooling_coil_load_30_kw, Qcoil40_kW_raw=raw_cooling_coil_load_40_kw, dQcoil_abs_kW=raw_cooling_coil_load_40_kw - raw_cooling_coil_load_30_kw, dQcoil_pct=(100 * (raw_cooling_coil_load_40_kw / raw_cooling_coil_load_30_kw - 1) if abs(raw_cooling_coil_load_30_kw) > 1e-9 else np.nan), dQevap_pct=percent_change('Qevap'), dQgen_pct=percent_change('Qgen'), dmr_pct=percent_change('mr'), dmfg_pct=percent_change('mfg'), dcond_pct=percent_change('cond'),
                   Qevap30=case_30_percent['Qevap'], Qevap40=case_40_percent['Qevap'], Qgen30=case_30_percent['Qgen'], Qgen40=case_40_percent['Qgen'], mfg40=case_40_percent['mfg'], Phigh40_bar=case_40_percent['P'], RH_room40=100 * ZONES_40_PERCENT_FRESH_AIR[condition_name]['RH'],
                   ADP40=ZONES_40_PERCENT_FRESH_AIR[condition_name]['ADP'], reheat40_kW_zone=ZONES_40_PERCENT_FRESH_AIR[condition_name]['Qreheat']))
fresh_air_comparison_results = pd.DataFrame(fresh_air_comparison_results); fresh_air_comparison_results.to_csv(f'{RESULTS_FOLDER}/forty_percent_fresh_air_case.csv', index=False); REPORT_DATA['fresh40'] = fresh_air_comparison_results.round(3).to_dict('records')
monsoon_night_uncontrolled_40_percent = air_handler_model.zone('Monsoon night', 0.40); monsoon_night_controlled_40_percent = ZONES_40_PERCENT_FRESH_AIR['Monsoon night']; monsoon_night_controlled_30_percent = ZONES_30_PERCENT_FRESH_AIR['Monsoon night']
uncontrolled_cooling_kw = monsoon_night_uncontrolled_40_percent['Qevap']; controlled_cooling_kw = monsoon_night_controlled_40_percent['Qevap']; uncontrolled_generator_heat_kw = condenser_model.controlled(SELECTED_CONDENSER_DESIGN, 'Monsoon night', uncontrolled_cooling_kw)['Qgen']; controlled_generator_heat_kw = condenser_model.controlled(SELECTED_CONDENSER_DESIGN, 'Monsoon night', controlled_cooling_kw)['Qgen']
REPORT_DATA['mn40_rhfix'] = dict(ADP_uncontrolled=monsoon_night_uncontrolled_40_percent['ADP'], ADP_controlled=monsoon_night_controlled_40_percent['ADP'], dADP_K=monsoon_night_controlled_40_percent['ADP'] - monsoon_night_uncontrolled_40_percent['ADP'], reheat_kW_district=monsoon_night_controlled_40_percent['Qreheat'] * ZONE_COUNT, Qevap_uncontrolled_kW=uncontrolled_cooling_kw, Qevap_controlled_kW=controlled_cooling_kw, extra_Qevap_kW=controlled_cooling_kw - uncontrolled_cooling_kw,
                       Qgen_uncontrolled_kW=uncontrolled_generator_heat_kw, Qgen_controlled_kW=controlled_generator_heat_kw, extra_Qgen_kW=controlled_generator_heat_kw - uncontrolled_generator_heat_kw, extra_Qgen_pct=100 * (controlled_generator_heat_kw / uncontrolled_generator_heat_kw - 1))
pd.DataFrame([REPORT_DATA['mn40_rhfix']]).to_csv(f'{RESULTS_FOLDER}/monsoon_night_humidity_adjustment.csv', index=False)
REPORT_DATA['mn40'] = dict(RH_uncontrolled=100 * monsoon_night_uncontrolled_40_percent['RH'], violates=bool(monsoon_night_uncontrolled_40_percent['viol']), ADP_ctl=monsoon_night_controlled_40_percent['ADP'], reheat_kW_district=monsoon_night_controlled_40_percent['Qreheat'] * ZONE_COUNT, Qevap30=monsoon_night_controlled_30_percent['Qevap'], Qevap40=monsoon_night_controlled_40_percent['Qevap'],
                 extra_Qgen_kW_total_load_growth_30to40=fresh_air_comparison_results.loc[fresh_air_comparison_results.Condition == 'Monsoon night', 'Qgen40'].item() - fresh_air_comparison_results.loc[fresh_air_comparison_results.Condition == 'Monsoon night', 'Qgen30'].item())
REPORT_DATA['limits40'] = dict(max_mfg=float(fresh_air_comparison_results.mfg40.max()), mfg_limit_active=bool((fresh_air_comparison_results.mfg40 > 80).any()), max_P_bar=float(np.nanmax(fresh_air_comparison_results.Phigh40_bar)), P_limit_active=bool((fresh_air_comparison_results.Phigh40_bar > 20).any()))

                                                                                                       
MAX_FLUE_GAS_KG_S = 80.0; TARGET_FLUE_GAS_KG_S = 76.0; flue_gas_remedy_rows = []; OPERATING_HOURS_PER_DAY = 12.0                                                                      
for condition_name in WEATHER_CONDITIONS:
    COOLING_LOAD_40_PERCENT = ZONES_40_PERCENT_FRESH_AIR[condition_name]['Qevap']
    if COOLING_LOAD_40_PERCENT <= 1e-6: continue
    if condenser_model.controlled(SELECTED_CONDENSER_DESIGN, condition_name, COOLING_LOAD_40_PERCENT)['mfg'] <= MAX_FLUE_GAS_KG_S: continue
    MAX_COOLING_AT_FLUE_LIMIT_KW = condenser_model.qevap_at_fg_limit(SELECTED_CONDENSER_DESIGN, condition_name); COOLING_SHORTFALL_KW = COOLING_LOAD_40_PERCENT - MAX_COOLING_AT_FLUE_LIMIT_KW; NIGHT_PARTNER_CONDITION = condition_name.replace('day', 'night'); NIGHT_COOLING_KW = ZONES_40_PERCENT_FRESH_AIR[NIGHT_PARTNER_CONDITION]['Qevap']; NIGHT_MAX_COOLING_KW = condenser_model.qevap_at_fg_limit(SELECTED_CONDENSER_DESIGN, NIGHT_PARTNER_CONDITION)
    CONDENSER_SIZE_FACTOR = condenser_model.size_for_fg_limit(SELECTED_CONDENSER_DESIGN, condition_name, COOLING_LOAD_40_PERCENT, maximum_flue_gas_mass_flow_kg_s=TARGET_FLUE_GAS_KG_S); RESIZED_OPERATING_RESULT = condenser_model.offdesign(SELECTED_CONDENSER_DESIGN, condition_name, COOLING_LOAD_40_PERCENT, size_factor=CONDENSER_SIZE_FACTOR); STORAGE_ENERGY_KWH = COOLING_SHORTFALL_KW * OPERATING_HOURS_PER_DAY
    flue_gas_remedy_rows.append(dict(Condition=condition_name, Qevap40_kW=COOLING_LOAD_40_PERCENT, Qevap_max_at_80kg_s_kW=MAX_COOLING_AT_FLUE_LIMIT_KW, shortfall_kW=COOLING_SHORTFALL_KW, shortfall_pct=100 * COOLING_SHORTFALL_KW / COOLING_LOAD_40_PERCENT,
                    optB_storage_MWh_th=STORAGE_ENERGY_KWH / 1e3, optB_chilledwater_m3_dT5K=STORAGE_ENERGY_KWH * 3600 / (4.186 * 5) / 1e3, night_partner=NIGHT_PARTNER_CONDITION, night_spare_flue_capacity_kW=NIGHT_MAX_COOLING_KW - NIGHT_COOLING_KW,
                    optB_charge_feasible=bool((NIGHT_MAX_COOLING_KW - NIGHT_COOLING_KW) * OPERATING_HOURS_PER_DAY >= 1.1 * STORAGE_ENERGY_KWH), optA_size_factor=CONDENSER_SIZE_FACTOR, optA_area_m2=SELECTED_CONDENSER_DESIGN['area'] * CONDENSER_SIZE_FACTOR, optA_extra_area_m2=SELECTED_CONDENSER_DESIGN['area'] * (CONDENSER_SIZE_FACTOR - 1),
                    optA_UA_kW_K=SELECTED_CONDENSER_DESIGN['UA'] / 1e3 * CONDENSER_SIZE_FACTOR, optA_sink_C_kW_K=SELECTED_CONDENSER_DESIGN['C'] * CONDENSER_SIZE_FACTOR, optA_sink_flow_scale=CONDENSER_SIZE_FACTOR, optA_Tcond_C=RESIZED_OPERATING_RESULT['Tcond'], optA_Phigh_bar=RESIZED_OPERATING_RESULT['Phigh_bar'], optA_mfg_kg_s=RESIZED_OPERATING_RESULT['mfg'], optA_target_mfg=TARGET_FLUE_GAS_KG_S))
flue_gas_remedy_results = pd.DataFrame(flue_gas_remedy_rows)
flue_gas_remedy_results.insert(0, 'basis', 'MITIGATED design (resized condenser) vs UNMITIGATED 30 %-fresh condenser; see resized_condenser_feasibility.csv for every condition')
flue_gas_remedy_results.to_csv(f'{RESULTS_FOLDER}/exhaust_gas_limit_solutions.csv', index=False); REPORT_DATA['flue_remedy'] = flue_gas_remedy_results.round(3).to_dict('records')
FINAL_CONDENSER_SIZE_FACTOR = float(flue_gas_remedy_results.optA_size_factor.max()) if len(flue_gas_remedy_results) else 1.0                                                                                              
mitigation_rows = []
for condition_name in WEATHER_CONDITIONS:
    COOLING_LOAD_40_PERCENT = ZONES_40_PERCENT_FRESH_AIR[condition_name]['Qevap']
    if COOLING_LOAD_40_PERCENT <= 1e-6: mitigation_rows.append(dict(Condition=condition_name, status='chiller OFF - zero cooling load')); continue
    UNMITIGATED_RESULT, MITIGATED_RESULT = condenser_model.controlled(SELECTED_CONDENSER_DESIGN, condition_name, COOLING_LOAD_40_PERCENT), condenser_model.controlled(SELECTED_CONDENSER_DESIGN, condition_name, COOLING_LOAD_40_PERCENT, size_factor=FINAL_CONDENSER_SIZE_FACTOR)
    mitigation_rows.append(dict(Condition=condition_name, status='operating', Qevap40_kW=COOLING_LOAD_40_PERCENT, size_factor=FINAL_CONDENSER_SIZE_FACTOR, condenser_area_m2=SELECTED_CONDENSER_DESIGN['area'] * FINAL_CONDENSER_SIZE_FACTOR, design_UA_kW_K=SELECTED_CONDENSER_DESIGN['UA'] / 1e3 * FINAL_CONDENSER_SIZE_FACTOR, design_sink_C_kW_K=SELECTED_CONDENSER_DESIGN['C'] * FINAL_CONDENSER_SIZE_FACTOR,
                    UNMIT_Tcond_C=UNMITIGATED_RESULT['Tcond'], UNMIT_Phigh_bar=UNMITIGATED_RESULT['Phigh_bar'], UNMIT_flue_kg_s=UNMITIGATED_RESULT['mfg'], UNMIT_feasible_80=bool(UNMITIGATED_RESULT['mfg'] <= MAX_FLUE_GAS_KG_S + 1e-3),
                    MIT_Tcond_C=MITIGATED_RESULT['Tcond'], MIT_Phigh_bar=MITIGATED_RESULT['Phigh_bar'], MIT_sink_sink_flow_fraction=MITIGATED_RESULT['sink_flow_fraction'], MIT_UA_kW_K=MITIGATED_RESULT['UA_kW_K'], MIT_Qgen_kW=MITIGATED_RESULT['Qgen'], MIT_COP=MITIGATED_RESULT['COP'], MIT_flue_kg_s=MITIGATED_RESULT['mfg'],
                    MIT_feasible_80=bool(MITIGATED_RESULT['mfg'] <= MAX_FLUE_GAS_KG_S + 1e-3), MIT_P_le_20bar=bool(MITIGATED_RESULT['Phigh_bar'] <= 20.0), MIT_Tcond_ge_25C=bool(MITIGATED_RESULT['Tcond'] >= 25.0 - 1e-6)))
resized_condenser_results = pd.DataFrame(mitigation_rows); resized_condenser_results.to_csv(f'{RESULTS_FOLDER}/resized_condenser_feasibility.csv', index=False); REPORT_DATA['mitigated'] = resized_condenser_results.round(4).to_dict('records'); REPORT_DATA['size_final'] = FINAL_CONDENSER_SIZE_FACTOR
                                                                             
                                                                                 
                                                                           
recommended_seasonal_rows = []
for fresh_fraction, zone_results in ((0.30, ZONES_30_PERCENT_FRESH_AIR), (0.40, ZONES_40_PERCENT_FRESH_AIR)):
    for condition_name in WEATHER_CONDITIONS:
        cooling_duty_kw = zone_results[condition_name]['Qevap']
        if cooling_duty_kw <= 1e-6:
            recommended_seasonal_rows.append(dict(fresh_air_fraction=fresh_fraction, Condition=condition_name, status='chiller OFF - zero cooling load', Qevap_kW=0.0))
            continue
        seasonal_result = condenser_model.controlled(SELECTED_CONDENSER_DESIGN, condition_name, cooling_duty_kw, size_factor=FINAL_CONDENSER_SIZE_FACTOR)
        recommended_seasonal_rows.append(dict(fresh_air_fraction=fresh_fraction, Condition=condition_name, status=seasonal_result['status'], condenser_type=SELECTED_CONDENSER_DESIGN['type'], condenser_size_factor=FINAL_CONDENSER_SIZE_FACTOR,
            condenser_area_m2=SELECTED_CONDENSER_DESIGN['area'] * FINAL_CONDENSER_SIZE_FACTOR, Qevap_kW=cooling_duty_kw, Tcond_C=seasonal_result['Tcond'], Phigh_bar=seasonal_result['Phigh_bar'], COP=seasonal_result['COP'], Qgen_kW=seasonal_result['Qgen'], Qcond_kW=seasonal_result['Qcond'], flue_gas_kg_s=seasonal_result['mfg'], sink_flow_fraction=seasonal_result['sink_flow_fraction'],
            flue_limit_ok=bool(seasonal_result['mfg'] <= MAX_FLUE_GAS_KG_S + 1e-3), pressure_limit_ok=bool(seasonal_result['Phigh_bar'] <= 20.0), generator_delivery_ok=bool(seasonal_result['b_gen_ok'])))
recommended_seasonal_frame = pd.DataFrame(recommended_seasonal_rows)
recommended_seasonal_frame.to_csv(f'{RESULTS_FOLDER}/recommended_condenser_seasonal_results.csv', index=False)
active_seasonal = recommended_seasonal_frame[recommended_seasonal_frame['status'] != 'chiller OFF - zero cooling load'].copy()
figure, plot_axes = plt.subplots(1, 2, figsize=(12, 4.5))
for fresh_fraction in (0.30, 0.40):
    fresh_rows = active_seasonal[active_seasonal.fresh_air_fraction == fresh_fraction]
    plot_axes[0].plot(fresh_rows.Condition, fresh_rows.COP, marker='o', label=f'{fresh_fraction:.0%} fresh air')
    plot_axes[1].plot(fresh_rows.Condition, fresh_rows.flue_gas_kg_s, marker='o', label=f'{fresh_fraction:.0%} fresh air')
plot_axes[0].set_ylabel('Thermal COP'); plot_axes[0].set_title('COP — recommended resized condenser')
plot_axes[1].axhline(MAX_FLUE_GAS_KG_S, linestyle='--', label='80 kg/s limit'); plot_axes[1].set_ylabel('Flue-gas flow (kg/s)'); plot_axes[1].set_title('Waste-heat feasibility')
for axis in plot_axes: axis.tick_params(axis='x', rotation=35); axis.grid(True, alpha=0.25); axis.legend()
figure.tight_layout(); figure.savefig(f'{RESULTS_FOLDER}/recommended_condenser_seasonal_performance.png', dpi=160, bbox_inches='tight'); plt.close(figure)
REPORT_DATA['recommended_seasonal'] = recommended_seasonal_frame.round(4).to_dict('records')
flue_gas_summary_rows = []
for condition_name in WEATHER_CONDITIONS:
    first_pass_flue_gas = first_pass_cycle_results.loc[first_pass_cycle_results.Condition == condition_name, 'mfg_kg_s']; first_pass_flue_gas = float(first_pass_flue_gas.iloc[0]) if len(first_pass_flue_gas) and 'mfg_kg_s' in first_pass_cycle_results and not pd.isna(first_pass_flue_gas.iloc[0]) else np.nan
    flue_gas_status_row = dict(Condition=condition_name, first_pass_Tdb_plus12=first_pass_flue_gas, final_30pct=OFF_DESIGN_BY_CONDITION[condition_name].get('mfg_ctl_kg_s', np.nan), final_40pct_asis=np.nan, final_40pct_upsized=np.nan)                                                                                             
    if ZONES_40_PERCENT_FRESH_AIR[condition_name]['Qevap'] > 1e-6:
        flue_gas_status_row['final_40pct_asis'] = condenser_model.controlled(SELECTED_CONDENSER_DESIGN, condition_name, ZONES_40_PERCENT_FRESH_AIR[condition_name]['Qevap'])['mfg']
        flue_gas_status_row['final_40pct_upsized'] = condenser_model.controlled(SELECTED_CONDENSER_DESIGN, condition_name, ZONES_40_PERCENT_FRESH_AIR[condition_name]['Qevap'], size_factor=FINAL_CONDENSER_SIZE_FACTOR)['mfg']                                                  
    flue_gas_summary_rows.append(flue_gas_status_row)
flue_gas_summary_results = pd.DataFrame(flue_gas_summary_rows); flue_gas_summary_results['limit_80_ok_final30'] = flue_gas_summary_results.final_30pct.fillna(0) <= MAX_FLUE_GAS_KG_S; flue_gas_summary_results['limit_80_ok_40_asis'] = flue_gas_summary_results.final_40pct_asis.fillna(0) <= MAX_FLUE_GAS_KG_S; flue_gas_summary_results['limit_80_ok_40_upsized'] = flue_gas_summary_results.final_40pct_upsized.fillna(0) <= MAX_FLUE_GAS_KG_S + 1e-3
flue_gas_summary_results['status_40pct'] = ['n/a - zero load' if pd.isna(unmitigated_flue_gas_mass_flow_kg_s) else ('INFEASIBLE before mitigation (%.1f kg/s > 80); FEASIBLE after condenser resizing (%.1f kg/s) - two different designs' % (unmitigated_flue_gas_mass_flow_kg_s, resized_flue_gas_mass_flow_kg_s) if unmitigated_flue_gas_mass_flow_kg_s > MAX_FLUE_GAS_KG_S else 'feasible as designed (%.1f kg/s <= 80), no resizing needed' % unmitigated_flue_gas_mass_flow_kg_s) for unmitigated_flue_gas_mass_flow_kg_s, resized_flue_gas_mass_flow_kg_s in zip(flue_gas_summary_results.final_40pct_asis, flue_gas_summary_results.final_40pct_upsized)]
flue_gas_summary_results = flue_gas_summary_results.rename(columns={'final_40pct_asis': '40pct_UNMITIGATED_30pct_design_condenser', 'final_40pct_upsized': '40pct_MITIGATED_resized_condenser', 'limit_80_ok_40_asis': 'limit_80_ok_40_unmitigated', 'limit_80_ok_40_upsized': 'limit_80_ok_40_mitigated'})
flue_gas_summary_results.to_csv(f'{RESULTS_FOLDER}/exhaust_gas_feasibility_summary.csv', index=False); REPORT_DATA['fg_status'] = flue_gas_summary_results.round(3).to_dict('records')

                                                                                                  
turndown_rows = []
for condition_name in WEATHER_CONDITIONS:
    for fresh_air_case_label, zone_results in (('30 % fresh', ZONES_30_PERCENT_FRESH_AIR), ('40 % fresh', ZONES_40_PERCENT_FRESH_AIR)):
        evaporator_cooling_kw = zone_results[condition_name]['Qevap']; design_load_fraction = evaporator_cooling_kw / DESIGN_COOLING_LOAD_KW
        turndown_rows.append(dict(Condition=condition_name, case=fresh_air_case_label, Qevap_kW=evaporator_cooling_kw, load_frac_of_design=design_load_fraction, state='OFF - zero cooling load (shutdown / free cooling)' if evaporator_cooling_kw <= 1e-6 else ('OPERATING at >= 25 % turndown' if design_load_fraction >= 0.25 else 'BELOW 25 % turndown (storage / staging needed)')))
pd.DataFrame(turndown_rows).to_csv(f'{RESULTS_FOLDER}/minimum_load_operating_status.csv', index=False); REPORT_DATA['turndown'] = pd.DataFrame(turndown_rows).round(3).to_dict('records')
minimum_load_check_rows = []
for condition_name in WEATHER_CONDITIONS:
    controlled_operating_result = condenser_model.controlled(SELECTED_CONDENSER_DESIGN, condition_name, 0.25 * DESIGN_COOLING_LOAD_KW); minimum_load_check_rows.append(dict(Condition=condition_name, Qevap_min_stable_kW=0.25 * DESIGN_COOLING_LOAD_KW, Tcond_ctl_C=controlled_operating_result['Tcond'], sink_sink_flow_fraction=controlled_operating_result.get('sink_flow_fraction', 1.0), COP=controlled_operating_result['COP'], Qgen_kW=controlled_operating_result['Qgen'], mfg_kg_s=controlled_operating_result['mfg'], action=controlled_operating_result['action']))
pd.DataFrame(minimum_load_check_rows).to_csv(f'{RESULTS_FOLDER}/hypothetical_minimum_load_check.csv', index=False); REPORT_DATA['turndown_floor_test'] = pd.DataFrame(minimum_load_check_rows).round(3).to_dict('records')

                                                                 
concentration_sensitivity_rows = []
for concentration_difference in (0.10, 0.12, 0.16, 0.20):
    for condenser_temperature_c in (25.0, 35.0, 40.0, 45.0):
        cycle_result = absorption_cycle_model.solve_cycle(5000.0, condenser_temperature_c + 273.15, concentration_difference=concentration_difference)
        concentration_sensitivity_rows.append(dict(concentration_difference=concentration_difference, Tcond_C=condenser_temperature_c, COP=cycle_result['COP'], Tgen_C=cycle_result['Tgen'] - 273.15, circulation_ratio=cycle_result['f'], Qgen_kW=cycle_result['Qgen']))
pd.DataFrame(concentration_sensitivity_rows).to_csv(f'{RESULTS_FOLDER}/solution_concentration_sensitivity.csv', index=False); REPORT_DATA['sens'] = pd.DataFrame(concentration_sensitivity_rows).round(3).to_dict('records')

                                           
PLOT_LINE_COLOURS = ['#2457A7', '#0F8B8D', '#D97706', '#C2414B', '#7656A8', '#475569']
figure, plot_axes = plt.subplots(figsize=(8, 6)); dry_bulb_temperature_values_c = np.linspace(8, 45, 200)
for relative_humidity_fraction in (.2, .4, .6, .8, 1.0): plot_axes.plot(dry_bulb_temperature_values_c, [1000 * air_handler_model.w_from_RH(dry_bulb_temperature_c, relative_humidity_fraction) for dry_bulb_temperature_c in dry_bulb_temperature_values_c], 'gray' if relative_humidity_fraction < 1 else 'k', lw=1 if relative_humidity_fraction < 1 else 2)
plot_axes.plot(dry_bulb_temperature_values_c, [1000 * air_handler_model.w_from_RH(dry_bulb_temperature_c, .75) if False else 1000 * air_handler_model.w_from_RH(22, .75) * 0 + 0 for dry_bulb_temperature_c in dry_bulb_temperature_values_c], alpha=0)
for condition_name, plot_color in zip(WEATHER_CONDITIONS, PLOT_LINE_COLOURS):
    zone_air_result = ZONES_30_PERCENT_FRESH_AIR[condition_name]; room_temperature_c = 22.0
    plot_axes.plot([zone_air_result['To'], zone_air_result['Tm'], room_temperature_c], [1000 * zone_air_result['wo'], 1000 * zone_air_result['wm'], 1000 * zone_air_result['wr']], ':', c=plot_color)
    plot_axes.plot(zone_air_result['To'], 1000 * zone_air_result['wo'], 'o', c=plot_color, label=condition_name); plot_axes.plot(zone_air_result['Tm'], 1000 * zone_air_result['wm'], 's', c=plot_color); plot_axes.plot(room_temperature_c, 1000 * zone_air_result['wr'], '*', c=plot_color, ms=10)
    if zone_air_result['wet']:
        plot_axes.plot([zone_air_result['Tm'], zone_air_result['ADP']], [1000 * zone_air_result['wm'], 1000 * air_handler_model.wsat(zone_air_result['ADP'])], '-', c=plot_color); plot_axes.plot(zone_air_result['ADP'], 1000 * air_handler_model.wsat(zone_air_result['ADP']), 'x', c=plot_color); plot_axes.plot(zone_air_result['Tcoil'], 1000 * zone_air_result['ws'], '^', c=plot_color)
plot_axes.set_xlim(8, 45); plot_axes.set_ylim(0, 34); plot_axes.set_xlabel('Dry-bulb T [C]'); plot_axes.set_ylabel('w [g/kg da]'); plot_axes.legend(fontsize=7)
plot_axes.set_title('Psychrometric chart: O(o) - M(s) - coil to ADP(x), S(^), R(*)  [30 % fresh air, ceiling-controlled]'); figure.savefig(f'{RESULTS_FOLDER}/air_moisture_process_chart.png', dpi=150, bbox_inches='tight'); plt.close(figure)
figure, plot_axes = plt.subplots(1, 2, figsize=(11, 4)); condition_index = np.arange(6)
plot_axes[0].bar(condition_index - .2, [COOLING_LOADS_30_PERCENT_FRESH_AIR[condition_name] for condition_name in WEATHER_CONDITIONS], .4, label='Q_evap'); plot_axes[0].bar(condition_index + .2, [OFF_DESIGN_BY_CONDITION[condition_name].get('Qgen_ctl_kW', 0) if 'Qgen_ctl_kW' in OFF_DESIGN_BY_CONDITION[condition_name] and OFF_DESIGN_BY_CONDITION[condition_name]['Qevap_kW'] > 0 and 'turndown' not in OFF_DESIGN_BY_CONDITION[condition_name]['status'] else 0 for condition_name in WEATHER_CONDITIONS], .4, label='Q_gen (off-design)')
plot_axes[0].set_xticks(condition_index); plot_axes[0].set_xticklabels(WEATHER_CONDITIONS, rotation=30, fontsize=8); plot_axes[0].set_ylabel('kW'); plot_axes[0].legend(); plot_axes[0].set_title('Seasonal duties (evaporative condenser)')
plot_axes[1].plot(condition_index, [OFF_DESIGN_BY_CONDITION[condition_name].get('Tcond_C', np.nan) for condition_name in WEATHER_CONDITIONS], 'o-', label='T_cond [C]'); plot_axes[1].axhline(25, c='r', ls='--', label='25 C valve limit'); secondary_axis = plot_axes[1].twinx(); secondary_axis.plot(condition_index, [OFF_DESIGN_BY_CONDITION[condition_name].get('COP', np.nan) for condition_name in WEATHER_CONDITIONS], 's--', c='#7656A8'); secondary_axis.set_ylabel('COP')
plot_axes[1].set_xticks(condition_index); plot_axes[1].set_xticklabels(WEATHER_CONDITIONS, rotation=30, fontsize=8); plot_axes[1].legend(); plot_axes[1].set_title('Off-design map (UA fixed)'); figure.savefig(f'{RESULTS_FOLDER}/seasonal_chiller_performance.png', dpi=150, bbox_inches='tight'); plt.close(figure)
figure, plot_axes = plt.subplots(figsize=(6, 4.2)); ammonia_mole_fraction_values = np.linspace(0.1, 0.93, 45); import equilibrium as equilibrium_model
for pressure_kpa, plot_color in ((500, '#2457A7'), (1100, '#0F8B8D'), (1700, '#C2414B')):
    bubble_temperature_values_c = [solution_properties.T_bub(pressure_kpa, ammonia_mole_fraction) - 273.15 for ammonia_mole_fraction in ammonia_mole_fraction_values]; plot_axes.plot(ammonia_mole_fraction_values, bubble_temperature_values_c, c=plot_color, label='bubble %d kPa (Eq.6)' % pressure_kpa)
    plot_axes.plot([equilibrium_model.y_of_px(pressure_kpa / 1e3, ammonia_mole_fraction) for ammonia_mole_fraction in ammonia_mole_fraction_values], bubble_temperature_values_c, '--', c=plot_color, label='dew %d kPa (Eq.6+7)' % pressure_kpa)
plot_axes.set_xlabel('NH3 mole fraction: x (solid) / y (dashed)'); plot_axes.set_ylabel('T [C]'); plot_axes.legend(fontsize=7); plot_axes.set_title('Patek-Klomfar T-x-y (own implementation)'); figure.savefig(f'{RESULTS_FOLDER}/vapour_liquid_equilibrium_check.png', dpi=150, bbox_inches='tight'); plt.close(figure)
                                                                         
                                                                             
                                                                             
figure, plot_axes = plt.subplots(figsize=(7.2, 4.6))
liquid_composition_grid = np.linspace(0.08, 0.90, 100)
for pressure_kpa, plot_color in ((500, '#2457A7'), (1000, '#0F8B8D'), (1500, '#C2414B')):
    pressure_mpa = pressure_kpa / 1000.0
    vapour_composition_eq8 = np.array([
        equilibrium_model.y_of_px_eq8(pressure_mpa, x) for x in liquid_composition_grid
    ])
    vapour_composition_eq67 = np.array([
        equilibrium_model.y_of_px(pressure_mpa, x) for x in liquid_composition_grid
    ])
    relative_vapour_composition_error_percent = 100.0 * (
        vapour_composition_eq8 - vapour_composition_eq67
    ) / vapour_composition_eq67
    plot_axes.plot(
        liquid_composition_grid, relative_vapour_composition_error_percent,
        color=plot_color, label=f'{pressure_kpa} kPa',
    )
plot_axes.axhline(0.0, color='0.25', linestyle='--', linewidth=1.0, label='Eq. (8) matches Eq. (6)+(7)')
plot_axes.axvline(0.2, color='0.45', linestyle=':', linewidth=1.0, label='x = 0.2')
plot_axes.set_xlabel('Liquid ammonia mole fraction, x')
plot_axes.set_ylabel('Relative y error of Eq. (8) vs Eqs. (6)+(7) [%]')
plot_axes.set_title('Direct Equation (8) correlation consistency check')
plot_axes.grid(True, alpha=0.3)
plot_axes.legend()
figure.tight_layout()
figure.savefig(f'{RESULTS_FOLDER}/equilibrium_eq8_consistency.png', dpi=160, bbox_inches='tight')
plt.close(figure)
REPORT_DATA['cycle_states_example'] = {result_key: (result_value if not isinstance(result_value, np.generic) else float(result_value)) for result_key, result_value in absorption_cycle_model.solve_cycle(DESIGN_COOLING_LOAD_KW, 313.65)['h'].items()}
json.dump(REPORT_DATA, open(f'{RESULTS_FOLDER}/results.json', 'w'), indent=1, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else bool(o) if isinstance(o, np.bool_) else str(o))
print('Governing chiller condition:', GOVERNING_CONDITION, '%.0f kW' % DESIGN_COOLING_LOAD_KW, '| decanter:', REPORT_DATA['gov_decanter'], '| chosen condenser:', CHOSEN_CONDENSER_TYPE)
for table_label, results_dataframe in (('Flue-gas remedy (40 % fresh)', flue_gas_remedy_results), ('Flue-gas status', flue_gas_summary_results), ('Turndown status', pd.DataFrame(turndown_rows)), ('Task1 AHU', air_handler_results), ('Task1 cycle (Tdb+12)', first_pass_cycle_results), ('Task2 types', condenser_comparison_results), ('Task2 off-design', off_design_results), ('Task3', condensate_results), ('Task4', fresh_air_comparison_results)):
    print('\n==', table_label); print(results_dataframe.round(3).to_string())
print('\nMonsoon night, 40 % fresh:', REPORT_DATA['mn40']); print('Limits 40 %:', REPORT_DATA['limits40'])
