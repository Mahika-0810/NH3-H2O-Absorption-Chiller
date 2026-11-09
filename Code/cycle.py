"""Ammonia-water absorption chiller calculations.

The cycle is calculated in stages: refrigerant flow, solution circulation,
pump work, solution heat exchange, and the generator/rectifier.  The same
property correlations and assumptions are used throughout the model.

Temperatures are in kelvin, pressure is in kPa, heat and power are in kW,
flow rates are in kg/s, and enthalpy is in kJ/kg.

State numbers used by the model:
1 rectifier vapour; 2 condenser liquid; 3 expansion-valve outlet; 4 evaporator
vapour; 5 strong solution from absorber; 6 pumped strong solution; 7 strong
solution after the solution heat exchanger; 8 weak solution from generator;
9 weak solution after the solution heat exchanger; 10 throttled weak solution.
"""
from __future__ import annotations
import numpy as np
import nh3_props as ammonia_properties
import solution_props as solution_properties

PROPERTY_MODEL_NAME = solution_properties.PROPERTY_MODEL_DESCRIPTION
DEFAULT_EVAPORATOR_TEMPERATURE_K = 275.15
RECTIFIER_OUTLET_AMMONIA_MASS_FRACTION = 0.995
PUMP_EFFICIENCY = 0.70
SOLUTION_HEAT_EXCHANGER_EFFECTIVENESS = 0.70
GENERATOR_HEAT_EXCHANGER_EFFECTIVENESS = 0.75
FLUE_GAS_SPECIFIC_HEAT_KJ_KG_K = 1.10
FLUE_GAS_INLET_TEMPERATURE_K = 573.15
MINIMUM_FLUE_GAS_OUTLET_TEMPERATURE_K = 413.15
MAXIMUM_FLUE_GAS_MASS_FLOW_KG_S = 80.0
MAXIMUM_HIGH_SIDE_PRESSURE_KPA = 2000.0
DEFAULT_CONCENTRATION_DIFFERENCE = 0.16
CONDENSER_TEMPERATURE_TOLERANCE_K = 1e-6


def _validate_inputs(cooling_duty_kw, condenser_temperature_k, concentration_difference, evaporator_temperature_k):
    named_values = {
        "evaporator_cooling_kw": cooling_duty_kw,
        "condenser_temperature_k": condenser_temperature_k,
        "concentration_difference": concentration_difference,
        "evaporator_temperature_k": evaporator_temperature_k,
    }
    for name, value in named_values.items():
        if not np.isfinite(value):
            raise ValueError(f"{name} must be a finite number; received {value!r}.")
    if cooling_duty_kw <= 0:
        raise ValueError(f"evaporator_cooling_kw must be greater than 0 kW; received {cooling_duty_kw}.")
    if evaporator_temperature_k <= 0 or condenser_temperature_k <= 0:
        raise ValueError("condenser_temperature_k and evaporator_temperature_k must be absolute temperatures in kelvin.")
    if condenser_temperature_k <= evaporator_temperature_k:
        raise ValueError(f"condenser_temperature_k ({condenser_temperature_k} K) must be higher than evaporator_temperature_k ({evaporator_temperature_k} K).")
    if not 0 < concentration_difference < 1:
        raise ValueError(f"concentration_difference must be between 0 and 1; received {concentration_difference}.")
    if condenser_temperature_k >= ammonia_properties.CRITICAL_TEMPERATURE_K:
        raise ValueError(f"condenser_temperature_k must be below ammonia's critical temperature ({ammonia_properties.CRITICAL_TEMPERATURE_K} K).")


def _refrigerant_side(cooling_duty_kw, condenser_temperature_k, evaporator_temperature_k):
    condenser_saturation = ammonia_properties.sat(condenser_temperature_k)
    evaporator_saturation = ammonia_properties.sat(evaporator_temperature_k)
    enthalpy_2 = solution_properties.h_NH3_liq(condenser_temperature_k)

    enthalpy_3 = enthalpy_2

    enthalpy_4 = solution_properties.h_NH3_vap(evaporator_temperature_k)

    specific_cooling = enthalpy_4 - enthalpy_3
    if not np.isfinite(specific_cooling) or specific_cooling <= 0:
        raise ValueError(f"The refrigerant specific cooling effect must be positive; got {specific_cooling:.6g} kJ/kg.")
    refrigerant_flow = cooling_duty_kw / specific_cooling
    return {
        "condenser_saturation": condenser_saturation,
        "evaporator_saturation": evaporator_saturation,
        "h2": enthalpy_2, "h3": enthalpy_3, "h4": enthalpy_4,
        "specific_cooling": specific_cooling, "refrigerant_flow": refrigerant_flow,
    }


def _solution_flow_and_pump_side(high_pressure_kpa, low_pressure_kpa, condenser_temperature_k, concentration_difference, refrigerant_flow):
    strong_fraction = solution_properties.mole_fraction_to_mass_fraction(
        solution_properties.x_bub(low_pressure_kpa, condenser_temperature_k)
    )
    weak_fraction = strong_fraction - concentration_difference
    if not 0 < weak_fraction < strong_fraction < 1:
        raise ValueError(
            "The selected degassing width gives a non-physical solution concentration: "
            f"x_strong={strong_fraction:.6g}, x_weak={weak_fraction:.6g}, "
            f"concentration_difference={concentration_difference:.6g}. Choose a smaller positive dx."
        )
    enthalpy_5 = solution_properties.h_liquid(low_pressure_kpa, strong_fraction, condenser_temperature_k)
    generator_temperature = solution_properties.T_bub(
        high_pressure_kpa, solution_properties.mass_fraction_to_mole_fraction(weak_fraction)
    )
    enthalpy_8 = solution_properties.h_liquid(high_pressure_kpa, weak_fraction, generator_temperature)

    circulation_ratio = (
        RECTIFIER_OUTLET_AMMONIA_MASS_FRACTION - weak_fraction
    ) / (strong_fraction - weak_fraction)
    strong_flow = circulation_ratio * refrigerant_flow
    weak_flow = (circulation_ratio - 1.0) * refrigerant_flow

    solution_specific_volume = 1.0 / solution_properties.rho_liquid(
        condenser_temperature_k, strong_fraction
    )
    specific_pump_work = (
        solution_specific_volume * (high_pressure_kpa - low_pressure_kpa) / PUMP_EFFICIENCY
    )
    enthalpy_6 = enthalpy_5 + specific_pump_work
    pump_outlet_temperature = solution_properties.T_from_hL(high_pressure_kpa, strong_fraction, enthalpy_6)
    return {
        "x_strong": strong_fraction, "x_weak": weak_fraction,
        "h5": enthalpy_5, "h8": enthalpy_8, "Tgen": generator_temperature,
        "circulation_ratio": circulation_ratio, "strong_flow": strong_flow, "weak_flow": weak_flow,
        "specific_volume": solution_specific_volume, "specific_pump_work": specific_pump_work,
        "h6": enthalpy_6, "T6": pump_outlet_temperature,
    }


def _solution_heat_exchanger(high_pressure_kpa, solution):
    strong_flow = solution["strong_flow"]
    weak_flow = solution["weak_flow"]
    strong_fraction = solution["x_strong"]
    weak_fraction = solution["x_weak"]
    generator_temperature = solution["Tgen"]
    weak_solution_enthalpy = solution["h8"]
    pumped_solution_enthalpy = solution["h6"]
    pump_outlet_temperature = solution["T6"]

    hot_side_limit = weak_flow * (
        weak_solution_enthalpy
        - solution_properties.h_liquid(high_pressure_kpa, weak_fraction, pump_outlet_temperature)
    ) if weak_flow > 0 else 0.0
    cold_side_limit = strong_flow * (
        solution_properties.h_liquid(high_pressure_kpa, strong_fraction, generator_temperature)
        - pumped_solution_enthalpy
    )
    heat_transfer = SOLUTION_HEAT_EXCHANGER_EFFECTIVENESS * min(hot_side_limit, cold_side_limit)
    enthalpy_7 = pumped_solution_enthalpy + heat_transfer / strong_flow
    enthalpy_9 = weak_solution_enthalpy - heat_transfer / weak_flow if weak_flow > 0 else weak_solution_enthalpy
    temperature_7 = solution_properties.T_from_hL(high_pressure_kpa, strong_fraction, enthalpy_7)
    temperature_9 = solution_properties.T_from_hL(high_pressure_kpa, weak_fraction, enthalpy_9)
    return {"Qshx": heat_transfer, "h7": enthalpy_7, "h9": enthalpy_9, "h10": enthalpy_9, "T7": temperature_7, "T9": temperature_9}


def _generator_and_rectifier(high_pressure_kpa, refrigerant_flow, solution_flow, heat_exchanger):
    strong_flow, weak_flow = solution_flow["strong_flow"], solution_flow["weak_flow"]
    strong_fraction = solution_flow["x_strong"]
    generator_temperature = solution_flow["Tgen"]
    enthalpy_7, enthalpy_8 = heat_exchanger["h7"], solution_flow["h8"]

    generator_liquid_mole_fraction = solution_properties.mass_fraction_to_mole_fraction(
        solution_flow["x_strong"]
    )
    generator_top_temperature = solution_properties.T_bub(
        high_pressure_kpa, generator_liquid_mole_fraction
    )
    vapour_fraction = solution_properties.mole_fraction_to_mass_fraction(
        solution_properties.y_eq(high_pressure_kpa, generator_liquid_mole_fraction)
    )
    feed_bubble_temperature = generator_top_temperature
    rectifier_temperature = solution_properties.T_dew(
        high_pressure_kpa, solution_properties.mass_fraction_to_mole_fraction(RECTIFIER_OUTLET_AMMONIA_MASS_FRACTION)
    )
    reflux_fraction = solution_properties.mole_fraction_to_mass_fraction(
        solution_properties.x_bub(high_pressure_kpa, rectifier_temperature)
    )
    enthalpy_1 = solution_properties.h_vapour(
        high_pressure_kpa, RECTIFIER_OUTLET_AMMONIA_MASS_FRACTION, rectifier_temperature
    )

    if vapour_fraction < RECTIFIER_OUTLET_AMMONIA_MASS_FRACTION:
        fraction_difference = vapour_fraction - reflux_fraction
        if abs(fraction_difference) < 1e-12:
            raise ValueError("Rectifier mass balance is singular: generator vapour and reflux concentrations are too close.")
        reflux_flow = refrigerant_flow * (RECTIFIER_OUTLET_AMMONIA_MASS_FRACTION - vapour_fraction) / fraction_difference
        generator_vapour_enthalpy = solution_properties.h_vapour(high_pressure_kpa, vapour_fraction, generator_top_temperature)
        generator_vapour_flow = refrigerant_flow + reflux_flow
        reflux_enthalpy = solution_properties.h_liquid(high_pressure_kpa, reflux_fraction, rectifier_temperature)
        rectifier_heat = generator_vapour_flow * generator_vapour_enthalpy - refrigerant_flow * enthalpy_1 - reflux_flow * reflux_enthalpy
    else:
        reflux_flow = 0.0
        rectifier_heat = 0.0
        vapour_fraction = RECTIFIER_OUTLET_AMMONIA_MASS_FRACTION
        generator_vapour_flow = refrigerant_flow

    generator_heat = (
        refrigerant_flow * enthalpy_1 + weak_flow * enthalpy_8
        + rectifier_heat - strong_flow * enthalpy_7
    )
    return {
        "h1": enthalpy_1, "T1": rectifier_temperature,
        "Ttop": generator_top_temperature, "Tfeed_bubble": feed_bubble_temperature,
        "yg": vapour_fraction, "generator_liquid_x": solution_flow["x_strong"],
        "generator_liquid_mole_fraction": generator_liquid_mole_fraction,
        "reflux_x": reflux_fraction, "reflux": reflux_flow,
        "Qrect": rectifier_heat, "Qgen": generator_heat, "generator_vapour_flow": generator_vapour_flow,
    }


def solve_cycle(evaporator_cooling_kw, condenser_temperature_k, concentration_difference=DEFAULT_CONCENTRATION_DIFFERENCE, evaporator_temperature_k=DEFAULT_EVAPORATOR_TEMPERATURE_K):
    _validate_inputs(evaporator_cooling_kw, condenser_temperature_k, concentration_difference, evaporator_temperature_k)
    high_pressure_kpa = float(ammonia_properties.psat(condenser_temperature_k))
    low_pressure_kpa = float(ammonia_properties.psat(evaporator_temperature_k))
    refrigerant = _refrigerant_side(
        evaporator_cooling_kw, condenser_temperature_k, evaporator_temperature_k
    )
    solution = _solution_flow_and_pump_side(high_pressure_kpa, low_pressure_kpa, condenser_temperature_k, concentration_difference, refrigerant["refrigerant_flow"])
    heat_exchanger = _solution_heat_exchanger(high_pressure_kpa, solution)
    generator = _generator_and_rectifier(high_pressure_kpa, refrigerant["refrigerant_flow"], solution, heat_exchanger)

    h1, h2, h3, h4 = generator["h1"], refrigerant["h2"], refrigerant["h3"], refrigerant["h4"]
    h5, h6, h7, h8, h9, h10 = solution["h5"], solution["h6"], heat_exchanger["h7"], solution["h8"], heat_exchanger["h9"], heat_exchanger["h10"]
    refrigerant_flow, strong_flow, weak_flow = refrigerant["refrigerant_flow"], solution["strong_flow"], solution["weak_flow"]
    generator_temperature = solution["Tgen"]
    pump_power = strong_flow * solution["specific_pump_work"]
    condenser_heat = refrigerant_flow * (h1 - h2)
    absorber_heat = refrigerant_flow * h4 + weak_flow * h10 - strong_flow * h5
    rectifier_heat = generator["Qrect"]
    generator_heat = generator["Qgen"]

    full_effectiveness_flue_outlet_k = FLUE_GAS_INLET_TEMPERATURE_K - GENERATOR_HEAT_EXCHANGER_EFFECTIVENESS * (FLUE_GAS_INLET_TEMPERATURE_K - generator_temperature)
    flue_outlet_k = full_effectiveness_flue_outlet_k
    flue_floor_active = flue_outlet_k < MINIMUM_FLUE_GAS_OUTLET_TEMPERATURE_K
    if flue_floor_active:
        flue_outlet_k = MINIMUM_FLUE_GAS_OUTLET_TEMPERATURE_K
    flue_mass_flow = generator_heat / (FLUE_GAS_SPECIFIC_HEAT_KJ_KG_K * (FLUE_GAS_INLET_TEMPERATURE_K - flue_outlet_k))

    total_input_power = generator_heat + pump_power
    cop_value = evaporator_cooling_kw / total_input_power if abs(total_input_power) > 1e-12 else np.nan
    reversible_cop_value = (evaporator_temperature_k / (condenser_temperature_k - evaporator_temperature_k)) * ((generator_temperature - condenser_temperature_k) / generator_temperature)
    total_energy_residual = (
        generator_heat + evaporator_cooling_kw + pump_power
        - condenser_heat - absorber_heat - rectifier_heat
    )
    total_mass_residual = strong_flow - (refrigerant_flow + weak_flow)
    ammonia_mass_residual = strong_flow * solution["x_strong"] - (weak_flow * solution["x_weak"] + refrigerant_flow * RECTIFIER_OUTLET_AMMONIA_MASS_FRACTION)
    shx_energy_residual = strong_flow * (h7 - h6) - weak_flow * (h8 - h9)
    pump_volume = solution["specific_volume"]
    condenser_saturation = refrigerant["condenser_saturation"]
    evaporator_saturation = refrigerant["evaporator_saturation"]

    cycle_result = {
        "Qevap": evaporator_cooling_kw, "Tcond": condenser_temperature_k,
        "Phigh": high_pressure_kpa, "Plow": low_pressure_kpa, "Tevap": evaporator_temperature_k,
        "props": PROPERTY_MODEL_NAME,
        "h": {"h1": h1, "h2": h2, "h3": h3, "h4": h4, "h5": h5, "h6": h6, "h7": h7, "h8": h8, "h9": h9, "h10": h10},
        "T": {"T1": generator["T1"], "T2": condenser_temperature_k, "T3": evaporator_temperature_k, "T4": evaporator_temperature_k, "T5": condenser_temperature_k, "T6": solution["T6"], "T7": heat_exchanger["T7"], "T8": generator_temperature, "T9": heat_exchanger["T9"], "T10": solution_properties.T_from_hL(low_pressure_kpa, solution["x_weak"], h10)},
        "mr": refrigerant_flow, "ms": strong_flow, "mw": weak_flow, "f": solution["circulation_ratio"],
        "x_strong": solution["x_strong"], "x_weak": solution["x_weak"], "Tgen": generator_temperature,
        "qe": refrigerant["specific_cooling"], "Wp": pump_power, "wp": solution["specific_pump_work"],
        "pump_m3h": strong_flow * pump_volume * 3600,
        "pump_head_kPa": high_pressure_kpa - low_pressure_kpa,
        "evap_vap_m3h": refrigerant_flow * evaporator_saturation["vg"] * 3600,
        "Qgen": generator_heat, "Qcond": condenser_heat, "Qabs": absorber_heat, "Qrect": rectifier_heat,
        "Qshx": heat_exchanger["Qshx"], "reflux": generator["reflux"], "yg": generator["yg"], "Ttop": generator["Ttop"],
        "Tfg_out": flue_outlet_k, "Tfg_out_eps": full_effectiveness_flue_outlet_k, "mfg": flue_mass_flow,
        "fg_floor_active": flue_floor_active,
        "eps_gen_eff": (FLUE_GAS_INLET_TEMPERATURE_K - flue_outlet_k) / (FLUE_GAS_INLET_TEMPERATURE_K - generator_temperature),
        "COP": cop_value, "COP_rev": reversible_cop_value, "Zg_cond": condenser_saturation["Zg"],
        "bal_E": total_energy_residual, "bal_E_rel": total_energy_residual / max(generator_heat, 1e-9),
        "bal_m": total_mass_residual, "bal_NH3": ammonia_mass_residual, "shx_E": shx_energy_residual,
        "Qreject_total": condenser_heat + absorber_heat + rectifier_heat,
        "chk": {
            "Phigh_ok": high_pressure_kpa <= MAXIMUM_HIGH_SIDE_PRESSURE_KPA,
            "Tcond_ok": condenser_temperature_k >= 298.15 - CONDENSER_TEMPERATURE_TOLERANCE_K,
            "mfg_ok": flue_mass_flow <= MAXIMUM_FLUE_GAS_MASS_FLOW_KG_S,
            "Tfg_ok": flue_outlet_k >= MINIMUM_FLUE_GAS_OUTLET_TEMPERATURE_K - 1e-6,
            "Tgen_lt_Tfgin": generator_temperature < FLUE_GAS_INLET_TEMPERATURE_K,
            "positive_Qgen": np.isfinite(generator_heat) and generator_heat > 0,
            "positive_pump_work": np.isfinite(pump_power) and pump_power >= 0,
            "Tgen_deliverable": (flue_outlet_k >= MINIMUM_FLUE_GAS_OUTLET_TEMPERATURE_K - 1e-6)
                and (flue_outlet_k >= generator_temperature + 5.0)
                and (generator_temperature <= FLUE_GAS_INLET_TEMPERATURE_K - 10.0),
            "T7_below_bubble": heat_exchanger["T7"] <= solution_properties.T_bub(high_pressure_kpa, solution_properties.mass_fraction_to_mole_fraction(solution["x_strong"])),
            "COP_lt_rev": cop_value < reversible_cop_value,
            "COP_band": 0.45 <= cop_value <= 0.65,
            "Tgen_below_Tc_NH3": generator_temperature < ammonia_properties.CRITICAL_TEMPERATURE_K,
        },
    }
    return cycle_result


if __name__ == "__main__":
    for condenser_temperature_k in (301.15, 313.15, 318.15):
        result = solve_cycle(5000.0, condenser_temperature_k)
        print(f"Tcond={condenser_temperature_k - 273.15:.1f} °C | Phigh={result['Phigh']:.0f} kPa | "
              f"Qgen={result['Qgen']:.0f} kW | COP={result['COP']:.3f} | "
              f"m_fg={result['mfg']:.1f} kg/s | energy residual={result['bal_E']:.3e} kW")
