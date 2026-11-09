"""Run the state-point/component-balance solution method.

Outputs (in results/):
- study_operating_point.csv: compact solved cycle summary
- study_state_points.csv: state-wise T, P, h, concentration and mass flow
- study_component_balances.csv: component and cycle residual checks
- study_temperature_sweep.csv: performance at user-configured condenser temperatures
- study_state_point_chart.png and study_method_performance.png
"""
from __future__ import annotations
import csv
import math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
import numpy as np
import cycle
import nh3_props
import solution_props
from user_inputs import (
    COOLING_DUTY_KW,
    EVAPORATOR_TEMPERATURE_C,
    CONDENSER_TEMPERATURE_C,
    DEGASSING_WIDTH,
    CONDENSER_SWEEP_OFFSETS_C,
)

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)
KELVIN_OFFSET = 273.15


def _write_csv(file_path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows to write to {file_path.name}.")
    with file_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _validate_study_inputs() -> list[float]:
    scalar_values = {
        "COOLING_DUTY_KW": COOLING_DUTY_KW,
        "EVAPORATOR_TEMPERATURE_C": EVAPORATOR_TEMPERATURE_C,
        "CONDENSER_TEMPERATURE_C": CONDENSER_TEMPERATURE_C,
        "DEGASSING_WIDTH": DEGASSING_WIDTH,
    }
    for name, value in scalar_values.items():
        if not math.isfinite(float(value)):
            raise ValueError(f"{name} must be a finite number.")
    if COOLING_DUTY_KW <= 0:
        raise ValueError("COOLING_DUTY_KW must be greater than zero.")
    if CONDENSER_TEMPERATURE_C <= EVAPORATOR_TEMPERATURE_C:
        raise ValueError("CONDENSER_TEMPERATURE_C must exceed EVAPORATOR_TEMPERATURE_C.")
    if not 0 < DEGASSING_WIDTH < 1:
        raise ValueError("DEGASSING_WIDTH must be between 0 and 1.")
    max_temp_c = nh3_props.CRITICAL_TEMPERATURE_K - KELVIN_OFFSET
    sweep_temperatures_c = sorted(set(round(CONDENSER_TEMPERATURE_C + offset, 6) for offset in CONDENSER_SWEEP_OFFSETS_C) | {float(CONDENSER_TEMPERATURE_C)})
    for temperature_c in sweep_temperatures_c:
        if not EVAPORATOR_TEMPERATURE_C < temperature_c < max_temp_c:
            raise ValueError(
                f"Every swept condenser temperature must be above the evaporator temperature "
                f"and below ammonia's critical temperature ({max_temp_c:.2f} °C); "
                f"received {temperature_c:.2f} °C. Edit user_inputs.py."
            )
                                                                                           
    return sweep_temperatures_c


def _state_point_rows(cycle_result: dict) -> list[dict]:
    high_pressure = cycle_result["Phigh"]
    low_pressure = cycle_result["Plow"]
    enthalpies = cycle_result["h"]
    temperatures = cycle_result["T"]
    strong_fraction = cycle_result["x_strong"]
    weak_fraction = cycle_result["x_weak"]
    refrigerant_flow = cycle_result["mr"]
    strong_flow = cycle_result["ms"]
    weak_flow = cycle_result["mw"]
    states = [
        (1, "Rectifier outlet vapour", high_pressure, temperatures["T1"], enthalpies["h1"], 0.995, refrigerant_flow, "Refrigerant"),
        (2, "Condenser outlet liquid", high_pressure, temperatures["T2"], enthalpies["h2"], 1.000, refrigerant_flow, "Refrigerant"),
        (3, "Expansion valve outlet", low_pressure, temperatures["T3"], enthalpies["h3"], 1.000, refrigerant_flow, "Refrigerant"),
        (4, "Evaporator outlet vapour", low_pressure, temperatures["T4"], enthalpies["h4"], 1.000, refrigerant_flow, "Refrigerant"),
        (5, "Absorber outlet strong solution", low_pressure, temperatures["T5"], enthalpies["h5"], strong_fraction, strong_flow, "Strong solution"),
        (6, "Pump outlet", high_pressure, temperatures["T6"], enthalpies["h6"], strong_fraction, strong_flow, "Strong solution"),
        (7, "Solution HX strong-side outlet", high_pressure, temperatures["T7"], enthalpies["h7"], strong_fraction, strong_flow, "Strong solution"),
        (8, "Generator outlet weak solution", high_pressure, temperatures["T8"], enthalpies["h8"], weak_fraction, weak_flow, "Weak solution"),
        (9, "Solution HX weak-side outlet", high_pressure, temperatures["T9"], enthalpies["h9"], weak_fraction, weak_flow, "Weak solution"),
        (10, "Weak-solution throttle outlet", low_pressure, temperatures["T10"], enthalpies["h10"], weak_fraction, weak_flow, "Weak solution"),
    ]
    return [
        {
            "state": state_number,
            "state_description": description,
            "stream": stream_name,
            "pressure_kPa": pressure,
            "temperature_C": temperature - KELVIN_OFFSET,
            "enthalpy_kJ_kg": enthalpy,
            "ammonia_mass_fraction": ammonia_fraction,
            "mass_flow_kg_s": mass_flow,
        }
        for state_number, description, pressure, temperature, enthalpy, ammonia_fraction, mass_flow, stream_name in states
    ]


def _component_balance_rows(cycle_result: dict) -> list[dict]:
    h = cycle_result["h"]
    mass_r = cycle_result["mr"]
    mass_strong = cycle_result["ms"]
    mass_weak = cycle_result["mw"]
    energy_checks = [
        ("Evaporator duty: Qe - mr*(h4-h3)", cycle_result["Qevap"] - mass_r * (h["h4"] - h["h3"]), "kW"),
        ("Condenser duty: Qcond - mr*(h1-h2)", cycle_result["Qcond"] - mass_r * (h["h1"] - h["h2"]), "kW"),
        ("Generator duty balance", cycle_result["Qgen"] - (mass_r * h["h1"] + mass_weak * h["h8"] + cycle_result["Qrect"] - mass_strong * h["h7"]), "kW"),
        ("Absorber duty balance", cycle_result["Qabs"] - (mass_r * h["h4"] + mass_weak * h["h10"] - mass_strong * h["h5"]), "kW"),
        ("Overall cycle energy balance", cycle_result["bal_E"], "kW"),
        ("Solution heat-exchanger energy balance", cycle_result["shx_E"], "kW"),
        ("Total mass balance", cycle_result["bal_m"], "kg/s"),
        ("Ammonia mass balance", cycle_result["bal_NH3"], "kg/s"),
    ]
    return [
        {"check": name, "residual": float(residual), "units": units, "absolute_tolerance": 1e-6,
         "result": "PASS" if math.isfinite(float(residual)) and abs(float(residual)) <= 1e-6 else "FAIL"}
        for name, residual, units in energy_checks
    ]


def main() -> None:
    sweep_temperatures_c = _validate_study_inputs()
    condenser_temperature_k = CONDENSER_TEMPERATURE_C + KELVIN_OFFSET
    evaporator_temperature_k = EVAPORATOR_TEMPERATURE_C + KELVIN_OFFSET
    solved_cycle = cycle.solve_cycle(
        COOLING_DUTY_KW, condenser_temperature_k,
        concentration_difference=DEGASSING_WIDTH,
        evaporator_temperature_k=evaporator_temperature_k,
    )

    summary_row = {
        "cooling_duty_kW": COOLING_DUTY_KW,
        "evaporator_temperature_C": EVAPORATOR_TEMPERATURE_C,
        "condenser_temperature_C": CONDENSER_TEMPERATURE_C,
        "degassing_width": DEGASSING_WIDTH,
        "high_pressure_kPa": solved_cycle["Phigh"],
        "low_pressure_kPa": solved_cycle["Plow"],
        "refrigerant_flow_kg_s": solved_cycle["mr"],
        "strong_solution_flow_kg_s": solved_cycle["ms"],
        "weak_solution_flow_kg_s": solved_cycle["mw"],
        "generator_heat_kW": solved_cycle["Qgen"],
        "pump_power_kW": solved_cycle["Wp"],
        "COP": solved_cycle["COP"],
        "condenser_heat_kW": solved_cycle["Qcond"],
        "absorber_heat_kW": solved_cycle["Qabs"],
        "rectifier_heat_kW": solved_cycle["Qrect"],
        "generator_temperature_C": solved_cycle["Tgen"] - KELVIN_OFFSET,
        "exhaust_gas_flow_kg_s": solved_cycle["mfg"],
        "energy_balance_residual_kW": solved_cycle["bal_E"],
        "ammonia_balance_residual_kg_s": solved_cycle["bal_NH3"],
    }
    _write_csv(RESULTS / "study_operating_point.csv", [summary_row])
    state_rows = _state_point_rows(solved_cycle)
    _write_csv(RESULTS / "study_state_points.csv", state_rows)
    balance_rows = _component_balance_rows(solved_cycle)
    _write_csv(RESULTS / "study_component_balances.csv", balance_rows)
    if not all(row["result"] == "PASS" for row in balance_rows):
        failed = [row for row in balance_rows if row["result"] != "PASS"]
        raise RuntimeError(f"Study balance check failed: {failed}")

    temperature_rows = []
    for temperature_c in sweep_temperatures_c:
        result = cycle.solve_cycle(
            COOLING_DUTY_KW, temperature_c + KELVIN_OFFSET,
            concentration_difference=DEGASSING_WIDTH,
            evaporator_temperature_k=evaporator_temperature_k,
        )
        temperature_rows.append({
            "cooling_duty_kW": COOLING_DUTY_KW,
            "evaporator_temperature_C": EVAPORATOR_TEMPERATURE_C,
            "condenser_temperature_C": temperature_c,
            "degassing_width": DEGASSING_WIDTH,
            "COP": result["COP"],
            "generator_heat_kW": result["Qgen"],
            "pump_power_kW": result["Wp"],
            "high_pressure_kPa": result["Phigh"],
            "generator_temperature_C": result["Tgen"] - KELVIN_OFFSET,
            "exhaust_gas_flow_kg_s": result["mfg"],
            "energy_balance_residual_kW": result["bal_E"],
            "checks_passed": all(result["chk"].values()),
        })
    _write_csv(RESULTS / "study_temperature_sweep.csv", temperature_rows)

                                                                                               
    figure, axes = plt.subplots(1, 3, figsize=(14.5, 4.8))
    state_numbers = [row["state"] for row in state_rows]
    state_temperatures = [row["temperature_C"] for row in state_rows]
    state_enthalpies = [row["enthalpy_kJ_kg"] for row in state_rows]
    axes[0].plot(state_numbers, state_temperatures, marker="o", linestyle="None", color="#2457A7")
    axes[0].set_xlabel("State number")
    axes[0].set_ylabel("Temperature (°C)")
    axes[0].set_title("Solved state temperatures")
    axes[0].set_xticks(state_numbers)
    axes[0].grid(True, alpha=0.3)
    axes[1].plot(state_numbers, state_enthalpies, marker="s", linestyle="None", color="#0F8B8D")
    axes[1].axhline(0.0, color="0.45", linestyle="--", linewidth=0.9)
    axes[1].text(0.02, 0.03, "Negative values depend on the selected enthalpy reference;\ncycle duties use enthalpy differences.", transform=axes[1].transAxes, fontsize=8, va="bottom")
    axes[1].set_xlabel("State number")
    axes[1].set_ylabel("Specific enthalpy, project reference (kJ/kg)")
    axes[1].set_title("Solved state enthalpies")
    axes[1].set_xticks(state_numbers)
    axes[1].grid(True, alpha=0.3)
    mass_flows = [row["mass_flow_kg_s"] for row in state_rows]
    axes[2].plot(state_numbers, mass_flows, marker="^", linestyle="None", color="#D97706")
    axes[2].set_xlabel("State number")
    axes[2].set_ylabel("Mass flow (kg/s)")
    axes[2].set_title("Stream mass flow rates")
    axes[2].set_xticks(state_numbers)
    axes[2].grid(True, alpha=0.3)
    figure.suptitle(
        f"State-point method | Qe={COOLING_DUTY_KW:g} kW, "
        f"Te={EVAPORATOR_TEMPERATURE_C:g} °C, Tc={CONDENSER_TEMPERATURE_C:g} °C, "
        f"dx={DEGASSING_WIDTH:g}"
    )
    figure.tight_layout()
    figure.savefig(RESULTS / "study_state_point_chart.png", dpi=180, bbox_inches="tight")
    plt.close(figure)

                                                                                          
    figure, cop_axis = plt.subplots(figsize=(8.3, 5.0))
    temperatures = [row["condenser_temperature_C"] for row in temperature_rows]
    cop_values = [row["COP"] for row in temperature_rows]
    generator_heat_values = [row["generator_heat_kW"] for row in temperature_rows]
    cop_axis.plot(temperatures, cop_values, marker="o", color="#2457A7", label="Thermal COP")
    cop_axis.set_xlabel("Condenser temperature (°C)")
    cop_axis.set_ylabel("Thermal COP (dimensionless)")
    cop_axis.grid(True, alpha=0.3)
    heat_axis = cop_axis.twinx()
    heat_axis.plot(temperatures, generator_heat_values, marker="s", linestyle="--", color="#0F8B8D", label="Generator heat input")
    heat_axis.set_ylabel("Generator heat input (kW)")
    cop_axis.axvline(CONDENSER_TEMPERATURE_C, color="#7656A8", linestyle=":", linewidth=1.5, label="Selected base temperature")
    lines = cop_axis.get_lines() + heat_axis.get_lines()
    cop_axis.legend(lines, [line.get_label() for line in lines], loc="best")
    cop_axis.set_title(
        f"Study performance | Qe={COOLING_DUTY_KW:g} kW, "
        f"Te={EVAPORATOR_TEMPERATURE_C:g} °C, dx={DEGASSING_WIDTH:g}"
    )
    figure.tight_layout()
    figure.savefig(RESULTS / "study_method_performance.png", dpi=180, bbox_inches="tight")
    plt.close(figure)

    print("Project state-point/component-balance method complete.")
    print(f"Base point: COP={solved_cycle['COP']:.5f}; Qgen={solved_cycle['Qgen']:.3f} kW; "
          f"P_high={solved_cycle['Phigh'] / 100:.3f} bar")
    print(f"Component and cycle balance checks: {sum(row['result'] == 'PASS' for row in balance_rows)}/{len(balance_rows)} passed.")
    print("Study CSVs and plots are saved in results/.")


if __name__ == "__main__":
    main()
