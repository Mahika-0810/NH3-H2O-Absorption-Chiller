"""Parametric studies driven by the editable values in user_inputs.py.

The selected condenser temperature anchors the temperature sweep, and the
selected concentration difference anchors the operating point. Outputs are
written to results/study_*.csv and results/study_*.png.
"""
from pathlib import Path
import csv
import math
import sys
import pandas as pd
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
import cycle
import nh3_props
import solution_props
import condenser
from user_inputs import (
    COOLING_DUTY_KW,
    EVAPORATOR_TEMPERATURE_C,
    CONDENSER_TEMPERATURE_C,
    DEGASSING_WIDTH,
    CONDENSER_SWEEP_OFFSETS_C,
    DEGASSING_WIDTHS_TO_COMPARE,
)

RESULTS_FOLDER = Path(__file__).resolve().parent / "results"
RESULTS_FOLDER.mkdir(exist_ok=True)
KELVIN_OFFSET = 273.15


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"Cannot write an empty study: {path.name}")
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    temps = sorted(set(
        round(CONDENSER_TEMPERATURE_C + offset, 6) for offset in CONDENSER_SWEEP_OFFSETS_C
    ) | {float(CONDENSER_TEMPERATURE_C)})
    for temp in temps:
        if temp <= EVAPORATOR_TEMPERATURE_C or temp + KELVIN_OFFSET >= nh3_props.CRITICAL_TEMPERATURE_K:
            raise ValueError("Temperature sweep contains a non-physical condenser temperature. Edit user_inputs.py.")

    temp_rows = []
    for t_cond in temps:
        cycle_data = cycle.solve_cycle(
            COOLING_DUTY_KW,
            t_cond + KELVIN_OFFSET,
            concentration_difference=DEGASSING_WIDTH,
            evaporator_temperature_k=EVAPORATOR_TEMPERATURE_C + KELVIN_OFFSET,
        )
        temp_rows.append({
            "Qevap_kW": COOLING_DUTY_KW,
            "Tevap_C": EVAPORATOR_TEMPERATURE_C,
            "Tcond_C": t_cond,
            "degassing_width_dx": DEGASSING_WIDTH,
            "COP": round(float(cycle_data["COP"]), 6),
            "Qgen_kW": round(float(cycle_data["Qgen"]), 3),
            "Tgen_C": round(float(cycle_data["Tgen"] - KELVIN_OFFSET), 3),
            "Phigh_bar": round(float(cycle_data["Phigh"] / 100.0), 4),
            "flue_gas_kg_s": round(float(cycle_data["mfg"]), 4),
            "pressure_below_20bar": bool(cycle_data["chk"]["Phigh_ok"]),
            "flue_gas_below_80kg_s": bool(cycle_data["chk"]["mfg_ok"]),
            "energy_balance_residual_kW": float(cycle_data["bal_E"]),
        })
    write_csv(RESULTS_FOLDER / "study_condenser_temperature_sweep.csv", temp_rows)

    conc_rows = []
    for dx in sorted(set(DEGASSING_WIDTHS_TO_COMPARE) | {float(DEGASSING_WIDTH)}):
        cycle_data = cycle.solve_cycle(
            COOLING_DUTY_KW,
            CONDENSER_TEMPERATURE_C + KELVIN_OFFSET,
            concentration_difference=dx,
            evaporator_temperature_k=EVAPORATOR_TEMPERATURE_C + KELVIN_OFFSET,
        )
        conc_rows.append({
            "Qevap_kW": COOLING_DUTY_KW,
            "Tevap_C": EVAPORATOR_TEMPERATURE_C,
            "Tcond_C": CONDENSER_TEMPERATURE_C,
            "degassing_width_dx": dx,
            "COP": round(float(cycle_data["COP"]), 6),
            "Qgen_kW": round(float(cycle_data["Qgen"]), 3),
            "Tgen_C": round(float(cycle_data["Tgen"] - KELVIN_OFFSET), 3),
            "Phigh_bar": round(float(cycle_data["Phigh"] / 100.0), 4),
            "flue_gas_kg_s": round(float(cycle_data["mfg"]), 4),
            "energy_balance_residual_kW": float(cycle_data["bal_E"]),
        })
    write_csv(RESULTS_FOLDER / "study_solution_concentration_sweep.csv", conc_rows)

                                                                                       
    figure, ax = plt.subplots(figsize=(8.2, 5.0))
    temps = [row["Tcond_C"] for row in temp_rows]
    ax.plot(temps, [row["COP"] for row in temp_rows], marker="o", color="#2457A7", label="Thermal COP")
    ax.axvline(CONDENSER_TEMPERATURE_C, color="#7656A8", linestyle=":", linewidth=1.5, label="Selected base temperature")
    ax.set_xlabel("Condenser temperature (°C)")
    ax.set_ylabel("Thermal COP (dimensionless)")
    ax.grid(True, alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(temps, [row["Qgen_kW"] for row in temp_rows], marker="s", color="#0F8B8D", linestyle="--", label="Generator heat input")
    ax2.set_ylabel("Generator heat input (kW)")
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [line.get_label() for line in lines], loc="best", fontsize=8)
    high_flow_rows = [row for row in temp_rows if row["flue_gas_kg_s"] > cycle.MAXIMUM_FLUE_GAS_MASS_FLOW_KG_S]
    if high_flow_rows:
        overloaded = max(high_flow_rows, key=lambda row: row["Tcond_C"])
        ax.annotate(
            f"Exhaust-gas limit exceeded\n({overloaded['flue_gas_kg_s']:.1f} kg/s)",
            xy=(overloaded["Tcond_C"], overloaded["COP"]),
            xytext=(0.54, 0.16), textcoords="axes fraction",
            arrowprops={"arrowstyle": "->", "linewidth": 0.8}, fontsize=8,
        )
    ax.set_title(
        f"Condenser-temperature sensitivity | Qe={COOLING_DUTY_KW:g} kW, "
        f"Te={EVAPORATOR_TEMPERATURE_C:g} °C, dx={DEGASSING_WIDTH:g}"
    )
    figure.tight_layout()
    figure.savefig(RESULTS_FOLDER / "study_condenser_temperature_plot.png", dpi=180, bbox_inches="tight")
    plt.close(figure)

                                                                                           
    figure, axis = plt.subplots(figsize=(7.5, 4.8))
    labels = [f"dx = {row['degassing_width_dx']:.3f}" for row in conc_rows]
    axis.bar(labels, [row["COP"] for row in conc_rows], color=["#2457A7", "#0F8B8D", "#D97706", "#C2414B", "#7656A8"][:len(labels)])
    axis.set_xlabel("Assumed degassing width (ammonia mass fraction)")
    axis.set_ylabel("Thermal COP (dimensionless)")
    axis.set_title(
        f"Degassing-width sensitivity | Qe={COOLING_DUTY_KW:g} kW, "
        f"Tc={CONDENSER_TEMPERATURE_C:g} °C"
    )
    axis.grid(axis="y", alpha=0.3)
    figure.tight_layout()
    figure.savefig(RESULTS_FOLDER / "study_solution_concentration_plot.png", dpi=180, bbox_inches="tight")
    plt.close(figure)

    print("Configured parametric study complete.")
    print("Condenser-temperature sweep:")
    for row in temp_rows:
        print(f"  Tc={row['Tcond_C']:>6.2f} °C | COP={row['COP']:.4f} | "
              f"Qgen={row['Qgen_kW']:.1f} kW | Tgen={row['Tgen_C']:.1f} °C | "
              f"Phigh={row['Phigh_bar']:.2f} bar")
    print("\nDegassing-width sensitivity:")
    for row in conc_rows:
        print(f"  dx={row['degassing_width_dx']:.3f} | COP={row['COP']:.4f} | "
              f"Qgen={row['Qgen_kW']:.1f} kW | Tgen={row['Tgen_C']:.1f} °C")
    print("CSV and PNG files saved in results/.")



def run_sensitivity():
    """Check how selected model assumptions affect the main outputs."""
    out_file = RESULTS_FOLDER / "assumption_sensitivity.csv"
    cooling_kw = 5000.0
    condenser_k = 313.15
    evaporator_k = 275.15
    dx = 0.16

    def solve():
        return cycle.solve_cycle(cooling_kw, condenser_k, dx, evaporator_k)

    rows = []
    base = solve()
    rows.append({
        "assumption": "baseline",
        "perturbation": "current stated assumptions",
        "COP": base["COP"],
        "Qgen_kW": base["Qgen"],
        "Wp_kW": base["Wp"],
        "mfg_kg_s": base["mfg"],
        "Qcond_kW": base["Qcond"],
    })

    original_density = solution_props.rho_liquid
    try:
        for scale in (0.90, 1.10):
            solution_props.rho_liquid = (
                lambda temperature_k, mass_fraction, factor=scale:
                original_density(temperature_k, mass_fraction) * factor
            )
            result = solve()
            rows.append({
                "assumption": "solution density",
                "perturbation": f"{(scale - 1) * 100:+.0f}% density",
                "COP": result["COP"],
                "Qgen_kW": result["Qgen"],
                "Wp_kW": result["Wp"],
                "mfg_kg_s": result["mfg"],
                "Qcond_kW": result["Qcond"],
            })
    finally:
        solution_props.rho_liquid = original_density

    original_cp = cycle.FLUE_GAS_SPECIFIC_HEAT_KJ_KG_K
    try:
        for scale in (0.90, 1.10):
            cycle.FLUE_GAS_SPECIFIC_HEAT_KJ_KG_K = original_cp * scale
            result = solve()
            rows.append({
                "assumption": "flue-gas specific heat",
                "perturbation": f"{(scale - 1) * 100:+.0f}% cp",
                "COP": result["COP"],
                "Qgen_kW": result["Qgen"],
                "Wp_kW": result["Wp"],
                "mfg_kg_s": result["mfg"],
                "Qcond_kW": result["Qcond"],
            })
    finally:
        cycle.FLUE_GAS_SPECIFIC_HEAT_KJ_KG_K = original_cp

    original_efficiency = cycle.GENERATOR_HEAT_EXCHANGER_EFFECTIVENESS
    try:
        for effectiveness in (0.65, 0.75, 0.85):
            cycle.GENERATOR_HEAT_EXCHANGER_EFFECTIVENESS = effectiveness
            result = solve()
            rows.append({
                "assumption": "generator/flue exchanger effectiveness",
                "perturbation": f"effectiveness={effectiveness:.2f}",
                "COP": result["COP"],
                "Qgen_kW": result["Qgen"],
                "Wp_kW": result["Wp"],
                "mfg_kg_s": result["mfg"],
                "Qcond_kW": result["Qcond"],
            })
    finally:
        cycle.GENERATOR_HEAT_EXCHANGER_EFFECTIVENESS = original_efficiency

    original_exponent = condenser.SINK_SIDE_UA_FLOW_EXPONENT
    try:
        for exponent in (0.6, 0.8, 1.0):
            condenser.SINK_SIDE_UA_FLOW_EXPONENT = exponent
            rows.append({
                "assumption": "condenser UA-flow exponent",
                "perturbation": f"exponent={exponent:.1f}",
                "COP": float("nan"),
                "Qgen_kW": float("nan"),
                "Wp_kW": float("nan"),
                "mfg_kg_s": float("nan"),
                "Qcond_kW": float("nan"),
                "UA_fraction_at_50pct_flow": 0.5 ** exponent,
            })
    finally:
        condenser.SINK_SIDE_UA_FLOW_EXPONENT = original_exponent

    pd.DataFrame(rows).to_csv(out_file, index=False)
    print(f"Wrote assumption sensitivity: {out_file.name} ({len(rows)} cases)")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--sensitivity":
        run_sensitivity()
    else:
        main()
