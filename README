# Waste-Heat-Driven Ammonia–Water Absorption Chiller for a District Cooling Network

ME3160 Power & Refrigeration Cycles · End-Semester Term Project 2

Group: _Mahika (ME24BTECH11018), Vaishnavi (ME24BTECH11045)_

This repository models the full chain **outdoor air → AHU psychrometrics → chilled-water load → NH₃–H₂O absorption cycle → condenser → waste-heat (flue-gas) demand**, and sizes the components that must survive a full year (six ambient conditions). COP is always an *output* of the solved cycle, never an input.

---

## 1. Folder layout

```
Project-02/
├── Code/                      all source code (run everything from here)
│   ├── main.py                single reproducible driver
│   ├── user_inputs.py         editable inputs for the reference operating point
│   ├── equilibrium.py         Pátek–Klomfar (1995) bubble/dew/enthalpy correlations (own implementation)
│   ├── solution_props.py      NH₃–H₂O solution properties built on equilibrium.py
│   ├── nh3_props.py           pure NH₃: Peng–Robinson EOS, Lee–Kesler/Pitzer P_sat, Cp₀ (B&S Table A.6)
│   ├── air_handler.py         moist-air relations, mixing, wet-bulb, coil ADP/BF, room-humidity iteration
│   ├── cycle.py               absorption-cycle solver (states 1–10, balances, COP, flue-gas demand)
│   ├── condenser.py           Task 2: air-cooled vs evaporative condenser, design and off-design
│   ├── run_cycle.py           seasonal loads, circuit size, condensate, 40 % fresh air (Tasks 1, 3, 4)
│   ├── cycle_study.py         state-point tables, component balances, state-point plots
│   ├── study.py               condenser-temperature and degassing-width studies (--sensitivity for more)
│   ├── check_balances.py      engineering checks and the Validation Table (CSV)
│   ├── test_validation.py     automatic assertions for Checks 1–3 (see section 4)
│   └── requirements.txt
├── tests/                     unit tests (equilibrium, cycle, cycle study)
├── results/                   generated tables, plots, results.json
└── Diagrams/                  static schematics: absorption cycle and one-zone AHU
```

`main.py` writes its output to a `results/` folder next to the scripts. In the submitted archive this folder has been placed beside `Code/`.

## 2. Setup

Python 3.11 or newer.

```bash
cd Code
python -m pip install -r requirements.txt
```

Packages: `numpy`, `scipy`, `pandas`, `matplotlib`, `pytest`, and `CoolProp`.
CoolProp is used **only for pure water** (saturation pressure and h_g inside the psychrometric relations). If it is not installed, `air_handler.py` falls back to a built-in Antoine/Bridgeman correlation; the validation checks still pass, but install CoolProp to reproduce the reported numbers exactly.

## 3. How to run

| Command (from `Code/`) | What it does |
|---|---|
| `python main.py` | Runs everything: cycle, seasonal analysis, condenser design, studies, plots, checks. Regenerates `results/`. |
| `python -m pytest -q test_validation.py ../tests` | Runs the validation test and all unit tests. |
| `python check_balances.py` | Prints and saves the Validation Table (Quantity / Hand value / Code output / % difference / Pass-Fail). |
| `python cycle_study.py` | Recomputes only the reference state-point tables and plots. |
| `python study.py` | Recomputes the condenser-temperature and degassing-width studies. |
| `python study.py --sensitivity` | One-at-a-time sensitivity to solution density, flue-gas c_p, ε_gen and condenser UA/flow exponent. |

To change the reference operating point, edit the four values at the top of `user_inputs.py` and run `python main.py` again:

```python
COOLING_DUTY_KW = 5000.0
EVAPORATOR_TEMPERATURE_C = 2.0
CONDENSER_TEMPERATURE_C = 40.0
DEGASSING_WIDTH = 0.16
```

## 4. Validation

The three required checks (Section 2.4 of the project statement) are asserted automatically in `test_validation.py`:

| Check | What is tested | Tolerance |
|---|---|---|
| 1. Psychrometrics | P_v, w, h of room and outdoor states; 70/30 mixing; T_wb of summer day (≈ 30.3 °C) and monsoon day (≈ 32.6 °C) | w, h: 0.5 %; T: 0.1 K; T_wb: 0.15 K |
| 2. Ammonia EOS vs Table B.2.1 | P_sat and h_fg at 0 °C and 40 °C; q_e = h_g(0 °C) − h_f(40 °C); ṁ_r at 5 MW | P_sat ±5 %, h_fg ±8 %, q_e 1050–1100 kJ/kg, ṁ_r 4.5–4.8 kg/s |
| 3. Cycle balances | mass, NH₃ and energy balances; COP < COP_rev; COP_rev(2/45/130 °C) ≈ 1.35; COP in 0.45–0.65; COP falls as T_cond rises | residuals < 1e-6; band as stated |

`check_balances.py` produces the full Validation Table, saved as `results/model_check_summary.csv`. Pátek–Klomfar is a separate tested module (`equilibrium.py`, tested in `tests/test_equilibrium.py`).

Note that the energy-balance residual in Check 3 closes by construction; it confirms bookkeeping, not property accuracy. The real evidence of property quality is Checks 1 and 2, the COP < COP_rev test and the COP band.

## 5. Models and assumptions

**Reference state.** T_ref = 273.15 K, P_ref = 101.325 kPa, ideal-gas state, h_ref = s_ref = 0, carried through all routines. The P&K pure-NH₃ saturated-vapour enthalpy is anchored to the Peng–Robinson value at 273.15 K so that the solution and refrigerant enthalpies share one basis.

**Ammonia–water solution.** Pátek & Klomfar (1995), implemented from scratch: bubble temperature (Eq. 6), dew temperature (Eq. 7), vapour composition and liquid/vapour enthalpies (Eqs. 8–10). Three coefficients of the handout tables differ from the physics and are corrected in `equilibrium.py` (each verified against known saturation data and covered by `tests/test_equilibrium.py`):

| Equation | Coefficient | Handout | Used | Test |
|---|---|---|---|---|
| Eq. 6 | a₁₃ (m=6, n=0) | 4.89596E−02 | 4.89596 | water boils at 2 MPa: 486 K vs 485.5 K true |
| Eq. 7 | a₃ (m=0, n=2) | −4.35624E−02 | +4.35624E−02 | NH₃ dew point at 0.1 MPa: 238.6 K vs 239.7 K true |
| Eq. 7 | a₁₀ (m=3, n=0) | −2.1078E+01 | −2.0178E+01 | H₂O dew point at 0.1 MPa: 373.0 K vs 372.8 K true |

**Pure ammonia.** Peng–Robinson EOS with T_c = 405.56 K, P_c = 11.333 MPa, ω = 0.256; ideal-gas part from the C_p⁰ polynomial (B&S Table A.6) plus the EOS departure function; P_sat from the Lee–Kesler/Pitzer correlation. A real-gas EOS is needed because Z_g of saturated vapour is about 0.95 at 0 °C and 0.88 at 40 °C. The EOS sets P_high and P_low and is what Check 2 validates. The refrigerant is treated as **pseudo-pure NH₃** (rectifier outlet is 99.5 % NH₃): enthalpies of states 2–4 use the P&K pure-NH₃ limit.

**Moist air.** Ideal-gas mixture (Z ≈ 0.999). P_sat,H₂O from CoolProp; w = 0.622 P_v/(P − P_v); h = 1.006 T + w(2501 + 1.86 T). Wet-bulb temperature from the adiabatic-saturation energy balance, solved iteratively. The room humidity ratio is found by iterating the moisture balance to |Δw| < 1e-10.

**Generator and rectifier closure.** The vapour leaving the top of the generator/rectifier column is taken in equilibrium with the strong feed (state 5); the weak solution (state 8) sets the generator outlet temperature and enthalpy. The weak-feed alternative is reported only as a sensitivity. The degassing width Δx = x_strong − x_weak = 0.16 is an assumed design parameter (it is not given in the project statement); varying it over 0.12–0.20 changes COP by less than 0.015.

**Other assumptions.** T_abs = T_cond; no pressure drops in vessels and pipes; flue gas c_p = 1.10 kJ/(kg·K), 300 °C inlet, ≥ 140 °C outlet; ε_gen = 0.75, ε_SHX = 0.70, η_pump = 0.70; solution density from a reciprocal mass-weighted estimate; condenser UA and sink capacity scale with an empirical flow exponent (0.8). These are screening assumptions, not vendor-certified values.

## 6. Where to find each result

| Project task | Main output files |
|---|---|
| Reference cycle (5 MW, 2 °C, 40 °C) | `results/study_operating_point.csv`, `study_state_points.csv`, `study_component_balances.csv` |
| Task 1: loads, Q_gen, circuit size | `results/details/air_handler_calculation_results.csv`, `seasonal_operating_point_results.csv`, `coupled_cycle_results.csv` |
| Task 1: psychrometric chart | `results/plots/air_moisture_process_chart.png` |
| Task 2: condenser comparison and off-design | `results/details/condenser_option_comparison.csv`, `resized_condenser_feasibility.csv`, `recommended_condenser_seasonal_results.csv` |
| Task 3: condensate | `results/details/condensate_recovery_estimates.csv` |
| Task 4: 40 % fresh air | `results/details/forty_percent_fresh_air_case.csv`, `monsoon_night_humidity_adjustment.csv`, `exhaust_gas_limit_solutions.csv` |
| Sensitivity | `results/details/assumption_sensitivity.csv`, `solution_concentration_sensitivity.csv` |
| Validation table | `results/model_check_summary.csv` |
| Schematics | `Diagrams/overview.png`, `Diagrams/AHU Airflow and Condensate Diagram.png` |

`results/details/` holds supporting CSVs and `results/plots/` holds all figures.

## 7. Headline results

| Quantity | Value |
|---|---|
| Reference COP (5 MW, T_evap 2 °C, T_cond 40 °C) | 0.522 |
| Reference ṁ_r / Q_gen / pump power | 4.65 kg/s / 9.53 MW / 42 kW |
| Governing condition for chiller, condenser and decanter | Monsoon day (Q_evap 4.91 MW, T_wb 32.6 °C, condensate 3,739 kg/h) |
| Flue-gas demand over the year (30 % fresh air) | 0 (winter night, chiller off) to 63.5 kg/s (limit 80) |
| Recommended condenser | Evaporative: UA 929 kW/K, 1,548 m², design T_cond 40.6 °C, P_high 15.9 bar |
| Air-cooled condenser | Rejected: P_high 21.5 bar exceeds the 20 bar vessel limit |
| 40 % fresh air | Coil load +10 to +23 %; monsoon day needs 6.0 MW and a 1.32× larger condenser (76 kg/s flue gas); only the 80 kg/s limit becomes active, the 20 bar limit does not |
| Biggest design lever | Heat-rejection temperature T_cond (COP −0.011 per K, Q_gen +200 kW per K) |

## 8. Limitations

- Results are model predictions from the correlations above; passing the checks shows consistency with the stated hand values, not agreement with real equipment.
- The ammonia EOS and the P&K enthalpy correlations are separate models. Their latent-heat values differ by a few percent at the checked points; this is a model-form uncertainty that internal balance closure cannot remove.
- The cubic EOS reproduces P_sat to about 1 % and h_fg to about 3–4 % of the table values, which is the expected accuracy for a polar, hydrogen-bonding fluid.
- Δx, solution density, condenser UA scaling and flue-gas properties are assumptions and would need plant or vendor data for a final design.
- Winter night has a negative raw coil balance; it is booked as a zero-load, chiller-off state with free cooling, and no percentage changes are reported for it.
