# Study Guide: State-Point Method for an Ammonia–Water Absorption Chiller

## Aim

This project evaluates a single-effect ammonia–water absorption chiller driven by waste heat. The project solver makes the state-point definitions and component mass/energy balances explicit. It calculates the thermal COP from the solved cycle rather than treating COP as a given input.

## Default operating point

Unless `user_inputs.py` is changed, the project calculation uses:

| Input | Default value |
|---|---:|
| Evaporator cooling duty | 5,000 kW |
| Evaporator temperature | 2 °C |
| Condenser/absorber temperature | 40 °C |
| Solution degassing width, `dx` | 0.16 mass fraction |
| Rectifier outlet ammonia mass fraction | 0.995 |

At that model operating point, the outputs are approximately:

| Calculated quantity | Result |
|---|---:|
| High-side pressure | 15.64 bar |
| Low-side pressure | 4.58 bar |
| Refrigerant mass flow | 4.649 kg/s |
| Strong solution flow | 19.762 kg/s |
| Weak solution flow | 15.113 kg/s |
| Generator heat input | 9,531.95 kW |
| Solution pump power | 42.08 kW |
| Thermal COP | 0.52225 |
| Generator solution temperature | 117.13 °C |
| Required exhaust-gas flow | 63.18 kg/s |

These values are model outputs, not measurements. The degassing width is a closure assumption because the assignment does not completely specify it. Do not present the selected width as a measured property.

## Calculation route and equations

1. **Pressures and refrigerant-side states.** Saturation pressures follow the selected condenser and evaporator temperatures. The refrigerant flow is found from the evaporator balance:

   `Qe = mr (h4 - h3)`, so `mr = Qe / (h4 - h3)`.

2. **Solution flow rates.** The circulation ratio is determined from total mass and ammonia balances using the strong/weak solution concentrations and rectifier outlet concentration. Then `ms = f*mr` and `mw = (f - 1)*mr`.

3. **Pump and solution heat exchanger.** Pump work is calculated from solution specific volume, pressure rise, and pump efficiency. The heat exchanger uses an enthalpy-based effectiveness limit on both streams.

4. **Generator and rectifier.** The model finds the generator temperature from the solution bubble-point relation, calculates vapour composition and rectifier reflux, and then applies the generator energy balance to calculate generator heat input.

5. **COP and balances.** The cycle COP is

   `COP = Qe / (Qgen + Wp)`.

   The model also checks the generator, condenser, absorber, evaporator and solution heat-exchanger balances, the overall cycle energy residual, total mass balance, ammonia balance, pressure limit and exhaust-gas constraints.

## Default condenser-temperature sensitivity

At a fixed 5,000 kW cooling duty, 2 °C evaporator temperature and `dx = 0.16`, the current configured sweep predicts:

| Condenser temperature | COP | Generator heat input | High-side pressure | Exhaust-gas flow |
|---:|---:|---:|---:|---:|
| 30 °C | 0.6385 | 7.81 MW | 11.72 bar | 45.09 kg/s |
| 34 °C | 0.5905 | 8.44 MW | 13.19 bar | 51.40 kg/s |
| 38 °C | 0.5445 | 9.14 MW | 14.79 bar | 58.90 kg/s |
| 40 °C | 0.5222 | 9.53 MW | 15.64 bar | 63.18 kg/s |
| 44 °C | 0.4789 | 10.39 MW | 17.45 bar | 73.10 kg/s |
| 48 °C | 0.4366 | 11.39 MW | 19.42 bar | 85.36 kg/s |

The 48 °C point exceeds the project's 80 kg/s exhaust-gas limit, so it must be flagged as infeasible under that constraint even though its high-side pressure remains below 20 bar. The calculation illustrates the expected direction of the trend: raising the condenser temperature increases the required pressure lift and generator heat input, while reducing predicted COP.

## Closure sensitivity

At the default 40 °C condenser temperature, the solution degassing-width sweep predicts approximately:

| `dx` | COP | Generator temperature |
|---:|---:|---:|
| 0.12 | 0.5126 | 108.43 °C |
| 0.16 | 0.5222 | 117.13 °C |
| 0.20 | 0.5254 | 126.22 °C |

This sensitivity is a reminder that the assumption influences the result. Choose and justify the assumed value based on the course handout or reference data, not just because it returns a preferred COP.

## Change inputs and regenerate the plots

Edit `user_inputs.py` and run `python main.py`. The following figures/tables are regenerated from the edited values:

- `study_operating_point.csv`
- `study_state_points.csv`
- `study_component_balances.csv`
- `study_temperature_sweep.csv`
- `study_state_point_chart.png`
- `study_method_performance.png`
- `study_condenser_temperature_sweep.csv` / `study_condenser_temperature_plot.png`
- `study_solution_concentration_sweep.csv` / `study_solution_concentration_plot.png`

Changing cooling duty changes mass flow rates and heat duties. When temperatures and concentrations remain fixed, COP can stay approximately the same because all flow rates and duties scale with load. Static schematics and the equilibrium-property chart represent the process layout and property correlations, so their underlying relationships do not depend on cooling duty. The plot styling and schematic layout have been redrawn for readability without changing the calculated result tables. The seasonal/air-handler plots are computed from the weather and HVAC inputs defined in `air_handler.py`.

## Verification and limitations

- Both implementations use the same thermodynamic property functions. Agreement verifies code consistency, not independent correctness of all property correlations.
- Internal residuals near zero confirm that equations close with the current inputs; they do not prove the property model accurately predicts a real chiller.
- Check the category and meaning of every programmed engineering check in `results/model_check_summary.csv`, and document the assignment assumptions and source data in any submission.

## Reproduce

Run from this project folder:

```bash
python -m pip install -r requirements.txt
python main.py
python -m pytest -q tests
```

## Updated model correction

To match the reference method, calculate generator-top equilibrium with the incoming strong solution: `Ttop = T_bub(P_high, x_strong)` and `ygen = y_eq(P_high, x_strong)`. The generator outlet temperature remains `Tgen = T_bub(P_high, x_weak)`. Enthalpy reference offsets are constant per species and only cancel when the corresponding species mass balances close. The pure-ammonia PR reference is an anchor for the P&K scale, not a claim that saturated-vapour enthalpy equals zero.

Run `python main.py` to regenerate the six seasonal cases, the 40% fresh-air comparison, recommended-resized-condenser seasonal table/plot, and assumption sensitivity table. Inspect `results/model_check_summary.csv`: `Review` means a screening/design assumption is outside its prior band; it is not a passed check.


## Modelling limitations to keep visible

- The generator vapour equilibrium basis is a modelling choice. Both the updated solver and retained routine use the incoming strong solution, matching the supplied reference method. The main state-point and performance plots are saved under `results/plots/`, while supporting analysis tables are in `results/details/`. The main numerical output tables are retained directly in `results/`.
- The refrigerant loop uses a pseudo-pure-ammonia property approximation downstream of the 99.5 mass-% rectifier outlet. The residual water is not represented as a separate condenser/evaporator stream, so the state-1 to state-2 composition change is a known simplification.
- The direct Equation (8) vapour-composition correlation is now checked against the equilibrium implied by Equations (6) and (7). It is still an empirical fit and should only be used within the paper's stated range.
- `Review` results in `model_check_summary.csv` are intentionally retained for the model-dependent COP screen; they are not silently relabelled as passing results.
