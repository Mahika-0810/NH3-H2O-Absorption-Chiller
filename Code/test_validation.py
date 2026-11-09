"""test_validation.py - automatic assertions for the three required validation checks (Project 2, Section 2.4).

Run from this folder:   python -m pytest -q test_validation.py     (or: python test_validation.py)

Check 1  psychrometric mixing + wet-bulb   : w, h within 0.5 %, T within 0.1 K, Twb within 0.15 K
Check 2  ammonia real-gas EOS vs Table B.2.1: Psat within +-5 %, hfg within +-8 %, q_e and m_r inside the stated bands
Check 3  absorption-cycle balances          : mass / NH3 / energy residuals ~ 0, COP < COP_rev, COP band, COP falls with T_cond

Hand values are those printed in the project statement. Only public routines of the project modules are called.
"""
import math

import air_handler as air
import nh3_props as nh3
import cycle

P_ATM = 101.325


def pct(code, hand):
    return 100.0 * (code - hand) / hand


def assert_rel(name, hand, code, tol_pct):
    assert abs(pct(code, hand)) <= tol_pct, f"{name}: hand {hand}, code {code:.6g}, diff {pct(code, hand):+.3f} % (> {tol_pct} %)"


def assert_abs(name, hand, code, tol):
    assert abs(code - hand) <= tol, f"{name}: hand {hand}, code {code:.6g}, |diff| {abs(code - hand):.4f} (> {tol})"


# ---------------------------------------------------------------- Check 1
def _states():
    w_room = air.w_from_RH(22.0, 0.75)
    w_out = air.w_from_RH(40.0, 0.50)
    t_mix, w_mix, h_mix = air.mix(22.0, w_room, 40.0, w_out, 0.30)
    return w_room, w_out, t_mix, w_mix, h_mix


def test_check1_room_and_outdoor_states():
    w_room, w_out, _, _, _ = _states()
    assert_rel("Pv room 22C/75%", 1.984, 0.75 * air.psat(22.0), 0.5)
    assert_rel("w room", 0.01242, w_room, 0.5)
    assert_rel("h room", 53.71, air.h_moist(22.0, w_room), 0.5)
    assert_rel("Pv outdoor 40C/50%", 3.692, 0.50 * air.psat(40.0), 0.5)
    assert_rel("w outdoor", 0.02352, w_out, 0.5)
    assert_rel("h outdoor", 100.82, air.h_moist(40.0, w_out), 0.5)


def test_check1_adiabatic_mixing_70_30():
    _, _, t_mix, w_mix, h_mix = _states()
    assert_rel("w mixed", 0.01575, w_mix, 0.5)
    assert_rel("h mixed", 67.84, h_mix, 0.5)
    assert_abs("T mixed [C]", 27.48, t_mix, 0.1)


def test_check1_wet_bulb_temperatures():
    assert_abs("Twb summer day [C]", 30.3, air.Twb(40.0, air.w_from_RH(40.0, 0.50)), 0.15)
    assert_abs("Twb monsoon day [C]", 32.6, air.Twb(35.0, air.w_from_RH(35.0, 0.85)), 0.15)


# ---------------------------------------------------------------- Check 2
def test_check2_saturation_pressure_vs_table_b21():
    assert_rel("Psat NH3 0 C [kPa]", 429.6, nh3.sat(273.15)["P"], 5.0)
    assert_rel("Psat NH3 40 C [kPa]", 1554.9, nh3.sat(313.15)["P"], 5.0)


def test_check2_latent_heat_vs_table_b21():
    assert_rel("hfg NH3 0 C [kJ/kg]", 1261.8, nh3.sat(273.15)["hfg"], 8.0)
    assert_rel("hfg NH3 40 C [kJ/kg]", 1098.8, nh3.sat(313.15)["hfg"], 8.0)


def test_check2_trend_with_temperature_is_physical():
    s0, s40 = nh3.sat(273.15), nh3.sat(313.15)
    assert s40["P"] > s0["P"]            # Psat rises with T
    assert s40["hfg"] < s0["hfg"]        # hfg falls with T
    assert s0["Zg"] > s40["Zg"]          # vapour departs more from ideal as T rises toward Tc
    assert s40["Zg"] < 0.97              # real-gas effect is several percent (justifies the cubic EOS)


def test_check2_refrigerating_effect_and_mass_flow():
    q_e = nh3.sat(273.15)["hg"] - nh3.sat(313.15)["hf"]
    assert 1050.0 <= q_e <= 1100.0, f"q_e = {q_e:.1f} kJ/kg outside 1050-1100"
    m_r = 5000.0 / q_e
    assert 4.5 <= m_r <= 4.8, f"m_r = {m_r:.3f} kg/s outside 4.5-4.8"


# ---------------------------------------------------------------- Check 3
def test_check3_cycle_balances_close():
    r = cycle.solve_cycle(5000.0, 313.15)
    assert abs(r["bal_E_rel"]) < 1e-6, "energy balance (Qgen+Qevap+W = Qcond+Qabs+Qrect)"
    assert abs(r["bal_m"]) < 1e-9, "mass balance ms = mr + mw"
    assert abs(r["bal_NH3"]) < 1e-9, "NH3 balance ms*xs = mw*xw + mr*xr"
    assert abs(r["shx_E"]) < 1e-6, "SHX energy balance"


def test_check3_reversible_cop_bound():
    t_e, t_c, t_g = 275.15, 318.15, 403.15
    cop_rev = (t_e / (t_c - t_e)) * ((t_g - t_c) / t_g)
    assert_rel("COP_rev(2/45/130 C)", 1.35, cop_rev, 1.0)
    r = cycle.solve_cycle(5000.0, t_c)
    assert r["COP"] < r["COP_rev"]
    assert r["COP"] < cop_rev + 1e-12


def test_check3_cop_band_and_trend_with_condensing_temperature():
    cops = [cycle.solve_cycle(5000.0, t + 273.15)["COP"] for t in (30, 35, 40, 45)]
    assert all(b < a for a, b in zip(cops, cops[1:])), f"COP must fall as T_cond rises: {cops}"
    assert 0.45 <= cops[2] <= 0.65, f"COP(40 C) = {cops[2]:.3f} outside 0.45-0.65"
    assert 0.45 <= cops[3] <= 0.65, f"COP(45 C) = {cops[3]:.3f} outside 0.45-0.65"


def test_check3_cop_is_output_not_input():
    r = cycle.solve_cycle(5000.0, 313.15)
    assert math.isclose(r["COP"], r["Qevap"] / (r["Qgen"] + r["Wp"]), rel_tol=1e-12)


if __name__ == "__main__":      # lightweight runner if pytest is not installed
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn(); print("PASS", name)
            except AssertionError as exc:
                failed += 1; print("FAIL", name, "->", exc)
    raise SystemExit(1 if failed else 0)
