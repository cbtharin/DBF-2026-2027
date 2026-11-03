#!/usr/bin/env python3
"""Physics and scoring invariants for the DBF 2026-27 model.

These are not unit tests of the code, they are sanity checks on the model:
things that must be true of any aeroplane, so that if you change a coefficient
and one of these trips, you broke the physics rather than the spreadsheet.

    python selftest.py
"""
from __future__ import annotations

import math
import sys
from dataclasses import replace

from dbf2027 import (AIRFOILS, PROPS, Aircraft, Config, Design, evaluate_scores,
                     mission_2, mission_3, turn_capability)
from dbf2027.config import Environment
from dbf2027.components import MOTORS, PACKS
from dbf2027.numerics import (PCHIP, bisect, brent, gauss_legendre_5,
                              golden_max, nelder_mead, observed_order,
                              pareto_front, richardson_extrapolate, rk4_step,
                              scan_then_refine)
from dbf2027.propulsion import PropulsionSystem
from dbf2027.scoring import _score_from_ratios, rescore
from dbf2027.units import KG2LB, N2LBF

PASS, FAIL = 0, 0


def check(name: str, cond: bool, detail: str = "") -> None:
    global PASS, FAIL
    if cond:
        PASS += 1
        print("  PASS  %s" % name)
    else:
        FAIL += 1
        print("  FAIL  %s   %s" % (name, detail))


def section(title: str) -> None:
    print("\n" + title)
    print("-" * len(title))


def _spar_at(cfg: Config, n_ult: float) -> float:
    """Spar mass for the baseline wing at a given ultimate load factor."""
    c2 = Config()
    c2.build.ultimate_load_factor_override = n_ult
    return Aircraft(baseline(), c2).mass_breakdown.items["wing spar"]


def _isa_cfg() -> Config:
    """A sea-level ISA config, for showing what the site costs."""
    c = Config()
    c.env = Environment.at_site("isa")
    return c


def baseline(**kw) -> Design:
    args = dict(span_m=1.55, chord_m=0.27, airfoil="SD7062",
                motor="SCORPION_SII_4025_440", prop="APC_12x8", n_motors=1,
                pack="GA_6S3300_60C", n_systems=1,
                sensor_mass_kg=0.9, n_simulators=3)
    args.update(kw)
    return Design(**args)


def main() -> int:
    cfg = Config()

    # ------------------------------------------------------------------
    section("Numerics")
    r = bisect(lambda x: x * x - 2.0, 0.0, 2.0)
    check("bisect finds sqrt(2)", abs(r - math.sqrt(2)) < 1e-6, "got %.9f" % r)
    x, f = nelder_mead(lambda v: (v[0] - 3.0) ** 2 + (v[1] + 1.0) ** 2,
                       [0.0, 0.0], [0.5, 0.5], [-10, -10], [10, 10])
    check("nelder-mead finds (3,-1)",
          abs(x[0] - 3) < 1e-3 and abs(x[1] + 1) < 1e-3, "got %s" % x)
    front = pareto_front([(1, 1), (2, 2), (0, 3), (3, 0)])
    check("pareto drops dominated point", 0 not in front and len(front) == 3,
          "got %s" % front)

    # --- the numerical methods the mission model actually relies on ---
    # Brent stops at tol1 = 2*eps*|x| + xtol/2, so the achievable error is
    # set by the xtol you ask for; at xtol -> 0 it is machine precision.
    rb = brent(lambda x: x * x - 2.0, 0.0, 2.0)                 # default 1e-12
    check("brent meets its default tolerance",
          abs(rb - math.sqrt(2)) < 1e-12, "error %.2e" % abs(rb - math.sqrt(2)))
    rb2 = brent(lambda x: x * x - 2.0, 0.0, 2.0, xtol=0.0)
    check("brent reaches machine precision when asked",
          abs(rb2 - math.sqrt(2)) < 1e-15, "error %.2e" % abs(rb2 - math.sqrt(2)))
    # and it must beat bisection badly at equal iteration count
    eb = abs(bisect(lambda x: x * x - 2.0, 0.0, 2.0, tol=1e-9, maxiter=12)
             - math.sqrt(2))
    er = abs(brent(lambda x: x * x - 2.0, 0.0, 2.0, xtol=0.0, maxiter=12)
             - math.sqrt(2))
    check("brent converges faster than bisection", er < eb / 1000.0,
          "brent %.2e vs bisect %.2e" % (er, eb))

    # 5-point Gauss-Legendre is exact for polynomials up to degree 9
    gl = gauss_legendre_5(lambda x: 7 * x ** 9 - 3 * x ** 4 + x, 0.0, 1.3)
    exact = 0.7 * 1.3 ** 10 - 0.6 * 1.3 ** 5 + 0.5 * 1.3 ** 2
    check("gauss-legendre exact to degree 9", abs(gl - exact) < 1e-10,
          "error %.2e" % abs(gl - exact))

    # RK4 on dy/dt = y should show 4th-order convergence
    def integrate(h):
        y = [1.0]
        t = 0.0
        while t < 1.0 - 1e-12:
            y = rk4_step(lambda _t, yy: [yy[0]], t, y, h)
            t += h
        return y[0]
    e1, e2, e3 = (abs(integrate(h) - math.e) for h in (0.1, 0.05, 0.025))
    order = math.log(e1 / e2) / math.log(2.0)
    check("rk4 converges at 4th order", 3.7 < order < 4.3, "order %.2f" % order)
    check("rk4 error tiny at h=0.025", e3 < 1e-8, "error %.2e" % e3)

    # PCHIP must not overshoot monotone data (provided for callers who
    # load tabulated propeller data via load_props_csv)
    xs = [0.0, 1.0, 2.0, 3.0, 4.0]
    ys = [0.0, 0.0, 0.0, 1.0, 1.0]
    sp = PCHIP(xs, ys)
    samples = [sp(i * 0.02) for i in range(201)]
    check("pchip does not overshoot monotone data",
          min(samples) >= -1e-12 and max(samples) <= 1.0 + 1e-12,
          "range %.3f..%.3f" % (min(samples), max(samples)))
    check("pchip is monotone on monotone data",
          all(samples[i + 1] >= samples[i] - 1e-12
              for i in range(len(samples) - 1)))

    # observed_order should recover a known convergence rate
    def synth(h):
        return 5.0 + 3.0 * h ** 2
    p_obs = observed_order(synth(0.4), synth(0.2), synth(0.1), r=2.0)
    check("observed_order recovers p = 2", abs(p_obs - 2.0) < 1e-6,
          "got %.6f" % p_obs)
    ex = richardson_extrapolate(synth(0.2), synth(0.1), 2.0, r=2.0)
    check("richardson recovers the exact value", abs(ex - 5.0) < 1e-9,
          "got %.9f" % ex)

    # scan_then_refine must beat a bare scan on a peak between samples
    def peaked(x):
        return -((x - 0.3141) ** 2)
    xr, fr = scan_then_refine(peaked, 0.0, 1.0, n_scan=8, tol=1e-9)
    check("scan_then_refine locates a peak between samples",
          abs(xr - 0.3141) < 1e-5, "got %.6f" % xr)

    # ------------------------------------------------------------------
    section("Seeding the search with known designs")
    import json as _json
    import os as _os
    import tempfile as _tf
    from dataclasses import asdict as _asdict
    from dbf2027.optimizer import (DesignSpace as _DSs, _seed_one,
                                   _seed_usable, default_seed_designs,
                                   load_seed_designs)
    _ex = default_seed_designs()
    check("every search is seeded with the saved example designs",
          len(_ex) >= 5 and any(d.motor == "NEU_6521_3P5_362" for d in _ex),
          "%d seeds" % len(_ex))
    _tdir = _tf.mkdtemp()
    _json.dump(_asdict(baseline()), open(_os.path.join(_tdir, "a.json"),
                                         "w"))
    _json.dump([{"not": "a design"}], open(_os.path.join(_tdir, "b.json"),
                                           "w"))
    check("a seed folder yields its designs and skips everything else",
          len(load_seed_designs([_tdir])) == 1)
    _sp = _DSs()
    check("a seed may name a retired part, but not an unknown one or one a "
          "restricted search left out",
          _seed_usable(baseline(pack="GA_6S4500_60C"), _sp)
          and not _seed_usable(baseline(motor="NO_SUCH_MOTOR"), _sp)
          and not _seed_usable(baseline(motor="KDE_4014XF_380"),
                               _DSs(motors=["NEU_6521_3P5_362"])))
    _rec = Design(**_json.load(open(_os.path.join(
        _os.path.dirname(_os.path.abspath(__file__)), "examples",
        "optimised_v1.json"))))
    # 25% over: the recommended design's sensor is set by the hot day,
    # so a few percent more still flies on the nominal one
    _heavy = replace(_rec, sensor_mass_kg=_rec.sensor_mass_kg * 1.25,
                     sensor_m3_kg=_rec.sensor_m3_kg * 1.25)
    _sr = _seed_one((_heavy, cfg, 12.0, None, None, None))
    check("a seed too heavy for today's rules is trimmed until it flies",
          _sr is not None and _sr["feasible"]
          and _sr["design"].sensor_mass_kg < _heavy.sensor_mass_kg,
          "" if _sr is None else "%.2f of %.2f kg"
          % (_sr["design"].sensor_mass_kg, _heavy.sensor_mass_kg))
    from dbf2027.optimizer import _dedupe_configs, _seeded_refine_limit
    _mk = list(MOTORS)[:4]
    _coarse = [{"design": baseline(motor=m)} for m in _mk[:3]]
    _sd = [{"design": baseline(motor=_mk[0])},     # already sampled
           {"design": baseline(motor=_mk[3])}]     # new to the search
    _lim = _seeded_refine_limit(_sd, _coarse, 2)
    _got = {r["design"].motor for r in _dedupe_configs(_sd + _coarse, _lim)}
    check("seeds add refinement slots rather than taking the sampled "
          "configurations' slots",
          _lim == 3 and _got == {_mk[0], _mk[1], _mk[3]},
          "%d slots: %s" % (_lim, sorted(_got)))

    # ------------------------------------------------------------------
    section("Numerical conditioning (each against a 60-digit reference)")
    from decimal import Decimal as _D, getcontext as _gctx
    _gctx().prec = 60
    _tiny = _D(10) ** -58

    def _dsin(x):
        term, total, k = x, x, 1
        while abs(term) > _tiny:
            term *= -x * x / ((2 * k) * (2 * k + 1))
            total += term
            k += 1
        return total

    def _dcos(x):
        term, total, k = _D(1), _D(1), 1
        while abs(term) > _tiny:
            term *= -x * x / ((2 * k - 1) * (2 * k))
            total += term
            k += 1
        return total

    def _dasin(x):                       # small arguments only
        term, total, k = x, x, 0
        while abs(term) > _tiny:
            term *= x * x * (2 * k + 1) ** 2 / ((2 * k + 2) * (2 * k + 3))
            total += term
            k += 1
        return total

    _DPI = _D("3.14159265358979323846264338327950288419716939937510582097")
    from dbf2027.weights import _elliptic_moment, _elliptic_shear
    _worst = 0.0
    for _e in (0.5, 0.9, 0.99, 0.999999, 1.0 - 1e-10):
        _th = 2 * _dasin(((_D(1) - _D(_e)) / 2).sqrt())
        _m_ref = (_D(9) / 4 * _dsin(_th) + _dsin(3 * _th) / 4
                  - 3 * _th * _dcos(_th)) / 6
        _s_ref = (2 * _th - _dsin(2 * _th)) / _DPI
        _worst = max(_worst,
                     abs(_elliptic_moment(_e) / float(_m_ref) - 1.0),
                     abs(_elliptic_shear(_e) / float(_s_ref) - 1.0))
    check("the elliptic spar load keeps full precision out to the tip",
          _worst < 1e-12, "worst relative error %.1e" % _worst)
    from dbf2027.numerics import (gci as _gci, observed_order as _oo,
                                  solve_linear as _sl)
    _x = _sl([[1e10, 1e30], [1.0, 1.0]], [1e30, 2.0])
    check("the linear solve pivots on scaled rows (a near-zero pivot hidden "
          "by a 1e30 row scale)",
          abs(_x[0] - 1.0) < 1e-12 and abs(_x[1] - 1.0) < 1e-12, str(_x))
    from dbf2027.propulsion import winding_step as _ws
    _dx = 2e-8                           # dt / tau
    _t = _ws(0.0, 0.0, 50.0, _dx * 0.5 * 100.0, 0.5, 100.0)
    _t_ref = float(25 * (1 - (-_D(_dx)).exp()))
    check("a short winding-temperature step keeps its precision",
          abs(_t / _t_ref - 1.0) < 1e-13, "%.3e vs %.3e" % (_t, _t_ref))
    check("convergence order is not read off rounding noise",
          _oo(1.0, 1.0 + 2.0 ** -52, 1.0 + 2.0 ** -51) is None)
    _g = _gci(1.0 + 1e-6, 1.0, 1e-12)
    _g_ref = float(_D("1.25") * (_D(1.0 + 1e-6) - 1)
                   / ((_D(1e-12) * _D(2).ln()).exp() - 1))
    check("the grid convergence index stays exact for a low order",
          abs(_g / _g_ref - 1.0) < 1e-12, "%.6e vs %.6e" % (_g, _g_ref))
    _env = replace(Config().env, rho=1.225 * (1.0 - 1e-9))
    _h_ref = float(-((_D(1.225 * (1.0 - 1e-9)) / _D("1.225")).ln()
                     / _D("4.2559")).exp() + 1) / 6.87535e-6
    check("density altitude keeps its precision near sea level",
          abs(_env.density_altitude_ft / _h_ref - 1.0) < 1e-6,
          "%.6e vs %.6e ft" % (_env.density_altitude_ft, _h_ref))
    from dbf2027.mission import turn_capability as _tc
    _tac = Aircraft(baseline(), cfg)
    _turn = _tc(_tac.aero_clean, _tac.prop_sys, _tac.w_m2_n, 0.70, cfg)
    _n = _turn.load_factor
    check("the bank angle follows tan(bank) = sqrt(n^2 - 1)",
          abs(math.tan(math.radians(_turn.bank_deg))
              - math.sqrt((_n - 1.0) * (_n + 1.0))) < 1e-9,
          "n %.4f, bank %.2f deg" % (_n, _turn.bank_deg))

    # ------------------------------------------------------------------
    section("Propeller model")
    for key in ("APC_12x6", "APC_13x8", "APC_16x10"):
        p = PROPS[key]
        n = 8000.0 / 60.0
        t0 = p.thrust(1.225, n, 0.0)
        p0 = p.power(1.225, n, 0.0)
        area = math.pi * (p.diameter_m / 2) ** 2
        fom = t0 ** 1.5 / math.sqrt(2 * 1.225 * area) / p0
        check("%s static FOM in 0.30..0.90 (measured range)" % key,
              0.30 <= fom <= 0.90, "FOM %.3f" % fom)
        etas = []
        for i in range(1, 100):
            j = i * 0.01 * p.j_zero
            v = j * n * p.diameter_m
            pw = p.power(1.225, n, v)
            etas.append(p.thrust(1.225, n, v) * v / max(pw, 1e-9))
        check("%s peak eta in 0.60..0.92" % key, 0.60 <= max(etas) <= 0.92,
              "eta_max %.3f" % max(etas))
        check("%s thrust falls monotonically with speed" % key,
              all(p.thrust(1.225, n, v) >= p.thrust(1.225, n, v + 1.0) - 1e-9
                  for v in range(0, 30)))
        check("%s makes no thrust past J_zero" % key,
              p.thrust(1.225, n, 1.05 * p.j_zero * n * p.diameter_m) == 0.0)

    # ------------------------------------------------------------------
    # --- measured propellers (UIUC) ---------------------------------
    _meas = [k for k, p_ in PROPS.items() if p_.measured]
    check("the catalogue propellers are UIUC wind-tunnel data",
          len(_meas) >= 25 and all(PROPS[k].table for k in PROPS),
          "%d of %d measured" % (len(_meas), len(PROPS)))
    _pm = PROPS["APC_14x10"]
    check("a measured table is read exactly at its points (CP)",
          all(abs(_pm.cp(j_) - cp_) < 1e-12 for j_, _c, cp_ in _pm.table))
    check("measured CT never rises with J",
          all(_pm.ct(0.02 * i) >= _pm.ct(0.02 * (i + 1)) - 1e-12
              for i in range(60)))
    check("and reaches zero at the zero-thrust advance ratio",
          _pm.ct(_pm.j_zero) == 0.0 and _pm.ct(0.9 * _pm.j_zero) > 0.0)

    section("Motor / propulsion")
    ps = PropulsionSystem(MOTORS["SCORPION_SII_4025_440"], PROPS["APC_12x8"], 1,
                          PACKS["GA_6S3300_60C"], cfg.env)
    static = ps.solve(0.0, 1.0, 1.0)
    fast = ps.solve(30.0, 1.0, 1.0)
    check("thrust decreases with airspeed",
          fast.thrust_total_n < static.thrust_total_n,
          "%.1f vs %.1f N" % (fast.thrust_total_n, static.thrust_total_n))
    check("current decreases as the prop unloads",
          fast.motor_current_a < static.motor_current_a,
          "%.1f vs %.1f A" % (fast.motor_current_a, static.motor_current_a))
    check("rpm rises as the prop unloads", fast.rpm > static.rpm,
          "%.0f vs %.0f" % (fast.rpm, static.rpm))
    half = ps.solve(20.0, 1.0, 0.55)
    full = ps.solve(20.0, 1.0, 1.0)
    check("less throttle gives less thrust", half.thrust_total_n < full.thrust_total_n)
    low_soc = ps.solve(20.0, 0.15, 1.0)
    check("a flat pack gives less thrust", low_soc.thrust_total_n < full.thrust_total_n,
          "%.1f vs %.1f N" % (low_soc.thrust_total_n, full.thrust_total_n))
    # Battery sag.  The loaded voltage must sit below open circuit by
    # exactly R * I (Kirchhoff on the pack side), and the lost volts must
    # cost thrust against a pack with no resistance.
    check("the pack sags under load: V = Voc - R I",
          abs(static.pack_voltage_v
              - (ps.pack.ocv(1.0) - ps.r_sys * static.pack_current_a)) < 1e-9
          and static.pack_voltage_v < ps.pack.ocv(1.0) - 0.1,
          "%.3f V loaded vs %.3f V open circuit at %.1f A"
          % (static.pack_voltage_v, ps.pack.ocv(1.0), static.pack_current_a))
    stiff = PropulsionSystem(MOTORS["SCORPION_SII_4025_440"], PROPS["APC_12x8"],
                             1, replace(PACKS["GA_6S3300_60C"],
                                        r_cell_mohm=1e-6),
                             cfg.env, r_harness_ohm=0.0)
    st_op = stiff.solve(0.0, 1.0, 1.0)
    check("sag costs static thrust", static.thrust_total_n < st_op.thrust_total_n,
          "%.2f vs %.2f N with no resistance"
          % (static.thrust_total_n, st_op.thrust_total_n))
    check("sag costs energy: power is drawn at open-circuit volts",
          abs(static.p_elec_w - ps.pack.ocv(1.0) * static.total_current_a
              / ps.esc.efficiency) < 1e-6,
          "%.1f W" % static.p_elec_w)
    check("with no resistance the loaded and open-circuit volts agree",
          abs(st_op.pack_voltage_v - stiff.pack.ocv(1.0)) < 1e-4,
          "%.4f vs %.4f V" % (st_op.pack_voltage_v, stiff.pack.ocv(1.0)))
    check("a measured cell resistance overrides the estimate",
          abs(replace(PACKS["GA_6S4500_60C"], r_cell_mohm=4.0).r_internal_ohm
              - 0.024) < 1e-12)
    # (up to ~170 mOhm for a 12S 2000 mAh pack with an honest, low 27C
    # label - MaxAmps - where the C-scaled estimate is conservative)
    check("pack resistance estimates are physical (1-200 mOhm)",
          all(0.001 <= pk.r_internal_ohm <= 0.200 for pk in PACKS.values()),
          "%.1f-%.1f mOhm" % (1e3 * min(pk.r_internal_ohm for pk in PACKS.values()),
                              1e3 * max(pk.r_internal_ohm for pk in PACKS.values())))
    check("higher C rating means lower estimated resistance",
          replace(PACKS["GA_6S4500_60C"], c_rating=25.0).r_internal_ohm
          > PACKS["GA_6S4500_60C"].r_internal_ohm)
    check("electrical power exceeds shaft power",
          static.p_elec_w > static.p_shaft_w,
          "%.0f vs %.0f W" % (static.p_elec_w, static.p_shaft_w))

    # The operating point is the root of motor torque minus propeller
    # torque.  Verify it really is a root, at every flight condition: a
    # bracketing or algebra mistake would be silent, shifting every thrust
    # value by a few percent without ever raising an error.  The current is
    # rebuilt from the circuit, not from the solver's closed form: the
    # motor sees throttle x the LOADED pack voltage it reported.
    worst_res = 0.0
    worst_scale = 1e-12
    worst_kvl = 0.0
    for v in (0.0, 9.0, 21.0, 34.0, 45.0):
        for soc_ in (1.0, 0.62, 0.21):
            for thr_ in (1.0, 0.74, 0.31):
                op = ps.solve(v, soc_, thr_)
                n_ = op.rpm / 60.0
                if n_ <= 1e-6:
                    continue
                worst_kvl = max(worst_kvl, abs(
                    op.pack_voltage_v
                    - (ps.pack.ocv(soc_) - ps.r_sys * op.pack_current_a)))
                vt = thr_ * op.pack_voltage_v
                cur = (vt - 60.0 * n_ / ps.motor.kv) / ps.rm
                # the no-load current at this speed (Motor.i0_at)
                q_motor = ps.kt * (cur - ps.motor.i0_at(n_, ps.i0_speed_share))
                jj = v / (n_ * ps.prop.diameter_m)
                q_prop = (ps.prop.cp(jj, 60.0 * n_) * ps.env.rho * n_ * n_
                          * ps.prop.diameter_m ** 5 / (2.0 * math.pi))
                worst_res = max(worst_res, abs(q_motor - q_prop))
                worst_scale = max(worst_scale, abs(q_prop))
    check("the torque balance is satisfied at every operating point",
          worst_res / worst_scale < 1e-9,
          "worst residual %.2e N*m (%.2e relative)"
          % (worst_res, worst_res / worst_scale))
    check("the pack circuit closes at every operating point",
          worst_kvl < 1e-9, "worst %.2e V" % worst_kvl)

    # Current is the thing the 100 A fuse is written against, so the two
    # current figures have to mean what the rules mean by them.
    two_sys = PropulsionSystem(MOTORS["SCORPION_SII_4025_440"], PROPS["APC_12x8"],
                               2, PACKS["GA_6S3300_60C"], cfg.env,
                               n_systems=2)
    one_sys = PropulsionSystem(MOTORS["SCORPION_SII_4025_440"], PROPS["APC_12x8"],
                               2, PACKS["GA_6S3300_60C"], cfg.env,
                               n_systems=1)
    r2_, r1_ = two_sys.rating_check(), one_sys.rating_check()
    check("each of two packs carries half the total current",
          abs(r2_["pack_current_a"] - 0.5 * r2_["total_current_a"]) < 1e-9,
          "%.1f of %.1f A" % (r2_["pack_current_a"], r2_["total_current_a"]))
    # Two packs each carry half the current, so each sags half as far and
    # the motors pull a little more in total - by a few percent, no more.
    # (8% on this pair since the propeller is read at its own rpm: the
    # faster-turning prop sits at a higher Reynolds number.)
    check("two packs sag less, so the total current is slightly higher",
          r1_["total_current_a"] < r2_["total_current_a"]
          < 1.10 * r1_["total_current_a"],
          "%.1f vs %.1f A" % (r2_["total_current_a"], r1_["total_current_a"]))
    _stiff_pk = replace(PACKS["GA_6S3300_60C"], r_cell_mohm=1e-6)
    _r2s = PropulsionSystem(MOTORS["SCORPION_SII_4025_440"], PROPS["APC_12x8"],
                            2, _stiff_pk, cfg.env, n_systems=2,
                            r_harness_ohm=0.0).rating_check()
    _r1s = PropulsionSystem(MOTORS["SCORPION_SII_4025_440"], PROPS["APC_12x8"],
                            2, _stiff_pk, cfg.env, n_systems=1,
                            r_harness_ohm=0.0).rating_check()
    check("with no resistance the pack count does not change the total",
          abs(_r2s["total_current_a"] - _r1s["total_current_a"]) < 1e-3,
          "%.4f vs %.4f A" % (_r2s["total_current_a"], _r1s["total_current_a"]))
    check("two systems carry twice the energy",
          abs(two_sys.total_watt_hours - 2.0 * one_sys.total_watt_hours) < 1e-9)
    check("the fuse is never rated above 100 A",
          all(PropulsionSystem(MOTORS["SCORPION_SII_4025_440"], PROPS["APC_12x8"],
                               1, pk, cfg.env).fuse_rating_a <= 100.0
              for pk in PACKS.values()))

    # and the returned operating point must satisfy the torque balance
    op_chk = ps.solve(18.0, 0.8, 1.0)
    n_chk = op_chk.rpm / 60.0
    q_motor = ps.kt * (op_chk.motor_current_a
                       - ps.motor.i0_at(n_chk, ps.i0_speed_share))
    q_prop = (ps.prop.cp(18.0 / (n_chk * ps.prop.diameter_m), 60.0 * n_chk)
              * ps.env.rho
              * n_chk * n_chk * ps.prop.diameter_m ** 5 / (2.0 * math.pi))
    check("torque balance residual is at floating-point noise",
          abs(q_motor - q_prop) < 1e-10,
          "residual %.2e N*m" % abs(q_motor - q_prop))

    # ------------------------------------------------------------------
    section("Aerodynamics")
    ac = Aircraft(baseline(), cfg)
    aero, w = ac.aero_clean, ac.w_m2_n
    # Not exactly sqrt(W) any more: CLmax falls with Reynolds number and the
    # trim download does not scale with weight.  Both are deliberate, so the
    # invariant is "close to sqrt(W), and never better than it".
    # With the airfoil data referenced at Re 500k, the baseline stalls
    # below its reference Reynolds number and at 4x the weight reaches the
    # CLmax cap, so CLmax rises with weight and the ratio sits a little
    # further under 2 than it did at a 300k reference.
    ratio = aero.v_stall(4 * w) / aero.v_stall(w)
    check("stall speed scales roughly as sqrt(W)", 1.85 < ratio < 2.10,
          "ratio %.4f" % ratio)
    check("heavier always stalls faster", aero.v_stall(2 * w) > aero.v_stall(w))
    check("CLmax falls at low Reynolds number",
          aero.cl_max_at_re(120e3) < aero.cl_max_at_re(400e3),
          "%.3f vs %.3f" % (aero.cl_max_at_re(120e3), aero.cl_max_at_re(400e3)))
    check("a cambered section needs a tail download",
          aero.tail_load(22.0, w) < 0.0, "%.2f N" % aero.tail_load(22.0, w))
    check("trim download makes the wing carry more than the weight",
          aero.cl_required(22.0, w) * aero.q(22.0) * aero.S > w)
    # The empennage is part of the lift balance: wing + tail = total lift.
    lt = aero.tail_load(22.0, w)
    lw = aero.cl_required(22.0, w) * aero.q(22.0) * aero.S
    check("wing lift + tail lift = total lift", abs((lw + lt) - w) < 1e-6,
          "%.4f + %.4f vs %.4f" % (lw, lt, w))
    check("tail load changes with load factor",
          abs(aero.tail_load(22.0, 3 * w) - aero.tail_load(22.0, w)) > 1e-6,
          "1g %.2f, 3g %.2f" % (aero.tail_load(22.0, w),
                                aero.tail_load(22.0, 3 * w)))
    check("tail lift coefficient stays inside the tail's own stall",
          abs(aero.tail_cl(22.0, w)) < cfg.build.tail_cl_max,
          "CL_t %.3f" % aero.tail_cl(22.0, w))
    check("trim costs drag",
          aero.drag(22.0, w) > aero.q(22.0) * aero.S
          * aero.cd_from_cl(w / (aero.q(22.0) * aero.S), 22.0) * 0.999)
    # With the CG exactly at the wing AC only the camber term survives, and
    # the tail load then does not scale with load factor.
    import copy as _copy
    a0 = _copy.copy(aero)
    a0.set_cg(0.25)
    w0 = w
    check("with CG at the wing AC the tail load is load-factor independent",
          abs(a0.tail_load(22.0, w0) - a0.tail_load(22.0, 3 * w0)) < 1e-9)
    # (a CG set aft of the AC: with the tails sized from their
    # requirements, m111, the baseline's own CG can sit just ahead of it)
    a1 = _copy.copy(aero)
    a1.set_cg(0.35)
    check("moving the CG aft unloads the tail",
          abs(a1.tail_load(22.0, w)) < abs(a0.tail_load(22.0, w)),
          "%.2f vs %.2f N" % (a1.tail_load(22.0, w), a0.tail_load(22.0, w)))
    check("a symmetric-ish section needs less trim than a cambered one",
          abs(Aircraft(baseline(airfoil="AG35"), cfg).aero_clean.tail_load(
              22.0, w))
          < abs(Aircraft(baseline(airfoil="S1223"), cfg).aero_clean.tail_load(
              22.0, w)))
    check("max_lift is consistent with the trim relation",
          abs((lambda v: (lambda Lm: (Lm - aero.tail_load(v, Lm))
                          - aero.q(v) * aero.S * aero.cl_max_at_v(v))
               (aero.max_lift(v)))(22.0)) < 1e-6)
    check("load factor raises stall speed",
          aero.v_stall(w, 2.0) > aero.v_stall(w))
    drags = [aero.drag(v, w) for v in (14, 18, 22, 26, 30, 36, 42)]
    check("drag polar has an interior minimum",
          min(drags) < drags[0] and min(drags) < drags[-1],
          "%s" % [round(d, 1) for d in drags])
    check("best L/D between 6 and 20",
          6.0 < aero.ld_ratio(aero.best_ld_speed(w), w) < 20.0,
          "%.1f" % aero.ld_ratio(aero.best_ld_speed(w), w))
    check("CD0 between 0.02 and 0.07", 0.02 < aero.cd0(20.0) < 0.07,
          "%.4f" % aero.cd0(20.0))
    check("induced drag grows with weight",
          aero.drag(25.0, 2 * w) > aero.drag(25.0, w))
    check("higher aspect ratio lowers induced drag",
          Aircraft(baseline(span_m=2.2, chord_m=0.19), cfg).aero_clean.k_induced
          < aero.k_induced)
    check("high-lift airfoil gives a higher CLmax",
          Aircraft(baseline(airfoil="S1223"), cfg).aero_clean.cl_max
          > Aircraft(baseline(airfoil="MH32"), cfg).aero_clean.cl_max)

    # ------------------------------------------------------------------
    section("Weights and structure")
    light = Aircraft(baseline(n_simulators=0), cfg)
    heavy = Aircraft(baseline(n_simulators=8), cfg)
    check("more payload means more gross weight",
          heavy.mtow_m2_kg > light.mtow_m2_kg)
    check("more payload means a heavier structure",
          heavy.empty_mass_kg > light.empty_mass_kg,
          "%.2f vs %.2f kg" % (heavy.empty_mass_kg, light.empty_mass_kg))
    # Same wing area, so at minimum gauge the two spars weigh the same
    # (the tube follows the section, and the section follows the chord).
    # Sized by LOAD, the longer, thinner wing needs far more spar.  (This
    # used to compare the gauge-limited masses, and passed only because
    # 0.2 x 2.1 rounds a hair above 0.35 x 1.2 in binary.)
    _ls = Config()
    _ls.build.spar_min_wall_m = 1e-9
    short = Aircraft(baseline(span_m=1.2, chord_m=0.35), _ls)
    long_ = Aircraft(baseline(span_m=2.1, chord_m=0.20), _ls)
    check("longer span means a heavier load-sized spar",
          long_.mass_breakdown.items["wing spar"]
          > 2.0 * short.mass_breakdown.items["wing spar"],
          "%.1f vs %.1f g" % (1e3 * long_.mass_breakdown.items["wing spar"],
                              1e3 * short.mass_breakdown.items["wing spar"]))
    check("mass buildup is self-consistent",
          abs(ac.mtow_m2_kg - (ac.empty_mass_kg + ac.m2_payload_kg)) < 1e-6)
    check("M3 gross is lighter than M2 gross", ac.mtow_m3_kg < ac.mtow_m2_kg)
    check("empty weight in a believable band (4..20 lb)",
          4.0 < ac.empty_mass_kg * KG2LB < 20.0,
          "%.2f lb" % (ac.empty_mass_kg * KG2LB))

    # ------------------------------------------------------------------
    section("Mission")
    r2 = mission_2(ac, cfg, dx=3.0)
    r3 = mission_3(ac, cfg, dx=3.0)
    check("Mission 2 flies the required laps",
          r2.ok and r2.laps == cfg.rules.m2_required_laps, r2.reason)
    check("Mission 2 fits the window", r2.time_s <= cfg.rules.mission_window_s,
          "%.1f s" % r2.time_s)
    check("Mission 3 fits the window", r3.time_s <= cfg.rules.mission_window_s,
          "%.1f s" % r3.time_s)
    check("lap times converge", len(r2.lap_times) >= 3
          and abs(r2.lap_times[-1] - r2.lap_times[-2]) < 1.0)
    check("corner speed sits between stall and Vmax",
          r2.v_stall < r2.v_turn < r2.v_max,
          "%.1f / %.1f / %.1f" % (r2.v_stall, r2.v_turn, r2.v_max))
    check("energy used does not exceed the pack",
          r2.energy_j < ac.pack.watt_hours * 3600.0,
          "%.1f Wh of %.1f" % (r2.energy_j / 3600, ac.pack.watt_hours))
    check("the deployed sensor costs drag",
          ac.aero_deployed.cd0(25.0) > ac.aero_clean.cd0(25.0))

    heavier = Aircraft(baseline(n_simulators=8), cfg)
    r2h = mission_2(heavier, cfg, dx=3.0)
    if r2h.ok:
        check("more payload means a longer Mission 2", r2h.time_s > r2.time_s,
              "%.1f vs %.1f s" % (r2h.time_s, r2.time_s))
        check("more payload means a longer takeoff roll",
              r2h.takeoff_m > r2.takeoff_m,
              "%.1f vs %.1f m" % (r2h.takeoff_m, r2.takeoff_m))
    else:
        check("8 simulators is correctly rejected as infeasible", True)

    windy = Config()
    windy.env.wind_mps = 12.0
    r2w = mission_2(Aircraft(baseline(), windy), windy, dx=3.0)
    check("wind costs lap time (or the mission)",
          (not r2w.ok) or r2w.time_s > r2.time_s,
          "%.1f vs %.1f s" % (r2w.time_s, r2.time_s))

    # ------------------------------------------------------------------
    section("Scoring")
    card = evaluate_scores(ac, cfg, dx=3.0)
    check("M1 is 1.0 when the aeroplane flies", card.m1 == 1.0)
    check("M2 is at least 1.0 on a successful flight", card.m2 >= 1.0)
    check("M3 is at least 2.0 on a successful flight", card.m3 >= 2.0)
    check("M2 never exceeds 2.0", card.m2 <= 2.0 + 1e-12, "%.4f" % card.m2)
    check("M3 never exceeds 3.0", card.m3 <= 3.0 + 1e-12, "%.4f" % card.m3)
    check("M2 raw equals scored weight over time",
          abs(card.m2_raw_lb_per_s
              - (ac.m2_scored_kg * KG2LB) / card.r2.time_s) < 1e-9)
    check("M3 raw equals laps times sensor weight",
          abs(card.m3_raw_lap_lb
              - card.r3.laps * ac.design.sensor_mass_kg * KG2LB) < 1e-9)
    sat = rescore(card, card.m2_raw_lb_per_s * 0.5,
                  card.m3_raw_lap_lb * 0.5, mode="field")
    check("in field mode, beating the field saturates the ratio at 1.0",
          abs(sat.m2 - 2.0) < 1e-9 and abs(sat.m3 - 3.0) < 1e-9)
    uncapped = rescore(card, card.m2_raw_lb_per_s * 0.5,
                       card.m3_raw_lap_lb * 0.5, mode="self")
    check("in self mode the ratio is not capped",
          uncapped.m2 > 2.0 and uncapped.m3 > 3.0)
    check("a weaker reference raises the score",
          rescore(card, 0.20, 80.0).total_mission < card.total_mission)

    broken = Aircraft(baseline(pack="TP_6S2200_25C",
                               n_simulators=12), cfg)
    bad = evaluate_scores(broken, cfg, dx=6.0)
    check("an infeasible aeroplane scores 0 on the missions it cannot fly",
          bad.m2 == 0.0 or bad.m3 == 0.0 or not bad.feasible, bad.reason)

    # --- lap counts and mission order (rules 3.3.3a, 3.3.3b, 3.3.3.1/2) --
    # Lap counts must come from the rules config, not literals: a knob in
    # the Advanced panel that changes nothing is worse than no knob.
    import copy
    from dbf2027.mission import mission_1
    fewer = copy.deepcopy(cfg)
    fewer.rules.m1_required_laps = 2
    fewer.rules.m2_required_laps = 4
    r1_3, r1_2 = mission_1(ac, cfg, dx=3.0), mission_1(ac, fewer, dx=3.0)
    check("Mission 1 flies the configured lap count, not a literal 3",
          r1_3.laps == cfg.rules.m1_required_laps and r1_2.laps == 2
          and r1_2.time_s < r1_3.time_s,
          "%d laps %.1f s vs %d laps %.1f s"
          % (r1_3.laps, r1_3.time_s, r1_2.laps, r1_2.time_s))
    r2_4 = mission_2(ac, fewer, dx=3.0)
    check("Mission 2 flies the configured lap count",
          r2_4.laps == 4 and r2_4.time_s < r2.time_s,
          "%d laps %.1f s vs %.1f s" % (r2_4.laps, r2_4.time_s, r2.time_s))
    check("the rules say 3 laps for Mission 1 and 5 for Mission 2",
          cfg.rules.m1_required_laps == 3 and cfg.rules.m2_required_laps == 5)

    # 40 laps cannot fit a 5 minute window, so Mission 2 fails while
    # Mission 3 - which flies laps for as long as the window lasts - is
    # untouched.  Mission 3 must then score nothing: it may not be flown
    # until Mission 2 has scored.
    no_m2 = copy.deepcopy(cfg)
    no_m2.rules.m2_required_laps = 40
    locked = evaluate_scores(ac, no_m2, dx=3.0)
    check("Mission 3 is locked until Mission 2 succeeds (rule 3.3.3a)",
          not locked.r2.ok and locked.r3.ok and locked.m3 == 0.0
          and locked.m3_raw_lap_lb == 0.0 and "3.3.3a" in locked.reason,
          "M2 ok=%s, M3 flew %d laps, scored %.3f"
          % (locked.r2.ok, locked.r3.laps, locked.m3))
    check("the Ground Mission is not gated by the flight missions",
          locked.gm == card.gm, "%.3f vs %.3f" % (locked.gm, card.gm))
    check("Mission 1 still scores when flown on its own",
          locked.m1 == 1.0 and locked.r1.ok)

    # A Mission 1 that cannot be completed, with a Mission 2 that can:
    # the waiver means the successful Mission 2 earns Mission 1 as well.
    waived = copy.deepcopy(cfg)
    waived.rules.m1_required_laps = 40
    w = evaluate_scores(ac, waived, dx=3.0)
    check("a successful Mission 2 also earns Mission 1 (rule 3.3.3b)",
          not w.r1.ok and w.r2.ok and w.m1 == 1.0,
          "M1 ok=%s M2 ok=%s M1 score %.1f" % (w.r1.ok, w.r2.ok, w.m1))

    # The score panel must not describe a failed Mission 2 as a flight:
    # it used to print "time for 5 laps" (the time the battery ran out)
    # and "M2 = 1 + 0.0000 = 0.000".
    from dbf2027.scoring import score_explanation
    _blocks = {b["title"]: b for b in score_explanation(ac, no_m2, locked)}
    _b2 = _blocks["Mission 2 - delivery flight"]
    check("a failed Mission 2 is reported as failed, with no lap time",
          "did not complete" in _b2["result"]
          and not any(l.startswith("time for") for l, _, _ in _b2["rows"]),
          _b2["result"])
    check("and Mission 3 is shown as locked, not as a score",
          "locked" in _blocks["Mission 3 - sensor flight"]["result"])

    # ------------------------------------------------------------------
    section("Weather - a propeller per condition (rule 3.2.1e)")
    from dbf2027.optimizer import weather_plan, OptimisationResult
    # A Hacker A40-14L V4 on a 16x8 and a 6S 4500: inside the motor's
    # power rating at the nominal 28 C, but over it in the denser 12 C
    # morning air with its own propeller (UIUC-measured data).  (A 14x12
    # until the propellers were read at their own rpm, 2026-10-03, which
    # left the 14x12 inside the rating on the morning too.  It was a
    # SunnySky X4120 until every flight was charged at full throttle,
    # 2026-10-02, which that motor's 14x12 cannot last five laps on; the
    # fixture before that broke a 100 A TOTAL across two systems, which is
    # not a rule and has been removed.)
    hot_draw = Design(span_m=1.8288, chord_m=0.356633, airfoil="SD7062",
                      motor="HACKER_A40_14L_V4_355", prop="APC_16x8",
                      n_motors=1, pack="GA_6S4500_60C", n_systems=1,
                      sensor_mass_kg=8.0, sensor_m3_kg=7.2,
                      n_simulators=0, taper=0.43396, drop_height_in=60.0)
    check("it is legal at the nominal 28 C",
          Aircraft(hot_draw, cfg).static_checks()["ok"])
    _morning = copy.deepcopy(cfg)
    _morning.env = Environment.at_site("tucson-morning")
    check("that design really is over the limit on its own propeller",
          not Aircraft(hot_draw, _morning).static_checks()["ok"])
    _wp = weather_plan(hot_draw, cfg, dx=12.0)
    _am = _wp["presets"]["tucson-morning"]
    check("every day's propeller keeps the peak inside each fuse",
          all(e["peak_pack_a"] <= e["fuse_a"] + 1e-9
              for e in _wp["presets"].values() if e["legal"]),
          ", ".join("%s %.0f/%.0f A" % (p, e["peak_pack_a"], e["fuse_a"])
                    for p, e in _wp["presets"].items()))
    check("on the morning a smaller propeller is fitted, legally",
          _am["swapped"] and _am["legal"] and _am["prop"] != hot_draw.prop
          and PROPS[_am["prop"]].diameter_m <= PROPS[hot_draw.prop].diameter_m,
          "%s, %s" % (_am["prop"], _am["reason"][:40]))
    check("and with it the design flies every preset",
          _wp["robust"] and _wp["code"] == "swap",
          "code %s" % _wp["code"])
    # (the plain baseline misses the window by 2 s on the hot day with its
    # own propeller now that sag is modelled; two simulators is legal
    # everywhere)
    _bp = weather_plan(baseline(n_simulators=2), cfg, dx=12.0)
    check("a design legal everywhere keeps its own propeller",
          _bp["code"] == "ok" and not _bp["swaps"], "code %s" % _bp["code"])
    _rows = [{"design": None, "robust": False}, {"design": None, "robust": True}]
    check("the recommendation skips a top score that cannot fly every preset",
          OptimisationResult(_rows, [], 0, 0).best_row is _rows[1])

    # ------------------------------------------------------------------
    section("Airfoil data")
    from dbf2027.components import AIRFOILS as _AF, Airfoil as _Airfoil
    check("the catalogue is 32 airfoiltools sections at Re 500k",
          len(_AF) == 32 and all(a.polar and a.re_ref == 500000.0
                                 and "airfoiltools" in a.note
                                 for a in _AF.values()), "%d" % len(_AF))
    _a = _AF["E205"]
    check("the model reads the polar itself - exact at every table point",
          all(abs(_a.cd_excess(cl) - (cd - _a.cd_min)) < 1e-12
              for cl, cd in _a.polar))
    _top = _a.polar[-1][0]
    check("drag keeps rising past the end of the polar",
          _a.cd_excess(_top + 0.2) > _a.cd_excess(_top)
          > _a.cd_excess(_a.cl_at_cd_min))
    _hand = _Airfoil("HAND", _a.cl_max_2d, _a.cd_min, _a.cl_at_cd_min,
                     _a.k_drag, _a.cl_break, _a.t_over_c, _a.cm_ac)
    check("a hand-entered airfoil falls back to the fitted bucket",
          not _hand.polar and _hand.cd_excess(_a.cl_at_cd_min) == 0.0
          and _hand.cd_excess(_a.cl_at_cd_min + 0.3) > 0.0)
    # --- motor data ----------------------------------------------------
    # The previous catalogue was unsourced and wrong where it mattered (a
    # "C-4120/12 465KV" that does not exist, a Scorpion 150 g too heavy).
    import re as _re
    from dbf2027.components import RETIRED as _RET, searchable as _srch
    _sm = [MOTORS[k] for k in _srch(MOTORS)]
    # pruned 2026-09-29 (results/catalogue_prune_2026-09-29) from 129 to
    # the ones that earn a place at some motor count; 59 low-Kv motors
    # added 2026-10-03 for the big high-pitch props (m65)
    check("the searched motors are pruned but still cover 1-4 motor "
          "designs, every motor sourced to a maker page",
          30 <= len(_sm) <= 130
          and all("http" in m.source and "VERIFY" not in m.source
                  for m in MOTORS.values())
          and min(m.mass_kg for m in _sm) < 0.20
          and max(m.mass_kg for m in _sm) > 0.6
          and max(m.max_cells for m in _sm) >= 12,
          "%d motors searched" % len(_sm))
    from dbf2027.optimizer import DesignSpace as _DSr
    _spr = _DSr()
    check("retired parts stay loadable but are left out of the search",
          _RET and all(k in MOTORS or k in PACKS for k in _RET)
          and not set(_RET) & (set(_spr.motors) | set(_spr.packs))
          and len(_spr.motors) + len(_RET) >= len(MOTORS))
    _kind = lambda p: (p.series, "LiHV" if p.v_full_cell > 4.3
                       else p.chemistry)
    check("every cell count and cell type in the catalogue keeps a "
          "searched pack",
          {_kind(p) for p in PACKS.values()}
          == {_kind(PACKS[k]) for k in _srch(PACKS)})
    # every researched motor passes the physics the research was checked
    # against: efficiency at its rated current, copper loss per gram, and
    # a motor constant on the makers' common Km-vs-mass trend
    _bad = []
    for _m in MOTORS.values():
        _v = _m.max_cells * 3.7
        _i = min(_m.i_max_a, _m.p_max_w / _v)
        _eta = (1 - _m.i0 / _i) * (1 - _i * _m.rm / _v)
        _km = 60.0 / (2 * math.pi * _m.kv) / math.sqrt(_m.rm)
        # (short-burst ratings - SunnySky's 30 s, Hacker's 15 s - reach
        # ~1.3 W/g of copper loss and ~77% efficiency at the rated current;
        # the winding-temperature model is what limits them in a mission)
        if not (0.75 <= _eta <= 0.96 and _i * _i * _m.rm / (_m.mass_kg * 1000)
                < 1.4 and 0.17 <= _km / _m.mass_kg ** 0.9 <= 0.55):
            _bad.append(_m.key)
    check("every motor passes the efficiency, copper-loss and Km checks",
          not _bad, str(_bad[:5]))
    _kv_bad = []
    for _m in MOTORS.values():
        _n = (_re.search(r"(?i)kv ?(\d{2,4})\b", _m.name)
              or _re.search(r"(?i)\b(\d{2,4}) ?kv\b", _m.name))
        if _n and abs(float(_n.group(1)) - _m.kv) > 0.05 * _m.kv:
            _kv_bad.append(_m.key)
    check("every motor's Kv matches the Kv in its own name", not _kv_bad,
          str(_kv_bad))
    # motor constant per kg tracks size closely; a resistance misread by
    # 10x or 1000x lands far outside this band
    _kmkg = [(60.0 / (2 * math.pi * m.kv)) / math.sqrt(m.rm) / m.mass_kg
             for m in MOTORS.values()]
    check("motor constants are physically consistent (Km/kg 0.1-1.5)",
          all(0.1 < v < 1.5 for v in _kmkg),
          "%.3f..%.3f" % (min(_kmkg), max(_kmkg)))
    check("no catalogue section is too thin for a spar (7%)",
          min(a.t_over_c for a in _AF.values()) >= 0.07,
          "%.3f" % min(a.t_over_c for a in _AF.values()))

    # ------------------------------------------------------------------
    section("Lifting line, stability and thermal state")
    from dbf2027.liftingline import (downwash_gradient, lift_curve_slope,
                                     span_efficiency)

    # span efficiency: bounded, and the classical taper optimum near 0.35
    check("span efficiency never exceeds 1",
          all(span_efficiency(ar, t) <= 1.0
              for ar in (4, 8, 14) for t in (0.3, 0.6, 1.0)))
    es = [(t, span_efficiency(8.0, t))
          for t in (1.0, 0.7, 0.5, 0.4, 0.35, 0.3, 0.2)]
    best = max(es, key=lambda kv: kv[1])[0]
    check("induced drag is minimised near taper 0.35 (classical result)",
          0.30 <= best <= 0.45, "peak at taper %.2f" % best)
    check("a rectangular wing loses span efficiency as AR grows",
          span_efficiency(12.0, 1.0) < span_efficiency(5.0, 1.0))
    check("taper improves a rectangular wing",
          span_efficiency(8.0, 0.5) > span_efficiency(8.0, 1.0))

    # lift-curve slope and downwash
    check("lift slope rises with AR and stays under 2 pi",
          lift_curve_slope(4) < lift_curve_slope(12) < 2 * math.pi)
    check("downwash is between 0 and 1", 0.0 < downwash_gradient(6) < 1.0)
    check("downwash weakens as aspect ratio grows",
          downwash_gradient(12) < downwash_gradient(5))

    # static margin is now an OUTPUT
    check("static margin is computed and the baseline is stable",
          ac.aero_clean.static_margin > 0.0,
          "%.3f c" % ac.aero_clean.static_margin)
    check("neutral point sits aft of the wing AC",
          ac.aero_clean.x_np_over_c > 0.25)
    aft = Config()
    aft.build.target_static_margin = 0.06
    check("a smaller target margin puts the CG further aft",
          Aircraft(baseline(), aft).aero_clean.static_margin
          < ac.aero_clean.static_margin)
    unstable = Config()
    unstable.build.target_static_margin = -0.05
    unstable.build.static_margin_band = 0.0
    check("a CG behind the neutral point is rejected as unstable",
          not Aircraft(baseline(), unstable).static_checks()["ok"])

    # --- where the wing stalls first: what prices taper ----------------
    from dbf2027.liftingline import spanwise_loading, stall_factor
    _sl = spanwise_loading(6.0, 0.45)
    _pts = [(0.0, _sl[0][1] * _sl[0][2])] + [(e, r * c) for e, r, c in _sl]
    _int = sum(0.5 * (a[1] + b_[1]) * (b_[0] - a[0])
               for a, b_ in zip(_pts, _pts[1:]))
    check("spanwise loading integrates back to the wing CL",
          abs(_int - 1.0) < 0.03, "%.4f" % _int)
    check("a rectangular wing stalls at the root",
          stall_factor(6.0, 1.0, 0.18)[1] < 0.05)
    check("a sharply tapered wing stalls near the tip",
          stall_factor(6.0, 0.20, 0.18)[1] > 0.8,
          "%.2f" % stall_factor(6.0, 0.20, 0.18)[1])
    check("taper below 0.45 costs CLmax (tip stall)",
          stall_factor(6.0, 0.20, 0.18)[0] < stall_factor(6.0, 0.30, 0.18)[0]
          < stall_factor(6.0, 0.45, 0.18)[0])
    from dbf2027.optimizer import DesignSpace as _DSt, _fit_taper
    check("the wing is rectangular for now: taper is held at 1.0, not "
          "searched (m101)", _DSt().taper_bounds == (1.0, 1.0)
          and _fit_taper(baseline(taper=0.45), _DSt()).taper == 1.0)
    from dbf2027.numerics import nelder_mead as _nm101
    _x101, _f101 = _nm101(lambda v: (v[0] - 1.0) ** 2 + (v[1] - 2.0) ** 2
                          + (v[2] + 0.5) ** 2, [0.0, 7.0, 0.0],
                          [0.5, 0.5, 0.5], [-5.0, 7.0, -5.0], [5.0, 7.0, 5.0])
    check("Nelder-Mead holds a variable with equal bounds and solves the rest",
          _x101[1] == 7.0 and abs(_x101[0] - 1.0) < 1e-3
          and abs(_x101[2] + 0.5) < 1e-3, str(_x101))

    # --- the layout: every station calculated --------------------------
    lay = ac.layout
    check("fuselage = nose + payload bay + tailcone",
          abs(ac.body.fuse_length_m - (lay.nose_length_m
                                       + (lay.bay_end_m - lay.bay_start_m)
                                       + lay.tailcone_length_m)) < 1e-6)
    check("the layout's stations reproduce its CG",
          abs(sum(q.mass_kg * q.x_m for q in lay.stations)
              / sum(q.mass_kg for q in lay.stations) - lay.x_cg_m) < 1e-9)
    check("every part of the aeroplane is placed (masses add up)",
          abs(sum(q.mass_kg for q in lay.stations) - ac.mtow_m2_kg) < 1e-6,
          "%.4f vs %.4f kg" % (sum(q.mass_kg for q in lay.stations),
                               ac.mtow_m2_kg))
    check("the static margin sits in the target band",
          cfg.build.target_static_margin - 1e-6 <= lay.static_margin
          <= cfg.build.target_static_margin + cfg.build.static_margin_band
          + 1e-6, "%.3f" % lay.static_margin)
    check("the aero model flies the layout's CG",
          abs(ac.aero_clean.static_margin - lay.static_margin) < 1e-9)
    check("the payload bay is centred on the CG target",
          abs(0.5 * (lay.bay_start_m + lay.bay_end_m) - lay.x_cg_target_m)
          < 1e-12)
    check("the wing and tail sit where the layout put them",
          abs(lay.x_tail_ac_m - lay.x_wing_ac_m - ac.body.tail_arm_m) < 1e-9)
    from dbf2027.layout import tail_arm_optimum_m as _arm_opt
    _bd = ac.body
    _arm = _arm_opt(ac.wing, _bd.fuse_width_m, _bd.fuse_height_m, cfg.build)
    check("the tail arm is the lightest one the tailcone taper allows",
          abs(_bd.tail_arm_m - max(_arm, _bd.tail_arm_min_m)) < 1e-9)
    # an exactly balanced design (twin motors, battery trims it)
    _tw = Aircraft(baseline(n_motors=2, n_systems=2, prop="APC_11x7",
                            pack="TP_6S2200_25C", n_simulators=0), cfg)
    check("a balanced layout has the same CG empty and loaded",
          (not _tw.layout.balanced)
          or abs(_tw.layout.x_cg_empty_m - _tw.layout.x_cg_m) < 0.02
          * _tw.wing.chord_m, "%.4f vs %.4f m" % (_tw.layout.x_cg_empty_m,
                                                 _tw.layout.x_cg_m))
    _heavy_tail = Config()
    _heavy_tail.build.h_tail_areal_density = 3.6
    check("heavier tails move the lightest tail arm aft",
          Aircraft(baseline(), _heavy_tail).body.tail_arm_optimum_m
          > ac.body.tail_arm_optimum_m)
    check("a bigger tail raises the neutral point",
          ac.aero_clean.v_h > 0.0)

    # --- the hard rules limits have to actually bite -------------------
    heavy = Aircraft(baseline(n_simulators=60), cfg)
    check("an aeroplane over 55 lb is rejected",
          any("55 lb" in f for f in heavy.static_checks()["failures"]),
          "%.1f lb: %s" % (heavy.mtow_m2_kg * KG2LB,
                           heavy.static_checks()["failures"][:2]))
    check("five motors are rejected (the team cap is 4)",
          any("limit is 4" in f for f in
              Aircraft(baseline(n_motors=5, prop="APC_10x5"),
                       cfg).static_checks()["failures"]))
    from dbf2027.layout import motor_stations
    _n3 = motor_stations(3, 1.8, 0.2)
    _n4 = motor_stations(4, 1.8, 0.2)
    check("an odd motor count puts one motor on the nose",
          _n3[0] == 1 and len(_n3[1]) == 2 and motor_stations(1, 1.8, 0.2)
          == (1, [], 0.0))
    check("wing motors sit mid-share, symmetric, clear of the fuselage",
          _n4[0] == 0 and len(_n4[1]) == 4
          and abs(sum(_n4[1])) < 1e-12
          and abs(_n4[1][2] - (0.1 + 0.5 * _n4[2])) < 1e-12
          and abs(_n4[2] - 0.4) < 1e-12)
    check("more propulsion systems than motors is rejected",
          any("systems than motors" in f for f in
              Aircraft(baseline(n_motors=1, n_systems=2),
                       cfg).static_checks()["failures"]))
    two_packs = Aircraft(baseline(n_motors=2, prop="APC_11x7",
                                  pack="GA_6S4500_60C", n_systems=2), cfg)
    check("two 99.9 Wh packs break the 100 Wh total",
          any("over the 100 Wh limit" in f
              for f in two_packs.static_checks()["failures"]),
          "%.1f Wh" % two_packs.total_watt_hours)
    check("a legal two-system aeroplane carries both packs' energy",
          abs(Aircraft(baseline(n_motors=2, prop="APC_11x7",
                                pack="TP_6S2200_25C", n_systems=2),
                       cfg).total_watt_hours - 2 * 48.84) < 0.1)

    # the three coupled loops must actually converge together
    from dbf2027.loads import design_load_factors
    w = ac.mtow_m2_kg * 9.80665
    again = design_load_factors(ac.aero_clean, w,
                                ac._max_level_speed(w), cfg.build)
    check("load case is consistent with the final propulsion state",
          abs(again.n_ultimate - ac.n_ultimate) < 1e-3,
          "%.5f vs %.5f" % (again.n_ultimate, ac.n_ultimate))

    # windmilling drag and pilot efficiency.  A measured propeller's
    # windmill drag comes from its own table, so switch THAT off.
    from unittest import mock
    from dbf2027.components import Propeller as _Prop
    r_yes = mission_2(ac, cfg, dx=6.0)
    with mock.patch.object(_Prop, "windmill_drag_coeff",
                           property(lambda self: 0.0)):
        r_no = mission_2(Aircraft(baseline(), cfg), cfg, dx=6.0)
    check("a freewheeling propeller changes the lap",
          abs(r_yes.time_s - r_no.time_s) > 1e-3,
          "%.2f vs %.2f s" % (r_yes.time_s, r_no.time_s))
    from dbf2027.mission import fly_straight as _fs, turn_capability as _tcap
    _tc = _tcap(ac.aero_clean, ac.prop_sys, ac.w_m2_n, 0.7, cfg, ac.n_limit)
    _t_, _e_, _vx = _fs(ac.aero_clean, ac.prop_sys, ac.mtow_m2_kg, ac.w_m2_n,
                        0.7, 152.4, _tc.v_turn, _tc.v_turn, 5.0, dx=6.0)
    check("the glide into a turn, planned with the freewheeling drag, "
          "arrives at corner speed",
          abs(_vx - _tc.v_turn) < 0.05,
          "%.3f vs %.3f m/s" % (_vx, _tc.v_turn))
    calm = Config()
    calm.build.windmill_drag_coeff = 0.0
    check("a measured propeller's windmill drag ignores the fallback",
          abs(mission_2(Aircraft(baseline(), calm), calm, dx=6.0).time_s
              - r_yes.time_s) < 1e-9)
    perfect = Config()
    perfect.build.pilot_efficiency = 1.0
    r_perf = mission_2(Aircraft(baseline(), perfect), perfect, dx=6.0)
    check("a real pilot is slower than a perfect one",
          r_yes.time_s > r_perf.time_s,
          "%.2f vs %.2f s" % (r_yes.time_s, r_perf.time_s))
    check("pilot efficiency scales lap time as expected",
          abs(r_yes.lap_times[-1] * cfg.build.pilot_efficiency
              - r_perf.lap_times[-1]) < 0.35,
          "%.3f vs %.3f" % (r_yes.lap_times[-1] * cfg.build.pilot_efficiency,
                            r_perf.lap_times[-1]))

    # ------------------------------------------------------------------
    section("Payload, landing and course")
    from dbf2027.payload import best_bay_arrangement, size_payload
    small = size_payload(0.4, cfg.rules)
    big = size_payload(3.0, cfg.rules)
    # The container is ONE fixed assumption, deliberately not resized per
    # sensor.  What must still vary is the sensor envelope inside it.
    check("the container is the same box for every sensor",
          big.box_l == small.box_l and big.box_volume_m3 == small.box_volume_m3)
    check("container tare is the fixed assumed value",
          abs(small.container_tare_kg - cfg.rules.container_tare_kg) < 1e-12)
    check("a heavier sensor still has a bigger envelope",
          big.sensor_volume_m3 > small.sensor_volume_m3)
    check("a small sensor fits the box",
          small.fits(cfg.rules.container_fill_fraction))
    check("a heavy sensor simply needs to be denser",
          size_payload(5.0, cfg.rules).required_density(
              5.0, cfg.rules.container_fill_fraction)
          > size_payload(1.0, cfg.rules).required_density(
              1.0, cfg.rules.container_fill_fraction))
    # No simulators here on purpose: every simulator now weighs a full
    # loaded container (rule 3.1.3), so a heavy sensor multiplies straight
    # into gross weight and the 55 lb limit would bind before density does.
    # This check is about the DENSITY ceiling, so it isolates it.
    check("there is no sensor mass ceiling short of impossible density",
          Aircraft(baseline(sensor_mass_kg=8.0, n_simulators=0),
                   cfg).static_checks()["ok"],
          "; ".join(Aircraft(baseline(sensor_mass_kg=8.0, n_simulators=0),
                             cfg).static_checks()["failures"])[:60])
    check("a sensor denser than lead is rejected",
          any("denser than lead" in f for f in
              Aircraft(baseline(sensor_mass_kg=60.0, n_simulators=0),
                       cfg).static_checks()["failures"]))
    # And the coupling itself is worth asserting: simulators are no longer
    # a cheap way to add scored weight.
    _one = Aircraft(baseline(sensor_mass_kg=1.0, n_simulators=1), cfg)
    check("a simulator weighs exactly one loaded container",
          abs(_one.sim_mass_kg - (1.0 + _one.payload_geom.container_tare_kg))
          < 1e-12, "%.4f kg" % _one.sim_mass_kg)
    check("Mission 2 scored weight is (1 + n) loaded containers",
          abs(_one.m2_scored_kg - 2.0 * _one.loaded_container_kg) < 1e-12)
    _two = Aircraft(baseline(sensor_mass_kg=2.0, n_simulators=1), cfg)
    # Every extra kilogram of sensor is carried (1 + n) times over, once
    # in the real container and once in each simulator that has to match
    # it.  The container tare is a constant and does not scale, so the
    # relationship is linear with slope (1 + n), not proportional.
    check("each extra kg of sensor is carried (1 + n) times",
          abs((_two.m2_scored_kg - _one.m2_scored_kg) - 2.0 * 1.0) < 1e-9,
          "%.4f kg for 1 kg more sensor" % (_two.m2_scored_kg
                                            - _one.m2_scored_kg))
    _three = Aircraft(baseline(sensor_mass_kg=2.0, n_simulators=3), cfg)
    check("and with three simulators, four times over",
          abs((_three.m2_scored_kg - Aircraft(
              baseline(sensor_mass_kg=1.0, n_simulators=3),
              cfg).m2_scored_kg) - 4.0 * 1.0) < 1e-9)
    check("drop padding is sized from the drop height",
          0.010 < small.pad_m < 0.100, "%.3f m" % small.pad_m)

    # --- derived load factors -----------------------------------------
    lc = ac.load_case
    check("limit load is at least the sustained turn it must fly",
          ac.n_limit >= 1.0 / math.cos(math.radians(cfg.build.max_bank_deg))
          - 1e-9 or lc.driver == "wing CLmax ceiling",
          "%.2f g, %s" % (ac.n_limit, lc.driver))
    check("ultimate is the safety factor times limit",
          abs(ac.n_ultimate - cfg.build.load_safety_factor * ac.n_limit) < 1e-9)
    check("limit load never exceeds what the wing can generate",
          ac.n_limit <= lc.n_aero_ceiling + 1e-9,
          "%.2f vs ceiling %.2f" % (ac.n_limit, lc.n_aero_ceiling))
    check("a lighter wing loading is more gust-critical",
          Aircraft(baseline(n_simulators=0), cfg).load_case.n_gust_cruise
          > Aircraft(baseline(n_simulators=8), cfg).load_case.n_gust_cruise)
    check("the gust case is milder at the dive speed than at cruise",
          lc.n_gust_dive < lc.n_gust_cruise)
    check("a derived load factor lands in a sane band",
          2.5 <= ac.n_ultimate <= 9.0, "%.2f" % ac.n_ultimate)
    ovr = Config()
    ovr.build.ultimate_load_factor_override = 4.5
    check("an override bypasses the V-n derivation",
          abs(Aircraft(baseline(), ovr).n_ultimate - 4.5) < 1e-9)
    # The spar is often GAUGE limited on this class of aeroplane: the
    # strength requirement falls below the minimum wall you can build, so
    # the load factor stops mattering.  Both behaviours must be right.
    check("spar mass never decreases with load factor",
          _spar_at(cfg, 7.0) >= _spar_at(cfg, 3.0) - 1e-12,
          "%.4f vs %.4f kg" % (_spar_at(cfg, 7.0), _spar_at(cfg, 3.0)))
    check("a load factor past minimum gauge does grow the spar",
          _spar_at(cfg, 40.0) > _spar_at(cfg, 3.0),
          "%.4f vs %.4f kg" % (_spar_at(cfg, 40.0), _spar_at(cfg, 3.0)))
    check("gauge-limited flag agrees with the required wall",
          ac.spar_gauge_limited == (ac.spar_t_req_m < cfg.build.spar_min_wall_m))
    check("the wall actually used is never below minimum gauge",
          ac.spar_t_used_m >= cfg.build.spar_min_wall_m - 1e-15)

    # --- factor of safety on every wing load case ---------------------
    from dbf2027.weights import (_elliptic_moment, _elliptic_shear,
                                 spar_design)
    check("elliptic root moment is 4/(3 pi) of (W/2)(b/2)",
          abs(_elliptic_moment(0.0) * 4.0 / math.pi - 4.0 / (3.0 * math.pi))
          < 1e-12)
    check("elliptic moment and shear vanish at the tip",
          abs(_elliptic_moment(1.0)) < 1e-12 and abs(_elliptic_shear(1.0))
          < 1e-12 and abs(_elliptic_shear(0.0) - 1.0) < 1e-12)
    check("every wing load case meets the factor of safety",
          ac.structure_ok and all(v >= cfg.build.load_safety_factor - 1e-9
                                  for v in ac.spar.fos.values()),
          ", ".join("%s %.2f" % kv for kv in ac.spar.fos.items()))
    # With no minimum gauge the load alone sizes the wall, so the factor
    # of safety must come out at EXACTLY 1.5 where each case governs:
    # flight bending at the root, the tip-lift test outboard.
    _thin = Config()
    _thin.build.spar_min_wall_m = 1e-9
    # (torsion left out here: carried by the tube alone, at 50 MPa, it
    # needs more wall outboard than either bending case)
    _sd = spar_design(ac.wing, ac.w_m2_n, _thin.build, ac.n_ultimate, 0.0)
    # (Every station where flight bending governs sits at exactly 1.5,
    # so which of them reports the minimum is rounding: what matters is
    # that it is the inboard half - outboard, the tip test governs.)
    check("load-sized spar: flight bending governs inboard at 1.5",
          abs(_sd.fos["flight bending"] - 1.5) < 1e-9
          and _sd.fos_at_eta["flight bending"] < 0.5,
          "%.4f at %.2f" % (_sd.fos["flight bending"],
                            _sd.fos_at_eta["flight bending"]))
    check("load-sized spar: the tip-lift test governs outboard at 1.5",
          abs(_sd.fos["wing-tip test"] - 1.5) < 1e-9
          and _sd.fos_at_eta["wing-tip test"] > 0.2,
          "%.4f at %.2f" % (_sd.fos["wing-tip test"],
                            _sd.fos_at_eta["wing-tip test"]))
    _fos2 = Config()
    _fos2.build.spar_min_wall_m = 1e-9
    _fos2.build.load_safety_factor = 2.0
    _sd2 = spar_design(ac.wing, ac.w_m2_n, _fos2.build, ac.n_ultimate,
                       0.0)
    check("a larger factor of safety makes a heavier load-sized spar",
          _sd2.mass_kg > _sd.mass_kg and min(_sd2.fos.values()) > 2.0 - 1e-9,
          "%.1f vs %.1f g" % (1e3 * _sd2.mass_kg, 1e3 * _sd.mass_kg))
    check("torsion adds to the spar",
          spar_design(ac.wing, ac.w_m2_n, _thin.build, ac.n_ultimate,
                      ac.q_dive_pa).mass_kg > _sd.mass_kg)
    _weak = Config()
    _weak.build.spar_sigma_allow_pa = 2e6
    check("a spar that cannot carry the loads is rejected",
          any("factor of safety" in f for f in
              Aircraft(baseline(), _weak).static_checks()["failures"]))
    rows, cols, layers, bl, bw, bh = best_bay_arrangement(
        9, small.box_l, small.box_w, small.box_h)
    check("bay packing holds every box", rows * cols * layers >= 9,
          "%dx%dx%d" % (rows, cols, layers))

    more = Aircraft(baseline(n_simulators=8), cfg)
    fewer = Aircraft(baseline(n_simulators=2), cfg)
    check("more simulators means more wetted area and more drag",
          more.aero_clean.cd0(25.0) * more.wing_area_m2
          > fewer.aero_clean.cd0(25.0) * fewer.wing_area_m2)

    lr = mission_2(ac, cfg, dx=6.0)
    # A taller pattern altitude means a longer climb at best-rate speed,
    # which is slower than lap speed.  It must therefore never make a
    # SUCCESSFUL mission finish sooner.  It used to: the climb-out ground
    # credit divided straight distance by whole-lap time, understating the
    # average ground speed and over-crediting the first lap, and an
    # arbitrary 0.35*lap floor let a tall climb cut it by 65 percent.
    from dbf2027.mission import lap_ground_distance as _lgd
    _times = []
    for _h in (5.0, 20.0, 35.0, 45.0):
        _c = Config()
        _c.course.pattern_altitude_m = _h
        _r = mission_2(Aircraft(baseline(), _c), _c, dx=6.0)
        if _r.ok:
            _times.append((_h, _r.time_s))
    check("a taller pattern altitude never speeds a mission up",
          len(_times) >= 2 and all(b[1] >= a[1] - 1e-9
                                   for a, b in zip(_times, _times[1:])),
          ", ".join("%.0fm:%.2fs" % t for t in _times))
    check("the climb credit cannot exceed the climb's own duration",
          len(_times) >= 2
          and (_times[-1][1] - _times[0][1])
          <= (_times[-1][0] - _times[0][0]) / 0.5,
          "%.2f s gained over %.0f m"
          % (_times[-1][1] - _times[0][1], _times[-1][0] - _times[0][0]))
    _c = Config()
    _turn = turn_capability(ac.aero_clean, ac.prop_sys, ac.w_m2_n, 0.7, _c,
                            ac.n_limit)
    check("lap ground distance exceeds the straights alone",
          _lgd(_c, _turn) > _c.course.straight_length_m,
          "%.0f vs %.0f m" % (_lgd(_c, _turn), _c.course.straight_length_m))

    check("landing is simulated and takes time", lr.landing_time_s > 0.0,
          "%.1f s" % lr.landing_time_s)
    check("landing is excluded from the scored window",
          abs(lr.total_flight_time_s - (lr.time_s + lr.landing_time_s)) < 1e-6)
    check("landing rollout is a sane distance", 3.0 < lr.landing_roll_m < 300.0,
          "%.1f m" % lr.landing_roll_m)
    # --- Ground Mission (rule 3.3.4) ----------------------------------
    from dbf2027.scoring import ground_mission
    from dbf2027.payload import crush_stroke_m, drop_impact_speed_mps
    gm_ok, gm_raw, gm_why = ground_mission(ac, cfg)
    check("the ground mission runs and scores", gm_ok, gm_why)
    check("GM raw is sensor weight times drop height",
          abs(gm_raw - ac.sensor_max_kg * KG2LB
              * ac.design.drop_height_in) < 1e-9)
    check("a 60 inch drop is 5 ft", abs(cfg.rules.max_drop_height_in
                                        - 60.0) < 1e-12)
    check("impact speed follows v = sqrt(2 g h)",
          abs(drop_impact_speed_mps(60.0)
              - math.sqrt(2 * 9.80665 * 60.0 * 0.0254)) < 1e-12)
    check("a taller drop needs more crush stroke",
          crush_stroke_m(60.0) > crush_stroke_m(30.0) > crush_stroke_m(12.0))
    check("crush stroke scales with height, not its square root",
          abs(crush_stroke_m(60.0) / crush_stroke_m(30.0) - 2.0) < 1e-9)
    _lowdrop = Aircraft(baseline(drop_height_in=12.0), cfg)
    check("a lower declared drop needs less foam",
          _lowdrop.payload_geom.pad_m < ac.payload_geom.pad_m,
          "%.1f vs %.1f mm" % (1000 * _lowdrop.payload_geom.pad_m,
                               1000 * ac.payload_geom.pad_m))
    check("a lower declared drop scores less",
          ground_mission(_lowdrop, cfg)[1] < gm_raw)
    check("GM is part of the total, so the ceiling is 7.5",
          abs(_score_from_ratios(1.0, 1.0, 1.0, "field", 1.0)[3] - 7.5)
          < 1e-12)
    check("a mission scores zero when it is not completed",
          _score_from_ratios(1.0, 0.0, 0.0, "field", 0.0)[3] == 1.0)

    # --- Mission 3 may fly a lighter sensor (rule 3.1.1) ---------------
    _light = Aircraft(baseline(sensor_m3_kg=0.5), cfg)
    check("Mission 3 can fly under the declared maximum",
          abs(_light.sensor_m3_kg - 0.5) < 1e-12
          and abs(_light.sensor_max_kg - 0.9) < 1e-12)
    check("Mission 2 still flies the declared maximum",
          abs(_light.loaded_container_kg
              - (0.9 + _light.payload_geom.container_tare_kg)) < 1e-12)
    check("the ground mission still uses the declared maximum",
          abs(ground_mission(_light, cfg)[1]
              - 0.9 * KG2LB * _light.design.drop_height_in) < 1e-9)
    check("a heavier M3 sensor than declared is rejected",
          any("heavier than the declared" in f for f in
              Aircraft(baseline(sensor_m3_kg=2.0),
                       cfg).static_checks()["failures"]))
    check("a lighter M3 sensor flies at lower weight",
          _light.mtow_m3_kg < ac.mtow_m3_kg)

    # --- towed sensor (rules 3.1.1, 3.1.2) ----------------------------
    check("the tow line is at least 1.5 x wingspan",
          ac.tow_line_length_m >= 1.5 * ac.wing.span_m - 1e-9,
          "%.2f m for a %.2f m span" % (ac.tow_line_length_m, ac.wing.span_m))
    check("the tow line is charged as drag",
          ac.deployed_drag_area_m2
          > cfg.rules.sensor_deployed_drag_area_m2 + 1e-9)
    _wide = Aircraft(baseline(span_m=1.80, chord_m=0.30), cfg)
    check("a wider span tows a longer line and more drag",
          _wide.deployed_drag_area_m2 > ac.deployed_drag_area_m2)
    check("deploying the sensor increases drag at a given speed",
          ac.aero_deployed.drag(22.0, ac.w_m3_n)
          > ac.aero_clean.drag(22.0, ac.w_m3_n))
    check("a sensor under six inches long is rejected",
          any("minimum is 6 in" in f for f in
              Aircraft(baseline(sensor_mass_kg=0.30, n_simulators=0),
                       cfg).static_checks()["failures"]))

    # --- state of charge is measured against EVERY pack ---------------
    # With two propulsion systems the packs discharge together.  Dividing
    # by one pack made the voltage sag twice as fast as it should, cost
    # thrust for the whole second half of a mission, and reported a
    # NEGATIVE energy margin for a flight that finished with energy left.
    _two = Aircraft(baseline(n_motors=2, n_systems=2, prop="APC_11x7",
                             pack="TP_6S2200_25C"), cfg)
    _r2 = mission_2(_two, cfg, dx=6.0)
    check("a completed mission never reports a negative energy margin",
          (not _r2.ok) or _r2.energy_margin >= -1e-9,
          "%.3f" % _r2.energy_margin)
    check("energy used never exceeds what every pack holds",
          _r2.energy_j <= _two.total_watt_hours * 3600.0 + 1e-6,
          "%.1f Wh of %.1f" % (_r2.energy_j / 3600.0,
                               _two.total_watt_hours))

    # --- Mission 2: energy margin and a continuous cruise speed ---------
    _m2ac = Aircraft(baseline(), cfg)
    _m2m = mission_2(_m2ac, cfg, dx=6.0)
    check("Mission 2 finishes inside its planning margin",
          _m2m.ok and _m2m.energy_j <= (_m2ac.depth_of_discharge
                                        - cfg.build.m2_energy_margin)
          * _m2ac.total_watt_hours * 3600.0 + 1e-6,
          "%.1f Wh of %.1f" % (_m2m.energy_j / 3600.0, _m2ac.total_watt_hours))
    # Heavier and heavier: once full throttle runs out of energy the cruise
    # speed eases off continuously, so the time grows smoothly instead of
    # jumping a whole probe step (which used to leave designs on a cliff).
    # (an energy-limited heavy lifter: the v1.16 thorough-search winner,
    # its sweep re-centred on 11.2 kg when the 88% landing point came in)
    _el = Design(span_m=1.8288, chord_m=0.32735, airfoil="DAE21",
                 motor="HACKER_A50_14L_TURNADO_V3_400", prop="APC_15x10",
                 n_motors=1, pack="TP_7S3850_70C", n_systems=1,
                 sensor_mass_kg=12.0, sensor_m3_kg=12.0 * 0.857,
                 n_simulators=0, taper=0.4894, drop_height_in=60.0,
                 tail_vol_h=0.55)        # its tail then (m112)
    _plan = copy.deepcopy(cfg)          # the planned-cruise machinery
    _plan.build.full_throttle_flights = False
    _pts = []
    for _i in range(7):
        _s = 11.2 + 0.2 * _i
        _r = mission_2(Aircraft(replace(_el, sensor_mass_kg=_s,
                                        sensor_m3_kg=_s * 0.857), _plan),
                       _plan, dx=6.0)
        if _r.ok:
            _pts.append((_s, _r.time_s))
    _jumps = [b[1] - a[1] for a, b in zip(_pts, _pts[1:])]
    check("with planned cruise, a heavier Mission 2 slows down smoothly, "
          "never jumps",
          len(_pts) >= 6 and 0.0 <= min(_jumps) and max(_jumps) < 8.0,
          "%d flyable, steps %s s" % (len(_pts),
                                      ", ".join("%.1f" % j for j in _jumps)))

    # --- every flight at full throttle (the team's call, 2026-10-02) ---
    from dbf2027.mission import (fly_mission as _fm, fly_straight as _fs2,
                                 simulate_lap as _sl2,
                                 turn_capability as _tc2)
    check("every flight is planned at full throttle by default",
          Config().build.full_throttle_flights is True)
    _fa = Aircraft(baseline(), cfg)
    _fw = _fa.w_m2_n
    _sa = (_fa.aero_clean, _fa.prop_sys, _fw / 9.80665, _fw, 0.8, 300.0,
           34.0, 24.0, 0.0, 6.0)
    _g0, _g1 = _fs2(*_sa), _fs2(*_sa, full_power=True)
    check("the glide into a turn is flown as a glide but charged at full "
          "throttle",
          abs(_g1[0] - _g0[0]) < 1e-9 and abs(_g1[2] - _g0[2]) < 1e-9
          and _g1[1] > _g0[1],
          "%.1f vs %.1f kJ over %.1f s" % (_g1[1] / 1e3, _g0[1] / 1e3,
                                           _g1[0]))
    _ft = _tc2(_fa.aero_clean, _fa.prop_sys, _fw, 0.8, cfg, _fa.n_limit)
    _l1 = _sl2(_fa.aero_clean, _fa.prop_sys, _fw / 9.80665, _fw, 0.8, cfg,
               _ft, _ft.v_turn, 6.0)
    _l0 = _sl2(_fa.aero_clean, _fa.prop_sys, _fw / 9.80665, _fw, 0.8,
               _plan, _ft, _ft.v_turn, 6.0)
    check("a full-throttle lap takes the same time and more energy than "
          "the throttled one",
          abs(_l1[0] - _l0[0]) < 1e-9 and _l1[1] > _l0[1],
          "%.1f vs %.1f Wh in %.1f s" % (_l1[1] / 3600, _l0[1] / 3600,
                                         _l1[0]))
    # An energy-limited design: with planned cruise it finishes by easing
    # off; at full throttle there is nothing to ease off to.
    _hv = Aircraft(replace(_el, sensor_mass_kg=12.0,
                           sensor_m3_kg=12.0 * 0.857), cfg)
    _hp = mission_2(Aircraft(_hv.design, _plan), _plan, dx=6.0)
    _hf = mission_2(_hv, cfg, dx=6.0)
    _hx = _fm(_hv, _hv.aero_clean, _hv.w_m2_n, cfg, required_laps=5,
              dx=6.0)
    check("Mission 2 never eases off: it is the full-throttle flight or "
          "nothing",
          _hf.ok == _hx.ok and abs(_hf.time_s - _hx.time_s) < 1e-9
          and abs(_hf.energy_j - _hx.energy_j) < 1e-6,
          "planned %s, full throttle %s" % (
              "%.0f s" % _hp.time_s if _hp.ok else "fails",
              "%.0f s" % _hf.time_s if _hf.ok else "fails: " + _hf.reason))
    from dbf2027.optimizer import _screen_one as _so
    _neu = ("NEU_6521_3P5_362", 1, "GA_6S4500_60C", 1)
    _sp = {pk: (_so((_neu[0], pk) + _neu[1:] + (cfg,)),
                _so((_neu[0], pk) + _neu[1:] + (_plan,)))
           for pk in ("APC_18x10", "APC_16x8")}
    check("the screen ranks a setup by how long full throttle lasts",
          _sp["APC_16x8"][0]["fom"] > _sp["APC_18x10"][0]["fom"]
          and _sp["APC_16x8"][1]["fom"] < _sp["APC_18x10"][1]["fom"],
          "18x10 lasts %.0f s, 16x8 %.0f s"
          % (_sp["APC_18x10"][0]["full_throttle_s"],
             _sp["APC_16x8"][0]["full_throttle_s"]))
    # --- every motor count gets its share of the search ----------------
    from dbf2027.optimizer import (_shortlist, _configs_by_count,
                                   SearchSettings as _SSs)
    _sv = ([{"n_motors": 1, "motor": "M%d" % i, "fom": 100.0 - i}
            for i in range(40)]
           + [{"n_motors": 6, "motor": "S%d" % i, "fom": 10.0 - i}
              for i in range(10)])
    _sl = _shortlist(_sv, _SSs(keep_propulsion=20, per_motor_keep=3), [1, 6])
    check("the shortlist gives every motor count its share",
          sum(r["n_motors"] == 6 for r in _sl) == 10
          and sum(r["n_motors"] == 1 for r in _sl) == 10,
          "%d one-motor, %d six-motor" % (
              sum(r["n_motors"] == 1 for r in _sl),
              sum(r["n_motors"] == 6 for r in _sl)))
    _rws = ([{"design": baseline(motor=m), "objective": 9.0 - 0.1 * i}
             for i, m in enumerate(list(MOTORS)[:6])]
            + [{"design": baseline(n_motors=2, motor=m),
                "objective": 5.0 - 0.1 * i}
               for i, m in enumerate(list(MOTORS)[:3])])
    _pk = _configs_by_count(_rws, 4, [1, 2])
    check("refinement and the finals keep every motor count's best",
          sum(r["design"].n_motors == 2 for r in _pk) == 2
          and _pk[0]["objective"] >= _pk[-1]["objective"],
          "%s" % [r["design"].n_motors for r in _pk])
    # --- wing-propeller layouts, swirl, roll trim (2026-10-02) ---------
    from dbf2027.layout import motor_stations as _mst
    _nz = Aircraft(baseline(), cfg)       # a nose tractor
    # RULE 3.1 limits the wing's span; read literally, a wingtip prop's
    # axis may sit on the tip.  Counting the discs, their edge stops there.
    _n, _ys, _sh = _mst(2, 1.8288, 0.185, "tip", 0.4064)
    _n, _ys2, _sh = _mst(2, 1.8288, 0.185, "tip", 0.4064, True)
    check("a wingtip prop sits on the tip (its disc's edge there if the "
          "discs count in the span)",
          abs(max(_ys) - 0.9144) < 1e-12
          and abs(max(_ys2) - (0.9144 - 0.2032)) < 1e-12
          and Config().rules.wingspan_counts_props is False,
          "%.4f / %.4f m" % (max(_ys), max(_ys2)))
    _st18 = Aircraft(baseline(span_m=1.8288, n_motors=4, prop="APC_24x12",
                              wing_prop_layout="stagger"), cfg)
    _st20 = Aircraft(baseline(span_m=1.8288, n_motors=4, prop="APC_27x13",
                              wing_prop_layout="stagger"), cfg)
    check("staggered props may overlap up to 30% - and pay for it",
          _st18.stagger_thrust_factor < 1.0
          and not any("staggered propellers overlap" in f
                      for f in _st18.static_checks()["failures"])
          and any("staggered propellers overlap" in f
                  for f in _st20.static_checks()["failures"]),
          "24x12 thrust x%.3f" % _st18.stagger_thrust_factor)
    _w2 = Aircraft(baseline(n_motors=2), cfg)
    _a2 = _w2.aero_clean
    _t2 = _w2.prop_sys.thrust(27.0, 0.7, 1.0)
    _cl2 = (_a2.cl_required(27.0, _w2.w_m2_n, _t2)
            * _a2.slipstream_gains(27.0, _t2)[0])
    check("wing props turning inboard-up give some induced drag back",
          0.8 < _a2.induced_factor(27.0, _t2, _cl2) < 1.0
          and _nz.aero_clean.induced_factor(
              27.0, _nz.prop_sys.thrust(27.0, 0.7, 1.0), 0.4) == 1.0,
          "x%.3f at %.1f deg swirl" % (_a2.induced_factor(27.0, _t2, _cl2),
                                      math.degrees(_a2.swirl_angle(27.0,
                                                                   _t2))))
    _tn = _nz.prop_sys.thrust(20.0, 0.9, 1.0)
    check("an unpaired nose motor costs roll-trim drag; a counter-rotating "
          "pair none",
          _nz.aero_clean.roll_trim_drag(20.0, _tn) > 0.0
          and _a2.roll_trim_drag(20.0, _t2) == 0.0)
    _big = Aircraft(baseline(n_motors=2, motor="KDE_6815XF_205"), cfg)
    _sml = Aircraft(baseline(n_motors=2, motor="KDE_2306XF_2050"), cfg)
    check("a nacelle's drag follows its motor's size",
          _big.aero_clean.nacelle_drag_area_m2
          > 4.0 * _sml.aero_clean.nacelle_drag_area_m2 > 0.0,
          "%.5f vs %.5f m2" % (_big.aero_clean.nacelle_drag_area_m2,
                               _sml.aero_clean.nacelle_drag_area_m2))
    from dbf2027.optimizer import _layouts_for, _swap_neighbours as _swn
    _nb2 = _swn(baseline(n_motors=2), _DSs(), cfg, _SSs(), {})
    check("the search tries every layout a motor count allows",
          _layouts_for(1) == ["mid"] and _layouts_for(2) == ["mid", "tip"]
          and _layouts_for(4) == ["mid", "tip", "stagger"]
          and any(n.wing_prop_layout == "tip" for n in _nb2)
          and baseline(n_motors=2).discrete_key()
          != baseline(n_motors=2, wing_prop_layout="tip").discrete_key())
    _m3f = mission_3(_fa, cfg, dx=6.0)
    _m3x = _fm(_fa, _fa.aero_deployed, _fa.w_m3_n, cfg, required_laps=None,
               dx=6.0)
    check("Mission 3 flies full throttle too, for as many laps as that gives",
          _m3f.laps == _m3x.laps and abs(_m3f.time_s - _m3x.time_s) < 1e-9,
          "%d laps" % _m3f.laps)

    # --- the swap test --------------------------------------------------
    from dbf2027.optimizer import (DesignSpace as _DSw, SearchSettings as _SSw,
                                   _swap_neighbours, _final_one)
    _base = baseline()
    _nbs = _swap_neighbours(_base, _DSw(), cfg, _SSw(), {})
    check("swaps try every other propeller and airfoil",
          {n.prop for n in _nbs} >= set(PROPS) - {_base.prop}
          and len({n.airfoil for n in _nbs}) == len(AIRFOILS))
    check("swaps try other motors and packs, each within its voltage rating",
          len({n.motor for n in _nbs}) > 5 and len({n.pack for n in _nbs}) > 5
          and all(PACKS[n.pack].series <= MOTORS[n.motor].max_cells
                  for n in _nbs))
    check("every swap changes exactly one component (plus its propeller)",
          all(sum((n.airfoil != _base.airfoil, n.motor != _base.motor,
                   (n.pack, n.n_systems) != (_base.pack, _base.n_systems)))
              <= 1 for n in _nbs))
    _fr = _final_one((_base, cfg, 12.0, None, None, None))
    check("final scoring re-checks a finalist at a 4x finer step",
          (not _fr.get("feasible")) or "hardened_steps" in _fr)
    if _fr.get("feasible"):
        from dbf2027.optimizer import evaluate_design as _ev103
        _fl = [_ev103(_fr["design"], cfg, 12.0 / k_).get("m3_laps") or 0
               for k_ in (4.0, 8.0)]
        check("a finalist flies at least the Mission 3 laps it is scored for at "
              "4x and 8x finer steps (m103)",
              min(_fl) >= (_fr.get("m3_laps") or 0),
              "%s laps scored, %s finer" % (_fr.get("m3_laps"), _fl))
    from dbf2027.optimizer import _prefetch_best_props
    from dbf2027.components import searchable as _srch
    _sp100 = _DSw()
    _sp100.motors = _srch(MOTORS)[:4]
    _sp100.props = sorted(PROPS)[:6]
    _sp100.packs = _srch(PACKS)[:3] + [_base.pack]
    _st100 = _SSw(processes=1)
    _plain100 = _swap_neighbours(_base, _sp100, cfg, _st100, {})
    _c100: dict = {}
    _n100 = _prefetch_best_props([_base], _sp100, cfg, _st100, _c100)
    _k100 = set(_c100)
    _pre100 = _swap_neighbours(_base, _sp100, cfg, _st100, _c100)
    check("swap test: the propeller ratings flown ahead on the pool give the "
          "same neighbours, and leave nothing to fly one at a time (m100)",
          _pre100 == _plain100 and set(_c100) == _k100 and _n100 > 0,
          "%d neighbours, %d flights" % (len(_pre100), _n100))

    # --- motor heating -------------------------------------------------
    from dbf2027.propulsion import motor_thermal, winding_step
    _r, _c = motor_thermal(MOTORS["SCORPION_SII_4025_440"], cfg.build)
    check("winding steady state is ambient + P R",
          abs(winding_step(25.0, 25.0, 100.0, 1e6, _r, _c) - (25.0 + 100.0 * _r))
          < 1e-9)
    check("after one time constant the winding is 63% of the way there",
          abs((winding_step(25.0, 25.0, 100.0, _r * _c, _r, _c) - 25.0)
              / (100.0 * _r) - (1.0 - math.exp(-1.0))) < 1e-9)
    check("a bigger motor sheds heat better (lower R)",
          motor_thermal(MOTORS["NEU_6530_3_295"], cfg.build)[0]
          < motor_thermal(MOTORS["SUNNYSKY_X3120_V3_585"], cfg.build)[0])
    _hot = Aircraft(baseline(), cfg).prop_sys_at(100.0)
    _cold = Aircraft(baseline(), cfg).prop_sys_at(20.0)
    check("hot windings cost thrust (copper resistance rises)",
          _hot.solve(0.0, 1.0, 1.0).thrust_total_n
          < _cold.solve(0.0, 1.0, 1.0).thrust_total_n
          and abs(_hot.rm / _cold.rm - (1.0 + cfg.build.copper_temp_coeff
                                        * 80.0)) < 1e-9)
    _m2h = mission_2(Aircraft(baseline(), cfg), cfg, dx=6.0)
    check("a mission reports its peak winding temperature",
          cfg.env.temp_c < _m2h.motor_temp_c < cfg.build.motor_temp_limit_c,
          "%.1f C" % _m2h.motor_temp_c)
    _tight = Config()
    _tight.build.motor_temp_limit_c = _tight.env.temp_c + 3.0
    _m2t = mission_2(Aircraft(baseline(), _tight), _tight, dx=6.0)
    check("a motor that would overheat fails the mission",
          not _m2t.ok and "motor windings" in _m2t.reason, _m2t.reason)
    check("the site presets carry their air temperature",
          Environment.at_site("tucson-hot").temp_c == 36.0
          and Environment.at_site("tucson-morning").temp_c == 12.0)

    # --- safe depth of discharge and the ESC low-voltage cutoff --------
    # Every completed mission, landing included, stays inside the pack's
    # safe depth of discharge (80% for LiPo) - and a battery-limited
    # Mission 3 should come close to it, not stop 25% short.
    _sd = Aircraft(baseline(), cfg)
    _sd_cap = _sd.usable_energy_j()
    check("a mission may use 95% of the pack by default",
          abs(_sd_cap - 0.95 * _sd.total_watt_hours * 3600.0) < 1e-6)
    _dod = Config()
    _dod.build.depth_of_discharge = 0.99
    check("the chemistry's 95% ceiling caps the depth of discharge",
          abs(Aircraft(baseline(), _dod).depth_of_discharge - 0.95) < 1e-12)
    _m2s, _m3s = mission_2(_sd, cfg, dx=6.0), mission_3(_sd, cfg, dx=6.0)
    check("M2 lands inside the safe depth of discharge",
          (not _m2s.ok) or _m2s.energy_j <= _sd_cap + 1e-6,
          "%.1f of %.1f Wh" % (_m2s.energy_j / 3600.0, _sd_cap / 3600.0))
    check("M3 lands inside the safe depth of discharge",
          (not _m3s.ok) or _m3s.energy_j <= _sd_cap + 1e-6,
          "%.1f of %.1f Wh" % (_m3s.energy_j / 3600.0, _sd_cap / 3600.0))
    check("every mission says what ended it",
          _m2s.stopped_by == "laps"
          and _m3s.stopped_by in ("battery", "window", "cap"),
          "M2 %s, M3 %s" % (_m2s.stopped_by, _m3s.stopped_by))
    if _m3s.ok and _m3s.stopped_by == "battery":
        # the landing is budgeted up front, so a battery-limited flight
        # can come up short of the cap by at most one lap plus landing
        _lap_e = _m3s.energy_j / max(_m3s.laps + 1, 1)
        check("a battery-limited M3 uses most of the safe budget",
              _m3s.energy_j >= _sd_cap - 2.0 * _lap_e,
              "%.1f of %.1f Wh" % (_m3s.energy_j / 3600.0, _sd_cap / 3600.0))
    _res_cfg = Config()
    _res_cfg.build.energy_reserve = 0.25
    _m3r = mission_3(Aircraft(baseline(), _res_cfg), _res_cfg, dx=6.0)
    check("an energy reserve in the build standard is honoured",
          (not _m3r.ok) or _m3r.energy_j <= 0.75 * _sd_cap + 1e-6,
          "%.1f Wh with 25%% held back" % (_m3r.energy_j / 3600.0))
    # The landing point (team's call, 2026-10-01): no lap is started that
    # would take any pack past land_at_discharge; the landing comes out of
    # what is left, and the 2% Mission 2 buffer it replaced is gone.
    _pj = _sd.total_energy_wh * 3600.0
    check("every pack turns for landing at 88%, with no Mission 2 buffer",
          abs(cfg.build.land_at_discharge - 0.88) < 1e-12
          and cfg.build.m2_energy_margin == 0.0)
    check("the default missions turn for landing by the 88% point",
          all((not r.ok) or r.energy_j <= (0.88 + 0.02) * _pj
              for r in (_m2s, _m3s)),
          "M2 %.1f, M3 %.1f of %.1f Wh" % (_m2s.energy_j / 3600.0,
                                            _m3s.energy_j / 3600.0,
                                            _pj / 3600.0))
    _lp = Config()
    _lp.build.land_at_discharge = 0.50
    _lp_ac = Aircraft(baseline(), _lp)
    _lp2, _lp3 = mission_2(_lp_ac, _lp, dx=6.0), mission_3(_lp_ac, _lp, dx=6.0)
    check("an early landing point ends Mission 2 and says so",
          not _lp2.ok and "landing point" in _lp2.reason, _lp2.reason)
    check("an early landing point stops Mission 3's laps there",
          _lp3.ok and _lp3.laps < _m3s.laps
          and _lp3.energy_j <= (0.50 + 0.02) * _pj,
          "%d laps, %.1f Wh" % (_lp3.laps, _lp3.energy_j / 3600.0))
    _hp = Config()
    _hp.build.land_at_discharge = 1.0
    _hp3 = mission_3(Aircraft(baseline(), _hp), _hp, dx=6.0)
    check("with the landing point at 100% the depth of discharge binds",
          (not _hp3.ok) or _hp3.energy_j <= _sd_cap + 1e-6)
    # A pack with ruinous resistance: full throttle drags it under the
    # ESC's cutoff, so it cannot fly however much energy it holds.
    PACKS["_SELFTEST_HIGH_R"] = replace(PACKS["GA_6S3300_60C"],
                                        key="_SELFTEST_HIGH_R",
                                        r_cell_mohm=80.0)
    try:
        _hr = mission_2(Aircraft(baseline(pack="_SELFTEST_HIGH_R"), cfg),
                        cfg, dx=6.0)
    finally:
        del PACKS["_SELFTEST_HIGH_R"]
    check("a pack that sags below the ESC cutoff cannot fly",
          (not _hr.ok) and "low-voltage cutoff" in _hr.reason, _hr.reason)

    # --- rows scored on different references must not be compared -----
    from dbf2027.scoring import apply_references
    _rows = [{"feasible": True, "m1": 1.0, "m2_raw_lb_per_s": 0.060,
              "m3_raw_lap_lb": 40.0, "gm_raw_lb_in": 500.0},
             {"feasible": True, "m1": 1.0, "m2_raw_lb_per_s": 0.080,
              "m3_raw_lap_lb": 60.0, "gm_raw_lb_in": 900.0}]
    apply_references(_rows, 0.060, 22.0, 120.0, "self")
    lenient = [r["objective"] for r in _rows]
    apply_references(_rows, 0.080, 60.0, 900.0, "self")
    strict = [r["objective"] for r in _rows]
    check("generous references inflate every objective",
          all(a > b for a, b in zip(lenient, strict)),
          "%s vs %s" % ([round(v, 2) for v in lenient],
                        [round(v, 2) for v in strict]))
    check("but the ORDER is unchanged by the reference scale",
          (lenient[1] > lenient[0]) == (strict[1] > strict[0]))
    check("a design better on every raw score wins under any reference",
          strict[1] > strict[0] and lenient[1] > lenient[0])

    # --- a design may throttle back rather than simply fail -----------
    from dbf2027.mission import fly_straight as _fs
    _ac = Aircraft(baseline(), cfg)
    _w = _ac.w_m2_n
    _args = (_ac.aero_clean, _ac.prop_sys, _w / 9.80665, _w, 0.8, 600.0,
             20.0, 1e9, 0.0, 6.0)
    t_open, e_open, v_open = _fs(*_args)
    t_cap, e_cap, v_cap = _fs(*_args, v_cruise=22.0)
    check("a commanded cruise speed is actually held",
          v_cap <= 22.0 + 0.5, "%.1f m/s" % v_cap)
    check("flat out is faster than a capped cruise", t_open < t_cap,
          "%.1f vs %.1f s" % (t_open, t_cap))
    check("and a capped cruise uses less energy", e_cap < e_open,
          "%.1f vs %.1f kJ" % (e_cap / 1000, e_open / 1000))
    check("an infinite cap is exactly the old flat-out behaviour",
          abs(_fs(*_args, v_cruise=float("inf"))[0] - t_open) < 1e-12)

    # --- the contest site --------------------------------------------
    from dbf2027.config import (CONDITIONS, DEFAULT_CONDITIONS, SITE_NAME,
                                SITE_ELEVATION_M)
    _dflt = Environment()
    _site = Environment.at_site(DEFAULT_CONDITIONS)
    # The dataclass defaults are hand-copied from at_site because a field
    # default cannot call the function that derives it.  If they drift,
    # every result quietly moves.
    check("the default environment IS the default site",
          abs(_dflt.rho - _site.rho) < 1e-5
          and abs(_dflt.nu - _site.nu) < 1e-8
          and abs(_dflt.wind_mps - _site.wind_mps) < 1e-9
          and abs(_dflt.field_elev_m - _site.field_elev_m) < 1e-9,
          "rho %.6f vs %.6f, nu %.4e vs %.4e"
          % (_dflt.rho, _site.rho, _dflt.nu, _site.nu))
    check("the site is Tucson", "Tucson" in SITE_NAME, SITE_NAME)
    check("and the fly-off dates are named", "2027" in SITE_NAME)
    check("every named condition is at the site elevation, except ISA",
          all(abs(v[0] - SITE_ELEVATION_M) < 1e-9
              for k, v in CONDITIONS.items() if k != "isa"))
    check("the site is well above sea level",
          600.0 < SITE_ELEVATION_M < 900.0, "%.0f m" % SITE_ELEVATION_M)
    _m = Environment.at_site("tucson-morning")
    _h = Environment.at_site("tucson-hot")
    check("morning air is denser than the nominal afternoon",
          _m.rho > _site.rho > _h.rho,
          "%.4f > %.4f > %.4f" % (_m.rho, _site.rho, _h.rho))
    check("and the hot afternoon is the windiest",
          _h.wind_mps > _site.wind_mps > _m.wind_mps)
    check("Tucson costs at least 3000 ft of density altitude",
          _site.density_altitude_ft > 3000.0,
          "%.0f ft" % _site.density_altitude_ft)
    check("sea-level ISA is still available for comparison",
          abs(Environment.at_site("isa").rho - 1.225) < 1e-3)
    check("an unknown condition falls back rather than raising",
          abs(Environment.at_site("no-such-place").rho - _site.rho) < 1e-9)
    check("the thinner site air raises the stall speed",
          Aircraft(baseline(), cfg).aero_clean.v_stall(100.0)
          > Aircraft(baseline(), _isa_cfg()).aero_clean.v_stall(100.0))

    # --- air density is the quantity the model actually flies on -----
    # Environment() is the CONTEST SITE now, not sea level, so ask for
    # ISA explicitly.
    _isa = Environment.at_site("isa")
    check("sea-level ISA is zero density altitude",
          abs(_isa.density_altitude_ft) < 2.0,
          "%.1f ft" % _isa.density_altitude_ft)
    check("thinner air is a higher density altitude",
          Environment(rho=0.95).density_altitude_ft
          > Environment(rho=1.10).density_altitude_ft > 0.0)
    # At the ISA temperature for 5000 ft (15 - 6.5 C per km) the density
    # must match the standard atmosphere.  Note the function takes the
    # temperature explicitly, so 1524 m at a sea-level 15 C is NOT ISA -
    # it is a hot day at altitude, and thinner.
    _isa5k = Environment.at_density_altitude(1524.0, 15.0 - 6.5 * 1.524)
    check("density matches the standard atmosphere at the ISA temperature",
          abs(_isa5k.rho - 1.225 * (1.0 - 6.87535e-6 * 5000.0) ** 4.2559)
          < 0.01, "%.4f" % _isa5k.rho)
    check("and that round-trips back to 5000 ft",
          abs(_isa5k.density_altitude_ft - 5000.0) < 120.0,
          "%.0f ft" % _isa5k.density_altitude_ft)
    check("a hot day at altitude is thinner than ISA there",
          Environment.at_density_altitude(1524.0, 15.0).rho < _isa5k.rho)
    _hot = Environment.at_density_altitude(0.0, 40.0)
    check("hot air is thinner than cold air at the same elevation",
          _hot.rho < Environment.at_density_altitude(0.0, 0.0).rho)
    check("high ground is thinner than low ground at the same temperature",
          Environment.at_density_altitude(2000.0, 15.0).rho
          < Environment.at_density_altitude(0.0, 15.0).rho)
    _thin = Config()
    _thin.env = Environment(rho=0.95, nu=cfg.env.nu)
    _r_thin = mission_2(Aircraft(baseline(), _thin), _thin, dx=6.0)
    _r_isa = mission_2(Aircraft(baseline(), cfg), cfg, dx=6.0)
    check("thin air never makes a mission quicker",
          (not _r_thin.ok) or (not _r_isa.ok)
          or _r_thin.time_s >= _r_isa.time_s - 1e-9,
          "%.1f vs %.1f s" % (_r_thin.time_s, _r_isa.time_s))
    check("and it raises the stall speed",
          Aircraft(baseline(), _thin).aero_clean.v_stall(100.0)
          > Aircraft(baseline(), cfg).aero_clean.v_stall(100.0))

    check("course has 2000 ft of straights",
          abs(cfg.course.straight_length_m - 2000 * 0.3048) < 1e-6,
          "%.1f m" % cfg.course.straight_length_m)
    check("course has 720 degrees of turning",
          abs(cfg.course.total_turn_deg - 720.0) < 1e-9)
    check("wingspan limit is 6 ft",
          abs(cfg.rules.max_wingspan_m - 6.0 * 0.3048) < 1e-9,
          "%.4f m" % cfg.rules.max_wingspan_m)
    over = Aircraft(baseline(span_m=2.5, chord_m=0.30), cfg)
    check("an over-span wing is rejected",
          any("span" in f for f in over.static_checks()["failures"]),
          "%s" % over.static_checks()["failures"])
    from dbf2027.optimizer import DesignSpace as _DS
    check("the search never samples beyond the span limit",
          _DS().continuous_bounds(cfg)[0][1] <= cfg.rules.max_wingspan_m + 1e-12)

    _lipo_s = sorted({pk.series for pk in PACKS.values()
                      if pk.chemistry == "LiPo"})
    check("LiPo packs span 3S to 13S", _lipo_s[0] <= 3 and _lipo_s[-1] == 13,
          str(_lipo_s))
    _hv = Aircraft(baseline(pack="MX_12S2000_27C",
                            motor="SUNNYSKY_X4120_V3_480"), cfg)
    check("a 12S pack on a 6S-rated motor is rejected",
          any("rated for 6S" in f for f in _hv.static_checks()["failures"]))
    check("a 12S pack gets a 12S-rated ESC",
          _hv.prop_sys.esc.max_cells >= 12, _hv.prop_sys.esc.key)
    from dbf2027.optimizer import _screen_one
    check("the propulsion screen drops a pack over the motor's rating",
          _screen_one(("SUNNYSKY_X4120_V3_480", "APC_10x7", 1,
                       "MX_12S2000_27C", 1, cfg)) is None)
    check("every catalogue motor states a cell rating",
          all(3 <= m.max_cells <= 14 for m in MOTORS.values()))
    for chem in ("LiPo", "LiIon"):
        check("%s packs are available" % chem,
              any(pk.chemistry == chem for pk in PACKS.values()))
    check("only lithium packs are catalogued",
          all(pk.chemistry in ("LiPo", "LiIon") for pk in PACKS.values()))
    check("every catalogued pack is under the 100 Wh per-pack limit",
          all(pk.watt_hours <= 100.0 for pk in PACKS.values()),
          "worst %.1f Wh" % max(pk.watt_hours for pk in PACKS.values()))
    check("Li-ion carries more energy per kilogram than LiPo",
          min(pk.watt_hours / pk.mass_kg
              for pk in PACKS.values() if pk.chemistry == "LiIon")
          > max(pk.watt_hours / pk.mass_kg
                for pk in PACKS.values() if pk.chemistry == "LiPo"))
    check("LiPo delivers more current than Li-ion",
          min(pk.c_rating for pk in PACKS.values() if pk.chemistry == "LiPo")
          > max(pk.c_rating for pk in PACKS.values()
                if pk.chemistry == "LiIon"))
    check("watt-hours are rated capacity times rated voltage",
          all(abs(pk.watt_hours - pk.capacity_ah * pk.v_nominal) < 1e-9
              for pk in PACKS.values()))
    _iac = Aircraft(baseline(), cfg)
    _vtop = _iac._max_level_speed(_iac.w_m2_n)
    _ks = [k_ for _v, k_ in _iac.install_thrust_factor]
    check("installed thrust is below free-air thrust at every speed, and "
          "the propulsion system carries the same table",
          all(0.80 <= k_ < 1.0 for k_ in _ks)
          and 0.85 < _iac.install_thrust_factor_at(_vtop) < 1.0
          and _iac.prop_sys.install_thrust_factor
          == _iac.install_thrust_factor,
          "%.4f at %.1f m/s" % (_iac.install_thrust_factor_at(_vtop), _vtop))
    check("installed-thrust table runs from standstill to the propellers' "
          "speed ceiling",
          _iac.install_thrust_factor[0][0] == 0.0
          and abs(_iac.install_thrust_factor[-1][0]
                  - _iac.prop_sys.speed_ceiling(1.0)) < 1e-3)

    # The current limits are checked at the WORST airspeed, not standstill:
    # a measured propeller's CP peaks above J = 0.
    _rc = _iac.prop_sys.rating_check()
    _fine = max(_iac.prop_sys.solve(0.25 * k, 1.0, 1.0).total_current_a
                for k in range(0, 241))
    check("the rating check finds the peak current over every airspeed",
          _rc["peak_current_speed_mps"] > 0.0
          and _rc["total_current_a"] > _rc["static_current_a"]
          and abs(_rc["total_current_a"] - _fine) < 0.05,
          "static %.1f A, peak %.2f A at %.1f m/s, scan %.2f A"
          % (_rc["static_current_a"], _rc["total_current_a"],
             _rc["peak_current_speed_mps"], _fine))
    _cap = Config()
    _cap.rules.max_fuse_current_a = 0.5 * (_rc["static_current_a"]
                                           + _rc["total_current_a"])
    _capped = Aircraft(baseline(), _cap).static_checks()
    check("a design over its fuse only in flight is rejected",
          not _capped["ok"]
          and any("blows its" in f for f in _capped["failures"]),
          "; ".join(_capped["failures"]))
    # The team's fuse: listed C x capacity, capped at 100 A, per system.
    _fz = [PropulsionSystem(MOTORS["SCORPION_SII_4025_440"], PROPS["APC_12x8"],
                            1, _p, cfg.env) for _p in PACKS.values()]
    check("every system's fuse is the pack's C x Ah, capped at 100 A",
          all(abs(_s.fuse_rating_a
                  - min(_s.pack.capacity_ah * _s.pack.c_rating, 100.0))
              < 1e-9 for _s in _fz))
    # A flight starts with the windings at air temperature: on a morning
    # colder than the catalogue's 20 C they draw more than the catalogue
    # figure, and that is what the fuse check holds the design to.
    _cold = copy.deepcopy(cfg)
    _cold.env = Environment.at_site("tucson-morning")
    _cac = Aircraft(baseline(), _cold)
    _cold_pk = _cac.prop_sys_at(_cold.env.temp_c).rating_check()[
        "pack_current_a"]
    check("the fuse check takes the windings at the morning's air "
          "temperature",
          _cold.env.temp_c < 20.0
          and _cold_pk > _cac.prop_sys.rating_check()["pack_current_a"]
          and abs(_cac.fuse_peak_a() - _cold_pk) < 1e-9,
          "%.2f vs %.2f A" % (_cold_pk,
                              _cac.prop_sys.rating_check()["pack_current_a"]))
    _cap2 = copy.deepcopy(_cold)
    _cap2.rules.max_fuse_current_a = 0.5 * (
        _cac.prop_sys.rating_check()["pack_current_a"] + _cold_pk)
    check("a design inside its fuse only with warm windings is rejected "
          "on a cold morning",
          any("blows its" in f
              for f in Aircraft(baseline(), _cap2).static_checks()["failures"]))
    # RULE 3.2.3.1d - each system has its own fuse; nothing caps the total.
    from dataclasses import fields as _rf
    check("no rule caps the current summed over every system",
          not any("total_current" in f.name for f in _rf(Config().rules)))
    _two = Aircraft(baseline(n_motors=2, n_systems=2,
                             pack="CNHL_6S2200_70C"), cfg)
    _trc = _two.prop_sys.rating_check()
    check("two systems are each held to their own fuse, not a shared one",
          abs(_trc["total_current_a"] - 2.0 * _trc["pack_current_a"]) < 1e-9
          and not any("total" in f for f in _two.static_checks()["failures"]),
          "%.0f A per pack, %.0f A total" % (_trc["pack_current_a"],
                                            _trc["total_current_a"]))

    # RULE 3.2.3.2d - the fuse is min(pack rating, the rules' fuse cap), and
    # the pack rating is the label's own current when it states one.
    _fc = Config()
    _fc.rules.max_fuse_current_a = 0.9 * _rc["pack_current_a"]
    _fac = Aircraft(baseline(), _fc)
    check("the fuse cap is the rules setting, not a hard-coded 100 A",
          abs(_fac.prop_sys.fuse_rating_a - _fc.rules.max_fuse_current_a)
          < 1e-9 and not _fac.static_checks()["ok"])
    _pk = PACKS[baseline().pack]
    check("a pack's labelled discharge current replaces capacity x C",
          replace(_pk, i_label_a=42.0).i_max_a == 42.0
          and _pk.i_max_a == _pk.capacity_ah * _pk.c_rating)
    check("watt-hours are rated capacity x rated voltage, with no voltage "
          "cap in the rules", all(
              abs(p.watt_hours - p.capacity_ah * p.v_nom_cell * p.series)
              < 1e-9 for p in PACKS.values())
          and not any("volt" in f.name for f in __import__(
              "dataclasses").fields(Config().rules)))

    # ------------------------------------------------------------------
    section("Researched ESCs and packs")
    from dbf2027.components import ESCS as _ESCS, pick_esc as _pick
    _sc = _ESCS["ESC_SCORPION_TRIBUNUS_III_14_120A"]
    _cc = _ESCS["ESC_CC_PHOENIX_EDGE_HV_120"]
    check("Scorpion's 60 s / 32 m/s rating holds 72% x 80% for 5 min at "
          "10 m/s (its own derating table)",
          abs(_sc.sustained_a(300, 10, 28) / _sc.i_cont_a - 0.72 * 0.80)
          < 1e-9)
    check("Castle's 5 min / 5 mph rating holds in full, and never above "
          "the published figure",
          abs(_cc.sustained_a(300, 10, 28) - _cc.i_cont_a) < 1e-9
          and all(e.sustained_a(60, 40, 20) <= e.i_cont_a + 1e-9
                  for e in _ESCS.values()))
    _pk = _pick(70.0, 6, 300, 10, 36)
    check("the ESC picked is the lightest that holds the current",
          _pk is not None and _pk.sustained_a(300, 10, 36) >= 84.0
          and all(e.mass_kg >= _pk.mass_kg for e in _ESCS.values()
                  if e.min_cells <= 6 <= e.max_cells
                  and e.sustained_a(300, 10, 36) >= 84.0), _pk.key)
    check("every ESC has no BEC or a documented way to disable it",
          all("no BEC" in e.source or "disabled" in e.source
              for e in _ESCS.values()))
    _escs = set()
    for _preset in ("tucson-morning", "tucson", "tucson-hot"):
        _cw = Config()
        _cw.env = Environment.at_site(_preset)
        _escs.add(Aircraft(baseline(), _cw).esc.key)
    check("one ESC for every weather preset (it is hardware)",
          len(_escs) == 1, str(_escs))
    _ph = Aircraft(baseline(), cfg)
    _r2 = mission_2(_ph, cfg, dx=6.0)
    check("the missions report the pack's temperature, above the air's",
          _r2.pack_temp_c > cfg.env.temp_c + 1.0, "%.1f C" % _r2.pack_temp_c)
    _cold = Config()
    _cold.build.pack_temp_limit_c = cfg.env.temp_c + 2.0
    _r2c = mission_2(Aircraft(baseline(), _cold), _cold, dx=6.0)
    check("a pack that would pass its temperature limit ends the mission",
          not _r2c.ok and "battery pack would reach" in _r2c.reason,
          _r2c.reason)
    _cool = Config()
    _cool.build.pack_specific_heat_j_kgk = 1e9
    check("pack heating scales with 1 / heat capacity",
          mission_2(Aircraft(baseline(), _cool), _cool, dx=6.0).pack_temp_c
          < cfg.env.temp_c + 0.01)
    # 90% of the pack held back: Mission 2 runs out of battery, but the
    # takeoff and climb it did fly still heated the pack and windings.
    _flat = Config()
    _flat.build.m2_energy_margin = 0.9
    _r2e = mission_2(Aircraft(baseline(), _flat), _flat, dx=6.0)
    check("a mission that runs out of battery still reports the pack's "
          "temperature, above the air's",
          not _r2e.ok and "out of battery" in _r2e.reason
          and _r2e.pack_temp_c > _flat.env.temp_c
          and _r2e.motor_temp_c > _flat.env.temp_c,
          "%s: pack %.1f C, windings %.1f C"
          % (_r2e.reason, _r2e.pack_temp_c, _r2e.motor_temp_c))
    check("a published AC resistance is entered x1.7 (the flight-average "
          "DC figure)",
          abs(PACKS["TURNIGY_PANTHER_6S4000_75C"].r_cell_ohm
              - 16.0 / 6 * 1.7 / 1000.0) < 1e-5)   # stored to 3 figures
    _gt = PACKS["AULINE_45A_6S4000_LI"]
    check("a measured true capacity sets the energy; the rules keep the "
          "label", abs(_gt.energy_wh - 3.65 * 6 * 3.6) < 1e-9
          and abs(_gt.watt_hours - 4.0 * 6 * 3.6) < 1e-9)
    _hv = PACKS["SMC_HV3_6S4200_90C"]
    check("a LiHV pack is 3.80 V nominal and 4.35 V full per cell",
          _hv.v_nom_cell == 3.80 and _hv.v_full_cell == 4.35
          and abs(_hv.watt_hours - 4.2 * 6 * 3.8) < 1e-9)

    # ------------------------------------------------------------------
    section("Motor count: what differs between 1 and 4 motors")
    from dbf2027.aerodynamics import junction_drag_area_m2 as _jda
    from dbf2027.layout import (COPPER_RESISTIVITY_OHM_M as _rho_cu,
                                motor_stations as _ms)
    _mk = MOTORS["SCORPION_SII_4025_440"]
    _nt = _mk.kv * _mk.i0_volts / 60.0
    check("the quoted no-load current holds at the maker's test speed",
          abs(_mk.i0_at(_nt, 0.5) - _mk.i0) < 1e-12)
    check("its eddy share grows in proportion to speed; share 0 is the old "
          "constant I0",
          abs(_mk.i0_at(2.0 * _nt, 0.5) - 1.5 * _mk.i0) < 1e-12
          and abs(_mk.i0_at(2.0 * _nt, 0.0) - _mk.i0) < 1e-12)
    check("Hacker's no-load current is quoted at 8.4 V, the rest at 10 V",
          all((m.i0_volts == 8.4) == k.startswith("HACKER")
              for k, m in MOTORS.items()))
    _pk = PACKS["GA_6S3300_60C"]
    _ps0 = PropulsionSystem(_mk, PROPS["APC_12x8"], 1, _pk, cfg.env,
                            i0_speed_share=0.0)
    _ps5 = PropulsionSystem(_mk, PROPS["APC_12x8"], 1, _pk, cfg.env,
                            i0_speed_share=0.5)
    _psl = PropulsionSystem(_mk, PROPS["APC_12x8"], 1, _pk, cfg.env,
                            i0_speed_share=0.5, r_lead_ohm=0.010)
    check("more no-load current at speed costs thrust on the same battery",
          _ps5.thrust(20.0, 1.0, 1.0) < _ps0.thrust(20.0, 1.0, 1.0))
    check("lead resistance costs thrust",
          _psl.thrust(20.0, 1.0, 1.0) < _ps5.thrust(20.0, 1.0, 1.0))
    _one = Aircraft(baseline(), cfg)
    _three = Aircraft(baseline(n_motors=3), cfg)
    _four = Aircraft(baseline(n_motors=4), cfg)
    check("a nose motor has no wing leads, nacelles or nacelle junctions",
          _one.n_wing_motors == 0 and _one.lead_r_ohm == 0.0
          and _one.lead_mass_kg == 0.0
          and _one.aero_clean.nacelle_junction_m2 == 0.0
          and not any(k.startswith(("wing nacelles", "wing power leads"))
                      for k in _one.mass_breakdown.items))
    _n3, _ys3, _ = _ms(3, _three.wing.span_m, _three.body.fuse_width_m)
    _a3 = _three.motor.i_max_a / cfg.build.lead_current_density_a_mm2 * 1e-6
    _l3 = abs(_ys3[-1]) + cfg.build.lead_route_m
    check("each wing motor's leads run its station plus the routing, in "
          "copper sized to its rated current; the mean is over every motor",
          _three.n_wing_motors == 2
          and all(abs(ln - _l3) < 1e-12 for ln in _three.lead_lengths_m)
          and abs(_three.lead_r_ohm
                  - 2 * 2.0 * _l3 * _rho_cu / _a3 / 3.0) < 1e-12)
    check("the lead resistance reaches the power solve",
          abs(_three.prop_sys.r_lead - _three.lead_r_ohm) < 1e-15
          and _one.prop_sys.r_lead == 0.0)
    check("an outboard wing motor runs a longer lead",
          max(_four.lead_lengths_m) > max(_three.lead_lengths_m)
          and len(set(round(x, 9) for x in _four.lead_lengths_m)) == 2)
    _it3 = _three.mass_breakdown.items
    check("the wing nacelles and leads are weighed, and placed",
          abs(_it3["wing nacelles (2)"] - 2 * cfg.build.nacelle_mass_frac
              * (_three.motor.mass_kg + _three.esc.mass_kg)) < 1e-12
          and abs(_it3["wing power leads (2)"] - _three.lead_mass_kg) < 1e-12
          and abs(sum(q.mass_kg for q in _three.layout.stations)
                  - _three.mtow_m2_kg) < 1e-6)
    _a = _three.aero_clean
    check("each side of a wing nacelle is a Hoerner junction at its chord",
          _a.nacelle_junction_m2 > 0.0
          and abs(_a.nacelle_junction_m2 - sum(
              2.0 * _jda(_three.wing.airfoil.t_over_c,
                         _three.wing.chord_at(y)) for y in _ys3)) < 1e-15)
    _v3 = _three._max_level_speed(_three.w_m2_n)
    _f1 = _a.slipstream_drag_area_m2(_three.prop.diameter_m, _v3)
    _j, _a.nacelle_junction_m2 = _a.nacelle_junction_m2, 0.0
    _f0 = _a.slipstream_drag_area_m2(_three.prop.diameter_m, _v3)
    _a.nacelle_junction_m2 = _j
    check("the nacelle junctions sit in their propellers' slipstream",
          abs((_f1 - _f0) - _j) < 1e-15)

    # ------------------------------------------------------------------
    section("Computed corrections (each replaced a fixed factor)")
    from dbf2027.aerodynamics import (AeroModel, junction_drag_area_m2,
                                      tail_form_factor)
    from dbf2027.layout import section_perimeters
    _a = _iac.aero_clean
    _db = _iac.body.fuse_width_m / _iac.wing.span_m
    check("span efficiency carries the fuselage factor 1 - 2 (d/b)^2",
          abs(_a.e - _a.e_wing * (1.0 - 2.0 * _db * _db)) < 1e-12
          and _a.e < _a.e_wing, "%.4f" % _a.e)
    _wide = AeroModel(_iac.wing, replace(_iac.body,
                                         fuse_width_m=2 * _iac.body.fuse_width_m),
                      cfg.env, cfg.build, 1)
    check("a wider fuselage costs span efficiency", _wide.e < _a.e)
    _c0 = Config()
    _c0.build.fuselage_lift_carryover = 0.0
    _a0 = AeroModel(_iac.wing, _iac.body, _c0.env, _c0.build, 1)
    check("with no lift carried over, CLmax loses the covered share",
          abs(_a0.cl_max / _a.cl_max - (1.0 - _a.covered_frac)) < 1e-12
          and _a.fuselage_clmax_factor == 1.0,
          "%.4f" % (_a0.cl_max / _a.cl_max))
    check("CLmax is section x build factor x lifting-line stall x fuselage"
          " (no second finite-span factor)",
          abs(_a.cl_max - _iac.airfoil.cl_max_2d * cfg.build.cl_max_3d_factor
              * _a.stall_factor * _a.fuselage_clmax_factor) < 1e-12)
    _ok = True
    for _lam in (0.2, 0.43, 0.7, 1.0):
        _w = replace(_iac.wing, taper=_lam)
        _ok &= abs(_w.planform_within_m2(_w.span_m) - _w.area_m2) < 1e-12
        _ok &= abs(_w.tip_chord_m - _lam * _w.root_chord_m) < 1e-12
    check("the drawn root and tip chords have the modelled area", _ok)
    check("the MAC of a rectangular wing is its chord",
          abs(replace(_iac.wing, taper=1.0).mac_m - _iac.wing.chord_m)
          < 1e-12)
    check("Hoerner junction drag: 0.75 t/c - 0.0003/(t/c)^2 on t^2",
          abs(junction_drag_area_m2(0.14, 1.0)
              - (0.105 - 0.0003 / 0.0196) * 0.0196) < 1e-12
          and junction_drag_area_m2(0.04, 1.0) == 0.0)
    check("tail form factor is Raymer's at 30% thickness",
          abs(tail_form_factor(0.09) - (1.18 + 100 * 0.09 ** 4)) < 1e-12)
    _parts = _a.drag_parts(20.0)
    check("the excrescence allowance applies to the skin only",
          abs(_parts["excrescence"] - cfg.build.excrescence_frac
              * (_parts["wing"] + _parts["fuselage"] + _parts["tails"]
                 + _parts["slipstream transition"]))
          < 1e-15 and abs(sum(_parts.values()) - _a.cd0(20.0)) < 1e-15)
    _cn = Config()
    _cn.build.wing_covered_surfaces = 0
    _an = AeroModel(_iac.wing, _iac.body, _cn.env, _cn.build, 1)
    check("skin the wing and fuselage cover carries no friction",
          _an.cd0(20.0) > _a.cd0(20.0))
    _cs = Config()
    _cs.build.fuse_corner_radius_frac = 0.0
    _p_sharp, _ = section_perimeters(0.2, 0.2, _cs.build)
    _p_round, _ = section_perimeters(0.2, 0.2, cfg.build)
    check("fuselage perimeter: sharp box 2(w + h), default rounding ~0.92",
          abs(_p_sharp - 0.8) < 1e-12
          and 0.91 < _p_round / _p_sharp < 0.93,
          "%.4f" % (_p_round / _p_sharp))
    _d = _iac.prop.diameter_m
    _adisc = math.pi * 0.25 * _d * _d
    _vn, _kn = _iac.install_thrust_factor[12]
    _ai = _iac.aero_clean
    check("installed thrust = 1 - (slipstream drag area) / (disc area), "
          "at each speed",
          abs(_kn - (1.0 - _ai.slipstream_drag_area_m2(_d, _vn) / _adisc))
          < 1e-5)
    check("a bigger disc loses a smaller share to scrubbing",
          _ai.installed_thrust_factor(1.5 * _d, _vn)
          > _ai.installed_thrust_factor(_d, _vn))
    check("scrubbing falls with speed (skin friction falls with Reynolds "
          "number)",
          _ai.installed_thrust_factor(_d, 2.0 * _vn)
          >= _ai.installed_thrust_factor(_d, _vn))
    # ---- the propeller slipstream: turbulence and lift -------------------
    from dbf2027.liftingline import blown_wing, jet_speed_ratio, stall_factor
    _k = _a._re_exp
    _c0, _f0 = blown_wing(_a.AR, _iac.wing.taper, (), _k)
    check("a wing with no jet is the plain lifting line",
          abs(_c0 - 2.0 * math.pi * _a.AR * 0 - _c0) < 1e-15
          and abs(_f0 - stall_factor(_a.AR, _iac.wing.taper, _k)[0]) < 2e-3
          and blown_wing(_a.AR, _iac.wing.taper, ((0.0, 0.2, 1.0),), _k)
          == (_c0, _f0))
    check("the jet's edge is a smooth ramp keeping its speed",
          jet_speed_ratio(0.0, ((0.0, 0.2, 1.5),)) == 1.5
          and jet_speed_ratio(0.3, ((0.0, 0.2, 1.5),)) == 1.0
          and 1.0 < jet_speed_ratio(0.2, ((0.0, 0.2, 1.5),)) < 1.5)
    _jet = ((0.0, 0.2, 1.6),)
    _g24 = blown_wing(_a.AR, _iac.wing.taper, _jet, _k)[0] / _c0
    _g48 = blown_wing(_a.AR, _iac.wing.taper, _jet, _k, n_terms=48)[0] / _c0
    _w = replace(_iac.wing, span_m=1.0, chord_m=1.0 / _a.AR)
    _naive = 1.0 + (1.6 ** 2 - 1.0) * _w.planform_within_m2(0.2) / _w.area_m2
    # (Not necessarily BELOW the strip estimate: a jet centred on the span
    # sheds vortices at its edges that raise the wing just outside it.)
    check("the blown lift has converged, close to the strip estimate",
          abs(_g24 - _g48) < 5e-3 and _g24 > 1.0
          and abs(_g24 - _naive) < 0.10 * (_naive - 1.0) + 0.02,
          "24 terms %.4f, 48 terms %.4f, strip %.4f" % (_g24, _g48, _naive))
    _T = _iac.prop_sys.thrust(20.0, 1.0, 1.0)
    check("full power lifts more before the stall, at a lower CL",
          _a.max_lift(20.0, _T) > _a.max_lift(20.0)
          and _a.cl_required(20.0, 50.0, _T) < _a.cl_required(20.0, 50.0))
    _thr = lambda v: _iac.prop_sys.thrust(v, 1.0, 1.0)
    # The fuselage takes the middle of a nose propeller's slipstream
    # (AeroModel._jet_bands): this 12 in prop is barely wider than the
    # 185 mm fuselage at the top wing's height, so it blows almost none of
    # the wing; a 16 in one blows a real strip either side.
    _a16 = Aircraft(baseline(prop="APC_16x8"), cfg).aero_clean
    _thr16 = Aircraft(baseline(prop="APC_16x8"), cfg).prop_sys
    check("a nose prop blows only the wing either side of the fuselage",
          _a.slipstream_gains(20.0, _T)[0] < 1.02
          < _a16.slipstream_gains(20.0, _thr16.thrust(20.0, 1.0, 1.0))[0],
          "12x8 x%.3f, 16x8 x%.3f" % (
              _a.slipstream_gains(20.0, _T)[0],
              _a16.slipstream_gains(20.0, _thr16.thrust(20.0, 1.0, 1.0))[0]))
    _w16 = Aircraft(baseline(prop="APC_16x8"), cfg).w_m2_n
    check("the power-on stall is below the power-off stall",
          _a16.v_stall_powered(_w16, lambda v: _thr16.thrust(v, 1.0, 1.0))
          < _a16.v_stall(_w16))
    check("power off, drag is exactly the unblown polar",
          abs(_a.drag(20.0, 60.0, 0.0) - _a.drag(20.0, 60.0)) < 1e-12)
    _p = _a.drag_parts(20.0)
    check("slipstream turbulence: the washed strip moved part-way to "
          "fully turbulent",
          abs(_p["slipstream transition"]
              - _a.transition_frac * max(_a._turbulent_section_cd(20.0)
                                         - _a._section_cd(20.0), 0.0)
              * _a.washed_strip_m2 / _a.S) < 1e-15
          and _p["slipstream transition"] > 0.0)
    _ct = Config()
    _ct.build.slipstream_transition_frac = 0.0
    check("with no transition the slipstream adds no skin drag",
          AeroModel(_iac.wing, _iac.body, _ct.env, _ct.build, 1,
                    prop_diameter_m=_iac.prop.diameter_m).drag_parts(20.0)
          ["slipstream transition"] == 0.0)
    _wm = [p.windmill_drag_coeff for p in PROPS.values() if p.measured]
    check("every measured propeller gives a windmill drag, 0.02-0.15",
          _wm and all(c is not None and 0.02 < c < 0.15 for c in _wm),
          "%d props" % len(_wm))

    # ------------------------------------------------------------------
    section("Search progress")
    from dbf2027.optimizer import _pool_map, _Progress, _STAGES
    _nums = [float(i) for i in range(40)]
    for _procs in (1, 2):
        _seen = []
        _out = _pool_map(math.sqrt, _nums, _procs,
                         lambda d, n: _seen.append((d, n)))
        check("results come back in order with progress on (%d process%s)"
              % (_procs, "" if _procs == 1 else "es"),
              _out == [math.sqrt(x) for x in _nums])
        check("progress counts up to the total (%d process%s)"
              % (_procs, "" if _procs == 1 else "es"),
              bool(_seen) and _seen[-1] == (40, 40)
              and all(b[0] > a[0] for a, b in zip(_seen, _seen[1:])),
              "%d calls, last %s" % (len(_seen), _seen[-1:] or None))
    _log = []
    _tr = _Progress(lambda f, m: _log.append((f, m)),
                    {"screen": (30000, 1.0), "coarse": (30000, 1.0),
                     "refine": (36, 90.0), "final": (15, 1.0),
                     "weather": (15, 1.0)}, 15)
    for _k, _, _ in _STAGES:
        _cb = _tr.stage(_k)
        for _i in range(1, 6):
            _cb(_i, 5)
    _tr.finish()
    _fr = [f for f, _ in _log]
    check("overall progress only ever moves forward",
          all(b >= a - 1e-12 for a, b in zip(_fr, _fr[1:])))
    check("and ends at exactly 100%", _fr[-1] == 1.0
          and abs(_fr[-2] - 1.0) < 1e-9, "%.6f" % _fr[-2])
    check("every stage is named in its message",
          all(any(lbl in m for _, m in _log) for _, lbl, _ in _STAGES))
    check("no progress callback means the fast path",
          _Progress(None, {}).stage("coarse") is None)

    # ---- m74 / m75 ----------------------------------------------------
    from dbf2027.liftingline import (span_efficiency as _se,
                                     spanwise_loading as _swl,
                                     lift_curve_slope as _lcs)
    check("lifting line: rectangular AR 6 wing has the textbook span "
          "efficiency (Anderson 5.3.2, delta ~ 0.05)",
          0.945 < _se(6.0, 1.0) < 0.962, "e %.4f" % _se(6.0, 1.0))
    # the lift slope implied by the solve (CL = 2 pi AR A_1) against the
    # independent Helmbold formula
    import dbf2027.liftingline as _llm
    _a1 = []
    _orig = _llm.solve_linear
    _llm.solve_linear = lambda m_, r_: _a1.append(_orig(m_, r_)) or _a1[-1]
    _llm._CACHE.clear()
    _se(6.0, 1.0)
    _llm.solve_linear = _orig
    _llm._CACHE.clear()
    _slope = 2.0 * math.pi * 6.0 * _a1[-1][0]
    check("lifting line: the solve's lift slope matches Helmbold's at AR 6",
          abs(_slope / _lcs(6.0) - 1.0) < 0.01,
          "%.3f vs %.3f /rad" % (_slope, _lcs(6.0)))
    check("lifting line: a bigger aspect ratio is not a free span "
          "efficiency (rectangular AR 12 < AR 6)",
          _se(12.0, 1.0) < _se(6.0, 1.0))
    _a = ac.aero_clean
    _h = _a.x_np_over_c
    _fus = (cfg.build.fuselage_np_shift_k * ac.body.fuse_width_m ** 2
            * ac.body.fuse_length_m / (_a.S * ac.wing.mac_m * _a.a_wing))
    # dCm/dalpha about a CG at x/c = h, per unit q S c: wing, tail (arm AC
    # to AC, less the CG's offset), fuselage
    _cma = (_a.a_wing * (_h - 0.25)
            - _a.eta_tail * _a.s_tail_h / _a.S * _a.a_tail
            * (1.0 - _a.downwash)
            * (_a.tail_arm / ac.wing.mac_m - (_h - 0.25))
            + _fus * _a.a_wing)
    check("neutral point: the pitching moment's slope is zero there, tail "
          "lift included", abs(_cma) < 1e-9, "%.2e /rad" % _cma)
    _st4 = Aircraft(baseline(n_motors=4, n_systems=1,
                             wing_prop_layout="stagger"), cfg)
    _placed = (sum(s_.mass_kg for s_ in _st4.layout.stations)
               - _st4.m2_payload_kg)
    check("staggered props: every mass item is placed in the layout "
          "(the forward booms included)",
          abs(_placed - sum(_st4.mass_breakdown.items.values())) < 1e-9
          and any(s_.name == "stagger extensions"
                  for s_ in _st4.layout.stations),
          "%.4f vs %.4f kg" % (_placed,
                               sum(_st4.mass_breakdown.items.values())))
    from dbf2027.mission import max_level_speed as _mls
    _v79 = _mls(ac.aero_clean, ac.prop_sys, ac.w_m2_n, 0.7)
    _t79 = ac.prop_sys.thrust(_v79, 0.7, 1.0)
    check("top speed is where full-throttle thrust meets the POWER-ON drag "
          "(m79)", abs(_t79 - ac.aero_clean.drag(_v79, ac.w_m2_n, _t79))
          < 1e-6 * _t79, "%.2f m/s" % _v79)
    _v79a = ac._max_level_speed(ac.w_m2_n)
    _t79a = ac.prop_sys.thrust(_v79a, 0.85, 1.0)
    check("Aircraft's own top-speed solve uses the power-on drag too",
          abs(_t79a - ac.aero_clean.drag(_v79a, ac.w_m2_n, _t79a))
          < 1e-4 * _t79a)
    from dbf2027.mission import takeoff as _to82
    _w82 = 3.0 * ac.w_m2_n
    _r82 = _to82(ac.aero_clean, ac.prop_sys, 3.0 * ac.mtow_m2_kg, _w82, cfg)
    _v82 = _r82.v_liftoff
    _t82 = ac.prop_sys.thrust(_v82, 1.0, 1.0)
    check("a roll that stalls out below rotation speed fails at once, where "
          "thrust no longer beats drag and friction (m82)",
          not _r82.ok and "never reached rotation speed" in _r82.reason
          and _r82.ground_roll_m == 0.0, _r82.reason)
    check("and the same aeroplane at its own weight still takes off",
          _to82(ac.aero_clean, ac.prop_sys, ac.mtow_m2_kg, ac.w_m2_n,
                cfg).ok)
    # ---- m86-m88 --------------------------------------------------------
    from dbf2027 import optimizer as _o86
    from dbf2027.components import MOTORS as _M86, PROPS as _P86
    for _n86 in (1, 2):
        _d86 = baseline(n_motors=_n86, n_systems=1)
        _a86 = Aircraft(_d86, cfg)
        _r86 = _o86._reference(_d86, cfg)
        _ps86 = _o86.PropulsionSystem(_a86.motor, _a86.prop, _n86, _a86.pack,
                                      cfg.env, n_systems=1, esc=_a86.esc)
        _own86 = (1 + cfg.build.contingency_frac) * sum(
            _a86.propulsion_items.values())
        check("screen: a reference airframe's own propulsion set weighs what "
              "the full model says (%d motor%s)" % (_n86, "s" if _n86 > 1 else ""),
              abs(_o86._installed_mass(_ps86, _r86, cfg) - _own86) < 1e-9,
              "%.4f vs %.4f kg" % (_o86._installed_mass(_ps86, _r86, cfg), _own86))
    _heavy = max((m for m in _M86.values() if m.max_cells >= 6
                  and 0.3 < m.mass_kg < 0.8), key=lambda m: m.mass_kg)
    _ps86h = _o86.PropulsionSystem(_heavy, _a86.prop, 2, _a86.pack, cfg.env,
                                   n_systems=1)
    check("screen: a heavier motor costs more than its own mass (mounts, "
          "nacelles, contingency)",
          _o86._installed_mass(_ps86h, _r86, cfg)
          - _o86._installed_mass(_ps86, _r86, cfg)
          > (1 + cfg.build.contingency_frac) * 2 * (
              (_heavy.mass_kg + _ps86h.esc.mass_kg)
              - (_a86.motor.mass_kg + _a86.esc.mass_kg)))
    check("CAM folding props carry their maker's rpm limits",
          sum(1 for p in _P86.values() if p.max_rpm > 0) == 20
          and all(p.max_rpm > 0 for k, p in _P86.items() if k.startswith("CAM_")))
    _fast = next(m for m in _M86.values()
                 if m.max_cells >= 6 and m.kv * 25.2 > 1.2 * _P86["CAM_14x12"].max_rpm)
    _rc86 = _o86.PropulsionSystem(_fast, _P86["CAM_14x12"], 1, _a86.pack,
                                  cfg.env).rating_check()
    check("a motor that could spin a CAM prop past its limit is caught",
          _rc86["prop_rpm_margin"] < 1.0,
          "%s: %.0f rpm unloaded, limit %.0f" % (_fast.key, _rc86["rpm_no_load"],
                                                 _P86["CAM_14x12"].max_rpm))
    check("no propeller under 8 in, and no KDE 1806, in the catalogue",
          min(p.diameter_m for p in _P86.values()) >= 8 * 0.0254 - 1e-9
          and "KDE_1806XF_2350" not in _M86)
    _gm = [(h_, ground_mission(Aircraft(baseline(drop_height_in=h_), cfg),
                               cfg)[0]) for h_ in (12.0, 36.0, 46.0, 60.0)]
    check("ground mission: foam sized for the declared drop always passes "
          "(no rounding failures)", all(ok_ for _, ok_ in _gm),
          str([h_ for h_, ok_ in _gm if not ok_]))
    # ---- m101: the stall limit ----------------------------------------
    _a101 = Aircraft(baseline(), cfg)
    _vt101 = _a101.aero_clean.v_stall_in_air(_a101.w_m2_n,
                                             *_a101.stall_air()[:2])
    _vf101 = _a101.stall_speed_factored()
    check("the stall limit is on the lift over 1.5: the factored stall is "
          "close to sqrt(1.5) x the true one",
          1.15 * _vt101 < _vf101 < 1.2248 * _vt101,
          "%.2f vs %.2f m/s" % (_vf101, _vt101))
    _c101 = Config()
    _c101.build.max_stall_mps = _vf101 + 0.05
    _ok101 = not any("mph limit" in f_ for f_ in
                     Aircraft(baseline(), _c101).static_checks()["failures"])
    _c101.build.max_stall_mps = _vf101 - 0.05
    _no101 = Aircraft(baseline(), _c101).static_checks()
    check("a design is rejected when its factored stall is over the limit, "
          "and passes just under it",
          _ok101 and not _no101["ok"]
          and any("mph limit" in f_ for f_ in _no101["failures"]),
          "; ".join(_no101["failures"])[:120])
    # ---- m104: low aspect ratios and the runway attitude ----------------
    from dbf2027.optimizer import DesignSpace as _DS104
    _lo104 = Aircraft(baseline(chord_m=baseline().span_m / 2.0), cfg)   # AR 2
    _al104 = _lo104.aero_clean
    _hi104 = Aircraft(baseline(), cfg)
    _ah104 = _hi104.aero_clean
    check("the search goes down to aspect ratio 1.5 (m104)",
          _DS104().ar_bounds[0] == 1.5
          and not any("aspect ratio" in f_ for f_ in _lo104.static_checks()["failures"]))
    check("an AR-2 wing cannot reach its CLmax at the runway attitude, so it "
          "touches down faster than it stalls (m104)",
          _al104.cl_at_ground_alpha() < _al104.cl_max_at_v(20.0)
          and _al104.v_stall(_lo104.w_m2_n, ground=True)
          > _al104.v_stall(_lo104.w_m2_n) + 0.1,
          "CL %.2f vs CLmax %.2f" % (_al104.cl_at_ground_alpha(),
                                     _al104.cl_max_at_v(20.0)))
    check("a normal wing is not attitude-limited: the same speeds with or "
          "without the runway limit",
          _ah104.v_stall(_hi104.w_m2_n, ground=True)
          == _ah104.v_stall(_hi104.w_m2_n))
    # ---- m106: the chord limit -------------------------------------------
    _c106 = Config()
    _d106 = baseline(chord_m=0.40)
    _c106.build.max_chord_m = 0.40 * 1.001
    _ok106 = not any("in limit" in f_ for f_ in
                     Aircraft(_d106, _c106).static_checks()["failures"])
    _c106.build.max_chord_m = 0.40 * 0.999
    _no106 = Aircraft(_d106, _c106).static_checks()["failures"]
    check("a wing whose root chord is over the chord limit is rejected, one "
          "just under it is not (m106)",
          _ok106 and any("in limit" in f_ for f_ in _no106))
    check("the model takes aspect ratios down to 1.0",
          not any("aspect ratio" in f_ for f_ in Aircraft(
              baseline(chord_m=baseline().span_m / 1.2), cfg)
              .static_checks()["failures"]))
    # ---- m107: coarse samples inside the chord cap and stall limit -------
    from dbf2027.optimizer import _fit_sample, DesignSpace as _DS107
    from dbf2027.components import AIRFOILS as _AF107
    _c107 = Config()
    _b107 = _DS107().continuous_bounds(_c107)
    _x107 = [1.8, 2.0, 15.0, 0.8, 1.0, 40.0]
    _same107 = _fit_sample(_x107, _AF107["SD7062"], _b107, _c107) == _x107
    _c107.build.max_chord_m = 0.45
    _c107.build.max_stall_mps = 30 * 0.44704
    _f107 = _fit_sample(_x107, _AF107["SD7062"], _b107, _c107)
    _m107 = (0.5 * _c107.env.rho * _c107.build.max_stall_mps ** 2
             * _f107[0] ** 2 / _f107[1] * 1.5 * _AF107["SD7062"].cl_max_2d
             / (1.5 * _c107.env.g))
    check("coarse samples are moved inside the chord cap and under what the "
          "wing could lift at the stall limit, and left alone without them "
          "(m107)", _same107 and abs(_f107[1] - 1.8 / 0.45) < 1e-12
          and _f107[2] <= _m107 + 1e-9 and _f107[2] < 15.0 and _f107[0] == 1.8,
          str(_f107[:3]))
    # ---- m108: a sample over the stall limit is brought down to it -------
    from dbf2027.optimizer import _worker_fit_stall
    _c108 = Config()
    _a108 = Aircraft(baseline(), _c108)
    _c108.build.max_stall_mps = 0.9 * _a108.stall_speed_factored()
    _r108 = _worker_fit_stall((baseline(), _c108, 12.0, None, None, None))
    _d108 = _r108.get("design")
    check("under a stall limit a coarse sample that is over it comes back with "
          "a lighter sensor that meets it (m108)",
          _d108 is not None and _d108.sensor_mass_kg < baseline().sensor_mass_kg
          and Aircraft(_d108, _c108).stall_speed_factored()
          <= _c108.build.max_stall_mps + 1e-9,
          "%s" % (_r108.get("reason") or "")[:80])
    # ---- m120: the stall fit by secant ------------------------------------
    from dbf2027 import optimizer as _op120
    _n120 = [0]
    _ac120 = _op120.Aircraft

    class _Count120(_ac120):
        def __init__(self, *a, **k):
            _n120[0] += 1
            super().__init__(*a, **k)
    _op120.Aircraft = _Count120
    try:
        _r120 = _worker_fit_stall((baseline(), _c108, 12.0, None, None, None))
    finally:
        _op120.Aircraft = _ac120
    _d120 = _r120.get("design")
    _v120 = (Aircraft(_d120, _c108).stall_speed_factored()
             if _d120 is not None else float("inf"))
    check("the stall fit lands on the limit, inside it, in a few aircraft "
          "builds (a secant, m120)",
          _d120 is not None and _v120 <= _c108.build.max_stall_mps + 1e-12
          and _v120 >= _c108.build.max_stall_mps * (1.0 - 2e-6)
          and _n120[0] <= 7,
          "%.7f of the limit in %d builds"
          % (_v120 / _c108.build.max_stall_mps, _n120[0]))
    # ---- m121: an empty coarse sweep goes on from the seeds --------------
    import inspect as _in121
    _src121 = _in121.getsource(_op120.optimise)
    check("a search whose coarse sweep flies nothing goes on from its seed "
          "designs, and stops only when none of them flies (m121)",
          "nothing flew in the coarse sweep" in _src121
          and "the seed designs - check" in _src121
          and "no feasible design in the coarse sweep - check" not in _src121
          and _src121.index("if not seeds:") > _src121.index("seed_rows("))
    # ---- m109: section lift data and the finite-wing slope -------------
    from dbf2027.components import AIRFOILS as _AF109
    from dbf2027.liftingline import lift_curve_slope as _lcs109
    check("every catalogue section carries its XFOIL zero-lift angle and "
          "lift slope (m109)",
          all(a.alpha_l0_deg is not None and 5.0 < a.cl_alpha_per_rad < 7.5
              and -12.0 < a.alpha_l0_deg <= 0.0 for a in _AF109.values()))
    check("the finite-wing lift slope tends to the section's at long aspect "
          "ratio and to pi AR / 2 at short (DATCOM/Polhamus)",
          abs(_lcs109(1e4, 6.0) - 6.0) < 1e-2
          and abs(_lcs109(0.5, 6.0) - math.pi * 0.25) < 0.02 * math.pi * 0.25
          and abs(_lcs109(6.0) - 2 * math.pi * 6 / (2 + math.sqrt(40))) < 1e-12)
    _s109 = Aircraft(baseline(airfoil="S1223"), cfg).aero_clean
    check("a high-lift section's zero-lift angle is its polar's, not the "
          "camber-moment estimate (S1223 -10.0 deg)",
          abs(math.degrees(_s109.alpha_zero_lift) + 10.008) < 1e-9)
    # ---- m110: winglets ---------------------------------------------------
    _b110 = baseline()
    _w110 = replace(_b110, winglet_h_m=0.10 * _b110.span_m)
    _a0_110, _a1_110 = Aircraft(_b110, cfg), Aircraft(_w110, cfg)
    check("winglets 10% of the span raise the effective aspect ratio by "
          "Hoerner's 1 + 1.9 h/b and cut the induced drag factor to match (m110)",
          abs(_a1_110.aero_clean.AR_eff / _a1_110.aero_clean.AR - 1.19) < 1e-9
          and abs(_a0_110.aero_clean.k_induced / _a1_110.aero_clean.k_induced
                  - 1.19) < 1e-6)
    check("winglets cost mass and skin drag, and steepen the lift slope",
          _a1_110.mass_breakdown.items.get("wing skin (winglets)", 0) > 0
          and _a1_110.aero_clean.drag_parts(20.0)["winglets"] > 0
          and _a0_110.aero_clean.drag_parts(20.0)["winglets"] == 0
          and _a1_110.aero_clean.a_wing > _a0_110.aero_clean.a_wing)
    from dbf2027.optimizer import DesignSpace as _DS110, _make_design as _md110
    check("the search carries winglet height as its seventh variable",
          len(_DS110().continuous_bounds(cfg)) >= 7
          and _DS110.VAR_NAMES[6] == "winglet_frac"
          and abs(_md110({"motor": _b110.motor, "prop": _b110.prop, "n_motors": 1,
                          "pack": _b110.pack, "n_systems": 1},
                         _b110.airfoil, [1.5, 5, 5, 1, 1, 60, 0.1], 0, 1.0)
                  .winglet_h_m - 0.15) < 1e-12)
    # ---- m111: tails sized from requirements ----------------------------
    _a111 = Aircraft(baseline(), cfg)
    _t111 = _a111.tail_sizing
    _ct111 = _a111.tail_trim_cl()[0]
    check("the horizontal tail is sized from its requirements and trims "
          "across the flight envelope within its headroom (m111)",
          _t111.get("mode") == "requirements"
          and _ct111 <= cfg.build.tail_trim_cl_frac * cfg.build.tail_cl_max + 0.01
          and _a111.body.s_horiz_m2 >= 0.99 * max(_t111["s_h_trim_m2"],
                                                  _t111["s_h_rotation_m2"],
                                                  _t111["s_h_floor_m2"]),
          "tail CL %.2f, sized by %s" % (_ct111, _t111.get("h_sized_by")))
    _2m111 = Aircraft(baseline(n_motors=2), cfg)
    check("a two-wing-motor fin holds one motor out; a single nose motor "
          "needs no such fin",
          _2m111.tail_sizing["s_v_motor_out_m2"] > 0
          and _2m111.body.s_vert_m2
          >= 0.99 * _2m111.tail_sizing["s_v_motor_out_m2"]
          and _a111.tail_sizing["s_v_motor_out_m2"] == 0.0,
          "%s" % _2m111.tail_sizing.get("v_sized_by"))
    _big112 = Aircraft(baseline(tail_vol_h=0.80), cfg)
    _vh112 = (_big112.body.s_horiz_m2 * _big112.body.tail_arm_m
              / (_big112.wing.area_m2 * _big112.wing.mac_m))
    check("a chosen tail volume above the requirements is flown, with less "
          "cruise download; one below them is not (m112)",
          abs(_vh112 - 0.80) < 1e-4
          and _big112.tail_sizing["h_sized_by"] == "choice"
          and abs(_big112.aero_clean.tail_load(22.0, _big112.w_m2_n))
          < abs(_a111.aero_clean.tail_load(22.0, _a111.w_m2_n))
          and abs(Aircraft(baseline(tail_vol_h=0.05), cfg).body.s_horiz_m2
                  / _a111.body.s_horiz_m2 - 1.0) < 1e-4, "Vh %.3f" % _vh112)
    from dbf2027.optimizer import DesignSpace as _DS112
    check("the tail volume is the 8th search variable",
          len(_DS112().continuous_bounds(cfg)) == 8
          and _DS112.VAR_NAMES[7] == "tail_vol_h")
    from dbf2027.optimizer import (SearchSettings as _SS114, _filled,
                                   _extra_starts)
    _s114 = _SS114(processes=64)
    _f114 = _SS114(processes=64, fill_cores=True)
    check("filling the cores rounds a stage up to whole rounds of the "
          "workers, and changes nothing unless asked (m114)",
          [_filled(n, _f114) for n in (0, 18, 32, 64, 100)]
          == [0, 64, 64, 64, 128]
          and [_filled(n, _s114) for n in (18, 100)] == [18, 100])
    _a114 = [{"design": baseline(span_m=1.50 + 0.01 * i)} for i in range(3)]
    _b114 = [{"design": baseline(airfoil="E205", span_m=1.50)}]
    _x114 = _extra_starts(_a114 + _b114, [_a114[0], _b114[0]], 5)
    check("too few configurations: the best get their next-best samples "
          "as more starts, never a duplicate (m114)",
          [round(r["design"].span_m, 4) for r in _x114] == [1.51, 1.52]
          and not _extra_starts(_a114[:1] + _a114[:1], _a114[:1], 3))
    _a116 = Aircraft(baseline(n_motors=2), cfg).aero_clean
    _dz116 = _a116.thrust_line_dz_m
    _lt1 = _a116.tail_load(18.0, 100.0, 20.0)
    _a116.thrust_line_dz_m = 0.0
    _lt0 = _a116.tail_load(18.0, 100.0, 20.0)
    _a116.thrust_line_dz_m = _dz116
    _lm116 = _a116.max_lift(18.0, 20.0)
    check("the trim carries the wing propellers' thrust line above the CG, "
          "in the tail load and the lift at the stall alike (m116)",
          _dz116 > 0.0
          and abs(_lt1 - _lt0 + 20.0 * _dz116 / _a116.tail_arm) < 1e-9
          and abs(_lm116 - _a116.tail_load(18.0, _lm116, 20.0)
                  - _a116.q(18.0) * _a116.S * _a116.cl_max_at_v(18.0)
                  * _a116.slipstream_gains(18.0, 20.0)[1]) < 1e-9 * _lm116
          and Aircraft(baseline(), cfg).aero_clean.thrust_line_dz_m == 0.0,
          "dz %.3f m" % _dz116)
    from dbf2027.weights import (_end_plated_moment, _end_plated_shear,
                                 _elliptic_moment, spar_design)
    _r117 = (1.0 + cfg.build.winglet_ar_gain * 0.15) ** 0.5
    _w117 = Aircraft(baseline(), cfg)
    # a load the minimum gauge does not cover
    _sp0 = spar_design(_w117.wing, 300.0, cfg.build, 8.0).mass_kg
    from dataclasses import replace as _rep117
    _sp1 = spar_design(_rep117(_w117.wing, winglet_h_m=0.15 * _w117.wing.span_m),
                       300.0, cfg.build, 8.0).mass_kg
    check("winglets load the tips: the spar's bending is the end-plated "
          "wing's, 6.6% more at the root for h/b 0.15, elliptic without "
          "them (m117)",
          abs(_end_plated_moment(0.0, _r117) / _elliptic_moment(0.0) - 1.066)
          < 0.002
          and max(abs(_end_plated_moment(e, 1.0 + 1e-9) - _elliptic_moment(e))
                  for e in (0.0, 0.3, 0.7, 0.95)) < 1e-6
          and abs(_end_plated_shear(0.0, _r117) - 1.0) < 1e-12
          and _end_plated_moment(1.0, _r117) == 0.0
          and _end_plated_shear(1.0, _r117) == 0.0 and _sp1 > _sp0,
          "root x%.4f" % (_end_plated_moment(0.0, _r117) / _elliptic_moment(0.0)))
    from dbf2027.config import Environment as _Env118
    _a118 = Aircraft(baseline(), cfg)
    _hot118 = _Env118.at_site("tucson-hot")
    _v_nom = _a118.aero_clean.v_stall(_a118.w_m2_n, load_factor=1.5, ground=True)
    _c118 = Config()
    _c118.env = _hot118
    _a118h = Aircraft(baseline(), _c118)
    _c118m = Config()
    _c118m.build.stall_check_conditions = ""
    check("the stall limit is checked in the hottest preset's air: the same "
          "factored stall the hot day's own aeroplane has, the mission's air "
          "when asked (m118)",
          cfg.build.stall_check_conditions == "tucson-hot"
          and abs(_a118.stall_speed_factored()
                  - _a118h.aero_clean.v_stall(_a118h.w_m2_n, load_factor=1.5,
                                              ground=True)) < 1e-9
          and _a118.stall_speed_factored() > _v_nom * 1.01
          and abs(Aircraft(baseline(), _c118m).stall_speed_factored() - _v_nom)
          < 1e-12,
          "%.3f vs %.3f m/s nominal" % (_a118.stall_speed_factored(), _v_nom))
    from dbf2027.mission import mission_2 as _m2_119, turn_capability as _tc119
    _a119 = Aircraft(baseline(), cfg)
    _r119 = _m2_119(_a119, cfg, 12.0)
    _t119 = _tc119(_a119.aero_clean, _a119.prop_sys, _a119.w_m2_n, 0.70, cfg,
                   _a119.n_limit)
    _n119 = _r119.turn_load_factor
    check("the results give the lap's turn: radius, load factor, bank and rate, "
          "the corner-speed turn the mission flies (m119)",
          _r119.ok and _r119.turn_radius_m > 0.0
          and abs(_r119.turn_radius_m - _t119.radius_m) < 1e-9
          and abs(_r119.v_turn - _t119.v_turn) < 1e-12
          and abs(_r119.turn_radius_m - _r119.v_turn ** 2
                  / (cfg.env.g * ((_n119 - 1.0) * (_n119 + 1.0)) ** 0.5))
          < 1e-6 * _r119.turn_radius_m
          and abs(math.tan(math.radians(_r119.turn_bank_deg))
                  - ((_n119 - 1.0) * (_n119 + 1.0)) ** 0.5) < 1e-9,
          "%.1f m at %.2f g" % (_r119.turn_radius_m, _n119))
    _ac115 = Aircraft(baseline(), cfg)
    _s115 = _ac115.summary()
    check("the results give the tail and fin areas, spans, chords and "
          "thickness, the tail arm, the fuselage and the stations, consistent "
          "with the tail volumes (m115)",
          abs(_s115["tail_area_m2"] * _s115["tail_arm_m"]
              / (_ac115.wing.area_m2 * _ac115.wing.mac_m)
              - _s115["vol_coeff_h"]) < 1e-12
          and abs(_s115["tail_span_m"] * _s115["tail_chord_m"]
                  - _s115["tail_area_m2"]) < 1e-12
          and abs(_s115["fin_height_m"] * _s115["fin_chord_m"]
                  - _s115["fin_area_m2"]) < 1e-12
          and abs(_s115["tail_thickness_m"] - cfg.build.tail_t_over_c
                  * _s115["tail_chord_m"]) < 1e-12
          and 0.0 < _s115["x_wing_le_m"] < _s115["x_cg_m"]
          < _s115["x_tail_le_m"] < _s115["fuse_length_m"]
          and _s115["fuse_length_m"] > _s115["tail_arm_m"] > 0.0)
    _c111 = Config()
    _c111.build.tail_sizing = "volume"
    _f111 = Aircraft(baseline(), _c111)
    check("tail_sizing 'volume' keeps the fixed tail volumes",
          abs(_f111.body.s_horiz_m2 * _f111.body.tail_arm_m
              / (_f111.wing.area_m2 * _f111.wing.mac_m)
              - _c111.build.vol_coeff_h) < 1e-9)
    # ---- m99 ------------------------------------------------------------
    # A near neighbour (chord 0.002% longer) shares every rounded table key
    # with the baseline; before m99 the baseline scored after it got its
    # numbers (5e-5 s on Mission 2).
    _b99 = baseline()
    _f99 = _o86.evaluate_design(_b99, cfg, 12.0)
    _o86.evaluate_design(replace(_b99, chord_m=_b99.chord_m * 1.00002), cfg, 12.0)
    _a99 = _o86.evaluate_design(_b99, cfg, 12.0)
    _k99 = ("objective", "m2_time_s", "m3_time_s", "m2_energy_margin")
    check("a design scores the same whatever was scored before it in the "
          "process (m99)", all(_f99.get(k_) == _a99.get(k_) for k_ in _k99),
          str([(k_, _f99.get(k_), _a99.get(k_)) for k_ in _k99
               if _f99.get(k_) != _a99.get(k_)]))

    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("%d passed, %d failed" % (PASS, FAIL))
    print("=" * 60)
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
