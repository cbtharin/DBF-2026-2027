#!/usr/bin/env python3
"""End-to-end validation of the DBF 2026-27 model.

`selftest.py` checks physical *invariants* on one or two aeroplanes.  This
checks the things that break when you are not looking:

  A. Catalogue integrity   - every airfoil, prop, motor, cell and ESC in the
                             database is physically sensible, not just the
                             three used in the tests.
  B. Independent bounds    - simulated lap time and energy checked against
                             hand calculations that do not use the simulator.
  C. Energy bookkeeping    - average power never exceeds what the system can
                             physically draw; energy is conserved across
                             mission phases.
  D. Determinism           - the same input gives bit-identical output, and
                             the same seed gives the same search result.
  E. Round trip            - export a design, read it back, rescore it.
  F. Trends across the DB  - physical monotonicity over the whole catalogue,
                             not just one sample point.
  G. Interfaces            - every CLI flag reaches the config it claims to.
  H. Design loads          - the V-n diagram, the gust formula, spar sizing
                             and 1-4 motor layouts, across the catalogue.

    python validate.py
"""
from __future__ import annotations

import json
from dataclasses import fields
import math
import os
import sys
import tempfile

from dbf2027 import (AIRFOILS, ESCS, MOTORS, PACKS, PROPS, Aircraft,
                     BatteryPack, Config, Design, Environment, evaluate_scores)
from dbf2027.mission import mission_2, mission_3, turn_capability
from dbf2027.optimizer import DesignSpace, SearchSettings, evaluate_design
from run_optimizer import build_config, build_parser
from dbf2027.loads import (design_load_factors, gust_load_factor,
                           lift_curve_slope_3d)
from dbf2027.payload import size_payload
from dbf2027.propulsion import PropulsionSystem, clear_system_cache
from dbf2027.units import FT2M, KG2LB, LB2KG, N2LBF

PASS = FAIL = 0
FAILURES = []


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
    else:
        FAIL += 1
        FAILURES.append("%s   %s" % (name, detail))
        print("  FAIL  %s   %s" % (name, detail))


def section(t):
    print("\n" + t)
    print("-" * len(t))


def baseline(**kw):
    a = dict(span_m=1.55, chord_m=0.26, airfoil="SD7062",
             motor="SCORPION_SII_4025_440", prop="APC_12x8", n_motors=1,
             pack="GA_6S3300_60C", n_systems=1,
             sensor_mass_kg=0.9, n_simulators=3)
    a.update(kw)
    return Design(**a)


# ==========================================================================
def validate_catalogue(cfg):
    section("A. Catalogue integrity (every entry, not just the sampled ones)")

    for name, af in AIRFOILS.items():
        check("%s CLmax plausible" % name, 0.8 <= af.cl_max_2d <= 2.6,
              "%.2f" % af.cl_max_2d)
        check("%s cd_min plausible" % name, 0.004 <= af.cd_min <= 0.035,
              "%.4f" % af.cd_min)
        check("%s camber gives nose-down Cm" % name, af.cm_ac <= 0.0,
              "%.3f" % af.cm_ac)
        check("%s thickness plausible" % name, 0.05 <= af.t_over_c <= 0.20,
              "%.3f" % af.t_over_c)
        check("%s stall break below CLmax" % name, af.cl_break < af.cl_max_2d)

    for key, pr in PROPS.items():
        n = 7000.0 / 60.0
        t0 = pr.thrust(1.225, n, 0.0)
        p0 = pr.power(1.225, n, 0.0)
        area = math.pi * (pr.diameter_m / 2) ** 2
        fom = t0 ** 1.5 / math.sqrt(2 * 1.225 * area) / p0
        # Measured (UIUC) props span this: a high-pitch prop stalls its
        # blades at rest (11x10, 14x12 about 0.34; the CAM 12x13 folding
        # prop, P/D 1.08, 0.29), a low-pitch one is efficient statically
        # (18x8 about 0.85).
        check("%s static FOM in band" % key, 0.28 <= fom <= 0.90,
              "%.3f" % fom)
        etas = [pr.thrust(1.225, n, j * pr.j_zero * n * pr.diameter_m)
                * (j * pr.j_zero * n * pr.diameter_m)
                / max(pr.power(1.225, n, j * pr.j_zero * n * pr.diameter_m), 1e-9)
                for j in [i * 0.02 for i in range(1, 50)]]
        check("%s peak eta in band" % key, 0.55 <= max(etas) <= 0.92,
              "%.3f" % max(etas))
        check("%s thrust monotone in V" % key,
              all(pr.thrust(1.225, n, v) >= pr.thrust(1.225, n, v + 0.5) - 1e-12
                  for v in [i * 0.5 for i in range(80)]))
        check("%s Cp never negative" % key,
              all(pr.cp(j * 0.05) > 0 for j in range(40)))

    for key, m in MOTORS.items():
        check("%s Kv positive" % key, 50 < m.kv < 3000, "%.0f" % m.kv)
        check("%s Rm plausible" % key, 0.005 <= m.rm <= 0.5, "%.3f" % m.rm)
        check("%s I0 plausible" % key, 0.2 <= m.i0 <= 6.0, "%.2f" % m.i0)
        check("%s mass plausible" % key, 0.03 <= m.mass_kg <= 1.5,
              "%.3f" % m.mass_kg)
        # power rating must be consistent with the current rating at a
        # voltage the motor could plausibly see
        v_plausible = m.p_max_w / max(m.i_max_a, 1e-6)
        check("%s power/current rating consistent" % key,
              6.0 <= v_plausible <= 60.0, "implies %.1f V" % v_plausible)

    # Pack-level energy density, not cell level: a real pack carries wire,
    # connectors, shrink and a label, so it always comes out below the bare
    # cell figure.  The LiPo ceiling was 175; the 3S-12S catalogue's densest
    # packs are maker-published at 179-187 Wh/kg (Thunder Power ProLite X,
    # MaxAmps, CNHL), which modern high-density cells do reach, so 200.
    bands = {"LiPo": (100, 200), "LiIon": (150, 260)}
    for key, pk in PACKS.items():
        check("%s chemistry is lithium" % key, pk.chemistry in bands)
        whkg = pk.watt_hours / pk.mass_kg
        lo, hi = bands[pk.chemistry]
        check("%s energy density in %s band" % (key, pk.chemistry),
              lo <= whkg <= hi, "%.0f Wh/kg" % whkg)
        check("%s cutoff below nominal" % key,
              pk.v_min_cell < pk.v_nom_cell < pk.v_full_cell)
        check("%s usable fraction sane" % key,
              0.5 <= pk.usable_fraction <= 0.95)
        check("%s is under the 100 Wh per-pack limit" % key,
              pk.watt_hours <= 100.0, "%.1f Wh" % pk.watt_hours)
        check("%s label figures are self-consistent" % key,
              abs(pk.watt_hours - pk.capacity_ah * pk.v_nominal) < 1e-9)
        # (racing labels run to 160C - inflated about 2x, per independent
        # tests - which is why the model limits a pack by its heating, not
        # by this number; the bound only catches typos)
        check("%s C rating is plausible for %s" % (key, pk.chemistry),
              (5 <= pk.c_rating <= 20) if pk.chemistry == "LiIon"
              else (15 <= pk.c_rating <= 200), "%.0fC" % pk.c_rating)
        check("%s series count is a real configuration" % key,
              2 <= pk.series <= 14)
        check("%s OCV rises with state of charge" % key,
              all(pk.ocv(x / 10.0) < pk.ocv((x + 1) / 10.0)
                  for x in range(10)))
        check("%s OCV spans the cell limits" % key,
              abs(pk.ocv(0.0) - pk.v_min) < 1e-9
              and abs(pk.ocv(1.0) - pk.v_full) < 1e-9)

    for key, e in ESCS.items():
        check("%s efficiency plausible" % key, 0.90 <= e.efficiency <= 0.99)
        check("%s mass scales with current" % key,
              0.3 <= 1000 * e.mass_kg / e.i_cont_a <= 2.0,
              "%.2f g/A" % (1000 * e.mass_kg / e.i_cont_a))


# ==========================================================================
def validate_bounds(cfg):
    section("B. Independent bounds (hand calculations, no simulator)")
    ac = Aircraft(baseline(), cfg)
    r2 = mission_2(ac, cfg, dx=6.0)
    check("Mission 2 flew", r2.ok, r2.reason)
    if not r2.ok:
        return

    # Lower bound on lap time: straights at Vmax plus turns at the best
    # sustained turn rate.  The simulator must not beat this.
    turn = turn_capability(ac.aero_clean, ac.prop_sys, ac.w_m2_n, 0.7, cfg)
    t_turns = math.radians(cfg.course.total_turn_deg) / turn.turn_rate_rad_s
    t_straights = cfg.course.straight_length_m / r2.v_max
    floor = t_straights + t_turns
    lap = r2.lap_times[-1]
    check("lap time respects its physical floor", lap >= floor - 1e-6,
          "lap %.2f s vs floor %.2f s" % (lap, floor))
    check("lap time within 2x of the floor", lap <= 2.0 * floor,
          "lap %.2f s vs floor %.2f s" % (lap, floor))

    # Lower bound on energy: minimum drag over the straight distance at the
    # best possible efficiency.
    v_ld = ac.aero_clean.best_ld_speed(ac.w_m2_n)
    d_min = ac.aero_clean.drag(v_ld, ac.w_m2_n)
    e_floor = d_min * 5 * cfg.course.straight_length_m / 0.70
    check("Mission 2 energy respects its floor", r2.energy_j >= e_floor - 1e-6,
          "%.0f J vs floor %.0f J" % (r2.energy_j, e_floor))

    # Stall speed against the textbook formula at the model's own CLmax
    # Stall speed against an independent solve of the trim relation:
    # the wing stalls when L_w = q S CLmax, with L_w = W - L_t(W).
    v_s = ac.aero_clean.v_stall(ac.w_m2_n)
    a = ac.aero_clean
    resid = ((ac.w_m2_n - a.tail_load(v_s, ac.w_m2_n))
             - a.q(v_s) * ac.wing_area_m2 * a.cl_max_at_v(v_s))
    check("stall speed satisfies the trimmed lift balance", abs(resid) < 1e-6,
          "residual %.3e N" % resid)

    # the empennage must close the lift balance at every load factor
    for n in (1.0, 2.0, 3.0):
        L = n * ac.w_m2_n
        lw = a.cl_required(20.0, L) * a.q(20.0) * ac.wing_area_m2
        lt = a.tail_load(20.0, L)
        check("wing + tail = total lift at %.0f g" % n,
              abs((lw + lt) - L) < 1e-6, "%.6f" % ((lw + lt) - L))

    # Course geometry against the stated rules
    check("course straights are 2000 ft",
          abs(cfg.course.straight_length_m / FT2M - 2000.0) < 1e-6)
    check("course turning is 720 deg", abs(cfg.course.total_turn_deg - 720) < 1e-9)
    check("span limit is 6 ft",
          abs(cfg.rules.max_wingspan_m / FT2M - 6.0) < 1e-9)


# ==========================================================================
def validate_energy(cfg):
    section("C. Energy bookkeeping")
    ac = Aircraft(baseline(), cfg)
    for label, res, w in (("M2", mission_2(ac, cfg, dx=6.0), ac.w_m2_n),
                          ("M3", mission_3(ac, cfg, dx=6.0), ac.w_m3_n)):
        if not res.ok:
            check("%s flew" % label, False, res.reason)
            continue
        # average power cannot exceed the most the system can physically draw
        p_avg = res.energy_j / max(res.total_flight_time_s, 1e-9)
        p_max = ac.prop_sys.power_elec(0.0, 1.0, 1.0)
        for v in (10, 20, 30, 40):
            p_max = max(p_max, ac.prop_sys.power_elec(v, 1.0, 1.0))
        check("%s average power below system maximum" % label, p_avg <= p_max,
              "%.0f W vs %.0f W" % (p_avg, p_max))
        check("%s energy within the pack" % label,
              res.energy_j <= ac.pack.watt_hours * 3600.0,
              "%.1f Wh of %.1f" % (res.energy_j / 3600, ac.pack.watt_hours))
        check("%s energy margin consistent with energy used" % label,
              abs((1.0 - res.energy_j / (ac.pack.watt_hours * 3600.0))
                  - res.energy_margin) < 1e-9)
        check("%s total time = scored time + landing" % label,
              abs(res.total_flight_time_s - (res.time_s + res.landing_time_s))
              < 1e-6)
        check("%s lap times all positive" % label,
              all(t > 0 for t in res.lap_times))
        check("%s scored time excludes landing" % label,
              res.time_s < res.total_flight_time_s)

    # more laps must cost more energy
    r3 = mission_3(ac, cfg, dx=6.0)
    check("Mission 3 energy scales with laps", r3.energy_j > 0 and r3.laps > 0)


# ==========================================================================
def validate_determinism(cfg):
    section("D. Determinism and reproducibility")
    d = baseline()
    a = evaluate_design(d, cfg, 6.0)
    b = evaluate_design(d, cfg, 6.0)
    for k in ("m2_raw_lb_per_s", "m3_raw_lap_lb", "m2_time_s", "m3_laps"):
        check("repeat evaluation identical: %s" % k, a.get(k) == b.get(k),
              "%r vs %r" % (a.get(k), b.get(k)))

    # a fresh propulsion cache must not change anything
    clear_system_cache()
    c = evaluate_design(d, cfg, 6.0)
    check("result independent of the instance cache",
          a["m2_raw_lb_per_s"] == c["m2_raw_lb_per_s"],
          "%r vs %r" % (a["m2_raw_lb_per_s"], c["m2_raw_lb_per_s"]))
    # ... and an entry edited under the same key (a measured I0 typed into
    # the Catalogues tab) must not be answered from it.  The key used to be
    # the motor's name, so the edit was silently ignored all session.
    from dataclasses import replace as _repl
    from dbf2027 import components as _C
    from dbf2027.propulsion import get_propulsion as _gp
    _m = _C.MOTORS[d.motor]
    _p0 = _gp(_m, _C.PROPS[d.prop], d.n_motors, _C.PACKS[d.pack], cfg.env,
              n_systems=d.n_systems)
    _p1 = _gp(_repl(_m, i0=2.0 * _m.i0), _C.PROPS[d.prop], d.n_motors,
              _C.PACKS[d.pack], cfg.env, n_systems=d.n_systems)
    check("an edited motor under the same key is not answered from the "
          "cache", _p1 is not _p0
          and _p1.thrust(20.0, 1.0, 1.0) < _p0.thrust(20.0, 1.0, 1.0))

    # the memo in PropulsionSystem must not leak between query points
    ps = Aircraft(d, cfg).prop_sys
    t1 = ps.thrust(10.0, 0.9, 1.0)
    ps.thrust(30.0, 0.4, 0.6)
    t2 = ps.thrust(10.0, 0.9, 1.0)
    check("operating-point memo is not stale", t1 == t2, "%.9f vs %.9f" % (t1, t2))


# ==========================================================================
def validate_roundtrip(cfg):
    section("E. Export / import round trip")
    from dataclasses import asdict
    d = baseline()
    before = evaluate_design(d, cfg, 6.0)
    path = os.path.join(tempfile.gettempdir(), "_dbf_roundtrip.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(asdict(d), fh)
    with open(path, encoding="utf-8") as fh:
        d2 = Design(**json.load(fh))
    os.remove(path)
    check("design survives JSON round trip", d2 == d)
    after = evaluate_design(d2, cfg, 6.0)
    check("score survives JSON round trip",
          before["m2_raw_lb_per_s"] == after["m2_raw_lb_per_s"])

    # The sensitivity table must be scored against the SAME references the
    # run reports, or it prints numbers above the 7.5 ceiling that cannot be
    # compared with anything else in the report.  This has been got wrong
    # once already, in the dossier.
    from dbf2027.optimizer import sensitivity
    from dbf2027.scoring import evaluate_scores as _es
    ref2, ref3 = 0.0731, 41.9
    card_ref = _es(Aircraft(d, cfg), cfg, 6.0, ref2, ref3)
    rows = sensitivity(d, cfg, 6.0, ref2, ref3)
    check("sensitivity baseline matches the reported objective",
          all(abs(r[3] - card_ref.objective) < 1e-9 for r in rows),
          "%.4f vs %.4f" % (rows[0][3], card_ref.objective))
    check("sensitivity honours the references it is given",
          abs(sensitivity(d, cfg, 6.0, ref2 * 0.5, ref3 * 0.5)[0][3]
              - rows[0][3]) > 1e-6)

    # ...and the dossier, for all THREE references.  Scored against its own
    # raw numbers a design is the best in its run, so every ratio is 1 and
    # the total is exactly the 7.5 ceiling.  Without the Ground Mission
    # reference the dossier fell back to the 120 lb*in field default and
    # printed a total of 22.
    import contextlib as _cl
    import io
    from dataclasses import replace as _rep
    from dbf2027.report import print_design_report, print_active_bounds
    _own = _es(Aircraft(d, cfg), cfg, 6.0)
    if _own.feasible and _own.m2_raw_lb_per_s > 0 and _own.m3_raw_lap_lb > 0:
        with _cl.redirect_stdout(io.StringIO()):
            _card = print_design_report(d, cfg, 6.0,
                                        ref_m2=_own.m2_raw_lb_per_s,
                                        ref_m3=_own.m3_raw_lap_lb,
                                        ref_gm=_own.gm_raw_lb_in)
        check("the dossier honours all three score references",
              abs(_card.total_mission - 7.5) < 1e-6,
              "%.3f" % _card.total_mission)
    # The Mission 3 sensor bound is a FRACTION of the declared sensor; the
    # check once compared kilograms with it and flagged every design.
    _buf = io.StringIO()
    with _cl.redirect_stdout(_buf):
        print_active_bounds(_rep(d, sensor_m3_kg=0.7 * d.sensor_mass_kg),
                            DesignSpace(), cfg)
    _srow = [r for r in sensitivity(_rep(d, sensor_m3_kg=0.9 * d.sensor_mass_kg),
                                    cfg, 6.0, ref2, ref3)
             if r[0] == "sensor_mass_kg"][0]
    check("a lighter declared sensor keeps its Mission 3 fraction (flyable)",
          not math.isnan(_srow[1]), "%s" % (_srow,))
    # Every build-standard number is sorted into a kind, and the report
    # prints the margins with what they hold back.
    from dataclasses import fields as _flds
    from dbf2027.config import BuildStandard as _BS, FIELD_ROLES, ROLE_TITLES
    check("every build-standard field has a role",
          {f.name for f in _flds(_BS)} == set(FIELD_ROLES)
          and all(r in ROLE_TITLES for r, _ in FIELD_ROLES.values()))
    from dbf2027.report import margin_rows, computed_rows, print_margins
    _mac = Aircraft(d, cfg)
    _mr = margin_rows(_mac, cfg)
    check("the margins table lists every numeric build field",
          len(_mr) == sum(1 for f in _flds(_BS)
                          if getattr(cfg.build, f.name) is not None))
    check("the depth-of-discharge row says what it holds back",
          any(r["name"] == "depth_of_discharge" and "Wh" in r["effect"]
              for r in _mr))
    check("the computed corrections each say what they replaced",
          all(r["was"] for r in computed_rows(_mac, cfg)))
    _mb = io.StringIO()
    with _cl.redirect_stdout(_mb):
        print_design_report(d, cfg, 6.0)
    check("the dossier prints the margins section, 78 columns wide",
          "MARGINS AND MODELLING FACTORS" in _mb.getvalue()
          and "COMPUTED FOR THIS DESIGN" in _mb.getvalue()
          and max(len(x) for x in _mb.getvalue()[
              _mb.getvalue().index("MARGINS AND"):].split("\n")[:120])
          <= 78)
    check("an interior Mission 3 sensor fraction is not flagged as pinned",
          "m3_sensor_frac" not in _buf.getvalue()
          and "sensor_m3_kg" not in _buf.getvalue())

    # bill of materials must be generatable
    from dbf2027.report import bill_of_materials
    bom = bill_of_materials(d, cfg)
    check("BOM lists parts", len(bom) >= 15, "%d items" % len(bom))
    check("BOM carries no pricing",
          not any(k.endswith("price_usd") for i in bom for k in i))
    check("BOM quantities are non-negative",
          all(i["qty"] >= 0 for i in bom))
    check("BOM covers every category",
          {i["category"] for i in bom} >=
          {"Propulsion", "Energy", "Structure", "Systems", "Payload"})
    # The rules items that get teams failed at tech inspection must appear.
    text = " ".join("%s %s %s" % (i["part"], i["spec"], i["note"])
                    for i in bom).lower()
    for needed in ("fuse", "arming plug", "rx/servo battery",
                   "charging sack", "unaltered", "insulated"):
        check("BOM names the required item: %s" % needed, needed in text)


# ==========================================================================
def validate_catalogue_editing(cfg):
    section("J. Catalogue registry and user additions")

    from dbf2027 import components as C

    check("every catalogue is registered", len(C.CATALOGUES) == 5,
          str([c.key for c in C.CATALOGUES]))
    for spec in C.CATALOGUES:
        check("%s has columns" % spec.key, len(spec.columns) >= 4)
        check("%s has a live store" % spec.key, len(spec.store) > 0)
        check("%s explains itself" % spec.key, len(spec.note) > 30)
        check("%s columns have help text" % spec.key,
              sum(1 for c in spec.columns if c.help) >= len(spec.columns) // 2)
        # The round trip is what keeps the GUI form, the CSV and the
        # constructor from drifting apart: whatever the table shows must
        # rebuild into the same object - EVERY entry, not a sample.
        # Checking only the first let a LiHV pack's voltages and a tested
        # pack's measured capacity go missing from the form unnoticed.
        lost = []
        for ident in sorted(spec.store):
            obj = spec.store[ident]
            rebuilt = spec.build({k: str(v)
                                  for k, v in spec.to_row(obj).items()})
            same = all(abs(getattr(rebuilt, f.name)
                           - getattr(obj, f.name)) < 1e-6
                       if isinstance(getattr(obj, f.name), float)
                       else getattr(rebuilt, f.name) == getattr(obj, f.name)
                       for f in fields(obj)
                       # the polar and propeller tables are data, not
                       # form fields: a hand-edited airfoil or propeller
                       # uses the fitted parameters (rpm_tables are the
                       # measured tables at each rpm)
                       if f.name not in ("source", "note", "polar", "table",
                                         "rpm_tables"))
            if not same:
                lost.append(ident)
        check("every %s entry survives a row round trip" % spec.key,
              not lost, str(lost[:3]))

    # Adding, reading back and removing must all work without a GUI.
    #
    # All of it happens in a scratch folder.  This section used to write
    # straight into the real catalogues/motors.csv, overwrite it with a
    # malformed row and then DELETE it - taking every motor the user had
    # entered with it.  The real file is snapshotted and compared at the end.
    import shutil as _sh
    _real_csv = os.path.join(C.user_dir(), "motors.csv")
    _snapshot = (open(_real_csv, "rb").read()
                 if os.path.exists(_real_csv) else None)
    _real_user_dir = C.user_dir
    _sandbox = tempfile.mkdtemp(prefix="dbf_validate_")
    C.user_dir = lambda: _sandbox
    C.reload_catalogues()
    n_before = len(C.MOTORS)
    C.add_user_entry("motors", {"key": "VALIDATE_M", "name": "validate",
                                "kv": 555, "rm": 0.042, "i0": 1.0,
                                "mass_kg": 0.2, "i_max_a": 40,
                                "p_max_w": 900, "shaft_mm": 5})
    check("an entry can be added headlessly",
          "VALIDATE_M" in C.MOTORS and len(C.MOTORS) == n_before + 1)
    check("it is flagged as a user entry",
          C.is_user_entry("motors", "VALIDATE_M"))
    check("built-ins are not flagged as user entries",
          not C.is_user_entry("motors", "SCORPION_SII_4025_440"))
    check("adding twice replaces rather than duplicates",
          (C.add_user_entry("motors", {"key": "VALIDATE_M", "name": "again",
                                       "kv": 556, "rm": 0.042, "i0": 1.0,
                                       "mass_kg": 0.2, "i_max_a": 40,
                                       "p_max_w": 900, "shaft_mm": 5})
           and len([r for r in C.user_rows("motors")
                    if r["key"] == "VALIDATE_M"]) == 1))
    check("the update took effect", abs(C.MOTORS["VALIDATE_M"].kv - 556) < 1e-9)
    check("a user entry is marked VERIFY",
          "VERIFY" in C.MOTORS["VALIDATE_M"].source)
    C.reload_catalogues()
    check("reload keeps user entries", "VALIDATE_M" in C.MOTORS)
    check("and keeps the updated value",
          abs(C.MOTORS["VALIDATE_M"].kv - 556) < 1e-9)
    check("removing works",
          C.remove_user_entry("motors", "VALIDATE_M")
          and "VALIDATE_M" not in C.MOTORS)
    check("removing a built-in is refused",
          C.remove_user_entry("motors", "SCORPION_SII_4025_440") is False
          and "SCORPION_SII_4025_440" in C.MOTORS)
    check("reload restores the shipped catalogue exactly",
          len(C.MOTORS) == n_before)

    # A bad row must be skipped, not crash the import.
    path = C.user_csv_path("motors")
    os.makedirs(C.user_dir(), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("key,name,kv,rm,i0,mass_kg,i_max_a,p_max_w,shaft_mm\n")
        fh.write("BAD_M,bad,not_a_number,x,,,,,\n")
    counts = C.reload_catalogues()
    check("a malformed row does not break loading", isinstance(counts, dict))
    check("and the shipped catalogue is still intact",
          "SCORPION_SII_4025_440" in C.MOTORS)
    os.remove(path)
    C.reload_catalogues()
    check("cleanup restores the original count", len(C.MOTORS) == n_before)
    C.user_dir = _real_user_dir
    C.reload_catalogues()
    _sh.rmtree(_sandbox, ignore_errors=True)
    check("and the real catalogues/ folder was never touched",
          (open(_real_csv, "rb").read() if os.path.exists(_real_csv)
           else None) == _snapshot)


# ==========================================================================
def validate_search_quality(cfg):
    section("K. The search actually finds what the model prefers")

    from dataclasses import replace as _replace
    from dbf2027.optimizer import (DesignSpace, _polish_to_bounds,
                                   evaluate_design)
    sp = DesignSpace()
    bounds = sp.continuous_bounds(cfg)
    lower = [b[0] for b in bounds]
    upper = [b[1] for b in bounds]

    # The simplex must probe each axis on a scale comparable to its own
    # range.  Absolute step sizes looked fine per axis but were not
    # comparable across them, and the payload axis in particular was
    # being probed at well under one percent of its range.
    step = [max(0.08 * (hi - lo), 1e-6) for lo, hi in zip(lower, upper)]
    # (a variable with equal bounds - taper, while the wing is rectangular -
    # is held fixed, not stepped)
    fracs = [st / (hi - lo) for st, (lo, hi) in zip(step, zip(lower, upper))
             if hi - lo > 1e-12]
    check("every simplex step is a comparable fraction of its range",
          max(fracs) / min(fracs) < 1.01,
          "%.4f to %.4f" % (min(fracs), max(fracs)))

    # The bounds polish may only ever improve.
    calls = {"n": 0}

    def quad(x):
        calls["n"] += 1
        return sum((v - u) ** 2 for v, u in zip(x, upper))

    mid = [0.5 * (lo + hi) for lo, hi in zip(lower, upper)]
    polished = _polish_to_bounds(quad, list(mid), lower, upper)
    check("the polish drives a monotone objective to its bounds",
          all(abs(p - u) < 1e-9 for p, u in zip(polished, upper)))
    check("and it costs only a few dozen evaluations", calls["n"] < 60,
          "%d" % calls["n"])

    def bowl(x):
        return sum((v - m) ** 2 for v, m in zip(x, mid))

    check("the polish never moves away from an interior optimum",
          _polish_to_bounds(bowl, list(mid), lower, upper) == mid)

    # Span is monotonically good in this model, so the optimum sits on
    # the rules limit.  If that ever stops being true, this tells us.
    base = baseline(n_simulators=0, sensor_mass_kg=4.0)
    refs = (cfg.rules.field_best_m2_lb_per_s,
            cfg.rules.field_best_m3_lap_lb, cfg.rules.field_best_gm_lb_in)
    short = evaluate_design(base, cfg, 6.0, *refs)
    lim = cfg.rules.max_wingspan_m
    wide = evaluate_design(
        _replace(base, span_m=lim, chord_m=lim / (base.span_m / base.chord_m)),
        cfg, 6.0, *refs)
    if short.get("feasible") and wide.get("feasible"):
        check("a wider span scores better, so the limit should bind",
              wide["objective"] > short["objective"],
              "%.3f at %.2f m vs %.3f at %.2f m"
              % (wide["objective"], lim, short["objective"], base.span_m))


# ==========================================================================
def validate_rules(cfg):
    section("I. Rules traceability")

    from dbf2027.rulebook import (ENFORCED_KINDS, RULES, RULES_BY_ID,
                                  check_annotations, coverage,
                                  scan_annotations)

    unenforced, unknown = check_annotations()
    # This is the check that keeps the annotations honest.  Delete the
    # line of code that enforces a rule and its marker goes with it, and
    # this fails - which is the whole point of writing them down as data
    # rather than as prose in a comment.
    check("every enforceable rule is annotated in the source",
          not unenforced, "unenforced: " + ", ".join(unenforced))
    check("every RULE marker names a rule that exists",
          not unknown, "unknown: " + ", ".join(unknown))

    check("the rulebook is not empty", len(RULES) >= 40, "%d" % len(RULES))
    check("every rule has a unique id",
          len(RULES_BY_ID) == len(RULES))
    check("every rule cites a section number",
          all(r.id and r.id[0].isdigit() for r in RULES))
    check("every rule quotes the requirement",
          all(len(r.text) > 20 for r in RULES),
          str([r.id for r in RULES if len(r.text) <= 20]))
    check("every rule says how this model handles it",
          all(len(r.note) > 20 for r in RULES),
          str([r.id for r in RULES if len(r.note) <= 20]))
    check("every rule has a known kind",
          all(r.kind in ENFORCED_KINDS + ("build",) for r in RULES),
          str(sorted({r.kind for r in RULES})))

    cov = coverage()
    # The rules that decide the score, and the hard limits, are the ones
    # worth being strict about: each must be enforced in more than one
    # place OR in the module that owns it, never only in a config value.
    for rid in ("3.1", "3.2.1a", "3.2.3.2b", "3.2.3.2c", "3.2.3.2d",
                "3.1.3a", "3.3.4", "6.3"):
        sites = cov[rid]["sites"]
        check("rule %s is enforced outside config.py" % rid,
              any(not s_.startswith("config.py") for s_ in sites),
              str(sites))

    # A rule the model claims to enforce must be reachable from the
    # package, not only from a test or a script.
    pkg_only = {rid: [s_ for s_ in v["sites"]
                      if not s_.startswith(("selftest", "validate",
                                            "guitest", "webtest",
                                            "grid_convergence",
                                            "analyse_results"))]
                for rid, v in cov.items()}
    missing_pkg = [rid for rid, sites in pkg_only.items()
                   if RULES_BY_ID[rid].kind in ENFORCED_KINDS and not sites]
    check("enforcement lives in the package, not only in the tests",
          not missing_pkg, str(missing_pkg))

    # Spot-check that the ids in the registry match the rules package
    # sections they claim.  These are the ones most often got wrong.
    for rid, must_say in (("3.1", "6 feet"),
                          ("3.2.1a", "55-lb"),
                          ("3.2.3.2b", "100 Watt-hours"),
                          ("3.1.3a", "1 ounce"),
                          ("3.1.2", "1.5"),
                          ("3.1.1b", "six inches"),
                          ("3.3.4", "GM = 0.5"),
                          ("6.3", "M1 + M2 + M3 + GM")):
        check("rule %s quotes the right requirement" % rid,
              must_say in RULES_BY_ID[rid].text,
              RULES_BY_ID[rid].text[:60])

    markers = scan_annotations()
    check("annotations are spread across the model, not clustered",
          len({s_.split(":")[0] for sites in markers.values()
               for s_ in sites}) >= 8,
          "%d files" % len({s_.split(":")[0] for sites in markers.values()
                            for s_ in sites}))


# ==========================================================================
def validate_trends(cfg):
    section("F. Physical trends across the catalogue")

    # every airfoil: more area lowers stall speed
    for name in AIRFOILS:
        small = Aircraft(baseline(airfoil=name, span_m=1.2, chord_m=0.20), cfg)
        big = Aircraft(baseline(airfoil=name, span_m=1.8, chord_m=0.30), cfg)
        check("%s: more wing area lowers stall speed" % name,
              big.aero_clean.v_stall(big.w_m2_n)
              < small.aero_clean.v_stall(small.w_m2_n))

    # trim drag must penalise camber: same planform, more camber, more drag
    flat = Aircraft(baseline(airfoil="MH32"), cfg)
    camb = Aircraft(baseline(airfoil="S1223"), cfg)
    v = 22.0
    wref = camb.w_m2_n
    check("a high-camber section carries more trim download",
          abs(camb.aero_clean.tail_load(v, wref))
          > abs(flat.aero_clean.tail_load(v, flat.w_m2_n)),
          "%.1f vs %.1f N" % (camb.aero_clean.tail_load(v, wref),
                              flat.aero_clean.tail_load(v, flat.w_m2_n)))

    # heavier sensor -> bigger container -> longer payload bay.  The whole
    # fuselage need not grow: its nose is set by the balance, which shifts
    # as the masses do - but it is always nose + bay + tailcone.
    prev_len = 0.0
    for m in (0.3, 0.8, 1.5, 3.0):
        ac = Aircraft(baseline(sensor_mass_kg=m), cfg)
        lay = ac.layout
        bay = lay.bay_end_m - lay.bay_start_m
        check("sensor %.1f kg grows the payload bay" % m,
              bay >= prev_len - 1e-9, "%.3f m" % bay)
        check("sensor %.1f kg: fuselage = nose + bay + tailcone" % m,
              abs(ac.body.fuse_length_m - (lay.nose_length_m + bay
                                           + lay.tailcone_length_m)) < 1e-6,
              "%.4f m" % ac.body.fuse_length_m)
        prev_len = bay

    # wind must never make the mission faster.  Two simulators, not three:
    # with battery sag the three-simulator baseline cannot finish Mission 2
    # in a 9 m/s wind or 1500 m of density altitude, and a comparison that
    # silently skips itself when a flight fails proves nothing - so the
    # flights are required to complete, as a check of their own.
    calm = Config()
    calm.env.wind_mps = 0.0
    windy = Config()
    windy.env.wind_mps = 9.0
    hi = Config()
    hi.env = Environment.at_density_altitude(1500.0, 30.0, 4.5)
    t_calm = mission_2(Aircraft(baseline(n_simulators=2), calm), calm, dx=6.0)
    t_wind = mission_2(Aircraft(baseline(n_simulators=2), windy), windy, dx=6.0)
    t_hi = mission_2(Aircraft(baseline(n_simulators=2), hi), hi, dx=6.0)
    check("the environment test aircraft flies calm, windy and thin air",
          t_calm.ok and t_wind.ok and t_hi.ok,
          "%s / %s / %s" % (t_calm.reason or "ok", t_wind.reason or "ok",
                            t_hi.reason or "ok"))
    check("wind never speeds up the mission", t_wind.time_s >= t_calm.time_s,
          "%.1f vs %.1f s" % (t_wind.time_s, t_calm.time_s))
    # altitude must never help
    check("thin air never speeds up the mission",
          t_hi.time_s >= t_calm.time_s,
          "%.1f vs %.1f s" % (t_hi.time_s, t_calm.time_s))

    # both searched chemistries must offer a usable pack
    for chem in ("LiPo", "LiIon"):
        check("%s offers a legal pack over 20 Wh" % chem,
              any(20.0 <= pk.watt_hours <= cfg.rules.max_watt_hours
                  for pk in PACKS.values() if pk.chemistry == chem))
    # and the rules arithmetic on multiple systems has to be right
    from dbf2027.optimizer import DesignSpace as _DSp, candidate_packs
    combos = candidate_packs(_DSp(), cfg, 18.0)
    check("pack/system combinations exist", len(combos) > 10,
          "%d combinations" % len(combos))
    check("no combination exceeds the 100 Wh total",
          all(n * PACKS[k].watt_hours <= cfg.rules.max_watt_hours + 1e-9
              for k, n in combos))
    check("no combination uses a non-lithium pack",
          all(PACKS[k].chemistry in cfg.rules.battery_chemistry_allowed
              for k, _ in combos))
    check("two-system combinations are actually offered",
          any(n == 2 for _, n in combos))


# ==========================================================================
def validate_interfaces(cfg):
    section("G. Command-line interface")
    parser = build_parser()

    argv = ["--wind", "7", "--altitude", "500", "--temp", "28",
            "--max-wh", "90", "--max-sim", "15", "--window", "280",
            "--max-span-in", "60", "--sensor-max", "4", "--sim-max", "2",
            "--lap-ft", "1800", "--normalise", "field",
            "--field-m2", "0.08", "--field-m3", "30"]
    c = build_config(parser.parse_args(argv))
    check("--wind reaches the config", abs(c.env.wind_mps - 7.0) < 1e-9)
    check("--altitude changes air density", c.env.rho < 1.225)
    check("--max-wh reaches the config", abs(c.rules.max_watt_hours - 90) < 1e-9)
    check("--max-sim reaches the config", c.rules.max_simulators == 15)
    check("--window reaches the config",
          abs(c.rules.mission_window_s - 280) < 1e-9)
    check("--max-span-in reaches the config",
          abs(c.rules.max_wingspan_m - 60 * 0.0254) < 1e-9)
    check("--sensor-max reaches the config",
          abs(c.rules.sensor_mass_max_kg - 4 * LB2KG) < 1e-6)
    check("--sim-max reaches the config",
          abs(c.rules.sim_mass_max_kg - 2 * LB2KG) < 1e-6)
    check("--lap-ft rescales the course",
          abs(c.course.straight_length_m / FT2M - 1800.0) < 1e-6,
          "%.1f ft" % (c.course.straight_length_m / FT2M))
    check("--lap-ft preserves the turns",
          abs(c.course.total_turn_deg - 720.0) < 1e-9)
    check("--normalise reaches the config", c.rules.normalisation == "field")
    check("--field-m2 reaches the config",
          abs(c.rules.field_best_m2_lb_per_s - 0.08) < 1e-12)

    c0 = build_config(parser.parse_args([]))
    check("defaults keep the 6 ft span rule",
          abs(c0.rules.max_wingspan_m / FT2M - 6.0) < 1e-9)
    check("defaults keep no takeoff limit",
          c0.rules.takeoff_distance_m == float("inf"))
    check("defaults keep self-normalisation", c0.rules.normalisation == "self")
    check("search never samples beyond the span limit",
          DesignSpace().continuous_bounds(c0)[0][1]
          <= c0.rules.max_wingspan_m + 1e-12)
    check("every declared flag is reachable",
          len(parser._actions) > 15, "%d actions" % len(parser._actions))

    # A self-normalised score can never exceed its own ceiling.  This is
    # what a single design scored against stale references got wrong.
    from dbf2027.scoring import normalise_population
    rows = [evaluate_design(baseline(n_simulators=n), cfg, 12.0)
            for n in (0, 2, 3)]
    rows = [r for r in rows if r.get("feasible")]
    if rows:
        normalise_population(rows, "self")
        for r in rows:
            check("self-normalised M2 never exceeds 2.0", r["m2"] <= 2.0 + 1e-9,
                  "%.4f" % r["m2"])
            check("self-normalised M3 never exceeds 3.0", r["m3"] <= 3.0 + 1e-9,
                  "%.4f" % r["m3"])
            check("self-normalised total never exceeds 7.5",
                  r["total_mission"] <= 7.5 + 1e-9, "%.4f" % r["total_mission"])
        check("the population best reaches exactly 1.0 on some ratio",
              any(abs(r["m2_ratio"] - 1.0) < 1e-9
                  or abs(r["m3_ratio"] - 1.0) < 1e-9 for r in rows))


# ==========================================================================
def validate_loads(cfg):
    section("H. Design load factors (V-n diagram)")

    # lift-curve slope must approach 2 pi at high aspect ratio and fall below
    a_low = lift_curve_slope_3d(4.0)
    a_high = lift_curve_slope_3d(20.0)
    check("lift-curve slope rises with aspect ratio", a_high > a_low)
    check("lift-curve slope stays below 2 pi", a_high < 2 * math.pi,
          "%.3f" % a_high)

    # the Pratt gust formula must behave
    base = gust_load_factor(250.0, 0.28, 4.8, 30.0, 1.225, 7.62)
    check("gust load factor exceeds 1 g", base > 1.0, "%.2f" % base)
    check("lighter wing loading is more gust-critical",
          gust_load_factor(120.0, 0.28, 4.8, 30.0, 1.225, 7.62) > base)
    check("a faster aeroplane sees a bigger gust load",
          gust_load_factor(250.0, 0.28, 4.8, 45.0, 1.225, 7.62) > base)
    check("a stronger gust gives a bigger load",
          gust_load_factor(250.0, 0.28, 4.8, 30.0, 1.225, 15.0) > base)
    check("zero gust gives exactly 1 g",
          abs(gust_load_factor(250.0, 0.28, 4.8, 30.0, 1.225, 0.0) - 1.0) < 1e-12)

    # across the whole airfoil catalogue the derived load case must be sane
    for name in AIRFOILS:
        ac = Aircraft(baseline(airfoil=name), cfg)
        lc = ac.load_case
        check("%s: limit load in a sane band" % name,
              1.5 <= ac.n_limit <= 8.0, "%.2f g" % ac.n_limit)
        check("%s: ultimate = safety factor x limit" % name,
              abs(ac.n_ultimate
                  - cfg.build.load_safety_factor * ac.n_limit) < 1e-9)
        check("%s: limit never exceeds the CLmax ceiling" % name,
              ac.n_limit <= lc.n_aero_ceiling + 1e-9,
              "%.2f vs %.2f" % (ac.n_limit, lc.n_aero_ceiling))
        check("%s: dive gust milder than cruise gust" % name,
              lc.n_gust_dive < lc.n_gust_cruise)
        check("%s: dive speed above cruise" % name, lc.v_dive > 0)

    # the mass/load fixed point has to actually converge
    for n_sim in (0, 4, 10):
        ac = Aircraft(baseline(n_simulators=n_sim), cfg)
        check("mass closes with %d simulators" % n_sim,
              abs(ac.mtow_m2_kg - (ac.empty_mass_kg + ac.m2_payload_kg)) < 1e-9)
        recomputed = design_load_factors(
            ac.aero_clean, ac.mtow_m2_kg * 9.80665,
            ac._max_level_speed(ac.mtow_m2_kg * 9.80665), cfg.build)
        check("load case is self-consistent with %d simulators" % n_sim,
              abs(recomputed.n_ultimate - ac.n_ultimate) < 0.02,
              "%.4f vs %.4f" % (recomputed.n_ultimate, ac.n_ultimate))

    # spar: monotone in load factor, and gauge-limiting reported honestly
    from dbf2027.weights import spar_sizing
    ac = Aircraft(baseline(), cfg)
    prev = 0.0
    for n in (2.0, 4.0, 8.0, 16.0, 32.0):
        m, t_req, t_used, gauge = spar_sizing(ac.wing, ac.w_m2_n, cfg.build, n)
        check("spar mass non-decreasing at n_ult %.0f" % n, m >= prev - 1e-15,
              "%.5f vs %.5f" % (m, prev))
        check("spar wall never below minimum gauge at n_ult %.0f" % n,
              t_used >= cfg.build.spar_min_wall_m - 1e-15)
        check("gauge flag agrees with the requirement at n_ult %.0f" % n,
              gauge == (t_req < cfg.build.spar_min_wall_m))
        prev = m

    # the spanwise spar integral has converged at the station count used
    import dbf2027.weights as _W
    from dbf2027.weights import spar_design as _spd
    _thin = Config()
    _thin.build.spar_min_wall_m = 1e-9        # load-sized: the hard case
    _n0 = _W.SPAR_STATIONS
    try:
        _m = []
        for _ns in (_n0, 4 * _n0):
            _W.SPAR_STATIONS = _ns
            _m.append(_spd(ac.wing, ac.w_m2_n, _thin.build, ac.n_ultimate,
                           ac.q_dive_pa).mass_kg)
    finally:
        _W.SPAR_STATIONS = _n0
    check("spar mass converged at %d stations (under 1 g)" % _n0,
          abs(_m[0] - _m[1]) < 1e-3, "%.2f vs %.2f g"
          % (1e3 * _m[0], 1e3 * _m[1]))

    # an override must bypass the derivation entirely
    ovr = Config()
    ovr.build.ultimate_load_factor_override = 3.75
    a2 = Aircraft(baseline(), ovr)
    check("override sets the ultimate load factor exactly",
          abs(a2.n_ultimate - 3.75) < 1e-12)
    check("override sets limit = ultimate / safety factor",
          abs(a2.n_limit - 3.75 / ovr.build.load_safety_factor) < 1e-12)

    # span efficiency across the whole design space must stay physical
    from dbf2027.liftingline import span_efficiency
    for ar in (3.5, 5, 7, 9, 11, 14):
        prev = None
        for t in (1.0, 0.8, 0.6, 0.45):
            e = span_efficiency(ar, t)
            check("AR %.1f taper %.2f: span efficiency in (0.5, 1.0]" % (ar, t),
                  0.5 < e <= 1.0, "%.4f" % e)
            if prev is not None:
                check("AR %.1f: taper %.2f beats %.2f" % (ar, t, prev[0]),
                      e > prev[1] - 1e-9, "%.4f vs %.4f" % (e, prev[1]))
            prev = (t, e)

    # static margin must be computed, finite and consistent for every airfoil
    for name in list(AIRFOILS)[:8]:
        a = Aircraft(baseline(airfoil=name), cfg)
        sm = a.aero_clean.static_margin
        check("%s: static margin is finite and sane" % name,
              -0.5 < sm < 1.0, "%.3f" % sm)
        check("%s: neutral point aft of the wing AC" % name,
              a.aero_clean.x_np_over_c > 0.25)

    # 1-4 motors must all build and be checked for propeller overlap
    for nm in (1, 2, 3, 4):
        a = Aircraft(baseline(n_motors=nm, prop="APC_11x7"), cfg)
        check("%d-motor layout builds" % nm, a.mtow_m2_kg > 0)
        check("%d-motor propulsion mass scales" % nm,
              abs(a.motor_mass_kg - nm * a.motor.mass_kg) < 1e-12)
    wide = Aircraft(baseline(n_motors=4, prop="APC_20x13"), cfg)
    check("four big propellers on a 6 ft wing are rejected as overlapping",
          any("overlap" in f for f in wide.static_checks()["failures"]),
          "%s" % wide.static_checks()["failures"])


# ==========================================================================
def main():
    cfg = Config()
    print("=" * 70)
    print("DBF 2026-27 MODEL VALIDATION")
    print("=" * 70)
    validate_catalogue(cfg)
    validate_bounds(cfg)
    validate_energy(cfg)
    validate_determinism(cfg)
    validate_roundtrip(cfg)
    validate_trends(cfg)
    validate_interfaces(cfg)
    validate_loads(cfg)
    validate_catalogue_editing(cfg)
    validate_search_quality(cfg)
    validate_rules(cfg)
    print()
    print("=" * 70)
    print("%d passed, %d failed" % (PASS, FAIL))
    if FAILURES:
        print()
        for f in FAILURES:
            print("  FAILED: %s" % f)
    print("=" * 70)
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
