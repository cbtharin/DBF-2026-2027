"""Human-readable output: rankings, a full design dossier, and a shopping list.

The bill of materials is the point of the whole exercise - it is what turns
"the optimiser likes AR 7.4" into something somebody can order on Monday.
"""
from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, is_dataclass
from typing import Any, Dict, List, Optional, Sequence

from .aircraft import Aircraft, Design
from .config import Config
from .mission import turn_capability
from .scoring import ScoreCard, evaluate_scores
from .units import (FT2M, IN2M, KG2LB, M2FT, M2IN, MS2MPH, N2LBF)


LINE = "=" * 78
THIN = "-" * 78


# ==========================================================================
# Tables
# ==========================================================================
def print_ranking(rows: Sequence[Dict[str, Any]], n: int = 12,
                  title: str = "RANKING") -> None:
    print()
    print(LINE)
    print(title)
    print(LINE)
    hdr = ("%-3s %-9s %6s %6s %5s %-22s %-10s %4s %5s %5s %6s %6s %6s %6s"
           % ("#", "airfoil", "span", "chord", "AR", "motor", "prop",
              "nM", "Wh", "M2t", "M3lap", "M2", "M3", "GM"))
    print(hdr)
    print(THIN)
    for i, r in enumerate(rows[:n], 1):
        d: Design = r["design"]
        print("%-3d %-9s %6.2f %6.3f %5.2f %-22s %-10s %4d %5.0f %5.1f %6d %6.3f %6.3f %6.3f"
              % (i, d.airfoil, d.span_m, d.chord_m, r["AR"],
                 r["motor"][:22], d.prop.replace("APC_", ""), d.n_motors,
                 r["total_Wh"], r.get("m2_time_s", 0.0), r.get("m3_laps", 0),
                 r["m2"], r["m3"], r.get("gm", 0.0)))
    print(THIN)
    print("M2t = Mission 2 elapsed time (s);  M3lap = laps flown in Mission 3")
    print("Scores are normalised against the best design in this run.")


def print_pareto(rows: Sequence[Dict[str, Any]]) -> None:
    print()
    print(LINE)
    print("PARETO FRONT  (raw Mission 2 weight/time  vs  raw Mission 3 laps x sensor)")
    print(LINE)
    print("%-3s %8s %9s  %6s %6s %5s  %-20s %-9s %4s %5s"
          % ("#", "M2 lb/s", "M3 lap*lb", "span", "chord", "sens", "motor",
             "prop", "nSim", "GMlbin"))
    print(THIN)
    for i, r in enumerate(rows, 1):
        d: Design = r["design"]
        print("%-3d %8.4f %9.2f  %6.2f %6.3f %5.2f  %-20s %-9s %4d %6.1f"
              % (i, r["m2_raw_lb_per_s"], r["m3_raw_lap_lb"], d.span_m,
                 d.chord_m, d.sensor_mass_kg * KG2LB, r["motor"][:20],
                 d.prop.replace("APC_", ""), d.n_simulators,
                 r.get("gm_raw_lb_in", 0.0)))
    print(THIN)
    print("Every row here is a legitimate answer.  Which one you build depends")
    print("entirely on how strong you think the field is in each mission.")


# ==========================================================================
# Dossier
# ==========================================================================
def print_design_report(design: Design, cfg: Config, dx: float = 2.0,
                        title: str = "OPTIMAL DESIGN",
                        ref_m2: Optional[float] = None,
                        ref_m3: Optional[float] = None) -> ScoreCard:
    """Full dossier for one design.

    `ref_m2` / `ref_m3` must be the normalisers the search actually used.
    Without them a self-normalised run would rescore this design against
    the field-mode defaults and print a total above its own 6.0 ceiling.
    """
    ac = Aircraft(design, cfg)
    card = evaluate_scores(ac, cfg, dx, ref_m2, ref_m3)
    s = ac.summary()
    checks = ac.static_checks()
    rc = checks["ratings"]

    print()
    print(LINE)
    print(title)
    print(LINE)

    # ---------------- geometry ----------------
    print("\nWING")
    print("  airfoil                %s" % ac.airfoil.name)
    print("    %s" % ac.airfoil.note)
    print("  span                   %.3f m   (%.1f in)"
          % (design.span_m, design.span_m * M2IN))
    print("  chord                  %.3f m   (%.2f in)"
          % (design.chord_m, design.chord_m * M2IN))
    print("  area                   %.4f m2  (%.0f in2)"
          % (s["area_m2"], s["area_m2"] / (IN2M ** 2)))
    print("  aspect ratio           %.2f" % s["AR"])
    print("  taper ratio            %.2f  (root %.0f mm, tip %.0f mm)"
          % (design.taper, 1000 * ac.wing.root_chord_m,
             1000 * ac.wing.root_chord_m * design.taper))
    print("  span efficiency e      %.3f  (lifting line)" % s["oswald_e"])
    print("  CL max (3D, built)     %.2f" % s["CLmax_3d"])
    print("  CD0 at 20 m/s          %.4f" % s["CD0_at_20ms"])
    print("  wing loading           %.1f N/m2  (%.2f oz/ft2)"
          % (s["wing_loading_Nm2"],
             s["wing_loading_Nm2"] * N2LBF * 16.0 / (M2FT ** 2)))

    b = ac.body
    print("\nFUSELAGE AND TAILS")
    print("  fuselage L x W x H     %.3f x %.3f x %.3f m"
          % (b.fuse_length_m, b.fuse_width_m, b.fuse_height_m))
    print("  payload bay holds      %d boxes (1 sensor container + %d simulators)"
          % (1 + design.n_simulators, design.n_simulators))
    print("  tail arm               %.3f m" % b.tail_arm_m)
    print("  horizontal tail area   %.4f m2 (%.0f in2)"
          % (b.s_horiz_m2, b.s_horiz_m2 / (IN2M ** 2)))
    print("  vertical tail area     %.4f m2 (%.0f in2)"
          % (b.s_vert_m2, b.s_vert_m2 / (IN2M ** 2)))

    # ---------------- weights ----------------
    print("\nWEIGHT  (W = m g)")
    for name, kg in ac.mass_breakdown.sorted_items():
        print("  %-24s %7.3f kg  (%6.2f lb)" % (name, kg, kg * KG2LB))
    print("  %-24s %7.3f kg  (%6.2f lb)  <- empty"
          % ("EMPTY", ac.empty_mass_kg, ac.empty_mass_kg * KG2LB))
    print("  %-24s %7.3f kg  (%6.2f lb)"
          % ("M2 payload (scored)", ac.m2_scored_kg, ac.m2_scored_kg * KG2LB))
    print("  %-24s %7.3f kg  (%6.2f lb)  <- W for Mission 2"
          % ("M2 gross", ac.mtow_m2_kg, ac.mtow_m2_kg * KG2LB))
    print("  %-24s %7.3f kg  (%6.2f lb)  <- W for Mission 3"
          % ("M3 gross", ac.mtow_m3_kg, ac.mtow_m3_kg * KG2LB))

    print("\nLONGITUDINAL STABILITY  (computed, not assumed)")
    a = ac.aero_clean
    print("  wing lift slope        %.3f /rad" % a.a_wing)
    print("  tail lift slope        %.3f /rad  (eta_t %.2f)"
          % (a.a_tail, a.eta_tail))
    print("  downwash de/dalpha     %.3f" % a.downwash)
    print("  tail volume ratio      %.3f" % a.v_h)
    print("  neutral point          %.3f c" % a.x_np_over_c)
    print("  CG                     %.3f c   (an input, not an output)"
          % a.x_cg_over_c)
    print("  STATIC MARGIN          %.3f c   %s"
          % (a.static_margin,
             "stable" if a.static_margin > 0 else "UNSTABLE"))
    print("    a positive margin means it will fly hands-off; verify in AVL")

    # ---------------- structure ----------------
    # The team asked for a verdict here, not a stress report: the load
    # factors and spar wall are still computed (the wing mass depends on
    # them, so they cannot be skipped), but what gets printed is whether
    # the wing holds together and what drove the case.
    lc = ac.load_case
    print("\nSTRUCTURE  (will it break?)")
    if ac.structure_ok:
        print("  VERDICT                NO - the spar carries the design")
        print("                         load with the working stress inside")
        print("                         the allowable, and the tube it needs")
        print("                         fits inside the aerofoil.")
    else:
        print("  VERDICT                YES, IT BREAKS - the spar this wing")
        print("                         needs is thicker than the section is")
        print("                         deep.  Shorten the span, deepen the")
        print("                         aerofoil, or carry less.")
    print("  sized for              %.2f g ultimate (%.2f g limit, %s)"
          % (ac.n_ultimate, ac.n_limit, lc.driver))
    if ac.spar_gauge_limited:
        print("  the spar is GAUGE limited: the minimum practical wall, not")
        print("  bending, sets its mass, so it is stronger than it needs to be.")

    # ---------------- propulsion ----------------
    print("\nPROPULSION  (T)")
    print("  motor                  %s  x%d" % (ac.motor.name, design.n_motors))
    print("    Kv %.0f rpm/V, Rm %.3f ohm, I0 %.1f A, %.0f g each"
          % (ac.motor.kv, ac.motor.rm, ac.motor.i0, ac.motor.mass_kg * 1000))
    print("  propeller              %s  (%.1f x %.1f in)"
          % (ac.prop.key, ac.prop.diameter_m * M2IN, ac.prop.pitch_m * M2IN))
    print("  ESC                    %s x%d"
          % (ac.prop_sys.esc.key, design.n_motors))
    print("  propulsion systems     %d  (one pack, fuse and arming plug"
          " each)" % ac.n_systems)
    print("  battery                %s" % ac.pack.name)
    print("    %s, %dS, %.1f V nominal, %.0f mAh, %.0fC, %.0f g"
          % (ac.pack.chemistry, ac.pack.series, ac.pack.v_nominal,
             ac.pack.capacity_ah * 1000, ac.pack.c_rating,
             ac.pack.mass_kg * 1000))
    print("    %.1f Wh per pack x %d = %.1f Wh total  (limit %.0f Wh)"
          % (ac.pack.watt_hours, ac.n_systems, ac.total_watt_hours,
             cfg.rules.max_watt_hours))
    print("    labelled max discharge %.0f A -> %.0f A blade fuse"
          % (ac.pack.i_max_a, ac.prop_sys.fuse_rating_a))
    print("  static thrust T0       %.1f N  (%.2f lbf)"
          % (rc["static_thrust_n"], rc["static_thrust_n"] * N2LBF))
    print("  static T/W  (M2)       %.2f" % s["static_T_over_W_m2"])
    print("  static rpm / tip Mach  %.0f / %.2f" % (rc["rpm"], rc["tip_mach"]))
    print("  worst-case currents    %.0f A per motor, %.0f A per pack,"
          " %.0f A total" % (rc["motor_current_a"], rc["pack_current_a"],
                             rc["total_current_a"]))
    print("    total limit %.0f A, fuse %.0f A -> %s"
          % (cfg.rules.max_total_current_a, rc["fuse_rating_a"],
             "inside both" if (rc["total_current_a"]
                               <= cfg.rules.max_total_current_a
                               and rc["fuse_margin"] >= 1.0) else "OVER"))
    print("  margins  motor I %.2fx  motor P %.2fx  ESC %.2fx  pack C %.2fx"
          "  fuse %.2fx"
          % (rc["motor_current_margin"], rc["motor_power_margin"],
             rc["esc_current_margin"], rc["pack_current_margin"],
             rc["fuse_margin"]))

    # ---------------- L / D / T / W at the design point ----------------
    print("\nFORCE BALANCE AT MISSION 2 CRUISE")
    v = card.r2.v_max if (card.r2 and card.r2.ok and card.r2.v_max > 0) else 25.0
    L = ac.w_m2_n
    D = ac.aero_clean.drag(v, L)
    T = ac.prop_sys.thrust(v, 0.7, 1.0)
    print("  V                      %.1f m/s  (%.0f mph)" % (v, v * MS2MPH))
    print("  L = W                  %.1f N  (%.2f lbf)" % (L, L * N2LBF))
    print("  D                      %.1f N  (%.2f lbf)" % (D, D * N2LBF))
    print("  T available            %.1f N  (%.2f lbf)" % (T, T * N2LBF))
    print("  L/D                    %.2f" % (L / max(D, 1e-9)))
    lt = ac.aero_clean.tail_load(v, L)
    print("  tail load (trim)       %+.1f N (%+.2f lbf), CL_t %+.2f"
          % (lt, lt * N2LBF, ac.aero_clean.tail_cl(v, L)))
    print("    %s; the wing therefore carries %.1f N, not %.1f N"
          % ("download" if lt < 0 else "upload", L - lt, L))
    lt_turn = ac.aero_clean.tail_load(v, ac.n_limit * L)
    print("  tail load at %.1f g       %+.1f N (%+.2f lbf)"
          % (ac.n_limit, lt_turn, lt_turn * N2LBF))
    print("  CL                     %.3f" % ac.aero_clean.cl_required(v, L))

    # ---------------- missions ----------------
    _print_mission("MISSION 1 - flight (no payload, feasibility)", card.r1)
    _print_mission("MISSION 2 - delivery flight", card.r2)
    _print_mission("MISSION 3 - sensor flight", card.r3)

    pg = ac.payload_geom
    print("\nPAYLOAD SELECTION")
    print("  sensor mass            %.3f kg  (%.2f lb)"
          % (design.sensor_mass_kg, design.sensor_mass_kg * KG2LB))
    print("  sensor envelope        %.0f x %.0f x %.0f mm  (%.2f L at %.0f kg/m3)"
          % (1000 * pg.sensor_l, 1000 * pg.sensor_w, 1000 * pg.sensor_h,
             1000 * pg.sensor_volume_m3,
             cfg.rules.sensor_packing_density_kg_m3))
    print("  shipping container     %.0f x %.0f x %.0f mm, %.3f kg (%.2f lb) tare"
          % (1000 * pg.box_l, 1000 * pg.box_w, 1000 * pg.box_h,
             pg.container_tare_kg, pg.container_tare_kg * KG2LB))
    print("    ONE fixed box for every design (assumed, not resized per")
    print("    sensor); %.0f mm foam per face holds %.0f g in a %.0f in drop"
          % (1000 * pg.pad_m, cfg.rules.sensor_survivable_g,
             design.drop_height_in))
    print("    usable volume %.2f L"
          % (1000 * cfg.rules.container_fill_fraction * pg.box_volume_m3))
    dens = ac.sensor_density_kg_m3
    print("  sensor density needed  %.0f kg/m3 to fit that box" % dens)
    if dens <= 1200:
        note = "easy - electronics and foam"
    elif dens <= 2800:
        note = "dense packing, or a little ballast"
    elif dens <= 5000:
        note = "needs real ballast (aluminium/steel shot)"
    elif dens <= 8000:
        note = "steel-shot ballast throughout"
    else:
        note = "lead ballast; check this is still a sensor"
    print("    %s  (water 1000, aluminium 2700, steel 7850, lead 11340)"
          % note)
    print("  simulators flown       %d at %.3f kg  (%.2f lb each)"
          % (design.n_simulators, ac.sim_mass_kg, ac.sim_mass_kg * KG2LB))
    print("    each simulator weighs exactly one loaded container "
          "(rule 3.1.3a),")     # RULE 3.1.3a
    print("    and is one container envelope of extra bay volume and drag")
    print("  Mission 3 sensor       %.3f kg  (%.2f lb)"
          % (ac.sensor_m3_kg, ac.sensor_m3_kg * KG2LB))
    print("    may be under the declared max; tow line %.2f m adds "
          "%.4f m2 of CD*A" % (ac.tow_line_length_m,
                               ac.tow_line_drag_area_m2))
    print("  Ground Mission         %.0f in declared drop, %.0f mm foam"
          % (design.drop_height_in, 1000 * ac.payload_geom.pad_m))
    print("  bay packing            %d rows x %d across x %d high"
          % (ac.bay_rows, ac.bay_cols, ac.bay_layers))
    print("  M2 scored weight       %.3f kg  (%.2f lb)"
          % (ac.m2_scored_kg, ac.m2_scored_kg * KG2LB))
    print("  -> declare at Tech Inspection a maximum of %d simulators"
          % max(design.n_simulators, 0))

    # ---------------- score ----------------
    print("\n" + LINE)
    print("SCORE")
    print(LINE)
    ref2 = (card.m2_raw_lb_per_s / card.m2_ratio) if card.m2_ratio > 0 else 0.0
    ref3 = (card.m3_raw_lap_lb / card.m3_ratio) if card.m3_ratio > 0 else 0.0
    refg = (card.gm_raw_lb_in / card.gm_ratio) if card.gm_ratio > 0 else 0.0
    print("  M2 raw  weight/time            %.4f lb/s" % card.m2_raw_lb_per_s)
    print("  M3 raw  laps x sensor weight   %.2f lap*lb" % card.m3_raw_lap_lb)
    print("  GM raw  weight x drop height   %.1f lb*in" % card.gm_raw_lb_in)
    if card.mode == "self":
        print("  normalised against             the best design in this run")
        print("    reference M2                 %.4f lb/s" % ref2)
        print("    reference M3                 %.2f lap*lb" % ref3)
        print("    reference GM                 %.1f lb*in" % refg)
    else:
        print("  assumed field best  M2         %.4f lb/s" % ref2)
        print("  assumed field best  M3         %.2f lap*lb" % ref3)
        print("  assumed field best  GM         %.1f lb*in" % refg)
    print(THIN)
    cap = "min(%s, 1)" if card.mode == "field" else "%s"
    print("  M1 = %.3f" % card.m1)
    print(("  M2 = 1 + " + cap + " = %.3f%s")
          % ("%.4f/%.4f" % (card.m2_raw_lb_per_s, ref2), card.m2,
             "   [SATURATED]" if (card.mode == "field" and card.m2_saturated)
             else ""))
    print(("  M3 = 2 + " + cap + " = %.3f%s")
          % ("%.2f/%.2f" % (card.m3_raw_lap_lb, ref3), card.m3,
             "   [SATURATED]" if (card.mode == "field" and card.m3_saturated)
             else ""))
    print(("  GM = 0.5 + " + cap + " = %.3f%s")
          % ("%.1f/%.1f" % (card.gm_raw_lb_in, refg), card.gm,
             "   [SATURATED]" if (card.mode == "field" and card.gm_saturated)
             else ""))
    print("  TOTAL MISSION SCORE = %.3f   (M1 + M2 + M3 + GM, ceiling 7.5)"
          % card.total_mission)
    if card.mode == "self":
        print()
        print("  Self-normalised: the best design in this run scores 1.0 on")
        print("  each ratio, so the ceiling is 7.5 and it is reached by")
        print("  whichever candidate leads ALL THREE scored events. These")
        print("  are relative numbers for choosing between your own")
        print("  candidates, not a prediction of your competition score.")
    elif card.m2_saturated or card.m3_saturated or card.gm_saturated:
        print()
        print("  NOTE: this design beats the assumed field best by %s."
              % ", ".join(
                  ([("%.0f%% on M2" % (100 * (card.m2_ratio - 1)))]
                   if card.m2_saturated else [])
                  + ([("%.0f%% on M3" % (100 * (card.m3_ratio - 1)))]
                     if card.m3_saturated else [])
                  + ([("%.0f%% on GM" % (100 * (card.gm_ratio - 1)))]
                     if card.gm_saturated else [])))
        print("  Max(...) in the rules includes you, so the ratio caps at 1.0.")
    if card.reason:
        print("  notes: %s" % card.reason)
    return card


def print_active_bounds(design: Design, space, cfg: Config) -> None:
    """Flag design variables sitting on a search bound.

    A variable pinned to its bound means the bound is the answer, not the
    physics - and most of these bounds are rules numbers you still have to
    verify.  Worth knowing before you believe an optimum.
    """
    b = space.continuous_bounds(cfg)
    checks = [("span_m", design.span_m, b[0]),
              ("aspect_ratio", design.span_m / design.chord_m, b[1]),
              ("sensor_mass_kg", design.sensor_mass_kg, b[2]),
              ("sensor_m3_kg",
               design.sensor_m3_kg or design.sensor_mass_kg, b[3]),
              ("drop_height_in", design.drop_height_in, b[5]),
              ("taper", design.taper, b[4])]
    hits = []
    for name, val, (lo, hi) in checks:
        span = max(hi - lo, 1e-9)
        if val >= hi - 0.01 * span:
            hits.append((name, val, "upper", hi))
        elif val <= lo + 0.01 * span:
            hits.append((name, val, "lower", lo))
    if design.n_simulators >= cfg.rules.max_simulators:
        hits.append(("n_simulators", float(design.n_simulators), "upper",
                     float(cfg.rules.max_simulators)))

    print()
    print(LINE)
    print("ACTIVE CONSTRAINTS")
    print(LINE)
    if not hits:
        print("  None - the optimum is interior, set by the physics.")
        return
    for name, val, side, bound in hits:
        print("  %-16s = %-8.3f  pinned to its %s bound (%.3f)"
              % (name, val, side, bound))
    print()
    print("  These are limits you imposed, not answers the aeroplane gave you.")
    print("  Check each against the actual 2026-27 rules and widen or fix it;")
    print("  the optimum will move.")
    if any(h[0] == "span_m" and h[2] == "upper" for h in hits):
        print()
        print("  SPAN IS AT ITS BOUND.  With no rules limit on size, bigger")
        print("  wings keep winning: area grows as span squared while the")
        print("  spar only grows linearly.  The model has no transport box,")
        print("  hangar door or pit table in it.  Decide a practical span")
        print("  yourself and set --max-span-in, or the optimiser will hand")
        print("  you an aeroplane you cannot get to the field.")


def _print_mission(title: str, r) -> None:
    print("\n" + title)
    if r is None:
        print("  not evaluated")
        return
    if not r.ok:
        print("  FAILED: %s" % (r.reason or "unknown"))
        if r.takeoff_m:
            print("  takeoff roll %.1f m (%.0f ft)" % (r.takeoff_m, r.takeoff_m * M2FT))
        return
    print("  laps                   %d" % r.laps)
    print("  elapsed time           %.1f s" % r.time_s)
    print("  takeoff ground roll    %.1f m  (%.0f ft)   [no rules limit]"
          % (r.takeoff_m, r.takeoff_m * M2FT))
    print("  stall speed            %.1f m/s (%.0f mph)"
          % (r.v_stall, r.v_stall * MS2MPH))
    print("  max level speed        %.1f m/s (%.0f mph)"
          % (r.v_max, r.v_max * MS2MPH))
    print("  corner speed in turns  %.1f m/s (%.0f mph), limited by %s"
          % (r.v_turn, r.v_turn * MS2MPH, r.turn_limited_by))
    if r.lap_times:
        shown = ", ".join("%.1f" % t for t in r.lap_times[:8])
        print("  lap times              %s%s s"
              % (shown, " ..." if len(r.lap_times) > 8 else ""))
    print("  energy used            %.1f Wh  (%.0f%% of the pack still in it)"
          % (r.energy_j / 3600.0, 100.0 * r.energy_margin))
    if r.landing_time_s > 0.0:
        print("  landing                %.1f s, %.0f m rollout  (NOT in the"
              " scored window)" % (r.landing_time_s, r.landing_roll_m))
        print("  total time on the clock %.1f s throttle-up to wheels stopped"
              % r.total_flight_time_s)


# ==========================================================================
# Bill of materials
# ==========================================================================
def bill_of_materials(design: Design, cfg: Config) -> List[Dict[str, Any]]:
    ac = Aircraft(design, cfg)
    b = cfg.build
    items: List[Dict[str, Any]] = []

    def add(category, part, spec, qty, note=""):
        items.append({"category": category, "part": part, "spec": spec,
                      "qty": qty, "note": note})

    # ---- propulsion -------------------------------------------------
    add("Propulsion", ac.motor.name,
        "Kv %.0f, Rm %.3f ohm, I0 %.1f A, %.0f g"
        % (ac.motor.kv, ac.motor.rm, ac.motor.i0, ac.motor.mass_kg * 1000),
        design.n_motors, ac.motor.source)
    add("Propulsion", "Propeller %s" % ac.prop.key,
        "%.1f x %.1f in, %.0f g" % (ac.prop.diameter_m * M2IN,
                                    ac.prop.pitch_m * M2IN,
                                    ac.prop.mass_kg * 1000),
        design.n_motors + 2, "buy spares - you will break them")
    add("Propulsion", "ESC %s" % ac.prop_sys.esc.key,
        "%.0f A continuous, up to %dS" % (ac.prop_sys.esc.i_cont_a,
                                          ac.prop_sys.esc.max_cells),
        design.n_motors, "BEC must be disabled (rule 3.2.3)")
    add("Propulsion", "Motor mount + spinner",
        "for %.0f mm shaft" % ac.motor.shaft_mm, design.n_motors, "")

    # ---- energy -----------------------------------------------------
    i_pack = ac.prop_sys.rating_check()["pack_current_a"]
    add("Energy", ac.pack.name,
        "%s %dS, %.1f V, %.0f mAh, %.0fC, %.1f Wh, %.0f g"
        % (ac.pack.chemistry, ac.pack.series, ac.pack.v_nominal,
           ac.pack.capacity_ah * 1000, ac.pack.c_rating,
           ac.pack.watt_hours, ac.pack.mass_kg * 1000),
        2 * ac.n_systems,
        "COTS, UNALTERED, original label must show chemistry/mAh/V/C. "
        "%d to fly + %d to charge. %s"
        % (ac.n_systems, ac.n_systems, ac.pack.source))
    if ac.n_systems > 1:
        add("Energy", "-> identical packs required", "all %d packs must be "
            "the same manufacturer and part number" % ac.n_systems,
            0, "rule 3.2.3.1; they may NOT be wired in series or parallel")
    add("Energy", "Blade fuse + holder",
        "%.0f A blade fuse, draws %.0f A per pack"
        % (ac.prop_sys.fuse_rating_a, i_pack),
        ac.n_systems,
        "own fuse harness per propulsion system, between battery and "
        "arming plug")
    add("Energy", "Arming plug",
        "externally accessible, rated above %.0f A" % i_pack,
        ac.n_systems, "one per propulsion system")
    add("Energy", "Rx/servo battery",
        "separate pack, %.0f g budget, any chemistry"
        % (b.rx_battery_mass_kg * 1000),
        1, "REQUIRED and independent of the propulsion packs")
    add("Energy", "Main connectors",
        "fully insulated style, XT90 or AS150 for %.0f A" % i_pack,
        4 * ac.n_systems, "no exposed conductors (rule 3.2.3)")
    add("Energy", "LiPo charging sack",
        "commercial, unaltered, manufacturer labelled", 2 * ac.n_systems,
        "batteries live in the sack except in the aeroplane or at tech "
        "inspection; the fuse harness goes in with them")

    # ---- structure --------------------------------------------------
    r_spar = max(0.40 * ac.wing.thickness_root_m, 0.006)
    m_root = (ac.n_ultimate * (ac.w_m2_n / 2.0) * 0.4244
              * (ac.wing.span_m / 2.0))
    t_req = m_root / (b.spar_sigma_allow_pa * math.pi * r_spar * r_spar)
    t_spec = max(t_req, b.spar_min_wall_m)
    add("Structure", "Carbon spar tube",
        "OD %.1f mm, wall %.2f mm, length %.0f mm"
        % (2000 * r_spar, 1000 * ac.spar_t_used_m, 1000 * ac.wing.span_m),
        1,
        "sized for %.0f N*m root bending at n_ult %.2f (%s)"
        % (m_root, ac.n_ultimate, ac.load_case.driver))
    add("Structure", "Wing cores / ribs",
        "%.0f in2 planform, %s section, t/c %.0f%%"
        % (ac.wing.area_m2 / (IN2M ** 2), ac.airfoil.name,
           100 * ac.airfoil.t_over_c),
        1, "")
    add("Structure", "Covering film",
        "%.1f m2 plus 30%% waste" % (2.1 * ac.wing.area_m2), 1, "")
    add("Structure", "Fuselage shell",
        "%.0f x %.0f x %.0f mm, %.2f m2 wetted"
        % (1000 * ac.body.fuse_length_m, 1000 * ac.body.fuse_width_m,
           1000 * ac.body.fuse_height_m, ac.body.fuse_wetted_m2),
        1, "carbon/ply box or moulded")
    add("Structure", "Tail surfaces",
        "H %.0f in2, V %.0f in2, arm %.0f mm"
        % (ac.body.s_horiz_m2 / (IN2M ** 2), ac.body.s_vert_m2 / (IN2M ** 2),
           1000 * ac.body.tail_arm_m),
        1, "")
    if design.n_motors == 1:
        gear_note = ("main gear %.0f mm legs, %.0f mm prop tip clearance"
                     % (1000 * ac.gear_height_m,
                        1000 * ac.static_checks()["prop_ground_clearance_m"]))
    else:
        gear_note = ("main gear %.0f mm legs (props are wing-mounted, "
                     "clearance is not the driver)" % (1000 * ac.gear_height_m))
    add("Structure", "Landing gear", gear_note, 1, "")

    # ---- systems ----------------------------------------------------
    add("Systems", "Servos", "%.0f g class, %d off"
        % (b.servo_mass_kg * 1000, b.n_servos_base), b.n_servos_base,
        "2 aileron, 1 elevator, 1 rudder, 1 payload door, 1 sensor deploy")
    add("Systems", "Receiver + telemetry", "6+ channel", 1,
        "powered from the separate Rx battery, not a BEC")
    add("Systems", "Sensor deploy mechanism",
        "winch/door, %.0f g budget" % (b.deploy_mech_kg * 1000), 1,
        "must fully stow before landing in M3")
    add("Systems", "Sensor light controller",
        "OFF / SOLID ON / FLASHING ON by remote command", 1,
        "M3 requires all three states on command")

    # ---- payload ----------------------------------------------------
    add("Payload", "Sensor article",
        "%.0f g target (%.2f lb)" % (design.sensor_mass_kg * 1000,
                                     design.sensor_mass_kg * KG2LB),
        1, "mass is a DESIGN VARIABLE - it multiplies the M3 score")
    pg = ac.payload_geom
    add("Payload", "Shipping container",
        "%.0f x %.0f x %.0f mm, %.0f mm foam, tare %.0f g"
        % (1000 * pg.box_l, 1000 * pg.box_w, 1000 * pg.box_h,
           1000 * pg.pad_m, pg.container_tare_kg * 1000),
        1, "fixed assumed box; foam sized for the declared %.0f in drop "
        "at %.0f g" % (design.drop_height_in, cfg.rules.sensor_survivable_g))
    add("Payload", "Container simulators",
        "%.0f g each (%.2f lb), same envelope as the container"
        % (ac.sim_mass_kg * 1000, ac.sim_mass_kg * KG2LB),
        design.n_simulators,
        "declare max %d at Tech Inspection" % design.n_simulators)
    return items


def print_bom(design: Design, cfg: Config) -> None:
    items = bill_of_materials(design, cfg)
    print()
    print(LINE)
    print("BILL OF MATERIALS  (parts to build this aeroplane)")
    print(LINE)
    cat = None
    for it in items:
        if it["category"] != cat:
            cat = it["category"]
            print("\n[%s]" % cat)
        qty = ("x%-3d" % it["qty"]) if it["qty"] else "    "
        print("  %-30s %s  %s" % (it["part"][:30], qty, it["spec"]))
        if it["note"]:
            print("  %-30s      %s" % ("", it["note"]))


# ==========================================================================
# Export
# ==========================================================================
def _row_to_flat(row: Dict[str, Any]) -> Dict[str, Any]:
    out = {}
    for k, v in row.items():
        if k == "design":
            d: Design = v
            out.update({"d_" + kk: vv for kk, vv in asdict(d).items()})
        elif isinstance(v, (int, float, str, bool)) or v is None:
            out[k] = v
    return out


def export_csv(rows: Sequence[Dict[str, Any]], path: str) -> None:
    flat = [_row_to_flat(r) for r in rows]
    if not flat:
        return
    keys: List[str] = []
    for f in flat:
        for k in f:
            if k not in keys:
                keys.append(k)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for f in flat:
            w.writerow(f)


def export_json(rows: Sequence[Dict[str, Any]], path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump([_row_to_flat(r) for r in rows], fh, indent=2)


def export_bom_csv(design: Design, cfg: Config, path: str) -> None:
    items = bill_of_materials(design, cfg)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(items[0].keys()))
        w.writeheader()
        w.writerows(items)


def _fmt_obj(v: float) -> str:
    if v is None or v != v or v < -1e6:      # nan or the infeasible sentinel
        return "infeasible"
    return "%10.3f" % v


def print_sensitivity(rows) -> None:
    print()
    print(LINE)
    print("SENSITIVITY  (objective = M2 + M3, one variable at a time)")
    print(LINE)
    print("%-22s %10s %10s %10s" % ("variable", "-10%", "baseline", "+10%"))
    print(THIN)
    for name, lo, hi, base in rows:
        print("%-22s %10s %10.3f %10s"
              % (name, _fmt_obj(lo).strip(), base, _fmt_obj(hi).strip()))
    print(THIN)
    print("Flat rows are dimensions you can trade away for buildability.")
    print("'infeasible' means that step breaks a constraint - usually a rules")
    print("limit or a power-system rating, not the aerodynamics.")
