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
from .aerodynamics import tail_form_factor
from .config import FIELD_ROLES, ROLE_TITLES, Config
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
    hdr = ("%-3s %-9s %6s %6s %5s %-22s %-10s %4s %5s %7s %7s %5s %6s %6s %6s %6s %4s"
           % ("#", "airfoil", "span", "chord", "AR", "motor", "prop",
              "nM", "Wh", "grossKg", "grossLb", "M2t", "M3lap", "M2", "M3",
              "GM", "wx"))
    print(hdr)
    print(THIN)
    for i, r in enumerate(rows[:n], 1):
        d: Design = r["design"]
        gross_kg = r.get("mtow_m2_kg", 0.0)
        print("%-3d %-9s %6.2f %6.3f %5.2f %-22s %-10s %4d %5.0f %7.2f %7.1f %5.1f %6d %6.3f %6.3f %6.3f %4s"
              % (i, d.airfoil, d.span_m, d.chord_m, r["AR"],
                 r["motor"][:22], d.prop.replace("APC_", ""), d.n_motors,
                 r["total_Wh"], gross_kg, gross_kg * KG2LB,
                 r.get("m2_time_s", 0.0), r.get("m3_laps", 0),
                 r["m2"], r["m3"], r.get("gm", 0.0), r.get("wx", "")))
    print(THIN)
    print("M2t = Mission 2 elapsed time (s);  M3lap = laps flown in Mission 3")
    if any("wx" in r for r in rows[:n]):
        print("wx = weather: ok = flies all three presets on its own propeller;"
              " swap = needs a")
        print("     different one on some day (legal, rule 3.2.1e); no = cannot"
              " fly every preset")
    # Display only - the 55 lb limit (rule 3.2.1a) is enforced in
    # aircraft.static_checks, not here.  It is checked at Mission 2 weight,
    # the heaviest of the three: Mission 2 must carry the sensor at its
    # maximum declared weight, and Mission 3 may fly a lighter one.
    print("gross = Mission 2 takeoff mass, the heaviest case; the rule limit "
          "is 55 lb")
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
                        ref_m3: Optional[float] = None,
                        weather: bool = False,
                        ref_gm: Optional[float] = None) -> ScoreCard:
    """Full dossier for one design.

    `ref_m2` / `ref_m3` / `ref_gm` must be the normalisers the search
    actually used.  Without them a self-normalised run would rescore this
    design against the field-mode defaults and print a total above its
    own 7.5 ceiling.
    """
    ac = Aircraft(design, cfg)
    card = evaluate_scores(ac, cfg, dx, ref_m2, ref_m3, ref_gm)
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
    print("  mean chord (S/b)       %.3f m   (%.2f in)   MAC %.3f m"
          % (design.chord_m, design.chord_m * M2IN, ac.wing.mac_m))
    print("  area                   %.4f m2  (%.0f in2)"
          % (s["area_m2"], s["area_m2"] / (IN2M ** 2)))
    print("  aspect ratio           %.2f" % s["AR"])
    print("  taper ratio            %.2f  (root %.0f mm, tip %.0f mm)"
          % (design.taper, 1000 * ac.wing.root_chord_m,
             1000 * ac.wing.tip_chord_m))
    print("  span efficiency e      %.3f  (lifting line %.3f x fuselage %.3f)"
          % (s["oswald_e"], ac.aero_clean.e_wing,
             ac.aero_clean.fuselage_e_factor))
    if design.winglet_h_m > 0.0:
        print("  winglets               %.0f mm high each (%.1f%% of span), "
              "root chord %.0f mm; effective AR %.2f (x%.3f)"
              % (1000 * design.winglet_h_m,
                 100 * design.winglet_h_m / design.span_m,
                 1000 * ac.wing.tip_chord_m, ac.aero_clean.AR_eff,
                 ac.aero_clean.winglet_ar_factor))
    print("  CL max (3D, built)     %.2f" % s["CLmax_3d"])
    _a = ac.aero_clean
    _thr = lambda v: ac.prop_sys.thrust(v, 1.0, 1.0)
    print("  M2 stall               %.1f m/s power off, %.1f m/s full power"
          % (_a.v_stall(ac.w_m2_n), _a.v_stall_powered(ac.w_m2_n, _thr)))
    print("  stall starts at        %.0f%% of the semi-span  (planform "
          "factor %.3f)" % (100 * ac.aero_clean.stall_eta,
                            ac.aero_clean.stall_factor))
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
    print("  tail arm               %.3f m  (lightest tails + tailcone "
          "%.3f m; the tailcone taper needs %.3f m)"
          % (b.tail_arm_m, b.tail_arm_optimum_m, b.tail_arm_min_m))
    print("  horizontal tail area   %.4f m2 (%.0f in2)"
          % (b.s_horiz_m2, b.s_horiz_m2 / (IN2M ** 2)))
    print("  vertical tail area     %.4f m2 (%.0f in2)"
          % (b.s_vert_m2, b.s_vert_m2 / (IN2M ** 2)))
    _vh = b.s_horiz_m2 * b.tail_arm_m / (ac.wing.area_m2 * ac.wing.mac_m)
    _vv = b.s_vert_m2 * b.tail_arm_m / (ac.wing.area_m2 * ac.wing.span_m)
    _ts = ac.tail_sizing
    if _ts.get("mode") == "requirements":
        print("  tail volumes           Vh %.3f, Vv %.4f - sized from the "
              "requirements%s" % (_vh, _vv,
                                  ", Vh chosen %.2f" % design.tail_vol_h
                                  if design.tail_vol_h > 0 else ""))
        print("    horizontal by %-9s trim %.0f in2, rotation %.0f in2 "
              "(%.1f m/s), floor %.0f in2"
              % (_ts["h_sized_by"].split(" ")[0],
                 _ts["s_h_trim_m2"] / IN2M ** 2,
                 _ts["s_h_rotation_m2"] / IN2M ** 2, _ts["v_rotate_mps"],
                 _ts["s_h_floor_m2"] / IN2M ** 2))
        print("    vertical by %-11s stability %.0f in2, motor out %.0f "
              "in2, floor %.0f in2"
              % (_ts["v_sized_by"], _ts["s_v_stability_m2"] / IN2M ** 2,
                 _ts["s_v_motor_out_m2"] / IN2M ** 2,
                 _ts["s_v_floor_m2"] / IN2M ** 2))
    else:
        print("  tail volumes           Vh %.3f, Vv %.4f - fixed" % (_vh, _vv))
    _ct, _cat = ac.tail_trim_cl()
    print("  tail CL to trim        %.2f at the %s (max %.2f)"
          % (_ct, _cat, ac.cfg.build.tail_cl_max))

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
    print("  CG                     %.3f c   (from where the parts sit - "
          "LAYOUT below)" % a.x_cg_over_c)
    print("  STATIC MARGIN          %.3f c   %s"
          % (a.static_margin,
             "stable" if a.static_margin > 0 else "UNSTABLE"))
    print("    a positive margin means it will fly hands-off; verify in AVL")
    print_layout(ac)

    # ---------------- structure ----------------
    # The team asked for a verdict here, not a stress report: the load
    # factors and spar wall are still computed (the wing mass depends on
    # them, so they cannot be skipped), but what gets printed is whether
    # the wing holds together and what drove the case.
    lc = ac.load_case
    sp = ac.spar
    print("\nSTRUCTURE  (will it break?)")
    if ac.structure_ok:
        print("  VERDICT                NO - the spar carries every wing load")
        print("                         case at the %.2fx factor of safety,"
              % sp.fos_required)
        print("                         and the tube it needs fits inside the")
        print("                         aerofoil all the way to the tip.")
    else:
        print("  VERDICT                YES, IT BREAKS - at %.0f%% span the"
              % (100 * sp.worst_eta))
        print("                         spar needs a wall thicker than the")
        print("                         section is deep.  Shorten the span,")
        print("                         deepen the aerofoil, or carry less.")
    print("  flight load            %.2f g limit (%s), x%.2f = %.2f g ultimate"
          % (ac.n_limit, lc.driver, sp.fos_required, ac.n_ultimate))
    print("  factor of safety       required %.2f; achieved, at the weakest"
          % sp.fos_required)
    print("                         station along the span:")
    for key, label in (("flight bending", "flight bending"),
                       ("wing-tip test", "wing-tip lift test (4.2.4)"),
                       ("shear + torsion", "shear + torsion")):
        got = sp.fos.get(key, float("inf"))
        print("    %-26s %s  at %3.0f%% span"
              % (label, ("%5.2f" % got) if got < 99.0 else " >99 ",
                 100 * sp.fos_at_eta.get(key, 0.0)))
    if ac.spar_gauge_limited:
        print("  the spar is GAUGE limited: the minimum practical wall, not")
        print("  load, sets it everywhere, so every factor is above %.2f."
              % sp.fos_required)

    # ---------------- propulsion ----------------
    print("\nPROPULSION  (T)")
    print("  motor                  %s  x%d" % (ac.motor.name, design.n_motors))
    print("    Kv %.0f rpm/V, Rm %.3f ohm, I0 %.2f A at %.1f V (%.2f A at "
          "the static rpm), %.0f g each"
          % (ac.motor.kv, ac.motor.rm, ac.motor.i0, ac.motor.i0_volts,
             ac.motor.i0_at(rc["rpm"] / 60.0, cfg.build.i0_speed_share),
             ac.motor.mass_kg * 1000))
    print("  propeller              %s  (%.1f x %.1f in)"
          % (ac.prop.key, ac.prop.diameter_m * M2IN, ac.prop.pitch_m * M2IN))
    _e = ac.prop_sys.esc
    print("  ESC                    %s x%d"
          % (_e.name or _e.key, design.n_motors))
    print("    %.0f A published; holds %.0f A for %.0f s at %.0f m/s cooling "
          "air, %.0f C" % (_e.i_cont_a, rc["esc_sustained_a"],
                           cfg.rules.mission_window_s,
                           cfg.build.esc_airflow_mps, cfg.env.temp_c))
    print("  propulsion systems     %d  (one pack, fuse and arming plug"
          " each)" % ac.n_systems)
    if ac.n_wing_motors:
        print("  wing power leads       %s" % wing_leads_text(ac, cfg.build))
        print("  wing nacelles          %d x %.0f g  (%.0f%% of motor + ESC)"
              % (ac.n_wing_motors,
                 1000 * ac.mass_breakdown.items.get(
                     "wing nacelles (%d)" % ac.n_wing_motors, 0.0)
                 / ac.n_wing_motors, 100 * cfg.build.nacelle_mass_frac))
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
    print("    peak at %s, full charge and throttle (static %.0f A)"
          % ("%.0f m/s" % rc["peak_current_speed_mps"]
             if rc["peak_current_speed_mps"] > 0 else "standstill",
             rc["static_current_a"]))
    _i_fuse = ac.fuse_peak_a()
    print("    fuse %.0f A (C x Ah, at most 100 A) on each of %d system(s);"
          " peak %.0f A per pack%s -> %s"
          % (rc["fuse_rating_a"], ac.n_systems, _i_fuse,
             " with %.0f C windings" % cfg.env.temp_c
             if cfg.env.temp_c < 20.0 else "",
             "inside it" if _i_fuse <= rc["fuse_rating_a"] + 1e-9
             else "OVER"))
    print("  margins  motor I %.2fx  motor P %.2fx  ESC %.2fx  pack C %.2fx"
          "  fuse %.2fx"
          % (rc["motor_current_margin"], rc["motor_power_margin"],
             rc["esc_current_margin"], rc["pack_current_margin"],
             rc["fuse_rating_a"] / max(_i_fuse, 1e-6)))

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
    _print_mission("MISSION 1 - flight (no payload, feasibility)", card.r1, cfg)
    _print_mission("MISSION 2 - delivery flight", card.r2, cfg)
    _print_mission("MISSION 3 - sensor flight", card.r3, cfg)
    if weather:
        print_weather_plan(design, cfg, dx)

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

    print_margins(ac, cfg)

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


def print_weather_plan(design: Design, cfg: Config,
                       dx: float = 6.0) -> Dict[str, Any]:
    """Which propeller to fly on each weather preset, and whether Mission 2
    completes on it.  RULE 3.2.1e allows a different propeller for each
    flight attempt; see optimizer.weather_plan."""
    from .optimizer import weather_plan
    plan = weather_plan(design, cfg, dx)
    print("\nWEATHER PLAN - the propeller to fly on each day (rule 3.2.1e: it")
    print("may change for each flight attempt; everything else stays as built)")
    for preset, e in plan["presets"].items():
        if e["m2_ok"]:
            out = "Mission 2 in %.1f s" % e["m2_time_s"]
        elif e["legal"]:
            out = "Mission 2 FAILS: " + e["reason"][:40]
        else:
            out = "NO LEGAL PROPELLER: " + e["reason"][:36]
        print("  %-15s %-11s %s%s" % (preset, e["prop"].replace("APC_", "APC "),
              out, "   (swapped from %s)" % design.prop.replace("APC_", "")
              if e["swapped"] else ""))
        if "peak_pack_a" in e:
            print("  %-15s %-11s   peak %.0f A per pack, fuse %.0f A%s"
                  % ("", "", e["peak_pack_a"], e["fuse_a"],
                     "" if e["peak_pack_a"] <= e["fuse_a"] + 1e-9
                     else "  OVER"))
        if e["swapped"] and e["reason"]:
            print("  %-15s %-11s   why: %s" % ("", "", e["reason"][:58]))
    print("  VERDICT                %s" % (
        "flyable on every preset" + (" with a propeller swap" if plan["swaps"] else "")
        if plan["robust"] else
        "NOT flyable on every preset - see above"))
    return plan


# ==========================================================================
# Margins and modelling factors
# ==========================================================================
def margin_rows(ac: Aircraft, cfg: Config) -> List[Dict[str, str]]:
    """Every build-standard number, with its role (config.FIELD_ROLES),
    what it is and - for the margins - what it holds back on this design.

    One list for the report and the GUI, so the two cannot disagree."""
    from dataclasses import fields as dc_fields
    b = cfg.build
    wh = ac.total_watt_hours
    contingency = dict(ac.mass_breakdown.sorted_items()).get("contingency",
                                                             0.0)
    effect = {
        "load_safety_factor": "ultimate %.2f g on a %.2f g limit"
                              % (ac.n_ultimate, ac.n_limit),
        "m2_energy_margin": ("%.1f Wh held back in Mission 2"
                             % (b.m2_energy_margin * wh)
                             if b.m2_energy_margin
                             else "none - the landing point replaced it"),
        "land_at_discharge": "every mission turns for landing with %.1f Wh "
                             "or more left" % ((1.0 - b.land_at_discharge)
                                               * ac.total_energy_wh),
        "depth_of_discharge": "every mission lands with %.1f Wh"
                              % ((1.0 - ac.depth_of_discharge) * wh),
        "energy_reserve": ("%.1f Wh" % (b.energy_reserve * wh)
                           if b.energy_reserve else "none"),
        "full_throttle_flights": ("no cruise speed to hold: every second "
                                  "of each flight charged at full power"
                                  if b.full_throttle_flights else
                                  "planned cruise speeds - the pilot must "
                                  "hold them"),
        "contingency_frac": "%.3f kg (%.2f lb) of the empty mass"
                            % (contingency, contingency * KG2LB),
        "target_static_margin": "built to %.3f c"
                                % ac.aero_clean.static_margin,
        "min_static_margin": "this design %.3f c"
                             % ac.aero_clean.static_margin,
        "prop_clearance_margin_m": "sets %.0f mm main gear legs"
                                   % (1000 * ac.gear_height_m),
        "max_stall_mps": ("this design %.1f mph on the factored lift in "
                          "%s air" % (ac.stall_speed_factored() * MS2MPH,
                                      ac.stall_air()[2])
                          + (" - no limit set" if b.max_stall_mps == float("inf")
                             else ", limit %.1f mph" % (b.max_stall_mps * MS2MPH))),
        "stall_lift_safety_factor": "true stall %.1f mph power off at the "
                                    "Mission 2 weight"
                                    % (ac.aero_clean.v_stall(ac.w_m2_n) * MS2MPH),
        "max_ground_alpha_deg": ("CL %.2f at it against CLmax %.2f - %s"
                                 % (ac.aero_clean.cl_at_ground_alpha(),
                                    ac.aero_clean.cl_max_at_v(
                                        ac.aero_clean.v_stall(ac.w_m2_n)),
                                    "it sets liftoff and touchdown"
                                    if ac.aero_clean.cl_at_ground_alpha()
                                    < ac.aero_clean.cl_max_at_v(
                                        ac.aero_clean.v_stall(ac.w_m2_n))
                                    else "the stall sets them")),
    }
    rows = []
    for f in dc_fields(b):
        role, note = FIELD_ROLES.get(f.name, ("design", ""))
        val = getattr(b, f.name)
        if val is None:
            continue
        rows.append({"role": role, "name": f.name,
                     "value": ("%.6g" % val if isinstance(val, float)
                               else str(val)),
                     "note": note, "effect": effect.get(f.name, "")})
    return rows


def awg_of(area_mm2: float) -> int:
    """The American wire gauge nearest a copper cross-section, mm^2."""
    d = math.sqrt(4.0 * max(area_mm2, 1e-6) / math.pi)
    return int(round(36.0 - 39.0 * math.log(d / 0.127) / math.log(92.0)))


def wing_leads_text(ac: Aircraft, b) -> str:
    """One line on a design's wing-motor power leads, or '' without."""
    if not ac.n_wing_motors:
        return ""
    area = ac.motor.i_max_a / max(b.lead_current_density_a_mm2, 0.1)
    r_each = ac.lead_r_ohm * ac.design.n_motors / ac.n_wing_motors
    return ("%d x %.2f m pairs of %.1f mm2 copper (%d AWG): %.1f mOhm "
            "each, %.0f g" % (ac.n_wing_motors, max(ac.lead_lengths_m),
                              area, awg_of(area), 1000 * r_each,
                              1000 * ac.lead_mass_kg))


def computed_rows(ac: Aircraft, cfg: Config) -> List[Dict[str, str]]:
    """The corrections the model now COMPUTES for this design, each with
    the fixed factor it replaced."""
    a, b = ac.aero_clean, cfg.build
    d = ac.prop.diameter_m
    a_disc = max(ac.design.n_motors, 1) * math.pi * 0.25 * d * d
    v_rep = ac._max_level_speed(ac.w_m2_n)       # the design's own top speed
    parts = a.drag_parts(v_rep)
    wm = ac.prop.windmill_drag_coeff
    thr = lambda v: ac.prop_sys.thrust(v, 1.0, 1.0)
    vs_off = a.v_stall(ac.w_m2_n)
    vs_on = a.v_stall_powered(ac.w_m2_n, thr)
    v_t = 1.3 * vs_off
    gl, gm = a.slipstream_gains(v_t, thr(v_t))
    sharp = 2.0 * (ac.body.fuse_width_m + ac.body.fuse_height_m)
    rc_rpm = ac.prop_sys.rating_check()["rpm"]
    return [
        {"name": "span efficiency e",
         "value": "%.3f" % a.e,
         "note": "lifting line %.3f x fuselage %.3f  (Kroo: 1 - 2 (d/b)^2, "
                 "d/b %.3f)" % (a.e_wing, a.fuselage_e_factor,
                                a.fuse_d_over_b),
         "was": "ignored"},
        {"name": "CLmax over the fuselage",
         "value": "x%.3f" % a.fuselage_clmax_factor,
         "note": "the fuselage covers %.1f%% of the planform; carry-over "
                 "%.2f" % (100 * a.covered_frac, b.fuselage_lift_carryover),
         "was": "ignored"},
        {"name": "wing / fuselage overlap",
         "value": "%.4f m2" % (2.0 * (1.0 - a.wing_wetted_frac) * a.S
                               + a.fuse_covered_m2),
         "note": "skin that is not wetted: the wing's covered surface and "
                 "the fuselage under it",
         "was": "counted as wetted"},
        {"name": "junction interference",
         "value": "%.5f m2" % a.junction_drag_area_m2,
         "note": "CD*A, Hoerner: 2 wing roots, 3 tail roots%s "
                 "(%.1f%% of CD0)"
                 % (", %d nacelle sides" % (2 * ac.n_wing_motors)
                    if ac.n_wing_motors else "",
                    100 * parts["junctions"] / max(a.cd0(v_rep), 1e-9)),
         "was": "inside a flat +10% (nacelle sides: ignored)"},
        {"name": "installed thrust",
         "value": "%.3f" % ac.install_thrust_factor_at(v_rep),
         "note": "at the %.1f m/s top speed: slipstream scrubs %.4f m2 of "
                 "drag area; disc %.3f m2"
                 % (v_rep, a.slipstream_drag_area_m2(d, v_rep), a_disc),
         "was": "flat 0.95"},
        {"name": "no-load current",
         "value": "%.2f A" % ac.motor.i0_at(rc_rpm / 60.0, b.i0_speed_share),
         "note": "at the static %.0f rpm, against the maker's %.2f A at "
                 "%.1f V (%.0f rpm); %.0f%% of it grows with speed (eddy "
                 "current)" % (rc_rpm, ac.motor.i0, ac.motor.i0_volts,
                               ac.motor.kv * ac.motor.i0_volts,
                               100 * b.i0_speed_share),
         "was": "the 10 V figure at every speed"},
        {"name": "wing-motor power leads",
         "value": ("%.1f mOhm" % (1000 * ac.lead_r_ohm)
                   if ac.n_wing_motors else "none"),
         "note": ("mean per motor, in the power solve: %s"
                  % wing_leads_text(ac, b)) if ac.n_wing_motors
                 else "no wing motors: every ESC sits by the battery",
         "was": "inside the fixed 3 mOhm harness"},
        {"name": "wing nacelles",
         "value": ("%.0f g" % (1000 * ac.mass_breakdown.items.get(
             "wing nacelles (%d)" % ac.n_wing_motors, 0.0)))
                  if ac.n_wing_motors else "none",
         "note": ("%d pods at %.0f%% of the motor + ESC they carry"
                  % (ac.n_wing_motors, 100 * b.nacelle_mass_frac))
                 if ac.n_wing_motors else "no wing motors",
         "was": "ignored (only the 25 g mount + spinner)"},
        {"name": "windmilling propeller",
         "value": ("%.3f" % wm) if wm is not None
                  else "%.3f" % b.windmill_drag_coeff,
         "note": ("Cd on disc area, from the %s table at zero torque"
                  % ac.prop.key) if wm is not None
                 else "no measured table: the build-standard fallback",
         "was": "flat 0.08"},
        {"name": "tail form factor",
         "value": "%.3f" % tail_form_factor(b.tail_t_over_c),
         "note": "Raymer, t/c %.2f" % b.tail_t_over_c,
         "was": "flat 1.25"},
        {"name": "fuselage perimeter",
         "value": "%.3f" % (ac.body.section_perimeter_m / max(sharp, 1e-9)),
         "note": "of the sharp-cornered box: rounded corners, radius %.2f "
                 "of the smaller side" % b.fuse_corner_radius_frac,
         "was": "flat 0.92"},
        {"name": "finite-span CLmax",
         "value": "x%.3f" % a.stall_factor,
         "note": "the lifting-line stall alone (first section to stall, at "
                 "%.0f%% semi-span)" % (100 * a.stall_eta),
         "was": "also x min(1, 0.80 + 0.04 AR)"},
        {"name": "slipstream turbulence",
         "value": "%.5f m2" % (parts["slipstream transition"] * a.S),
         "note": "CD*A: the %.3f m2 of wing behind the propeller moves %.0f%% "
                 "of the way to fully turbulent (%.1f%% of CD0)"
                 % (a.washed_strip_m2, 100 * a.transition_frac,
                    100 * parts["slipstream transition"]
                    / max(a.cd0(v_rep), 1e-9)),
         "was": "ignored"},
        {"name": "slipstream lift",
         "value": "x%.3f" % gl,
         "note": "lift at the same angle at full throttle, %.0f m/s; CLmax "
                 "x%.3f; M2 stall %.1f m/s at full power against %.1f "
                 "power off (lifting line with the jet's speed)"
                 % (v_t, gm, vs_on, vs_off),
         "was": "ignored"},
    ]


def print_margins(ac: Aircraft, cfg: Config) -> None:
    """What this answer rests on: the deliberate margins with what they
    hold back, the corrections now computed, and the estimates still
    waiting for a measurement."""
    import textwrap
    pad = " " * 40

    def row(name: str, value: str, text: str) -> None:
        body = textwrap.wrap(text, 78 - 40) or [""]
        print("  %-26s %-10s %s" % (name, value, body[0]))
        for line in body[1:]:
            print(pad + line)

    print()
    print(LINE)
    print("MARGINS AND MODELLING FACTORS  (what this answer rests on)")
    print(LINE)
    rows = margin_rows(ac, cfg)
    for role in ("margin", "calibrate"):
        print()
        print(ROLE_TITLES[role].upper())
        for r in rows:
            if r["role"] == role:
                row(r["name"], r["value"], r["note"]
                    + ((" -> " + r["effect"]) if r["effect"] else ""))
    print()
    print("COMPUTED FOR THIS DESIGN  (each replaced a fixed factor)")
    for r in computed_rows(ac, cfg):
        row(r["name"], r["value"], "%s. Was: %s." % (r["note"], r["was"]))
    print()
    print("  Design choices, operating settings and constants are in the")
    print("  GUI's Advanced tab, grouped the same way.")


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
              ("m3_sensor_frac",
               (design.sensor_m3_kg / design.sensor_mass_kg
                if design.sensor_m3_kg and design.sensor_mass_kg > 0
                else 1.0), b[3]),
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


def _print_mission(title: str, r, cfg: Config) -> None:
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
    print("  stall speed, power off %.1f m/s (%.0f mph)"
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
    if r.motor_temp_c > 0.0:
        print("  motor windings         %.0f C peak  (limit %.0f C, air %.0f C)"
              % (r.motor_temp_c, cfg.build.motor_temp_limit_c,
                 cfg.env.temp_c))
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
        "Kv %.0f, Rm %.3f ohm, I0 %.2f A at %.1f V, %.0f g"
        % (ac.motor.kv, ac.motor.rm, ac.motor.i0, ac.motor.i0_volts,
           ac.motor.mass_kg * 1000),
        design.n_motors, ac.motor.source)
    add("Propulsion", "Propeller %s" % ac.prop.key,
        "%.1f x %.1f in, %.0f g" % (ac.prop.diameter_m * M2IN,
                                    ac.prop.pitch_m * M2IN,
                                    ac.prop.mass_kg * 1000),
        design.n_motors + 2, "buy spares - you will break them")
    _e = ac.prop_sys.esc
    add("Propulsion", "ESC %s" % (_e.name or _e.key),
        "%.0f A published (%.0f A for the %.0f s mission at %.0f m/s "
        "cooling air), %d-%dS, %.0f g"
        % (_e.i_cont_a, _e.sustained_a(cfg.rules.mission_window_s,
                                       cfg.build.esc_airflow_mps,
                                       cfg.env.temp_c),
           cfg.rules.mission_window_s, cfg.build.esc_airflow_mps,
           _e.min_cells, _e.max_cells, _e.mass_kg * 1000),
        design.n_motors, "BEC must be disabled (rule 3.2.3b); " + _e.source)
    # 0 = the maker does not publish the shaft size; say so, don't guess
    if ac.n_wing_motors:
        _a = ac.motor.i_max_a / max(b.lead_current_density_a_mm2, 0.1)
        add("Propulsion", "Wing power leads",
            "%.1f mm2 copper (%d AWG), 2 x %.2f m per wing motor"
            % (_a, awg_of(_a), max(ac.lead_lengths_m)), ac.n_wing_motors,
            "battery bay to the ESC in each nacelle; %.1f mOhm and %.0f g "
            "a pair" % (1000 * ac.lead_r_ohm * design.n_motors
                        / ac.n_wing_motors,
                        1000 * ac.lead_mass_kg / ac.n_wing_motors))
        add("Propulsion", "Wing nacelle", "motor + ESC pod, firewall and "
            "pylon", ac.n_wing_motors, "budgeted at %.0f%% of the motor + "
            "ESC mass - weigh a built one" % (100 * b.nacelle_mass_frac))
    add("Propulsion", "Motor mount + spinner",
        ("for %.0f mm shaft" % ac.motor.shaft_mm if ac.motor.shaft_mm > 0
         else "shaft size not published - check with the maker"),
        design.n_motors, "")

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
        "arming plug. A pack with two output leads: combine them in this "
        "harness BEFORE the fuse, and it feeds this system only "
        "(rule 3.2.3.1e)")
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
    sp = ac.spar
    add("Structure", "Carbon spar tube",
        "OD %.1f mm root to %.1f mm tip, wall %.2f mm max, length %.0f mm"
        % (2000 * sp.r_root_m, 2000 * sp.r_tip_m, 1000 * sp.t_used_m,
           1000 * ac.wing.span_m),
        1,
        "factor of safety %.2f on flight (n_ult %.2f, %s), tip-lift test, "
        "shear and torsion"
        % (sp.fos_required, ac.n_ultimate, ac.load_case.driver))
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
    print("SENSITIVITY  (objective = M2 + M3 + GM, one variable at a time)")
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


def print_weights_and_speeds(design: Design, cfg: Config,
                             dx: float = 6.0) -> None:
    """Just the weights and the estimated speeds of one design - what
    `run_optimizer.py --no-report` prints instead of the full dossier."""
    from .mission import (lap_ground_distance, mission_1, mission_2,
                          mission_3, turn_capability)
    ac = Aircraft(design, cfg)
    print("\nWEIGHTS AND SPEEDS  (%s, span %.2f m, AR %.1f, taper %.2f, "
          "%s, %s)" % (design.airfoil, ac.wing.span_m, ac.aspect_ratio,
                       design.taper, ac.motor.name, design.prop))
    print("  sensor %.1f lb (%.2f kg) | empty %.1f lb | M2 gross %.1f lb "
          "(%.2f kg) | M3 gross %.1f lb"
          % (ac.sensor_max_kg * KG2LB, ac.sensor_max_kg,
             ac.empty_mass_kg * KG2LB, ac.mtow_m2_kg * KG2LB,
             ac.mtow_m2_kg, ac.mtow_m3_kg * KG2LB))

    def spd(v: float) -> str:
        return "%.1f m/s (%.0f mph)" % (v, v * MS2MPH)

    for name, fn, aero, w in (("M1", mission_1, ac.aero_clean, ac.w_m1_n),
                              ("M2", mission_2, ac.aero_clean, ac.w_m2_n),
                              ("M3", mission_3, ac.aero_deployed, ac.w_m3_n)):
        r = fn(ac, cfg, dx=dx)
        if not r.ok:
            print("  %s  did not complete: %s" % (name, r.reason))
            continue
        turn = turn_capability(aero, ac.prop_sys, w, 0.70, cfg, ac.n_limit)
        laps = r.lap_times[1:] or r.lap_times     # lap 1 includes the climb
        avg = lap_ground_distance(cfg, turn) / (sum(laps) / len(laps))
        print("  %s  stall %s | top %s | turn %s | average lap %s | "
              "%d laps in %.0f s | windings %.0f C"
              % (name, spd(r.v_stall), spd(r.v_max), spd(r.v_turn),
                 spd(avg), r.laps, r.time_s, r.motor_temp_c))


def print_layout(ac) -> None:
    """Where every part sits, from the nose, and the balance it gives."""
    lay = ac.layout
    b = ac.cfg.build
    print("\nLAYOUT  (every station calculated - dbf2027/layout.py)")
    print("  nose                   %.3f m  (motor section %.3f, battery "
          "%.3f, avionics %.3f)" % (lay.nose_length_m, lay.motor_section_m,
                                     lay.battery_length_m,
                                     lay.avionics_length_m))
    print("  payload bay            %.3f - %.3f m  (centred on the CG)"
          % (lay.bay_start_m, lay.bay_end_m))
    print("  tailcone               %.3f m  to the tail post"
          % lay.tailcone_length_m)
    print("  fuselage length        %.3f m  (%.1f in)"
          % (lay.fuse_length_m, lay.fuse_length_m / IN2M))
    print("  wing leading edge      %.3f m at the root, %.3f m at the mean "
          "chord" % (lay.x_wing_le_root_m, lay.x_wing_le_mac_m))
    print("  tail aerodynamic ctr   %.3f m  (arm %.3f m)"
          % (lay.x_tail_ac_m, lay.tail_arm_m))
    print("  battery centre         %.3f m  %s"
          % (lay.x_battery_m, "(behind the payload bay)"
             if lay.x_battery_m > lay.bay_start_m else "(in the nose)"))
    print("  gear                   nose %.3f m, mains %.3f m (%.0f deg "
          "tip-back)" % (lay.x_nose_gear_m, lay.x_main_gear_m,
                         b.gear_tipback_deg))
    print("  CG                     M2 %.3f m, M1 (empty) %.3f m, M3 %.3f m;"
          " target %.3f m" % (lay.x_cg_m, lay.x_cg_empty_m, lay.x_cg_m3_m,
                              lay.x_cg_target_m))
    print("  static margin          M2 %.3f c, M1 %.3f c  (target %.2f, "
          "band to %.2f)" % (lay.static_margin, lay.static_margin_empty,
                             b.target_static_margin,
                             b.target_static_margin + b.static_margin_band))
    print("  balance                %s" % lay.note)
    print("  %-40s %8s %8s %9s" % ("station", "kg", "x m", "moment"))
    for st in sorted(lay.stations, key=lambda q: q.x_m):
        print("  %-40s %8.3f %8.3f %9.3f"
              % (st.name[:40], st.mass_kg, st.x_m,
                 st.mass_kg * (st.x_m - lay.x_cg_target_m)))


def print_top_configurations(rows, cfg: Config, n: int = 10,
                             out_dir: str = "", dx: float = 1.5) -> None:
    """A breakout of each of the top `n` ranked designs: configuration,
    weights, speeds, mission results, motor temperature, the propeller for
    each weather preset and the layout.  With `out_dir`, each design is
    also written to out_dir/top<n>/NN_<airfoil>_<motor>.json so any of
    them can be loaded into the GUI or --evaluate.

    The missions are flown at `dx`, the dossier's step: a design on the
    energy limit can lose its last Mission 3 lap at a coarser one, and the
    breakout then disagrees with the dossier for the same aeroplane."""
    top = list(rows)[:n]
    sub = "top%d" % n
    print("\n" + "=" * 78)
    print("TOP %d CONFIGURATIONS" % len(top))
    print("=" * 78)
    if out_dir:
        os.makedirs(os.path.join(out_dir, sub), exist_ok=True)
    for i, r in enumerate(top, 1):
        d = r["design"]
        ac = Aircraft(d, cfg)
        wx = r.get("weather") or {}
        code = wx.get("code", r.get("wx", "?"))
        print("\n#%d  score %.3f  (M1 %.0f + M2 %.3f + M3 %.3f + GM %.3f)"
              "   weather: %s"
              % (i, r.get("total_mission", 0.0), r.get("m1", 0.0),
                 r.get("m2", 0.0), r.get("m3", 0.0), r.get("gm", 0.0),
                 {"ok": "flies every preset on its own propeller",
                  "swap": "flies every preset with a propeller swap",
                  "no": "FAILS a weather preset"}.get(code, code)))
        print("  %s, span %.3f m, chord %.3f m, AR %.2f, taper %.2f"
              % (d.airfoil, ac.wing.span_m, ac.wing.chord_m,
                 ac.aspect_ratio, d.taper))
        print("  %d x %s on %s, %d x %s (%.1f Wh, %dS)"
              % (d.n_motors, ac.motor.name, d.prop, ac.n_systems,
                 ac.pack.name, ac.total_watt_hours, ac.pack.series))
        if wx.get("presets"):
            print("  propellers: %s" % ", ".join(
                "%s %s%s" % (p.replace("tucson-", "").replace(
                    "tucson", "nominal"), e["prop"].replace("APC_", ""),
                    "" if e["m2_ok"] else " (M2 FAILS)")
                for p, e in wx["presets"].items()))
        lay = ac.layout
        print("  fuselage %.2f m, static margin %.3f c, battery %s"
              % (lay.fuse_length_m, lay.static_margin,
                 "behind the bay" if lay.x_battery_m > lay.bay_start_m
                 else "in the nose"))
        print_weights_and_speeds(d, cfg, dx)
        if out_dir:
            name = "%02d_%s_%s.json" % (i, d.airfoil, d.motor)
            with open(os.path.join(out_dir, sub, name), "w",
                      encoding="utf-8") as fh:
                json.dump(asdict(d), fh, indent=2)
    if out_dir:
        print("\nwrote the %d designs to %s" % (len(top),
                                                os.path.join(out_dir, sub)))
