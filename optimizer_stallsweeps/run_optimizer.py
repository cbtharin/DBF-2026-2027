#!/usr/bin/env python3
"""DBF 2026-27 design optimiser - command line entry point.

    python run_optimizer.py --quick
    python run_optimizer.py --wind 6 --altitude 400 --max-wh 100
    python run_optimizer.py --thorough --field-m2 0.07 --field-m3 26
    python run_optimizer.py --evaluate examples/baseline.json

Pure standard library: no pip install required.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import asdict, replace

from dbf2027.config import CONDITIONS, DEFAULT_CONDITIONS, SITE_NAME
from dbf2027 import (Aircraft, Config, Design, DesignSpace, Environment,
                     SearchSettings, evaluate_design, load_motors_csv,
                     load_packs_csv, load_props_csv, optimise, sensitivity)
from dbf2027.report import (export_bom_csv, export_csv, export_json,
                            print_active_bounds, print_bom,
                            print_design_report, print_pareto,
                            print_ranking, print_sensitivity,
                            print_top_configurations,
                            print_weights_and_speeds)
from dbf2027.scoring import rescore
from dbf2027.units import FT2M


def build_config(args) -> Config:
    cfg = Config()
    # A named condition sets all three together; an explicit --altitude,
    # --temp or --wind then overrides just that one.
    alt, temp, wind = CONDITIONS.get(args.conditions,
                                     CONDITIONS[DEFAULT_CONDITIONS])
    if args.altitude is not None:
        alt = args.altitude
    if args.temp is not None:
        temp = args.temp
    if args.wind is not None:
        wind = args.wind
    cfg.env = Environment.at_density_altitude(alt, temp, wind)
    if args.air_density > 0:
        # Viscosity still comes from the temperature; only the density is
        # replaced, because that is the part a weather station gives you.
        cfg.env = Environment(rho=args.air_density, nu=cfg.env.nu,
                              wind_mps=cfg.env.wind_mps,
                              field_elev_m=cfg.env.field_elev_m,
                              temp_c=cfg.env.temp_c)
    if args.takeoff_ft > 0:
        cfg.rules.takeoff_distance_m = args.takeoff_ft * FT2M
    cfg.rules.max_watt_hours = args.max_wh
    cfg.rules.normalisation = args.normalise
    cfg.rules.max_simulators = args.max_sim
    cfg.rules.mission_window_s = args.window
    cfg.rules.field_best_m2_lb_per_s = args.field_m2
    cfg.rules.field_best_m3_lap_lb = args.field_m3
    if args.max_span_in > 0:
        cfg.rules.max_wingspan_m = args.max_span_in * 0.0254
    if args.sensor_max > 0:
        cfg.rules.sensor_mass_max_kg = args.sensor_max * 0.45359237
    if args.sim_max > 0:
        cfg.rules.sim_mass_max_kg = args.sim_max * 0.45359237
    if args.lap_ft > 0:
        # scale the default pattern to a measured lap length, keeping the
        # proportions of the classic course
        scale = (args.lap_ft * FT2M) / cfg.course.straight_length_m
        cfg.course.segments = [(k, v * scale) if k == "straight" else (k, v)
                               for k, v in cfg.course.segments]
    return cfg


def build_parser() -> argparse.ArgumentParser:
    """All command-line options, in one place.

    Separated from main() so the validation suite can check that
    every flag reaches the config value it claims to, without
    running a search or reaching into module globals.
    """
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    # environment
    # These default to None, not to numbers: --conditions supplies all
    # three, and any one given explicitly overrides just that one.
    p.add_argument("--wind", type=float, default=None,
                   help="steady wind along the runway, m/s "
                        "(default: from --conditions)")
    p.add_argument("--altitude", type=float, default=None,
                   help="field elevation, m (default: from --conditions)")
    p.add_argument("--temp", type=float, default=None,
                   help="air temperature, C (default: from --conditions)")
    # rules
    p.add_argument("--takeoff-ft", type=float, default=0.0,
                   help="takeoff distance limit in ft; 0 = no limit (2026-27)")
    p.add_argument("--max-wh", type=float, default=100.0,
                   help="battery energy limit, Wh (VERIFY)")
    p.add_argument("--max-sim", type=int, default=40,
                   help="cap on simulators the search will try (no rules max; "
                        "drag and weight are what actually limit it)")
    p.add_argument("--window", type=float, default=300.0,
                   help="mission flight window, s (default 300)")
    p.add_argument("--max-span-in", type=float, default=0.0,
                   help="wingspan limit in inches; 0 keeps the 2026-27 rule "
                        "of 6 ft (72 in)")
    p.add_argument("--sensor-max", type=float, default=0.0,
                   help="upper search bound on sensor mass in lb (no rules "
                        "cap; 0 keeps the wide default)")
    p.add_argument("--sim-max", type=float, default=0.0,
                   help="upper search bound on mass per simulator in lb "
                        "(no rules cap; 0 keeps the wide default)")
    p.add_argument("--lap-ft", type=float, default=0.0,
                   help="total straight-leg length per lap in ft; 0 keeps the "
                        "2026-27 course (2 x 1000 ft, 180 at each end, 360 at "
                        "midfield)")
    p.add_argument("--normalise", default="self",
                   choices=("self", "field"),
                   help="score against the best design in this run (self, "
                        "default) or against an estimate of the competition")
    # scoring normalisers
    p.add_argument("--field-m2", type=float, default=0.060,
                   help="estimated best weight/time in the field, lb/s")
    p.add_argument("--field-m3", type=float, default=22.0,
                   help="estimated best laps x sensor weight in the field, lap*lb")
    # search
    p.add_argument("--quick", action="store_true", help="fast, coarse search")
    p.add_argument("--thorough", action="store_true", help="slow, wide search")
    p.add_argument("--medium", action="store_true",
                   help="between quick and the default; for a search held "
                        "to one motor count")
    p.add_argument("--max-stall-mph", type=float, default=0.0,
                   help="fastest stall allowed, mph: power off at the Mission "
                        "2 weight, on the lift divided by --stall-lift-fos "
                        "(BuildStandard.max_stall_mps; 0 = no limit)")
    p.add_argument("--stall-air", default="",
                   help="the air --max-stall-mph is checked in: a --conditions "
                        "name, or 'mission' for the mission's own (default "
                        "tucson-hot, the hottest weather preset)")
    p.add_argument("--stall-lift-fos", type=float, default=0.0,
                   help="factor of safety on lift for --max-stall-mph "
                        "(default 1.5)")
    p.add_argument("--max-ground-alpha", type=float, default=0.0,
                   help="most wing angle of attack on the runway, deg: "
                        "tail-strike pitch plus wing incidence "
                        "(BuildStandard.max_ground_alpha_deg, default 15)")
    p.add_argument("--max-chord-in", type=float, default=0.0,
                   help="longest wing chord allowed, inches (the root chord; "
                        "BuildStandard.max_chord_m; 0 = no limit)")
    p.add_argument("--ar", nargs=2, type=float, default=None,
                   metavar=("MIN", "MAX"),
                   help="aspect-ratio range to search (default 1.5 14; the "
                        "model goes down to 1.0)")
    p.add_argument("--tail-vh", nargs=2, type=float, default=None,
                   metavar=("MIN", "MAX"),
                   help="horizontal tail volume range searched (default 0.30 "
                        "0.90); the tail is never smaller than its "
                        "requirements call for")
    p.add_argument("--winglets", nargs=2, type=float, default=None,
                   metavar=("MIN", "MAX"),
                   help="winglet height range searched, as a share of the span "
                        "(default 0 0.15; 0 0 for none)")
    p.add_argument("--taper", nargs=2, type=float, default=None,
                   metavar=("MIN", "MAX"),
                   help="taper ratio range to search (default 1 1: a "
                        "rectangular wing; 0.2 1.0 searches it)")
    p.add_argument("--nacelle-cd", type=float, default=None,
                   help="wing nacelle drag coefficient on frontal area "
                        "(BuildStandard.nacelle_cd, default 0.20) - for a "
                        "sensitivity run")
    p.add_argument("--motor-counts", nargs="*", type=int, default=[],
                   help="motor counts to search, e.g. 2, or 1 2 3 "
                        "(default 1-4; above 4 also lifts the team's "
                        "max_motors cap to match)")
    p.add_argument("--processes", type=int, default=0,
                   help="worker processes, 0 = cpu_count-1")
    p.add_argument("--seed", type=int, default=20262027)
    p.add_argument("--top", type=int, default=10,
                   help="designs broken out in full and written one per "
                        "file to OUT/topN (default 10)")
    p.add_argument("--seeds", nargs="*", default=[],
                   help="extra known designs every search also refines: "
                        "JSON files, or folders of them (the examples/ "
                        "designs are always seeded unless --no-seeds)")
    p.add_argument("--no-seed-refine", action="store_true",
                   help="use the seeds only as the screen's reference "
                        "airframes, not as designs to refine (for comparing "
                        "restricted searches on equal terms)")
    p.add_argument("--no-seeds", action="store_true",
                   help="refine only what the search samples")
    p.add_argument("--keep-propulsion", type=int, default=0,
                   help="motor/prop/pack sets carried past the screen, 0 = "
                        "the mode's own (260 thorough); raise it with the "
                        "catalogue")
    p.add_argument("--refine", type=int, default=0,
                   help="configurations locally refined, 0 = the mode's "
                        "own (80 thorough)")
    p.add_argument("--swap-parents", type=int, default=0,
                   help="refined leaders given every single-part swap, 0 = "
                        "the mode's own (10 thorough)")
    p.add_argument("--swap-refine", type=int, default=0,
                   help="swaps that beat their parent then refined in full, "
                        "0 = 1.5 x --swap-parents or the mode's own (15 "
                        "thorough)")
    p.add_argument("--fill-cores", action="store_true",
                   help="round the refinements, the swaps' refinements and "
                        "the finalists up to whole rounds of the workers, "
                        "and refine the best configurations from more than "
                        "one start when there are too few (for the cluster)")
    p.add_argument("--finalists", type=int, default=0,
                   help="configurations carried into final scoring, 0 = "
                        "the mode's own (8 quick, 15 default, 25 "
                        "thorough); a few fail the final checks, so ask "
                        "for more than --top")
    # data overrides
    p.add_argument("--motors-csv", default="", help="measured motor data")
    p.add_argument("--props-csv", default="", help="measured propeller data")
    p.add_argument("--conditions", default=DEFAULT_CONDITIONS,
                   choices=sorted(CONDITIONS),
                   help="named weather for the fly-off site; sets "
                        "altitude, temperature and wind together")
    p.add_argument("--air-density", type=float, default=0.0,
                   help="kg/m3; overrides the value derived from "
                        "--altitude and --temp, for when you have "
                        "measured it at the field")
    p.add_argument("--packs-csv", default="",
                   help="your own COTS pack shortlist, straight off the "
                        "labels: key,name,chemistry,series,capacity_mah,"
                        "c_rating,mass_g")
    # modes
    p.add_argument("--evaluate", default="",
                   help="evaluate a single design from a JSON file instead of searching")
    p.add_argument("--airfoils", default="",
                   help="comma-separated airfoil keys to search, e.g. "
                        "E214,NACA4412 (default: the whole catalogue)")
    p.add_argument("--motors", default="",
                   help="comma-separated motor keys to search, e.g. "
                        "NEU_6521_3P5_362 (default: the whole catalogue)")
    p.add_argument("--props", default="",
                   help="comma-separated propeller keys to search, e.g. "
                        "APC_20x22.5CD (default: the whole catalogue)")
    p.add_argument("--packs", default="",
                   help="comma-separated battery-pack keys to search, e.g. "
                        "MX_12S2000_27C (default: the searchable catalogue)")
    p.add_argument("--out", default="results", help="output directory")
    p.add_argument("--no-export", action="store_true")
    p.add_argument("--no-report", action="store_true",
                   help="after the search print only the best design's "
                        "weights and speeds: skips the full report, its "
                        "fine-step weather re-check, the parts list and "
                        "the sensitivity study (minutes, mostly on one "
                        "core)")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.motors_csv:
        print("loaded %d motors from %s"
              % (load_motors_csv(args.motors_csv), args.motors_csv))
    if args.props_csv:
        print("loaded %d propellers from %s"
              % (load_props_csv(args.props_csv), args.props_csv))
    if args.packs_csv:
        print("loaded %d battery packs from %s"
              % (load_packs_csv(args.packs_csv), args.packs_csv))

    cfg = build_config(args)
    if args.nacelle_cd is not None:
        cfg.build.nacelle_cd = args.nacelle_cd
        print("  nacelle Cd         %.2f on frontal area" % args.nacelle_cd)
    if args.max_chord_in > 0:
        cfg.build.max_chord_m = args.max_chord_in * 0.0254
        print("  chord limit        %.1f in (root chord)" % args.max_chord_in)
    if args.max_ground_alpha > 0:
        cfg.build.max_ground_alpha_deg = args.max_ground_alpha
        print("  runway attitude    %.1f deg wing angle of attack at most"
              % args.max_ground_alpha)
    if args.stall_lift_fos > 0:
        cfg.build.stall_lift_safety_factor = args.stall_lift_fos
    if args.stall_air:
        if args.stall_air != "mission" and args.stall_air not in CONDITIONS:
            raise SystemExit("--stall-air: %s is not a --conditions name or "
                             "'mission'" % args.stall_air)
        cfg.build.stall_check_conditions = ("" if args.stall_air == "mission"
                                            else args.stall_air)
    if args.max_stall_mph > 0:
        cfg.build.max_stall_mps = args.max_stall_mph * 0.44704
        print("  stall limit        %.1f mph power off at the Mission 2 weight, "
              "on lift / %.2f (true stall about %.1f mph)"
              % (args.max_stall_mph, cfg.build.stall_lift_safety_factor,
                 args.max_stall_mph / cfg.build.stall_lift_safety_factor ** 0.5))
        _sa = cfg.build.stall_check_conditions
        print("  stall checked in   %s" % (
            "the mission's air" if not _sa else "%s air (%.0f C, %.4f kg/m3)"
            % (_sa, CONDITIONS[_sa][1], Environment.at_site(_sa).rho)))

    print("=" * 78)
    print("DBF 2026-27 MISSION OPTIMISER")
    print("=" * 78)
    print("  density altitude   %.0f ft  (%.0f%% of sea-level ISA)"
          % (cfg.env.density_altitude_ft, 100.0 * cfg.env.rho / 1.225))
    print("  air density        %.4f kg/m3 %s"
          % (cfg.env.rho,
             "(measured, overrides elevation and temperature)"
             if args.air_density > 0
             else "(from %.0f m, %.0f C)"
             % (cfg.env.field_elev_m, args.temp
                if args.temp is not None
                else CONDITIONS.get(args.conditions,
                                    CONDITIONS[DEFAULT_CONDITIONS])[1])))
    print("  conditions         %s  [%s]" % (SITE_NAME, args.conditions))
    print("  wind               %.1f m/s" % cfg.env.wind_mps)
    tl = cfg.rules.takeoff_distance_m
    print("  takeoff limit      %s"
          % ("none" if tl == float("inf") else "%.1f m (%.0f ft)"
             % (tl, tl / FT2M)))
    print("  battery limit      %.0f Wh" % cfg.rules.max_watt_hours)
    sp_lim = cfg.rules.max_wingspan_m
    print("  wingspan limit     %s"
          % ("none" if sp_lim == float("inf")
             else "%.3f m (%.0f in)" % (sp_lim, sp_lim / 0.0254)))
    print("  lap straights      %.0f m, turns %.0f deg total"
          % (cfg.course.straight_length_m, cfg.course.total_turn_deg))
    print("  window             %.0f s, Mission 2 needs %d laps"
          % (cfg.rules.mission_window_s, cfg.rules.m2_required_laps))
    print("  scoring            %s-normalised" % cfg.rules.normalisation)

    # ---------------- single design ----------------
    if args.evaluate:
        with open(args.evaluate, encoding="utf-8") as fh:
            d = Design(**json.load(fh))
        # A single design self-normalises to 1.0 by definition, which says
        # nothing.  Score it against the field estimate instead and say so.
        if cfg.rules.normalisation == "self":
            cfg.rules.normalisation = "field"
            print()
            print("  (one design cannot be self-normalised - it would score")
            print("   1.0 on both ratios by definition, so this is scored")
            print("   against the --field-m2 / --field-m3 estimates)")
        print_design_report(d, cfg, dx=1.5, title="DESIGN EVALUATION",
                            weather=True)
        print_active_bounds(d, DesignSpace(), cfg)
        print_bom(d, cfg)
        print_sensitivity(sensitivity(d, cfg, dx=3.0))
        return 0

    # ---------------- search ----------------
    settings = (SearchSettings.quick() if args.quick else
                SearchSettings.thorough() if args.thorough else
                SearchSettings.medium() if args.medium else
                SearchSettings())
    settings.processes = args.processes
    settings.seed = args.seed
    if args.finalists > 0:
        settings.top_for_final = args.finalists
    if args.keep_propulsion > 0:
        settings.keep_propulsion = args.keep_propulsion
        settings.per_motor_keep = max(settings.per_motor_keep,
                                      args.keep_propulsion // 10)
    if args.refine > 0:
        settings.top_for_refine = args.refine
    if args.swap_parents > 0:
        settings.swap_parents = args.swap_parents
        settings.swap_refine = max(settings.swap_refine,
                                   (3 * args.swap_parents) // 2)
    if args.swap_refine > 0:
        settings.swap_refine = args.swap_refine
    settings.fill_cores = args.fill_cores
    from dbf2027.optimizer import default_seed_designs, load_seed_designs
    settings.seed_designs = ([] if args.no_seeds else
                             default_seed_designs()
                             + load_seed_designs(args.seeds))
    space = DesignSpace()
    if args.taper:
        space.taper_bounds = (min(args.taper), max(args.taper))
    if args.tail_vh:
        space.tail_vh_bounds = (max(min(args.tail_vh), 0.0), max(args.tail_vh))
    print("  tail volume Vh     %.2f to %.2f (never below the requirements)"
          % space.tail_vh_bounds)
    if args.winglets:
        space.winglet_frac_bounds = (max(min(args.winglets), 0.0),
                                     max(args.winglets))
    print("  winglets           %s" % ("none" if space.winglet_frac_bounds[1] <= 0
                                       else "%.2f to %.2f of the span"
                                       % space.winglet_frac_bounds))
    if args.ar:
        space.ar_bounds = (max(min(args.ar), 1.0), max(args.ar))
        print("  aspect ratio       %.2f to %.2f" % space.ar_bounds)
    print("  taper              %s" % ("%.2f (fixed)" % space.taper_bounds[0]
                                       if space.taper_bounds[0] == space.taper_bounds[1]
                                       else "%.2f to %.2f" % space.taper_bounds))
    if settings.fill_cores:
        from dbf2027 import optimizer as _opt
        print("  rounds filled      refinements, swap refinements and "
              "finalists in whole rounds of %d workers"
              % _opt._n_processes(settings))
    if args.airfoils:
        from dbf2027.components import AIRFOILS
        wanted = [a.strip() for a in args.airfoils.split(",") if a.strip()]
        unknown = [a for a in wanted if a not in AIRFOILS]
        if unknown:
            print("unknown airfoil(s): %s\nthe catalogue has: %s"
                  % (", ".join(unknown), ", ".join(sorted(AIRFOILS))))
            return 2
        space.airfoils = wanted
        # The coarse sweep samples every propulsion set on every airfoil,
        # so restricting the airfoils would shrink it by the same factor -
        # 990 designs for one airfoil instead of 29,700, of which one was
        # feasible.  Keep the sample budget of a full search instead, so a
        # one-airfoil search explores the rest of the space as widely.
        scale = max(1, int(round(len(AIRFOILS) / float(len(wanted)))))
        settings.lhs_points *= scale
        print("  airfoils           %s only  (%d samples per propulsion "
              "set, x%d, to keep the full search's coverage)"
              % (", ".join(wanted), settings.lhs_points, scale))
    if args.motors:
        from dbf2027.components import MOTORS
        wanted_m = [m.strip() for m in args.motors.split(",") if m.strip()]
        unknown = [m for m in wanted_m if m not in MOTORS]
        if unknown:
            print("unknown motor(s): %s\nthe catalogue has: %s"
                  % (", ".join(unknown), ", ".join(sorted(MOTORS))))
            return 2
        space.motors = wanted_m
        print("  motors             %s only" % ", ".join(wanted_m))
    if args.props:
        from dbf2027.components import PROPS
        wanted_p = [q.strip() for q in args.props.split(",") if q.strip()]
        unknown = [q for q in wanted_p if q not in PROPS]
        if unknown:
            print("unknown propeller(s): %s\nthe catalogue has: %s"
                  % (", ".join(unknown), ", ".join(sorted(PROPS))))
            return 2
        space.props = wanted_p
        print("  propellers         %s only" % ", ".join(wanted_p))
    if args.no_seed_refine:
        settings.seeds_per_count = 0
    if args.packs:
        from dbf2027.components import PACKS
        wanted_b = [q.strip() for q in args.packs.split(",") if q.strip()]
        unknown = [q for q in wanted_b if q not in PACKS]
        if unknown:
            print("unknown pack(s): %s\nthe catalogue has: %s"
                  % (", ".join(unknown), ", ".join(sorted(PACKS))))
            return 2
        space.packs = wanted_b
        print("  battery packs      %s only" % ", ".join(wanted_b))
    if args.motor_counts:
        space.motor_counts = sorted(set(args.motor_counts))
        cfg.rules.max_motors = max(cfg.rules.max_motors,
                                   max(space.motor_counts))
        print("  motor counts       %s" % ", ".join(
            str(n) for n in space.motor_counts))

    t0 = time.time()
    res = optimise(cfg, space, settings, verbose=True)
    dt = time.time() - t0
    print("\nsearch finished in %.1f s (%d designs evaluated, %d propulsion "
          "combinations shortlisted)" % (dt, res.n_evaluated, res.shortlist_size))

    if not res.ranked:
        print("\nNothing feasible.  Loosen the takeoff limit, raise the Wh cap, "
              "or widen the payload bounds.")
        return 1

    print_ranking(res.ranked, n=12,
                  title="TOP DESIGNS (%s-normalised score)" % res.mode)
    by_worst = sorted((r for r in res.ranked if "worst_total" in r),
                      key=lambda r: -r["worst_total"])
    if by_worst:
        print("\nBY THE WORST WEATHER PRESET (the score on each design's "
              "worst day; 0 = cannot fly one)")
        for r in by_worst[:12]:
            d = r.get("design_robust") or r["design"]
            print("  #%-3d worst %.3f  nominal %.3f  %d x %s on %s (%s), "
                  "sensor %.1f lb%s"
                  % (res.ranked.index(r) + 1, r["worst_total"],
                     r["total_mission"], d.n_motors, d.motor, d.prop,
                     d.wing_prop_layout, d.sensor_mass_kg * 2.20462,
                     "" if not r.get("robust_trim_pct") else
                     " (trimmed %d%% to fly every preset)"
                     % r["robust_trim_pct"]))
    print_pareto(res.pareto)
    print_top_configurations(res.ranked, cfg, args.top,
                             "" if args.no_export else args.out)

    if args.no_report:
        best = res.best_row.get("design_robust") or res.best_row["design"]
        print_weights_and_speeds(best, cfg)
        if not args.no_export:
            os.makedirs(args.out, exist_ok=True)
            export_csv(res.ranked, os.path.join(args.out,
                                                "ranked_designs.csv"))
            with open(os.path.join(args.out, "best_design.json"), "w",
                      encoding="utf-8") as fh:
                json.dump(asdict(best), fh, indent=2)
            print("\nwrote %s/ranked_designs.csv, best_design.json"
                  % args.out)
        return 0

    # The recommendation is the best design that can be flown on every
    # weather preset (with a propeller swap where needed, rule 3.2.1e) -
    # not simply the top score, which may be illegal on a cool morning.
    best_row = res.best_row
    # its weather-robust version, where it had to be trimmed to fly every
    # preset
    best = best_row.get("design_robust") or best_row["design"]
    if best_row is not res.ranked[0]:
        print("\nrecommending #%d: the best score on its worst weather "
              "preset (%.3f), where #1 scores %.3f"
              % (res.ranked.index(best_row) + 1,
                 best_row.get("worst_total", 0.0),
                 res.ranked[0].get("worst_total", 0.0)))
    if best_row.get("robust_trim_pct"):
        print("  (its sensor trimmed %d%%, to %.1f lb, so it flies every "
              "weather preset)" % (best_row["robust_trim_pct"],
                                   best.sensor_mass_kg * 2.20462))
    elif not best_row.get("robust"):
        print("\nNo finalist can be flown on every weather preset - the "
              "design below is the top score, flagged 'no' in the wx column.")
    print_design_report(best, cfg, dx=1.5, title="RECOMMENDED DESIGN",
                        ref_m2=res.ref_m2, ref_m3=res.ref_m3,
                        ref_gm=res.ref_gm, weather=True)
    print_active_bounds(best, space, cfg)
    print_bom(best, cfg)
    print_sensitivity(sensitivity(best, cfg, dx=3.0,
                                  ref_m2=res.ref_m2,
                                  ref_m3=res.ref_m3,
                                  ref_gm=res.ref_gm))

    # mission-specialist picks, useful when you know the field
    by_m2 = max(res.ranked, key=lambda r: r["m2_raw_lb_per_s"])
    by_m3 = max(res.ranked, key=lambda r: r["m3_raw_lap_lb"])
    print("\n" + "=" * 78)
    print("SPECIALIST PICKS")
    print("=" * 78)
    print("  best Mission 2 raw score: %.4f lb/s  ->  %s span %.2f m, %s, %d sims"
          % (by_m2["m2_raw_lb_per_s"], by_m2["design"].airfoil,
             by_m2["design"].span_m, by_m2["design"].prop,
             by_m2["design"].n_simulators))
    print("  best Mission 3 raw score: %.2f lap*lb -> %s span %.2f m, %s, sensor %.2f lb"
          % (by_m3["m3_raw_lap_lb"], by_m3["design"].airfoil,
             by_m3["design"].span_m, by_m3["design"].prop,
             by_m3["design"].sensor_mass_kg * 2.20462))

    if not args.no_export:
        os.makedirs(args.out, exist_ok=True)
        export_csv(res.ranked, os.path.join(args.out, "ranked_designs.csv"))
        export_json(res.ranked, os.path.join(args.out, "ranked_designs.json"))
        export_csv(res.pareto, os.path.join(args.out, "pareto_front.csv"))
        export_bom_csv(best, cfg, os.path.join(args.out, "bill_of_materials.csv"))
        with open(os.path.join(args.out, "best_design.json"), "w",
                  encoding="utf-8") as fh:
            json.dump(asdict(best), fh, indent=2)
        print("\nwrote %s/ranked_designs.csv, pareto_front.csv, "
              "bill_of_materials.csv, best_design.json" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
