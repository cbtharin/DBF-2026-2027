"""Staged search over the design space.

The space is mixed integer / categorical / continuous and the objective needs
a flight simulation, so a single global method is the wrong tool.  Six stages
instead, each cheap enough to feed the next:

  1. Propulsion screen  - every motor x prop x count x pack: throws out
     everything that melts a motor, over-Cs a pack or cannot turn a prop,
     then flies each survivor around Mission 2 on a reference airframe of
     its motor count (takeoff, climb, top speed, turns, laps) and rates it
     on that time and whether the energy lasts.  Keeps a diverse shortlist.
  2. Coarse sweep       - shortlist x airfoils x Latin-hypercube over the
     continuous variables, coarse integration step.
  3. Refinement         - bounded Nelder-Mead on the continuous variables
     (span, aspect ratio, sensor mass, Mission 3 sensor fraction, taper -
     held at 1.0, a rectangular wing, for now - and drop height) for the
     best configurations, then an exhaustive sweep of
     the integer simulator count, then a short re-refinement.  Known good
     designs (the saved examples, past winners) are seeded in first, so a
     search never finishes below a design it was handed.
  4. Swap test          - every single-component swap (propeller, airfoil,
     motor with its best propellers, pack) on the best configurations; the
     swaps that beat their parent are refined and join the finalists.
  5. Final scoring      - re-scored, and re-checked at a 4x finer step;
     a finalist on the energy limit has its sensor trimmed until it passes
     both, so the report never shows it failing.
  6. Weather check      - every finalist on all three presets, with a
     propeller swap where needed.

Everything after stage 1 is embarrassingly parallel and runs on a Pool.
"""
from __future__ import annotations

import copy
import glob
import json
import math
import os
import random
import time
from dataclasses import dataclass, field, replace
from multiprocessing import Pool, cpu_count
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from .aircraft import Aircraft, Design
from .components import AIRFOILS, MOTORS, PACKS, PROPS, searchable
from .config import Config, Environment
from .mission import mission_2
from .numerics import lhs, nelder_mead, pareto_front
from .propulsion import PropulsionSystem
from .payload import size_payload
from .scoring import (ScoreCard, apply_references, evaluate_scores,
                      normalise_population, population_references)
from .units import KG2LB


# ==========================================================================
# Search space
# ==========================================================================
@dataclass
class DesignSpace:
    airfoils: List[str] = field(default_factory=lambda: list(AIRFOILS.keys()))
    # Every catalogue motor and pack but the retired ones
    # (components.RETIRED: kept so saved designs load, never searched
    # unless named here).
    motors: List[str] = field(default_factory=lambda: searchable(MOTORS))
    props: List[str] = field(default_factory=lambda: list(PROPS.keys()))
    # Up to four motors (BuildStandard.max_motors, a team cap).
    motor_counts: List[int] = field(default_factory=lambda: [1, 2, 3, 4])
    # Whole COTS packs, picked off the shelf.  Nothing is built from cells.
    packs: List[str] = field(default_factory=lambda: searchable(PACKS))
    # Independent propulsion systems, each carrying one identical pack.
    # Two motors can share one pack (1 system) or have one each (2 systems);
    # they may not be paralleled onto a shared pack pair.
    system_counts: List[int] = field(default_factory=lambda: [1, 2, 3, 4])

    # Upper bound is clamped to the rules span limit in continuous_bounds(),
    # so the search never spends samples on aeroplanes that cannot pass tech
    # inspection.
    span_bounds: Tuple[float, float] = (0.80, 3.20)
    # Taper is a real variable: the lifting-line solve prices its induced
    # drag, and liftingline.stall_factor prices its tip stall - a sharply
    # tapered wing loads its small, low-Reynolds tips hardest and stalls
    # there first, losing CLmax.  That is the term that used to be missing
    # and held the floor at 0.45.  0.20 leaves a tip chord that is still a
    # buildable, spar-carrying section.
    #
    # For now the wing is RECTANGULAR (m101): equal bounds hold taper at
    # 1.0 and Nelder-Mead leaves it out of the simplex.  Set (0.20, 1.00)
    # - run_optimizer.py --taper 0.2 1.0, or the GUI's taper fields - to
    # search it again; the model prices any taper either way.
    taper_bounds: Tuple[float, float] = (1.00, 1.00)
    # Aspect ratio is sampled directly and the chord derived from it, so
    # every sample lands inside the range the aero model is valid over.
    # The lower bound was pinning at 5.0 - with span capped by the rules and
    # the payload heavy, the model wants chord - so it is opened to the
    # bottom of the modelled range: 1.5, where the lifting-line tables stop
    # (m104; it was 3.5).  What keeps a stubby wing honest is the runway
    # attitude (BuildStandard.max_ground_alpha_deg): it makes its lift
    # slowly with angle, so it lifts off and lands faster than its stall.
    ar_bounds: Tuple[float, float] = (1.5, 14.0)
    sensor_bounds: Optional[Tuple[float, float]] = None      # default: rules
    # Fraction of the declared maximum sensor weight flown in Mission 3.
    # 1.0 is "fly the declared maximum"; lower trades weight for laps.
    m3_frac_bounds: Tuple[float, float] = (0.35, 1.00)
    # Winglet height at each tip as a share of the span (m110): 0 is none.
    # Hoerner's end-plate relation holds to about 0.2.
    winglet_frac_bounds: Tuple[float, float] = (0.0, 0.15)
    # Horizontal tail volume chosen (m112): the tail is never smaller than
    # its requirements (m111), but a bigger one trims with less download.
    # Below the requirement the variable has no effect.
    tail_vh_bounds: Tuple[float, float] = (0.30, 0.90)
    # Declared Ground Mission drop height, inches.  Higher scores more but
    # needs thicker foam, which makes the container heavier to fly.
    drop_height_bounds: Tuple[float, float] = (12.0, 60.0)
    taper: float = 1.0

    def continuous_bounds(self, cfg: Config) -> List[Tuple[float, float]]:
        r = cfg.rules
        span_lo, span_hi = self.span_bounds
        span_hi = min(span_hi, r.max_wingspan_m)          # RULE 3.1
        # Only ensure the range is not inverted.  This used to be
        # `min(span_lo, span_hi * 0.5)`, which silently dragged any lower
        # bound above half the span limit back down - harmless while the
        # bounds were hard-coded, but a lie once they became a control
        # the user sets: asking for 1.40 to 1.55 m quietly searched from
        # 0.775 m.
        span_lo = min(span_lo, span_hi)
        sb = self.sensor_bounds or (r.sensor_mass_min_kg, r.sensor_mass_max_kg)
        # RULE 3.1.3a - simulator mass is NOT a variable: it is fixed at
        # the weight of the loaded container.  Its slot in the vector is
        # taken by the Mission 3 sensor fraction, which IS free - rule
        # 3.1.1 lets M3 fly a lighter sensor than the declared maximum, and
        # M3 scores laps x weight, so there is a real optimum in between.
        return [(span_lo, span_hi), self.ar_bounds, sb, self.m3_frac_bounds,
                self.taper_bounds, self.drop_height_bounds,
                self.winglet_frac_bounds, self.tail_vh_bounds]

    VAR_NAMES = ("span_m", "aspect_ratio", "sensor_mass_kg", "m3_sensor_frac",
                 "taper", "drop_height_in", "winglet_frac", "tail_vol_h")


@dataclass
class SearchSettings:
    keep_propulsion: int = 110
    per_motor_keep: int = 12
    lhs_points: int = 9
    top_for_refine: int = 36
    top_for_final: int = 15
    nm_iters: int = 90
    # With RK4 plus event location the mission answer is independent of dx
    # to ~1e-5 s over a 32x range (see grid_convergence.py), so all three
    # stages use the same step and it is chosen for speed, not accuracy.
    # Staging different step sizes used to bias which designs survived
    # screening, because the coarse stage had ~20 percent error.
    dx_coarse: float = 12.0
    dx_refine: float = 12.0
    dx_final: float = 12.0
    seed: int = 20262027
    processes: int = 0            # 0 -> cpu_count()-1
    min_pack_wh: float = 18.0
    # Swap test (stage 4): how many of the best refined configurations get
    # every single-component swap tried, and how many of the swaps that
    # beat their parent are then fully refined.
    swap_parents: int = 6
    swap_refine: int = 6
    # Fill the cores (m114): the refinements, the swaps' refinements and
    # the finalists rounded up to whole rounds of the workers, since a
    # stage of long tasks lasts as long as its last round; with fewer
    # configurations than refinement slots, the best are refined from more
    # than one start.  For the cluster, where a part-filled round of 64
    # leaves most of a node idle.
    fill_cores: bool = False
    # Designs every search also refines, ahead of the coarse sweep's best
    # (seed_rows).  None = the saved designs in examples/; [] = none.
    seed_designs: Optional[List[Design]] = None
    # at most this many seeds per motor count, the best of each
    seeds_per_count: int = 6
    # The airframe each motor count's propulsion sets are flown on in the
    # screen (stage 1 sets it; the swap test's prop choice uses it too).
    reference_designs: Optional[Dict[int, Design]] = None

    @staticmethod
    def quick() -> "SearchSettings":
        return SearchSettings(keep_propulsion=40, per_motor_keep=5,
                              lhs_points=6, top_for_refine=12,
                              top_for_final=8, nm_iters=45,
                              swap_parents=3, swap_refine=3)

    @staticmethod
    def medium() -> "SearchSettings":
        """For a search held to one motor count: every refinement slot
        goes to that count, so fewer of them cover it as well."""
        return SearchSettings(keep_propulsion=60, per_motor_keep=8,
                              lhs_points=9, top_for_refine=16,
                              top_for_final=10, nm_iters=90,
                              swap_parents=3, swap_refine=4)

    @staticmethod
    def thorough() -> "SearchSettings":
        return SearchSettings(keep_propulsion=260, per_motor_keep=26,
                              lhs_points=16, top_for_refine=80,
                              top_for_final=25, nm_iters=160,
                              swap_parents=10, swap_refine=15)


# ==========================================================================
# Single-design evaluation (the parallel unit of work)
# ==========================================================================
def quick_reject(ac: Aircraft, cfg: Config) -> Optional[str]:
    """Cheap screen before paying for three mission simulations.

    A full evaluation costs ~100 ms; these checks cost well under one.  With
    payload unbounded by any rule, most of the sampled space cannot fly, so
    making rejection cheap matters more than making the sampler clever.
    Everything here is a necessary condition, so nothing flyable is lost.
    """
    checks = ac.static_checks()
    if not checks["ok"]:
        return "; ".join(checks["failures"])

    t0 = checks["ratings"]["static_thrust_n"]
    if t0 / ac.w_m2_n < 0.18:
        return "static thrust/weight %.2f - cannot climb out and turn" % (
            t0 / ac.w_m2_n)

    aero = ac.aero_clean
    v_s = aero.v_stall(ac.w_m2_n)
    # The propellers make no thrust past their zero-thrust advance ratio at
    # the fastest a full pack can spin them: no level flight above that.
    v_ceiling = ac.prop_sys.speed_ceiling(1.0)
    if v_s * 1.02 >= v_ceiling:
        return ("stall speed %.1f m/s is above the %.1f m/s the propellers "
                "make thrust to" % (v_s, v_ceiling))

    # Energy is the binding constraint once there is no takeoff box.  The
    # screen must be a NECESSARY condition - rejecting something the
    # simulator would have accepted silently deletes good designs - so this
    # is a genuine physical LOWER bound, not an estimate:
    #
    #   * the aeroplane cannot do better than its minimum drag, W/(L/D)max;
    #   * it cannot fly less far than the straight-leg distance (turns only
    #     add arc length);
    #   * it cannot beat a generous best-case total efficiency.
    #
    # An earlier version guessed at full-throttle cruise power instead and
    # over-predicted energy by about 2x, because a real lap spends a large
    # fraction gliding power-off into the turns.  That threw away a third of
    # the flyable designs.
    # Minimum drag in closed form.  For a parabolic polar CD = CD0 + k CL^2,
    # (L/D)max = 1/(2 sqrt(CD0 k)), so D_min = 2 W sqrt(CD0 k).  The real
    # polar here adds camber drag, a stall rise and trim drag on top, so the
    # true minimum drag is always at least this - which is what makes it a
    # valid bound.  It costs one cd0 evaluation instead of a 20-point scan.
    # Drag can be no lower than at the lowest CD0 the aeroplane can reach -
    # at its highest possible speed, skin friction falling with Reynolds
    # number - so this stays a lower bound.
    # The slipstream can take the drag below the power-off minimum (swirl
    # recovery and the blown wing: measured down to 0.92 of it on 1-2 motor
    # designs), so 0.80 of it keeps this a lower bound.
    cd0 = aero.cd0(v_ceiling)
    drag_min = (0.80 * 2.0 * ac.w_m2_n
                * math.sqrt(max(cd0 * aero.k_induced, 1e-12)))
    dist_min = cfg.rules.m2_required_laps * cfg.course.straight_length_m
    e_floor = drag_min * dist_min / 0.62      # 0.62 = optimistic prop*motor*esc
    if e_floor > ac.usable_energy_j():
        return ("Mission 2 needs at least %.0f Wh even flown perfectly, "
                "pack holds %.0f Wh usable"
                % (e_floor / 3600.0, ac.usable_energy_j() / 3600.0))

    # Speed ceiling, also in closed form: a propeller makes no thrust beyond
    # its zero-thrust advance ratio, so V can never exceed J_zero * n * D at
    # the highest rpm the pack can spin it.  Turns only add time on top.
    t_floor = dist_min / max(v_ceiling, 1.0)
    if t_floor > cfg.rules.mission_window_s:
        return ("five laps take at least %.0f s even at the propeller's "
                "speed ceiling, window is %.0f s"
                % (t_floor, cfg.rules.mission_window_s))
    return None


def evaluate_design(design: Design, cfg: Config, dx: float = 6.0,
                    ref_m2: Optional[float] = None,
                    ref_m3: Optional[float] = None,
                    ref_gm: Optional[float] = None,
                    ac: Optional[Aircraft] = None) -> Dict[str, Any]:
    """Score one design; `ac` is its Aircraft if one is already built.

    The stall fit (m120) passes the aeroplane it built for the stall test,
    which is flown instead of being built again."""
    try:
        ac = Aircraft(design, cfg) if ac is None else ac
        why = quick_reject(ac, cfg)
        if why is not None:
            return {"design": design, "feasible": False, "objective": -1e9,
                    "reason": why, "m2_raw_lb_per_s": 0.0,
                    "m3_raw_lap_lb": 0.0, "gm_raw_lb_in": 0.0,
                    "m1": 0.0, "m2": 0.0, "m3": 0.0, "gm": 0.0}
        card: ScoreCard = evaluate_scores(ac, cfg, dx, ref_m2, ref_m3,
                                          ref_gm)
    except Exception as exc:                       # never let one design kill a sweep
        return {"design": design, "feasible": False, "objective": -1e9,
                "reason": "%s: %s" % (type(exc).__name__, exc)}

    row: Dict[str, Any] = {
        "design": design,
        "feasible": card.feasible,
        "objective": card.objective if card.feasible else -1e9,
        "m1": card.m1, "m2": card.m2, "m3": card.m3, "gm": card.gm,
        "total_mission": card.total_mission,
        "m2_raw_lb_per_s": card.m2_raw_lb_per_s,
        "m3_raw_lap_lb": card.m3_raw_lap_lb,
        "gm_raw_lb_in": card.gm_raw_lb_in,
        "reason": card.reason,
    }
    if card.r2 is not None:
        row.update({"m2_time_s": card.r2.time_s, "m2_laps": card.r2.laps,
                    "m2_landing_s": card.r2.landing_time_s,
                    "m2_total_s": card.r2.total_flight_time_s,
                    "m2_takeoff_m": card.r2.takeoff_m,
                    "m2_vmax": card.r2.v_max, "m2_vstall": card.r2.v_stall,
                    "m2_vturn": card.r2.v_turn,
                    "m2_energy_margin": card.r2.energy_margin,
                    "m2_turn_limit": card.r2.turn_limited_by,
                    # the lap's turn (m119)
                    "m2_turn_radius_m": card.r2.turn_radius_m,
                    "m2_turn_load_factor": card.r2.turn_load_factor,
                    "m2_turn_bank_deg": card.r2.turn_bank_deg,
                    "m2_turn_rate_deg_s": card.r2.turn_rate_deg_s})
    if card.r3 is not None:
        row.update({"m3_laps": card.r3.laps, "m3_time_s": card.r3.time_s,
                    "m3_energy_margin": card.r3.energy_margin,
                    "m3_vturn": card.r3.v_turn,
                    "m3_turn_radius_m": card.r3.turn_radius_m,
                    "m3_turn_load_factor": card.r3.turn_load_factor,
                    "m3_turn_bank_deg": card.r3.turn_bank_deg})
    row.update(ac.summary())
    return row


# ==========================================================================
# Weather: a propeller per condition
# ==========================================================================
#: The weather a design must meet, densest first.  A propeller absorbs the
#: most power - and the current and fuse limits bind - in the coolest air.
WEATHER_PRESETS = ("tucson-morning", "tucson", "tucson-hot")
#: Where Mission 2 must complete for a design to count as robust: all
#: three.  With one propeller for every condition the hot afternoon failed
#: almost everything; with a propeller per condition it usually does not,
#: and a failed Mission 2 also locks Mission 3 (rule 3.3.3a), so it is too
#: expensive a day to write off.
MUST_COMPLETE = WEATHER_PRESETS


def _prop_distance(a, b) -> float:
    return (abs(a.diameter_m - b.diameter_m) / b.diameter_m
            + abs(a.pitch_m - b.pitch_m) / b.pitch_m)


def _fly_preset(design: Design, c: Config, dx: float,
                alts: Sequence[str], max_alternatives: int) -> Dict[str, Any]:
    out = _fly_preset_prop(design, c, dx, alts, max_alternatives)
    # The day's peak against each system's fuse, on the propeller chosen
    # for it (RULE 3.2.3.2d) - the static checks already refused any
    # propeller that would blow it; this is what the report shows.
    a = Aircraft(replace(design, prop=out["prop"]), c)
    out["peak_pack_a"] = a.fuse_peak_a()
    out["fuse_a"] = a.prop_sys.fuse_rating_a
    return out


def _fly_preset_prop(design: Design, c: Config, dx: float,
                     alts: Sequence[str],
                     max_alternatives: int) -> Dict[str, Any]:
    ac = Aircraft(design, c)
    chk = ac.static_checks()
    if chk["ok"]:
        r2 = mission_2(ac, c, dx)
        if r2.ok:
            return {"prop": design.prop, "swapped": False, "legal": True,
                    "m2_ok": True, "m2_time_s": r2.time_s, "reason": ""}
        own_reason = r2.reason
    else:
        own_reason = chk["failures"][0]
    best, first_legal, tried = None, None, 0
    for k in alts:
        if tried >= max_alternatives:
            break
        d = replace(design, prop=k)
        a = Aircraft(d, c)
        if not a.static_checks()["ok"]:
            continue
        tried += 1
        first_legal = first_legal or k
        r = mission_2(a, c, dx)
        if r.ok and (best is None or r.time_s < best[1]):
            best = (k, r.time_s)
    why = "own propeller: " + own_reason
    if best:
        return {"prop": best[0], "swapped": True, "legal": True, "m2_ok": True,
                "m2_time_s": best[1], "reason": why}
    if chk["ok"]:
        return {"prop": design.prop, "swapped": False, "legal": True,
                "m2_ok": False, "m2_time_s": None, "reason": own_reason}
    return {"prop": first_legal or design.prop, "swapped": bool(first_legal),
            "legal": bool(first_legal), "m2_ok": False, "m2_time_s": None,
            "reason": why}


def weather_plan(design: Design, cfg: Config, dx: float = 6.0,
                 presets: Sequence[str] = WEATHER_PRESETS,
                 max_alternatives: int = 8) -> Dict[str, Any]:
    """Can this design fly every day it might meet - and on which propeller?

    RULE 3.2.1e - the propeller may be changed for each flight attempt.  A
    design whose propeller over-draws in cool morning air is therefore not
    dead: the team fits a smaller one for that flight.  For each preset
    this flies the design's own propeller if it is legal there and
    completes Mission 2; otherwise it tries the nearest catalogue
    propellers (by diameter and pitch) and keeps the fastest one that
    completes.  Everything else stays as built (rule 3.3.3c).

    Checking every limit at one air density was a real gap: at 28 C the
    search drove nearly every design to 98-99 A, and in the denser morning
    air the same aeroplanes drew 101-104 A - over the 100 A rule.
    """
    own = PROPS[design.prop]
    alts = sorted((k for k in PROPS if k != design.prop),
                  key=lambda k: _prop_distance(PROPS[k], own))
    days: Dict[str, Dict[str, Any]] = {}
    for preset in presets:
        c = copy.deepcopy(cfg)
        c.env = Environment.at_site(preset)
        days[preset] = _fly_preset(design, c, dx, alts, max_alternatives)
    legal_all = all(e["legal"] for e in days.values())
    robust = legal_all and all(days[p]["m2_ok"] for p in MUST_COMPLETE
                               if p in days)
    swaps = [p for p, e in days.items() if e["swapped"]]
    return {"presets": days, "legal_everywhere": legal_all, "robust": robust,
            "all_days": legal_all and all(e["m2_ok"] for e in days.values()),
            "swaps": swaps,
            "code": ("no" if not robust else "swap" if swaps else "ok")}


def _weather_worker(task) -> Dict[str, Any]:
    design, cfg, dx = task
    return weather_plan(design, cfg, dx)


def _worst_preset_one(task) -> float:
    """A finalist's score on its worst weather preset, on the run's own
    scale: each day flown with the propeller its weather plan fits that
    day, and a day it cannot fly counts 0."""
    design, cfg, dx, plan, ref2, ref3, refg = task
    return _worst_preset_score(design, cfg, dx, plan, ref2, ref3, refg)


def _weather_robust_one(task):
    """A finalist made to fly every weather preset: its sensor trimmed 1%
    at a time (to about 18% lighter) until it does, then scored on its
    worst day.  The search sizes each design for the nominal day, right
    on the 300 s Mission 2 window, so the hot day - slower - is usually
    missed; this is the version of it that is not.

    Returns (the design's own weather plan, the robust design or None,
    its trim in percent, its worst-preset score)."""
    design, cfg, dx, ref2, ref3, refg = task
    plan = weather_plan(design, cfg, dx)
    frac = (design.sensor_m3_kg / design.sensor_mass_kg
            if design.sensor_mass_kg > 0 and design.sensor_m3_kg > 0 else 1.0)
    d, p = design, plan
    for k in range(0, 21):
        if k:
            s_ = design.sensor_mass_kg * 0.99 ** k
            d = replace(design, sensor_mass_kg=s_, sensor_m3_kg=s_ * frac)
            p = weather_plan(d, cfg, dx)
        if p["robust"]:
            w = _worst_preset_score(d, cfg, dx, p, ref2, ref3, refg)
            if w > 0.0:
                return plan, d, k, w
    return plan, None, None, 0.0


def _worst_preset_score(design, cfg, dx, plan, ref2, ref3, refg) -> float:
    worst = float("inf")
    for preset, e in plan["presets"].items():
        if not (e.get("legal") and e.get("m2_ok")):
            return 0.0
        c = copy.deepcopy(cfg)
        c.env = Environment.at_site(preset)
        r = evaluate_design(replace(design, prop=e["prop"]), c, dx,
                            ref2, ref3, refg)
        worst = min(worst, r["total_mission"] if r.get("feasible") else 0.0)
    return 0.0 if worst == float("inf") else worst


def _worker(task) -> Dict[str, Any]:
    design, cfg, dx = task[0], task[1], task[2]
    ref2 = task[3] if len(task) > 3 else None
    ref3 = task[4] if len(task) > 4 else None
    refg = task[5] if len(task) > 5 else None
    return evaluate_design(design, cfg, dx, ref2, ref3, refg)


def _worker_fit_stall(task) -> Dict[str, Any]:
    """A coarse sample scored as _worker does, except that one over the
    stall limit (BuildStandard.max_stall_mps) has its sensor brought down -
    the heaviest that meets the limit - and that design is scored instead
    (m108).  The samples land on the stall boundary, where a stall-limited
    optimum sits.

    The limit is found by a secant, not a bisection (m120): the stall goes
    very nearly as the square root of the Mission 2 weight, and the weight
    as the sensor, so v^2 is close to linear in the sensor scale.  A first
    step from that law and a bracketed secant on v^2 - lim^2 (Illinois)
    take 2-4 Aircraft builds where the bisection took 15, and stop closer
    to the limit (1e-6 against 3e-5).  The answer is always a sensor the
    limit passes; a secant that has not converged in 8 builds hands its
    bracket to the bisection.  The aeroplane built for the stall test is
    the one flown, not built again."""
    design, cfg, dx = task[0], task[1], task[2]
    r2, r3, rg = (tuple(task[3:]) + (None, None, None))[:3]

    def built(d: Design) -> Optional[Aircraft]:
        try:
            return Aircraft(d, cfg)
        except Exception:                      # noqa: BLE001
            return None

    a1 = built(design)
    try:
        why = quick_reject(a1, cfg) if a1 is not None else None
    except Exception:                          # noqa: BLE001
        a1 = None
    if a1 is None:
        return _worker(task)                   # evaluate_design says why

    def score(d: Design, ac: Aircraft) -> Dict[str, Any]:
        return evaluate_design(d, cfg, dx, r2, r3, rg, ac=ac)

    if why is None or "mph limit" not in why:
        return score(design, a1)
    lim = cfg.build.max_stall_mps
    s_min = cfg.rules.sensor_mass_min_kg / max(design.sensor_mass_kg, 1e-9)
    if s_min >= 1.0:
        return score(design, a1)               # its rejection

    def at(s: float) -> Design:
        return replace(design, sensor_mass_kg=design.sensor_mass_kg * s,
                       sensor_m3_kg=design.sensor_m3_kg * s)

    def f(ac: Optional[Aircraft]) -> float:
        """v^2 - lim^2: > 0 over the limit (or no aeroplane at all)."""
        if ac is None:
            return float("inf")
        try:
            return ac.stall_speed_factored() ** 2 - lim * lim
        except Exception:                      # noqa: BLE001
            return float("inf")

    # s = 1 is over the limit (that is why we are here); a first step from
    # v^2 proportional to the weight, the sensor carried 1 + n_sim times
    f_hi, hi = f(a1), 1.0
    lo, f_lo, a_lo = None, 0.0, None
    if math.isfinite(f_hi):
        dw_ds = (1 + design.n_simulators) * design.sensor_mass_kg * cfg.env.g
        s1 = 1.0 - a1.w_m2_n * f_hi / (f_hi + lim * lim) / max(dw_ds, 1e-9)
        s1 = min(max(s1, s_min), 1.0 - 1e-9)
        a_s = built(at(s1))
        f1 = f(a_s)
        if f1 > 0.0:
            hi, f_hi = s1, f1
        else:
            lo, f_lo, a_lo = s1, f1, a_s
    if lo is None:
        a_s = built(at(s_min))
        f_min = f(a_s)
        if f_min > 0.0:
            return score(design, a1)           # not even the lightest sensor
        lo, f_lo, a_lo = s_min, f_min, a_s
    tol = (lim * 1e-6) * 2.0 * lim             # 1e-6 of the limit, on v^2
    w_lo, w_hi, side = f_lo, f_hi, 0           # Illinois weights
    for _ in range(8):
        if -f_lo <= tol:
            break
        if math.isfinite(w_hi):
            s = lo + (hi - lo) * (-w_lo) / (w_hi - w_lo)
            s = min(max(s, lo + 1e-12), hi - 1e-12)
        else:
            s = 0.5 * (lo + hi)                # no aeroplane at hi: halve
        a_s = built(at(s))
        fs = f(a_s)
        if fs > 0.0:
            hi, f_hi, w_hi = s, fs, fs
            if side == -1:                     # lo kept twice: halve its weight
                w_lo *= 0.5
            side = -1
        else:
            lo, f_lo, w_lo, a_lo = s, fs, fs, a_s
            if side == 1:
                w_hi *= 0.5
            side = 1
    else:
        if -f_lo > tol:
            # not converged: finish the bracket as the bisection did
            for _ in range(14):
                mid = 0.5 * (lo + hi)
                a_s = built(at(mid))
                if f(a_s) > 0.0:
                    hi = mid
                else:
                    lo, a_lo = mid, a_s
    return score(at(lo), a_lo)


def _indexed_call(job):
    i, func, task = job
    return i, func(task)


def _pool_map(func, tasks, processes: int,
              on_done: Optional[Callable[[int, int], None]] = None,
              min_parallel: int = 8):
    """map() over a process pool, in order.

    With `on_done`, results are collected as they finish (out of order,
    then put back in place) and on_done(done, total) is called after each
    chunk - that is what drives a progress bar.  It runs in THIS process,
    never in a worker.

    Fewer than `min_parallel` tasks run here, serially, because starting a
    pool costs more than a handful of cheap tasks.  Expensive tasks pass 2:
    a refinement takes minutes, and a search restricted to one airfoil can
    have only a few of them - which used to run one after another on one
    core while the rest sat idle.
    """
    n = len(tasks)
    if processes == 1 or n < max(min_parallel, 2):
        out = []
        for i, t in enumerate(tasks):
            out.append(func(t))
            if on_done is not None:
                on_done(i + 1, n)
        return out
    with Pool(processes=processes) as pool:
        chunk = max(1, n // (processes * 8))
        if on_done is None:
            return pool.map(func, tasks, chunksize=chunk)
        out: List[Any] = [None] * n
        done = 0
        for i, r in pool.imap_unordered(
                _indexed_call, [(i, func, t) for i, t in enumerate(tasks)],
                chunksize=chunk):
            out[i] = r
            done += 1
            on_done(done, n)
        return out


#: progress(fraction_done, message) - called from the thread running the
#: search, never from a worker process.  A GUI must not touch widgets from
#: it: hand the two values to its main loop instead.
ProgressFn = Callable[[float, str], None]

_STAGES = (("screen", "propulsion screen", "combinations"),
           ("coarse", "coarse sweep", "designs"),
           ("refine", "local refinement", "configurations"),
           ("swap", "swap test", "swapped designs"),
           ("swap_refine", "refining the best swaps", "configurations"),
           ("final", "final scoring", "designs"),
           ("weather", "weather check", "finalists"))

#: Seconds of one worker's time per task in each stage (refinement: per
#: task per Nelder-Mead iteration), measured on v1.16 searches.  Only the
#: RATIOS matter: they weight the stages so the bar moves at a roughly
#: even pace through all five.
_UNIT_COST = {"screen": 0.0003, "coarse": 0.28, "refine": 13.5,
              "swap": 0.8, "swap_refine": 13.5, "final": 7.0,
              "weather": 10.0}


class _Progress:
    """Turns per-stage counts into one overall fraction.

    `tasks` maps each stage to (number of tasks, cost scale per task).  A
    stage's weight is its estimated wall-clock time: the per-task cost times
    the number of rounds the worker pool needs, ceil(tasks / processes).
    """

    def __init__(self, fn: Optional[ProgressFn],
                 tasks: Dict[str, Tuple[float, float]], processes: int = 1):
        self.fn = fn
        procs = max(int(processes), 1)
        cost = {}
        for k, _, _ in _STAGES:
            n, scale = tasks.get(k, (0.0, 1.0))
            rounds = math.ceil(max(n, 0.0) / procs)
            cost[k] = _UNIT_COST[k] * scale * rounds
        total = sum(cost.values()) or 1.0
        self.share, self.base, acc = {}, {}, 0.0
        for k, _, _ in _STAGES:
            self.base[k] = acc
            self.share[k] = cost[k] / total
            acc += self.share[k]
        self.t0 = time.time()

    def stage(self, key: str) -> Optional[Callable[[int, int], None]]:
        if self.fn is None:
            return None
        idx = [k for k, _, _ in _STAGES].index(key)
        label, unit = _STAGES[idx][1], _STAGES[idx][2]

        def on_done(done: int, total: int) -> None:
            f = self.base[key] + self.share[key] * done / max(total, 1)
            self.fn(min(f, 1.0), "Stage %d of %d, %s: %d of %d %s"
                    % (idx + 1, len(_STAGES), label, done, total, unit))
        self.fn(min(self.base[key], 1.0), "Stage %d of %d, %s: starting"
                % (idx + 1, len(_STAGES), label))
        return on_done

    def finish(self) -> None:
        if self.fn is not None:
            self.fn(1.0, "search finished in %.0f s" % (time.time() - self.t0))


def _n_processes(settings: SearchSettings) -> int:
    if settings.processes > 0:
        return settings.processes
    try:
        return max(1, cpu_count() - 1)
    except NotImplementedError:
        return 1


def _filled(n: int, settings: SearchSettings) -> int:
    """`n` tasks rounded up to whole rounds of the workers when the search
    fills its cores (m114); `n` otherwise."""
    w = _n_processes(settings)
    if not settings.fill_cores or w <= 1 or n <= 0:
        return n
    return -(-n // w) * w


# ==========================================================================
# Stage 1 - propulsion screen
# ==========================================================================
def candidate_packs(space: DesignSpace, cfg: Config,
                    min_wh: float) -> List[Tuple[str, int]]:
    """Legal (pack, number of propulsion systems) combinations.

    RULE 3.2.3.1a, 3.2.3.1c, 3.2.3.2b, 3.2.3.2c.  Every system carries
    one identical pack, so N systems means N times the pack energy, and
    that total is what the 100 Wh limit applies to.  Each individual pack
    must also clear the FAA 100 Wh hand-carry limit.
    """
    r = cfg.rules
    out = []
    for pk_key in space.packs:
        pack = PACKS[pk_key]
        if pack.chemistry not in r.battery_chemistry_allowed:
            continue
        if pack.watt_hours > r.max_pack_watt_hours + 1e-9:
            continue
        for n_sys in space.system_counts:
            total_wh = n_sys * pack.watt_hours
            if total_wh > r.max_watt_hours + 1e-9:
                continue
            if total_wh < min_wh:
                continue
            out.append((pk_key, n_sys))
    return out


_V_MIN_CACHE: Dict[Tuple, float] = {}
_REF_CACHE: Dict[Tuple, Optional[Dict[str, Any]]] = {}


def _min_lap_airspeed(cfg: Config) -> float:
    """The slowest constant airspeed that can fly Mission 2's laps inside
    the window: each straight at V into or with the wind (course headings),
    every turn at the steepest bank and load factor allowed, slowed by the
    pilot efficiency, and no time for the takeoff or climb - a lower bound.
    At a fixed bank a turn takes longer the faster it is flown, so the lap
    time falls and then rises with speed: the limit is the root on the
    falling side."""
    c, b, r, env = cfg.course, cfg.build, cfg.rules, cfg.env
    key = (tuple(c.segments), tuple(c.headings), b.max_bank_deg,
           b.max_turn_load_factor, b.pilot_efficiency, env.wind_mps,
           r.m2_required_laps, r.mission_window_s, env.g)
    if key in _V_MIN_CACHE:
        return _V_MIN_CACHE[key]
    n_max = min(1.0 / math.cos(math.radians(b.max_bank_deg)),
                b.max_turn_load_factor)
    tan_phi = math.sqrt(max(n_max * n_max - 1.0, 1e-9))
    w = env.wind_mps
    turn_rad = sum(math.radians(v) for kind, v in c.segments if kind == "turn")
    legs = [(v, h) for (kind, v), h in zip(c.segments, c.headings)
            if kind == "straight"]

    def lap_time(v: float) -> float:
        t = sum(length / max(v - w * h, 1e-6) for length, h in legs)
        return t + turn_rad * v / (env.g * tan_phi) / max(b.pilot_efficiency, 0.5)

    budget = r.mission_window_s / max(r.m2_required_laps, 1)
    vs = [w + 0.5 + 0.25 * k for k in range(800)]
    v_fast = min(vs, key=lap_time)
    if lap_time(v_fast) > budget:
        _V_MIN_CACHE[key] = v_fast
        return v_fast
    lo, hi = w + 0.5, v_fast
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if lap_time(mid) > budget:
            lo = mid
        else:
            hi = mid
    _V_MIN_CACHE[key] = hi
    return hi


def _min_course_distance(cfg: Config) -> float:
    """The shortest Mission 2 can be flown, m: the laps' straights plus
    every turn at the tightest radius the bank and load-factor limits allow
    at _min_lap_airspeed, R = V^2 / (g tan phi)."""
    c, b = cfg.course, cfg.build
    v = _min_lap_airspeed(cfg)
    n_max = min(1.0 / math.cos(math.radians(b.max_bank_deg)),
                b.max_turn_load_factor)
    tan_phi = math.sqrt(max(n_max * n_max - 1.0, 1e-9))
    r_min = v * v / (cfg.env.g * tan_phi)
    lap = sum(v_ if kind == "straight" else math.radians(v_) * r_min
              for kind, v_ in c.segments)
    return max(cfg.rules.m2_required_laps, 1) * lap


def _reference(design: Design, cfg: Config) -> Optional[Dict[str, Any]]:
    """The reference airframe: its clean aerodynamics, its Mission 2
    weight less its own propulsion, its load limit and gear height - built
    once per worker."""
    key = (design, round(cfg.env.rho, 6), round(cfg.env.wind_mps, 4),
           cfg.build.full_throttle_flights)
    if key not in _REF_CACHE:
        try:
            ac = Aircraft(design, cfg)
            # Everything but the propulsion installation (with its share
            # of the contingency) stays; each candidate brings its own
            # (_installed_mass), worked out on this airframe's geometry.
            ref = dict(aero=ac.aero_clean, n_limit=ac.n_limit,
                       gear_m=ac.gear_height_m, span_m=ac.wing.span_m,
                       fuse_width_m=ac.body.fuse_width_m,
                       prop_axis_m=ac.prop_axis_m,
                       layout=design.wing_prop_layout)
            own = (1.0 + cfg.build.contingency_frac) * sum(
                ac.propulsion_items.values())
            ref["w_rest_n"] = ac.w_m2_n - own * cfg.env.g
            _REF_CACHE[key] = ref
        except Exception:                      # noqa: BLE001
            _REF_CACHE[key] = None
    return _REF_CACHE[key]


def _wing_leads(motor, prop, n_motors: int, ref: Dict[str, Any],
                cfg: Config) -> Tuple[float, float, int]:
    """(lead resistance per motor, ohm; lead mass, kg; wing motors) for
    this motor and propeller on the reference airframe (layout.wing_power_leads)."""
    from .layout import wing_power_leads
    r_l, m_l, lengths = wing_power_leads(
        n_motors, ref["span_m"], ref["fuse_width_m"], motor.i_max_a,
        cfg.build, ref["layout"], prop.diameter_m,
        cfg.rules.wingspan_counts_props)
    return r_l, m_l, len(lengths)


def _installed_mass(ps: PropulsionSystem, ref: Dict[str, Any],
                    cfg: Config) -> float:
    """The propulsion installation's mass, contingency included, kg, as
    the full model would charge it on the reference airframe: the same
    items (aircraft.propulsion_mass_items) with this set's motor, propeller,
    ESC, packs, leads and gear."""
    from .aircraft import (gear_extra_mass_kg, gear_height_for,
                           propulsion_mass_items)
    from .layout import motor_stations
    b, nm = cfg.build, ps.n_motors
    _r, m_leads, n_wing = _wing_leads(ps.motor, ps.prop, nm, ref, cfg)
    gear_x = gear_extra_mass_kg(
        gear_height_for(ps.prop.diameter_m, ref["prop_axis_m"], b), b)
    n_ext = 0
    if ref["layout"] == "stagger":
        _n0, ys, _s = motor_stations(nm, ref["span_m"], ref["fuse_width_m"],
                                     "stagger", ps.prop.diameter_m,
                                     cfg.rules.wingspan_counts_props)
        n_ext = len(ys) // 2 if len(ys) >= 4 else 0
    items = propulsion_mass_items(ps.motor, ps.prop, ps.esc, ps.pack, nm,
                                  ps.n_systems, n_wing, m_leads, gear_x,
                                  n_ext, b)
    return (1.0 + b.contingency_frac) * sum(items.values())


def _flown_merit(ps: PropulsionSystem, ref: Dict[str, Any],
                 cfg: Config) -> Optional[Tuple[float, float, float]]:
    """(merit, top speed, Mission 2 time) of this propulsion on the
    reference airframe, or None if it cannot take off, hold level flight or
    turn.  Everything at full throttle, as the missions are flown."""
    from .mission import takeoff, climb, turn_capability, max_level_speed
    b, r, env = cfg.build, cfg.rules, cfg.env
    aero = ref["aero"]
    w = ref["w_rest_n"] + _installed_mass(ps, ref, cfg) * env.g
    mass = w / env.g
    to = takeoff(aero, ps, mass, w, cfg, soc=1.0, gear_height_m=ref["gear_m"])
    if not to.ok:
        return None
    soc = 0.7
    cl = climb(aero, ps, mass, w, cfg, soc, to.v_liftoff)
    if cl is None:
        return None
    t_cl, e_cl, _v_cl, _x = cl
    v_top = max_level_speed(aero, ps, w, soc)
    t_top = ps.thrust(v_top, soc, 1.0)
    if (t_top < aero.drag(v_top, w, t_top) * 0.999
            or v_top <= env.wind_mps + 0.5):
        return None
    tc = turn_capability(aero, ps, w, soc, cfg, ref["n_limit"])
    if tc.turn_rate_rad_s <= 1e-3:
        return None
    c = cfg.course
    turn_rad = sum(math.radians(v) for kind, v in c.segments if kind == "turn")
    t_turns = turn_rad / tc.turn_rate_rad_s / max(b.pilot_efficiency, 0.5)
    t_str = sum(v / (v_top - env.wind_mps * h)
                for (kind, v), h in zip(c.segments, c.headings)
                if kind == "straight")
    t_laps = max(r.m2_required_laps, 1) * (t_str + t_turns)
    t_m2 = to.time_s + t_cl + t_laps
    p_full = ps.solve(v_top, soc, 1.0).p_elec_w
    e_need = to.energy_j + e_cl + p_full * t_laps
    e_have = b.land_at_discharge * ps.total_watt_hours * 3600.0
    merit = 1000.0 / t_m2
    merit *= min(1.0, e_have / max(e_need, 1.0)) ** 2
    merit *= min(1.0, r.mission_window_s / t_m2) ** 4
    return merit, v_top, t_m2


def _screen_one(task) -> Optional[Dict[str, Any]]:
    # Fresh shared tables for each set, as for each Aircraft (m99): the
    # screen flies many sets on one cached reference airframe.
    from .liftingline import clear_caches
    from .propulsion import clear_system_cache
    clear_caches()
    clear_system_cache()
    mk, pk, nm, pack_key, n_sys, cfg = task[:6]
    ref_design = task[6] if len(task) > 6 else None
    if nm % n_sys != 0:          # motors must divide evenly across systems
        return None
    motor, prop = MOTORS[mk], PROPS[pk]
    pack = PACKS[pack_key]
    if pack.series > motor.max_cells:      # over the motor's voltage rating
        return None
    # Wing propellers side by side must fit the largest legal span, each in
    # its share of it with 3 cm clear a side (layout.motor_stations, the
    # overlap check in Aircraft.static_checks).  A necessary condition -
    # it allows no fuselage at all - so it only removes what can never fly.
    # (staggered props may overlap their neighbours by up to 30% of a
    # diameter: from the fuselage to the tip that fits a larger prop)
    per_side = nm // 2
    if per_side:
        half = 0.5 * cfg.rules.max_wingspan_m
        d_mid = half / per_side - 0.06
        d_stag = (half - 0.03) / (1.0 + 0.7 * (per_side - 1))
        if prop.diameter_m > max(d_mid, d_stag):
            return None
    # On a reference airframe each wing motor has its own power leads,
    # whose resistance the power solve carries, as the full model's does.
    ref = _reference(ref_design, cfg) if ref_design is not None else None
    r_lead = (_wing_leads(motor, prop, nm, ref, cfg)[0]
              if ref is not None else 0.0)
    ps = PropulsionSystem(motor, prop, nm, pack, cfg.env,
                          n_systems=n_sys,
                          r_harness_ohm=cfg.build.harness_resistance_ohm,
                          max_fuse_a=cfg.rules.max_fuse_current_a,
                          esc_airflow_mps=cfg.build.esc_airflow_mps,
                          mission_s=cfg.rules.mission_window_s,
                          r_lead_ohm=r_lead,
                          i0_speed_share=cfg.build.i0_speed_share)
    rc = ps.rating_check()
    if (rc["motor_current_margin"] < 1.0 or rc["motor_power_margin"] < 1.0
            or rc["esc_current_margin"] < 1.0 or rc["pack_current_margin"] < 1.0
            or rc["tip_mach"] > 0.82 or rc["prop_rpm_margin"] < 1.0):
        return None
    # RULE 3.2.3.2d - each system's blade fuse (100 A at most) is a hard
    # stop, checked at the peak current over every airspeed (rating_check).
    # There is no rule on the total across systems.
    if rc["fuse_margin"] < 1.0:
        return None
    if rc["static_thrust_n"] < 4.0:
        return None
    # No lap speed is assumed.  The set is flown on the reference airframe
    # of its motor count (_flown_merit): takeoff, climb, top speed where
    # thrust meets drag, the fastest sustained turn, Mission 2's laps - and
    # rated on that Mission 2 time and whether the energy lasts.  (It was
    # rated on its thrust at exactly 26 and 34 m/s, speeds picked by hand.)
    v_rated, t_m2 = 0.0, 0.0
    if ref is not None:
        fm = _flown_merit(ps, ref, cfg)
        fom = fm[0] if fm else 0.0
        if fm:
            v_rated, t_m2 = fm[1], fm[2]
        p_full = ps.solve(max(v_rated, _min_lap_airspeed(cfg)), 0.70,
                          1.0).p_elec_w
    else:
        # No reference airframe (no seeds): thrust power at the set's own
        # best speed above the slowest airspeed that could finish in the
        # window.
        v = v_lo = _min_lap_airspeed(cfg)
        op_lo = ps.solve(v_lo, 0.70, 1.0)
        if op_lo.thrust_total_n <= 0.5:
            return None
        best_op, best_tv = op_lo, op_lo.thrust_total_n * v_lo
        v_rated = v_lo
        for _ in range(60):
            v *= 1.08
            op = ps.solve(v, 0.70, 1.0)
            tv = op.thrust_total_n * v
            if tv > best_tv:
                best_op, best_tv, v_rated = op, tv, v
            elif op.thrust_total_n <= 0.0 or tv < 0.85 * best_tv:
                break
        fom = ((rc["static_thrust_n"] ** 0.35) * (best_tv ** 0.65)
               / (ps.mass_kg ** 0.30))
        p_full = best_op.p_elec_w
    t_full = (cfg.build.land_at_discharge * ps.total_watt_hours * 3600.0
              / max(p_full, 1.0))
    if ref is None and cfg.build.full_throttle_flights:
        fom *= min(1.0, t_full / cfg.rules.mission_window_s) ** 2
    return {"motor": mk, "prop": pk, "n_motors": nm, "pack": pack_key,
            "n_systems": n_sys, "fom": fom, "full_throttle_s": t_full,
            "v_rated": v_rated, "m2_time_s_ref": t_m2,
            "static_thrust_n": rc["static_thrust_n"],
            "mass_kg": ps.mass_kg, "wh": ps.total_watt_hours}


def _fit_taper(d: Design, space: DesignSpace) -> Design:
    """A saved design with its taper put inside the search's bounds (a
    tapered seed in a rectangular-wing search flies rectangular)."""
    lo, hi = space.taper_bounds
    t = min(max(d.taper, lo), hi)
    return d if t == d.taper else replace(d, taper=t)


def _reference_designs(space: DesignSpace, cfg: Config,
                       settings: SearchSettings,
                       verbose: bool = True) -> Dict[int, Design]:
    """The airframe each motor count's sets are flown on in the screen:
    the best seed of that count, made flyable and scored as stage 3 scores
    it; a count with no seed borrows the best seed's airframe with its
    motor count changed."""
    seeds = [_fit_taper(d, space) for d in
             (settings.seed_designs if settings.seed_designs is not None
              else default_seed_designs())]
    seen, usable = set(), []
    for d in seeds:
        if (d.motor in MOTORS and d.prop in PROPS and d.pack in PACKS
                and d.airfoil in AIRFOILS and d.discrete_key() not in seen):
            seen.add(d.discrete_key())
            usable.append(d)
    if not usable:
        if verbose:
            print("  no seeds: sets rated on thrust power at their own best "
                  "speed")
        return {}
    rows = _pool_map(_seed_one, [(d, cfg, settings.dx_coarse, None, None,
                                  None) for d in usable],
                     _n_processes(settings), None, min_parallel=2)
    rows = [r for r in rows if r]
    if not rows:
        return {}
    rows.sort(key=lambda r: -r["objective"])
    out: Dict[int, Design] = {}
    for nm in space.motor_counts:
        best = next((r for r in rows if r["design"].n_motors == nm), None)
        if best is not None:
            out[nm] = best["design"]
        else:
            d = rows[0]["design"]
            out[nm] = replace(d, n_motors=nm, n_systems=1,
                              wing_prop_layout=("mid" if nm > 1 else
                                                d.wing_prop_layout))
    if verbose:
        for nm, d in sorted(out.items()):
            print("  %d-motor sets flown on: %s, span %.2f m, chord %.3f m, "
                  "sensor %.1f lb (the best seed: %s on %s)"
                  % (nm, d.airfoil, d.span_m, d.chord_m,
                     d.sensor_mass_kg * 2.20462, d.motor, d.prop))
    return out


def screen_propulsion(space: DesignSpace, cfg: Config,
                      settings: SearchSettings,
                      verbose: bool = True,
                      on_done: Optional[Callable[[int, int], None]] = None
                      ) -> List[Dict[str, Any]]:
    packs = candidate_packs(space, cfg, settings.min_pack_wh)
    refs = _reference_designs(space, cfg, settings, verbose)
    settings.reference_designs = refs
    tasks = [(mk, pk, nm, pack_key, n_sys, cfg, refs.get(nm))
             for mk in space.motors
             for pk in space.props
             for nm in space.motor_counts
             for (pack_key, n_sys) in packs
             if nm % n_sys == 0]
    if verbose:
        print("  screening %d propulsion combinations "
              "(%d pack/system options)..." % (len(tasks), len(packs)))
    raw = _pool_map(_screen_one, tasks, _n_processes(settings), on_done)
    survivors = [r for r in raw if r]
    survivors.sort(key=lambda r: -r["fom"])
    keep = _shortlist(survivors, settings, space.motor_counts)
    if verbose:
        print("  %d passed ratings, shortlisted %d" % (len(survivors), len(keep)))
    return keep


def _shortlist(survivors: List[Dict[str, Any]], settings: SearchSettings,
               counts: Sequence[int]) -> List[Dict[str, Any]]:
    """The propulsion shortlist from the screened setups, best first.

    Every motor count the search allows gets an equal share of it first
    (the best of that count, at most per_motor_keep per motor): ranked on
    one merit, the setups of one count can take every place - a 260-place
    shortlist over 1-12 motors held none with more than four, so no
    aeroplane was ever designed around them.  Then the best few per motor,
    then the best of what is left."""
    limit = settings.keep_propulsion
    keep: List[Dict[str, Any]] = []
    taken = set()
    per_motor: Dict[str, int] = {}

    def take(r):
        keep.append(r)
        taken.add(id(r))
        per_motor[r["motor"]] = per_motor.get(r["motor"], 0) + 1

    counts = sorted(set(counts))
    if len(counts) > 1:
        share = max(1, limit // len(counts))
        for nm in counts:
            got = 0
            for r in survivors:
                if got >= share:
                    break
                if (r["n_motors"] == nm and id(r) not in taken
                        and per_motor.get(r["motor"], 0)
                        < settings.per_motor_keep):
                    take(r)
                    got += 1
    for r in survivors:
        if len(keep) >= limit:
            break
        if (id(r) not in taken
                and per_motor.get(r["motor"], 0) < settings.per_motor_keep):
            take(r)
    for r in survivors:
        if len(keep) >= limit:
            break
        if id(r) not in taken:
            take(r)
    return keep


def _configs_by_count(rows: List[Dict[str, Any]], limit: int,
                      counts: Sequence[int],
                      flyable_only: bool = False) -> List[Dict[str, Any]]:
    """One row per configuration, best first (_dedupe_configs) - but each
    motor count in `counts` first gets its share of the `limit`, so the
    refinement and the final scoring always include every count's best,
    not just the leading count's.  A count's share goes to its designs
    that fly first; with `flyable_only` (the finals) only to those - a
    place given to a design that cannot fly is a place wasted."""
    counts = sorted(set(counts))
    if len(counts) <= 1:
        return _dedupe_configs(rows, limit)
    share = max(1, limit // len(counts))
    seen, out = set(), []
    for nm in counts:
        got = 0
        for want in ((True,) if flyable_only else (True, False)):
            for r in rows:
                if got >= share or len(out) >= limit:
                    break
                d = r["design"]
                k = d.discrete_key()
                if (d.n_motors == nm and k not in seen
                        and bool(r.get("feasible", True)) == want):
                    seen.add(k)
                    out.append(r)
                    got += 1
    for r in rows:
        if len(out) >= limit:
            break
        k = r["design"].discrete_key()
        if k not in seen:
            seen.add(k)
            out.append(r)
    out.sort(key=lambda r: -r["objective"])
    return out


# ==========================================================================
# Stage 2/3 helpers
# ==========================================================================
# ==========================================================================
# Seeds: known good designs, refined alongside the sampled ones
# ==========================================================================
EXAMPLES_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")

#: How far a seed's sensor is stepped down, 3% at a time, to make a design
#: saved under older rules or an older model fly again (to ~40% of it).
_SEED_SCALES = tuple(0.97 ** k for k in range(31))


def load_seed_designs(paths: Sequence[str]) -> List[Design]:
    """Designs from JSON files, or from every *.json in a folder.  Anything
    that is not a saved design (a ranked list, a parts list) is skipped."""
    files: List[str] = []
    for p in paths:
        if os.path.isdir(p):
            files += sorted(glob.glob(os.path.join(p, "*.json")))
        elif os.path.isfile(p):
            files.append(p)
    out: List[Design] = []
    for f in files:
        try:
            with open(f, encoding="utf-8") as fh:
                data = json.load(fh)
            out.append(Design(**data))
        except (OSError, ValueError, TypeError):
            continue
    return out


def default_seed_designs() -> List[Design]:
    """The saved designs in examples/ - the recommended design among them."""
    return load_seed_designs([EXAMPLES_DIR])


def _seed_usable(d: Design, space: DesignSpace) -> bool:
    """A seed's parts must exist, and a restricted search must not have
    ruled them out: each is in the search space, or retired (a seed may
    name a retired part - that is what keeps a saved design loadable)."""
    from .components import RETIRED
    for key, allowed, store in ((d.motor, space.motors, MOTORS),
                                (d.prop, space.props, PROPS),
                                (d.pack, space.packs, PACKS),
                                (d.airfoil, space.airfoils, AIRFOILS)):
        if key not in store:
            return False
        if key not in allowed and key not in RETIRED:
            return False
    return (d.n_motors in space.motor_counts
            and d.n_systems in space.system_counts)


def _seed_one(task) -> Optional[Dict[str, Any]]:
    """A seed scored as the search scores it - its sensor stepped down 3%
    at a time until it flies under today's rules and model, since a design
    saved under older ones may now be too heavy.  None if it never does."""
    design, cfg, dx, r2, r3, rg = task
    frac = (design.sensor_m3_kg / design.sensor_mass_kg
            if design.sensor_mass_kg > 0 and design.sensor_m3_kg > 0 else 1.0)
    for f in _SEED_SCALES:
        s = design.sensor_mass_kg * f
        row = evaluate_design(replace(design, sensor_mass_kg=s,
                                      sensor_m3_kg=s * frac),
                              cfg, dx, r2, r3, rg)
        if row.get("feasible"):
            row["seeded"] = True
            return row
    return None


def seed_rows(space: DesignSpace, cfg: Config, settings: SearchSettings,
              refs: Tuple[Optional[float], Optional[float], Optional[float]],
              verbose: bool = True) -> List[Dict[str, Any]]:
    """The known good designs, made flyable and scored - best first.

    Random sampling decides which configurations get refined, and a
    configuration none of whose samples flies is never refined at all.
    That is how three searches running missed the single-motor design that
    beat all of their winners.  A seed takes a refinement slot whatever
    its samples did, so the search starts from the best it already knows.
    """
    if settings.seeds_per_count <= 0:
        # --no-seed-refine: the seeds gave the reference airframes only
        if verbose:
            print("  known designs not refined (--no-seed-refine)")
        return []
    seeds = [_fit_taper(d, space) for d in
             (settings.seed_designs if settings.seed_designs is not None
              else default_seed_designs())]
    seen, usable = set(), []
    for d in seeds:
        if not _seed_usable(d, space) or d.discrete_key() in seen:
            continue
        seen.add(d.discrete_key())
        usable.append(d)
    if not usable:
        return []
    rows = _pool_map(_seed_one, [(d, cfg, settings.dx_coarse, refs[0],
                                  refs[1], refs[2]) for d in usable],
                     _n_processes(settings), None, min_parallel=2)
    rows = [r for r in rows if r]
    rows.sort(key=lambda r: -r["objective"])
    # The same number of seeds for every motor count: the saved designs
    # are mostly single-motor winners of earlier searches, and seeding one
    # count with forty polished designs and another with none would refine
    # the first far harder.
    per: Dict[int, int] = {}
    kept = []
    for r in rows:
        nm = r["design"].n_motors
        if per.get(nm, 0) < settings.seeds_per_count:
            per[nm] = per.get(nm, 0) + 1
            kept.append(r)
    rows = kept
    if verbose:
        print("  seeded with %d known designs, %d of them flying (sensors "
              "trimmed where today's rules needed it)" % (len(usable),
                                                          len(rows)))
    return rows


def _layouts_for(n_motors: int) -> List[str]:
    """The wing-propeller layouts open to a motor count: none to choose
    for a nose motor alone; mid or tip with one a side; staggered too
    with two or more a side (layout.motor_stations)."""
    per_side = int(n_motors) // 2
    if per_side == 0:
        return ["mid"]
    return ["mid", "tip"] + (["stagger"] if per_side >= 2 else [])


def _make_design(base: Dict[str, Any], airfoil: str, x: Sequence[float],
                 n_sim: int, taper: float) -> Design:
    """x = (span, AR, declared sensor mass, M3 sensor fraction, taper,
    drop height, winglet height / span, horizontal tail volume)."""
    span = max(x[0], 0.3)
    ar = min(max(x[1], 1.0), 15.5)
    taper = min(max(x[4], 0.20), 1.0) if len(x) > 4 else taper
    sensor = max(x[2], 1e-3)
    frac = min(max(x[3], 0.05), 1.0) if len(x) > 3 else 1.0
    drop = min(max(x[5], 1.0), 60.0) if len(x) > 5 else 60.0
    return Design(span_m=span, chord_m=span / ar, airfoil=airfoil,
                  motor=base["motor"], prop=base["prop"],
                  n_motors=base["n_motors"], pack=base["pack"],
                  n_systems=base["n_systems"],
                  sensor_mass_kg=sensor, n_simulators=n_sim,
                  sensor_m3_kg=sensor * frac,
                  drop_height_in=round(drop),      # RULE 3.3.4, whole in
                  taper=taper,
                  wing_prop_layout=base.get("layout", "mid"),
                  winglet_h_m=(max(x[6], 0.0) * span if len(x) > 6
                               else base.get("winglet_h_m", 0.0)),
                  tail_vol_h=(max(x[7], 0.0) if len(x) > 7
                              else base.get("tail_vol_h", 0.0)))


def _plausible_sim_count(base: Dict[str, Any], span: float, chord: float,
                         sensor_kg: float, sim_kg: float,
                         rng: random.Random, cfg: Config) -> int:
    """How many simulators could this power system plausibly lift?

    RULE 3.1.3a: `sim_kg` is no longer free - every simulator weighs the
    same as the loaded container, so a heavy sensor makes every extra box
    heavy too.  That coupling is why the count has to be sampled against
    the sensor actually chosen rather than independently.

    Sampling the count independently of the aeroplane wastes almost the whole
    budget on designs that cannot leave the ground.  Bound it by requiring a
    static thrust-to-weight of at least 0.30, subtract a rough empty mass,
    and sample a fraction of what is left.
    """
    t0 = base.get("static_thrust_n", 0.0)
    m_thrust = (t0 / 0.30) / 9.80665
    # Energy is the binding limit now that there is no takeoff box: five laps
    # of cruise on a 100 Wh pack.  W = E * (L/D) * eta / (V * t).
    usable_j = base.get("wh", 100.0) * 3600.0 * 0.78
    # W = E (L/D) eta / d over the course: d is the shortest Mission 2 can
    # be flown - every turn at the tightest radius the bank and load-factor
    # limits allow at the slowest speed that still finishes in the window
    # (_min_course_distance), so this stays an upper bound on the mass.
    m_energy = usable_j * 9.0 * 0.60 / (9.80665 * _min_course_distance(cfg))
    m_empty = 1.4 + 3.0 * (span * chord) + base.get("mass_kg", 1.0)
    # RULE 3.2.1a - 55 lb gross is a hard ceiling, so it is a bound on how
    # much cargo is worth sampling at all.
    # the loaded container - sensor and its container - flies too
    m_cap = min(m_thrust, m_energy, cfg.rules.max_gross_mass_kg)
    if cfg.build.max_stall_mps < float("inf"):
        # what the wing could lift at the stall limit on a CL of 3, more than
        # any built wing makes (m107)
        m_cap = min(m_cap, 0.5 * cfg.env.rho * cfg.build.max_stall_mps ** 2
                    * span * chord * 3.0
                    / (cfg.build.stall_lift_safety_factor * cfg.env.g))
    capacity = m_cap - m_empty - sensor_kg - cfg.rules.container_tare_kg
    if capacity <= 0.0:
        return 0
    n = int(round(rng.random() * capacity / max(sim_kg, 0.05)))
    return max(0, min(n, cfg.rules.max_simulators))


def _fit_sample(x: Sequence[float], airfoil, bounds, cfg: Config) -> List[float]:
    """A coarse sample moved inside two limits it could not otherwise meet
    (m107): the chord cap (aspect ratio at least span / max_chord_m) and the
    stall limit (the sensor scaled into the range a wing of this area could
    possibly lift at max_stall_mps - the whole aeroplane's mass on 1.5 x the
    section's clmax, a bound no built wing reaches).  Necessary conditions
    only, so nothing that could fly is lost; with neither limit set the
    sample comes back unchanged."""
    b = cfg.build
    x = list(x)
    if b.max_chord_m < float("inf"):
        x[1] = max(x[1], x[0] / b.max_chord_m)
    if b.max_stall_mps < float("inf"):
        lo, hi = bounds[2]
        area = x[0] * x[0] / x[1]
        m_max = (0.5 * cfg.env.rho * b.max_stall_mps ** 2 * area
                 * 1.5 * airfoil.cl_max_2d
                 / (b.stall_lift_safety_factor * cfg.env.g))
        top = min(hi, max(m_max, lo))
        if hi > lo:
            x[2] = lo + (x[2] - lo) / (hi - lo) * (top - lo)
    return x


def coarse_sweep(shortlist: List[Dict[str, Any]], space: DesignSpace,
                 cfg: Config, settings: SearchSettings,
                 verbose: bool = True,
                 refs: Tuple[Optional[float], Optional[float], Optional[float]] = (None, None, None),
                 on_done: Optional[Callable[[int, int], None]] = None
                 ) -> List[Dict[str, Any]]:
    rng = random.Random(settings.seed)
    bounds = space.continuous_bounds(cfg)
    max_sim = cfg.rules.max_simulators
    tasks: List[Tuple[Design, Config, float]] = []
    for base in shortlist:
        pts = lhs(settings.lhs_points, bounds, rng)
        lays = _layouts_for(base["n_motors"])
        for af in space.airfoils:
            for kx, x in enumerate(pts):
                x = _fit_sample(x, AIRFOILS[af], bounds, cfg)
                span, ar = x[0], x[1]
                # each simulator weighs one loaded container (rule 3.1.3a)
                n_sim = _plausible_sim_count(
                    base, span, span / ar, x[2],
                    x[2] + cfg.rules.container_tare_kg, rng, cfg)
                # the wing-propeller layouts take the samples in turn
                b_l = dict(base, layout=lays[kx % len(lays)])
                tasks.append((_make_design(b_l, af, x, n_sim, space.taper),
                              cfg, settings.dx_coarse, refs[0], refs[1],
                              refs[2]))
    if verbose:
        print("  coarse sweep: %d designs..." % len(tasks))
    work = (_worker_fit_stall if cfg.build.max_stall_mps < float("inf")
            else _worker)
    rows = _pool_map(work, tasks, _n_processes(settings), on_done)
    rows = [r for r in rows if r.get("feasible")]
    rows.sort(key=lambda r: -r["objective"])
    if verbose:
        print("  %d feasible%s" % (len(rows), " (samples over the stall limit "
                                   "had their sensor brought down to it)"
                                   if work is _worker_fit_stall else ""))
    return rows


def _polish_to_bounds(negobj, x: List[float], lower: List[float],
                      upper: List[float], passes: int = 2) -> List[float]:
    """Try each variable at its bounds, keep what helps.

    Nelder-Mead is a local method, so where it lands depends on where it
    started: the same configuration refined from two different coarse
    seeds converged to span 1.60 m and to 1.83 m, and only the second is
    right.  This is the cheap insurance - a dozen evaluations per pass
    against ninety for the simplex - and it catches exactly the failure
    that matters here, a variable that wanted to sit at its limit and was
    left in the middle.

    It only ever accepts an improvement, so it cannot make a design
    worse.
    """
    best = list(x)
    best_f = negobj(best)
    for _ in range(passes):
        improved = False
        for i in range(len(best)):
            for edge in (lower[i], upper[i]):
                if abs(best[i] - edge) < 1e-12:
                    continue
                cand = list(best)
                cand[i] = edge
                f = negobj(cand)
                if f < best_f - 1e-12:
                    best_f, best, improved = f, cand, True
        if not improved:
            break
    return best


def _refine_one(task) -> Dict[str, Any]:
    row, cfg, space, settings, ref2, ref3, refg = task
    d: Design = row["design"]
    bounds = space.continuous_bounds(cfg)
    lower = [b[0] for b in bounds]
    upper = [b[1] for b in bounds]
    base = {"motor": d.motor, "prop": d.prop, "n_motors": d.n_motors,
            "pack": d.pack, "n_systems": d.n_systems,
            "layout": d.wing_prop_layout}

    def negobj(x: List[float], n_sim: int) -> float:
        des = _make_design(base, d.airfoil, x, n_sim, d.taper)
        return -evaluate_design(des, cfg, settings.dx_refine,
                                ref2, ref3, refg)["objective"]

    frac = (d.sensor_m3_kg / d.sensor_mass_kg
            if d.sensor_mass_kg > 0 and d.sensor_m3_kg > 0 else 1.0)
    x0 = [d.span_m, d.span_m / d.chord_m, d.sensor_mass_kg, frac,
          d.taper, d.drop_height_in, d.winglet_h_m / max(d.span_m, 1e-9),
          d.tail_vol_h]
    # The initial simplex is scaled to each variable's RANGE, not written
    # out as absolute numbers.  Hard-coded steps looked reasonable per
    # axis but were not comparable across them: 0.12 kg of sensor mass
    # against a 0.15-20 kg range is a 0.6 percent probe, so the simplex
    # had to expand many times to reach the heavy end and, inside the
    # iteration budget, usually did not.  From a typical seed this cost
    # about 6 percent of objective - the search returned a 6.7 kg sensor
    # where 9.2 kg scored better - and it was worse than simply tripling
    # the iteration count, which is the signature of a badly shaped
    # simplex rather than too little effort.
    step = [max(0.08 * (hi - lo), 1e-6)
            for lo, hi in zip(lower, upper)]

    n_sim = d.n_simulators
    x, f = nelder_mead(lambda v: negobj(v, n_sim), x0, step, lower, upper,
                       maxiter=settings.nm_iters)

    # Integer sweep on the simulator count with the refined geometry frozen.
    # The objective is genuinely non-smooth here: the bay repacks in discrete
    # steps as boxes are added, so a scan beats any gradient method.  Walk
    # upward and stop once the aeroplane stops being able to fly.
    best_n, best_f = n_sim, f
    misses = 0
    for k in range(0, cfg.rules.max_simulators + 1):
        if k == n_sim:
            continue
        fk = negobj(x, k)
        if fk < best_f:
            best_f, best_n = fk, k
            misses = 0
        elif fk > 1e8:              # infeasible
            misses += 1
            if misses >= 3 and k > best_n:
                break
    if best_n != n_sim:
        x, best_f = nelder_mead(lambda v: negobj(v, best_n), x, step,
                                lower, upper, maxiter=settings.nm_iters // 2)

    # Then check whether any variable simply wanted to be at a limit.
    x = _polish_to_bounds(lambda v: negobj(v, best_n), x, lower, upper)

    des = _make_design(base, d.airfoil, x, best_n, d.taper)
    out = evaluate_design(des, cfg, settings.dx_refine, ref2, ref3, refg)
    # Never hand back something worse than the seed.  Nelder-Mead with
    # clipped bounds is not guaranteed to improve on its starting point,
    # and a refinement stage that can make a design worse turns "more
    # effort" into "a worse answer" - which is exactly what it looked
    # like from the outside.
    if out.get("objective", -1e9) < row.get("objective", -1e9):
        return row
    return out


def refine(rows: List[Dict[str, Any]], space: DesignSpace, cfg: Config,
           settings: SearchSettings, verbose: bool = True,
           refs: Tuple[Optional[float], Optional[float], Optional[float]] = (None, None, None),
           on_done: Optional[Callable[[int, int], None]] = None
           ) -> List[Dict[str, Any]]:
    limit = _filled(settings.top_for_refine, settings)
    seeds = _configs_by_count(rows, limit, space.motor_counts)
    extra: List[Dict[str, Any]] = []
    if settings.fill_cores and len(seeds) < limit:
        # too few configurations for the slots: the best get a second
        # start, filling the round they need and never adding one (m114)
        extra = _extra_starts(rows, seeds,
                              _filled(len(seeds), settings) - len(seeds))
    if verbose:
        print("  refining %d configurations%s..."
              % (len(seeds), ", and %d more start%s of the best of them"
                 % (len(extra), "s" if len(extra) > 1 else "")
                 if extra else ""))
    tasks = [(r, cfg, space, settings, refs[0], refs[1], refs[2])
             for r in seeds + extra]
    out = _pool_map(_refine_one, tasks, _n_processes(settings), on_done,
                    min_parallel=2)
    out = [r for r in out if r.get("feasible")]
    out.sort(key=lambda r: -r["objective"])
    return out


# ==========================================================================
# Stage 4 - swap test
# ==========================================================================
#: Sensor masses each swapped design is scored at, as fractions of its
#: parent's: a better component usually carries a little more.
_SWAP_SCAN = (0.97, 1.0, 1.03, 1.06)


def _best_props(mk: str, nm: int, pack_key: str, n_sys: int, cfg: Config,
                props: Sequence[str], k: int,
                cache: Dict[Tuple, List[str]],
                ref: Optional[Design] = None) -> List[str]:
    """The k propellers the propulsion screen rates best for this motor
    and pack - the same figure of merit stage 1 ranks on, flown on the same
    reference airframe (`ref`, SearchSettings.reference_designs)."""
    key = (mk, nm, pack_key, n_sys, k)
    if key not in cache:
        scored = []
        for pk in props:
            r = _screen_one((mk, pk, nm, pack_key, n_sys, cfg, ref))
            if r:
                scored.append((r["fom"], pk))
        scored.sort(reverse=True)
        cache[key] = [pk for _, pk in scored[:k]]
    return cache[key]


def _prefetch_best_props(designs: Sequence[Design], space: DesignSpace,
                         cfg: Config, settings: SearchSettings,
                         cache: Dict[Tuple, List[str]]) -> int:
    """Fill _best_props' cache for every question _swap_neighbours will
    ask about these designs, flying the ratings on the worker pool.

    The same screen flights and the same sort as _best_props, so the
    neighbours come out unchanged; done one at a time in the main process
    they were ~20,000 flights and a quarter of the swap test's time on 64
    cores, every worker idle (m100).  Returns the number of flights.
    """
    want: Dict[Tuple, int] = {}            # cache key -> motor count
    for d in designs:
        motor, pack = MOTORS[d.motor], PACKS[d.pack]
        for mk in space.motors:
            if mk == d.motor or pack.series > MOTORS[mk].max_cells:
                continue
            want[(mk, d.n_motors, d.pack, d.n_systems, 2)] = d.n_motors
        for pack_key, n_sys in candidate_packs(space, cfg, settings.min_pack_wh):
            if ((pack_key, n_sys) == (d.pack, d.n_systems)
                    or d.n_motors % n_sys
                    or PACKS[pack_key].series > motor.max_cells):
                continue
            want[(d.motor, d.n_motors, pack_key, n_sys, 1)] = d.n_motors
    want = {k: v for k, v in want.items() if k not in cache}
    sets = sorted({k[:4] for k in want})   # (motor, motors, pack, systems)
    refs = settings.reference_designs or {}
    tasks = [(mk, pk, nm, pack_key, n_sys, cfg, refs.get(nm))
             for mk, nm, pack_key, n_sys in sets for pk in space.props]
    rated = _pool_map(_screen_one, tasks, _n_processes(settings))
    scored: Dict[Tuple, List[Tuple[float, str]]] = {s: [] for s in sets}
    for t, r in zip(tasks, rated):
        if r:
            scored[(t[0], t[2], t[3], t[4])].append((r["fom"], t[1]))
    for key in want:
        best = sorted(scored[key[:4]], reverse=True)
        cache[key] = [pk for _, pk in best[:key[4]]]
    return len(tasks)


def _swap_neighbours(d: Design, space: DesignSpace, cfg: Config,
                     settings: SearchSettings,
                     cache: Dict[Tuple, List[str]]) -> List[Design]:
    """Every design one component away: another propeller, another
    airfoil, another motor (with the two propellers that suit it best),
    another pack or system count (with its own and the best propeller)."""
    out: List[Design] = []
    motor, pack = MOTORS[d.motor], PACKS[d.pack]
    ref = (settings.reference_designs or {}).get(d.n_motors)
    for pk in space.props:
        if pk != d.prop:
            out.append(replace(d, prop=pk))
    for af in space.airfoils:
        if af != d.airfoil:
            out.append(replace(d, airfoil=af))
    for lay in _layouts_for(d.n_motors):
        if lay != d.wing_prop_layout:
            out.append(replace(d, wing_prop_layout=lay))
    for mk in space.motors:
        if mk == d.motor or pack.series > MOTORS[mk].max_cells:
            continue
        for pk in _best_props(mk, d.n_motors, d.pack, d.n_systems, cfg,
                              space.props, 2, cache, ref):
            out.append(replace(d, motor=mk, prop=pk))
    for pack_key, n_sys in candidate_packs(space, cfg, settings.min_pack_wh):
        if ((pack_key, n_sys) == (d.pack, d.n_systems)
                or d.n_motors % n_sys
                or PACKS[pack_key].series > motor.max_cells):
            continue
        props = {d.prop} | set(_best_props(d.motor, d.n_motors, pack_key,
                                           n_sys, cfg, space.props, 1, cache, ref))
        for pk in sorted(props):
            out.append(replace(d, pack=pack_key, n_systems=n_sys, prop=pk))
    return out


def swap_test(parents: List[Dict[str, Any]], space: DesignSpace, cfg: Config,
              settings: SearchSettings, refs, verbose: bool = True,
              on_scan=None, make_on_refine=None) -> List[Dict[str, Any]]:
    """Try every single-component swap on the best refined designs.

    Refinement polishes the continuous variables of a configuration; it
    never changes the motor, propeller, pack or airfoil.  So a design can
    be one swap away from a better one and never find it - which is what
    a manual swap test turned up.  This does it automatically: each
    swapped design is scored at a few sensor masses around its parent's,
    and the swaps that beat their parent are refined in full and join the
    finalists.
    """
    cache: Dict[Tuple, List[str]] = {}
    # every propeller rating the neighbours need, on the pool (m100)
    _prefetch_best_props([row["design"] for row in parents], space, cfg,
                         settings, cache)
    tasks, owner = [], []
    for pi, row in enumerate(parents):
        d = row["design"]
        frac = (d.sensor_m3_kg / d.sensor_mass_kg
                if d.sensor_mass_kg > 0 and d.sensor_m3_kg > 0 else 1.0)
        for nb in _swap_neighbours(d, space, cfg, settings, cache):
            for f in _SWAP_SCAN:
                s_ = d.sensor_mass_kg * f
                tasks.append((replace(nb, sensor_mass_kg=s_,
                                      sensor_m3_kg=s_ * frac),
                              cfg, settings.dx_refine, refs[0], refs[1],
                              refs[2]))
                owner.append(pi)
    if verbose:
        print("  %d swapped designs around %d finalists (%d scorings)..."
              % (len(tasks) // len(_SWAP_SCAN), len(parents), len(tasks)))
    rows = _pool_map(_worker, tasks, _n_processes(settings), on_scan)
    # the best scoring of each swapped configuration, kept only if it beats
    # the parent it came from
    best: Dict[Tuple, Dict[str, Any]] = {}
    for pi, r in zip(owner, rows):
        if not r.get("feasible"):
            continue
        if r["objective"] <= parents[pi].get("objective", -1e9) + 1e-9:
            continue
        k = r["design"].discrete_key()
        if k not in best or r["objective"] > best[k]["objective"]:
            best[k] = r
    have = {p["design"].discrete_key() for p in parents}
    winners = _configs_by_count(
        sorted((r for k, r in best.items() if k not in have),
               key=lambda r: -r["objective"]),
        _filled(settings.swap_refine, settings), space.motor_counts)
    if verbose:
        print("  %d swaps beat their parent; refining the best %d"
              % (len(best), len(winners)))
    if not winners:
        return []
    # the refinement's progress stage starts only now that it is reached
    on_refine = make_on_refine() if make_on_refine else None
    out = _pool_map(_refine_one, [(r, cfg, space, settings, refs[0], refs[1],
                                   refs[2]) for r in winners],
                    _n_processes(settings), on_refine, min_parallel=2)
    return [r for r in out if r.get("feasible")]


# ==========================================================================
# Stage 5 - final scoring, hardened against the energy limit
# ==========================================================================
#: The finer steps a finalist must hold at, as divisors of the final step
#: (8: the step the report flies).
_FINE_STEPS = (4.0, 8.0)
#: The most 0.25% trims of the Mission 3 sensor spent holding a lap (~10%).
_M3_LAP_TRIMS = 40


def _final_one(task) -> Dict[str, Any]:
    """Score a finalist at the final step, and make sure it holds at four
    and eight times finer ones (eight: the step the report flies).

    Designs sit on the energy limit, where the last 0.01 Wh between two
    integration steps can decide whether Mission 2 finishes - and whether
    Mission 3 completes its last lap.  At each check:
      * Mission 2 fails at any step: the declared sensor is trimmed 0.25
        percent (the Mission 3 sensor with it);
      * Mission 3 flies a lap short at a finer step: only the Mission 3
        sensor is trimmed 0.25 percent, until the lap holds at every step
        (m103) - five laps at a hair less weight score far more than four,
        and Mission 2 is untouched.  If _M3_LAP_TRIMS trims cannot hold
        it, the design is scored as the finest step flies it.
    So the report never shows the recommended design failing, or flying a
    lap fewer than its score counted.
    """
    # (A finalist that fails at the final step itself is trimmed the same
    # way rather than dropped: refined right onto the energy limit, it can
    # come back a hair over it - the final step is not the step it was
    # refined at - and dropping it threw away a third of every run's
    # finalists, every count's best among them.  Before m99 the shared
    # blown-wing tables could also round a jet differently in another
    # worker; every aircraft now starts from empty tables.)
    design, cfg, dx, ref2, ref3, refg = task
    d = design
    r: Dict[str, Any] = {}
    m3_trims = 0
    for i in range(24 + _M3_LAP_TRIMS):
        r = evaluate_design(d, cfg, dx, ref2, ref3, refg)
        if r.get("feasible"):
            fine = [evaluate_design(d, cfg, dx / k, ref2, ref3, refg)
                    for k in _FINE_STEPS]
            if all(f.get("feasible") for f in fine):
                laps = r.get("m3_laps") or 0
                if all((f.get("m3_laps") or 0) >= laps for f in fine):
                    r["hardened_steps"] = i
                    r["m3_lap_trims"] = m3_trims
                    return r
                if m3_trims < _M3_LAP_TRIMS:
                    m3_trims += 1
                    d = replace(d, sensor_m3_kg=d.sensor_m3_kg * 0.9975)
                    continue
                r = fine[-1]                 # scored as the finest step flies it
                r["hardened_steps"] = i
                r["m3_lap_trims"] = m3_trims
                r["scored_at_finest_step"] = True
                return r
        frac = (d.sensor_m3_kg / d.sensor_mass_kg
                if d.sensor_mass_kg > 0 and d.sensor_m3_kg > 0 else 1.0)
        s_ = d.sensor_mass_kg * 0.9975
        d = replace(d, sensor_mass_kg=s_, sensor_m3_kg=s_ * frac)
    return r


def _dedupe_configs(rows: List[Dict[str, Any]], limit: int) -> List[Dict[str, Any]]:
    """One seed per discrete configuration, best first."""
    seen = set()
    out = []
    for r in rows:
        k = r["design"].discrete_key()
        if k in seen:
            continue
        seen.add(k)
        out.append(r)
        if len(out) >= limit:
            break
    return out


def _extra_starts(rows: List[Dict[str, Any]], chosen: List[Dict[str, Any]],
                  k: int) -> List[Dict[str, Any]]:
    """Up to `k` more refinement starts for the configurations `chosen`
    (m114): each one's next-best sample in `rows`, best configuration
    first, then each one's third-best, and so on.  Nelder-Mead is local, so
    a second start is a second chance at the configuration's best."""
    seen = {r["design"].key() for r in chosen}
    spare: Dict[Tuple, List[Dict[str, Any]]] = {}
    for r in rows:
        if r["design"].key() not in seen:
            seen.add(r["design"].key())
            spare.setdefault(r["design"].discrete_key(), []).append(r)
    order = [r["design"].discrete_key() for r in chosen]
    out: List[Dict[str, Any]] = []
    depth = 0
    while len(out) < k:
        more = False
        for key in order:
            q = spare.get(key, ())
            if depth < len(q):
                out.append(q[depth])
                more = True
                if len(out) >= k:
                    break
        if not more:
            break
        depth += 1
    return out


def _seeded_refine_limit(seeds: List[Dict[str, Any]],
                         coarse: List[Dict[str, Any]], limit: int) -> int:
    """Refinement slots once the seeds are in: the sampled configurations'
    usual `limit`, plus one per seed whose configuration is not already
    among them - a seed adds a slot, it never takes one from the search."""
    sampled = {r["design"].discrete_key()
               for r in _dedupe_configs(coarse, limit)}
    return limit + len({r["design"].discrete_key() for r in seeds} - sampled)


# ==========================================================================
# Driver
# ==========================================================================
@dataclass
class OptimisationResult:
    ranked: List[Dict[str, Any]]
    pareto: List[Dict[str, Any]]
    shortlist_size: int
    n_evaluated: int
    ref_m2: float = 0.0          # normaliser actually used in the report
    ref_m3: float = 0.0
    ref_gm: float = 0.0
    mode: str = "self"

    @property
    def best_row(self) -> Optional[Dict[str, Any]]:
        """The best design on its worst weather preset (worst_total) -
        the design that scores most whatever day it meets, not the top
        nominal score, which may not fly a hot afternoon.

        Falls back to the best robust design, then to the top of the
        ranking, so a result is never empty - but a non-robust row is
        flagged "no" in the ranking's weather column.
        """
        scored = [r for r in self.ranked if r.get("worst_total", 0.0) > 0.0]
        if scored:
            return max(scored, key=lambda r: r["worst_total"])
        for r in self.ranked:
            if r.get("robust"):
                return r
        return self.ranked[0] if self.ranked else None


def optimise(cfg: Config, space: Optional[DesignSpace] = None,
             settings: Optional[SearchSettings] = None,
             verbose: bool = True,
             progress: Optional[ProgressFn] = None) -> OptimisationResult:
    """The staged search.  `progress(fraction, message)`, if given, is
    called as work completes - see ProgressFn."""
    space = space or DesignSpace()
    settings = settings or SearchSettings()
    n_eval = 0
    n_packs = len(candidate_packs(space, cfg, settings.min_pack_wh))
    track = _Progress(progress, {
        "screen": (len(space.motors) * len(space.props)
                   * len(space.motor_counts) * n_packs, 1.0),
        "coarse": (settings.keep_propulsion * settings.lhs_points
                   * len(space.airfoils), 1.0),
        "refine": (_filled(settings.top_for_refine, settings),
                   settings.nm_iters),
        "swap": (settings.swap_parents * 160 * len(_SWAP_SCAN), 1.0),
        "swap_refine": (_filled(settings.swap_refine, settings),
                        settings.nm_iters),
        "final": (_filled(settings.top_for_final, settings), 1.0),
        "weather": (_filled(settings.top_for_final, settings), 1.0)},
        _n_processes(settings))
    t_stage = time.time()

    def stage_time() -> None:
        nonlocal t_stage
        if verbose:
            print("  (%.0f s)" % (time.time() - t_stage))
        t_stage = time.time()

    if verbose:
        print("Stage 1: propulsion screen")
    shortlist = screen_propulsion(space, cfg, settings, verbose,
                                  track.stage("screen"))
    stage_time()
    if not shortlist:
        raise RuntimeError("no propulsion combination passed its ratings - "
                           "loosen the catalogue or the Wh limit")

    if verbose:
        print("Stage 2: coarse design sweep")
    refs = (cfg.rules.field_best_m2_lb_per_s,
            cfg.rules.field_best_m3_lap_lb,
            cfg.rules.field_best_gm_lb_in)
    coarse = coarse_sweep(shortlist, space, cfg, settings, verbose, refs,
                          track.stage("coarse"))
    stage_time()
    n_eval += len(shortlist) * settings.lhs_points * len(space.airfoils)
    if not coarse and verbose:
        # the seeds may still fly: a limit near the edge of what flies can
        # leave every random sample out and a known design in (m121)
        print("  nothing flew in the coarse sweep: the search goes on from "
              "the seed designs")

    # In "self" mode the normaliser is the population, so update it from
    # the coarse sweep before refining.  Refining against a stale reference
    # would weight M2 against M3 wrongly.
    if cfg.rules.normalisation == "self" and coarse:
        refs = population_references(coarse)
        # Put the coarse rows on the NEW scale before anything compares
        # them with refined rows.  Without this the pool mixes objectives
        # computed against two different sets of references: the coarse
        # ones keep the generous field-estimate denominators and outrank
        # refined designs that are better on every raw score.  That is
        # how a 1.60 m span coarse sample beat the 1.83 m refinement of
        # itself - it was not a better aeroplane, it was scored on a
        # different scale.
        apply_references(coarse, refs[0], refs[1], refs[2],
                         cfg.rules.normalisation)
        coarse.sort(key=lambda r: -r["objective"])
        if verbose:
            print("  self-normalising against the coarse best: "
                  "M2 %.4f lb/s, M3 %.2f lap*lb, GM %.0f lb*in" % refs)

    if verbose:
        print("Stage 3: local refinement")
    # Known good designs first, so they always get a refinement slot
    # (refine keeps one row per configuration, in this order) - on top of
    # the sampled configurations' own slots, not in place of any of them.
    seeds = seed_rows(space, cfg, settings, refs, verbose)
    n_eval += len(seeds)
    if not coarse:
        if not seeds:
            raise RuntimeError("no feasible design in the coarse sweep or among "
                               "the seed designs - check the takeoff distance, "
                               "Wh limit, stall limit and payload bounds")
        # the seeds are the population now: normalise against them, as the
        # coarse sweep would have been (m121)
        if cfg.rules.normalisation == "self":
            refs = population_references(seeds)
            apply_references(seeds, refs[0], refs[1], refs[2],
                             cfg.rules.normalisation)
            seeds.sort(key=lambda r: -r["objective"])
            if verbose:
                print("  self-normalising against the seeds' best: "
                      "M2 %.4f lb/s, M3 %.2f lap*lb, GM %.0f lb*in" % refs)
    refined = refine(seeds + coarse, space, cfg,
                     replace(settings, top_for_refine=_seeded_refine_limit(
                         seeds, coarse, settings.top_for_refine)),
                     verbose, refs, track.stage("refine"))
    stage_time()
    n_eval += len(refined) * (settings.nm_iters + cfg.rules.max_simulators)

    pool_rows = refined + seeds + coarse[: settings.top_for_refine]
    pool_rows.sort(key=lambda r: -r["objective"])

    # Stage 4: every single-component swap on the best configurations; the
    # swaps that beat their parent are refined and join the pool.
    if verbose:
        print("Stage 4: swap test")
    # every count's leaders get their swaps tried, not just the leading
    # count's (at least one parent per count)
    parents = _configs_by_count(
        pool_rows, max(settings.swap_parents, len(set(space.motor_counts))),
        space.motor_counts, flyable_only=True)
    swapped = swap_test(parents, space, cfg, settings, refs, verbose,
                        track.stage("swap"),
                        lambda: track.stage("swap_refine"))
    stage_time()
    n_eval += len(parents) * 160 * len(_SWAP_SCAN)
    pool_rows = pool_rows + swapped
    pool_rows.sort(key=lambda r: -r["objective"])
    finalists = _configs_by_count(pool_rows,
                                  _filled(settings.top_for_final, settings),
                                  space.motor_counts, flyable_only=True)

    if verbose:
        print("Stage 5: final high-resolution scoring of %d designs "
              "(checked again at 4x and 8x finer steps)" % len(finalists))
    tasks = [(r["design"], cfg, settings.dx_final, refs[0], refs[1], refs[2])
             for r in finalists]
    final_rows = _pool_map(_final_one, tasks, _n_processes(settings),
                           track.stage("final"), min_parallel=2)
    stage_time()
    dropped = [r for r in final_rows if not r.get("feasible")]
    if verbose and dropped:
        print("  %d finalists do not fly at the final step:" % len(dropped))
        for r in dropped[:16]:
            d_ = r.get("design")
            print("    %s: %s" % (
                ("%d x %s, %s" % (d_.n_motors, d_.motor, d_.prop))
                if d_ is not None else "?", str(r.get("reason", ""))[:90]))
    final_rows = [r for r in final_rows if r.get("feasible")]
    held = [r for r in final_rows if r.get("m3_lap_trims")]
    if verbose and held:
        print("  %d finalists flew Mission 3 a lap short at a finer step; "
              "their Mission 3 sensor was trimmed until the lap held "
              "(at most %.1f%%)" % (len(held), 100 * (1 - 0.9975 ** max(
                  r["m3_lap_trims"] for r in held))))
    final_rows.sort(key=lambda r: -r["objective"])
    n_eval += len(tasks)

    # Final scores are normalised against the finished population, so the
    # best design sits at exactly 1.0 on each ratio and everything else is
    # reported as a fraction of it.
    if cfg.rules.normalisation == "self":
        ref2, ref3, refg = normalise_population(final_rows, "self")
    else:
        ref2, ref3, refg = refs
    final_rows.sort(key=lambda r: -r["objective"])

    # Stage 6: every finalist flown on each weather preset, swapping the
    # propeller where its own over-draws or fails (rule 3.2.1e).
    if verbose:
        print("Stage 6: weather check of %d finalists (%s)"
              % (len(final_rows), ", ".join(WEATHER_PRESETS)))
    # Each finalist's own weather plan, and its weather-robust version:
    # trimmed until it flies every preset, and scored on its worst day on
    # this run's scale (_weather_robust_one).
    out = _pool_map(_weather_robust_one,
                    # the same resolution the finalists were scored at:
                    # they sit on the energy limit, where a coarser step
                    # can flip a pass into a fail
                    [(r["design"], cfg, settings.dx_final, ref2, ref3, refg)
                     for r in final_rows],
                    _n_processes(settings), track.stage("weather"),
                    min_parallel=2)
    stage_time()
    for r, (plan, d_rob, trim, w) in zip(final_rows, out):
        r["weather"] = plan
        r["robust"] = plan["robust"]
        r["wx"] = plan["code"]
        r["worst_total"] = w
        r["robust_trim_pct"] = trim
        if d_rob is not None:
            r["design_robust"] = d_rob
            r["robust_sensor_kg"] = d_rob.sensor_mass_kg

    pts = [(r["m2_raw_lb_per_s"], r["m3_raw_lap_lb"]) for r in final_rows]
    front_idx = pareto_front(pts) if pts else []
    pareto = [final_rows[i] for i in front_idx]
    pareto.sort(key=lambda r: -r["m2_raw_lb_per_s"])

    track.finish()
    return OptimisationResult(final_rows, pareto, len(shortlist), n_eval,
                              ref2, ref3, refg, cfg.rules.normalisation)


# ==========================================================================
# Sensitivity
# ==========================================================================
def sensitivity(design: Design, cfg: Config, dx: float = 6.0,
                ref_m2: Optional[float] = None,
                ref_m3: Optional[float] = None,
                ref_gm: Optional[float] = None
                ) -> List[Tuple[str, float, float, float]]:
    """One-at-a-time +/-10 percent sweep on the continuous variables.

    Returns (variable, low value, high value, d(objective)) rows so you can see
    which dimension is actually worth arguing about in the design review.

    `ref_m2`/`ref_m3`/`ref_gm` must be the SAME references the run
    reported its score against.  Without them a self-normalised run falls
    back to the default field estimates, and because self mode does not cap
    the ratio at 1 the table prints scores above the 7.5 ceiling - numbers
    that cannot be compared with anything else in the report.
    """
    base = evaluate_design(design, cfg, dx, ref_m2, ref_m3,
                           ref_gm)["objective"]
    out = []
    fields = [("span_m", design.span_m), ("chord_m", design.chord_m),
              ("sensor_mass_kg", design.sensor_mass_kg),
              ("sensor_m3_kg", design.sensor_m3_kg or design.sensor_mass_kg),
              ("drop_height_in", design.drop_height_in)]
    for name, val in fields:
        lo_d = _replace(design, name, val * 0.90)
        hi_d = _replace(design, name, val * 1.10)
        if name == "sensor_mass_kg" and design.sensor_m3_kg:
            # Keep the Mission 3 fraction, as the search does: a lighter
            # declared sensor with the M3 sensor left alone is illegal
            # (rule 3.1.1), so the step would always read "infeasible".
            lo_d = _replace(lo_d, "sensor_m3_kg", design.sensor_m3_kg * 0.90)
            hi_d = _replace(hi_d, "sensor_m3_kg", design.sensor_m3_kg * 1.10)
        lo = evaluate_design(lo_d, cfg, dx, ref_m2, ref_m3,
                             ref_gm)["objective"]
        hi = evaluate_design(hi_d, cfg, dx, ref_m2, ref_m3,
                             ref_gm)["objective"]
        lo = lo if lo > -1e8 else float("nan")
        hi = hi if hi > -1e8 else float("nan")
        out.append((name, lo, hi, base))
    # integer: simulator count
    for delta in (-1, +1):
        k = max(0, min(cfg.rules.max_simulators, design.n_simulators + delta))
        if k == design.n_simulators:
            continue
        d2 = _replace(design, "n_simulators", k)
        v = evaluate_design(d2, cfg, dx, ref_m2, ref_m3,
                            ref_gm)["objective"]
        out.append(("n_simulators=%d" % k, v, v, base))
    return out


def _replace(design: Design, field_name: str, value) -> Design:
    from dataclasses import replace as dc_replace
    return dc_replace(design, **{field_name: value})
