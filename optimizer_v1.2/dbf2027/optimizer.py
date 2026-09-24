"""Staged search over the design space.

The space is mixed integer / categorical / continuous and the objective needs
a flight simulation, so a single global method is the wrong tool.  Four stages
instead, each cheap enough to feed the next:

  1. Propulsion screen  - every motor x prop x count x pack, two operating
     points each, no lap simulation.  Throws out everything that melts a
     motor, over-Cs a pack or cannot turn a prop.  Keeps a diverse shortlist.
  2. Coarse sweep       - shortlist x airfoils x Latin-hypercube over the
     continuous variables, coarse integration step.
  3. Refinement         - bounded Nelder-Mead on (span, chord, sensor mass,
     simulator mass) for the best configurations, then an exhaustive sweep of
     the integer simulator count, then a short re-refinement.
  4. Final scoring      - fine integration step, Pareto front, report.

Everything after stage 1 is embarrassingly parallel and runs on a Pool.
"""
from __future__ import annotations

import math
import os
import random
from dataclasses import dataclass, field
from multiprocessing import Pool, cpu_count
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .aircraft import Aircraft, Design
from .components import AIRFOILS, MOTORS, PACKS, PROPS
from .config import Config
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
    motors: List[str] = field(default_factory=lambda: list(MOTORS.keys()))
    props: List[str] = field(default_factory=lambda: list(PROPS.keys()))
    # Two motors maximum.
    motor_counts: List[int] = field(default_factory=lambda: [1, 2])
    # Whole COTS packs, picked off the shelf.  Nothing is built from cells.
    packs: List[str] = field(default_factory=lambda: list(PACKS.keys()))
    # Independent propulsion systems, each carrying one identical pack.
    # Two motors can share one pack (1 system) or have one each (2 systems);
    # they may not be paralleled onto a shared pack pair.
    system_counts: List[int] = field(default_factory=lambda: [1, 2])

    # Upper bound is clamped to the rules span limit in continuous_bounds(),
    # so the search never spends samples on aeroplanes that cannot pass tech
    # inspection.
    span_bounds: Tuple[float, float] = (0.80, 3.20)
    # Taper is a real variable now that the lifting-line solve prices it.
    # Floored at 0.45: the model has no tip-stall penalty, and a sharply
    # tapered wing tip-stalls, so letting it run to 0.2 would be the
    # optimiser exploiting a missing term.
    taper_bounds: Tuple[float, float] = (0.45, 1.00)
    # Aspect ratio is sampled directly and the chord derived from it, so
    # every sample lands inside the range the aero model is valid over.
    # The lower bound was pinning at 5.0 - with span capped by the rules and
    # the payload heavy, the model wants chord - so it is opened to the
    # bottom of the modelled range.
    ar_bounds: Tuple[float, float] = (3.5, 14.0)
    sensor_bounds: Optional[Tuple[float, float]] = None      # default: rules
    # Fraction of the declared maximum sensor weight flown in Mission 3.
    # 1.0 is "fly the declared maximum"; lower trades weight for laps.
    m3_frac_bounds: Tuple[float, float] = (0.35, 1.00)
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
                self.taper_bounds, self.drop_height_bounds]

    VAR_NAMES = ("span_m", "aspect_ratio", "sensor_mass_kg", "m3_sensor_frac",
                 "taper", "drop_height_in")


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

    @staticmethod
    def quick() -> "SearchSettings":
        return SearchSettings(keep_propulsion=40, per_motor_keep=5,
                              lhs_points=6, top_for_refine=12,
                              top_for_final=8, nm_iters=45)

    @staticmethod
    def thorough() -> "SearchSettings":
        return SearchSettings(keep_propulsion=260, per_motor_keep=26,
                              lhs_points=16, top_for_refine=80,
                              top_for_final=25, nm_iters=160)


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
    if v_s > 32.0:
        return "stall speed %.0f m/s at Mission 2 weight - unflyable" % v_s

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
    cd0 = aero.cd0(max(1.3 * v_s, 12.0))
    drag_min = 2.0 * ac.w_m2_n * math.sqrt(max(cd0 * aero.k_induced, 1e-12))
    dist_min = cfg.rules.m2_required_laps * cfg.course.straight_length_m
    e_floor = drag_min * dist_min / 0.62      # 0.62 = optimistic prop*motor*esc
    if e_floor > ac.usable_energy_j():
        return ("Mission 2 needs at least %.0f Wh even flown perfectly, "
                "pack holds %.0f Wh usable"
                % (e_floor / 3600.0, ac.usable_energy_j() / 3600.0))

    # Speed ceiling, also in closed form: a propeller makes no thrust beyond
    # its zero-thrust advance ratio, so V can never exceed J_zero * n * D at
    # the highest rpm the pack can spin it.  Turns only add time on top.
    n_ceiling = ac.motor.kv * ac.pack.v_full / 60.0
    v_ceiling = ac.prop.j_zero * n_ceiling * ac.prop.diameter_m
    t_floor = dist_min / max(v_ceiling, 1.0)
    if t_floor > cfg.rules.mission_window_s:
        return ("five laps take at least %.0f s even at the propeller's "
                "speed ceiling, window is %.0f s"
                % (t_floor, cfg.rules.mission_window_s))
    return None


def evaluate_design(design: Design, cfg: Config, dx: float = 6.0,
                    ref_m2: Optional[float] = None,
                    ref_m3: Optional[float] = None,
                    ref_gm: Optional[float] = None) -> Dict[str, Any]:
    try:
        ac = Aircraft(design, cfg)
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
                    "m2_turn_limit": card.r2.turn_limited_by})
    if card.r3 is not None:
        row.update({"m3_laps": card.r3.laps, "m3_time_s": card.r3.time_s,
                    "m3_energy_margin": card.r3.energy_margin})
    row.update(ac.summary())
    return row


def _worker(task) -> Dict[str, Any]:
    design, cfg, dx = task[0], task[1], task[2]
    ref2 = task[3] if len(task) > 3 else None
    ref3 = task[4] if len(task) > 4 else None
    refg = task[5] if len(task) > 5 else None
    return evaluate_design(design, cfg, dx, ref2, ref3, refg)


def _pool_map(func, tasks, processes: int):
    if processes == 1 or len(tasks) < 8:
        return [func(t) for t in tasks]
    with Pool(processes=processes) as pool:
        chunk = max(1, len(tasks) // (processes * 8))
        return pool.map(func, tasks, chunksize=chunk)


def _n_processes(settings: SearchSettings) -> int:
    if settings.processes > 0:
        return settings.processes
    try:
        return max(1, cpu_count() - 1)
    except NotImplementedError:
        return 1


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


def _screen_one(task) -> Optional[Dict[str, Any]]:
    mk, pk, nm, pack_key, n_sys, cfg = task
    if nm % n_sys != 0:          # motors must divide evenly across systems
        return None
    motor, prop = MOTORS[mk], PROPS[pk]
    pack = PACKS[pack_key]
    ps = PropulsionSystem(motor, prop, nm, pack, cfg.env,
                          n_systems=n_sys)
    rc = ps.rating_check()
    if (rc["motor_current_margin"] < 1.0 or rc["motor_power_margin"] < 1.0
            or rc["esc_current_margin"] < 1.0 or rc["pack_current_margin"] < 1.0
            or rc["tip_mach"] > 0.82):
        return None
    # RULE 3.2.3.2d - the 100 A blade fuse and the total budget are hard stops,
    # and they bind at static full throttle, which is exactly this point.
    if (rc["total_current_a"] > cfg.rules.max_total_current_a
            or rc["fuse_margin"] < 1.0):
        return None
    if rc["static_thrust_n"] < 4.0:
        return None
    op_cruise = ps.solve(26.0, 0.70, 1.0)
    op_fast = ps.solve(34.0, 0.70, 1.0)
    if op_cruise.thrust_total_n <= 0.5:
        return None
    # A DBF airframe needs static thrust for the takeoff box AND thrust that
    # survives to lap speed.  Reward both, penalise dead mass.
    fom = ((rc["static_thrust_n"] ** 0.35)
           * ((op_cruise.thrust_total_n * 26.0 + op_fast.thrust_total_n * 34.0)
              ** 0.65)
           / (ps.mass_kg ** 0.30))
    return {"motor": mk, "prop": pk, "n_motors": nm, "pack": pack_key,
            "n_systems": n_sys, "fom": fom,
            "static_thrust_n": rc["static_thrust_n"],
            "mass_kg": ps.mass_kg, "wh": ps.total_watt_hours}


def screen_propulsion(space: DesignSpace, cfg: Config,
                      settings: SearchSettings,
                      verbose: bool = True) -> List[Dict[str, Any]]:
    packs = candidate_packs(space, cfg, settings.min_pack_wh)
    tasks = [(mk, pk, nm, pack_key, n_sys, cfg)
             for mk in space.motors
             for pk in space.props
             for nm in space.motor_counts
             for (pack_key, n_sys) in packs
             if nm % n_sys == 0]
    if verbose:
        print("  screening %d propulsion combinations "
              "(%d pack/system options)..." % (len(tasks), len(packs)))
    raw = _pool_map(_screen_one, tasks, _n_processes(settings))
    survivors = [r for r in raw if r]
    survivors.sort(key=lambda r: -r["fom"])

    # Keep a diverse shortlist: the best few per motor, then fill globally.
    per_motor: Dict[str, int] = {}
    keep: List[Dict[str, Any]] = []
    for r in survivors:
        c = per_motor.get(r["motor"], 0)
        if c < settings.per_motor_keep:
            keep.append(r)
            per_motor[r["motor"]] = c + 1
        if len(keep) >= settings.keep_propulsion:
            break
    for r in survivors:
        if len(keep) >= settings.keep_propulsion:
            break
        if r not in keep:
            keep.append(r)
    if verbose:
        print("  %d passed ratings, shortlisted %d" % (len(survivors), len(keep)))
    return keep


# ==========================================================================
# Stage 2/3 helpers
# ==========================================================================
def _make_design(base: Dict[str, Any], airfoil: str, x: Sequence[float],
                 n_sim: int, taper: float) -> Design:
    """x = (span, AR, declared sensor mass, M3 sensor fraction, taper,
    drop height)."""
    span = max(x[0], 0.3)
    ar = min(max(x[1], 3.2), 15.5)
    taper = min(max(x[4], 0.30), 1.0) if len(x) > 4 else taper
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
                  taper=taper)


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
    m_energy = usable_j * 9.0 * 0.60 / (9.80665 * 24.0 * 280.0)
    m_empty = 1.4 + 3.0 * (span * chord) + base.get("mass_kg", 1.0)
    # RULE 3.2.1a - 55 lb gross is a hard ceiling, so it is a bound on how
    # much cargo is worth sampling at all.
    capacity = (min(m_thrust, m_energy, cfg.rules.max_gross_mass_kg)
                - m_empty - sensor_kg)
    if capacity <= 0.0:
        return 0
    n = int(round(rng.random() * capacity / max(sim_kg, 0.05)))
    return max(0, min(n, cfg.rules.max_simulators))


def coarse_sweep(shortlist: List[Dict[str, Any]], space: DesignSpace,
                 cfg: Config, settings: SearchSettings,
                 verbose: bool = True,
                 refs: Tuple[Optional[float], Optional[float], Optional[float]] = (None, None, None)
                 ) -> List[Dict[str, Any]]:
    rng = random.Random(settings.seed)
    bounds = space.continuous_bounds(cfg)
    max_sim = cfg.rules.max_simulators
    tasks: List[Tuple[Design, Config, float]] = []
    for base in shortlist:
        pts = lhs(settings.lhs_points, bounds, rng)
        for af in space.airfoils:
            for x in pts:
                span, ar = x[0], x[1]
                n_sim = _plausible_sim_count(base, span, span / ar, x[2],
                                             x[3], rng, cfg)
                tasks.append((_make_design(base, af, x, n_sim, space.taper),
                              cfg, settings.dx_coarse, refs[0], refs[1],
                              refs[2]))
    if verbose:
        print("  coarse sweep: %d designs..." % len(tasks))
    rows = _pool_map(_worker, tasks, _n_processes(settings))
    rows = [r for r in rows if r.get("feasible")]
    rows.sort(key=lambda r: -r["objective"])
    if verbose:
        print("  %d feasible" % len(rows))
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
            "pack": d.pack, "n_systems": d.n_systems}

    def negobj(x: List[float], n_sim: int) -> float:
        des = _make_design(base, d.airfoil, x, n_sim, d.taper)
        return -evaluate_design(des, cfg, settings.dx_refine,
                                ref2, ref3, refg)["objective"]

    frac = (d.sensor_m3_kg / d.sensor_mass_kg
            if d.sensor_mass_kg > 0 and d.sensor_m3_kg > 0 else 1.0)
    x0 = [d.span_m, d.span_m / d.chord_m, d.sensor_mass_kg, frac,
          d.taper, d.drop_height_in]
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
           refs: Tuple[Optional[float], Optional[float], Optional[float]] = (None, None, None)
           ) -> List[Dict[str, Any]]:
    seeds = _dedupe_configs(rows, settings.top_for_refine)
    if verbose:
        print("  refining %d configurations..." % len(seeds))
    tasks = [(r, cfg, space, settings, refs[0], refs[1], refs[2])
             for r in seeds]
    out = _pool_map(_refine_one, tasks, _n_processes(settings))
    out = [r for r in out if r.get("feasible")]
    out.sort(key=lambda r: -r["objective"])
    return out


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


def optimise(cfg: Config, space: Optional[DesignSpace] = None,
             settings: Optional[SearchSettings] = None,
             verbose: bool = True) -> OptimisationResult:
    space = space or DesignSpace()
    settings = settings or SearchSettings()
    n_eval = 0

    if verbose:
        print("Stage 1: propulsion screen")
    shortlist = screen_propulsion(space, cfg, settings, verbose)
    if not shortlist:
        raise RuntimeError("no propulsion combination passed its ratings - "
                           "loosen the catalogue or the Wh limit")

    if verbose:
        print("Stage 2: coarse design sweep")
    refs = (cfg.rules.field_best_m2_lb_per_s,
            cfg.rules.field_best_m3_lap_lb,
            cfg.rules.field_best_gm_lb_in)
    coarse = coarse_sweep(shortlist, space, cfg, settings, verbose, refs)
    n_eval += len(shortlist) * settings.lhs_points * len(space.airfoils)
    if not coarse:
        raise RuntimeError("no feasible design in the coarse sweep - check the "
                           "takeoff distance, Wh limit and payload bounds")

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
    refined = refine(coarse, space, cfg, settings, verbose, refs)
    n_eval += len(refined) * (settings.nm_iters + cfg.rules.max_simulators)

    pool_rows = refined + coarse[: settings.top_for_refine]
    pool_rows.sort(key=lambda r: -r["objective"])
    finalists = _dedupe_configs(pool_rows, settings.top_for_final)

    if verbose:
        print("Stage 4: final high-resolution scoring of %d designs"
              % len(finalists))
    tasks = [(r["design"], cfg, settings.dx_final, refs[0], refs[1], refs[2])
             for r in finalists]
    final_rows = _pool_map(_worker, tasks, _n_processes(settings))
    final_rows = [r for r in final_rows if r.get("feasible")]
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

    pts = [(r["m2_raw_lb_per_s"], r["m3_raw_lap_lb"]) for r in final_rows]
    front_idx = pareto_front(pts) if pts else []
    pareto = [final_rows[i] for i in front_idx]
    pareto.sort(key=lambda r: -r["m2_raw_lb_per_s"])

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

    `ref_m2`/`ref_m3` must be the SAME references the run reported its score
    against.  Without them a self-normalised run falls back to the default
    field estimates, and because self mode does not cap the ratio at 1 the
    table prints scores above the 6.0 ceiling - numbers that cannot be
    compared with anything else in the report.
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
