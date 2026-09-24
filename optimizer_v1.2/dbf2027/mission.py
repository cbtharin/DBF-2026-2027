"""Mission simulation: takeoff roll, climb, laps, energy budget.

Force balance used everywhere:
    level flight    L = W,          T = D + m dV/dt
    ground roll     m dV/dt = T - D - mu (W - L)
    steady turn     L = n W,  n = 1/cos(phi),  turn rate = g sqrt(n^2 - 1) / V
    climb           ROC = V (T - D) / W

Numerical scheme
----------------
Each straight and the ground roll are integrated as an ODE in time on the
state y = [airspeed, ground distance, energy] with classical RK4, so lap time,
distance and energy are all 4th-order accurate and consistent with each other
(energy is a state, not a running sum tacked on afterwards).

Two events are found by root-finding rather than by stepping past them:
  * the point on a straight where the pilot must cut power to arrive at the
    turn at corner speed;
  * the end of the leg itself.
Both are located by bisection on the length of the final RK4 step, so the
integrator lands exactly on the event instead of overshooting by up to one
step.  That is what removes the step-size sensitivity.

Corner speed and best rate of climb are found by a coarse scan followed by
golden-section refinement inside the winning bracket - a bare scan returns
the best *sample*, which is a discretisation error of order (range / n_scan).

Wind is a steady component along the runway: it changes ground speed on the
straights and therefore lap time, and it does not change turn rate.

Landing is simulated - descent, approach, flare and braked rollout - because
it consumes battery and because you must land successfully to score at all.
Its duration is reported as `landing_time_s` and rolled into
`total_flight_time_s`, but it is NOT added to the scored mission time: the
rules stop the clock as the aircraft crosses the line at the end of the last
lap.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Sequence, Tuple

from .aerodynamics import AeroModel
from .aircraft import Aircraft
from .config import Config
from .numerics import (brent, gauss_legendre_5, rk4_step, scan_then_refine)
from .propulsion import PropulsionSystem

# Re-simulate a lap once the pack has drained this much since the last
# simulated lap, instead of reusing the converged lap.  Bounds the error
# introduced by freezing the lap while the battery sags.
SOC_RESIM_THRESHOLD = 0.08

_MIN_GROUND_SPEED = 0.5      # m/s, keeps dt finite in a near-gale headwind


# ==========================================================================
# Results
# ==========================================================================
@dataclass
class TakeoffResult:
    ok: bool
    ground_roll_m: float
    time_s: float
    energy_j: float
    v_liftoff: float
    reason: str = ""


@dataclass
class TurnCapability:
    v_turn: float
    turn_rate_rad_s: float
    load_factor: float
    bank_deg: float
    radius_m: float
    limited_by: str


@dataclass
class MissionResult:
    ok: bool
    laps: int
    time_s: float                 # throttle-up to crossing the line on the last lap
    energy_j: float
    energy_margin: float          # fraction of the pack still unused
    lap_times: List[float] = field(default_factory=list)
    v_max: float = 0.0
    v_stall: float = 0.0
    v_turn: float = 0.0
    turn_limited_by: str = ""
    takeoff_m: float = 0.0
    landing_time_s: float = 0.0     # not in the scored window, see below
    landing_roll_m: float = 0.0
    total_flight_time_s: float = 0.0   # throttle-up to wheels stopped
    reason: str = ""


# ==========================================================================
# RK4 driver with event location
# ==========================================================================
def _rk4_until(deriv: Callable[[float, Sequence[float]], List[float]],
               t0: float, y0: Sequence[float], h: float,
               phi: Callable[[Sequence[float]], float],
               max_steps: int = 4000,
               refine_iters: int = 18
               ) -> Tuple[float, List[float], bool]:
    """Integrate with RK4 steps of size h until phi(y) crosses zero.

    The crossing step is shortened by bisection so the returned state sits on
    the event to ~h/2^refine_iters, instead of up to a full step past it.
    Returns (t, y, event_reached).
    """
    t = t0
    y = list(y0)
    if phi(y) <= 0.0:
        return t, y, True
    for _ in range(max_steps):
        y_new = rk4_step(deriv, t, y, h)
        if phi(y_new) <= 0.0:
            lo, hi = 0.0, h
            for _ in range(refine_iters):
                mid = 0.5 * (lo + hi)
                if phi(rk4_step(deriv, t, y, mid)) > 0.0:
                    lo = mid
                else:
                    hi = mid
            return t + hi, rk4_step(deriv, t, y, hi), True
        t += h
        y = y_new
    return t, y, False


# ==========================================================================
# Steady-state performance points
# ==========================================================================
def max_level_speed(aero: AeroModel, ps: PropulsionSystem, weight_n: float,
                    soc: float) -> float:
    v_lo = aero.v_stall(weight_n) * 1.02
    v_hi = 60.0

    def excess(v: float) -> float:
        return ps.thrust(v, soc, 1.0) - aero.drag(v, weight_n)

    if excess(v_lo) <= 0.0:
        return v_lo
    if excess(v_hi) > 0.0:
        return v_hi
    return brent(excess, v_lo, v_hi, xtol=1e-9)


def turn_capability(aero: AeroModel, ps: PropulsionSystem, weight_n: float,
                    soc: float, cfg: Config,
                    n_struct_limit: Optional[float] = None) -> TurnCapability:
    """Corner speed: the airspeed that maximises sustained turn rate.

    At each speed the load factor is the smallest of the CLmax limit, the
    thrust limit (found by Brent, not by scanning) and the structural/bank
    limit.  The outer maximisation over speed is a scan plus golden section.
    """
    g = cfg.env.g
    # Normal flight must stay inside LIMIT load, not ultimate.  That limit is
    # derived per design from the V-n diagram rather than assumed.
    n_struct = (n_struct_limit if n_struct_limit is not None
                else cfg.build.max_turn_load_factor)
    n_bank = 1.0 / math.cos(math.radians(cfg.build.max_bank_deg))
    v_s1 = aero.v_stall(weight_n)
    v_max = max_level_speed(aero, ps, weight_n, soc)
    if v_max <= v_s1 * 1.02:
        return TurnCapability(v_s1, 1e-3, 1.0, 0.0, 1e6,
                              "cannot sustain level flight")

    def load_factor_at(v: float) -> Tuple[float, str]:
        n_cl = aero.max_lift(v) / weight_n
        if n_cl <= 1.0:
            return 1.0, "CLmax"
        t_avail = ps.thrust(v, soc, 1.0)
        n_hi = min(n_cl, n_struct, n_bank)

        def excess(n: float) -> float:
            return t_avail - aero.drag(v, n * weight_n)

        if excess(1.0) <= 0.0:
            return 1.0, "thrust"
        if excess(n_hi) >= 0.0:
            lim = ("CLmax" if n_hi == n_cl else
                   ("structure" if n_hi == n_struct else "bank limit"))
            return n_hi, lim
        return brent(excess, 1.0, n_hi, xtol=1e-10), "thrust"

    def rate(v: float) -> float:
        n, _ = load_factor_at(v)
        if n <= 1.0 + 1e-9:
            return 0.0
        return g * math.sqrt(n * n - 1.0) / v

    v_best, rate_best = scan_then_refine(rate, v_s1 * 1.03, v_max,
                                         n_scan=14, tol=1e-4)
    if rate_best <= 1e-6:
        return TurnCapability(v_s1, 1e-3, 1.0, 0.0, 1e6, "no sustained turn")
    n, lim = load_factor_at(v_best)
    return TurnCapability(v_best, rate_best, n,
                          math.degrees(math.acos(min(1.0, 1.0 / n))),
                          v_best / rate_best, lim)


def _decel_distance(aero: AeroModel, mass_kg: float, weight_n: float,
                    v_from: float, v_to: float) -> float:
    """Power-off distance to slow from v_from to v_to.

    From m v dv = -D dx,  x = integral_{v_to}^{v_from} m v / D(v) dv,
    evaluated with 5-point Gauss-Legendre (exact to degree 9) rather than a
    midpoint sum.
    """
    if v_to >= v_from:
        return 0.0
    return gauss_legendre_5(
        lambda v: mass_kg * v / max(aero.drag(v, weight_n), 1e-6),
        v_to, v_from)


# ==========================================================================
# Straight legs
# ==========================================================================
def fly_straight(aero: AeroModel, ps: PropulsionSystem, mass_kg: float,
                 weight_n: float, soc: float, ground_len_m: float,
                 v_in: float, v_out_target: float, wind_along: float,
                 dx: float = 6.0,
                 v_cruise: float = float("inf")) -> Tuple[float, float, float]:
    """Integrate one straight leg.  Returns (time_s, energy_J, v_exit).

    `wind_along` is the wind component opposing motion (+ for a headwind).
    `dx` is a nominal ground-track step used to size the RK4 time step; with
    RK4 plus event location the answer is nearly independent of it.

    Flown as: throttle up towards `v_cruise`, hold it on part throttle,
    then power off at the last moment that still allows arrival at the
    turn at corner speed.

    `v_cruise` is what lets a design trade speed for energy.  With it at
    infinity this is the old behaviour - flat out everywhere - and a
    design whose battery could not sustain that simply failed, even when
    it would have completed the mission comfortably a few m/s slower.
    That is not how anyone flies, and it was costing real payload: at the
    optimum the aeroplane was burning 84 percent more energy than a
    best-L/D cruise over the same ground.
    """
    v_floor = 0.90 * aero.v_stall(weight_n)
    # Once the power is cut the propellers keep turning, and a freewheeling
    # disc is a real drag item - it shortens the glide and so costs lap
    # time.  Leaving it out flattered the deceleration.
    disc_area = (math.pi * (ps.prop.diameter_m * 0.5) ** 2) * ps.n_motors
    windmill_area = disc_area * aero.build.windmill_drag_coeff

    def make_deriv(powered: bool):
        def deriv(_t: float, y: Sequence[float]) -> List[float]:
            v = max(y[0], v_floor)
            drag = aero.drag(v, weight_n)
            if powered:
                if v >= v_cruise - 1e-9:
                    # At the commanded speed: hold it on whatever throttle
                    # that takes, rather than accelerating past it.
                    thr = ps.throttle_for_thrust(v, soc, drag)
                    thrust = ps.thrust(v, soc, thr)
                    power = ps.power_elec(v, soc, thr)
                else:
                    thrust = ps.thrust(v, soc, 1.0)
                    power = ps.power_elec(v, soc, 1.0)
            else:
                thrust, power = 0.0, 0.0
                drag += 0.5 * aero.env.rho * v * v * windmill_area
            return [(thrust - drag) / mass_kg,
                    max(v - wind_along, _MIN_GROUND_SPEED),
                    power]
        return deriv

    y = [max(v_in, 1.0), 0.0, 0.0]
    t = 0.0
    h = min(0.5, max(0.01, dx / max(y[0] - wind_along, 5.0)))

    need_decel = v_out_target < 1e8
    if need_decel:
        def phi_switch(s: Sequence[float]) -> float:
            v = max(s[0], v_floor)
            return ((ground_len_m - s[1])
                    - _decel_distance(aero, mass_kg, weight_n, v, v_out_target))

        t, y, _ = _rk4_until(make_deriv(True), t, y, h, phi_switch)

    # remaining distance, powered if no turn follows, otherwise gliding down
    powered_tail = not need_decel
    deriv = make_deriv(powered_tail)
    t, y, _ = _rk4_until(deriv, t, y, h,
                         lambda s: ground_len_m - s[1])
    return t, y[2], max(y[0], v_floor)


# ==========================================================================
# Takeoff and climb
# ==========================================================================
def takeoff(aero: AeroModel, ps: PropulsionSystem, mass_kg: float,
            weight_n: float, cfg: Config, soc: float = 1.0,
            gear_height_m: Optional[float] = None) -> TakeoffResult:
    """Ground roll to rotation, integrated with RK4 and stopped exactly at
    the liftoff speed."""
    b, env, r = cfg.build, cfg.env, cfg.rules
    v_s = aero.v_stall(weight_n)
    v_lof = b.rotate_margin * v_s
    gh = b.gear_height_m if gear_height_m is None else gear_height_m
    h_wheel = gh + 0.5 * aero.body.fuse_height_m
    hb = 16.0 * h_wheel / aero.wing.span_m
    ground_effect = (hb * hb) / (1.0 + hb * hb)     # induced drag multiplier

    cl_g = b.takeoff_cl_ground
    af = aero.wing.airfoil
    cd_extra = (af.k_drag * (cl_g - af.cl_at_cd_min) ** 2
                + ground_effect * aero.k_induced * cl_g * cl_g)

    def accel(v: float) -> float:
        q = 0.5 * env.rho * v * v
        lift = q * aero.S * cl_g
        drag = q * aero.S * (aero.cd0(max(v, 3.0)) + cd_extra)
        fric = b.rolling_friction * max(weight_n - lift, 0.0)
        return (ps.thrust(v, soc, 1.0) - drag - fric) / mass_kg

    # Static check first: if it will not roll, say so rather than integrating.
    if accel(max(env.wind_mps, 0.5)) <= 0.02:
        return TakeoffResult(False, 0.0, 0.0, 0.0, env.wind_mps,
                             "cannot accelerate from rest")

    def deriv(_t: float, y: Sequence[float]) -> List[float]:
        v = max(y[0], 0.1)
        return [accel(v),
                max(v - env.wind_mps, 0.0),
                ps.power_elec(v, soc, 1.0)]

    y0 = [max(env.wind_mps, 0.5), 0.0, 0.0]   # airspeed at brake release = headwind
    t, y, reached = _rk4_until(deriv, 0.0, y0, 0.05,
                               lambda s: v_lof - s[0], max_steps=3000)
    x, e = y[1], y[2]
    if not reached:
        return TakeoffResult(False, x, t, e, y[0],
                             "never reached rotation speed")
    ok = x <= r.takeoff_distance_m
    return TakeoffResult(ok, x, t, e, y[0],
                         "" if ok else "ground roll %.1f m exceeds the %.1f m limit"
                         % (x, r.takeoff_distance_m))


def climb(aero: AeroModel, ps: PropulsionSystem, mass_kg: float,
          weight_n: float, cfg: Config, soc: float,
          v_start: float) -> Tuple[float, float, float, float]:
    """Climb to pattern altitude at best rate of climb.

    Returns (time_s, energy_J, v_exit, ground_distance_m).
    """
    h = cfg.course.pattern_altitude_m
    v_s = aero.v_stall(weight_n)

    def roc(v: float) -> float:
        return v * (ps.thrust(v, soc, 1.0) - aero.drag(v, weight_n)) / weight_n

    best_v, best_roc = scan_then_refine(roc, v_s * 1.10, v_s * 2.8,
                                        n_scan=12, tol=1e-4)
    if best_roc <= 0.2:
        return 60.0, 0.0, v_start, 0.0          # effectively cannot climb
    t = h / best_roc
    e = ps.power_elec(best_v, soc, 1.0) * t
    v_horiz = math.sqrt(max(best_v * best_v - best_roc * best_roc, 1.0))
    return t, e, best_v, v_horiz * t


def landing(aero: AeroModel, ps: PropulsionSystem, mass_kg: float,
            weight_n: float, cfg: Config, soc: float,
            gear_height_m: Optional[float] = None
            ) -> Tuple[float, float, float]:
    """Descend, approach, flare and roll out.

    Returns (time_s, energy_J, ground_roll_m).

    RULE 3.3.3d.  "Landing is not part of the
    5 minute time window" and the clock stops when the aircraft crosses the
    line at the end of the final lap.  So this time is reported and its
    energy is charged against the battery, but it is deliberately NOT added
    to the scored mission time.  A successful landing is still required to
    score at all, so this phase also decides feasibility.
    """
    b, env = cfg.build, cfg.env
    v_s = aero.v_stall(weight_n)
    v_app = b.approach_margin * v_s

    # descent from pattern altitude at a steady rate, partial power
    t_desc = cfg.course.pattern_altitude_m / max(b.descent_rate_mps, 0.5)
    thr_desc = ps.throttle_for_thrust(
        v_app, soc, max(0.35 * aero.drag(v_app, weight_n), 0.0))
    e_desc = ps.power_elec(v_app, soc, thr_desc) * t_desc

    # flare and touchdown, engine at idle
    t_flare = b.flare_time_s
    v_td = 1.10 * v_s

    # ground roll with brakes, integrated with RK4 down to taxi speed
    gh = b.gear_height_m if gear_height_m is None else gear_height_m
    cl_g = b.takeoff_cl_ground
    af = aero.wing.airfoil
    h_wheel = gh + 0.5 * aero.body.fuse_height_m
    hb = 16.0 * h_wheel / aero.wing.span_m
    ground_effect = (hb * hb) / (1.0 + hb * hb)
    cd_extra = (af.k_drag * (cl_g - af.cl_at_cd_min) ** 2
                + ground_effect * aero.k_induced * cl_g * cl_g)

    def decel(v: float) -> float:
        q = 0.5 * env.rho * v * v
        lift = q * aero.S * cl_g
        drag = q * aero.S * (aero.cd0(max(v, 3.0)) + cd_extra)
        fric = b.braking_friction * max(weight_n - lift, 0.0)
        return -(drag + fric) / mass_kg

    def deriv(_t: float, y: Sequence[float]) -> List[float]:
        v = max(y[0], 0.5)
        return [decel(v), max(v - env.wind_mps, 0.0), 0.0]

    y0 = [v_td, 0.0, 0.0]
    t_roll, y, _ = _rk4_until(deriv, 0.0, y0, 0.05,
                              lambda st: st[0] - 2.0, max_steps=4000)
    return t_desc + t_flare + t_roll, e_desc, y[1]


# ==========================================================================
# Lap
# ==========================================================================
# RULE 3.3.1 - the course segments and headings come from cfg.course,
# which lays the start/finish line between markers 500 ft each way.
def simulate_lap(aero: AeroModel, ps: PropulsionSystem, mass_kg: float,
                 weight_n: float, soc: float, cfg: Config,
                 turn: TurnCapability, v_entry: float,
                 dx: float = 6.0,
                 v_cruise: float = float("inf")
                 ) -> Tuple[float, float, float]:
    """One full lap of the course.  Returns (time_s, energy_J, v_exit)."""
    segs = cfg.course.segments
    heads = cfg.course.headings
    wind = cfg.env.wind_mps
    t_tot = 0.0
    e_tot = 0.0
    v = v_entry

    # throttle and power needed to hold the sustained turn
    thrust_turn = aero.drag(turn.v_turn, turn.load_factor * weight_n)
    thr = ps.throttle_for_thrust(turn.v_turn, soc, thrust_turn)
    p_turn = ps.power_elec(turn.v_turn, soc, thr)

    for i, (kind, val) in enumerate(segs):
        if kind == "straight":
            wind_along = wind * heads[i]        # +1 upwind, -1 downwind, 0 cross
            nxt = segs[(i + 1) % len(segs)]
            v_target = turn.v_turn if nxt[0] == "turn" else 1e9
            dt_, de_, v = fly_straight(aero, ps, mass_kg, weight_n, soc,
                                       val, v, v_target, wind_along, dx,
                                       v_cruise)
            t_tot += dt_
            e_tot += de_
        else:
            dt_ = math.radians(val) / max(turn.turn_rate_rad_s, 1e-6)
            t_tot += dt_
            e_tot += p_turn * dt_
            v = turn.v_turn
    return t_tot, e_tot, v


def lap_ground_distance(cfg: Config, turn: TurnCapability) -> float:
    """Ground track of one lap, in metres: straights plus turn arcs.

    Turn *time* is angle / turn rate and never needs a radius, which is why
    the rest of the model does not carry one.  Average ground speed does
    need it, so it is recovered here from the sustained turn:

        r = v_turn / omega,        arc = theta * r

    Used only to convert a climb-out ground credit into a time credit.
    """
    d = cfg.course.straight_length_m
    omega = max(turn.turn_rate_rad_s, 1e-6)
    radius = turn.v_turn / omega
    for kind, val in cfg.course.segments:
        if kind == "turn":
            d += math.radians(val) * radius
    return d


# ==========================================================================
# Missions
# ==========================================================================
def fly_mission(ac: Aircraft, aero: AeroModel, weight_n: float, cfg: Config,
                required_laps: Optional[int] = None,
                max_laps: int = 40, dx: float = 6.0,
                energy_reserve: float = 0.08,
                v_cruise: float = float("inf")) -> MissionResult:
    """Common engine for M1/M2/M3.

    required_laps set  -> fly exactly that many, report the time (Mission 2).
    required_laps None -> fly as many as fit in the window (Mission 3).
    """
    ps = ac.prop_sys
    mass = weight_n / cfg.env.g
    # Every pack on the aeroplane, not one of them.  With two propulsion
    # systems the packs discharge together, so state of charge must be
    # measured against the total; dividing by a single pack made the
    # voltage sag twice as fast as it should, understated thrust for the
    # whole second half of a mission, and reported a NEGATIVE energy
    # margin for a flight that had in fact finished with energy to spare.
    pack_j = max(ac.total_watt_hours * 3600.0, 1e-6)
    usable = ac.usable_energy_j() * (1.0 - energy_reserve)
    window = cfg.rules.mission_window_s

    v_s = aero.v_stall(weight_n)
    to = takeoff(aero, ps, mass, weight_n, cfg, soc=1.0,
                 gear_height_m=ac.gear_height_m)
    if not to.ok:
        return MissionResult(False, 0, 0.0, to.energy_j, 1.0,
                             v_stall=v_s, takeoff_m=to.ground_roll_m,
                             reason=to.reason or "takeoff")

    # RULE 3.2.1c - every joule of the takeoff comes out of the same
    # propulsion pack budget the laps draw on.  There is no external
    # assist and no separate launch energy.
    energy = to.energy_j
    soc = max(0.03, 1.0 - energy / pack_j)
    t_cl, e_cl, v_cl, x_cl = climb(aero, ps, mass, weight_n, cfg, soc,
                                   to.v_liftoff)
    energy += e_cl
    t_now = to.time_s + t_cl

    turn = turn_capability(aero, ps, weight_n, 0.70, cfg,
                           getattr(ac, 'n_limit', None))
    if turn.turn_rate_rad_s <= 1e-3:
        return MissionResult(False, 0, t_now, energy, 0.0, v_stall=v_s,
                             takeoff_m=to.ground_roll_m,
                             reason="cannot sustain a level turn")
    v_max = max_level_speed(aero, ps, weight_n, 0.70)

    lap_times: List[float] = []
    v = v_cl
    laps = 0
    cached: Optional[Tuple[float, float]] = None    # (time, energy)
    cached_soc = 1.0
    credit = x_cl                                   # ground covered during climb

    while True:
        if required_laps is not None and laps >= required_laps:
            break
        if laps >= max_laps:
            break
        soc = max(0.03, 1.0 - energy / pack_j)

        # Reuse the converged lap only while the pack state is close to the
        # state it was simulated at; otherwise re-simulate.
        if cached is not None and abs(cached_soc - soc) <= SOC_RESIM_THRESHOLD:
            lt, le = cached
        else:
            lt, le, v = simulate_lap(aero, ps, mass, weight_n, soc, cfg,
                                     turn, v, dx, v_cruise)
            cached, cached_soc = (lt, le), soc

        # No pilot flies the theoretical lap.  One factor, applied to time
        # and to the energy that time costs.
        eff = max(min(cfg.build.pilot_efficiency, 1.0), 0.5)
        lt, le = lt / eff, le / eff

        if laps == 0 and credit > 0.0:
            # The climb-out is flown over the course, so the ground it
            # covers counts against the first lap.  The time that saves is
            # that distance flown at the lap's average GROUND speed - which
            # needs the true lap ground distance, turn arcs included, not
            # just the straights.
            #
            # It is then capped at the climb's own duration.  Climbing is
            # done at best-rate speed, which is slower than lap speed, so
            # covering ground while climbing can never save more time than
            # the climb took; without the cap a taller pattern altitude
            # came out FASTER, which it must not.
            v_avg = lap_ground_distance(cfg, turn) / max(lt, 1e-3)
            saved = min(credit / max(v_avg, 1.0), t_cl)
            lt = max(lt - saved, 0.0)

        if energy + le > usable:
            if required_laps is not None:
                return MissionResult(False, laps, t_now, energy, 0.0,
                                     lap_times, v_max, v_s, turn.v_turn,
                                     turn.limited_by, to.ground_roll_m,
                                     reason="out of battery before %d laps"
                                     % required_laps)
            break
        if required_laps is None and t_now + lt > window:
            break

        t_now += lt
        energy += le
        laps += 1
        lap_times.append(lt)

    if required_laps is not None and laps < required_laps:
        return MissionResult(False, laps, t_now, energy, 0.0, lap_times,
                             v_max, v_s, turn.v_turn, turn.limited_by,
                             to.ground_roll_m,
                             reason="could not complete required laps")
    if required_laps is not None and t_now > window:
        return MissionResult(False, laps, t_now, energy, 0.0, lap_times,
                             v_max, v_s, turn.v_turn, turn.limited_by,
                             to.ground_roll_m,
                             reason="%.1f s exceeds the %.0f s window"
                             % (t_now, window))
    if laps == 0:
        return MissionResult(False, 0, t_now, energy, 0.0, lap_times,
                             v_max, v_s, turn.v_turn, turn.limited_by,
                             to.ground_roll_m,
                             reason="no laps completed in the window")

    # Landing: required to score, charged to the battery, reported
    # separately, and excluded from the scored window per the rules.
    soc = max(0.03, 1.0 - energy / pack_j)
    t_land, e_land, roll = landing(aero, ps, mass, weight_n, cfg, soc,
                                   ac.gear_height_m)
    if energy + e_land > ac.usable_energy_j():
        return MissionResult(False, laps, t_now, energy + e_land, 0.0,
                             lap_times, v_max, v_s, turn.v_turn,
                             turn.limited_by, to.ground_roll_m,
                             t_land, roll, t_now + t_land,
                             reason="not enough battery left to land safely")
    energy += e_land
    margin = 1.0 - energy / pack_j
    return MissionResult(True, laps, t_now, energy, margin, lap_times,
                         v_max, v_s, turn.v_turn, turn.limited_by,
                         to.ground_roll_m, t_land, roll, t_now + t_land)


# ==========================================================================
# Choosing a cruise speed
# ==========================================================================
# Full throttle is the right way to fly Mission 2 if the battery can
# stand it, because M2 is scored on weight over TIME.  It is the wrong
# way to fly it if the battery cannot: the aeroplane then simply failed,
# where a real pilot would have eased off and finished a few seconds
# slower for a lower but non-zero score.  Mission 3 is different again -
# it is scored on LAPS, so the right speed is the one that gets the most
# laps out of the pack, which is rarely full throttle.
# Three probes, not seven.  Each one is a whole extra mission
# simulation, and near the optimum almost every design IS energy limited
# - that is the binding constraint - so the search fires on most of
# them.  At seven probes a full run went from about five minutes to well
# over twenty for no measured change in the answer; three keeps the
# capability at a cost the search can absorb.
_CRUISE_STEPS = 3


def _fastest_that_finishes(ac, aero, weight_n, cfg, required_laps, dx,
                           v_max) -> MissionResult:
    """Fly as fast as the energy budget allows.

    Tries full throttle first, because when it works it is optimal and
    costs one simulation.  Only if that runs the pack dry does it search
    downwards for the quickest speed that still completes the mission.
    """
    best = fly_mission(ac, aero, weight_n, cfg, required_laps=required_laps,
                       dx=dx)
    if best.ok:
        return best
    v_s = aero.v_stall(weight_n)
    lo, hi = 1.15 * v_s, max(v_max, 1.2 * v_s)
    found = None
    for i in range(_CRUISE_STEPS):
        # Walk down from the top: the first speed that completes is the
        # fastest one that does, which is what M2 wants.
        v = hi - (hi - lo) * (i + 1) / (_CRUISE_STEPS + 1)
        r = fly_mission(ac, aero, weight_n, cfg,
                        required_laps=required_laps, dx=dx, v_cruise=v)
        if r.ok:
            found = r
            break
    return found if found is not None else best


#: Above this much charge left, the flight ended because the WINDOW ran
#: out, not the battery - and then flying slower can only lose laps, so
#: there is nothing to search for and the extra simulations are waste.
_ENERGY_LIMITED_MARGIN = 0.15


def _most_laps(ac, aero, weight_n, cfg, dx, v_max) -> MissionResult:
    """Fly for laps, not for speed.

    Mission 3 scores laps times sensor weight, so the cruise speed worth
    flying is whichever yields the most laps inside the window - slower
    is often better, because the pack lasts longer.

    The search only runs when the full-throttle flight was stopped by the
    BATTERY.  If it was stopped by the clock, full throttle already got
    every lap available and searching would cost seven simulations to
    confirm it.  That distinction is the difference between a 0.1 s
    evaluation and a 2 s one, and at 24,000 designs per search it is the
    difference between minutes and hours.
    """
    best = fly_mission(ac, aero, weight_n, cfg, required_laps=None, dx=dx)
    if best.ok and best.energy_margin > _ENERGY_LIMITED_MARGIN:
        return best
    v_s = aero.v_stall(weight_n)
    lo, hi = 1.15 * v_s, max(v_max, 1.2 * v_s)
    for i in range(_CRUISE_STEPS):
        v = lo + (hi - lo) * (i + 1) / (_CRUISE_STEPS + 1)
        r = fly_mission(ac, aero, weight_n, cfg, required_laps=None, dx=dx,
                        v_cruise=v)
        if not r.ok:
            continue
        if (r.laps > best.laps
                or (r.laps == best.laps and best.ok
                    and r.time_s < best.time_s)
                or not best.ok):
            best = r
    return best


def mission_1(ac: Aircraft, cfg: Config, dx: float = 6.0) -> MissionResult:
    """Flight mission with no payload - pass/fail (3 laps in the window)."""
    return _fastest_that_finishes(ac, ac.aero_clean, ac.w_m1_n, cfg, 3, dx,
                                  ac._max_level_speed(ac.w_m1_n))


def mission_2(ac: Aircraft, cfg: Config, dx: float = 6.0) -> MissionResult:
    """Delivery flight: sensor in shipping container + optional simulators,
    5 laps, scored on weight / time."""
    return _fastest_that_finishes(ac, ac.aero_clean, ac.w_m2_n, cfg,
                                  cfg.rules.m2_required_laps, dx,
                                  ac._max_level_speed(ac.w_m2_n))


def mission_3(ac: Aircraft, cfg: Config, dx: float = 6.0) -> MissionResult:
    """Sensor flight: sensor deployed for the whole pattern, maximise laps."""
    return _most_laps(ac, ac.aero_deployed, ac.w_m3_n, cfg, dx,
                      ac._max_level_speed(ac.w_m3_n))
