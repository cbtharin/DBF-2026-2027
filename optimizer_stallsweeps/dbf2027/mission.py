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
from .propulsion import PropulsionSystem, motor_thermal, winding_step

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
    # Peak motor winding temperature over the flight, C, and the pack's
    # temperature at its end.  A failed flight reports both as far as it
    # got, so neither is ever left at zero.
    motor_temp_c: float = 0.0
    pack_temp_c: float = 0.0
    # What ended the lap loop: "laps" (the required count was flown),
    # "window" (the clock), "battery" (the safe depth of discharge or the
    # ESC low-voltage cutoff) or "cap" (max_laps).  Mission 3's slower-
    # cruise search runs only when this is "battery".
    stopped_by: str = ""
    # The corner-speed turn every lap flies (m119): radius, m; load factor;
    # bank, deg; rate, deg/s.  0 if the flight ended before it was found.
    turn_radius_m: float = 0.0
    turn_load_factor: float = 0.0
    turn_bank_deg: float = 0.0
    turn_rate_deg_s: float = 0.0


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
    # no thrust, so no level flight, above the propellers' speed ceiling
    v_hi = max(ps.speed_ceiling(soc), v_lo * 1.001)

    # full throttle, so the drag is the power-on drag at that thrust
    def excess(v: float) -> float:
        t = ps.thrust(v, soc, 1.0)
        return t - aero.drag(v, weight_n, t)

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
    # A sustained turn is flown at full throttle, so the slipstream's lift
    # is there: the wing reaches more lift before it stalls, and makes a
    # given lift at a lower CL (AeroModel.slipstream_gains).
    v_s1 = aero.v_stall_powered(weight_n, lambda v: ps.thrust(v, soc, 1.0))
    v_max = max_level_speed(aero, ps, weight_n, soc)
    if v_max <= v_s1 * 1.02:
        return TurnCapability(v_s1, 1e-3, 1.0, 0.0, 1e6,
                              "cannot sustain level flight")

    def load_factor_at(v: float) -> Tuple[float, str]:
        t_avail = ps.thrust(v, soc, 1.0)
        n_cl = aero.max_lift(v, t_avail) / weight_n
        if n_cl <= 1.0:
            return 1.0, "CLmax"
        n_hi = min(n_cl, n_struct, n_bank)

        def excess(n: float) -> float:
            return t_avail - aero.drag(v, n * weight_n, t_avail)

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
        # (n - 1)(n + 1), not n^2 - 1: no cancellation in a gentle turn
        return g * math.sqrt((n - 1.0) * (n + 1.0)) / v

    v_best, rate_best = scan_then_refine(rate, v_s1 * 1.03, v_max,
                                         n_scan=14, tol=1e-4)
    if rate_best <= 1e-6:
        return TurnCapability(v_s1, 1e-3, 1.0, 0.0, 1e6, "no sustained turn")
    n, lim = load_factor_at(v_best)
    # bank = atan(sqrt(n^2 - 1)): acos(1/n) has an infinite slope as n
    # nears 1 and amplifies every rounding error in n
    return TurnCapability(v_best, rate_best, n,
                          math.degrees(math.atan(math.sqrt(
                              max((n - 1.0) * (n + 1.0), 0.0)))),
                          v_best / rate_best, lim)


def _decel_distance(aero: AeroModel, mass_kg: float, weight_n: float,
                    v_from: float, v_to: float,
                    windmill_area_m2: float = 0.0,
                    wind_along: float = 0.0) -> float:
    """Power-off GROUND distance to slow from airspeed v_from to v_to.

    m dv/dt = -D and dx/dt = v - w (w the headwind component), so
        x = integral_{v_to}^{v_from} m (v - w) / D(v) dv,
    D the airframe's drag plus the freewheeling propellers',
    0.5 rho v^2 A_disc Cd_wm (windmill_area_m2 = A_disc Cd_wm) - the drag
    and ground speed the glide is then flown with - evaluated with 5-point
    Gauss-Legendre (exact to degree 9).
    """
    if v_to >= v_from:
        return 0.0
    rho = aero.env.rho
    return gauss_legendre_5(
        lambda v: mass_kg * max(v - wind_along, _MIN_GROUND_SPEED)
        / max(aero.drag(v, weight_n) + 0.5 * rho * v * v * windmill_area_m2,
              1e-6),
        v_to, v_from)


# ==========================================================================
# Straight legs
# ==========================================================================
def fly_straight(aero: AeroModel, ps: PropulsionSystem, mass_kg: float,
                 weight_n: float, soc: float, ground_len_m: float,
                 v_in: float, v_out_target: float, wind_along: float,
                 dx: float = 6.0,
                 v_cruise: float = float("inf"),
                 loss_out: Optional[List[float]] = None,
                 pack_out: Optional[List[float]] = None,
                 full_power: bool = False
                 ) -> Tuple[float, float, float]:
    """Integrate one straight leg.  Returns (time_s, energy_J, v_exit).
    If `loss_out` is given, the motor heat of the leg (J) is appended; if
    `pack_out` is, the Joule heat in ONE pack (J), I_pack^2 R_pack.

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

    `full_power` (BuildStandard.full_throttle_flights) charges the whole
    leg at full throttle: no cruise speed is held, and the glide into the
    turn is flown as a glide - the time and exit speed do not change - but
    costs the energy and heat full throttle would at that airspeed.
    """
    v_floor = 0.90 * aero.v_stall(weight_n)
    # Once the power is cut the propellers keep turning, and a freewheeling
    # disc is a real drag item - it shortens the glide and so costs lap
    # time.  Leaving it out flattered the deceleration.
    # The coefficient comes from the propeller's measured table
    # (Propeller.windmill_drag_coeff); the build standard's is the
    # fallback for a propeller without one.
    disc_area = (math.pi * (ps.prop.diameter_m * 0.5) ** 2) * ps.n_motors
    cd_wm = ps.prop.windmill_drag_coeff
    if cd_wm is None:
        cd_wm = aero.build.windmill_drag_coeff
    windmill_area = disc_area * cd_wm

    def make_deriv(powered: bool):
        def deriv(_t: float, y: Sequence[float]) -> List[float]:
            v = max(y[0], v_floor)
            if powered:
                if v >= v_cruise - 1e-9 and not full_power:
                    # At the commanded speed: hold it on whatever throttle
                    # that takes, rather than accelerating past it - the
                    # thrust that equals the power-on drag it causes (two
                    # fixed-point passes, as the turn does).
                    t_hold = aero.drag(v, weight_n)
                    for _ in range(2):
                        t_hold = aero.drag(v, weight_n, t_hold)
                    thr = ps.throttle_for_thrust(v, soc, t_hold)
                else:
                    thr = 1.0
                thrust = ps.thrust(v, soc, thr)
                # power on: the slipstream's lift and swirl, the blown
                # wing's section drag, the roll trim (AeroModel.drag)
                drag = aero.drag(v, weight_n, thrust)
                power = ps.power_elec(v, soc, thr)
                heat = ps.motor_loss(v, soc, thr)
                pack_heat = ps.pack_heat(v, soc, thr)
            else:
                thrust = 0.0
                drag = (aero.drag(v, weight_n)
                        + 0.5 * aero.env.rho * v * v * windmill_area)
                if full_power:
                    # flown as a glide, charged as full throttle
                    power = ps.power_elec(v, soc, 1.0)
                    heat = ps.motor_loss(v, soc, 1.0)
                    pack_heat = ps.pack_heat(v, soc, 1.0)
                else:
                    power, heat, pack_heat = 0.0, 0.0, 0.0
            return [(thrust - drag) / mass_kg,
                    max(v - wind_along, _MIN_GROUND_SPEED),
                    power, heat, pack_heat]
        return deriv

    y = [max(v_in, 1.0), 0.0, 0.0, 0.0, 0.0]
    t = 0.0
    h = min(0.5, max(0.01, dx / max(y[0] - wind_along, 5.0)))

    need_decel = v_out_target < 1e8
    if need_decel:
        def phi_switch(s: Sequence[float]) -> float:
            v = max(s[0], v_floor)
            return ((ground_len_m - s[1])
                    - _decel_distance(aero, mass_kg, weight_n, v, v_out_target,
                                      windmill_area, wind_along))

        t, y, _ = _rk4_until(make_deriv(True), t, y, h, phi_switch)

    # remaining distance, powered if no turn follows, otherwise gliding down
    powered_tail = not need_decel
    deriv = make_deriv(powered_tail)
    t, y, _ = _rk4_until(deriv, t, y, h,
                         lambda s: ground_len_m - s[1])
    if loss_out is not None:
        loss_out.append(y[3])
    if pack_out is not None:
        pack_out.append(y[4])
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
    # rotation is flown at full power, so on the power-on stall speed
    v_s = aero.v_stall_powered(weight_n, lambda v: ps.thrust(v, soc, 1.0))
    # ...unless the wing cannot make that lift at the attitude the
    # aeroplane can take on its wheels (a low-aspect-ratio wing needs a
    # large angle of attack): then it lifts off where its runway attitude
    # holds its weight (m104).  The same speed when the attitude does not
    # limit it.
    v_att = aero.v_stall_powered(weight_n, lambda v: ps.thrust(v, soc, 1.0),
                                 ground=True)
    v_lof = max(b.rotate_margin * v_s, v_att)
    gh = b.gear_height_m if gear_height_m is None else gear_height_m
    h_wheel = gh + 0.5 * aero.body.fuse_height_m
    hb = 16.0 * h_wheel / aero.wing.span_m
    ground_effect = (hb * hb) / (1.0 + hb * hb)     # induced drag multiplier

    cl_g = b.takeoff_cl_ground
    af = aero.wing.airfoil
    cd_extra = (af.cd_excess(cl_g)
                + ground_effect * aero.k_induced * cl_g * cl_g)

    def accel(v: float) -> float:
        q = 0.5 * env.rho * v * v
        t = ps.thrust(v, soc, 1.0)
        # the slipstream lifts the blown wing in the roll too
        lift = q * aero.S * cl_g * aero.slipstream_gains(v, t)[0]
        drag = q * aero.S * (aero.cd0(max(v, 3.0)) + cd_extra)
        fric = b.rolling_friction * max(weight_n - lift, 0.0)
        return (t - drag - fric) / mass_kg

    # Static check first: if it will not roll, say so rather than integrating.
    if accel(max(env.wind_mps, 0.5)) <= 0.02:
        return TakeoffResult(False, 0.0, 0.0, 0.0, env.wind_mps,
                             "cannot accelerate from rest")
    # dv/dt = a(v) alone, so the speed rises towards the first zero of a(v)
    # and never passes it: rotation is reached only if a(v) > 0 all the way
    # there.  Checked first - integrating a roll that stalls out ran 3000
    # steps before giving up.
    v0 = max(env.wind_mps, 0.5)
    for i in range(1, 41):
        v = v0 + (v_lof - v0) * i / 41.0
        if accel(v) <= 0.0:
            return TakeoffResult(False, 0.0, 0.0, 0.0, v,
                                 "never reached rotation speed: thrust "
                                 "falls below drag and rolling friction "
                                 "at %.1f m/s" % v)

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
          v_start: float) -> Optional[Tuple[float, float, float, float]]:
    """From liftoff to pattern altitude at full throttle, by energy height.

    The aeroplane must gain specific energy
        dE_s = h + (V_c^2 - V_0^2) / (2 g)
    - the height plus the speed change from liftoff (V_0 = v_start) to the
    speed of best climb V_c - and at full throttle gains it at its specific
    excess power P_s(V) = V (T - D) / W.

    1. The speed change, flown level first:
           t_1 = integral_{V_0}^{V_c} m dV / (T - D),
           x_1 = integral m V dV / (T - D),
           E_1 = integral P_elec m dV / (T - D)        (Gauss-Legendre)
       If V_c < V_0 the surplus kinetic energy climbs instead:
           h_left = h - (V_0^2 - V_c^2) / (2 g).
    2. The climb at V_c:  t_2 = h_left / P_s(V_c),  E_2 = P_elec(V_c) t_2,
       ground distance sqrt(V_c^2 - P_s(V_c)^2) t_2.

    Energy is always charged.  Returns (time_s, energy_J, v_exit,
    ground_distance_m), or None if the aeroplane has no excess power at its
    best-climb speed or anywhere between liftoff and it - it cannot climb.
    """
    g = cfg.env.g
    h = cfg.course.pattern_altitude_m
    v_s = aero.v_stall_powered(weight_n, lambda v: ps.thrust(v, soc, 1.0))

    def excess(v: float) -> float:
        t = ps.thrust(v, soc, 1.0)
        return t - aero.drag(v, weight_n, t)

    def p_s(v: float) -> float:
        return v * excess(v) / weight_n

    best_v, best_ps = scan_then_refine(p_s, v_s * 1.10, v_s * 2.8,
                                       n_scan=12, tol=1e-4)
    if best_ps <= 0.0:
        return None
    v0 = max(v_start, 0.5)
    t1 = e1 = x1 = 0.0
    h_left = h
    if best_v > v0:
        # the level speed-up needs excess thrust all the way
        for i in range(7):
            if excess(v0 + (best_v - v0) * i / 6.0) <= 0.0:
                return None
        t1 = gauss_legendre_5(lambda v: mass_kg / excess(v), v0, best_v)
        x1 = gauss_legendre_5(lambda v: mass_kg * v / excess(v), v0, best_v)
        e1 = gauss_legendre_5(
            lambda v: ps.power_elec(v, soc, 1.0) * mass_kg / excess(v),
            v0, best_v)
    else:
        h_left = max(h - (v0 * v0 - best_v * best_v) / (2.0 * g), 0.0)
    t2 = h_left / best_ps
    e2 = ps.power_elec(best_v, soc, 1.0) * t2
    v_horiz = math.sqrt(max(best_v * best_v - best_ps * best_ps, 1.0))
    return t1 + t2, e1 + e2, best_v, x1 + v_horiz * t2


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
    # the slowest it can touch down: its stall, or where its runway
    # attitude holds its weight if that is faster (m104)
    v_s = aero.v_stall(weight_n, ground=True)
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
    cd_extra = (af.cd_excess(cl_g)
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
                 v_cruise: float = float("inf"),
                 loss_out: Optional[List[float]] = None,
                 pack_out: Optional[List[float]] = None
                 ) -> Tuple[float, float, float]:
    """One full lap of the course.  Returns (time_s, energy_J, v_exit).
    If `loss_out` is given, the lap's motor heat (J) is appended; if
    `pack_out` is, the Joule heat in one pack (J)."""
    segs = cfg.course.segments
    heads = cfg.course.headings
    wind = cfg.env.wind_mps
    full = cfg.build.full_throttle_flights
    if full:
        v_cruise = float("inf")
    t_tot = 0.0
    e_tot = 0.0
    v = v_entry

    # throttle and power needed to hold the sustained turn: the thrust that
    # equals the drag, which itself depends (a little) on that thrust
    # through the slipstream's lift - two fixed-point passes settle it
    thrust_turn = aero.drag(turn.v_turn, turn.load_factor * weight_n)
    for _ in range(2):
        thrust_turn = aero.drag(turn.v_turn, turn.load_factor * weight_n,
                                thrust_turn)
    # (or, flying at full throttle throughout, full throttle: the turn is
    # flown as the sustained turn, but charged as full power)
    thr = (1.0 if full else
           ps.throttle_for_thrust(turn.v_turn, soc, thrust_turn))
    p_turn = ps.power_elec(turn.v_turn, soc, thr)
    q_turn = ps.motor_loss(turn.v_turn, soc, thr)
    qp_turn = ps.pack_heat(turn.v_turn, soc, thr)
    heat: List[float] = []
    pheat: List[float] = []

    for i, (kind, val) in enumerate(segs):
        if kind == "straight":
            wind_along = wind * heads[i]        # +1 upwind, -1 downwind, 0 cross
            nxt = segs[(i + 1) % len(segs)]
            v_target = turn.v_turn if nxt[0] == "turn" else 1e9
            dt_, de_, v = fly_straight(aero, ps, mass_kg, weight_n, soc,
                                       val, v, v_target, wind_along, dx,
                                       v_cruise, heat, pheat,
                                       full_power=full)
            t_tot += dt_
            e_tot += de_
        else:
            dt_ = math.radians(val) / max(turn.turn_rate_rad_s, 1e-6)
            t_tot += dt_
            e_tot += p_turn * dt_
            heat.append(q_turn * dt_)
            pheat.append(qp_turn * dt_)
            v = turn.v_turn
    if loss_out is not None:
        loss_out.append(sum(heat))
    if pack_out is not None:
        pack_out.append(sum(pheat))
    return t_tot, e_tot, v


def _with_turn(turn: TurnCapability, res: MissionResult) -> MissionResult:
    """A mission's result with its turn filled in (m119)."""
    if turn.turn_rate_rad_s > 1e-3:
        res.turn_radius_m = turn.radius_m
        res.turn_load_factor = turn.load_factor
        res.turn_bank_deg = turn.bank_deg
        res.turn_rate_deg_s = math.degrees(turn.turn_rate_rad_s)
    return res


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
                energy_reserve: Optional[float] = None,
                v_cruise: float = float("inf"),
                hold_back: float = 0.0) -> MissionResult:
    """Common engine for M1/M2/M3.

    required_laps set  -> fly exactly that many, report the time (Mission 2).
    required_laps None -> fly as many as fit in the window (Mission 3).

    Battery limits, all hard:
      * the depth of discharge - takeoff, laps and landing together may
        take at most cfg.build.depth_of_discharge of the rated energy (95%:
        every mission flies a fresh pack), less cfg.build.energy_reserve
        if one is set;
      * the landing point - no lap is started that would take any pack
        past cfg.build.land_at_discharge of its energy (88%): the flight
        turns for landing there and lands out of what is left;
      * the landing is budgeted BEFORE the last lap is started, so a
        flight never flies a lap it then cannot land from;
      * the ESC low-voltage cutoff - once the LOADED pack voltage at full
        throttle in the turn would fall below pack.v_cutoff_cell per cell,
        the ESC would pull the power, so the flight ends there even if
        the energy budget is not yet spent.
    """
    ps = ac.prop_sys
    mass = weight_n / cfg.env.g
    # Every pack on the aeroplane, not one of them.  With two propulsion
    # systems the packs discharge together, so state of charge must be
    # measured against the total; dividing by a single pack made the
    # voltage sag twice as fast as it should, understated thrust for the
    # whole second half of a mission, and reported a NEGATIVE energy
    # margin for a flight that had in fact finished with energy to spare.
    pack_j = max(ac.total_energy_wh * 3600.0, 1e-6)
    if energy_reserve is None:
        energy_reserve = cfg.build.energy_reserve
    usable = ac.usable_energy_j() * (1.0 - max(0.0, min(energy_reserve, 0.9)))
    # a planning margin held back on top (Mission 2: m2_energy_margin)
    usable -= max(hold_back, 0.0) * pack_j
    window = cfg.rules.mission_window_s
    v_cut = ac.pack.v_cutoff_cell * ac.pack.series

    # Full throttle from a full pack is the heaviest load the flight ever
    # puts on the battery.  If that alone sags the pack below the ESC
    # cutoff, the pack cannot fly this aeroplane.  Until the aeroplane is
    # off the ground no heat is tracked, so the pack and windings are
    # reported at air temperature.
    if ps.operating(0.0, 1.0, 1.0)[3] < v_cut:
        return MissionResult(False, 0, 0.0, 0.0, 1.0,
                             reason="pack sags below the ESC low-voltage "
                                    "cutoff at full throttle on a full "
                                    "charge", stopped_by="battery",
                             motor_temp_c=cfg.env.temp_c,
                             pack_temp_c=cfg.env.temp_c)

    v_s = aero.v_stall(weight_n)
    to = takeoff(aero, ps, mass, weight_n, cfg, soc=1.0,
                 gear_height_m=ac.gear_height_m)
    if not to.ok:
        return MissionResult(False, 0, 0.0, to.energy_j, 1.0,
                             v_stall=v_s, takeoff_m=to.ground_roll_m,
                             reason=to.reason or "takeoff",
                             motor_temp_c=cfg.env.temp_c,
                             pack_temp_c=cfg.env.temp_c)

    # RULE 3.2.1c - every joule of the takeoff comes out of the same
    # propulsion pack budget the laps draw on.  There is no external
    # assist and no separate launch energy.
    energy = to.energy_j
    soc = max(0.03, 1.0 - energy / pack_j)
    cl = climb(aero, ps, mass, weight_n, cfg, soc, to.v_liftoff)
    if cl is None:
        return MissionResult(False, 0, to.time_s, energy, 1.0,
                             v_stall=v_s, takeoff_m=to.ground_roll_m,
                             reason="cannot climb to pattern altitude at "
                                    "full throttle",
                             motor_temp_c=cfg.env.temp_c,
                             pack_temp_c=cfg.env.temp_c)
    t_cl, e_cl, v_cl, x_cl = cl
    energy += e_cl
    t_now = to.time_s + t_cl

    # Motor heating (propulsion.motor_thermal): the windings start at air
    # temperature, heat with their losses, cool to the prop wash, and every
    # lap is flown on copper at the temperature it has reached.  Takeoff
    # and climb are at full throttle, so their heat is their energy times
    # the full-throttle loss fraction.
    b = cfg.build
    n_mot = max(ac.design.n_motors, 1)
    r_th, c_th = motor_thermal(ac.motor, b)
    t_air = cfg.env.temp_c

    def loss_frac(v_: float) -> float:
        op_ = ps.solve(v_, 1.0, 1.0)
        return op_.p_motor_loss_w / max(op_.p_elec_w, 1e-9)

    q_init = (to.energy_j * loss_frac(0.5 * to.v_liftoff)
              + e_cl * loss_frac(v_cl))
    t_wind = winding_step(t_air, t_air, q_init / n_mot / max(t_now, 1e-3),
                          t_now, r_th, c_th)
    t_peak = t_wind

    # Pack heating: each pack heats with its own I^2 R (the DC resistance,
    # polarisation included) and, over a five-minute flight inside the
    # fuselage, sheds little of it - so none is credited.  Takeoff and
    # climb are at full throttle.
    c_pack = max(ac.pack.mass_kg * b.pack_specific_heat_j_kgk, 1e-6)

    def pack_rate(v_: float) -> float:
        return ps.pack_heat(v_, 1.0, 1.0)

    t_pack = (t_air + (pack_rate(0.5 * to.v_liftoff) * to.time_s
                       + pack_rate(v_cl) * t_cl) / c_pack)

    turn = turn_capability(aero, ps, weight_n, 0.70, cfg,
                           getattr(ac, 'n_limit', None))
    if turn.turn_rate_rad_s <= 1e-3:
        return _with_turn(turn, MissionResult(False, 0, t_now, energy, 0.0, v_stall=v_s,
                                              takeoff_m=to.ground_roll_m,
                                              reason="cannot sustain a level turn",
                                              motor_temp_c=t_peak, pack_temp_c=t_pack))
    v_max = max_level_speed(aero, ps, weight_n, 0.70)

    # The landing is paid for out of the same budget, so hold it back
    # from the laps up front.  Estimated at the deepest charge the flight
    # can reach, where the pack voltage is lowest and the energy for the
    # same descent power is the largest - so the estimate is conservative.
    soc_floor = max(0.03, 1.0 - usable / pack_j)
    e_land_budget = landing(aero, ps, mass, weight_n, cfg, soc_floor,
                            ac.gear_height_m)[1]
    lap_budget = usable - e_land_budget
    # The landing point: turn for landing once any pack has used
    # land_at_discharge of its energy.  Every system carries one identical
    # pack at the same current, so each pack is at the same fraction as
    # the whole - the first pack to reach the point is all of them.
    land_point_j = min(max(cfg.build.land_at_discharge, 0.0), 1.0) * pack_j
    land_point_binds = land_point_j < lap_budget
    lap_budget = min(lap_budget, land_point_j)

    lap_times: List[float] = []
    v = v_cl
    laps = 0
    cached: Optional[Tuple[float, float, float, float]] = None
    # (time, energy, motor heat, pack heat)
    cached_soc = 1.0
    cached_ps = None
    credit = x_cl                                   # ground covered during climb
    stopped_by = "cap"

    while True:
        if required_laps is not None and laps >= required_laps:
            stopped_by = "laps"
            break
        if laps >= max_laps:
            break
        soc = max(0.03, 1.0 - energy / pack_j)

        # ESC low-voltage cutoff, checked where the current is highest:
        # full throttle at the turn speed, the slowest part of the lap.
        if ps.operating(turn.v_turn, soc, 1.0)[3] < v_cut:
            if required_laps is not None:
                return _with_turn(turn, MissionResult(False, laps, t_now, energy, 0.0,
                                                      lap_times, v_max, v_s, turn.v_turn,
                                                      turn.limited_by, to.ground_roll_m,
                                                      reason="ESC low-voltage cutoff after "
                                                      "%d of %d laps" % (laps, required_laps),
                                                      motor_temp_c=t_peak,
                                                      pack_temp_c=t_pack,
                                                      stopped_by="battery"))
            stopped_by = "battery"
            break

        # Reuse the converged lap only while the pack state is close to the
        # state it was simulated at, and the windings are still in the
        # same 2 C band; otherwise re-simulate.
        ps_hot = ac.prop_sys_at(t_wind)
        if (cached is not None and abs(cached_soc - soc) <= SOC_RESIM_THRESHOLD
                and cached_ps is ps_hot):
            lt, le, lq, lqp = cached
        else:
            heat_lap: List[float] = []
            pheat_lap: List[float] = []
            lt, le, v = simulate_lap(aero, ps_hot, mass, weight_n, soc, cfg,
                                     turn, v, dx, v_cruise, heat_lap,
                                     pheat_lap)
            lq, lqp = heat_lap[0], pheat_lap[0]
            cached, cached_soc, cached_ps = (lt, le, lq, lqp), soc, ps_hot

        # No pilot flies the theoretical lap.  One factor, applied to time
        # and to the energy that time costs.
        eff = max(min(cfg.build.pilot_efficiency, 1.0), 0.5)
        lt, le, lq, lqp = lt / eff, le / eff, lq / eff, lqp / eff

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

        if energy + le > lap_budget:
            if required_laps is not None:
                return _with_turn(turn, MissionResult(False, laps, t_now, energy, 0.0,
                                                      lap_times, v_max, v_s, turn.v_turn,
                                                      turn.limited_by, to.ground_roll_m,
                                                      reason=("a pack would pass the %.0f%% "
                                                              "landing point before %d laps"
                                                              % (100.0
                                                                 * cfg.build.land_at_discharge,
                                                                 required_laps)
                                                              if land_point_binds else
                                                              "out of battery before %d laps"
                                                              % required_laps),
                                                      motor_temp_c=t_peak,
                                                      pack_temp_c=t_pack,
                                                      stopped_by="battery"))
            stopped_by = "battery"
            break
        if required_laps is None and t_now + lt > window:
            stopped_by = "window"
            break
        # the windings after this lap: past the limit, the lap is not flown
        t_next = winding_step(t_wind, t_air, lq / n_mot / max(lt, 1e-3),
                              lt, r_th, c_th)
        if t_next > b.motor_temp_limit_c:
            if required_laps is not None:
                return _with_turn(turn, MissionResult(False, laps, t_now, energy, 0.0,
                                                      lap_times, v_max, v_s, turn.v_turn,
                                                      turn.limited_by, to.ground_roll_m,
                                                      reason="motor windings would reach "
                                                      "%.0f C on lap %d (limit %.0f C)"
                                                      % (t_next, laps + 1,
                                                         b.motor_temp_limit_c),
                                                      motor_temp_c=t_next,
                                                      pack_temp_c=t_pack,
                                                      stopped_by="motor heat"))
            stopped_by = "motor heat"
            break
        t_wind = t_next
        t_peak = max(t_peak, t_wind)
        # the pack after this lap: past its limit, the lap is not flown
        t_pack_next = t_pack + lqp / c_pack
        if t_pack_next > b.pack_temp_limit_c:
            if required_laps is not None:
                return _with_turn(turn, MissionResult(False, laps, t_now, energy, 0.0,
                                                      lap_times, v_max, v_s, turn.v_turn,
                                                      turn.limited_by, to.ground_roll_m,
                                                      reason="battery pack would reach "
                                                      "%.0f C on lap %d (limit %.0f C)"
                                                      % (t_pack_next, laps + 1,
                                                         b.pack_temp_limit_c),
                                                      motor_temp_c=t_peak,
                                                      pack_temp_c=t_pack_next,
                                                      stopped_by="pack heat"))
            stopped_by = "pack heat"
            break
        t_pack = t_pack_next

        t_now += lt
        energy += le
        laps += 1
        lap_times.append(lt)

    if required_laps is not None and laps < required_laps:
        return _with_turn(turn, MissionResult(False, laps, t_now, energy, 0.0, lap_times,
                                              v_max, v_s, turn.v_turn, turn.limited_by,
                                              to.ground_roll_m,
                                              reason="could not complete required laps",
                                              motor_temp_c=t_peak, pack_temp_c=t_pack))
    if required_laps is not None and t_now > window:
        return _with_turn(turn, MissionResult(False, laps, t_now, energy, 0.0, lap_times,
                                              v_max, v_s, turn.v_turn, turn.limited_by,
                                              to.ground_roll_m,
                                              reason="%.1f s exceeds the %.0f s window"
                                              % (t_now, window),
                                              motor_temp_c=t_peak, pack_temp_c=t_pack))
    if laps == 0:
        return _with_turn(turn, MissionResult(False, 0, t_now, energy, 0.0, lap_times,
                                              v_max, v_s, turn.v_turn, turn.limited_by,
                                              to.ground_roll_m,
                                              reason="no laps completed in the window",
                                              motor_temp_c=t_peak, pack_temp_c=t_pack))

    # Landing: required to score, charged to the battery, reported
    # separately, and excluded from the scored window per the rules.
    soc = max(0.03, 1.0 - energy / pack_j)
    t_land, e_land, roll = landing(aero, ps, mass, weight_n, cfg, soc,
                                   ac.gear_height_m)
    if energy + e_land > usable * (1.0 + 1e-9):
        return _with_turn(turn, MissionResult(False, laps, t_now, energy + e_land, 0.0,
                                              lap_times, v_max, v_s, turn.v_turn,
                                              turn.limited_by, to.ground_roll_m,
                                              t_land, roll, t_now + t_land,
                                              reason="not enough battery left to land safely",
                                              stopped_by=stopped_by, motor_temp_c=t_peak,
                                              pack_temp_c=t_pack))
    energy += e_land
    margin = 1.0 - energy / pack_j
    return _with_turn(turn, MissionResult(True, laps, t_now, energy, margin, lap_times,
                                          v_max, v_s, turn.v_turn, turn.limited_by,
                                          to.ground_roll_m, t_land, roll, t_now + t_land,
                                          stopped_by=stopped_by, motor_temp_c=t_peak,
                                          pack_temp_c=t_pack))


# ==========================================================================
# Choosing a cruise speed
# ==========================================================================
# Since 2026-10-02 every flight is planned at full throttle end to end
# (BuildStandard.full_throttle_flights), a factor of safety: both searches
# below then fly the full-throttle flight and nothing else.  What follows
# is why they exist, for when that is switched off.
#
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
#: Bisections between the first completing probe and the failing one above
#: it.  Four take the cruise speed to 1/64 of the speed range (~0.3 m/s):
#: the fastest the budget allows, found continuously rather than dropping
#: a whole probe step.
_CRUISE_REFINE = 4


def _fastest_that_finishes(ac, aero, weight_n, cfg, required_laps, dx,
                           v_max, hold_back: float = 0.0) -> MissionResult:
    """Fly as fast as the energy budget allows.

    Tries full throttle first, because when it works it is optimal and
    costs one simulation.  Only if that runs the pack dry does it search
    downwards: coarse probes to find a cruise speed that completes, then
    bisection towards the failing speed above it, so the answer is the
    fastest speed that still finishes - continuous in the design.

    With only the probes, a design one gram too heavy for the first probe
    dropped straight to the second, which on a heavy aeroplane misses the
    time window: the score fell off a cliff, and the search parked designs
    on the lip of it with ~0.001 Wh to spare.  Now a heavier design simply
    cruises a little slower.
    """
    def fly(v):
        return fly_mission(ac, aero, weight_n, cfg,
                           required_laps=required_laps, dx=dx, v_cruise=v,
                           hold_back=hold_back)

    def too_slow(r) -> bool:
        return "exceeds the" in r.reason        # missed the time window

    best = fly(float("inf"))
    if best.ok or cfg.build.full_throttle_flights:
        return best
    v_s = aero.v_stall(weight_n)
    lo, hi = 1.15 * v_s, max(v_max, 1.2 * v_s)
    found, v_ok, v_bad, v_slow = None, None, hi, None
    for i in range(_CRUISE_STEPS):
        # Walk down from the top: the first speed that completes.
        v = hi - (hi - lo) * (i + 1) / (_CRUISE_STEPS + 1)
        r = fly(v)
        if r.ok:
            found, v_ok = r, v
            break
        if too_slow(r):
            v_slow = v
            break
        v_bad = v
    if found is None and v_slow is not None:
        # One probe ran out of energy, the next ran out of time: any speed
        # that finishes lies between them.  Bisect for one.
        a, b = v_slow, v_bad
        for _ in range(_CRUISE_REFINE + 2):
            mid = 0.5 * (a + b)
            r = fly(mid)
            if r.ok:
                found, v_ok = r, mid
                break
            if too_slow(r):
                a = mid
            else:
                b = mid
        v_bad = b
    if found is None:
        return best
    for _ in range(_CRUISE_REFINE):
        mid = 0.5 * (v_ok + v_bad)
        r = fly(mid)
        if r.ok:
            found, v_ok = r, mid
        else:
            v_bad = mid
    return found


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
    # Stopped by the clock (or the lap cap): full throttle already got
    # every lap there was.  This used to test energy_margin > 0.15, but
    # energy_margin counts the whole pack, and a flight that has spent its
    # whole safe budget still has 20% of the pack left - so the test was
    # always true and the slower-cruise search never ran at all.
    if cfg.build.full_throttle_flights or (
            best.ok and best.stopped_by not in ("battery", "motor heat",
                                                "pack heat")):
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
    # RULE 3.3.3.1 - the lap count comes from the rules config like
    # Mission 2's does.  It used to be a literal 3 here, which left
    # `m1_required_laps` - editable in the GUI's Advanced panel - with no
    # effect at all.
    return _fastest_that_finishes(ac, ac.aero_clean, ac.w_m1_n, cfg,
                                  cfg.rules.m1_required_laps, dx,
                                  ac._max_level_speed(ac.w_m1_n))


def mission_2(ac: Aircraft, cfg: Config, dx: float = 6.0) -> MissionResult:
    """Delivery flight: sensor in shipping container + optional simulators,
    5 laps, scored on weight / time."""
    # RULE 3.3.3a - an unfinished Mission 2 also locks Mission 3.  Its
    # margin for model error is the landing point (land_at_discharge);
    # m2_energy_margin, 0 by default, can add more.
    return _fastest_that_finishes(ac, ac.aero_clean, ac.w_m2_n, cfg,
                                  cfg.rules.m2_required_laps, dx,
                                  ac._max_level_speed(ac.w_m2_n),
                                  hold_back=cfg.build.m2_energy_margin)


def mission_3(ac: Aircraft, cfg: Config, dx: float = 6.0) -> MissionResult:
    """Sensor flight: sensor deployed for the whole pattern, maximise laps."""
    return _most_laps(ac, ac.aero_deployed, ac.w_m3_n, cfg, dx,
                      ac._max_level_speed(ac.w_m3_n))
