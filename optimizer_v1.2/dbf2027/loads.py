"""Design load factors, derived rather than assumed.

The spar used to be sized to a hand-waved 4.5 g.  This builds a V-n diagram
for each aeroplane instead and takes the limit load factor from whichever of
three things actually governs it:

1.  **Manoeuvre.**  The mission itself demands a sustained turn at the bank
    limit, n = 1/cos(phi).  A pilot flying for time overshoots, so that is
    multiplied by an overshoot factor.

2.  **Gust.**  The classic Pratt discrete-gust formula.  Light wing loading
    makes model aircraft *more* gust-sensitive than full-scale, so this term
    frequently governs a lightly loaded DBF wing:

        mu  = 2 (W/S) / (rho c a g)              mass ratio
        Kg  = 0.88 mu / (5.3 + mu)               gust alleviation
        n   = 1 + Kg rho V a U / (2 W/S)

    evaluated at cruise with the rough-air gust and at the design dive speed
    with half of it, taking whichever is worse.

3.  **Aerodynamic ceiling.**  You cannot pull more g than the wing can
    generate.  At the design dive speed the wing stalls at

        n_aero = q_D S CLmax / W

    so the limit load factor is capped there - designing a spar for a load
    the wing physically cannot produce is wasted weight.

Limit load is the largest of (1) and (2), capped by (3) and floored for
sanity.  Ultimate = 1.5 x limit, the standard factor of safety.

If tech inspection specifies a wing test, set `BuildStandard.
ultimate_load_factor_override` and this whole module is bypassed.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from .units import FT2M


@dataclass
class LoadCase:
    n_manoeuvre: float       # what the flown turns demand, with overshoot
    n_gust_cruise: float     # Pratt gust at cruise
    n_gust_dive: float       # Pratt gust at the design dive speed
    n_aero_ceiling: float    # most the wing can generate at V_D
    n_limit: float           # governing limit load factor
    n_ultimate: float        # 1.5 x limit
    v_dive: float
    driver: str              # which case governs

    @property
    def summary(self) -> str:
        return ("limit %.2f g (%s), ultimate %.2f g"
                % (self.n_limit, self.driver, self.n_ultimate))


def lift_curve_slope_3d(aspect_ratio: float) -> float:
    """Finite-wing lift-curve slope, per radian (Helmbold / Prandtl)."""
    a0 = 2.0 * math.pi
    return a0 * aspect_ratio / (aspect_ratio + 2.0)


def gust_load_factor(w_over_s: float, chord_m: float, cl_alpha: float,
                     v_mps: float, rho: float, u_gust_mps: float) -> float:
    """Pratt discrete-gust load factor."""
    w_over_s = max(w_over_s, 1e-6)
    mu = 2.0 * w_over_s / (rho * chord_m * cl_alpha * 9.80665)
    kg = 0.88 * mu / (5.3 + mu)
    return 1.0 + kg * rho * v_mps * cl_alpha * u_gust_mps / (2.0 * w_over_s)


def design_load_factors(aero, weight_n: float, v_max_level: float,
                        build) -> LoadCase:
    """Build the V-n diagram and return the governing load case."""
    rho = aero.env.rho
    s = aero.S
    chord = aero.wing.chord_m
    w_over_s = weight_n / max(s, 1e-9)
    a = lift_curve_slope_3d(aero.AR)

    v_c = max(v_max_level, 5.0)
    v_d = build.dive_speed_factor * v_c

    # 1. manoeuvre: the sustained turn the mission asks for, plus overshoot
    n_bank = 1.0 / math.cos(math.radians(build.max_bank_deg))
    n_man = n_bank * build.pilot_overshoot_factor

    # 2. gust, at cruise and at the dive speed
    n_gc = gust_load_factor(w_over_s, chord, a, v_c, rho,
                            build.gust_cruise_mps)
    n_gd = gust_load_factor(w_over_s, chord, a, v_d, rho,
                            build.gust_dive_mps)

    # 3. aerodynamic ceiling at the dive speed
    cl_max = aero.cl_max_at_v(v_d)
    n_aero = 0.5 * rho * v_d * v_d * s * cl_max / max(weight_n, 1e-9)

    candidates = [(n_man, "manoeuvre"), (n_gc, "gust at cruise"),
                  (n_gd, "gust at dive speed"),
                  (build.load_factor_floor, "floor")]
    n_limit, driver = max(candidates, key=lambda t: t[0])
    if n_aero < n_limit:
        n_limit, driver = n_aero, "wing CLmax ceiling"

    return LoadCase(n_man, n_gc, n_gd, n_aero, n_limit,
                    build.load_safety_factor * n_limit, v_d, driver)
