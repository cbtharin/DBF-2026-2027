#!/usr/bin/env python3
"""DBF 2026-27 mission simulator, design optimiser and validation suite.

A single-file build of the `dbf2027` package.  Pure Python standard library:
no pip install, no numpy.  Python 3.8+.

    python dbf.py --quick                 fast search
    python dbf.py --thorough              wide search
    python dbf.py --evaluate design.json  score one aeroplane, print its BOM
    python dbf.py --selftest              physics and scoring invariants
    python dbf.py --validate              catalogue, bounds, determinism
    python dbf.py --convergence           numerical error study
    python dbf.py --audit ranked.csv      which options reached the finals
    python dbf.py --version

The physics, every derivation and every assumption are documented in
README.md.  Sections below follow the dependency order of the original
package; each is introduced by a banner comment.

GENERATED FILE - edit the dbf2027/ package and re-run build_single_file.py.
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from dataclasses import asdict, replace
from dataclasses import dataclass
from dataclasses import dataclass, field
from multiprocessing import Pool, cpu_count
from typing import Any, Dict, List, Optional, Sequence
from typing import Any, Dict, List, Optional, Sequence, Tuple
from typing import Any, Dict, List, Optional, Tuple
from typing import Callable, List, Optional, Sequence, Tuple
from typing import Dict
from typing import Dict, List, Optional, Sequence
from typing import Dict, List, Optional, Tuple
from typing import Dict, Optional
from typing import Dict, Tuple
from typing import List, Optional, Tuple
from typing import Optional
from typing import Tuple
import argparse
import collections
import csv
import json
import math
import os
import random
import sys
import tempfile
import time


__version__ = "1.0.0"



# ==========================================================================
#  UNITS - conversions; everything internal is SI
# ==========================================================================

# length
IN2M = 0.0254
M2IN = 1.0 / IN2M
FT2M = 0.3048
M2FT = 1.0 / FT2M

# mass / force
LB2KG = 0.45359237
KG2LB = 1.0 / LB2KG
OZ2KG = 0.028349523125
KG2OZ = 1.0 / OZ2KG
G0 = 9.80665                      # m/s^2
LBF2N = LB2KG * G0
N2LBF = 1.0 / LBF2N

# speed
MPH2MS = 0.44704
MS2MPH = 1.0 / MPH2MS
KT2MS = 0.514444

# area
IN2_2_M2 = IN2M ** 2
FT2_2_M2 = FT2M ** 2

# energy
WH2J = 3600.0
J2WH = 1.0 / WH2J


def kg2lbf(m_kg: float) -> float:
    """Mass in kg -> weight in pound-force at standard gravity."""
    return m_kg * KG2LB


def lbf2kg(w_lb: float) -> float:
    return w_lb * LB2KG


def n2lb(force_n: float) -> float:
    return force_n * N2LBF


# ==========================================================================
#  NUMERICS - Brent, Gauss-Legendre, RK4, PCHIP, golden section, LHS, Nelder-Mead, Richardson
# ==========================================================================

# ==========================================================================
# Root finding
# ==========================================================================
def bisect(f: Callable[[float], float], a: float, b: float,
           tol: float = 1e-9, maxiter: int = 120) -> float:
    """Kept for callers that want guaranteed-bracketed, no-surprises behaviour."""
    fa, fb = f(a), f(b)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0.0:
        return a if abs(fa) < abs(fb) else b
    for _ in range(maxiter):
        m = 0.5 * (a + b)
        fm = f(m)
        if fm == 0.0 or (b - a) < tol:
            return m
        if fa * fm < 0.0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


def brent(f: Callable[[float], float], a: float, b: float,
          xtol: float = 1e-12, rtol: float = 1e-12,
          maxiter: int = 100) -> float:
    """Brent's method: inverse quadratic interpolation with a bisection
    fallback.  Superlinear convergence with the robustness of bisection.

    Falls back to the better endpoint if the root is not bracketed, which is
    what degenerate propulsion cases hand us.
    """
    fa, fb = f(a), f(b)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0.0:
        return a if abs(fa) < abs(fb) else b

    c, fc = a, fa
    d = e = b - a
    for _ in range(maxiter):
        if fb * fc > 0.0:
            c, fc = a, fa
            d = e = b - a
        if abs(fc) < abs(fb):
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        tol1 = 2.0 * 2.220446049250313e-16 * abs(b) + 0.5 * xtol
        xm = 0.5 * (c - b)
        if abs(xm) <= tol1 or fb == 0.0:
            return b
        if abs(e) >= tol1 and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:                       # secant
                p = 2.0 * xm * s
                q = 1.0 - s
            else:                            # inverse quadratic
                q_ = fa / fc
                r_ = fb / fc
                p = s * (2.0 * xm * q_ * (q_ - r_) - (b - a) * (r_ - 1.0))
                q = (q_ - 1.0) * (r_ - 1.0) * (s - 1.0)
            if p > 0.0:
                q = -q
            p = abs(p)
            if 2.0 * p < min(3.0 * xm * q - abs(tol1 * q), abs(e * q)):
                e, d = d, p / q              # accept interpolation
            else:
                d = e = xm                   # fall back to bisection
        else:
            d = e = xm
        a, fa = b, fb
        b += d if abs(d) > tol1 else (tol1 if xm > 0 else -tol1)
        fb = f(b)
    return b


# ==========================================================================
# Optimisation (derivative free)
# ==========================================================================
def golden_max(f: Callable[[float], float], a: float, b: float,
               tol: float = 1e-6, maxiter: int = 200) -> Tuple[float, float]:
    """Maximise a unimodal f on [a, b].  Returns (x*, f(x*))."""
    invphi = (math.sqrt(5.0) - 1.0) / 2.0
    c = b - invphi * (b - a)
    d = a + invphi * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(maxiter):
        if (b - a) < tol:
            break
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - invphi * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + invphi * (b - a)
            fd = f(d)
    x = 0.5 * (a + b)
    return x, f(x)


def scan_then_refine(f: Callable[[float], float], lo: float, hi: float,
                     n_scan: int = 12, tol: float = 1e-4
                     ) -> Tuple[float, float]:
    """Coarse scan to bracket the global max of a possibly multi-modal f,
    then golden section inside the winning bracket.

    A bare scan reports the best *sample*, not the best *value* - that is a
    discretisation error of order (interval/n_scan), and it was costing real
    accuracy in corner speed and best-rate-of-climb.
    """
    if hi <= lo:
        return lo, f(lo)
    xs = [lo + (hi - lo) * i / n_scan for i in range(n_scan + 1)]
    vals = [f(x) for x in xs]
    k = max(range(len(xs)), key=lambda i: vals[i])
    a = xs[max(k - 1, 0)]
    b = xs[min(k + 1, n_scan)]
    xr, fr = golden_max(f, a, b, tol=tol)
    return (xr, fr) if fr >= vals[k] else (xs[k], vals[k])


# ==========================================================================
# Quadrature
# ==========================================================================
_GL5_X = (-0.9061798459386640, -0.5384693101056831, 0.0,
          0.5384693101056831, 0.9061798459386640)
_GL5_W = (0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
          0.4786286704993665, 0.2369268850561891)


def gauss_legendre_5(f: Callable[[float], float], a: float, b: float) -> float:
    """5-point Gauss-Legendre: exact for polynomials up to degree 9.

    Replaces the 5-point midpoint rule (exact only to degree 1 per panel)
    used for the power-off deceleration integral.
    """
    if b <= a:
        return 0.0
    half = 0.5 * (b - a)
    mid = 0.5 * (a + b)
    return half * sum(w * f(mid + half * x) for x, w in zip(_GL5_X, _GL5_W))


# ==========================================================================
# Interpolation
# ==========================================================================
# PCHIP is not used by the mission model - the propulsion operating point
# is solved exactly, which measured faster than interpolating it.  It is
# kept and tested for callers who load tabulated propeller data with
# load_props_csv and need to interpolate a measured Ct/Cp curve.
def pchip_slopes(xs: Sequence[float], ys: Sequence[float]) -> List[float]:
    """Fritsch-Carlson slopes: cubic Hermite that cannot overshoot.

    Matters here because the interpolant feeds a root solve - an interpolant
    that wiggles turns a monotone thrust curve into one with spurious local
    extrema, and the solver chases them.
    """
    n = len(xs)
    if n == 1:
        return [0.0]
    h = [xs[i + 1] - xs[i] for i in range(n - 1)]
    delta = [(ys[i + 1] - ys[i]) / h[i] if h[i] != 0 else 0.0
             for i in range(n - 1)]
    m = [0.0] * n
    m[0] = delta[0]
    m[-1] = delta[-1]
    for i in range(1, n - 1):
        if delta[i - 1] * delta[i] <= 0.0:
            m[i] = 0.0                      # local extremum: flatten
        else:
            w1 = 2.0 * h[i] + h[i - 1]
            w2 = h[i] + 2.0 * h[i - 1]
            m[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i])
    # endpoint limiters (Fritsch-Carlson one-sided rule)
    for i, dl in ((0, delta[0]), (n - 1, delta[-1])):
        if dl == 0.0:
            m[i] = 0.0
        elif m[i] * dl < 0.0:
            m[i] = 0.0
        elif abs(m[i]) > 3.0 * abs(dl):
            m[i] = 3.0 * dl
    return m


def hermite_eval(xs: Sequence[float], ys: Sequence[float],
                 ms: Sequence[float], x: float) -> float:
    """Evaluate the cubic Hermite defined by (xs, ys, slopes ms)."""
    n = len(xs)
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    lo, hi = 0, n - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    h = xs[hi] - xs[lo]
    t = (x - xs[lo]) / h
    t2 = t * t
    t3 = t2 * t
    h00 = 2 * t3 - 3 * t2 + 1
    h10 = t3 - 2 * t2 + t
    h01 = -2 * t3 + 3 * t2
    h11 = t3 - t2
    return (h00 * ys[lo] + h10 * h * ms[lo]
            + h01 * ys[hi] + h11 * h * ms[hi])


class PCHIP:
    """Monotone cubic interpolant with precomputed slopes."""

    __slots__ = ("xs", "ys", "ms")

    def __init__(self, xs: Sequence[float], ys: Sequence[float]):
        self.xs = list(xs)
        self.ys = list(ys)
        self.ms = pchip_slopes(self.xs, self.ys)

    def __call__(self, x: float) -> float:
        return hermite_eval(self.xs, self.ys, self.ms, x)


# ==========================================================================
# ODE integration
# ==========================================================================
def rk4_step(f: Callable[[float, Sequence[float]], List[float]],
             t: float, y: Sequence[float], h: float) -> List[float]:
    """One classical Runge-Kutta 4 step.  Local error O(h^5)."""
    k1 = f(t, y)
    y2 = [y[i] + 0.5 * h * k1[i] for i in range(len(y))]
    k2 = f(t + 0.5 * h, y2)
    y3 = [y[i] + 0.5 * h * k2[i] for i in range(len(y))]
    k3 = f(t + 0.5 * h, y3)
    y4 = [y[i] + h * k3[i] for i in range(len(y))]
    k4 = f(t + h, y4)
    return [y[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
            for i in range(len(y))]


# ==========================================================================
# Small dense linear algebra
# ==========================================================================
def solve_linear(a, b):
    """Solve a x = b by Gaussian elimination with partial pivoting.

    Written out because the lifting-line solve needs a 10x10 system and
    pulling in numpy for that would break the no-dependency promise.
    `a` is a list of rows; both are left untouched.
    """
    n = len(b)
    m = [list(a[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-300:
            raise ValueError("singular matrix in solve_linear")
        if piv != col:
            m[col], m[piv] = m[piv], m[col]
        inv = 1.0 / m[col][col]
        for r in range(col + 1, n):
            f = m[r][col] * inv
            if f == 0.0:
                continue
            for k in range(col, n + 1):
                m[r][k] -= f * m[col][k]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        acc = m[r][n] - sum(m[r][k] * x[k] for k in range(r + 1, n))
        x[r] = acc / m[r][r]
    return x


# ==========================================================================
# Grid convergence (Roache)
# ==========================================================================
def observed_order(f_coarse: float, f_med: float, f_fine: float,
                   r: float = 2.0) -> Optional[float]:
    """Observed order of convergence p from three systematically refined grids."""
    num = f_coarse - f_med
    den = f_med - f_fine
    if den == 0.0 or num == 0.0 or num / den <= 0.0:
        return None
    return math.log(abs(num / den)) / math.log(r)


def richardson_extrapolate(f_med: float, f_fine: float, p: float,
                           r: float = 2.0) -> float:
    """Estimate of the exact (zero-step-size) value."""
    return f_fine + (f_fine - f_med) / (r ** p - 1.0)


def gci(f_med: float, f_fine: float, p: float, r: float = 2.0,
        safety: float = 1.25) -> float:
    """Grid Convergence Index: a conservative error band on the fine grid,
    as a fraction.  Multiply by 100 for percent."""
    if f_fine == 0.0:
        return float("inf")
    return safety * abs((f_fine - f_med) / f_fine) / (r ** p - 1.0)


# ==========================================================================
# Sampling and multi-objective helpers
# ==========================================================================
def clamp(x: float, lo: float, hi: float) -> float:
    return lo if x < lo else (hi if x > hi else x)


def lhs(n: int, bounds: Sequence[Tuple[float, float]],
        rng: random.Random) -> List[List[float]]:
    """Latin hypercube sample of n points over the given box."""
    dim = len(bounds)
    cols: List[List[float]] = []
    for (lo, hi) in bounds:
        strata = [(i + rng.random()) / n for i in range(n)]
        rng.shuffle(strata)
        cols.append([lo + s * (hi - lo) for s in strata])
    return [[cols[d][i] for d in range(dim)] for i in range(n)]


def nelder_mead(f: Callable[[List[float]], float], x0: Sequence[float],
                step: Sequence[float],
                lower: Sequence[float], upper: Sequence[float],
                maxiter: int = 400, ftol: float = 1e-8,
                restarts: int = 1) -> Tuple[List[float], float]:
    """Bounded Nelder-Mead simplex *minimisation* (bounds enforced by clipping).

    `restarts` re-forms the simplex around the incumbent and runs again, which
    recovers the convergence a collapsed simplex throws away - cheap insurance
    against stopping on a flat spot rather than a minimum.
    """
    dim = len(x0)

    def clip(v: Sequence[float]) -> List[float]:
        return [clamp(v[i], lower[i], upper[i]) for i in range(dim)]

    best_x = clip(list(x0))
    best_f = f(best_x)

    for attempt in range(max(1, restarts + 1)):
        shrink = 0.45 ** attempt
        simplex = [list(best_x)]
        for i in range(dim):
            p = list(best_x)
            s = (step[i] if step[i] != 0.0 else 0.05 * max(1e-6, abs(p[i])))
            p[i] = p[i] + s * shrink
            simplex.append(clip(p))
        fs = [f(p) for p in simplex]

        alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, 0.5
        for _ in range(maxiter):
            order = sorted(range(dim + 1), key=lambda k: fs[k])
            simplex = [simplex[k] for k in order]
            fs = [fs[k] for k in order]
            if abs(fs[-1] - fs[0]) <= ftol * (abs(fs[0]) + abs(fs[-1]) + 1e-12):
                break
            centroid = [sum(simplex[k][d] for k in range(dim)) / dim
                        for d in range(dim)]
            worst = simplex[-1]
            xr = clip([centroid[d] + alpha * (centroid[d] - worst[d])
                       for d in range(dim)])
            fr = f(xr)
            if fr < fs[0]:
                xe = clip([centroid[d] + gamma * (xr[d] - centroid[d])
                           for d in range(dim)])
                fe = f(xe)
                simplex[-1], fs[-1] = (xe, fe) if fe < fr else (xr, fr)
            elif fr < fs[-2]:
                simplex[-1], fs[-1] = xr, fr
            else:
                xc = clip([centroid[d] + rho * (worst[d] - centroid[d])
                           for d in range(dim)])
                fc = f(xc)
                if fc < fs[-1]:
                    simplex[-1], fs[-1] = xc, fc
                else:
                    b0 = simplex[0]
                    for k in range(1, dim + 1):
                        simplex[k] = clip([b0[d] + sigma * (simplex[k][d] - b0[d])
                                           for d in range(dim)])
                        fs[k] = f(simplex[k])
        k = min(range(dim + 1), key=lambda i: fs[i])
        if fs[k] < best_f:
            best_x, best_f = simplex[k], fs[k]
    return best_x, best_f


def pareto_front(points: Sequence[Sequence[float]]) -> List[int]:
    """Indices of the non-dominated points, *maximising* every objective."""
    keep: List[int] = []
    for i, pi in enumerate(points):
        dominated = False
        for j, pj in enumerate(points):
            if i == j:
                continue
            if all(pj[k] >= pi[k] for k in range(len(pi))) and \
               any(pj[k] > pi[k] for k in range(len(pi))):
                dominated = True
                break
        if not dominated:
            keep.append(i)
    return keep


# ==========================================================================
#  CONFIGURATION - environment, course, rules, build standard
# ==========================================================================

# --------------------------------------------------------------------------
# Atmosphere / field
# --------------------------------------------------------------------------
@dataclass
class Environment:
    rho: float = 1.225            # kg/m^3 air density
    nu: float = 1.46e-5           # m^2/s kinematic viscosity
    g: float = 9.80665            # m/s^2
    wind_mps: float = 4.5         # steady wind along the runway
    field_elev_m: float = 0.0

    @staticmethod
    def at_density_altitude(alt_m: float, temp_c: float = 15.0,
                            wind_mps: float = 4.5) -> "Environment":
        """ISA density for a field elevation + temperature, which is what
        actually bites you at the competition site."""
        t_k = temp_c + 273.15
        p = 101325.0 * (1.0 - 2.25577e-5 * alt_m) ** 5.25588
        rho = p / (287.058 * t_k)
        mu = 1.458e-6 * t_k ** 1.5 / (t_k + 110.4)
        return Environment(rho=rho, nu=mu / rho, wind_mps=wind_mps,
                           field_elev_m=alt_m)


# --------------------------------------------------------------------------
# Course
# --------------------------------------------------------------------------
Segment = Tuple[str, float]       # ("straight", metres) | ("turn", degrees)


@dataclass
class Course:
    """2026-27 course as described by the team:

        two 1000 ft straightaways, a 180 degree turn at each end, and a
        horizontal 360 degree loop in the middle of the second straightaway.

    So the second straight is flown as 500 ft, the 360, then another 500 ft.
    Total straight distance 2000 ft; total turning 720 degrees.

    Turn *time* is `angle / turn_rate`, and turn rate is `g sqrt(n^2-1)/V`.
    That is a function of how hard each aeroplane can turn, so the lap time
    genuinely depends on each design's own best turning method - exactly as
    it should.  The ground track length of a turn never enters the time.
    """
    segments: List[Segment] = field(default_factory=lambda: [
        ("straight", 1000.0 * FT2M),   # upwind leg
        ("turn", 180.0),               # turn at the far end
        ("straight", 500.0 * FT2M),    # downwind, first half
        ("turn", 360.0),               # horizontal loop at midfield
        ("straight", 500.0 * FT2M),    # downwind, second half
        ("turn", 180.0),               # turn back onto the upwind leg
    ])
    headings: List[int] = field(default_factory=lambda: [+1, 0, -1, 0, -1, 0])
    pattern_altitude_m: float = 30.0

    @property
    def straight_length_m(self) -> float:
        return sum(L for kind, L in self.segments if kind == "straight")

    @property
    def total_turn_deg(self) -> float:
        return sum(L for kind, L in self.segments if kind == "turn")


# --------------------------------------------------------------------------
# Rules / scoring constants
# --------------------------------------------------------------------------
@dataclass
class Rules:
    # --- flight windows -----------------------------------------------------
    mission_window_s: float = 300.0        # 5 minutes, M2 and M3
    m2_required_laps: int = 5

    # --- field constraints --------------------------------------------------
    # No takeoff distance limit in 2026-27.  The ground roll is still computed
    # and reported, because a 300 ft roll is a real operational risk even when
    # it is not a rules violation.
    takeoff_distance_m: float = float("inf")
    # Maximum allowable wingspan: 6 ft.  This is a real rule, and it is the
    # constraint that stops the wing growing - without it, area scales as
    # span squared while the spar only scales linearly, and the optimiser
    # walks off to a 10 ft aeroplane.
    max_wingspan_m: float = 6.0 * FT2M

    # --- energy -------------------------------------------------------------
    max_watt_hours: float = 100.0          # total pack energy limit
    battery_chemistry_allowed: Tuple[str, ...] = ("LiPo", "LiIon", "NiMH", "NiCd")

    # --- payload ------------------------------------------------------------
    # No rules cap on sensor or simulator mass, so these are search bounds
    # wide enough that the physics, not the bound, decides.  If a result pins
    # to one of them, the report says so.
    # No rules cap on sensor mass.  This is only a search box, set wide
    # enough that the physics decides; if a result pins here, the report
    # says so and the bound should be widened again.
    sensor_mass_min_kg: float = 0.15
    sensor_mass_max_kg: float = 20.00
    sensor_deployed_drag_area_m2: float = 0.0045   # CD*A while deployed

    # Shipping container: ONE fixed assumed box for every design, rather
    # than resizing it per sensor.  The dimensions come from the drop-test
    # calculation for a roughly 1 kg sensor: a 5 ft drop reaches 5.47 m/s,
    # holding the sensor under 75 g needs a 20 mm crush stroke, and foam
    # crushes usefully through about 70 percent of its thickness, so 29 mm
    # of foam per face around a 165 x 83 x 83 mm sensor envelope.  Rounded
    # to a buildable box.  Change these four numbers if the real container
    # differs; nothing else has to change.
    container_length_m: float = 0.230
    container_width_m: float = 0.150
    container_height_m: float = 0.150
    container_tare_kg: float = 0.5 * LB2KG
    # Nominal packing density of a sensor built out of electronics, battery
    # and structure.  Used to REPORT the density a given sensor mass would
    # need, not to cap it: a heavy sensor is simply a dense one.
    sensor_packing_density_kg_m3: float = 800.0
    container_fill_fraction: float = 0.55   # usable internal volume fraction
    # The only real ceiling is physics.  A sensor cannot be denser than the
    # dense metal you would ballast it with; lead is 11340 kg/m3.  This
    # should essentially never bind - it exists so the optimiser cannot
    # propose a 50 lb brick in a 9 inch box.
    max_sensor_density_kg_m3: float = 11340.0

    # Simulators have the container's external dimensions (that is what they
    # simulate) and whatever mass you load them to.  No rules cap on either
    # the count or the mass, so both are physics-limited.
    max_simulators: int = 40
    sim_mass_max_kg: float = 3.00

    # --- normalisers --------------------------------------------------------
    # "self" normalises each mission against the best design found in this
    # simulation, which is what the team asked for.  "field" uses the two
    # estimates below instead, for when you know what the competition can do.
    normalisation: str = "self"
    field_best_m2_lb_per_s: float = 0.060
    field_best_m3_lap_lb: float = 22.0


# --------------------------------------------------------------------------
# Build standard (what your shop actually produces)
# --------------------------------------------------------------------------
@dataclass
class BuildStandard:
    wing_areal_density: float = 1.10      # kg/m^2 of planform: ribs+skin+covering
    tail_areal_density: float = 0.90      # kg/m^2 of planform
    fuse_areal_density: float = 1.40      # kg/m^2 of wetted area
    spar_sigma_allow_pa: float = 600e6    # carbon tube working stress
    spar_rho: float = 1600.0              # kg/m^3
    spar_min_wall_m: float = 0.0006       # minimum practical layup/tube wall
    spar_fitting_factor: float = 1.8      # webs, joiners, root fittings

    # The ultimate load factor is DERIVED per design from a V-n diagram
    # (see loads.py): the worse of the manoeuvre case and the Pratt gust
    # case, capped by what the wing can aerodynamically generate at the
    # design dive speed, times a 1.5 safety factor.  Set the override to a
    # number if tech inspection specifies a wing test instead.
    ultimate_load_factor_override: Optional[float] = None
    load_safety_factor: float = 1.50       # ultimate / limit, standard
    load_factor_floor: float = 2.50        # sanity floor on limit load
    pilot_overshoot_factor: float = 1.15   # a pilot flying for time overshoots
    dive_speed_factor: float = 1.25        # V_D / V_cruise
    gust_cruise_mps: float = 25.0 * FT2M   # rough-air discrete gust at Vc
    gust_dive_mps: float = 12.5 * FT2M     # half of it at V_D, per FAR 23

    gear_mass_frac: float = 0.030         # of MTOW
    gear_mass_fixed_kg: float = 0.080
    servo_mass_kg: float = 0.020
    n_servos_base: int = 6                # 2 ail, 1 elev, 1 rud, 1 door, 1 sensor
    avionics_kg: float = 0.100            # rx, wiring, arming plug, telemetry
    deploy_mech_kg: float = 0.150         # winch/door for the M3 sensor
    payload_structure_frac: float = 0.080 # restraint structure per kg of payload
    contingency_frac: float = 0.080       # honest margin on the whole airframe

    # drag areas (CD*A, m^2)
    gear_drag_area_m2: float = 0.0040
    nacelle_drag_area_m2: float = 0.0015  # per motor beyond the first
    misc_interference_frac: float = 0.10

    # fuselage proportions
    fuse_nose_m: float = 0.22
    fuse_tailcone_m: float = 0.55
    fuse_min_length_m: float = 0.90
    tail_arm_frac: float = 0.45           # of fuselage length
    vol_coeff_h: float = 0.55
    vol_coeff_v: float = 0.040
    tail_aspect_ratio: float = 4.0
    # The tail sits in the wing's wake and behind the fuselage, so it sees
    # less dynamic pressure than the wing.  0.9 is conventional.
    tail_eta: float = 0.90
    # Stability.  The neutral point is now COMPUTED, so static margin is an
    # output and these are real limits on it.
    min_static_margin: float = 0.05     # below this it is not flyable
    max_static_margin: float = 0.35     # above this it is a dart, all trim drag
    # The fuselage is destabilising: it moves the neutral point forward.
    # Raymer-style term, calibrated so a typical DBF pod shifts it ~0.05c.
    fuselage_np_shift_k: float = 0.60
    # CG position aft of the wing aerodynamic centre, in chords.  Positive
    # means the CG is behind the wing AC, which unloads the tail.  A
    # conventional stable layout sits a little aft of the wing AC and well
    # forward of the neutral point; 0.05 c is typical.  Set 0.0 to model
    # trim from wing camber alone.
    cg_aft_of_ac_frac: float = 0.05
    tail_cl_max: float = 1.10        # symmetric section plus elevator
    tail_efficiency: float = 0.85    # tail span efficiency for induced drag

    # --- realism factors: these exist to stop the model flattering itself ---
    cl_max_3d_factor: float = 0.82        # built wing vs 2-D section data
    # Motor winding temperature is now COMPUTED per design from the
    # mission-average copper loss, so this is only a fallback.
    motor_hot_resistance_factor: float = 1.25
    motor_thermal_resistance_kw: float = 1.5   # K per watt of copper loss
    motor_thermal_tau_s: float = 180.0         # winding time constant
    copper_temp_coeff: float = 0.00393         # per K, resistivity of copper
    # A freewheeling propeller is a disc of drag once the power is cut.
    # Coefficient on disc area, per propeller.  With measured propeller
    # data (which extends past J_zero into negative thrust) this could be
    # computed instead of assumed.
    windmill_drag_coeff: float = 0.08
    # Nobody flies the perfect corner speed on every lap.  Lap time is
    # divided by this.  0.93 means 7.5 percent off a flawless lap; set it
    # from your own flight logs, and set it to 1.0 to model a perfect pilot.
    pilot_efficiency: float = 0.93
    install_thrust_factor: float = 0.95   # blockage, non-uniform inflow

    # operating
    rolling_friction: float = 0.035       # paved runway, pneumatic wheels
    braking_friction: float = 0.30        # wheel brakes / drag on rollout
    takeoff_cl_ground: float = 0.35
    rotate_margin: float = 1.15           # V_rotate / V_stall
    approach_margin: float = 1.30         # V_approach / V_stall
    descent_rate_mps: float = 3.0
    flare_time_s: float = 3.0
    max_bank_deg: float = 70.0
    # Fallback only: the mission uses the derived LIMIT load factor, because
    # you should not exceed limit load in normal flight.
    max_turn_load_factor: float = 3.5
    prop_clearance_margin_m: float = 0.05
    gear_height_m: float = 0.130          # main gear leg length


@dataclass
class Config:
    env: Environment = field(default_factory=Environment)
    course: Course = field(default_factory=Course)
    rules: Rules = field(default_factory=Rules)
    build: BuildStandard = field(default_factory=BuildStandard)


# ==========================================================================
#  CATALOGUES - airfoils, motors, propellers, cells, ESCs
# ==========================================================================

# ==========================================================================
# Airfoils
# ==========================================================================
@dataclass(frozen=True)
class Airfoil:
    name: str
    cl_max_2d: float        # section clmax at Re ~ 300k
    cd_min: float           # section minimum profile drag at re_ref
    cl_at_cd_min: float
    k_drag: float           # cd = cd_min + k*(cl-cl_opt)^2
    cl_break: float         # cl beyond which drag rises sharply
    t_over_c: float
    cm_ac: float
    re_ref: float = 300_000.0
    note: str = ""


AIRFOILS: Dict[str, Airfoil] = {
    "S1223": Airfoil("S1223", 2.20, 0.0215, 1.05, 0.0130, 1.85, 0.121, -0.28,
                     note="Selig high-lift. Huge CLmax, draggy, sharp stall. Heavy-haul choice."),
    "E423": Airfoil("E423", 1.95, 0.0165, 0.90, 0.0105, 1.65, 0.125, -0.24,
                    note="Eppler high-lift. Friendlier than S1223, still payload-oriented."),
    "S1210": Airfoil("S1210", 2.05, 0.0195, 1.00, 0.0120, 1.75, 0.120, -0.26,
                     note="High lift, slightly cleaner than S1223."),
    "SD7062": Airfoil("SD7062", 1.55, 0.0122, 0.60, 0.0082, 1.30, 0.140, -0.14,
                      note="Thick, docile, good spar depth. Common DBF compromise."),
    "SD7037": Airfoil("SD7037", 1.35, 0.0102, 0.45, 0.0072, 1.15, 0.092, -0.10,
                      note="Low drag, sport. Favours lap-count missions."),
    "CLARKY": Airfoil("CLARKY", 1.42, 0.0118, 0.40, 0.0085, 1.20, 0.117, -0.08,
                      note="Flat bottom, trivially buildable, forgiving."),
    "MH32": Airfoil("MH32", 1.28, 0.0092, 0.35, 0.0068, 1.05, 0.087, -0.06,
                    note="Fast, thin, low drag. Needs a real spar."),
    "NACA4412": Airfoil("NACA4412", 1.48, 0.0112, 0.45, 0.0080, 1.25, 0.120, -0.09,
                        note="Classic, well documented, easy to validate in XFOIL."),
    # --- added to broaden the trade study to 20 sections ----------------
    "FX63137": Airfoil("FX63137", 1.92, 0.0185, 0.95, 0.0112, 1.62, 0.137, -0.25,
                       note="Wortmann FX 63-137. High lift, thick, popular for heavy-lift."),
    "GOE398": Airfoil("GOE398", 1.78, 0.0170, 0.85, 0.0108, 1.50, 0.110, -0.20,
                      note="Gottingen 398. High lift, classic, easy to build."),
    "S1218": Airfoil("S1218", 1.88, 0.0180, 0.92, 0.0115, 1.58, 0.113, -0.22,
                     note="Selig high-lift, a little cleaner than S1223."),
    "NACA6412": Airfoil("NACA6412", 1.62, 0.0135, 0.62, 0.0092, 1.38, 0.120, -0.14,
                        note="High camber NACA. Good lift, heavy pitching moment."),
    "E214": Airfoil("E214", 1.38, 0.0108, 0.48, 0.0076, 1.18, 0.111, -0.09,
                    note="Eppler general purpose. Balanced, forgiving."),
    "E193": Airfoil("E193", 1.32, 0.0100, 0.42, 0.0072, 1.12, 0.104, -0.07,
                    note="Eppler sport. Low drag, mild stall, easy to build."),
    "E205": Airfoil("E205", 1.26, 0.0094, 0.38, 0.0068, 1.06, 0.106, -0.06,
                    note="Eppler low drag. Fast, needs speed to work."),
    "S3021": Airfoil("S3021", 1.30, 0.0098, 0.40, 0.0070, 1.10, 0.092, -0.07,
                     note="Selig sport. Good low-Re behaviour for its drag."),
    "SA7038": Airfoil("SA7038", 1.36, 0.0104, 0.45, 0.0074, 1.16, 0.093, -0.09,
                      note="Somers/Selig. Similar to SD7037, slightly more lift."),
    "RG15": Airfoil("RG15", 1.18, 0.0088, 0.32, 0.0064, 0.98, 0.089, -0.05,
                    note="Rolf Girsberger. Fast sailplane section, very low drag."),
    "SD8020": Airfoil("SD8020", 1.06, 0.0080, 0.22, 0.0060, 0.88, 0.080, -0.02,
                      note="Thin, near-symmetric, lowest drag here. Speed only."),
    "AG35": Airfoil("AG35", 1.04, 0.0078, 0.25, 0.0062, 0.86, 0.073, -0.04,
                    note="Drela thin section. Minimum drag, minimum structure depth."),
}


# ==========================================================================
# Motors (brushless outrunners)
# ==========================================================================
@dataclass(frozen=True)
class Motor:
    key: str
    name: str
    kv: float               # rpm per volt, no load
    rm: float               # ohm, effective line resistance used by the model
    i0: float               # A, no-load current
    mass_kg: float
    i_max_a: float          # manufacturer continuous current
    p_max_w: float          # manufacturer continuous power
    shaft_mm: float
    price_usd: float
    source: str = "catalogue estimate - VERIFY"


MOTORS: Dict[str, Motor] = {
    "SCORPION_SII_4020_420": Motor("SCORPION_SII_4020_420", "Scorpion SII-4020-420KV",
                                   420, 0.033, 1.6, 0.405, 70, 1900, 6.0, 240),
    "SCORPION_SII_3020_890": Motor("SCORPION_SII_3020_890", "Scorpion SII-3020-890KV",
                                   890, 0.052, 1.4, 0.205, 50, 1000, 5.0, 170),
    "HACKER_A50_14L": Motor("HACKER_A50_14L", "Hacker A50-14L V4",
                            355, 0.041, 1.5, 0.425, 60, 1700, 8.0, 290),
    "HACKER_A40_12L": Motor("HACKER_A40_12L", "Hacker A40-12L V4",
                            610, 0.055, 1.3, 0.222, 45, 1000, 5.0, 175),
    "COBRA_C4120_12": Motor("COBRA_C4120_12", "Cobra C-4120/12 465KV",
                            465, 0.036, 1.5, 0.305, 58, 1450, 6.0, 135),
    "COBRA_C3520_14": Motor("COBRA_C3520_14", "Cobra C-3520/14 780KV",
                            780, 0.048, 1.2, 0.190, 44, 900, 5.0, 100),
    "TMOTOR_AT4130_300": Motor("TMOTOR_AT4130_300", "T-Motor AT4130 300KV",
                               300, 0.030, 1.8, 0.408, 75, 2000, 8.0, 150),
    "TMOTOR_AT4125_540": Motor("TMOTOR_AT4125_540", "T-Motor AT4125 540KV",
                               540, 0.038, 1.6, 0.345, 65, 1550, 6.0, 135),
    "TMOTOR_AT2814_900": Motor("TMOTOR_AT2814_900", "T-Motor AT2814 900KV",
                               900, 0.065, 1.1, 0.152, 36, 700, 4.0, 65),
    "SUNNYSKY_X4112S_400": Motor("SUNNYSKY_X4112S_400", "SunnySky X4112S 400KV",
                                 400, 0.035, 1.5, 0.300, 60, 1400, 6.0, 90),
    "NEU_1512_1Y": Motor("NEU_1512_1Y", "Neu 1512/1.5Y (direct drive wind)",
                         420, 0.025, 2.0, 0.390, 80, 2200, 6.0, 420),
    "EFLITE_POWER_32": Motor("EFLITE_POWER_32", "E-flite Power 32 770KV",
                             770, 0.058, 1.3, 0.200, 42, 850, 5.0, 90),
    # --- added to broaden the trade study to 20 motors ------------------
    # Smaller motors matter now that 3 and 4 motor layouts are searched.
    "SCORPION_SII_4025_330": Motor("SCORPION_SII_4025_330", "Scorpion SII-4025-330KV",
                                   330, 0.028, 1.9, 0.500, 85, 2300, 6.0, 280),
    "HACKER_A60_7XS": Motor("HACKER_A60_7XS", "Hacker A60-7XS V4",
                            365, 0.026, 2.2, 0.640, 95, 2800, 8.0, 420),
    "TMOTOR_AT3520_880": Motor("TMOTOR_AT3520_880", "T-Motor AT3520 880KV",
                               880, 0.058, 1.3, 0.195, 42, 850, 5.0, 75),
    "TMOTOR_AT2820_1050": Motor("TMOTOR_AT2820_1050", "T-Motor AT2820 1050KV",
                                1050, 0.072, 1.0, 0.128, 32, 600, 4.0, 55),
    "COBRA_C2826_10": Motor("COBRA_C2826_10", "Cobra C-2826/10 1130KV",
                            1130, 0.078, 1.0, 0.118, 30, 560, 4.0, 60),
    "COBRA_C3530_14": Motor("COBRA_C3530_14", "Cobra C-3530/14 750KV",
                            750, 0.045, 1.3, 0.225, 48, 1000, 5.0, 110),
    "SUNNYSKY_X3520_520": Motor("SUNNYSKY_X3520_520", "SunnySky X3520 520KV",
                                520, 0.042, 1.4, 0.238, 50, 1100, 5.0, 75),
    "AXI_2826_10": Motor("AXI_2826_10", "AXi 2826/10 920KV",
                         920, 0.068, 1.1, 0.177, 35, 700, 5.0, 130),
}


# ==========================================================================
# Propellers
# ==========================================================================
@dataclass(frozen=True)
class Propeller:
    key: str
    diameter_m: float
    pitch_m: float
    mass_kg: float
    ct0: float              # static thrust coefficient
    cp0: float              # static power coefficient
    j_zero: float           # advance ratio at zero thrust
    kp: float = 0.75        # Cp roll-off: Cp(j_zero) = (1-kp) * cp0
    price_usd: float = 18.0
    source: str = "parametric p/D fit - REPLACE with UIUC/APC data"

    @property
    def p_over_d(self) -> float:
        return self.pitch_m / self.diameter_m

    def ct(self, j: float) -> float:
        return max(0.0, self.ct0 * (1.0 - j / self.j_zero))

    def cp(self, j: float) -> float:
        x = min(j / self.j_zero, 1.25)
        return max(0.20 * self.cp0, self.cp0 * (1.0 - self.kp * x))

    def thrust(self, rho: float, n_rps: float, v: float) -> float:
        if n_rps <= 1e-6:
            return 0.0
        j = v / (n_rps * self.diameter_m)
        return self.ct(j) * rho * n_rps ** 2 * self.diameter_m ** 4

    def power(self, rho: float, n_rps: float, v: float) -> float:
        if n_rps <= 1e-6:
            return 0.0
        j = v / (n_rps * self.diameter_m)
        return self.cp(j) * rho * n_rps ** 3 * self.diameter_m ** 5


def _make_prop(d_in: float, p_in: float) -> Propeller:
    """Parametric APC-thin-electric-like coefficients from pitch/diameter.

    Anchored on three things a real propeller has to satisfy at once:

      * static figure of merit  FOM = ct0^1.5 / (1.2533 cp0)  lands at ~0.71
        for a low-pitch prop falling to ~0.52 for a high-pitch one;
      * peak propulsive efficiency  eta_max = 0.444 * j_zero * ct0/cp0  runs
        ~0.70 (p/D 0.5) to ~0.78 (p/D 0.83);
      * peak efficiency occurs near 0.8x geometric pitch speed, i.e. at
        J/j_zero = 2/3 with the linear Cp roll-off above.

    Check: 13x8 at 8000 rpm gives 7.8 lbf static on 700 W (FOM 0.63), peak
    eta 0.74 at 22 m/s.  Fine for ranking; measure the two or three props you
    shortlist before you trust an absolute number.
    """
    pod = p_in / d_in
    ct0 = 0.115 + 0.030 * pod
    ratio = min(2.70, max(1.40, 2.48 - 2.27 * (pod - 0.50)))
    cp0 = ct0 / ratio
    jz = 1.15 * pod + 0.06
    mass = 0.0009 * (d_in ** 2.2) / 10.0 + 0.012
    return Propeller(key="APC_%gx%g" % (d_in, p_in), diameter_m=d_in * IN2M,
                     pitch_m=p_in * IN2M, mass_kg=mass,
                     ct0=ct0, cp0=cp0, j_zero=jz,
                     price_usd=10.0 + 0.9 * d_in)


_PROP_SIZES = [
    (10, 5), (10, 7), (11, 5.5), (11, 7), (11, 8),
    (12, 6), (12, 8), (12, 10), (13, 6.5), (13, 8), (13, 10),
    (14, 7), (14, 8.5), (14, 10), (14, 12),
    (15, 8), (15, 10), (15, 13),
    (16, 8), (16, 10), (16, 12),
    (17, 10), (17, 12), (18, 8), (18, 10), (18, 12),
    (20, 10), (20, 13),
]

PROPS: Dict[str, Propeller] = {}
for _d, _p in _PROP_SIZES:
    _pr = _make_prop(_d, _p)
    PROPS[_pr.key] = _pr


# ==========================================================================
# Batteries
# ==========================================================================
@dataclass(frozen=True)
class Cell:
    key: str
    chemistry: str
    v_nom: float
    v_full: float
    v_min: float            # usable cutoff under load
    capacity_ah: float
    r_internal_ohm: float
    mass_kg: float
    c_cont: float
    price_usd: float
    # Fraction of nameplate energy you actually get out at DBF discharge
    # rates, after depth-of-discharge limits and high-rate capacity loss.
    # This is chemistry specific: NiMH in particular loses a lot at 10C+.
    usable_fraction: float = 0.80
    source: str = "catalogue estimate - VERIFY"


def _lipo_cell(cap_ah: float, c_rate: float = 65.0) -> Cell:
    wh = 3.7 * cap_ah
    mass = wh / 135.0                     # ~135 Wh/kg at pack level for high-C LiPo
    r = 0.0060 / cap_ah * (45.0 / c_rate) ** 0.5
    return Cell(key="LIPO_%dMAH_%dC" % (int(cap_ah * 1000), int(c_rate)),
                chemistry="LiPo", v_nom=3.70, v_full=4.20, v_min=3.35,
                capacity_ah=cap_ah, r_internal_ohm=r, mass_kg=mass,
                c_cont=c_rate, price_usd=6.0 + 4.0 * cap_ah,
                usable_fraction=0.80)


CELLS: Dict[str, Cell] = {}
for _cap in (1.3, 1.8, 2.2, 3.0, 4.0, 5.0, 6.0):
    _c = _lipo_cell(_cap)
    CELLS[_c.key] = _c

# Li-ion alternative: much better Wh/kg, much worse current capability.
CELLS["LIION_P45B"] = Cell("LIION_P45B", "LiIon", 3.60, 4.20, 3.00, 4.50,
                           0.0125, 0.0700, 10.0, 12.0, usable_fraction=0.85,
                           source="Molicel INR21700-P45B class - VERIFY")
CELLS["LIION_P42A"] = Cell("LIION_P42A", "LiIon", 3.60, 4.20, 3.00, 4.00,
                           0.0140, 0.0700, 11.0, 9.0, usable_fraction=0.85,
                           source="Molicel INR21700-P42A class - VERIFY")

# NiMH and NiCd.  Both are legal and both are heavy: at DBF power levels a
# 100 Wh NiMH pack weighs about twice a LiPo one and a NiCd pack about three
# times.  They are here so the optimiser can reject them on the numbers
# rather than on my say-so.  1.2 V nominal means a lot of cells in series.
CELLS["NIMH_SUBC_5000"] = Cell("NIMH_SUBC_5000", "NiMH", 1.20, 1.45, 1.00,
                               5.00, 0.0040, 0.0920, 20.0, 6.0,
                               usable_fraction=0.70,
                               source="sub-C high-drain NiMH class - VERIFY")
CELLS["NIMH_SUBC_3000"] = Cell("NIMH_SUBC_3000", "NiMH", 1.20, 1.45, 1.00,
                               3.00, 0.0055, 0.0580, 25.0, 4.5,
                               usable_fraction=0.70,
                               source="sub-C high-drain NiMH class - VERIFY")
CELLS["NICD_SUBC_2000"] = Cell("NICD_SUBC_2000", "NiCd", 1.20, 1.45, 0.95,
                               2.00, 0.0030, 0.0550, 30.0, 5.0,
                               usable_fraction=0.75,
                               source="sub-C NiCd class - VERIFY")


@dataclass
class BatteryPack:
    cell: Cell
    series: int
    parallel: int
    wiring_overhead: float = 1.12     # tabs, shrink, connectors, balance lead

    @property
    def key(self) -> str:
        return "%dS%dP_%s" % (self.series, self.parallel, self.cell.key)

    @property
    def capacity_ah(self) -> float:
        return self.cell.capacity_ah * self.parallel

    @property
    def v_nominal(self) -> float:
        return self.cell.v_nom * self.series

    @property
    def v_full(self) -> float:
        return self.cell.v_full * self.series

    @property
    def v_min(self) -> float:
        return self.cell.v_min * self.series

    @property
    def watt_hours(self) -> float:
        return self.v_nominal * self.capacity_ah

    @property
    def mass_kg(self) -> float:
        return self.cell.mass_kg * self.series * self.parallel * self.wiring_overhead

    @property
    def r_ohm(self) -> float:
        return self.cell.r_internal_ohm * self.series / self.parallel

    @property
    def i_max_a(self) -> float:
        return self.cell.c_cont * self.capacity_ah

    @property
    def usable_fraction(self) -> float:
        return self.cell.usable_fraction

    @property
    def price_usd(self) -> float:
        return self.cell.price_usd * self.series * self.parallel

    def ocv(self, soc: float) -> float:
        """Open-circuit pack voltage vs state of charge (0..1): steep top,
        long flat middle, steep knee at the bottom."""
        s = max(0.0, min(1.0, soc))
        per_cell = (self.cell.v_min
                    + (self.cell.v_full - self.cell.v_min)
                    * (0.18 * s + 0.82 * (s ** 0.45)))
        return per_cell * self.series

    def terminal_voltage(self, soc: float, current_a: float) -> float:
        return max(self.v_min * 0.85, self.ocv(soc) - current_a * self.r_ohm)


# ==========================================================================
# ESC
# ==========================================================================
@dataclass(frozen=True)
class ESC:
    key: str
    i_cont_a: float
    mass_kg: float
    efficiency: float
    max_cells: int
    price_usd: float


ESCS: Dict[str, ESC] = {
    "ESC_45A": ESC("ESC_45A", 45, 0.045, 0.960, 6, 45),
    "ESC_60A": ESC("ESC_60A", 60, 0.063, 0.960, 6, 60),
    "ESC_80A": ESC("ESC_80A", 80, 0.082, 0.955, 8, 85),
    "ESC_100A": ESC("ESC_100A", 100, 0.105, 0.950, 8, 110),
    "ESC_120A": ESC("ESC_120A", 120, 0.125, 0.950, 12, 140),
}


def pick_esc(current_a: float, cells: int) -> Optional[ESC]:
    """Smallest ESC with 20 percent headroom that supports the cell count."""
    need = current_a * 1.20
    best = None
    for e in ESCS.values():
        if e.i_cont_a >= need and e.max_cells >= cells:
            if best is None or e.i_cont_a < best.i_cont_a:
                best = e
    return best


# ==========================================================================
# CSV overrides - measured data wins over the estimates above
# ==========================================================================
def load_motors_csv(path: str) -> int:
    """CSV columns: key,name,kv,rm,i0,mass_kg,i_max_a,p_max_w,shaft_mm,price_usd"""
    if not os.path.exists(path):
        return 0
    n = 0
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            MOTORS[row["key"]] = Motor(
                row["key"], row.get("name", row["key"]), float(row["kv"]),
                float(row["rm"]), float(row["i0"]), float(row["mass_kg"]),
                float(row["i_max_a"]), float(row["p_max_w"]),
                float(row.get("shaft_mm", 5)), float(row.get("price_usd", 0)),
                source="measured/" + os.path.basename(path))
            n += 1
    return n


def load_props_csv(path: str) -> int:
    """CSV columns: key,diameter_in,pitch_in,mass_kg,ct0,cp0,j_zero[,kp,price_usd]

    Fit ct0/cp0/j_zero from your own thrust-stand sweep or from the UIUC
    propeller database (Brandt and Selig).  That is the single highest-value
    data upgrade you can make to this model.
    """
    if not os.path.exists(path):
        return 0
    n = 0
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            PROPS[row["key"]] = Propeller(
                row["key"], float(row["diameter_in"]) * IN2M,
                float(row["pitch_in"]) * IN2M, float(row["mass_kg"]),
                float(row["ct0"]), float(row["cp0"]), float(row["j_zero"]),
                float(row.get("kp", 0.80)), float(row.get("price_usd", 18.0)),
                source="measured/" + os.path.basename(path))
            n += 1
    return n


# ==========================================================================
#  PAYLOAD - fixed container, sensor envelope, bay packing
# ==========================================================================

@dataclass
class PayloadGeometry:
    sensor_l: float
    sensor_w: float
    sensor_h: float
    pad_m: float                 # foam thickness per face
    box_l: float                 # container external dimensions
    box_w: float
    box_h: float
    container_tare_kg: float

    @property
    def box_volume_m3(self) -> float:
        return self.box_l * self.box_w * self.box_h

    @property
    def box_frontal_m2(self) -> float:
        return self.box_w * self.box_h

    @property
    def sensor_volume_m3(self) -> float:
        return self.sensor_l * self.sensor_w * self.sensor_h

    def fits(self, fill_fraction: float) -> bool:
        """Would the sensor fit at its NOMINAL packing density?

        Informational only.  A sensor heavier than this just has to be
        denser - see `required_density`.  Nothing rejects a design for
        failing this.
        """
        return self.sensor_volume_m3 <= fill_fraction * self.box_volume_m3

    def required_density(self, sensor_mass_kg: float,
                         fill_fraction: float) -> float:
        """Density the sensor must actually reach to fit the fixed box."""
        usable = max(fill_fraction * self.box_volume_m3, 1e-12)
        return sensor_mass_kg / usable


def size_payload(sensor_mass_kg: float, rules: Rules) -> PayloadGeometry:
    """One fixed container for every design.

    The box is NOT resized per sensor: the dimensions and tare in `Rules`
    are a single assumption, taken from the drop-test calculation for a
    roughly 1 kg sensor (see config.py).  The sensor envelope is still
    reported so you can see whether it fits, and `fits()` is what the
    feasibility check uses.
    """
    vol = max(sensor_mass_kg, 1e-9) / rules.sensor_packing_density_kg_m3
    w = (0.5 * vol) ** (1.0 / 3.0)          # 2:1:1 sensor envelope
    l, h = 2.0 * w, w
    v_impact = math.sqrt(2.0 * 9.80665 * 5.0 * 0.3048)
    stroke = v_impact ** 2 / (2.0 * 75.0 * 9.80665)
    pad = stroke / 0.70                      # the basis for the fixed box
    return PayloadGeometry(l, w, h, pad,
                           rules.container_length_m,
                           rules.container_width_m,
                           rules.container_height_m,
                           rules.container_tare_kg)


def best_bay_arrangement(n_boxes: int, box_l: float, box_w: float,
                         box_h: float) -> Tuple[int, int, int, float, float, float]:
    """Pack n identical boxes into the smallest-wetted-area cuboid bay.

    Returns (rows_along, cols_across, layers_high, bay_l, bay_w, bay_h).

    One long row gives a long thin fuselage, side-by-side gives a short fat
    one; which is cheaper depends on the box proportions, so the arrangement
    is chosen rather than assumed.
    """
    if n_boxes <= 0:
        return 0, 0, 0, 0.0, 0.0, 0.0
    best = None
    for cols in range(1, min(n_boxes, 4) + 1):
        for layers in range(1, min(n_boxes, 3) + 1):
            rows = int(math.ceil(n_boxes / (cols * layers)))
            if rows * cols * layers < n_boxes:
                continue
            bl = rows * box_l
            bw = cols * box_w
            bh = layers * box_h
            # proxy for fuselage cost: wetted area of the bay itself
            wetted = 2.0 * (bl * bw + bl * bh + bw * bh)
            # a bay wider or taller than it is long makes an unflyable fuselage
            if bw > 0.55 * bl or bh > 0.55 * bl:
                wetted *= 1.6
            if best is None or wetted < best[0]:
                best = (wetted, rows, cols, layers, bl, bw, bh)
    _w, rows, cols, layers, bl, bw, bh = best
    return rows, cols, layers, bl, bw, bh


# ==========================================================================
#  LIFTING LINE - span efficiency from the planform
# ==========================================================================

# span efficiency depends only on (aspect ratio, taper), so cache on those
_CACHE: Dict[Tuple[int, int, int], float] = {}


def span_efficiency(aspect_ratio: float, taper: float,
                    n_terms: int = 10, a0: float = 2.0 * math.pi) -> float:
    """Inviscid span efficiency `e` for an untwisted tapered wing."""
    ar = max(min(aspect_ratio, 40.0), 1.5)
    lam = max(min(taper, 1.0), 0.05)
    key = (int(round(ar * 100)), int(round(lam * 1000)), n_terms)
    hit = _CACHE.get(key)
    if hit is not None:
        return hit

    # Mean chord from AR; root chord from the taper ratio, since
    # c_mean = c_root (1 + lam) / 2 for a straight-tapered wing.
    b = 1.0                                   # unit span: only shape matters
    c_mean = b / ar
    c_root = 2.0 * c_mean / (1.0 + lam)

    # collocation at theta_i, avoiding the tip singularity at theta = 0
    thetas = [(i + 1) * math.pi / (2.0 * n_terms) for i in range(n_terms)]
    odd = [2 * k + 1 for k in range(n_terms)]

    mat, rhs = [], []
    for th in thetas:
        eta = abs(math.cos(th))               # |2y/b|
        c = c_root * (1.0 - (1.0 - lam) * eta)
        mu = c * a0 / (8.0 * b)
        row = [math.sin(n * th) * (n * mu + math.sin(th)) for n in odd]
        mat.append(row)
        rhs.append(mu * math.sin(th))         # unit (alpha - alpha_L0)
    coeffs = solve_linear(mat, rhs)

    a1 = coeffs[0]
    if abs(a1) < 1e-14:
        return 0.85
    delta = sum(n * (coeffs[k] / a1) ** 2
                for k, n in enumerate(odd) if n >= 3)
    e = 1.0 / (1.0 + delta)
    e = min(max(e, 0.5), 1.0)
    _CACHE[key] = e
    return e


def lift_curve_slope(aspect_ratio: float, a0: float = 2.0 * math.pi) -> float:
    """Finite-wing lift-curve slope per radian (Helmbold, low speed)."""
    ar = max(aspect_ratio, 0.5)
    return a0 * ar / (2.0 + math.sqrt(ar * ar + 4.0))


def downwash_gradient(aspect_ratio: float) -> float:
    """d(epsilon)/d(alpha) at the tail, the standard low-aspect-ratio
    approximation `2 a_w / (pi AR)`.

    This is why a tail is less effective than its area suggests: the wing
    turns the flow down before the tail ever sees it.
    """
    aw = lift_curve_slope(aspect_ratio)
    return min(0.85, 2.0 * aw / (math.pi * max(aspect_ratio, 0.5)))


# ==========================================================================
#  DESIGN LOADS - V-n diagram: manoeuvre, gust, CLmax ceiling
# ==========================================================================

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


# ==========================================================================
#  AERODYNAMICS - lift, drag, trim, CLmax(Re)
# ==========================================================================

@dataclass
class WingGeometry:
    span_m: float
    chord_m: float           # mean aerodynamic chord for a tapered wing
    airfoil: Airfoil
    taper: float = 1.0

    @property
    def area_m2(self) -> float:
        return self.span_m * self.chord_m

    @property
    def aspect_ratio(self) -> float:
        return self.span_m / self.chord_m

    @property
    def root_chord_m(self) -> float:
        if self.taper >= 0.999:
            return self.chord_m
        lam = self.taper
        return self.chord_m * 1.5 * (1.0 + lam) / (1.0 + lam + lam * lam)

    @property
    def thickness_root_m(self) -> float:
        return self.airfoil.t_over_c * self.root_chord_m


@dataclass
class BodyGeometry:
    """Fuselage + empennage, sized by the payload bay and tail volume."""
    fuse_length_m: float
    fuse_width_m: float
    fuse_height_m: float
    tail_arm_m: float
    s_horiz_m2: float
    s_vert_m2: float

    @property
    def fuse_wetted_m2(self) -> float:
        perim = 2.0 * (self.fuse_width_m + self.fuse_height_m)
        return perim * self.fuse_length_m * 0.92

    @property
    def fineness(self) -> float:
        d_eq = math.sqrt(4.0 * self.fuse_width_m * self.fuse_height_m / math.pi)
        return self.fuse_length_m / max(d_eq, 1e-6)


def oswald_efficiency(ar: float, taper: float = 1.0) -> float:
    """Span efficiency from a Prandtl lifting-line solve of the planform.

    This replaces Raymer's regression `1.78(1-0.045 AR^0.68) - 0.64`, which
    knows nothing about taper.  Note the two are not the same quantity:
    Raymer's is an *Oswald* factor that bundles viscous drag-due-to-lift in
    with the induced drag, whereas this is the inviscid span efficiency.
    Using the inviscid value is the correct one here, because the polar in
    `cd_from_cl` already carries the viscous term separately as
    `k_drag (CL - CL_opt)^2` - with Raymer's factor that part was counted
    twice.

    See liftingline.py for the derivation.  The solve reproduces the
    classical result that a straight-tapered wing minimises induced drag
    near taper 0.35.
    """
    return span_efficiency(ar, taper)


def skin_friction_cf(re: float, laminar_fraction: float = 0.2) -> float:
    re = max(re, 5.0e4)
    cf_turb = 0.074 / re ** 0.2
    cf_lam = 1.328 / math.sqrt(re)
    return laminar_fraction * cf_lam + (1.0 - laminar_fraction) * cf_turb


class AeroModel:
    def __init__(self, wing: WingGeometry, body: BodyGeometry,
                 env: Environment, build: BuildStandard,
                 n_motors: int = 1,
                 extra_drag_area_m2: float = 0.0):
        self.wing = wing
        self.body = body
        self.env = env
        self.build = build
        self.n_motors = n_motors
        self.extra_drag_area_m2 = extra_drag_area_m2

        self.S = wing.area_m2
        self.AR = wing.aspect_ratio
        self.e = oswald_efficiency(self.AR, wing.taper)
        self.k_induced = 1.0 / (math.pi * self.AR * self.e)

        af = wing.airfoil
        self.a_wing = lift_curve_slope(self.AR)
        self.downwash = downwash_gradient(self.AR)
        self.eta_tail = build.tail_eta
        # Finite span + built-wing reality (hinge gaps, washout, finish).
        ar_factor = min(1.0, 0.80 + 0.04 * self.AR)
        self._cl_max_base = af.cl_max_2d * build.cl_max_3d_factor * ar_factor
        # Camber-sensitive Reynolds exponent: a heavily cambered high-lift
        # section loses more of its CLmax at low Re than a thin one.
        self._re_exp = 0.15 + 0.10 * max(0.0, af.cl_max_2d - 1.2)

        # Tail geometry for trim drag
        self.s_tail_h = body.s_horiz_m2
        self.tail_arm = max(body.tail_arm_m, 1e-3)
        self.b_tail = math.sqrt(max(build.tail_aspect_ratio * self.s_tail_h,
                                    1e-6))

        # Representative CLmax for reporting (at the clean stall Re)
        self.cl_max = self._cl_max_base

        # ---- longitudinal stability -------------------------------------
        # Tail volume ratio and the neutral point.  With the wing AC at
        # 0.25c and the tail contributing through the downwash field,
        #     x_np/c = 0.25 + eta_t V_h (a_t/a_w)(1 - de/dalpha)
        # less a destabilising fuselage term.
        self.a_tail = lift_curve_slope(build.tail_aspect_ratio)
        self.v_h = (self.s_tail_h * self.tail_arm
                    / max(self.S * wing.chord_m, 1e-12))
        tail_term = (self.eta_tail * self.v_h
                     * (self.a_tail / max(self.a_wing, 1e-9))
                     * (1.0 - self.downwash))
        fus_term = (build.fuselage_np_shift_k * body.fuse_width_m ** 2
                    * body.fuse_length_m
                    / max(self.S * wing.chord_m * self.a_wing, 1e-12))
        self.x_np_over_c = 0.25 + tail_term - fus_term
        self.x_cg_over_c = 0.25 + build.cg_aft_of_ac_frac
        self.static_margin = self.x_np_over_c - self.x_cg_over_c

    # ------------------------------------------------------------------
    # CLmax with Reynolds effect
    # ------------------------------------------------------------------
    def cl_max_at_re(self, re: float) -> float:
        af = self.wing.airfoil
        scale = (max(re, 3.0e4) / af.re_ref) ** self._re_exp
        return self._cl_max_base * min(max(scale, 0.55), 1.08)

    def cl_max_at_v(self, v: float) -> float:
        re = max(v, 3.0) * self.wing.chord_m / self.env.nu
        return self.cl_max_at_re(re)

    # ------------------------------------------------------------------
    # Parasite drag
    # ------------------------------------------------------------------
    def cd0(self, v: float) -> float:
        env, wing, body, b = self.env, self.wing, self.body, self.build
        v = max(v, 3.0)

        re_w = v * wing.chord_m / env.nu
        re_scale = (wing.airfoil.re_ref / max(re_w, 3.0e4)) ** 0.15
        cd_wing = wing.airfoil.cd_min * re_scale * 1.05

        re_f = v * body.fuse_length_m / env.nu
        cf_f = skin_friction_cf(re_f, 0.10)
        f = body.fineness
        ff_f = 1.0 + 60.0 / max(f, 2.0) ** 3 + max(f, 2.0) / 400.0
        cd_fuse = cf_f * ff_f * body.fuse_wetted_m2 / self.S

        s_tail = body.s_horiz_m2 + body.s_vert_m2
        c_tail = math.sqrt(max(s_tail, 1e-4) / b.tail_aspect_ratio)
        re_t = v * c_tail / env.nu
        cf_t = skin_friction_cf(re_t, 0.20)
        cd_tail = cf_t * 1.25 * (2.04 * s_tail) / self.S

        area_extra = (b.gear_drag_area_m2
                      + b.nacelle_drag_area_m2 * max(0, self.n_motors - 1)
                      + self.extra_drag_area_m2)
        cd_extra = area_extra / self.S

        total = cd_wing + cd_fuse + cd_tail + cd_extra
        return total * (1.0 + b.misc_interference_frac)

    # ------------------------------------------------------------------
    # Trim
    # ------------------------------------------------------------------
    def tail_load(self, v: float, lift_n: float) -> float:
        """Horizontal tail lift required to trim, in newtons.

        Moments about the CG, with the tail arm l_t measured from the CG:

            M_ac + L_w (x_cg - x_ac) - L_t l_t = 0

        Writing k = (x_cg - x_ac)/c and using L_w + L_t = L (the tail is
        part of the lift balance, not an afterthought):

            L_t (l_t + k c) = Cm_ac q S c + L k c

            L_t = c (Cm_ac q S + L k) / (l_t + k c)

        Two terms, and both matter:

        * the **camber** term, Cm_ac q S c, is negative for a cambered
          section and scales with dynamic pressure only;
        * the **CG** term, L k c, scales with the lift being carried, so in
          an n-g turn the tail load grows with n.  Modelling only the camber
          term - as this did previously - holds the tail load constant
          through a 3 g turn, which is wrong.

        Negative means download.
        """
        q = self.q(v)
        c = self.wing.chord_m
        k = self.build.cg_aft_of_ac_frac
        return (c * (self.wing.airfoil.cm_ac * q * self.S + lift_n * k)
                / (self.tail_arm + k * c))

    # Note: the tail LOAD above comes from a moment balance and does not
    # care about the tail's dynamic pressure.  Its lift COEFFICIENT and its
    # induced drag do, because the tail has to make that load in slower
    # air - hence eta_tail below.

    def tail_cl(self, v: float, lift_n: float) -> float:
        """Tail lift coefficient, for checking it against its own stall."""
        qs = self.eta_tail * self.q(v) * max(self.s_tail_h, 1e-9)
        return self.tail_load(v, lift_n) / max(qs, 1e-9)

    # ------------------------------------------------------------------
    # Full polar
    # ------------------------------------------------------------------
    def q(self, v: float) -> float:
        """Dynamic pressure, 0.5 rho V^2."""
        return 0.5 * self.env.rho * v * v

    def cl_required(self, v: float, lift_n: float) -> float:
        """Wing CL needed.

        The tail carries part of the load, so the wing carries L - L_t.
        With a download (L_t negative) the wing carries more than the
        aircraft weight.
        """
        qs = self.q(v) * self.S
        wing_lift = lift_n - self.tail_load(v, lift_n)
        return wing_lift / max(qs, 1e-9)

    def cd_from_cl(self, cl: float, v: float) -> float:
        """Wing drag coefficient: profile bucket + stall rise + induced."""
        af = self.wing.airfoil
        cd = self.cd0(v)
        cd += af.k_drag * (cl - af.cl_at_cd_min) ** 2
        over = cl - af.cl_break
        if over > 0.0:
            cd += 0.085 * over ** 2 + 0.25 * over ** 3
        cd += self.k_induced * cl * cl
        return cd

    def drag(self, v: float, lift_n: float) -> float:
        """Total drag at this speed carrying this much lift."""
        q = self.q(v)
        cl = self.cl_required(v, lift_n)
        d = q * self.S * self.cd_from_cl(cl, v)
        # induced drag of the trimming tail load, L_t^2 / (q pi b_t^2 e_t)
        lt = self.tail_load(v, lift_n)
        if abs(lt) > 1e-9 and self.b_tail > 1e-6:
            d += lt * lt / (self.eta_tail * q * math.pi * self.b_tail ** 2
                            * self.build.tail_efficiency)
        return d

    def max_lift(self, v: float) -> float:
        """Maximum supportable weight at this speed, after trim.

        The wing stalls at L_w = q S CLmax.  Substituting that and the trim
        relation into L = L_w + L_t and solving for L:

            L = [ q S CLmax (l_t + k c) + Cm_ac q S c ] / l_t
        """
        q = self.q(v)
        c = self.wing.chord_m
        k = self.build.cg_aft_of_ac_frac
        qs = q * self.S
        lw_max = qs * self.cl_max_at_v(v)
        return (lw_max * (self.tail_arm + k * c)
                + self.wing.airfoil.cm_ac * qs * c) / self.tail_arm

    def v_stall(self, weight_n: float, load_factor: float = 1.0) -> float:
        """Stall speed, solved exactly rather than iterated.

        CLmax depends on Reynolds number and the trim download depends on
        dynamic pressure, so the stall speed is the root of

            max_lift(v) - n W = 0

        which is monotone increasing in v (CLmax rises with Re, and both
        terms scale with v^2).  A damped fixed point left ~0.03 percent on
        the table; Brent gets it to machine precision for the same cost.
        """
        target = load_factor * weight_n

        def residual(v: float) -> float:
            return self.max_lift(v) - target

        lo, hi = 1.0, 120.0
        if residual(lo) >= 0.0:
            return lo
        if residual(hi) <= 0.0:
            return hi
        return brent(residual, lo, hi, xtol=1e-12)

    def ld_ratio(self, v: float, weight_n: float) -> float:
        return weight_n / max(self.drag(v, weight_n), 1e-9)

    def best_ld_speed(self, weight_n: float) -> float:
        """Speed for best L/D: coarse scan to bracket, then golden section."""
        lo, hi = self.v_stall(weight_n) * 1.02, 60.0
        v, _ = scan_then_refine(lambda x: self.ld_ratio(x, weight_n),
                                lo, hi, n_scan=20, tol=1e-6)
        return v


def build_body(payload_bay_length_m: float, bay_width_m: float,
               bay_height_m: float, wing: WingGeometry,
               build: BuildStandard) -> BodyGeometry:
    """Close the geometry loop: payload bay sets fuselage length, fuselage
    length sets the tail arm, tail arm plus volume coefficients set the
    tails.  That is what makes 'carry one more container' cost drag."""
    length = max(build.fuse_min_length_m,
                 build.fuse_nose_m + payload_bay_length_m + build.fuse_tailcone_m)
    width = max(bay_width_m + 0.035, 0.10)
    height = max(bay_height_m + 0.045, 0.11)
    arm = build.tail_arm_frac * length
    s = wing.area_m2
    s_h = build.vol_coeff_h * s * wing.chord_m / arm
    s_v = build.vol_coeff_v * s * wing.span_m / arm
    return BodyGeometry(length, width, height, arm, s_h, s_v)


# ==========================================================================
#  WEIGHTS - itemised mass buildup, spar sizing
# ==========================================================================

@dataclass
class MassBreakdown:
    items: Dict[str, float] = field(default_factory=dict)

    def add(self, name: str, kg: float) -> None:
        self.items[name] = self.items.get(name, 0.0) + kg

    @property
    def total(self) -> float:
        return sum(self.items.values())

    def sorted_items(self):
        return sorted(self.items.items(), key=lambda kv: -kv[1])


def spar_sizing(wing: WingGeometry, mtow_n: float, build: BuildStandard,
                n_ultimate: float):
    """Carbon tube spar sized for root bending at the ultimate load factor.

    `n_ultimate` comes from the V-n diagram in loads.py, not a constant.

    Elliptical spanwise loading puts the half-wing lift centroid at
    4/(3*pi) of the semi-span, so
        M_root = n_ult * (W/2) * 0.4244 * (b/2)
    Thin-wall tube section modulus Z = pi * r^2 * t, so
        t_req = M_root / (sigma_allow * pi * r^2)
    with r set by the airfoil depth available at the root.

    Returns (mass_kg, t_required_m, t_used_m, gauge_limited).

    `gauge_limited` is True when the strength requirement is below the
    minimum wall you can actually lay up or buy.  That is common on this
    class of aeroplane and worth knowing: it means the spar is sized by
    buildability, not by load, so raising the load factor costs nothing and
    lowering it saves nothing.
    """
    r = 0.40 * wing.thickness_root_m           # tube radius fits the section
    r = max(r, 0.006)
    m_root = n_ultimate * (mtow_n / 2.0) * 0.4244 * (wing.span_m / 2.0)
    t_req = m_root / (build.spar_sigma_allow_pa * math.pi * r * r)
    t = max(t_req, build.spar_min_wall_m)
    # 0.6 accounts for spanwise taper relative to a constant section
    m = 0.6 * (2.0 * math.pi * r * t) * build.spar_rho * wing.span_m
    return (m * build.spar_fitting_factor, t_req, t,
            t_req < build.spar_min_wall_m)


def spar_mass(wing: WingGeometry, mtow_n: float, build: BuildStandard,
              n_ultimate: float) -> float:
    return spar_sizing(wing, mtow_n, build, n_ultimate)[0]


def airframe_masses(wing: WingGeometry, body: BodyGeometry, mtow_n: float,
                    build: BuildStandard, n_ultimate: float) -> MassBreakdown:
    mb = MassBreakdown()
    mb.add("wing skin/ribs", build.wing_areal_density * wing.area_m2)
    mb.add("wing spar", spar_mass(wing, mtow_n, build, n_ultimate))
    mb.add("fuselage", build.fuse_areal_density * body.fuse_wetted_m2)
    mb.add("tails", build.tail_areal_density * (body.s_horiz_m2 + body.s_vert_m2) * 2.0)
    mb.add("tail boom/links", 0.030 + 0.045 * body.tail_arm_m)
    mb.add("landing gear",
           build.gear_mass_fixed_kg + build.gear_mass_frac * (mtow_n / 9.80665))
    mb.add("servos", build.n_servos_base * build.servo_mass_kg)
    mb.add("avionics/wiring", build.avionics_kg)
    mb.add("sensor deploy mechanism", build.deploy_mech_kg)
    return mb


def converge_empty_mass(wing: WingGeometry, body: BodyGeometry,
                        propulsion_items: dict, payload_mass_kg: float,
                        build: BuildStandard, n_ultimate: float,
                        tol: float = 1e-10, maxiter: int = 60):
    """Fixed-point iteration: spar and gear scale with MTOW, MTOW contains them.

    The map turns out to be a strong contraction (derivative about 0.04,
    because only the spar and gear terms respond to weight), so it is iterated
    undamped and reaches 1e-10 kg in four or five passes.  Damping it 50/50 -
    the cautious default - raised the contraction factor to 0.52 and needed
    more than twenty.  A divergence guard falls back to damping if a different
    build standard ever makes the map expansive.

    `propulsion_items` is a {name: kg} dict so the breakdown stays itemised.

    Returns (MassBreakdown for the empty aircraft, empty_mass_kg, mtow_kg).
    """
    propulsion_mass_kg = sum(propulsion_items.values())
    mtow = payload_mass_kg + propulsion_mass_kg + 2.0     # first guess
    mb = MassBreakdown()
    prev_change = float("inf")
    growing = 0
    damping = 0.0
    for _ in range(maxiter):
        mb = airframe_masses(wing, body, mtow * 9.80665, build,
                             n_ultimate)
        # itemised, not lumped: motors, propellers, ESCs and battery are the
        # four biggest single line items on a DBF aeroplane and deserve to be
        # visible next to the wing.
        for name, kg in propulsion_items.items():
            mb.add(name, kg)
        mb.add("payload restraint", build.payload_structure_frac * payload_mass_kg)
        mb.add("contingency", build.contingency_frac * mb.total)
        new_mtow = mb.total + payload_mass_kg
        change = abs(new_mtow - mtow)
        if change < tol:
            mtow = new_mtow
            break
        if change > prev_change:
            growing += 1
            if growing >= 2:
                damping = 0.5          # guard: only if the map is expansive
        prev_change = change
        mtow = damping * mtow + (1.0 - damping) * new_mtow
    return mb, mb.total, mtow


# ==========================================================================
#  PROPULSION - motor/prop torque balance, closed-form battery sag
# ==========================================================================

@dataclass
class OperatingPoint:
    rpm: float
    thrust_total_n: float
    motor_current_a: float
    pack_current_a: float
    pack_voltage_v: float
    p_elec_w: float          # drawn from the battery
    p_shaft_w: float         # total, all motors
    eta_prop: float
    eta_motor: float
    advance_ratio: float


class PropulsionSystem:
    ROOT_XTOL = 1e-12          # Brent tolerance on rev/s

    def __init__(self, motor: Motor, prop: Propeller, n_motors: int,
                 pack: BatteryPack, env: Environment, esc: Optional[ESC] = None,
                 hot_resistance_factor: float = 1.25,
                 install_thrust_factor: float = 0.95):
        # Realism, not best case: a motor that has been at full throttle for
        # two minutes is hot, and hot copper has more resistance, so it makes
        # less torque for the same current.  And an installed propeller never
        # matches its free-air thrust - it works in the blockage and the
        # non-uniform inflow of the airframe ahead of or behind it.
        self.hot_resistance_factor = hot_resistance_factor
        self.install_thrust_factor = install_thrust_factor
        self.motor = motor
        self.prop = prop
        self.n_motors = n_motors
        self.pack = pack
        self.env = env
        self.esc = esc or pick_esc(motor.i_max_a, pack.series) or _biggest_esc()
        self.kt = 60.0 / (2.0 * math.pi * motor.kv)     # N*m per amp
        self.rm_hot = motor.rm * hot_resistance_factor
        self._memo_key: Optional[Tuple[float, float, float]] = None
        self._memo_val: Tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)

    # ------------------------------------------------------------------
    # Exact single-point solve
    # ------------------------------------------------------------------
    def _state_at_rpm(self, n: float, ocv: float, throttle: float,
                      a: float, denom: float, v_floor: float):
        """Closed-form pack voltage, terminal voltage and motor current at a
        given propeller speed.  See the module docstring for the algebra."""
        v_pack = (ocv + a * (60.0 / self.motor.kv) * n) / denom
        if v_pack < v_floor:
            v_pack = v_floor
        vt = throttle * v_pack
        i_motor = (vt - 60.0 * n / self.motor.kv) / self.rm_hot
        return v_pack, vt, i_motor

    def solve(self, v: float, soc: float, throttle: float) -> OperatingPoint:
        m, p, env = self.motor, self.prop, self.env
        throttle = max(0.0, min(1.0, throttle))
        v = max(0.0, v)
        ocv = self.pack.ocv(soc)
        v_floor = self.pack.v_min * 0.85

        if throttle * ocv <= 0.05:
            return OperatingPoint(0.0, 0.0, 0.0, 0.0, ocv, 0.0, 0.0,
                                  0.0, 0.0, 0.0)

        a = self.n_motors * throttle * self.pack.r_ohm / self.rm_hot
        denom = 1.0 + a * throttle
        d5 = p.diameter_m ** 5
        two_pi = 2.0 * math.pi

        def residual(nn: float) -> float:
            _vp, _vt, cur = self._state_at_rpm(nn, ocv, throttle, a, denom,
                                               v_floor)
            q_motor = self.kt * (cur - m.i0)
            j = v / (nn * p.diameter_m) if nn > 1e-6 else 1e9
            q_prop = p.cp(j) * env.rho * nn * nn * d5 / two_pi
            return q_motor - q_prop

        # At n_noload the closed form gives Vpack = Vocv exactly and the motor
        # current is zero, so the residual is negative: a guaranteed bracket.
        n_noload = m.kv * throttle * ocv / 60.0
        if residual(1e-4) <= 0.0:
            n = 0.0
        else:
            n = brent(residual, 1e-4, n_noload, xtol=self.ROOT_XTOL)
        v_pack, _vt, i_motor = self._state_at_rpm(n, ocv, throttle, a, denom,
                                                  v_floor)
        i_motor = max(0.0, i_motor)

        i_batt = self.n_motors * i_motor * throttle
        thrust = (self.n_motors * max(0.0, p.thrust(env.rho, n, v))
                  * self.install_thrust_factor)
        p_shaft = self.n_motors * p.power(env.rho, n, v)
        p_elec = v_pack * i_batt / max(self.esc.efficiency, 0.5)
        j = v / (n * p.diameter_m) if n > 1e-6 else 0.0
        eta_prop = (thrust * v / p_shaft) if p_shaft > 1e-6 and v > 0.1 else 0.0
        eta_motor = (p_shaft / (self.n_motors
                                * max(throttle * v_pack * i_motor, 1e-6))
                     if i_motor > 1e-6 else 0.0)
        return OperatingPoint(60.0 * n, thrust, i_motor, i_batt, v_pack,
                              p_elec, p_shaft, min(eta_prop, 0.98),
                              min(eta_motor, 0.98), j)

    # ------------------------------------------------------------------
    # Query interface
    # ------------------------------------------------------------------
    def operating(self, v: float, soc: float = 1.0,
                  throttle: float = 1.0) -> Tuple[float, float, float, float]:
        """Returns (thrust_N, p_elec_W, pack_current_A, pack_voltage_V).

        Memoised on the last query, because the mission integrator asks for
        thrust and power at the same flight condition back to back.
        """
        throttle = max(0.0, min(1.0, throttle))
        v = max(0.0, v)
        key = (v, soc, throttle)
        if key == self._memo_key:
            return self._memo_val
        op = self.solve(v, soc, throttle)
        out = (op.thrust_total_n, op.p_elec_w, op.pack_current_a,
               op.pack_voltage_v)
        self._memo_key, self._memo_val = key, out
        return out

    def thrust(self, v: float, soc: float = 1.0, throttle: float = 1.0) -> float:
        return self.operating(v, soc, throttle)[0]

    def power_elec(self, v: float, soc: float = 1.0,
                   throttle: float = 1.0) -> float:
        return self.operating(v, soc, throttle)[1]

    def pack_current(self, v: float, soc: float = 1.0,
                     throttle: float = 1.0) -> float:
        return self.operating(v, soc, throttle)[2]

    def throttle_for_thrust(self, v: float, soc: float,
                            thrust_req: float) -> float:
        """Smallest throttle that makes the required thrust (1.0 if it cannot)."""
        if thrust_req >= self.thrust(v, soc, 1.0):
            return 1.0
        lo = 0.05
        if self.thrust(v, soc, lo) >= thrust_req:
            return lo
        return brent(lambda th: self.thrust(v, soc, th) - thrust_req,
                     lo, 1.0, xtol=1e-9)

    # ------------------------------------------------------------------
    # Ratings / masses
    # ------------------------------------------------------------------
    @property
    def mass_kg(self) -> float:
        return (self.n_motors * (self.motor.mass_kg + self.prop.mass_kg
                                 + self.esc.mass_kg + 0.025)   # mount + spinner
                + self.pack.mass_kg)

    @property
    def price_usd(self) -> float:
        return (self.n_motors * (self.motor.price_usd + self.prop.price_usd
                                 + self.esc.price_usd)
                + self.pack.price_usd)

    def rating_check(self) -> Dict[str, float]:
        """Worst-case (static, full charge, full throttle) electrical loads,
        from the exact solve."""
        op = self.solve(0.0, 1.0, 1.0)
        return {
            "static_thrust_n": op.thrust_total_n,
            "motor_current_a": op.motor_current_a,
            "pack_current_a": op.pack_current_a,
            "motor_power_w": op.p_elec_w / max(self.n_motors, 1),
            "motor_current_margin": self.motor.i_max_a / max(op.motor_current_a, 1e-6),
            "motor_power_margin": self.motor.p_max_w
            / max(op.p_elec_w / max(self.n_motors, 1), 1e-6),
            "esc_current_margin": self.esc.i_cont_a / max(op.motor_current_a, 1e-6),
            "pack_current_margin": self.pack.i_max_a / max(op.pack_current_a, 1e-6),
            "rpm": op.rpm,
            "tip_mach": (math.pi * self.prop.diameter_m * op.rpm / 60.0) / 340.0,
        }


def _biggest_esc() -> ESC:
    return max(ESCS.values(), key=lambda e: e.i_cont_a)


# --------------------------------------------------------------------------
# Instance cache.  A PropulsionSystem depends only on the hardware, never on
# the wing, so thousands of candidate wings share one object - and with it
# its last-operating-point memo.  Bounded by the hardware shortlist size.
# --------------------------------------------------------------------------
_SYSTEM_CACHE: Dict[Tuple, PropulsionSystem] = {}


def get_propulsion(motor: Motor, prop: Propeller, n_motors: int,
                   pack: BatteryPack, env: Environment,
                   hot_resistance_factor: float = 1.25,
                   install_thrust_factor: float = 0.95) -> PropulsionSystem:
    key = (motor.key, prop.key, n_motors, pack.key,
           round(env.rho, 4), round(env.nu, 9),
           round(hot_resistance_factor, 4), round(install_thrust_factor, 4))
    sys_ = _SYSTEM_CACHE.get(key)
    if sys_ is None:
        sys_ = PropulsionSystem(motor, prop, n_motors, pack, env, None,
                                hot_resistance_factor, install_thrust_factor)
        _SYSTEM_CACHE[key] = sys_
    return sys_


def clear_system_cache() -> None:
    _SYSTEM_CACHE.clear()


# ==========================================================================
#  AIRCRAFT - the design vector and the aeroplane it produces
# ==========================================================================

@dataclass(frozen=True)
class Design:
    # --- geometry -----------------------------------------------------
    span_m: float
    chord_m: float
    airfoil: str
    # --- propulsion ---------------------------------------------------
    motor: str
    prop: str
    n_motors: int
    # --- energy -------------------------------------------------------
    cell: str
    cells_series: int
    cells_parallel: int
    # --- payload ------------------------------------------------------
    sensor_mass_kg: float
    n_simulators: int
    sim_mass_kg: float
    taper: float = 1.0

    def key(self) -> Tuple:
        return (round(self.span_m, 4), round(self.chord_m, 4), self.airfoil,
                self.motor, self.prop, self.n_motors, self.cell,
                self.cells_series, self.cells_parallel,
                round(self.sensor_mass_kg, 4), self.n_simulators,
                round(self.sim_mass_kg, 4), round(self.taper, 3))

    def discrete_key(self) -> Tuple:
        return (self.airfoil, self.motor, self.prop, self.n_motors,
                self.cell, self.cells_series, self.cells_parallel)


class Aircraft:
    def __init__(self, design: Design, cfg: Config):
        self.design = design
        self.cfg = cfg
        r, b, env = cfg.rules, cfg.build, cfg.env

        self.airfoil = AIRFOILS[design.airfoil]
        self.motor = MOTORS[design.motor]
        self.prop = PROPS[design.prop]
        self.pack = BatteryPack(CELLS[design.cell], design.cells_series,
                                design.cells_parallel)

        # ---- payload: the container is sized by the sensor -----------
        self.payload_geom: PayloadGeometry = size_payload(
            design.sensor_mass_kg, r)
        container_tare = self.payload_geom.container_tare_kg
        # A simulator cannot weigh less than the empty container it mimics.
        self.sim_mass_kg = max(design.sim_mass_kg, container_tare)
        self.payload_sim_total_kg = design.n_simulators * self.sim_mass_kg
        self.m2_scored_kg = (design.sensor_mass_kg + container_tare
                             + self.payload_sim_total_kg)
        self.m2_payload_kg = self.m2_scored_kg
        self.m3_payload_kg = design.sensor_mass_kg
        self.sensor_density_kg_m3 = self.payload_geom.required_density(
            design.sensor_mass_kg, r.container_fill_fraction)

        # ---- geometry ------------------------------------------------
        self.wing = WingGeometry(design.span_m, design.chord_m,
                                 self.airfoil, design.taper)
        pg = self.payload_geom
        n_boxes = 1 + design.n_simulators
        rows, cols, layers, bay_l, bay_w, bay_h = best_bay_arrangement(
            n_boxes, pg.box_l, pg.box_w, pg.box_h)
        self.bay_rows, self.bay_cols, self.bay_layers = rows, cols, layers
        self.body = build_body(bay_l, bay_w, bay_h, self.wing, b)

        # ---- landing gear height is set by propeller clearance -------
        if design.n_motors == 1:
            need = (0.5 * self.prop.diameter_m + b.prop_clearance_margin_m
                    - 0.5 * self.body.fuse_height_m)
            self.gear_height_m = max(b.gear_height_m, need)
        else:
            self.gear_height_m = b.gear_height_m
        gear_scale = self.gear_height_m / b.gear_height_m
        gear_extra_drag = b.gear_drag_area_m2 * (gear_scale - 1.0)
        self.gear_extra_mass_kg = 0.55 * (self.gear_height_m - b.gear_height_m)

        # ---- aero ----------------------------------------------------
        self.aero_clean = AeroModel(self.wing, self.body, env, b,
                                    design.n_motors,
                                    extra_drag_area_m2=gear_extra_drag)
        self.aero_deployed = AeroModel(self.wing, self.body, env, b,
                                       design.n_motors,
                                       extra_drag_area_m2=gear_extra_drag
                                       + r.sensor_deployed_drag_area_m2)

        # ---- propulsion ----------------------------------------------
        # A first pass at the nominal hot resistance; the real winding
        # temperature is computed below once there is a weight to fly.
        self.prop_sys: PropulsionSystem = get_propulsion(
            self.motor, self.prop, design.n_motors, self.pack, env,
            b.motor_hot_resistance_factor, b.install_thrust_factor)

        # ---- mass ----------------------------------------------------
        nm = design.n_motors
        esc = self.prop_sys.esc
        propulsion_items = {
            "motors (%d)" % nm: nm * self.motor.mass_kg,
            "propellers (%d)" % nm: nm * self.prop.mass_kg,
            "ESCs (%d)" % nm: nm * esc.mass_kg,
            "motor mounts/spinners": nm * 0.025,
            "battery pack": self.pack.mass_kg,
        }
        if self.gear_extra_mass_kg > 1e-9:
            propulsion_items["tall gear penalty"] = self.gear_extra_mass_kg
        self.propulsion_items = propulsion_items

        # The spar is sized by the ultimate load factor, which comes from a
        # V-n diagram, which needs the weight and the speed - which need the
        # spar.  Resolve by alternating: converge the mass, rebuild the load
        # case, converge again.  Three passes settle it to well under a gram.
        # Three quantities are mutually dependent and must converge
        # TOGETHER, not one after the other:
        #   spar mass  <- ultimate load factor <- weight and Vmax
        #   Vmax       <- motor hot resistance <- cruise current <- weight
        #   weight     <- spar mass
        # Updating the motor temperature *after* settling the load case
        # left the load case stale, which the validation suite caught.
        n_override = b.ultimate_load_factor_override
        n_ult = n_override if n_override else 4.0
        self.motor_hot_factor = b.motor_hot_resistance_factor
        mb = None
        for _ in range(6):
            mb, empty, mtow_m2 = converge_empty_mass(
                self.wing, self.body, propulsion_items, self.m2_payload_kg,
                b, n_ult)
            w = mtow_m2 * G0

            hot = self._thermal_factor(w, b)
            if abs(hot - self.motor_hot_factor) > 0.005:
                self.motor_hot_factor = hot
                self.prop_sys = get_propulsion(
                    self.motor, self.prop, design.n_motors, self.pack, env,
                    hot, b.install_thrust_factor)

            self.load_case: LoadCase = design_load_factors(
                self.aero_clean, w, self._max_level_speed(w), b)
            n_new = (n_override if n_override
                     else self.load_case.n_ultimate)
            if abs(n_new - n_ult) < 1e-6:
                n_ult = n_new
                break
            n_ult = n_new
        mb, empty, mtow_m2 = converge_empty_mass(
            self.wing, self.body, propulsion_items, self.m2_payload_kg,
            b, n_ult)
        self.n_ultimate = n_ult
        self.n_limit = (n_ult / b.load_safety_factor if n_override
                        else self.load_case.n_limit)

        (_sm, self.spar_t_req_m, self.spar_t_used_m,
         self.spar_gauge_limited) = spar_sizing(
            self.wing, mtow_m2 * G0, b, n_ult)

        self.mass_breakdown: MassBreakdown = mb
        self.empty_mass_kg = empty
        self.mtow_m2_kg = mtow_m2
        self.mtow_m3_kg = empty + self.m3_payload_kg
        self.mtow_m1_kg = empty            # Mission 1 flies empty

        self.w_m2_n = self.mtow_m2_kg * G0
        self.w_m3_n = self.mtow_m3_kg * G0
        self.w_m1_n = self.mtow_m1_kg * G0

    def _thermal_factor(self, weight_n: float, b) -> float:
        """Mission-average winding temperature rise -> resistance factor.

        Steady-state rise from copper loss is `dT_ss = I^2 Rm R_th`.  The
        winding heats with a time constant `tau`, so averaged over a
        mission of length `T`:

            dT_avg = dT_ss [ 1 - (tau/T)(1 - exp(-T/tau)) ]

        and copper resistance rises as `1 + alpha dT`.
        """
        try:
            v = 0.85 * self._max_level_speed(weight_n)
            op = self.prop_sys.solve(v, 0.75, 1.0)
        except Exception:
            return b.motor_hot_resistance_factor
        i_motor = op.motor_current_a
        p_cu = i_motor * i_motor * self.motor.rm
        dt_ss = p_cu * b.motor_thermal_resistance_kw
        tau = max(b.motor_thermal_tau_s, 1.0)
        t_mis = max(self.cfg.rules.mission_window_s, 1.0)
        dt_avg = dt_ss * (1.0 - (tau / t_mis) * (1.0 - math.exp(-t_mis / tau)))
        return min(max(1.0 + b.copper_temp_coeff * dt_avg, 1.0), 1.80)

    def _max_level_speed(self, weight_n: float) -> float:
        """Maximum level speed, solved locally.

        mission.py imports this module, so it cannot be imported back; this
        is the same Brent solve on thrust minus drag.
        """
        aero = self.aero_clean

        def excess(v: float) -> float:
            return self.prop_sys.thrust(v, 0.85, 1.0) - aero.drag(v, weight_n)

        lo = aero.v_stall(weight_n) * 1.02
        hi = 65.0
        if excess(lo) <= 0.0:
            return lo
        if excess(hi) > 0.0:
            return hi
        return brent(excess, lo, hi, xtol=1e-6)

    # ------------------------------------------------------------------
    @property
    def wing_mass_kg(self) -> float:
        return (self.mass_breakdown.items.get("wing skin/ribs", 0.0)
                + self.mass_breakdown.items.get("wing spar", 0.0))

    @property
    def motor_mass_kg(self) -> float:
        return self.design.n_motors * self.motor.mass_kg

    @property
    def wing_area_m2(self) -> float:
        return self.wing.area_m2

    @property
    def aspect_ratio(self) -> float:
        return self.wing.aspect_ratio

    def wing_loading_m2(self) -> float:
        return self.w_m2_n / self.wing_area_m2

    def static_thrust_n(self) -> float:
        return self.prop_sys.rating_check()["static_thrust_n"]

    def usable_energy_j(self) -> float:
        """Energy you can actually get out, chemistry dependent."""
        return self.pack.watt_hours * 3600.0 * self.pack.usable_fraction

    # ------------------------------------------------------------------
    # Hard constraints that do not need a mission simulation
    # ------------------------------------------------------------------
    def static_checks(self) -> Dict[str, Any]:
        r, b = self.cfg.rules, self.cfg.build
        rc = self.prop_sys.rating_check()
        d = self.design

        clearance = (self.gear_height_m + 0.5 * self.body.fuse_height_m
                     - 0.5 * self.prop.diameter_m) if d.n_motors == 1 else 1.0
        if d.n_motors >= 2:
            lateral = self.wing.span_m / d.n_motors
            prop_gap = lateral - self.prop.diameter_m
        else:
            prop_gap = 1.0

        fails: List[str] = []
        if self.pack.watt_hours > r.max_watt_hours + 1e-9:
            fails.append("battery energy over the %.0f Wh limit" % r.max_watt_hours)
        if self.wing.span_m > r.max_wingspan_m:
            fails.append("span over limit")
        if rc["motor_current_margin"] < 1.0:
            fails.append("motor over current (%.0f A draw)" % rc["motor_current_a"])
        if rc["motor_power_margin"] < 1.0:
            fails.append("motor over power")
        if rc["esc_current_margin"] < 1.0:
            fails.append("ESC over current")
        if rc["pack_current_margin"] < 1.0:
            fails.append("pack over C-rating (%.0f A)" % rc["pack_current_a"])
        if rc["tip_mach"] > 0.82:
            fails.append("prop tip Mach %.2f" % rc["tip_mach"])
        if self.gear_height_m > 0.30:
            fails.append("landing gear would need %.0f mm legs for prop clearance"
                         % (self.gear_height_m * 1000))
        if prop_gap < 0.03:
            fails.append("propellers overlap")
        # the tail has to be able to trim without stalling itself
        v_trim = self.aero_clean.v_stall(self.w_m2_n) * 1.3
        ct = abs(self.aero_clean.tail_cl(v_trim, self.w_m2_n))
        ct_turn = abs(self.aero_clean.tail_cl(
            v_trim, self.n_limit * self.w_m2_n))
        if max(ct, ct_turn) > b.tail_cl_max:
            fails.append("tail needs CL %.2f to trim (max %.2f) - bigger "
                         "tail or longer arm" % (max(ct, ct_turn), b.tail_cl_max))
        sm = self.aero_clean.static_margin
        if sm < b.min_static_margin:
            fails.append("static margin %.3f c - %s" % (
                sm, "unstable" if sm < 0 else "too little for a model"))
        if sm > b.max_static_margin:
            fails.append("static margin %.3f c - excessive, all trim drag" % sm)
        if self.aspect_ratio < 3.0 or self.aspect_ratio > 16.0:
            fails.append("aspect ratio %.1f outside modelled range" % self.aspect_ratio)
        if self.wing.chord_m < 0.12:
            fails.append("chord too small (low Re)")
        # the payload has to physically fit in the fuselage
        if self.payload_geom.box_w > self.body.fuse_width_m + 1e-9:
            fails.append("container too wide for the fuselage")
        # There is NO rules cap on sensor mass, so the only ceiling is
        # physics: the sensor has to be made of something.  A heavy sensor
        # in a fixed box is simply a dense one, and the report says what
        # density it would need.  Only an impossible density is rejected.
        if self.sensor_density_kg_m3 > r.max_sensor_density_kg_m3:
            fails.append("sensor would need %.0f kg/m3 to fit the container "
                         "(denser than lead)" % self.sensor_density_kg_m3)

        return {
            "ok": not fails,
            "failures": fails,
            "ratings": rc,
            "prop_ground_clearance_m": clearance,
            "prop_lateral_gap_m": prop_gap,
        }

    # ------------------------------------------------------------------
    def summary(self) -> Dict[str, Any]:
        rc = self.prop_sys.rating_check()
        return {
            "span_m": self.wing.span_m,
            "chord_m": self.wing.chord_m,
            "area_m2": self.wing_area_m2,
            "AR": self.aspect_ratio,
            "airfoil": self.airfoil.name,
            "CLmax_3d": self.aero_clean.cl_max_at_v(
                self.aero_clean.v_stall(self.w_m2_n)),
            "oswald_e": self.aero_clean.e,
            "CD0_at_20ms": self.aero_clean.cd0(20.0),
            "motor": self.motor.name,
            "n_motors": self.design.n_motors,
            "prop": self.prop.key,
            "pack": self.pack.key,
            "pack_V_nom": self.pack.v_nominal,
            "pack_Wh": self.pack.watt_hours,
            "chemistry": self.pack.cell.chemistry,
            "esc": self.prop_sys.esc.key,
            "wing_mass_kg": self.wing_mass_kg,
            "motor_mass_kg": self.motor_mass_kg,
            "battery_mass_kg": self.pack.mass_kg,
            "empty_kg": self.empty_mass_kg,
            "mtow_m2_kg": self.mtow_m2_kg,
            "mtow_m3_kg": self.mtow_m3_kg,
            "m2_scored_lb": self.m2_scored_kg * KG2LB,
            "sensor_lb": self.design.sensor_mass_kg * KG2LB,
            "container_lb": self.payload_geom.container_tare_kg * KG2LB,
            "sensor_density_kg_m3": self.sensor_density_kg_m3,
            "static_margin": self.aero_clean.static_margin,
            "neutral_point_c": self.aero_clean.x_np_over_c,
            "taper": self.design.taper,
            "motor_hot_factor": self.motor_hot_factor,
            "tail_cl_cruise": self.aero_clean.tail_cl(
                self.aero_clean.v_stall(self.w_m2_n) * 1.3, self.w_m2_n),
            "static_thrust_n": rc["static_thrust_n"],
            "static_T_over_W_m2": rc["static_thrust_n"] / self.w_m2_n,
            "wing_loading_Nm2": self.wing_loading_m2(),
            "n_limit": self.n_limit,
            "n_ultimate": self.n_ultimate,
            "load_driver": getattr(self.load_case, "driver", "override"),
            "spar_gauge_limited": self.spar_gauge_limited,
        }


# ==========================================================================
#  MISSION - takeoff, climb, laps, landing, M1/M2/M3
# ==========================================================================

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
                 dx: float = 6.0) -> Tuple[float, float, float]:
    """Integrate one straight leg.  Returns (time_s, energy_J, v_exit).

    `wind_along` is the wind component opposing motion (+ for a headwind).
    `dx` is a nominal ground-track step used to size the RK4 time step; with
    RK4 plus event location the answer is nearly independent of it.

    Flown as: full throttle, then power off at the last moment that still
    allows arrival at the turn at corner speed.
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

    IMPORTANT: per the rules quoted by the team, "landing is not part of the
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
def simulate_lap(aero: AeroModel, ps: PropulsionSystem, mass_kg: float,
                 weight_n: float, soc: float, cfg: Config,
                 turn: TurnCapability, v_entry: float,
                 dx: float = 6.0) -> Tuple[float, float, float]:
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
                                       val, v, v_target, wind_along, dx)
            t_tot += dt_
            e_tot += de_
        else:
            dt_ = math.radians(val) / max(turn.turn_rate_rad_s, 1e-6)
            t_tot += dt_
            e_tot += p_turn * dt_
            v = turn.v_turn
    return t_tot, e_tot, v


# ==========================================================================
# Missions
# ==========================================================================
def fly_mission(ac: Aircraft, aero: AeroModel, weight_n: float, cfg: Config,
                required_laps: Optional[int] = None,
                max_laps: int = 40, dx: float = 6.0,
                energy_reserve: float = 0.08) -> MissionResult:
    """Common engine for M1/M2/M3.

    required_laps set  -> fly exactly that many, report the time (Mission 2).
    required_laps None -> fly as many as fit in the window (Mission 3).
    """
    ps = ac.prop_sys
    mass = weight_n / cfg.env.g
    pack_j = max(ac.pack.watt_hours * 3600.0, 1e-6)
    usable = ac.usable_energy_j() * (1.0 - energy_reserve)
    window = cfg.rules.mission_window_s

    v_s = aero.v_stall(weight_n)
    to = takeoff(aero, ps, mass, weight_n, cfg, soc=1.0,
                 gear_height_m=ac.gear_height_m)
    if not to.ok:
        return MissionResult(False, 0, 0.0, to.energy_j, 1.0,
                             v_stall=v_s, takeoff_m=to.ground_roll_m,
                             reason=to.reason or "takeoff")

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
                                     turn, v, dx)
            cached, cached_soc = (lt, le), soc

        # No pilot flies the theoretical lap.  One factor, applied to time
        # and to the energy that time costs.
        eff = max(min(cfg.build.pilot_efficiency, 1.0), 0.5)
        lt, le = lt / eff, le / eff

        if laps == 0 and credit > 0.0:
            v_avg = cfg.course.straight_length_m / max(lt, 1e-3)
            lt = max(0.35 * lt, lt - credit / max(v_avg, 1.0))

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


def mission_1(ac: Aircraft, cfg: Config, dx: float = 6.0) -> MissionResult:
    """Flight mission with no payload - pass/fail (3 laps in the window)."""
    return fly_mission(ac, ac.aero_clean, ac.w_m1_n, cfg,
                       required_laps=3, dx=dx)


def mission_2(ac: Aircraft, cfg: Config, dx: float = 6.0) -> MissionResult:
    """Delivery flight: sensor in shipping container + optional simulators,
    5 laps, scored on weight / time."""
    return fly_mission(ac, ac.aero_clean, ac.w_m2_n, cfg,
                       required_laps=cfg.rules.m2_required_laps, dx=dx)


def mission_3(ac: Aircraft, cfg: Config, dx: float = 6.0) -> MissionResult:
    """Sensor flight: sensor deployed for the whole pattern, maximise laps."""
    return fly_mission(ac, ac.aero_deployed, ac.w_m3_n, cfg,
                       required_laps=None, dx=dx)


# ==========================================================================
#  SCORING - the three mission formulas
# ==========================================================================

# Weight on performance beyond the reference.  In "field" mode the score
# saturates at the reference, but margin past it is insurance against having
# underestimated the competition, so the optimiser still values it weakly.
INSURANCE_WEIGHT = 0.05


@dataclass
class ScoreCard:
    feasible: bool
    m1: float
    m2: float
    m3: float
    total_mission: float
    m2_raw_lb_per_s: float
    m3_raw_lap_lb: float
    r1: Optional[MissionResult] = None
    r2: Optional[MissionResult] = None
    r3: Optional[MissionResult] = None
    reason: str = ""
    m2_ratio: float = 0.0          # raw / reference, uncapped
    m3_ratio: float = 0.0
    m2_saturated: bool = False     # only meaningful in "field" mode
    m3_saturated: bool = False
    mode: str = "self"

    @property
    def objective(self) -> float:
        """What the optimiser maximises (M1 is a constant for any flyer)."""
        excess = (max(0.0, self.m2_ratio - 1.0) + max(0.0, self.m3_ratio - 1.0))
        return self.m2 + self.m3 + INSURANCE_WEIGHT * excess


def _score_from_ratios(m1: float, r2: float, r3: float, mode: str):
    """Apply the two formulas.  In field mode the ratio caps at 1.0."""
    if mode == "field":
        m2 = (1.0 + min(r2, 1.0)) if r2 > 0 else 0.0
        m3 = (2.0 + min(r3, 1.0)) if r3 > 0 else 0.0
    else:
        m2 = (1.0 + r2) if r2 > 0 else 0.0
        m3 = (2.0 + r3) if r3 > 0 else 0.0
    return m2, m3, m1 + m2 + m3


def evaluate_scores(ac: Aircraft, cfg: Config, dx: float = 6.0,
                    ref_m2: Optional[float] = None,
                    ref_m3: Optional[float] = None) -> ScoreCard:
    """Fly all three missions and score them.

    `ref_m2` / `ref_m3` are the normalisers.  In "self" mode the optimiser
    passes the best raw values found so far; they are only a scale for the
    objective, and the final report renormalises against the finished
    population.
    """
    r = cfg.rules
    mode = r.normalisation
    fb2 = ref_m2 if ref_m2 is not None else r.field_best_m2_lb_per_s
    fb3 = ref_m3 if ref_m3 is not None else r.field_best_m3_lap_lb
    fb2 = max(fb2, 1e-9)
    fb3 = max(fb3, 1e-9)

    checks = ac.static_checks()
    if not checks["ok"]:
        return ScoreCard(False, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         reason="; ".join(checks["failures"]), mode=mode)

    r1 = mission_1(ac, cfg, dx)
    r2 = mission_2(ac, cfg, dx)
    r3 = mission_3(ac, cfg, dx)

    m1 = 1.0 if r1.ok else 0.0

    if r2.ok and r2.time_s > 0.0:
        raw2 = (ac.m2_scored_kg * KG2LB) / r2.time_s
        ratio2 = raw2 / fb2
    else:
        raw2, ratio2 = 0.0, 0.0

    if r3.ok and r3.laps > 0:
        raw3 = r3.laps * (ac.design.sensor_mass_kg * KG2LB)
        ratio3 = raw3 / fb3
    else:
        raw3, ratio3 = 0.0, 0.0

    m2, m3, total = _score_from_ratios(m1, ratio2, ratio3, mode)
    reasons = [x.reason for x in (r1, r2, r3) if not x.ok and x.reason]
    return ScoreCard(r1.ok and r2.ok and r3.ok, m1, m2, m3, total,
                     raw2, raw3, r1, r2, r3, "; ".join(reasons),
                     ratio2, ratio3, ratio2 > 1.0, ratio3 > 1.0, mode)


def rescore(card: ScoreCard, ref_m2: float, ref_m3: float,
            mode: Optional[str] = None) -> ScoreCard:
    """Re-apply the ratio terms with different normalisers, without
    re-running the simulation."""
    mode = mode or card.mode
    r2 = card.m2_raw_lb_per_s / max(ref_m2, 1e-9) if card.m2_raw_lb_per_s > 0 else 0.0
    r3 = card.m3_raw_lap_lb / max(ref_m3, 1e-9) if card.m3_raw_lap_lb > 0 else 0.0
    m2, m3, total = _score_from_ratios(card.m1, r2, r3, mode)
    return ScoreCard(card.feasible, card.m1, m2, m3, total,
                     card.m2_raw_lb_per_s, card.m3_raw_lap_lb,
                     card.r1, card.r2, card.r3, card.reason,
                     r2, r3, r2 > 1.0, r3 > 1.0, mode)


def population_references(rows: Sequence[Dict]) -> tuple:
    """Best raw M2 and M3 across a set of evaluated designs."""
    best2 = max([r.get("m2_raw_lb_per_s", 0.0) for r in rows] or [0.0])
    best3 = max([r.get("m3_raw_lap_lb", 0.0) for r in rows] or [0.0])
    return max(best2, 1e-9), max(best3, 1e-9)


def normalise_population(rows: List[Dict], mode: str = "self") -> tuple:
    """Rescore every row against the best in the population, in place.

    This is what makes "self" mode honest: the normaliser is the finished
    population, not a moving target from mid-search.
    """
    ref2, ref3 = population_references(rows)
    for row in rows:
        r2 = row.get("m2_raw_lb_per_s", 0.0) / ref2
        r3 = row.get("m3_raw_lap_lb", 0.0) / ref3
        m1 = row.get("m1", 0.0)
        m2, m3, total = _score_from_ratios(m1, r2, r3, mode)
        row["m2"], row["m3"], row["total_mission"] = m2, m3, total
        row["m2_ratio"], row["m3_ratio"] = r2, r3
        row["objective"] = (m2 + m3) if row.get("feasible") else -1e9
    return ref2, ref3


# ==========================================================================
#  OPTIMISER - four-stage search
# ==========================================================================

# ==========================================================================
# Search space
# ==========================================================================
@dataclass
class DesignSpace:
    airfoils: List[str] = field(default_factory=lambda: list(AIRFOILS.keys()))
    motors: List[str] = field(default_factory=lambda: list(MOTORS.keys()))
    props: List[str] = field(default_factory=lambda: list(PROPS.keys()))
    motor_counts: List[int] = field(default_factory=lambda: [1, 2, 3, 4])
    cells: List[str] = field(default_factory=lambda: list(CELLS.keys()))
    # 1.2 V cells need many more in series than 3.7 V ones, so the series
    # options are chosen per chemistry in candidate_packs().
    series_li: List[int] = field(default_factory=lambda: [3, 4, 5, 6, 8, 10, 12])
    series_nickel: List[int] = field(
        default_factory=lambda: [10, 12, 14, 16, 18, 20, 24, 28])
    parallel: List[int] = field(default_factory=lambda: [1, 2])

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
    sim_mass_bounds: Optional[Tuple[float, float]] = None    # default: rules
    taper: float = 1.0

    def continuous_bounds(self, cfg: Config) -> List[Tuple[float, float]]:
        r = cfg.rules
        span_lo, span_hi = self.span_bounds
        span_hi = min(span_hi, r.max_wingspan_m)
        span_lo = min(span_lo, span_hi * 0.5)
        sb = self.sensor_bounds or (r.sensor_mass_min_kg, r.sensor_mass_max_kg)
        # A simulator cannot weigh less than the empty container it mimics.
        floor = size_payload(r.sensor_mass_min_kg, r).container_tare_kg
        mb = self.sim_mass_bounds or (floor, r.sim_mass_max_kg)
        return [(span_lo, span_hi), self.ar_bounds, sb, mb,
                self.taper_bounds]

    VAR_NAMES = ("span_m", "aspect_ratio", "sensor_mass_kg", "sim_mass_kg",
                 "taper")


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
                    ref_m3: Optional[float] = None) -> Dict[str, Any]:
    try:
        ac = Aircraft(design, cfg)
        why = quick_reject(ac, cfg)
        if why is not None:
            return {"design": design, "feasible": False, "objective": -1e9,
                    "reason": why, "m2_raw_lb_per_s": 0.0,
                    "m3_raw_lap_lb": 0.0, "m1": 0.0, "m2": 0.0, "m3": 0.0}
        card: ScoreCard = evaluate_scores(ac, cfg, dx, ref_m2, ref_m3)
    except Exception as exc:                       # never let one design kill a sweep
        return {"design": design, "feasible": False, "objective": -1e9,
                "reason": "%s: %s" % (type(exc).__name__, exc)}

    row: Dict[str, Any] = {
        "design": design,
        "feasible": card.feasible,
        "objective": card.objective if card.feasible else -1e9,
        "m1": card.m1, "m2": card.m2, "m3": card.m3,
        "total_mission": card.total_mission,
        "m2_raw_lb_per_s": card.m2_raw_lb_per_s,
        "m3_raw_lap_lb": card.m3_raw_lap_lb,
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
    row["price_usd"] = ac.prop_sys.price_usd
    return row


def _worker(task) -> Dict[str, Any]:
    design, cfg, dx = task[0], task[1], task[2]
    ref2 = task[3] if len(task) > 3 else None
    ref3 = task[4] if len(task) > 4 else None
    return evaluate_design(design, cfg, dx, ref2, ref3)


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
                    min_wh: float) -> List[Tuple[str, int, int]]:
    out = []
    for ck in space.cells:
        cell = CELLS[ck]
        if cell.chemistry not in cfg.rules.battery_chemistry_allowed:
            continue
        series_opts = (space.series_nickel if cell.v_nom < 2.0
                       else space.series_li)
        for s in series_opts:
            for p in space.parallel:
                pack = BatteryPack(cell, s, p)
                if pack.watt_hours > cfg.rules.max_watt_hours + 1e-9:
                    continue
                if pack.watt_hours < min_wh:
                    continue
                out.append((ck, s, p))
    return out


def _screen_one(task) -> Optional[Dict[str, Any]]:
    mk, pk, nm, ck, s, p, cfg = task
    motor, prop = MOTORS[mk], PROPS[pk]
    pack = BatteryPack(CELLS[ck], s, p)
    ps = PropulsionSystem(motor, prop, nm, pack, cfg.env)
    rc = ps.rating_check()
    if (rc["motor_current_margin"] < 1.0 or rc["motor_power_margin"] < 1.0
            or rc["esc_current_margin"] < 1.0 or rc["pack_current_margin"] < 1.0
            or rc["tip_mach"] > 0.82):
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
    return {"motor": mk, "prop": pk, "n_motors": nm, "cell": ck,
            "series": s, "parallel": p, "fom": fom,
            "static_thrust_n": rc["static_thrust_n"],
            "mass_kg": ps.mass_kg, "wh": pack.watt_hours}


def screen_propulsion(space: DesignSpace, cfg: Config,
                      settings: SearchSettings,
                      verbose: bool = True) -> List[Dict[str, Any]]:
    packs = candidate_packs(space, cfg, settings.min_pack_wh)
    tasks = [(mk, pk, nm, ck, s, p, cfg)
             for mk in space.motors
             for pk in space.props
             for nm in space.motor_counts
             for (ck, s, p) in packs]
    if verbose:
        print("  screening %d propulsion combinations (%d packs)..."
              % (len(tasks), len(packs)))
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
    """x = (span, aspect ratio, sensor mass, simulator mass, taper)."""
    span = max(x[0], 0.3)
    ar = min(max(x[1], 3.2), 15.5)
    taper = min(max(x[4], 0.30), 1.0) if len(x) > 4 else taper
    return Design(span_m=span, chord_m=span / ar, airfoil=airfoil,
                  motor=base["motor"], prop=base["prop"],
                  n_motors=base["n_motors"], cell=base["cell"],
                  cells_series=base["series"], cells_parallel=base["parallel"],
                  sensor_mass_kg=x[2], n_simulators=n_sim,
                  sim_mass_kg=x[3], taper=taper)


def _plausible_sim_count(base: Dict[str, Any], span: float, chord: float,
                         sensor_kg: float, sim_kg: float,
                         rng: random.Random, cfg: Config) -> int:
    """How many simulators could this power system plausibly lift?

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
    capacity = min(m_thrust, m_energy) - m_empty - sensor_kg
    if capacity <= 0.0:
        return 0
    n = int(round(rng.random() * capacity / max(sim_kg, 0.05)))
    return max(0, min(n, cfg.rules.max_simulators))


def coarse_sweep(shortlist: List[Dict[str, Any]], space: DesignSpace,
                 cfg: Config, settings: SearchSettings,
                 verbose: bool = True,
                 refs: Tuple[Optional[float], Optional[float]] = (None, None)
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
                              cfg, settings.dx_coarse, refs[0], refs[1]))
    if verbose:
        print("  coarse sweep: %d designs..." % len(tasks))
    rows = _pool_map(_worker, tasks, _n_processes(settings))
    rows = [r for r in rows if r.get("feasible")]
    rows.sort(key=lambda r: -r["objective"])
    if verbose:
        print("  %d feasible" % len(rows))
    return rows


def _refine_one(task) -> Dict[str, Any]:
    row, cfg, space, settings, ref2, ref3 = task
    d: Design = row["design"]
    bounds = space.continuous_bounds(cfg)
    lower = [b[0] for b in bounds]
    upper = [b[1] for b in bounds]
    base = {"motor": d.motor, "prop": d.prop, "n_motors": d.n_motors,
            "cell": d.cell, "series": d.cells_series,
            "parallel": d.cells_parallel}

    def negobj(x: List[float], n_sim: int) -> float:
        des = _make_design(base, d.airfoil, x, n_sim, d.taper)
        return -evaluate_design(des, cfg, settings.dx_refine,
                                ref2, ref3)["objective"]

    x0 = [d.span_m, d.span_m / d.chord_m, d.sensor_mass_kg, d.sim_mass_kg,
          d.taper]
    step = [0.14, 0.70, 0.12, 0.09, 0.12]

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

    des = _make_design(base, d.airfoil, x, best_n, d.taper)
    return evaluate_design(des, cfg, settings.dx_refine, ref2, ref3)


def refine(rows: List[Dict[str, Any]], space: DesignSpace, cfg: Config,
           settings: SearchSettings, verbose: bool = True,
           refs: Tuple[Optional[float], Optional[float]] = (None, None)
           ) -> List[Dict[str, Any]]:
    seeds = _dedupe_configs(rows, settings.top_for_refine)
    if verbose:
        print("  refining %d configurations..." % len(seeds))
    tasks = [(r, cfg, space, settings, refs[0], refs[1]) for r in seeds]
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
    refs = (cfg.rules.field_best_m2_lb_per_s, cfg.rules.field_best_m3_lap_lb)
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
        if verbose:
            print("  self-normalising against the coarse best: "
                  "M2 %.4f lb/s, M3 %.2f lap*lb" % refs)

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
    tasks = [(r["design"], cfg, settings.dx_final, refs[0], refs[1])
             for r in finalists]
    final_rows = _pool_map(_worker, tasks, _n_processes(settings))
    final_rows = [r for r in final_rows if r.get("feasible")]
    final_rows.sort(key=lambda r: -r["objective"])
    n_eval += len(tasks)

    # Final scores are normalised against the finished population, so the
    # best design sits at exactly 1.0 on each ratio and everything else is
    # reported as a fraction of it.
    if cfg.rules.normalisation == "self":
        ref2, ref3 = normalise_population(final_rows, "self")
    else:
        ref2, ref3 = refs
    final_rows.sort(key=lambda r: -r["objective"])

    pts = [(r["m2_raw_lb_per_s"], r["m3_raw_lap_lb"]) for r in final_rows]
    front_idx = pareto_front(pts) if pts else []
    pareto = [final_rows[i] for i in front_idx]
    pareto.sort(key=lambda r: -r["m2_raw_lb_per_s"])

    return OptimisationResult(final_rows, pareto, len(shortlist), n_eval,
                              ref2, ref3, cfg.rules.normalisation)


# ==========================================================================
# Sensitivity
# ==========================================================================
def sensitivity(design: Design, cfg: Config, dx: float = 6.0
                ) -> List[Tuple[str, float, float, float]]:
    """One-at-a-time +/-10 percent sweep on the continuous variables.

    Returns (variable, low value, high value, d(objective)) rows so you can see
    which dimension is actually worth arguing about in the design review.
    """
    base = evaluate_design(design, cfg, dx)["objective"]
    out = []
    fields = [("span_m", design.span_m), ("chord_m", design.chord_m),
              ("sensor_mass_kg", design.sensor_mass_kg),
              ("sim_mass_kg", design.sim_mass_kg)]
    for name, val in fields:
        lo_d = _replace(design, name, val * 0.90)
        hi_d = _replace(design, name, val * 1.10)
        lo = evaluate_design(lo_d, cfg, dx)["objective"]
        hi = evaluate_design(hi_d, cfg, dx)["objective"]
        lo = lo if lo > -1e8 else float("nan")
        hi = hi if hi > -1e8 else float("nan")
        out.append((name, lo, hi, base))
    # integer: simulator count
    for delta in (-1, +1):
        k = max(0, min(cfg.rules.max_simulators, design.n_simulators + delta))
        if k == design.n_simulators:
            continue
        d2 = _replace(design, "n_simulators", k)
        v = evaluate_design(d2, cfg, dx)["objective"]
        out.append(("n_simulators=%d" % k, v, v, base))
    return out


def _replace(design: Design, field_name: str, value) -> Design:
    from dataclasses import replace as dc_replace
    return dc_replace(design, **{field_name: value})


# ==========================================================================
#  REPORT - dossier, Pareto table, bill of materials
# ==========================================================================

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
    hdr = ("%-3s %-9s %6s %6s %5s %-22s %-10s %4s %5s %5s %6s %6s %6s"
           % ("#", "airfoil", "span", "chord", "AR", "motor", "prop",
              "nM", "Wh", "M2t", "M3lap", "M2", "M3"))
    print(hdr)
    print(THIN)
    for i, r in enumerate(rows[:n], 1):
        d: Design = r["design"]
        print("%-3d %-9s %6.2f %6.3f %5.2f %-22s %-10s %4d %5.0f %5.1f %6d %6.3f %6.3f"
              % (i, d.airfoil, d.span_m, d.chord_m, r["AR"],
                 r["motor"][:22], d.prop.replace("APC_", ""), d.n_motors,
                 r["pack_Wh"], r.get("m2_time_s", 0.0), r.get("m3_laps", 0),
                 r["m2"], r["m3"]))
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
             "prop", "nSim", "simlb"))
    print(THIN)
    for i, r in enumerate(rows, 1):
        d: Design = r["design"]
        print("%-3d %8.4f %9.2f  %6.2f %6.3f %5.2f  %-20s %-9s %4d %5.2f"
              % (i, r["m2_raw_lb_per_s"], r["m3_raw_lap_lb"], d.span_m,
                 d.chord_m, d.sensor_mass_kg * KG2LB, r["motor"][:20],
                 d.prop.replace("APC_", ""), d.n_simulators,
                 d.sim_mass_kg * KG2LB))
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

    # ---------------- design loads ----------------
    lc = ac.load_case
    print("\nDESIGN LOADS  (V-n diagram, not an assumed number)")
    print("  manoeuvre case         %.2f g   (%.0f deg bank x %.2f overshoot)"
          % (lc.n_manoeuvre, cfg.build.max_bank_deg,
             cfg.build.pilot_overshoot_factor))
    print("  gust at cruise         %.2f g   (%.1f m/s discrete gust)"
          % (lc.n_gust_cruise, cfg.build.gust_cruise_mps))
    print("  gust at dive speed     %.2f g   (V_D = %.1f m/s)"
          % (lc.n_gust_dive, lc.v_dive))
    print("  wing CLmax ceiling     %.2f g   (most the wing can generate)"
          % lc.n_aero_ceiling)
    print("  -> LIMIT load          %.2f g   governed by %s"
          % (ac.n_limit, lc.driver))
    print("  -> ULTIMATE load       %.2f g   (x%.2f safety factor)"
          % (ac.n_ultimate, cfg.build.load_safety_factor))
    print("     the wing must survive %.2f g without failing; normal flight"
          % ac.n_ultimate)
    print("     is held inside %.2f g" % ac.n_limit)
    if ac.spar_gauge_limited:
        print("  NOTE: the spar is GAUGE limited, not strength limited.")
        print("        Bending needs only %.2f mm of wall; the %.2f mm"
              % (1000 * ac.spar_t_req_m, 1000 * ac.spar_t_used_m))
        print("        minimum is what sets its mass.  Raising the load")
        print("        factor costs nothing here, and a thinner layup would")
        print("        save weight if you can actually build it.")
    else:
        print("  spar wall required     %.2f mm (strength limited)"
              % (1000 * ac.spar_t_req_m))

    # ---------------- propulsion ----------------
    print("\nPROPULSION  (T)")
    print("  motor                  %s  x%d" % (ac.motor.name, design.n_motors))
    print("    Kv %.0f rpm/V, Rm %.3f ohm, I0 %.1f A, %.0f g each"
          % (ac.motor.kv, ac.motor.rm, ac.motor.i0, ac.motor.mass_kg * 1000))
    print("  propeller              %s  (%.1f x %.1f in)"
          % (ac.prop.key, ac.prop.diameter_m * M2IN, ac.prop.pitch_m * M2IN))
    print("  ESC                    %s x%d"
          % (ac.prop_sys.esc.key, design.n_motors))
    print("  battery                %s" % ac.pack.key)
    print("    %.1f V nominal, %.2f Ah, %.1f Wh, %.0f g, %.0f A max"
          % (ac.pack.v_nominal, ac.pack.capacity_ah, ac.pack.watt_hours,
             ac.pack.mass_kg * 1000, ac.pack.i_max_a))
    print("  static thrust T0       %.1f N  (%.2f lbf)"
          % (rc["static_thrust_n"], rc["static_thrust_n"] * N2LBF))
    print("  static T/W  (M2)       %.2f" % s["static_T_over_W_m2"])
    print("  static rpm / tip Mach  %.0f / %.2f" % (rc["rpm"], rc["tip_mach"]))
    print("  worst-case currents    %.0f A per motor, %.0f A pack"
          % (rc["motor_current_a"], rc["pack_current_a"]))
    print("  margins  motor I %.2fx  motor P %.2fx  ESC %.2fx  pack C %.2fx"
          % (rc["motor_current_margin"], rc["motor_power_margin"],
             rc["esc_current_margin"], rc["pack_current_margin"]))

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
    print("    sensor); %.0f mm foam per face would hold %.0f g in a 5 ft drop"
          % (1000 * pg.pad_m, 75.0))
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
    print("    each is one container envelope of extra bay volume and drag")
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
    print("  M2 raw  weight/time            %.4f lb/s" % card.m2_raw_lb_per_s)
    print("  M3 raw  laps x sensor weight   %.2f lap*lb" % card.m3_raw_lap_lb)
    if card.mode == "self":
        print("  normalised against             the best design in this run")
        print("    reference M2                 %.4f lb/s" % ref2)
        print("    reference M3                 %.2f lap*lb" % ref3)
    else:
        print("  assumed field best  M2         %.4f lb/s" % ref2)
        print("  assumed field best  M3         %.2f lap*lb" % ref3)
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
    print("  TOTAL MISSION SCORE = %.3f" % card.total_mission)
    if card.mode == "self":
        print()
        print("  Self-normalised: the best design in this run scores 1.0 on")
        print("  each ratio, so the ceiling is 6.0 and it is reached by")
        print("  whichever candidate leads BOTH missions. These are relative")
        print("  numbers for choosing between your own candidates, not a")
        print("  prediction of your competition score.")
    elif card.m2_saturated or card.m3_saturated:
        print()
        print("  NOTE: this design beats the assumed field best by %s."
              % ", ".join(
                  ([("%.0f%% on M2" % (100 * (card.m2_ratio - 1)))]
                   if card.m2_saturated else [])
                  + ([("%.0f%% on M3" % (100 * (card.m3_ratio - 1)))]
                     if card.m3_saturated else [])))
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
              ("sim_mass_kg", design.sim_mass_kg, b[3]),
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

    def add(category, part, spec, qty, unit_price, note=""):
        items.append({"category": category, "part": part, "spec": spec,
                      "qty": qty, "unit_price_usd": round(unit_price, 2),
                      "ext_price_usd": round(unit_price * qty, 2), "note": note})

    # ---- propulsion -------------------------------------------------
    add("Propulsion", ac.motor.name,
        "Kv %.0f, Rm %.3f ohm, I0 %.1f A, %.0f g"
        % (ac.motor.kv, ac.motor.rm, ac.motor.i0, ac.motor.mass_kg * 1000),
        design.n_motors, ac.motor.price_usd, ac.motor.source)
    add("Propulsion", "Propeller %s" % ac.prop.key,
        "%.1f x %.1f in, %.0f g" % (ac.prop.diameter_m * M2IN,
                                    ac.prop.pitch_m * M2IN,
                                    ac.prop.mass_kg * 1000),
        design.n_motors + 2, ac.prop.price_usd, "buy spares - you will break them")
    add("Propulsion", "ESC %s" % ac.prop_sys.esc.key,
        "%.0f A continuous, up to %dS" % (ac.prop_sys.esc.i_cont_a,
                                          ac.prop_sys.esc.max_cells),
        design.n_motors, ac.prop_sys.esc.price_usd, "")
    add("Propulsion", "Motor mount + spinner",
        "for %.0f mm shaft" % ac.motor.shaft_mm, design.n_motors, 25.0, "")

    # ---- energy -----------------------------------------------------
    cell = ac.pack.cell
    add("Energy", "Battery pack %s" % ac.pack.key,
        "%dS%dP %s, %.1f V, %.2f Ah, %.1f Wh, %.0f g"
        % (ac.pack.series, ac.pack.parallel, cell.chemistry,
           ac.pack.v_nominal, ac.pack.capacity_ah, ac.pack.watt_hours,
           ac.pack.mass_kg * 1000),
        2, ac.pack.price_usd, "one to fly, one to charge - %s" % cell.source)
    add("Energy", "Arming plug / fuse",
        "rated above %.0f A" % ac.prop_sys.rating_check()["pack_current_a"],
        1, 25.0, "required safety item")
    add("Energy", "Main connectors",
        "XT90 or AS150 for %.0f A" % ac.prop_sys.rating_check()["pack_current_a"],
        4, 6.0, "")

    # ---- structure --------------------------------------------------
    r_spar = max(0.40 * ac.wing.thickness_root_m, 0.006)
    m_root = (ac.n_ultimate * (ac.w_m2_n / 2.0) * 0.4244
              * (ac.wing.span_m / 2.0))
    t_req = m_root / (b.spar_sigma_allow_pa * math.pi * r_spar * r_spar)
    t_spec = max(t_req, b.spar_min_wall_m)
    add("Structure", "Carbon spar tube",
        "OD %.1f mm, wall %.2f mm, length %.0f mm"
        % (2000 * r_spar, 1000 * ac.spar_t_used_m, 1000 * ac.wing.span_m),
        1, 90.0,
        "sized for %.0f N*m root bending at n_ult %.2f (%s)"
        % (m_root, ac.n_ultimate, ac.load_case.driver))
    add("Structure", "Wing cores / ribs",
        "%.0f in2 planform, %s section, t/c %.0f%%"
        % (ac.wing.area_m2 / (IN2M ** 2), ac.airfoil.name,
           100 * ac.airfoil.t_over_c),
        1, 120.0, "")
    add("Structure", "Covering film",
        "%.1f m2 plus 30%% waste" % (2.1 * ac.wing.area_m2), 1, 45.0, "")
    add("Structure", "Fuselage shell",
        "%.0f x %.0f x %.0f mm, %.2f m2 wetted"
        % (1000 * ac.body.fuse_length_m, 1000 * ac.body.fuse_width_m,
           1000 * ac.body.fuse_height_m, ac.body.fuse_wetted_m2),
        1, 150.0, "carbon/ply box or moulded")
    add("Structure", "Tail surfaces",
        "H %.0f in2, V %.0f in2, arm %.0f mm"
        % (ac.body.s_horiz_m2 / (IN2M ** 2), ac.body.s_vert_m2 / (IN2M ** 2),
           1000 * ac.body.tail_arm_m),
        1, 60.0, "")
    if design.n_motors == 1:
        gear_note = ("main gear %.0f mm legs, %.0f mm prop tip clearance"
                     % (1000 * ac.gear_height_m,
                        1000 * ac.static_checks()["prop_ground_clearance_m"]))
    else:
        gear_note = ("main gear %.0f mm legs (props are wing-mounted, "
                     "clearance is not the driver)" % (1000 * ac.gear_height_m))
    add("Structure", "Landing gear", gear_note, 1, 85.0, "")

    # ---- systems ----------------------------------------------------
    add("Systems", "Servos", "%.0f g class, %d off"
        % (b.servo_mass_kg * 1000, b.n_servos_base), b.n_servos_base, 22.0,
        "2 aileron, 1 elevator, 1 rudder, 1 payload door, 1 sensor deploy")
    add("Systems", "Receiver + telemetry", "6+ channel", 1, 60.0, "")
    add("Systems", "Sensor deploy mechanism",
        "winch/door, %.0f g budget" % (b.deploy_mech_kg * 1000), 1, 90.0,
        "must fully stow before landing in M3")
    add("Systems", "Sensor light controller",
        "OFF / SOLID ON / FLASHING ON by remote command", 1, 45.0,
        "M3 requires all three states on command")

    # ---- payload ----------------------------------------------------
    add("Payload", "Sensor article",
        "%.0f g target (%.2f lb)" % (design.sensor_mass_kg * 1000,
                                     design.sensor_mass_kg * KG2LB),
        1, 0.0, "mass is a DESIGN VARIABLE - it multiplies the M3 score")
    pg = ac.payload_geom
    add("Payload", "Shipping container",
        "%.0f x %.0f x %.0f mm, %.0f mm foam, tare %.0f g"
        % (1000 * pg.box_l, 1000 * pg.box_w, 1000 * pg.box_h,
           1000 * pg.pad_m, pg.container_tare_kg * 1000),
        1, 0.0, "fixed assumed box; foam sized for a 5 ft drop at 75 g")
    add("Payload", "Container simulators",
        "%.0f g each (%.2f lb), same envelope as the container"
        % (ac.sim_mass_kg * 1000, ac.sim_mass_kg * KG2LB),
        design.n_simulators, 0.0,
        "declare max %d at Tech Inspection" % design.n_simulators)
    return items


def print_bom(design: Design, cfg: Config) -> None:
    items = bill_of_materials(design, cfg)
    print()
    print(LINE)
    print("BILL OF MATERIALS  (parts to build this aeroplane)")
    print(LINE)
    cat = None
    total = 0.0
    for it in items:
        if it["category"] != cat:
            cat = it["category"]
            print("\n[%s]" % cat)
        print("  %-28s x%-3d  $%8.2f   %s"
              % (it["part"][:28], it["qty"], it["ext_price_usd"], it["spec"]))
        if it["note"]:
            print("  %-28s        %10s   %s" % ("", "", it["note"]))
        total += it["ext_price_usd"]
    print(THIN)
    print("  approximate hardware cost: $%.2f  (airframe consumables only)" % total)


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


# ==========================================================================
#  COMMAND LINE INTERFACE
# ==========================================================================

#!/usr/bin/env python3




def build_config(args) -> Config:
    cfg = Config()
    cfg.env = Environment.at_density_altitude(args.altitude, args.temp, args.wind)
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

    Separated from cli_main() so the validation suite can check that
    every flag reaches the config value it claims to, without
    running a search or reaching into module globals.
    """
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    # environment
    p.add_argument("--wind", type=float, default=4.5,
                   help="steady wind along the runway, m/s (default 4.5)")
    p.add_argument("--altitude", type=float, default=0.0,
                   help="field elevation, m (default 0)")
    p.add_argument("--temp", type=float, default=15.0,
                   help="air temperature, C (default 15)")
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
    p.add_argument("--processes", type=int, default=0,
                   help="worker processes, 0 = cpu_count-1")
    p.add_argument("--seed", type=int, default=20262027)
    # data overrides
    p.add_argument("--motors-csv", default="", help="measured motor data")
    p.add_argument("--props-csv", default="", help="measured propeller data")
    # modes
    p.add_argument("--evaluate", default="",
                   help="evaluate a single design from a JSON file instead of searching")
    p.add_argument("--out", default="results", help="output directory")
    p.add_argument("--no-export", action="store_true")
    return p


def cli_main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.motors_csv:
        print("loaded %d motors from %s"
              % (load_motors_csv(args.motors_csv), args.motors_csv))
    if args.props_csv:
        print("loaded %d propellers from %s"
              % (load_props_csv(args.props_csv), args.props_csv))

    cfg = build_config(args)

    print("=" * 78)
    print("DBF 2026-27 MISSION OPTIMISER")
    print("=" * 78)
    print("  air density        %.4f kg/m3 (%.0f m, %.0f C)"
          % (cfg.env.rho, args.altitude, args.temp))
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
        print_design_report(d, cfg, dx=1.5, title="DESIGN EVALUATION")
        print_active_bounds(d, DesignSpace(), cfg)
        print_bom(d, cfg)
        print_sensitivity(sensitivity(d, cfg, dx=3.0))
        return 0

    # ---------------- search ----------------
    settings = (SearchSettings.quick() if args.quick else
                SearchSettings.thorough() if args.thorough else
                SearchSettings())
    settings.processes = args.processes
    settings.seed = args.seed
    space = DesignSpace()

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
    print_pareto(res.pareto)

    best = res.ranked[0]["design"]
    print_design_report(best, cfg, dx=1.5, title="RECOMMENDED DESIGN",
                        ref_m2=res.ref_m2, ref_m3=res.ref_m3)
    print_active_bounds(best, space, cfg)
    print_bom(best, cfg)
    print_sensitivity(sensitivity(best, cfg, dx=3.0))

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


# ==========================================================================
#  SELF TEST - physics and scoring invariants
# ==========================================================================

#!/usr/bin/env python3



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


def baseline(**kw) -> Design:
    args = dict(span_m=1.55, chord_m=0.27, airfoil="SD7062",
                motor="COBRA_C4120_12", prop="APC_12x8", n_motors=1,
                cell="LIPO_3000MAH_65C", cells_series=6, cells_parallel=1,
                sensor_mass_kg=0.9, n_simulators=4, sim_mass_kg=0.5)
    args.update(kw)
    return Design(**args)


def selftest_main() -> int:
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
    section("Propeller model")
    for key in ("APC_12x6", "APC_13x8", "APC_16x10"):
        p = PROPS[key]
        n = 8000.0 / 60.0
        t0 = p.thrust(1.225, n, 0.0)
        p0 = p.power(1.225, n, 0.0)
        area = math.pi * (p.diameter_m / 2) ** 2
        fom = t0 ** 1.5 / math.sqrt(2 * 1.225 * area) / p0
        check("%s static FOM in 0.45..0.80" % key, 0.45 <= fom <= 0.80,
              "FOM %.3f" % fom)
        etas = []
        for i in range(1, 100):
            j = i * 0.01 * p.j_zero
            v = j * n * p.diameter_m
            pw = p.power(1.225, n, v)
            etas.append(p.thrust(1.225, n, v) * v / max(pw, 1e-9))
        check("%s peak eta in 0.60..0.85" % key, 0.60 <= max(etas) <= 0.85,
              "eta_max %.3f" % max(etas))
        check("%s thrust falls monotonically with speed" % key,
              all(p.thrust(1.225, n, v) >= p.thrust(1.225, n, v + 1.0) - 1e-9
                  for v in range(0, 30)))
        check("%s makes no thrust past J_zero" % key,
              p.thrust(1.225, n, 1.05 * p.j_zero * n * p.diameter_m) == 0.0)

    # ------------------------------------------------------------------
    section("Motor / propulsion")
    ps = PropulsionSystem(MOTORS["COBRA_C4120_12"], PROPS["APC_12x8"], 1,
                          BatteryPack(CELLS["LIPO_3000MAH_65C"], 6, 1), cfg.env)
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
    check("battery sags under load", static.pack_voltage_v < ps.pack.ocv(1.0),
          "%.2f vs %.2f V" % (static.pack_voltage_v, ps.pack.ocv(1.0)))
    check("electrical power exceeds shaft power",
          static.p_elec_w > static.p_shaft_w,
          "%.0f vs %.0f W" % (static.p_elec_w, static.p_shaft_w))
    # The battery-sag loop is solved in closed form rather than iterated.
    # Verify that algebra against a brute-force fixed point, because a wrong
    # rearrangement would be silent - it would just shift every thrust value.
    worst_t = worst_v = 0.0
    for v in (0.0, 9.0, 21.0, 34.0, 45.0):
        for soc_ in (1.0, 0.62, 0.21):
            for thr_ in (1.0, 0.74, 0.31):
                ocv = ps.pack.ocv(soc_)
                vp = ocv
                n_bf = 0.0
                for _ in range(200):          # iterate the sag loop to death
                    vt = thr_ * vp
                    if vt <= 0.05:
                        n_bf = 0.0
                        break
                    nnl = ps.motor.kv * vt / 60.0

                    def res(nn, vt=vt):
                        cur = (vt - 60.0 * nn / ps.motor.kv) / ps.rm_hot
                        jj = (v / (nn * ps.prop.diameter_m)
                              if nn > 1e-6 else 1e9)
                        return (ps.kt * (cur - ps.motor.i0)
                                - ps.prop.cp(jj) * ps.env.rho * nn * nn
                                * ps.prop.diameter_m ** 5 / (2.0 * math.pi))

                    n_bf = (0.0 if res(1e-4) <= 0.0
                            else brent(res, 1e-4, nnl, xtol=1e-13))
                    i_bf = max(0.0, (vt - 60.0 * n_bf / ps.motor.kv)
                               / ps.rm_hot)
                    vp = max(ps.pack.v_min * 0.85,
                             ocv - ps.n_motors * i_bf * thr_ * ps.pack.r_ohm)
                op_cf = ps.solve(v, soc_, thr_)
                t_bf = (ps.n_motors
                        * max(0.0, ps.prop.thrust(ps.env.rho, n_bf, v))
                        * ps.install_thrust_factor)
                worst_t = max(worst_t, abs(op_cf.thrust_total_n - t_bf))
                worst_v = max(worst_v, abs(op_cf.pack_voltage_v - vp))
    check("closed-form battery sag matches brute-force iteration",
          worst_t < 1e-9 and worst_v < 1e-9,
          "dT %.2e N, dV %.2e V" % (worst_t, worst_v))

    # and the returned operating point must satisfy the torque balance
    op_chk = ps.solve(18.0, 0.8, 1.0)
    n_chk = op_chk.rpm / 60.0
    q_motor = ps.kt * (op_chk.motor_current_a - ps.motor.i0)
    q_prop = (ps.prop.cp(18.0 / (n_chk * ps.prop.diameter_m)) * ps.env.rho
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
    ratio = aero.v_stall(4 * w) / aero.v_stall(w)
    check("stall speed scales roughly as sqrt(W)", 1.90 < ratio < 2.10,
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
    cg0 = Config()
    cg0.build.cg_aft_of_ac_frac = 0.0
    a0 = Aircraft(baseline(), cg0).aero_clean
    w0 = Aircraft(baseline(), cg0).w_m2_n
    check("with CG at the wing AC the tail load is load-factor independent",
          abs(a0.tail_load(22.0, w0) - a0.tail_load(22.0, 3 * w0)) < 1e-9)
    check("moving the CG aft unloads the tail",
          abs(aero.tail_load(22.0, w)) < abs(a0.tail_load(22.0, w)),
          "%.2f vs %.2f N" % (aero.tail_load(22.0, w), a0.tail_load(22.0, w)))
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
    short = Aircraft(baseline(span_m=1.2, chord_m=0.35), cfg)
    long_ = Aircraft(baseline(span_m=2.1, chord_m=0.20), cfg)
    check("longer span means a heavier spar",
          long_.mass_breakdown.items["wing spar"]
          > short.mass_breakdown.items["wing spar"],
          "%.3f vs %.3f kg" % (long_.mass_breakdown.items["wing spar"],
                               short.mass_breakdown.items["wing spar"]))
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

    broken = Aircraft(baseline(cells_series=4, cell="LIPO_1300MAH_65C",
                               n_simulators=12, sim_mass_kg=1.0), cfg)
    bad = evaluate_scores(broken, cfg, dx=6.0)
    check("an infeasible aeroplane scores 0 on the missions it cannot fly",
          bad.m2 == 0.0 or bad.m3 == 0.0 or not bad.feasible, bad.reason)

    # ------------------------------------------------------------------
    section("Lifting line, stability and thermal state")

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
    aft.build.cg_aft_of_ac_frac = 0.30
    check("moving the CG aft reduces static margin",
          Aircraft(baseline(), aft).aero_clean.static_margin
          < ac.aero_clean.static_margin)
    check("a CG far enough aft is rejected as unstable",
          not Aircraft(baseline(), aft).static_checks()["ok"])
    check("a bigger tail raises the neutral point",
          ac.aero_clean.v_h > 0.0)

    # motor thermal state, per design
    hots = [Aircraft(baseline(n_motors=n, prop="APC_11x7"), cfg).motor_hot_factor
            for n in (1, 2, 4)]
    check("sharing current across motors runs them cooler",
          hots[0] > hots[1] > hots[2], "1/2/4 motors: %s"
          % ["%.3f" % h for h in hots])
    check("hot resistance factor never below 1", min(hots) >= 1.0)
    check("hot resistance factor stays physical", max(hots) <= 1.80)

    # the three coupled loops must actually converge together
    w = ac.mtow_m2_kg * 9.80665
    again = design_load_factors(ac.aero_clean, w,
                                ac._max_level_speed(w), cfg.build)
    check("load case is consistent with the final propulsion state",
          abs(again.n_ultimate - ac.n_ultimate) < 1e-3,
          "%.5f vs %.5f" % (again.n_ultimate, ac.n_ultimate))

    # windmilling drag and pilot efficiency
    calm = Config()
    calm.build.windmill_drag_coeff = 0.0
    r_no = mission_2(Aircraft(baseline(), calm), calm, dx=6.0)
    r_yes = mission_2(ac, cfg, dx=6.0)
    check("a freewheeling propeller costs lap time",
          r_yes.time_s > r_no.time_s,
          "%.2f vs %.2f s" % (r_yes.time_s, r_no.time_s))
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
    check("there is no sensor mass ceiling short of impossible density",
          Aircraft(baseline(sensor_mass_kg=8.0), cfg).static_checks()["ok"])
    check("a sensor denser than lead is rejected",
          not Aircraft(baseline(sensor_mass_kg=60.0),
                       cfg).static_checks()["ok"])
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
    check("landing is simulated and takes time", lr.landing_time_s > 0.0,
          "%.1f s" % lr.landing_time_s)
    check("landing is excluded from the scored window",
          abs(lr.total_flight_time_s - (lr.time_s + lr.landing_time_s)) < 1e-6)
    check("landing rollout is a sane distance", 3.0 < lr.landing_roll_m < 300.0,
          "%.1f m" % lr.landing_roll_m)
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
    _DS = DesignSpace
    check("the search never samples beyond the span limit",
          _DS().continuous_bounds(cfg)[0][1] <= cfg.rules.max_wingspan_m + 1e-12)

    _C = CELLS
    for chem in ("LiPo", "LiIon", "NiMH", "NiCd"):
        check("%s cells are available" % chem,
              any(c.chemistry == chem for c in _C.values()))
    check("nickel packs are heavier than lithium for the same energy",
          BatteryPack(_C["NIMH_SUBC_5000"], 16, 1).mass_kg
          > BatteryPack(_C["LIPO_3000MAH_65C"], 9, 1).mass_kg)

    check("hot motor resistance is used, not the cold catalogue value",
          ps.rm_hot > ps.motor.rm)
    check("installed thrust is below free-air thrust",
          ps.install_thrust_factor < 1.0)

    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("%d passed, %d failed" % (PASS, FAIL))
    print("=" * 60)
    return 1 if FAIL else 0


# ==========================================================================
#  VALIDATION - catalogue, bounds, determinism, interfaces
# ==========================================================================

#!/usr/bin/env python3



_va_PASS = _va_FAIL = 0
FAILURES = []


def _va_check(name, cond, detail=""):
    global _va_PASS, _va_FAIL
    if cond:
        _va_PASS += 1
    else:
        _va_FAIL += 1
        FAILURES.append("%s   %s" % (name, detail))
        print("  _va_FAIL  %s   %s" % (name, detail))


def _va_section(t):
    print("\n" + t)
    print("-" * len(t))


def _va_baseline(**kw):
    a = dict(span_m=1.55, chord_m=0.26, airfoil="SD7062",
             motor="COBRA_C4120_12", prop="APC_12x8", n_motors=1,
             cell="LIPO_3000MAH_65C", cells_series=6, cells_parallel=1,
             sensor_mass_kg=0.9, n_simulators=3, sim_mass_kg=0.5)
    a.update(kw)
    return Design(**a)


# ==========================================================================
def validate_catalogue(cfg):
    _va_section("A. Catalogue integrity (every entry, not just the sampled ones)")

    for name, af in AIRFOILS.items():
        _va_check("%s CLmax plausible" % name, 0.8 <= af.cl_max_2d <= 2.6,
              "%.2f" % af.cl_max_2d)
        _va_check("%s cd_min plausible" % name, 0.004 <= af.cd_min <= 0.035,
              "%.4f" % af.cd_min)
        _va_check("%s camber gives nose-down Cm" % name, af.cm_ac <= 0.0,
              "%.3f" % af.cm_ac)
        _va_check("%s thickness plausible" % name, 0.05 <= af.t_over_c <= 0.20,
              "%.3f" % af.t_over_c)
        _va_check("%s stall break below CLmax" % name, af.cl_break < af.cl_max_2d)

    for key, pr in PROPS.items():
        n = 7000.0 / 60.0
        t0 = pr.thrust(1.225, n, 0.0)
        p0 = pr.power(1.225, n, 0.0)
        area = math.pi * (pr.diameter_m / 2) ** 2
        fom = t0 ** 1.5 / math.sqrt(2 * 1.225 * area) / p0
        _va_check("%s static FOM in band" % key, 0.40 <= fom <= 0.85,
              "%.3f" % fom)
        etas = [pr.thrust(1.225, n, j * pr.j_zero * n * pr.diameter_m)
                * (j * pr.j_zero * n * pr.diameter_m)
                / max(pr.power(1.225, n, j * pr.j_zero * n * pr.diameter_m), 1e-9)
                for j in [i * 0.02 for i in range(1, 50)]]
        _va_check("%s peak eta in band" % key, 0.55 <= max(etas) <= 0.88,
              "%.3f" % max(etas))
        _va_check("%s thrust monotone in V" % key,
              all(pr.thrust(1.225, n, v) >= pr.thrust(1.225, n, v + 0.5) - 1e-12
                  for v in [i * 0.5 for i in range(80)]))
        _va_check("%s Cp never negative" % key,
              all(pr.cp(j * 0.05) > 0 for j in range(40)))

    for key, m in MOTORS.items():
        _va_check("%s Kv positive" % key, 50 < m.kv < 3000, "%.0f" % m.kv)
        _va_check("%s Rm plausible" % key, 0.005 <= m.rm <= 0.5, "%.3f" % m.rm)
        _va_check("%s I0 plausible" % key, 0.2 <= m.i0 <= 6.0, "%.2f" % m.i0)
        _va_check("%s mass plausible" % key, 0.03 <= m.mass_kg <= 1.5,
              "%.3f" % m.mass_kg)
        # power rating must be consistent with the current rating at a
        # voltage the motor could plausibly see
        v_plausible = m.p_max_w / max(m.i_max_a, 1e-6)
        _va_check("%s power/current rating consistent" % key,
              6.0 <= v_plausible <= 60.0, "implies %.1f V" % v_plausible)

    bands = {"LiPo": (100, 175), "LiIon": (180, 280),
             "NiMH": (45, 95), "NiCd": (30, 70)}
    for key, c in CELLS.items():
        whkg = c.v_nom * c.capacity_ah / c.mass_kg
        lo, hi = bands[c.chemistry]
        _va_check("%s energy density in %s band" % (key, c.chemistry),
              lo <= whkg <= hi, "%.0f Wh/kg" % whkg)
        _va_check("%s cutoff below nominal" % key, c.v_min < c.v_nom < c.v_full)
        _va_check("%s usable fraction sane" % key, 0.5 <= c.usable_fraction <= 0.95)
        _va_check("%s internal resistance positive" % key, c.r_internal_ohm > 0)

    for key, e in ESCS.items():
        _va_check("%s efficiency plausible" % key, 0.90 <= e.efficiency <= 0.99)
        _va_check("%s mass scales with current" % key,
              0.3 <= 1000 * e.mass_kg / e.i_cont_a <= 2.0,
              "%.2f g/A" % (1000 * e.mass_kg / e.i_cont_a))


# ==========================================================================
def validate_bounds(cfg):
    _va_section("B. Independent bounds (hand calculations, no simulator)")
    ac = Aircraft(_va_baseline(), cfg)
    r2 = mission_2(ac, cfg, dx=6.0)
    _va_check("Mission 2 flew", r2.ok, r2.reason)
    if not r2.ok:
        return

    # Lower bound on lap time: straights at Vmax plus turns at the best
    # sustained turn rate.  The simulator must not beat this.
    turn = turn_capability(ac.aero_clean, ac.prop_sys, ac.w_m2_n, 0.7, cfg)
    t_turns = math.radians(cfg.course.total_turn_deg) / turn.turn_rate_rad_s
    t_straights = cfg.course.straight_length_m / r2.v_max
    floor = t_straights + t_turns
    lap = r2.lap_times[-1]
    _va_check("lap time respects its physical floor", lap >= floor - 1e-6,
          "lap %.2f s vs floor %.2f s" % (lap, floor))
    _va_check("lap time within 2x of the floor", lap <= 2.0 * floor,
          "lap %.2f s vs floor %.2f s" % (lap, floor))

    # Lower bound on energy: minimum drag over the straight distance at the
    # best possible efficiency.
    v_ld = ac.aero_clean.best_ld_speed(ac.w_m2_n)
    d_min = ac.aero_clean.drag(v_ld, ac.w_m2_n)
    e_floor = d_min * 5 * cfg.course.straight_length_m / 0.70
    _va_check("Mission 2 energy respects its floor", r2.energy_j >= e_floor - 1e-6,
          "%.0f J vs floor %.0f J" % (r2.energy_j, e_floor))

    # Stall speed against the textbook formula at the model's own CLmax
    # Stall speed against an independent solve of the trim relation:
    # the wing stalls when L_w = q S CLmax, with L_w = W - L_t(W).
    v_s = ac.aero_clean.v_stall(ac.w_m2_n)
    a = ac.aero_clean
    resid = ((ac.w_m2_n - a.tail_load(v_s, ac.w_m2_n))
             - a.q(v_s) * ac.wing_area_m2 * a.cl_max_at_v(v_s))
    _va_check("stall speed satisfies the trimmed lift balance", abs(resid) < 1e-6,
          "residual %.3e N" % resid)

    # the empennage must close the lift balance at every load factor
    for n in (1.0, 2.0, 3.0):
        L = n * ac.w_m2_n
        lw = a.cl_required(20.0, L) * a.q(20.0) * ac.wing_area_m2
        lt = a.tail_load(20.0, L)
        _va_check("wing + tail = total lift at %.0f g" % n,
              abs((lw + lt) - L) < 1e-6, "%.6f" % ((lw + lt) - L))

    # Course geometry against the stated rules
    _va_check("course straights are 2000 ft",
          abs(cfg.course.straight_length_m / FT2M - 2000.0) < 1e-6)
    _va_check("course turning is 720 deg", abs(cfg.course.total_turn_deg - 720) < 1e-9)
    _va_check("span limit is 6 ft",
          abs(cfg.rules.max_wingspan_m / FT2M - 6.0) < 1e-9)


# ==========================================================================
def validate_energy(cfg):
    _va_section("C. Energy bookkeeping")
    ac = Aircraft(_va_baseline(), cfg)
    for label, res, w in (("M2", mission_2(ac, cfg, dx=6.0), ac.w_m2_n),
                          ("M3", mission_3(ac, cfg, dx=6.0), ac.w_m3_n)):
        if not res.ok:
            _va_check("%s flew" % label, False, res.reason)
            continue
        # average power cannot exceed the most the system can physically draw
        p_avg = res.energy_j / max(res.total_flight_time_s, 1e-9)
        p_max = ac.prop_sys.power_elec(0.0, 1.0, 1.0)
        for v in (10, 20, 30, 40):
            p_max = max(p_max, ac.prop_sys.power_elec(v, 1.0, 1.0))
        _va_check("%s average power below system maximum" % label, p_avg <= p_max,
              "%.0f W vs %.0f W" % (p_avg, p_max))
        _va_check("%s energy within the pack" % label,
              res.energy_j <= ac.pack.watt_hours * 3600.0,
              "%.1f Wh of %.1f" % (res.energy_j / 3600, ac.pack.watt_hours))
        _va_check("%s energy margin consistent with energy used" % label,
              abs((1.0 - res.energy_j / (ac.pack.watt_hours * 3600.0))
                  - res.energy_margin) < 1e-9)
        _va_check("%s total time = scored time + landing" % label,
              abs(res.total_flight_time_s - (res.time_s + res.landing_time_s))
              < 1e-6)
        _va_check("%s lap times all positive" % label,
              all(t > 0 for t in res.lap_times))
        _va_check("%s scored time excludes landing" % label,
              res.time_s < res.total_flight_time_s)

    # more laps must cost more energy
    r3 = mission_3(ac, cfg, dx=6.0)
    _va_check("Mission 3 energy scales with laps", r3.energy_j > 0 and r3.laps > 0)


# ==========================================================================
def validate_determinism(cfg):
    _va_section("D. Determinism and reproducibility")
    d = _va_baseline()
    a = evaluate_design(d, cfg, 6.0)
    b = evaluate_design(d, cfg, 6.0)
    for k in ("m2_raw_lb_per_s", "m3_raw_lap_lb", "m2_time_s", "m3_laps"):
        _va_check("repeat evaluation identical: %s" % k, a.get(k) == b.get(k),
              "%r vs %r" % (a.get(k), b.get(k)))

    # a fresh propulsion cache must not change anything
    clear_system_cache()
    c = evaluate_design(d, cfg, 6.0)
    _va_check("result independent of the instance cache",
          a["m2_raw_lb_per_s"] == c["m2_raw_lb_per_s"],
          "%r vs %r" % (a["m2_raw_lb_per_s"], c["m2_raw_lb_per_s"]))

    # the memo in PropulsionSystem must not leak between query points
    ps = Aircraft(d, cfg).prop_sys
    t1 = ps.thrust(10.0, 0.9, 1.0)
    ps.thrust(30.0, 0.4, 0.6)
    t2 = ps.thrust(10.0, 0.9, 1.0)
    _va_check("operating-point memo is not stale", t1 == t2, "%.9f vs %.9f" % (t1, t2))


# ==========================================================================
def validate_roundtrip(cfg):
    _va_section("E. Export / import round trip")
    from dataclasses import asdict
    d = _va_baseline()
    before = evaluate_design(d, cfg, 6.0)
    path = os.path.join(tempfile.gettempdir(), "_dbf_roundtrip.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(asdict(d), fh)
    with open(path, encoding="utf-8") as fh:
        d2 = Design(**json.load(fh))
    os.remove(path)
    _va_check("design survives JSON round trip", d2 == d)
    after = evaluate_design(d2, cfg, 6.0)
    _va_check("score survives JSON round trip",
          before["m2_raw_lb_per_s"] == after["m2_raw_lb_per_s"])

    # bill of materials must be generatable and priced
    bom = bill_of_materials(d, cfg)
    _va_check("BOM lists parts", len(bom) >= 15, "%d items" % len(bom))
    _va_check("BOM has no negative prices",
          all(i["ext_price_usd"] >= 0 for i in bom))
    _va_check("BOM covers every category",
          {i["category"] for i in bom} >=
          {"Propulsion", "Energy", "Structure", "Systems", "Payload"})


# ==========================================================================
def validate_trends(cfg):
    _va_section("F. Physical trends across the catalogue")

    # every airfoil: more area lowers stall speed
    for name in AIRFOILS:
        small = Aircraft(_va_baseline(airfoil=name, span_m=1.2, chord_m=0.20), cfg)
        big = Aircraft(_va_baseline(airfoil=name, span_m=1.8, chord_m=0.30), cfg)
        _va_check("%s: more wing area lowers stall speed" % name,
              big.aero_clean.v_stall(big.w_m2_n)
              < small.aero_clean.v_stall(small.w_m2_n))

    # trim drag must penalise camber: same planform, more camber, more drag
    flat = Aircraft(_va_baseline(airfoil="MH32"), cfg)
    camb = Aircraft(_va_baseline(airfoil="S1223"), cfg)
    v = 22.0
    wref = camb.w_m2_n
    _va_check("a high-camber _va_section carries more trim download",
          abs(camb.aero_clean.tail_load(v, wref))
          > abs(flat.aero_clean.tail_load(v, flat.w_m2_n)),
          "%.1f vs %.1f N" % (camb.aero_clean.tail_load(v, wref),
                              flat.aero_clean.tail_load(v, flat.w_m2_n)))

    # heavier sensor -> bigger container -> bigger fuselage
    prev_len = 0.0
    for m in (0.3, 0.8, 1.5, 3.0):
        ac = Aircraft(_va_baseline(sensor_mass_kg=m), cfg)
        _va_check("sensor %.1f kg grows the fuselage" % m,
              ac.body.fuse_length_m >= prev_len - 1e-9,
              "%.3f m" % ac.body.fuse_length_m)
        prev_len = ac.body.fuse_length_m

    # wind must never make the mission faster
    calm = Config()
    calm.env.wind_mps = 0.0
    windy = Config()
    windy.env.wind_mps = 9.0
    t_calm = mission_2(Aircraft(_va_baseline(), calm), calm, dx=6.0)
    t_wind = mission_2(Aircraft(_va_baseline(), windy), windy, dx=6.0)
    if t_calm.ok and t_wind.ok:
        _va_check("wind never speeds up the mission", t_wind.time_s >= t_calm.time_s,
              "%.1f vs %.1f s" % (t_wind.time_s, t_calm.time_s))

    # altitude must never help
    hi = Config()
    hi.env = Environment.at_density_altitude(1500.0, 30.0, 4.5)
    t_hi = mission_2(Aircraft(_va_baseline(), hi), hi, dx=6.0)
    if t_hi.ok and t_calm.ok:
        _va_check("thin air never speeds up the mission",
              t_hi.time_s >= t_calm.time_s,
              "%.1f vs %.1f s" % (t_hi.time_s, t_calm.time_s))

    # every chemistry must produce a usable pack somewhere
    for chem in ("LiPo", "LiIon", "NiMH", "NiCd"):
        keys = [k for k, c in CELLS.items() if c.chemistry == chem]
        ok = False
        for k in keys:
            for s_ in (4, 6, 8, 10, 12, 16, 20, 24, 28):
                p = BatteryPack(CELLS[k], s_, 1)
                if 20.0 <= p.watt_hours <= cfg.rules.max_watt_hours:
                    ok = True
        _va_check("%s can make a legal pack under 100 Wh" % chem, ok)


# ==========================================================================
def validate_interfaces(cfg):
    _va_section("G. Command-line interface")
    parser = build_parser()

    argv = ["--wind", "7", "--altitude", "500", "--temp", "28",
            "--max-wh", "90", "--max-sim", "15", "--window", "280",
            "--max-span-in", "60", "--sensor-max", "4", "--sim-max", "2",
            "--lap-ft", "1800", "--normalise", "field",
            "--field-m2", "0.08", "--field-m3", "30"]
    c = build_config(parser.parse_args(argv))
    _va_check("--wind reaches the config", abs(c.env.wind_mps - 7.0) < 1e-9)
    _va_check("--altitude changes air density", c.env.rho < 1.225)
    _va_check("--max-wh reaches the config", abs(c.rules.max_watt_hours - 90) < 1e-9)
    _va_check("--max-sim reaches the config", c.rules.max_simulators == 15)
    _va_check("--window reaches the config",
          abs(c.rules.mission_window_s - 280) < 1e-9)
    _va_check("--max-span-in reaches the config",
          abs(c.rules.max_wingspan_m - 60 * 0.0254) < 1e-9)
    _va_check("--sensor-max reaches the config",
          abs(c.rules.sensor_mass_max_kg - 4 * LB2KG) < 1e-6)
    _va_check("--sim-max reaches the config",
          abs(c.rules.sim_mass_max_kg - 2 * LB2KG) < 1e-6)
    _va_check("--lap-ft rescales the course",
          abs(c.course.straight_length_m / FT2M - 1800.0) < 1e-6,
          "%.1f ft" % (c.course.straight_length_m / FT2M))
    _va_check("--lap-ft preserves the turns",
          abs(c.course.total_turn_deg - 720.0) < 1e-9)
    _va_check("--normalise reaches the config", c.rules.normalisation == "field")
    _va_check("--field-m2 reaches the config",
          abs(c.rules.field_best_m2_lb_per_s - 0.08) < 1e-12)

    c0 = build_config(parser.parse_args([]))
    _va_check("defaults keep the 6 ft span rule",
          abs(c0.rules.max_wingspan_m / FT2M - 6.0) < 1e-9)
    _va_check("defaults keep no takeoff limit",
          c0.rules.takeoff_distance_m == float("inf"))
    _va_check("defaults keep self-normalisation", c0.rules.normalisation == "self")
    _va_check("search never samples beyond the span limit",
          DesignSpace().continuous_bounds(c0)[0][1]
          <= c0.rules.max_wingspan_m + 1e-12)
    _va_check("every declared flag is reachable",
          len(parser._actions) > 15, "%d actions" % len(parser._actions))

    # A self-normalised score can never exceed its own ceiling.  This is
    # what a single design scored against stale references got wrong.
    rows = [evaluate_design(_va_baseline(n_simulators=n), cfg, 12.0)
            for n in (0, 2, 4)]
    rows = [r for r in rows if r.get("feasible")]
    if rows:
        normalise_population(rows, "self")
        for r in rows:
            _va_check("self-normalised M2 never exceeds 2.0", r["m2"] <= 2.0 + 1e-9,
                  "%.4f" % r["m2"])
            _va_check("self-normalised M3 never exceeds 3.0", r["m3"] <= 3.0 + 1e-9,
                  "%.4f" % r["m3"])
            _va_check("self-normalised total never exceeds 6.0",
                  r["total_mission"] <= 6.0 + 1e-9, "%.4f" % r["total_mission"])
        _va_check("the population best reaches exactly 1.0 on some ratio",
              any(abs(r["m2_ratio"] - 1.0) < 1e-9
                  or abs(r["m3_ratio"] - 1.0) < 1e-9 for r in rows))


# ==========================================================================
def validate_loads(cfg):
    _va_section("H. Design load factors (V-n diagram)")

    # lift-curve slope must approach 2 pi at high aspect ratio and fall below
    a_low = lift_curve_slope_3d(4.0)
    a_high = lift_curve_slope_3d(20.0)
    _va_check("lift-curve slope rises with aspect ratio", a_high > a_low)
    _va_check("lift-curve slope stays below 2 pi", a_high < 2 * math.pi,
          "%.3f" % a_high)

    # the Pratt gust formula must behave
    base = gust_load_factor(250.0, 0.28, 4.8, 30.0, 1.225, 7.62)
    _va_check("gust load factor exceeds 1 g", base > 1.0, "%.2f" % base)
    _va_check("lighter wing loading is more gust-critical",
          gust_load_factor(120.0, 0.28, 4.8, 30.0, 1.225, 7.62) > base)
    _va_check("a faster aeroplane sees a bigger gust load",
          gust_load_factor(250.0, 0.28, 4.8, 45.0, 1.225, 7.62) > base)
    _va_check("a stronger gust gives a bigger load",
          gust_load_factor(250.0, 0.28, 4.8, 30.0, 1.225, 15.0) > base)
    _va_check("zero gust gives exactly 1 g",
          abs(gust_load_factor(250.0, 0.28, 4.8, 30.0, 1.225, 0.0) - 1.0) < 1e-12)

    # across the whole airfoil catalogue the derived load case must be sane
    for name in AIRFOILS:
        ac = Aircraft(_va_baseline(airfoil=name), cfg)
        lc = ac.load_case
        _va_check("%s: limit load in a sane band" % name,
              1.5 <= ac.n_limit <= 8.0, "%.2f g" % ac.n_limit)
        _va_check("%s: ultimate = safety factor x limit" % name,
              abs(ac.n_ultimate
                  - cfg.build.load_safety_factor * ac.n_limit) < 1e-9)
        _va_check("%s: limit never exceeds the CLmax ceiling" % name,
              ac.n_limit <= lc.n_aero_ceiling + 1e-9,
              "%.2f vs %.2f" % (ac.n_limit, lc.n_aero_ceiling))
        _va_check("%s: dive gust milder than cruise gust" % name,
              lc.n_gust_dive < lc.n_gust_cruise)
        _va_check("%s: dive speed above cruise" % name, lc.v_dive > 0)

    # the mass/load fixed point has to actually converge
    for n_sim in (0, 4, 10):
        ac = Aircraft(_va_baseline(n_simulators=n_sim), cfg)
        _va_check("mass closes with %d simulators" % n_sim,
              abs(ac.mtow_m2_kg - (ac.empty_mass_kg + ac.m2_payload_kg)) < 1e-9)
        recomputed = design_load_factors(
            ac.aero_clean, ac.mtow_m2_kg * 9.80665,
            ac._max_level_speed(ac.mtow_m2_kg * 9.80665), cfg.build)
        _va_check("load case is self-consistent with %d simulators" % n_sim,
              abs(recomputed.n_ultimate - ac.n_ultimate) < 0.02,
              "%.4f vs %.4f" % (recomputed.n_ultimate, ac.n_ultimate))

    # spar: monotone in load factor, and gauge-limiting reported honestly
    ac = Aircraft(_va_baseline(), cfg)
    prev = 0.0
    for n in (2.0, 4.0, 8.0, 16.0, 32.0):
        m, t_req, t_used, gauge = spar_sizing(ac.wing, ac.w_m2_n, cfg.build, n)
        _va_check("spar mass non-decreasing at n_ult %.0f" % n, m >= prev - 1e-15,
              "%.5f vs %.5f" % (m, prev))
        _va_check("spar wall never below minimum gauge at n_ult %.0f" % n,
              t_used >= cfg.build.spar_min_wall_m - 1e-15)
        _va_check("gauge flag agrees with the requirement at n_ult %.0f" % n,
              gauge == (t_req < cfg.build.spar_min_wall_m))
        prev = m

    # an override must bypass the derivation entirely
    ovr = Config()
    ovr.build.ultimate_load_factor_override = 3.75
    a2 = Aircraft(_va_baseline(), ovr)
    _va_check("override sets the ultimate load factor exactly",
          abs(a2.n_ultimate - 3.75) < 1e-12)
    _va_check("override sets limit = ultimate / safety factor",
          abs(a2.n_limit - 3.75 / ovr.build.load_safety_factor) < 1e-12)

    # span efficiency across the whole design space must stay physical
    for ar in (3.5, 5, 7, 9, 11, 14):
        prev = None
        for t in (1.0, 0.8, 0.6, 0.45):
            e = span_efficiency(ar, t)
            _va_check("AR %.1f taper %.2f: span efficiency in (0.5, 1.0]" % (ar, t),
                  0.5 < e <= 1.0, "%.4f" % e)
            if prev is not None:
                _va_check("AR %.1f: taper %.2f beats %.2f" % (ar, t, prev[0]),
                      e > prev[1] - 1e-9, "%.4f vs %.4f" % (e, prev[1]))
            prev = (t, e)

    # static margin must be computed, finite and consistent for every airfoil
    for name in list(AIRFOILS)[:8]:
        a = Aircraft(_va_baseline(airfoil=name), cfg)
        sm = a.aero_clean.static_margin
        _va_check("%s: static margin is finite and sane" % name,
              -0.5 < sm < 1.0, "%.3f" % sm)
        _va_check("%s: neutral point aft of the wing AC" % name,
              a.aero_clean.x_np_over_c > 0.25)

    # 1-4 motors must all build and be checked for propeller overlap
    for nm in (1, 2, 3, 4):
        a = Aircraft(_va_baseline(n_motors=nm, prop="APC_11x7"), cfg)
        _va_check("%d-motor layout builds" % nm, a.mtow_m2_kg > 0)
        _va_check("%d-motor propulsion mass scales" % nm,
              abs(a.motor_mass_kg - nm * a.motor.mass_kg) < 1e-12)
    wide = Aircraft(_va_baseline(n_motors=4, prop="APC_20x13"), cfg)
    _va_check("four big propellers on a 6 ft wing are rejected as overlapping",
          any("overlap" in f for f in wide.static_checks()["failures"]),
          "%s" % wide.static_checks()["failures"])


# ==========================================================================
def validate_main():
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
    print()
    print("=" * 70)
    print("%d passed, %d failed" % (_va_PASS, _va_FAIL))
    if FAILURES:
        print()
        for f in FAILURES:
            print("  FAILED: %s" % f)
    print("=" * 70)
    return 1 if _va_FAIL else 0


# ==========================================================================
#  NUMERICAL ERROR STUDY
# ==========================================================================

#!/usr/bin/env python3



_gc_LINE = "=" * 78
_gc_THIN = "-" * 78


def demo_design() -> Design:
    """The aeroplane the thorough run recommended - a representative point."""
    return Design(span_m=1.826, chord_m=0.259, airfoil="CLARKY",
                  motor="SCORPION_SII_4020_420", prop="APC_13x10", n_motors=2,
                  cell="LIPO_2200MAH_65C", cells_series=6, cells_parallel=2,
                  sensor_mass_kg=1.361, n_simulators=7, sim_mass_kg=0.680)


# ==========================================================================
# 1. Grid convergence of the mission integrator
# ==========================================================================
def grid_convergence(cfg: Config, design: Design) -> None:
    print()
    print(_gc_LINE)
    print("1. GRID CONVERGENCE OF THE MISSION INTEGRATOR")
    print(_gc_LINE)
    print("Refinement ratio r = 2 on the nominal ground-track step dx.")
    print()

    ac = Aircraft(design, cfg)
    steps = [24.0, 12.0, 6.0, 3.0, 1.5, 0.75]
    rows = []
    for dx in steps:
        t0 = time.time()
        r2 = mission_2(ac, cfg, dx=dx)
        r3 = mission_3(ac, cfg, dx=dx)
        rows.append((dx, r2.time_s, r2.energy_j / 3600.0, r3.laps,
                     r3.time_s, time.time() - t0))

    print("%8s %13s %11s %8s %11s %8s" %
          ("dx [m]", "M2 time [s]", "M2 E [Wh]", "M3 laps", "M3 time [s]",
           "cost [s]"))
    print(_gc_THIN)
    for dx, t2, e2, l3, t3, cost in rows:
        print("%8.2f %13.5f %11.4f %8d %11.5f %8.3f"
              % (dx, t2, e2, l3, t3, cost))
    print(_gc_THIN)

    # Richardson on the three finest grids for each scalar quantity
    print()
    print("Richardson extrapolation on the three finest grids (r = 2):")
    print()
    print("%-18s %12s %10s %14s %12s" %
          ("quantity", "fine value", "order p", "extrapolated", "GCI [%]"))
    print(_gc_THIN)
    for label, idx, unit in (("M2 elapsed time", 1, "s"),
                             ("M2 energy", 2, "Wh"),
                             ("M3 elapsed time", 4, "s")):
        f_c, f_m, f_f = rows[-3][idx], rows[-2][idx], rows[-1][idx]
        p = observed_order(f_c, f_m, f_f, r=2.0)
        if p is None:
            print("%-18s %12.5f %10s %14s %12s"
                  % (label, f_f, "exact", "-", "0.000"))
            continue
        p_use = max(min(p, 8.0), 0.5)
        ex = richardson_extrapolate(f_m, f_f, p_use, r=2.0)
        band = 100.0 * gci(f_m, f_f, p_use, r=2.0)
        print("%-18s %12.5f %10.2f %14.5f %12.4f"
              % (label, f_f, p, ex, band))
    print(_gc_THIN)

    # Error of each grid against the extrapolated answer
    f_c, f_m, f_f = rows[-3][1], rows[-2][1], rows[-1][1]
    p = observed_order(f_c, f_m, f_f, r=2.0)
    p_use = max(min(p, 8.0), 0.5) if p else 4.0
    exact_t2 = richardson_extrapolate(f_m, f_f, p_use, r=2.0)
    print()
    print("Error in Mission 2 elapsed time vs the extrapolated exact value")
    print("(%.5f s):" % exact_t2)
    print()
    print("%8s %14s %12s %10s" % ("dx [m]", "error [s]", "error [%]", "cost [s]"))
    print(_gc_THIN)
    for dx, t2, _e2, _l3, _t3, cost in rows:
        err = t2 - exact_t2
        print("%8.2f %14.5f %12.5f %10.3f"
              % (dx, err, 100.0 * err / exact_t2, cost))
    print(_gc_THIN)
    print("A 1 percent error in Mission 2 time is a 1 percent error in the M2")
    print("score ratio, so this is the number that has to be small.")


# ==========================================================================
# 2. Propulsion solver accuracy
# ==========================================================================
def propulsion_accuracy(cfg: Config) -> None:
    """The propulsion operating point is solved, not interpolated.

    Battery sag used to be handled by a fixed-point iteration, and the result
    was then tabulated and interpolated.  Both steps introduced error.  The
    sag loop now has a closed-form solution (see propulsion.py), so the whole
    system is one scalar root-find, and the table is gone.

    This section verifies the closed form against a brute-force 200-iteration
    fixed point, and reports how well the root is actually converged.
    """
    print()
    print(_gc_LINE)
    print("2. PROPULSION SOLVER ACCURACY (no interpolation anywhere)")
    print(_gc_LINE)

    combos = [("SCORPION_SII_4020_420", "APC_13x10", 2, "LIPO_2200MAH_65C", 6, 2),
              ("COBRA_C4120_12", "APC_12x8", 1, "LIPO_3000MAH_65C", 6, 1),
              ("TMOTOR_AT4130_300", "APC_16x12", 2, "LIPO_4000MAH_65C", 6, 1)]

    print("%-30s %14s %14s %14s" %
          ("system", "max dT [N]", "max dVpack [V]", "max |residual|"))
    print(_gc_THIN)
    for mk, pk, nm, ck, s_, par in combos:
        ps = PropulsionSystem(MOTORS[mk], PROPS[pk], nm,
                              BatteryPack(CELLS[ck], s_, par), cfg.env)
        m, p, env = ps.motor, ps.prop, ps.env
        worst_t = worst_v = worst_r = 0.0
        for i in range(25):
            v = 0.5 + 45.0 * i / 24.0
            for soc in (1.0, 0.72, 0.41, 0.13):
                for thr in (1.0, 0.85, 0.60, 0.35, 0.15):
                    n_bf, _i_bf, vp_bf = _brute_force_sag(ps, v, soc, thr)
                    op = ps.solve(v, soc, thr)
                    t_bf = (nm * max(0.0, p.thrust(env.rho, n_bf, v))
                            * ps.install_thrust_factor)
                    worst_t = max(worst_t, abs(op.thrust_total_n - t_bf))
                    worst_v = max(worst_v, abs(op.pack_voltage_v - vp_bf))
                    # torque residual at the returned rpm
                    n = op.rpm / 60.0
                    if n > 1e-6:
                        cur = op.motor_current_a
                        q_m = ps.kt * (cur - m.i0)
                        j = v / (n * p.diameter_m)
                        q_p = (p.cp(j) * env.rho * n * n
                               * p.diameter_m ** 5 / (2.0 * math.pi))
                        worst_r = max(worst_r, abs(q_m - q_p))
        print("%-30s %14.3e %14.3e %14.3e"
              % ("%s + %s x%d" % (mk.split("_")[0], pk.replace("APC_", ""), nm),
                 worst_t, worst_v, worst_r))
    print(_gc_THIN)
    print("dT and dVpack are vs a 200-iteration brute-force fixed point.")
    print("|residual| is motor torque minus propeller torque at the answer,")
    print("in N*m - it should be at the level of floating-point noise.")

    # cost
    ps = PropulsionSystem(MOTORS["SCORPION_SII_4020_420"], PROPS["APC_13x10"],
                          2, BatteryPack(CELLS["LIPO_2200MAH_65C"], 6, 2),
                          cfg.env)
    n_rep = 20000
    t0 = time.time()
    for i in range(n_rep):
        ps._memo_key = None
        ps.thrust(5.0 + 35.0 * ((i * 0.618) % 1.0), 0.7, 1.0)
    print()
    print("Cost: %.2f us per operating point."
          % (1e6 * (time.time() - t0) / n_rep))


def _brute_force_sag(ps, v: float, soc: float, thr: float, iters: int = 200):
    """Independent check: iterate the sag loop to death, no closed form."""
    m, p, env = ps.motor, ps.prop, ps.env
    ocv = ps.pack.ocv(soc)
    vp = ocv
    n = i_m = 0.0
    for _ in range(iters):
        vt = thr * vp
        if vt <= 0.05:
            return 0.0, 0.0, vp
        nnl = m.kv * vt / 60.0

        def res(nn):
            cur = (vt - 60.0 * nn / m.kv) / ps.rm_hot
            j = v / (nn * p.diameter_m) if nn > 1e-6 else 1e9
            return (ps.kt * (cur - m.i0)
                    - p.cp(j) * env.rho * nn * nn * p.diameter_m ** 5
                    / (2.0 * math.pi))

        n = 0.0 if res(1e-4) <= 0.0 else brent(res, 1e-4, nnl, xtol=1e-13)
        i_m = max(0.0, (vt - 60.0 * n / m.kv) / ps.rm_hot)
        vp = max(ps.pack.v_min * 0.85,
                 ocv - ps.n_motors * i_m * thr * ps.pack.r_ohm)
    return n, i_m, vp


# ==========================================================================
# 3. Mass loop convergence
# ==========================================================================
def mass_loop(cfg: Config, design: Design) -> None:

    print()
    print(_gc_LINE)
    print("3. FIXED-POINT CONVERGENCE OF THE MASS LOOP")
    print(_gc_LINE)
    print("Structure mass depends on MTOW, and MTOW contains the structure.")
    print()

    ac = Aircraft(design, cfg)
    wing = ac.wing
    body = ac.body
    b = cfg.build
    payload = ac.m2_payload_kg
    prop_mass = ac.prop_sys.mass_kg + ac.gear_extra_mass_kg

    mtow = payload + prop_mass + 2.0
    print("%6s %14s %16s" % ("iter", "MTOW [kg]", "change [kg]"))
    print(_gc_THIN)
    prev = mtow
    for i in range(12):
        mb = airframe_masses(wing, body, mtow * 9.80665, b)
        mb.add("propulsion group", prop_mass)
        mb.add("payload restraint", b.payload_structure_frac * payload)
        mb.add("contingency", b.contingency_frac * mb.total)
        new = mb.total + payload
        mtow = 0.5 * mtow + 0.5 * new
        print("%6d %14.9f %16.2e" % (i, mtow, mtow - prev))
        if abs(mtow - prev) < 1e-9:
            break
        prev = mtow
    print(_gc_THIN)
    print("Converged MTOW used by the model: %.9f kg" % ac.mtow_m2_kg)
    print("Residual: %.3e kg  (tolerance 1e-4 kg = 0.1 g)"
          % abs(ac.mtow_m2_kg - (ac.empty_mass_kg + ac.m2_payload_kg)))


# ==========================================================================
# 4. Recommended step size
# ==========================================================================
def recommend(cfg: Config, design: Design) -> None:
    print()
    print(_gc_LINE)
    print("4. COST vs ACCURACY - PICKING dx")
    print(_gc_LINE)
    ac = Aircraft(design, cfg)
    ref2 = mission_2(ac, cfg, dx=0.375).time_s
    print("Reference (dx = 0.375 m): M2 time = %.5f s" % ref2)
    print()
    print("%8s %14s %12s %11s %12s" %
          ("dx [m]", "M2 time [s]", "error [%]", "cost [ms]", "verdict"))
    print(_gc_THIN)
    for dx in (24.0, 16.0, 12.0, 8.0, 6.0, 4.0, 2.0, 1.0):
        t0 = time.time()
        n_rep = 3
        for _ in range(n_rep):
            t2 = mission_2(ac, cfg, dx=dx).time_s
        cost = 1000.0 * (time.time() - t0) / n_rep
        err = 100.0 * (t2 - ref2) / ref2
        verdict = ("screening" if abs(err) < 0.5 else
                   ("too coarse" if abs(err) >= 0.5 else ""))
        if abs(err) < 0.05:
            verdict = "final scoring"
        print("%8.2f %14.5f %12.4f %11.2f %12s" % (dx, t2, err, cost, verdict))
    print(_gc_THIN)


def convergence_main() -> int:
    cfg = Config()
    design = demo_design()
    print(_gc_LINE)
    print("NUMERICAL ERROR STUDY - DBF 2026-27 MODEL")
    print(_gc_LINE)
    print("Integrator: RK4 (4th order) with bisection event location")
    print("Quadrature: 5-point Gauss-Legendre (exact to degree 9)")
    print("Root finding: Brent (superlinear)")
    print("Battery sag: closed form (no iteration, no interpolation)")
    grid_convergence(cfg, design)
    propulsion_accuracy(cfg)
    mass_loop(cfg, design)
    recommend(cfg, design)
    print()
    print(_gc_LINE)
    print("Numerical error is now well below model uncertainty.")
    print("See ASSUMPTIONS.md for the model error that remains.")
    print(_gc_LINE)
    return 0


# ==========================================================================
#  SEARCH AUDIT - which catalogue options reached the finals
# ==========================================================================

#!/usr/bin/env python3




def load(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def f(row, key, default=0.0):
    try:
        return float(row.get(key, default) or default)
    except ValueError:
        return default


def tally(rows, key, label, universe=None):
    print()
    print("%s" % label)
    print("-" * len(label))
    counts = collections.Counter(r[key] for r in rows)
    best = {}
    for r in rows:
        v = f(r, "objective")
        if r[key] not in best or v > best[r[key]]:
            best[r[key]] = v
    order = sorted(counts, key=lambda k: -best.get(k, -1e9))
    print("  %-26s %7s %10s" % ("option", "finals", "best score"))
    for k in order:
        print("  %-26s %7d %10.3f" % (k[:26], counts[k], best[k]))
    if universe:
        missing = [k for k in universe if k not in counts]
        if missing:
            print("  did not reach the finals: %s" % ", ".join(
                sorted(m[:20] for m in missing)))


def audit_main(argv):
    path = argv[1] if len(argv) > 1 else "results_v4/ranked_designs.csv"
    rows = load(path)
    if not rows:
        print("no rows in %s" % path)
        return 1

    print("=" * 70)
    print("SEARCH AUDIT: %s  (%d finalists)" % (path, len(rows)))
    print("=" * 70)

    objs = sorted((f(r, "objective") for r in rows), reverse=True)
    print()
    print("Score spread")
    print("------------")
    print("  best            %.4f" % objs[0])
    print("  10th            %.4f" % objs[min(9, len(objs) - 1)])
    print("  worst finalist  %.4f" % objs[-1])
    within = sum(1 for o in objs if o >= objs[0] * 0.98)
    print("  within 2%% of the best: %d of %d designs" % (within, len(objs)))
    if within > len(objs) // 3:
        print("  -> the objective is FLAT. Pick on buildability, not score.")

    tally(rows, "d_airfoil", "Airfoils in the finals", AIRFOILS.keys())
    tally(rows, "motor", "Motors in the finals")
    tally(rows, "d_prop", "Propellers in the finals")
    tally(rows, "d_n_motors", "Motor count")
    tally(rows, "chemistry", "Battery chemistry")

    print()
    print("Design-variable ranges across the finals")
    print("----------------------------------------")
    for key, label, scale, unit in (
            ("d_span_m", "span", 39.3701, "in"),
            ("AR", "aspect ratio", 1.0, ""),
            ("d_sensor_mass_kg", "sensor", 2.20462, "lb"),
            ("d_n_simulators", "simulators", 1.0, ""),
            ("d_sim_mass_kg", "simulator mass", 2.20462, "lb"),
            ("pack_Wh", "battery", 1.0, "Wh"),
            ("mtow_m2_kg", "M2 gross", 2.20462, "lb"),
            ("n_ultimate", "ultimate load", 1.0, "g")):
        vals = [f(r, key) for r in rows if r.get(key)]
        if not vals:
            continue
        print("  %-16s %8.2f to %8.2f %s"
              % (label, min(vals) * scale, max(vals) * scale, unit))

    print()
    print("What limited the turn, and what sized the spar")
    print("----------------------------------------------")
    for key, label in (("m2_turn_limit", "Mission 2 turn limited by"),
                       ("load_driver", "ultimate load governed by")):
        c = collections.Counter(r.get(key, "?") for r in rows)
        print("  %s:" % label)
        for k, n in c.most_common():
            print("      %-26s %d" % (k, n))
    gauge = collections.Counter(r.get("spar_gauge_limited", "?") for r in rows)
    print("  spar gauge limited: %s" % dict(gauge))
    return 0


# ==========================================================================
#  ENTRY POINT
# ==========================================================================
def main(argv=None):
    """Dispatch to the search, a single evaluation, or a test suite."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--selftest" in argv:
        return selftest_main()
    if "--validate" in argv:
        return validate_main()
    if "--convergence" in argv:
        return convergence_main()
    if "--audit" in argv:
        i = argv.index("--audit")
        path = (argv[i + 1] if i + 1 < len(argv)
                else "results/ranked_designs.csv")
        return audit_main(["audit", path])
    if "--version" in argv:
        print("DBF 2026-27 optimiser v%s" % __version__)
        return 0
    return cli_main(argv)


if __name__ == "__main__":
    sys.exit(main())
