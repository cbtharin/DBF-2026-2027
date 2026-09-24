"""Dependency-free numerical methods.

The whole toolchain is pure-stdlib on purpose: it has to run on whatever laptop
is open in the shop at 2am, without a working `pip install`.

Accuracy notes
--------------
Every routine here states its order of accuracy, because the mission model
integrates an ODE thousands of times and first-order error shows up directly
in lap time:

    brent               root finding, superlinear, machine-precision capable
    gauss_legendre_5    quadrature, exact for polynomials up to degree 9
    rk4_step            ODE, local error O(h^5), global O(h^4)
    PCHIP               monotone cubic interpolation, O(h^3), no overshoot
    scan_then_refine    bracket a maximum by scan, then golden section
    golden_max          unimodal maximisation, no derivatives
    richardson_*        grid-convergence estimates (observed order, GCI)
    solve_linear        dense Gaussian elimination, partial pivoting
"""
from __future__ import annotations

import math
import random
from typing import Callable, List, Optional, Sequence, Tuple


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
