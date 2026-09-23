"""Prandtl lifting-line: span efficiency from the actual planform.

Replaces Raymer's curve fit `e = 1.78(1 - 0.045 AR^0.68) - 0.64`, which is a
regression over historical aircraft and knows nothing about taper.  This
solves the monoplane equation instead, so span efficiency comes out of the
geometry and taper becomes a design variable the optimiser can actually use.

Derivation
----------
Represent the circulation as a Fourier sine series in the spanwise angle
`theta`, where `y = -(b/2) cos(theta)`:

    Gamma(theta) = 4 b V  sum_n  A_n sin(n theta)

Substituting into the lifting-line integro-differential equation and
collocating gives the monoplane equation at each station:

    sum_n  A_n sin(n theta) [ n mu + sin(theta) ]  =  mu (alpha - alpha_L0) sin(theta)

with the local geometry parameter

    mu(theta) = c(theta) a_0 / (8 b)

and `a_0` the section lift-curve slope (2 pi per radian).  For an untwisted
wing the right-hand side is a constant times sin(theta), so the solution
scales linearly with `(alpha - alpha_L0)` and the *shape* of the loading -
which is all span efficiency depends on - is independent of angle of attack.

Only odd harmonics survive for a symmetric wing.  From the solution:

    CL = pi AR A_1
    delta = sum_{n = 3,5,...}  n (A_n / A_1)^2
    e = 1 / (1 + delta)

An elliptical wing gives A_3 = A_5 = ... = 0, hence delta = 0 and e = 1,
which is the standard check the self-test makes.
"""
from __future__ import annotations

import math
from typing import Dict, Tuple

from .numerics import solve_linear

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
