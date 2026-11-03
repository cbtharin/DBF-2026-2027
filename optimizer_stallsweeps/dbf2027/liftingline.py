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

(b the full span).  The induced angle is then
alpha_i = 2 sum_n n A_n sin(n theta) / sin(theta), and the section law
Gamma = 1/2 V c a_0 (alpha - alpha_L0 - alpha_i) collocated at each station
gives the monoplane equation

    sum_n  A_n sin(n theta) [ 2 n mu + sin(theta) ]  =  mu (alpha - alpha_L0) sin(theta)

with the local geometry parameter

    mu(theta) = c(theta) a_0 / (8 b)

and `a_0` the section lift-curve slope (2 pi per radian).  (The textbook
form - Anderson, Fundamentals of Aerodynamics 5.3.2 - writes
Gamma = 2 b V sum A'_n sin(n theta), mu' = c a_0 / (4 b) and the row
[n mu' + sin theta]: the same equation, A'_n = 2 A_n.  The solve reproduces
the textbook values: a rectangular wing of AR 6 has e = 0.954, delta ~ 0.05,
and a lift slope of 4.53 /rad, Helmbold's 4.53.)  For an untwisted
wing the right-hand side is a constant times sin(theta), so the solution
scales linearly with `(alpha - alpha_L0)` and the *shape* of the loading -
which is all span efficiency depends on - is independent of angle of attack.

Only odd harmonics survive for a symmetric wing.  From the solution:

    CL = 2 pi AR A_1
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


def clear_caches() -> None:
    """Empty this module's shared tables.  They are keyed on rounded
    inputs but hold what the first caller's exact inputs gave, so every
    aircraft starts with them empty (Aircraft.__init__): a design's score
    must not depend on what was scored before it in the same process."""
    for table in (_CACHE, _STALL_CACHE, _BLOWN_CACHE, _TERMS_CACHE):
        table.clear()


def span_efficiency(aspect_ratio: float, taper: float,
                    n_terms: int = 10, a0: float = 2.0 * math.pi) -> float:
    """Inviscid span efficiency `e` for an untwisted tapered wing."""
    ar = max(min(aspect_ratio, 40.0), 1.0)
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
        row = [math.sin(n * th) * (2.0 * n * mu + math.sin(th)) for n in odd]
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
    """Finite-wing lift-curve slope per radian, unswept, low speed - the
    DATCOM (Polhamus) form 2 pi AR / (2 + sqrt((AR/k)^2 + 4)), k = a0 / 2 pi
    the section's slope over the thin-airfoil value.  It tends to a0 as the
    aspect ratio grows and to the slender-wing pi AR / 2 as it shrinks.
    (It was a0 AR / (2 + sqrt(AR^2 + 4)): the same for a0 = 2 pi, but for
    any other section slope wrong at long aspect ratio, m109.)"""
    ar = max(aspect_ratio, 0.5)
    k = max(a0, 1e-6) / (2.0 * math.pi)
    return 2.0 * math.pi * ar / (2.0 + math.sqrt((ar / k) ** 2 + 4.0))


def downwash_gradient(aspect_ratio: float, a0: float = 2.0 * math.pi) -> float:
    """d(epsilon)/d(alpha) at the tail, the standard low-aspect-ratio
    approximation `2 a_w / (pi AR)`.

    This is why a tail is less effective than its area suggests: the wing
    turns the flow down before the tail ever sees it.
    """
    aw = lift_curve_slope(aspect_ratio, a0)
    return min(0.85, 2.0 * aw / (math.pi * max(aspect_ratio, 0.5)))


# ==========================================================================
# Where the wing stalls first
# ==========================================================================
_STALL_CACHE: Dict[Tuple[int, int, int], Tuple[float, float]] = {}


def spanwise_loading(aspect_ratio: float, taper: float, n_terms: int = 24,
                     a0: float = 2.0 * math.pi):
    """Local lift coefficient along the span, as a ratio to the wing CL.

    Returns [(eta, cl/CL, c/c_mean), ...] from root (eta 0) to tip, for an
    untwisted straight-tapered wing - the same monoplane-equation solve as
    span_efficiency.  With Gamma = 4 b V sum A_n sin(n theta), the local
    section cl = 2 Gamma / (V c) and the wing CL = 2 pi AR A_1.
    """
    ar = max(min(aspect_ratio, 40.0), 1.0)
    lam = max(min(taper, 1.0), 0.05)
    b = 1.0
    c_mean = b / ar
    c_root = 2.0 * c_mean / (1.0 + lam)
    thetas = [(i + 1) * math.pi / (2.0 * n_terms) for i in range(n_terms)]
    odd = [2 * k + 1 for k in range(n_terms)]
    mat, rhs = [], []
    for th in thetas:
        eta = abs(math.cos(th))
        c = c_root * (1.0 - (1.0 - lam) * eta)
        mu = c * a0 / (8.0 * b)
        mat.append([math.sin(n * th) * (2.0 * n * mu + math.sin(th)) for n in odd])
        rhs.append(mu * math.sin(th))
    coeffs = solve_linear(mat, rhs)
    cl_wing = 2.0 * math.pi * ar * coeffs[0]
    out = []
    for th in thetas:
        eta = abs(math.cos(th))
        c = c_root * (1.0 - (1.0 - lam) * eta)
        gam = 4.0 * b * sum(coeffs[k] * math.sin(n * th)
                            for k, n in enumerate(odd))
        out.append((eta, 2.0 * gam / c / cl_wing, c / c_mean))
    out.sort()
    return out


def stall_factor(aspect_ratio: float, taper: float,
                 re_exponent: float) -> Tuple[float, float]:
    """(factor, eta) - the wing reaches its maximum lift when its FIRST
    section reaches the section's own CLmax, not when the average does.

    Each section works at cl = (cl/CL)(eta) x CL and can reach
    cl_max_2d x (Re_local / Re_ref)^k, and its Reynolds number scales with
    its own chord.  So the wing CLmax is the airfoil CLmax at the mean
    chord's Reynolds number times

        F = min over eta of  (c/c_mean)^k / (cl/CL)

    and eta is where that minimum sits - where stall begins.  A rectangular
    wing stalls at the root, where its loading peaks (F = 0.88 at AR 6.5);
    moderate taper evens the loading out (0.92 at taper 0.45, stalling at
    eta 0.66), but sharp taper moves the peak outboard onto smaller,
    lower-Reynolds chords, so the wing stalls near the tip first and F
    falls (0.83-0.85 at taper 0.25, AR 6.5, for k 0.21-0.15).  That is the
    tip-stall penalty that makes low taper cost something.
    """
    key = (int(round(aspect_ratio * 100)), int(round(taper * 1000)),
           int(round(re_exponent * 1000)))
    hit = _STALL_CACHE.get(key)
    if hit is not None:
        return hit
    best = (float("inf"), 0.0)
    for eta, ratio, c_rel in spanwise_loading(aspect_ratio, taper):
        if eta > 0.98:              # the tip vortex core, not a section
            continue
        f = (c_rel ** re_exponent) / max(ratio, 1e-9)
        if f < best[0]:
            best = (f, eta)
    _STALL_CACHE[key] = best
    return best


# ==========================================================================
# A wing in a propeller slipstream
# ==========================================================================
_BLOWN_CACHE: Dict[Tuple, Tuple[float, float]] = {}


def jet_speed_ratio(eta: float, jets) -> float:
    """Airspeed / free stream at half-span station eta, inside any of the
    `jets` - each (eta_centre, eta_half_width, u[, eta_edge]).  A
    slipstream's edge is a mixing layer, not a step: the speed eases from
    u to 1 over a raised cosine either side of the nominal edge, which
    keeps the jet's integrated speed increase.  The layer is a quarter of
    the half-width, or `eta_edge` when given - a strip cut from a wider
    slipstream (beside the fuselage) keeps the whole slipstream's mixing
    layer, not one scaled to the strip.  (A hard step is also what a
    Fourier series resolves worst.)"""
    out = 1.0
    for jet in jets:
        centre, half, u = jet[0], jet[1], jet[2]
        eps = min(jet[3], half) if len(jet) > 3 else 0.25 * half
        d = abs(eta - centre)
        if d <= half - eps:
            s = 1.0
        elif d >= half + eps:
            s = 0.0
        else:
            s = 0.5 * (1.0 + math.cos(math.pi * (d - (half - eps))
                                      / (2.0 * eps)))
        out = max(out, 1.0 + (u - 1.0) * s)
    return out


def jet_swirl(eta: float, jets) -> float:
    """The swirl pattern at half-span station eta: each jet's 5th element
    is +1 where its swirl is upwash (the inboard side of an inboard-up
    propeller), -1 where it is downwash, 0 for none; eased at the edges as
    the speed is (jet_speed_ratio)."""
    out = 0.0
    for jet in jets:
        if len(jet) < 5 or not jet[4]:
            continue
        centre, half = jet[0], jet[1]
        eps = min(jet[3], half)
        d = abs(eta - centre)
        if d <= half - eps:
            s = 1.0
        elif d >= half + eps:
            s = 0.0
        else:
            s = 0.5 * (1.0 + math.cos(math.pi * (d - (half - eps))
                                      / (2.0 * eps)))
        out += jet[4] * s
    return out


_TERMS_CACHE: Dict = {}


def blown_wing_terms(aspect_ratio: float, taper: float, jets,
                     n_terms: int = 24,
                     a0: float = 2.0 * math.pi) -> Tuple[float, ...]:
    """The blown wing's lift and induced drag, split so the angle of
    attack and the swirl angle can be put in afterwards.

    The swirl in a wing propeller's slipstream changes the angle the
    sections behind it see by +eps on one side, -eps on the other
    (jet_swirl); the lifting line is linear, so its solution is
    alpha A + eps P, A for unit angle and P for unit swirl.  Returns
    (CL_a, CL_p, S_aa, S_ap, S_pp, W_a, W_p) such that

        CL  = alpha CL_a + eps CL_p
        CDi = 4 pi AR (alpha^2 S_aa + 2 alpha eps S_ap + eps^2 S_pp)
              - eps (alpha W_a + eps W_p)

    The first term is the trailing vortices' drag (on the free stream,
    with these coefficients' normalisation); the second is swirl
    recovery - where the swirl is upwash the local lift tilts forward,
    and that forward component is thrust taken back from the swirl.
    """
    ar = max(min(aspect_ratio, 40.0), 1.0)
    lam = max(min(taper, 1.0), 0.05)
    key = (int(round(ar * 100)), int(round(lam * 1000)), n_terms,
           tuple(tuple(round(v, 4) for v in j) for j in jets))
    hit = _TERMS_CACHE.get(key)
    if hit is not None:
        return hit
    b = 1.0
    c_mean = b / ar
    c_root = 2.0 * c_mean / (1.0 + lam)
    thetas = [(i + 1) * math.pi / (2.0 * n_terms) for i in range(n_terms)]
    odd = [2 * k + 1 for k in range(n_terms)]
    mat, rhs_a, rhs_p = [], [], []
    for th in thetas:
        eta = abs(math.cos(th))
        c = c_root * (1.0 - (1.0 - lam) * eta)
        mu = c * a0 / (8.0 * b)
        u = jet_speed_ratio(eta, jets)
        p = jet_swirl(eta, jets)
        mat.append([math.sin(n * th) * (2.0 * n * mu + math.sin(th)) for n in odd])
        rhs_a.append(mu * u * math.sin(th))
        rhs_p.append(mu * u * p * math.sin(th))
    ca = solve_linear(mat, rhs_a)
    cp = (solve_linear(mat, rhs_p) if any(abs(x) > 0.0 for x in rhs_p)
          else [0.0] * n_terms)
    m = 6 * n_terms
    cla = clp = wa = wp = 0.0
    for i in range(m):
        th = (i + 0.5) * (0.5 * math.pi) / m
        eta = abs(math.cos(th))
        u = jet_speed_ratio(eta, jets)
        p = jet_swirl(eta, jets)
        sa = sum(ca[k] * math.sin(n * th) for k, n in enumerate(odd))
        sp = sum(cp[k] * math.sin(n * th) for k, n in enumerate(odd))
        cla += u * sa * math.sin(th)
        clp += u * sp * math.sin(th)
        wa += u * sa * math.sin(th) * p
        wp += u * sp * math.sin(th) * p
    scale = 4.0 * ar * 2.0 * (0.5 * math.pi) / m
    out = (scale * cla, scale * clp,
           sum(n * ca[k] * ca[k] for k, n in enumerate(odd)),
           sum(n * ca[k] * cp[k] for k, n in enumerate(odd)),
           sum(n * cp[k] * cp[k] for k, n in enumerate(odd)),
           scale * wa, scale * wp)
    _TERMS_CACHE[key] = out
    return out


def blown_wing(aspect_ratio: float, taper: float, jets,
               re_exponent: float, n_terms: int = 24,
               a0: float = 2.0 * math.pi,
               eta_body: float = 0.0) -> Tuple[float, float]:
    """(CL per unit angle, stall factor) of a wing with part of its span in
    faster air - a propeller slipstream.

    `jets` is a sequence of (eta_centre, eta_half_width, u) on the half
    span (eta = |2y/b|): the air there arrives at u times the free stream
    (jet_speed_ratio).

    The lifting line carries over with one change.  A section's lift is
    set by its OWN airspeed, Gamma = 1/2 V_loc c a0 (alpha - w/V_loc), and
    the trailing vorticity induces the same downwash w as before, so the
    monoplane equation becomes

        sum_n A_n sin(n theta) [2 n mu + sin theta] = mu u(theta) alpha sin theta

    - the right-hand side scaled by the local speed ratio.  The lift is
    then Kutta-Joukowski at the local speed,

        CL = 4 AR  integral_0^pi  u  sum_n A_n sin(n theta) sin theta  d theta

    (referenced to the free stream), and each section can reach its CLmax
    at its own Reynolds number, which the faster air raises by u too, so
    the stall factor becomes

        F = min over eta of  (u c / c_mean)^k / (cl_local / CL)

    The stall search skips the span the fuselage covers (eta below
    `eta_body`): a wing section sitting on the body is not where a wing-
    body combination stalls first, and with a nose propeller's slipstream
    blowing the exposed root right beside it, its edge upwash would
    otherwise make that covered strip look like the first to stall.

    What this gets right that a strip estimate (lift x u^2 over the blown
    area) does not: a jet only a propeller wide cannot carry its full
    dynamic-pressure lift, because the extra circulation it sheds at its
    edges induces a downwash on it - about 15% of the strip estimate's
    gain at u = 1.6 across a nose tractor's slipstream.  What it leaves
    out: swirl (the rotating slipstream raises the angle of attack on one
    side and lowers it on the other) and the jet's deflection by the wing.
    """
    ar = max(min(aspect_ratio, 40.0), 1.0)
    lam = max(min(taper, 1.0), 0.05)
    key = (int(round(ar * 100)), int(round(lam * 1000)),
           int(round(re_exponent * 1000)), n_terms,
           int(round(eta_body * 10000)),
           tuple(tuple(round(v, 4) for v in j) for j in jets))
    hit = _BLOWN_CACHE.get(key)
    if hit is not None:
        return hit

    b = 1.0
    c_mean = b / ar
    c_root = 2.0 * c_mean / (1.0 + lam)
    thetas = [(i + 1) * math.pi / (2.0 * n_terms) for i in range(n_terms)]
    odd = [2 * k + 1 for k in range(n_terms)]
    mat, rhs, st = [], [], []
    for th in thetas:
        eta = abs(math.cos(th))
        c = c_root * (1.0 - (1.0 - lam) * eta)
        mu = c * a0 / (8.0 * b)
        u = jet_speed_ratio(eta, jets)
        mat.append([math.sin(n * th) * (2.0 * n * mu + math.sin(th)) for n in odd])
        rhs.append(mu * u * math.sin(th))
        st.append((th, eta, c, u))
    coeffs = solve_linear(mat, rhs)

    def series(th: float) -> float:
        return sum(coeffs[k] * math.sin(n * th) for k, n in enumerate(odd))

    # the lift integral: midpoint rule on the half-span (symmetric wing)
    m = 6 * n_terms
    acc = 0.0
    for i in range(m):
        th = (i + 0.5) * (0.5 * math.pi) / m
        acc += (jet_speed_ratio(abs(math.cos(th)), jets) * series(th)
                * math.sin(th))
    cl_wing = 4.0 * ar * 2.0 * acc * (0.5 * math.pi) / m
    # where stall starts: the series is continuous, so search a fine grid
    # rather than the collocation points alone
    f_best = float("inf")
    for i in range(m):
        th = (i + 0.5) * (0.5 * math.pi) / m
        eta = abs(math.cos(th))
        if eta > 0.98:              # the tip vortex core, not a section
            continue
        if eta < eta_body:          # on the fuselage, not free wing
            continue
        c = c_root * (1.0 - (1.0 - lam) * eta)
        u = jet_speed_ratio(eta, jets)
        cl_local = 2.0 * 4.0 * b * series(th) / (u * c)
        ratio = cl_local / max(cl_wing, 1e-12)
        f_best = min(f_best, (u * c / c_mean) ** re_exponent
                     / max(ratio, 1e-9))
    out = (cl_wing, f_best)
    _BLOWN_CACHE[key] = out
    return out
