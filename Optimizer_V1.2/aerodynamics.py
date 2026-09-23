"""Aerodynamic model: L, D and the drag polar.

Level of fidelity: component parasite-drag buildup + airfoil section polar +
lifting-line induced drag + trim.  Right for concept selection; validate the
winner in XFLR5/AVL before committing.

Giving every airfoil a fair hearing
-----------------------------------
Two effects are included specifically so that high-lift sections are judged
on what they actually deliver rather than on their catalogue CLmax:

1.  **Reynolds-dependent CLmax.**  A DBF wing runs at Re 150k-500k, well
    below where section data is usually quoted, and high-camber sections
    lose more at low Re than thin ones do.  So the Re penalty scales with
    the section's own camber.  Quoting S1223 at CLmax 2.2 on a 9 inch chord
    is exactly the kind of optimism this model is supposed to avoid.

2.  **Trim drag, with the empennage in the lift balance.**  The horizontal
    tail is not a passive drag item: it carries part of the load.  Its
    share comes from moments about the CG and has a camber term (from
    Cm_ac) and a CG term (from where the centre of gravity sits relative
    to the wing AC).  The wing then carries L - L_t, which for a cambered
    section is *more* than the aircraft weight, and the tail makes induced
    drag of its own.  S1223 has Cm = -0.28 against MH32's -0.06; without
    this the high-lift sections win on lift and pay nothing for it.

    Because the CG term scales with the lift being carried, the tail load
    grows in a turn - which the earlier camber-only model missed.

Sign/unit conventions
---------------------
    L, D, T, W   newtons
    V            true airspeed, m/s
    S            reference (wing planform) area, m^2
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

from .components import Airfoil
from .config import BuildStandard, Environment
from .liftingline import (downwash_gradient, lift_curve_slope,
                          span_efficiency)
from .numerics import brent, scan_then_refine


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
