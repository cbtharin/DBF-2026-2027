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
from dataclasses import dataclass, replace
from typing import List, Optional

from .components import Airfoil
from .config import BuildStandard, Environment
from .liftingline import (blown_wing, blown_wing_terms, downwash_gradient, lift_curve_slope,
                          span_efficiency, stall_factor)
from .numerics import brent, scan_then_refine


@dataclass
class WingGeometry:
    span_m: float
    # Mean chord, S / b - the lifting line's c_mean, so the area is b c.
    chord_m: float
    airfoil: Airfoil
    taper: float = 1.0
    # height of the vertical winglet at each tip, m (0 = none, m110)
    winglet_h_m: float = 0.0

    @property
    def area_m2(self) -> float:
        return self.span_m * self.chord_m

    @property
    def winglet_area_m2(self) -> float:
        """Planform of ONE winglet: root chord the wing's tip chord, taper
        0.5."""
        return max(self.winglet_h_m, 0.0) * 0.75 * self.tip_chord_m

    @property
    def aspect_ratio(self) -> float:
        return self.span_m / self.chord_m

    @property
    def root_chord_m(self) -> float:
        """Root chord of the straight-tapered wing with this area:
        S = b c_r (1 + lam) / 2, so c_r = 2 c / (1 + lam).

        This used the mean-AERODYNAMIC-chord formula, which drew a wing
        smaller than the one the model flew - 5% less area at taper 0.43,
        13% at 0.20 - while the lifting line (and the area) used the mean
        chord.  The printed root and tip chords are what gets built."""
        lam = min(max(self.taper, 0.0), 1.0)
        return 2.0 * self.chord_m / (1.0 + lam)

    @property
    def tip_chord_m(self) -> float:
        return self.root_chord_m * min(max(self.taper, 0.0), 1.0)

    @property
    def mac_m(self) -> float:
        """Mean aerodynamic chord, (2/3) c_r (1 + lam + lam^2)/(1 + lam)."""
        lam = min(max(self.taper, 0.0), 1.0)
        return (2.0 / 3.0 * self.root_chord_m * (1.0 + lam + lam * lam)
                / (1.0 + lam))

    @property
    def thickness_root_m(self) -> float:
        return self.airfoil.t_over_c * self.root_chord_m

    def chord_at(self, y: float) -> float:
        """Local chord `y` metres out from the centreline."""
        eta = min(abs(y) / max(0.5 * self.span_m, 1e-9), 1.0)
        lam = min(max(self.taper, 0.0), 1.0)
        return self.root_chord_m * (1.0 - (1.0 - lam) * eta)

    def planform_within_m2(self, width: float) -> float:
        """Planform inside a centred strip `width` wide (the chord is linear,
        so its mean over the strip is the chord a quarter-width out)."""
        w = min(max(width, 0.0), self.span_m)
        return w * self.chord_at(0.25 * w)


@dataclass
class BodyGeometry:
    """Fuselage + empennage.  The length is the nose (sized by the parts in
    it and the balance, layout.py) plus the payload bay plus a tailcone that
    tapers to the tail post; the tail arm is the lightest one that fits."""
    fuse_length_m: float
    fuse_width_m: float
    fuse_height_m: float
    tail_arm_m: float
    s_horiz_m2: float
    s_vert_m2: float
    nose_length_m: float = 0.0
    tailcone_length_m: float = 0.0
    tail_post_perimeter_m: float = 0.0
    tail_arm_optimum_m: float = 0.0     # lightest tails + tailcone
    tail_arm_min_m: float = 0.0         # shortest the tailcone taper allows
    # Full-section perimeter, rounded corners included (layout.
    # section_perimeters); 0 means a sharp-cornered box.
    section_perimeter_m: float = 0.0

    @property
    def fuse_wetted_m2(self) -> float:
        """Skin area of the whole shell - what the fuselage WEIGHS.  The
        drag uses this less whatever the wing covers (AeroModel)."""
        perim = (self.section_perimeter_m
                 or 2.0 * (self.fuse_width_m + self.fuse_height_m))
        l_tc = min(max(self.tailcone_length_m, 0.0), self.fuse_length_m)
        # constant section to the end of the payload bay, then a linear
        # taper to the tail post
        return (perim * (self.fuse_length_m - l_tc)
                + 0.5 * (perim + self.tail_post_perimeter_m) * l_tc)

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


def junction_drag_area_m2(t_over_c: float, chord_m: float) -> float:
    """Interference drag of one wing (or tail) root meeting a fuselage
    wall, as CD*A in m^2, at zero lift.  Hoerner, Fluid-Dynamic Drag
    ch. 8, on the section thickness squared:

        D_int / (q t^2) = 0.75 (t/c) - 0.0003 / (t/c)^2

    Thin sections come out at or below zero (the junction behaves like a
    fillet); that is floored at zero rather than counted as a thrust.
    """
    tc = max(t_over_c, 1e-3)
    t = tc * max(chord_m, 0.0)
    return max(0.0, 0.75 * tc - 0.0003 / (tc * tc)) * t * t


def tail_form_factor(t_over_c: float) -> float:
    """Raymer's thickness form factor for a tail surface with its maximum
    thickness at 30% chord, at low Mach: 1 + (0.6/0.3) t/c + 100 (t/c)^4."""
    return 1.0 + 2.0 * t_over_c + 100.0 * t_over_c ** 4


def skin_friction_cf(re: float, laminar_fraction: float = 0.2) -> float:
    re = max(re, 5.0e4)
    cf_turb = 0.074 / re ** 0.2
    cf_lam = 1.328 / math.sqrt(re)
    return laminar_fraction * cf_lam + (1.0 - laminar_fraction) * cf_turb


def _hough_ordway_mean(hub: float, n: int = 2000) -> float:
    """Disc-area mean of r* sqrt(1 - r*), r* = (x - hub) / (1 - hub), so
    that the shape divided by it averages 1 over the disc."""
    acc = 0.0
    for k in range(n):
        x = (k + 0.5) / n
        if x > hub:
            rs = (x - hub) / (1.0 - hub)
            acc += rs * math.sqrt(max(1.0 - rs, 0.0)) * 2.0 * x / n
    return max(acc, 1e-12)


class AeroModel:
    def __init__(self, wing: WingGeometry, body: BodyGeometry,
                 env: Environment, build: BuildStandard,
                 n_motors: int = 1,
                 extra_drag_area_m2: float = 0.0,
                 prop_diameter_m: float = 0.0,
                 nacelle_width_m: float = 0.0,
                 wing_prop_layout: str = "mid",
                 props_in_span: bool = False):
        # each wing nacelle's drag area: measured if given, else its
        # coefficient on the frontal area of the can and cowl
        self.nacelle_drag_area_m2 = (
            build.nacelle_drag_area_m2
            if build.nacelle_drag_area_m2 is not None
            else build.nacelle_cd * 0.25 * math.pi * nacelle_width_m ** 2)
        self.wing = wing
        self.body = body
        self.env = env
        self.build = build
        self.n_motors = n_motors
        self.extra_drag_area_m2 = extra_drag_area_m2

        self.S = wing.area_m2
        self.AR = wing.aspect_ratio
        # motors in wing nacelles: all of them but the nose one of an odd
        # count (layout.motor_stations)
        self.n_wing_motors = max(n_motors, 1) - max(n_motors, 1) % 2

        # ---- the fuselage on the wing (computed) --------------------------
        # The fuselage covers a strip of the planform its own width wide.
        self.fuse_d_over_b = body.fuse_width_m / max(wing.span_m, 1e-9)
        self.covered_frac = (wing.planform_within_m2(body.fuse_width_m)
                             / max(wing.planform_within_m2(wing.span_m),
                                   1e-12))
        self.s_covered = self.covered_frac * self.S
        # Span efficiency: the lift distribution dips where the fuselage
        # interrupts the span.  Kroo (Aircraft Design: Synthesis and
        # Analysis): e = e_wing (1 - 2 (d/b)^2).
        self.e_wing = oswald_efficiency(self.AR, wing.taper)
        self.fuselage_e_factor = max(0.5, 1.0 - 2.0 * self.fuse_d_over_b ** 2)
        self.e = self.e_wing * self.fuselage_e_factor
        # Winglets raise the effective aspect ratio (m110): Hoerner's
        # end-plate relation AR (1 + k h/b).  It sets the induced drag, the
        # lift slope and the downwash; CLmax stays the planform's.
        self.winglet_ar_factor = 1.0 + (build.winglet_ar_gain
                                        * max(wing.winglet_h_m, 0.0)
                                        / max(wing.span_m, 1e-9))
        self.AR_eff = self.AR * self.winglet_ar_factor
        self.k_induced = 1.0 / (math.pi * self.AR_eff * self.e)
        # CLmax: the share of the covered planform's lift that survives at
        # the stall (BuildStandard.fuselage_lift_carryover).
        carry = min(max(build.fuselage_lift_carryover, 0.0), 1.0)
        self.fuselage_clmax_factor = 1.0 - (1.0 - carry) * self.covered_frac
        # Skin the wing loses to the fuselage: one surface of the covered
        # strip for a top wing, both for a mid wing.
        n_cov = min(max(build.wing_covered_surfaces, 0), 2)
        self.wing_wetted_frac = max(0.0, 1.0 - 0.5 * n_cov * self.covered_frac)
        # ...and the fuselage top (or bottom) the wing sits on is not wetted
        # either.  A mid wing only takes the root sections' footprint out of
        # the fuselage sides, which is left in.
        self.fuse_covered_m2 = self.s_covered if n_cov == 1 else 0.0

        af = wing.airfoil
        # the section's own lift slope (XFOIL, Components._LIFT_XFOIL_500K;
        # 2 pi for a hand-entered section), m109
        self.a0_section = af.cl_alpha_per_rad or 2.0 * math.pi
        self.a_wing = lift_curve_slope(self.AR_eff, self.a0_section)
        self.downwash = downwash_gradient(self.AR_eff, self.a0_section)
        self.eta_tail = build.tail_eta
        # Camber-sensitive Reynolds exponent: a heavily cambered high-lift
        # section loses more of its CLmax at low Re than a thin one.
        self._re_exp = 0.15 + 0.10 * max(0.0, af.cl_max_2d - 1.2)
        # Where the wing stalls first (liftingline.stall_factor): the
        # planform's peak local loading against each section's own,
        # chord-dependent Reynolds number.  This is what prices taper - a
        # sharply tapered wing tip-stalls - and it IS the finite-span
        # correction: a separate min(1, 0.80 + 0.04 AR) factor that also sat
        # here counted it twice for AR under 5.  cl_max_3d_factor is the
        # rest: the built wing's gaps, finish and washout.
        self.stall_factor, self.stall_eta = stall_factor(
            self.AR, wing.taper, self._re_exp)
        self._cl_max_base = (af.cl_max_2d * build.cl_max_3d_factor
                             * self.stall_factor * self.fuselage_clmax_factor)

        # ---- junction interference (computed) -----------------------------
        # Hoerner per junction.  The wing roots meet the fuselage at the
        # chord a half-width out; a top wing has the junction on its lower
        # surface only, so half the through-wing value.  The tails meet the
        # tail post: two for the stabiliser, one for the fin.
        tc_t = build.tail_t_over_c
        c_h = math.sqrt(max(body.s_horiz_m2, 1e-6) / build.tail_aspect_ratio)
        c_v = math.sqrt(max(body.s_vert_m2, 1e-6)
                        / (0.5 * build.tail_aspect_ratio))
        self.wing_junction_m2 = (2.0 * 0.5 * max(n_cov, 1)
                                 * junction_drag_area_m2(
                                     af.t_over_c,
                                     wing.chord_at(0.5 * body.fuse_width_m)))
        self.junction_drag_area_m2 = (
            self.wing_junction_m2
            + 2.0 * junction_drag_area_m2(tc_t, c_h)
            + junction_drag_area_m2(tc_t, c_v)
            + (2.0 * junction_drag_area_m2(build.winglet_t_over_c,
                                           wing.tip_chord_m)
               if wing.winglet_h_m > 0.0 else 0.0))

        # ---- the propeller slipstream (computed) --------------------------
        # Two effects on top of the scrubbing drag (installed_thrust_factor):
        # the slipstream is turbulent, so the wing behind the propeller
        # loses part of its laminar run (drag_parts, "slipstream
        # transition"); and it is fast, so that wing lifts more at the same
        # angle of attack (slipstream_gains - lift and stall).
        from .layout import motor_stations
        self.prop_diameter_m = max(prop_diameter_m, 0.0)
        self.transition_frac = min(max(build.slipstream_transition_frac,
                                       0.0), 1.0)
        n_nose, ys, share = motor_stations(n_motors, wing.span_m,
                                           body.fuse_width_m,
                                           wing_prop_layout,
                                           self.prop_diameter_m,
                                           props_in_span)
        self._n_nose, self._wing_ys, self._share = n_nose, ys, share
        # The thrust line's height above the CG, which the trim carries
        # (m116): the nose propeller on the fuselage centreline, the wing
        # propellers on the wing's line - the fuselage top for a top wing,
        # its middle for a mid wing, as Aircraft sets the gear - each motor
        # alike, the CG at mid-height as the rotation case takes it.
        n_wing = 2 * len(ys)
        self.thrust_line_dz_m = (
            n_wing / float(max(n_nose + n_wing, 1))
            * (0.5 * body.fuse_height_m if build.wing_covered_surfaces == 1
               else 0.0))
        # full-throttle thrust and torque against airspeed, per propeller
        # (set_power_table); without it, no swirl and no roll trim
        self._pt_v: List[float] = []
        self._pt_t: List[float] = []
        self._pt_q: List[float] = []
        self._terms_cache: dict = {}
        # Each wing nacelle is a body the wing passes through: a Hoerner
        # junction on each side of it, at the local section thickness.
        self.nacelle_junction_m2 = sum(
            2.0 * junction_drag_area_m2(af.t_over_c, wing.chord_at(y))
            for y in ys)
        self.junction_drag_area_m2 += self.nacelle_junction_m2
        self._fuse_washed = bool(n_nose) and self.prop_diameter_m > 0.0
        self.washed_strip_m2 = self._washed_strip_m2(self.prop_diameter_m)
        # Where each slipstream meets the wing: (spanwise centre, height of
        # the wing above the thrust line, axial distance from the disc to
        # the quarter chord, the most span it can cover, and the half-width
        # of the body behind the propeller).  A top wing sits half a
        # fuselage height above a nose propeller's axis, so it cuts the
        # slipstream below its full width; a wing propeller is on the
        # wing's own line.  The body - the fuselage behind a nose
        # propeller, the nacelle behind a wing one - takes the middle of
        # the slipstream: only the wing either side of it is blown.
        r = 0.5 * self.prop_diameter_m
        self._jets = []
        if self.prop_diameter_m > 0.0:
            if n_nose:
                h = (0.5 * body.fuse_height_m if n_cov == 1 else 0.0)
                x = (max(body.nose_length_m, 0.3 * self.prop_diameter_m)
                     + 0.25 * wing.root_chord_m)
                self._jets.append((0.0, h, x, 0.5 * wing.span_m,
                                   0.5 * body.fuse_width_m))
            for y in ys:
                if y > 0.0:
                    self._jets.append((y, 0.0,
                                       0.55 * r + 0.25 * wing.chord_at(y),
                                       0.5 * share, 0.5 * nacelle_width_m))
        self._ho_mean = _hough_ordway_mean(build.prop_hub_frac)
        self._gain_cache: dict = {}

        # Tail geometry for trim drag
        self.s_tail_h = body.s_horiz_m2
        self.tail_arm = max(body.tail_arm_m, 1e-3)
        self.b_tail = math.sqrt(max(build.tail_aspect_ratio * self.s_tail_h,
                                    1e-6))

        # Representative CLmax for reporting: the base value, before the
        # Reynolds scaling cl_max_at_re applies
        self.cl_max = self._cl_max_base

        # ---- longitudinal stability -------------------------------------
        # Tail volume ratio and the neutral point, on the mean aerodynamic
        # chord (x/c below means x/MAC from the MAC's leading edge).  With
        # the wing AC at 0.25 MAC, the tail arm l_t measured AC to AC and
        # the tail's lift counted in the aeroplane's (as tail_load counts
        # it), dM_cg/dalpha = 0 puts the neutral point at
        #     x_np/c = 0.25 + [eta_t V_h (a_t/a_w)(1 - de/dalpha) - fus]
        #                     / (1 + eta_t (S_h/S)(a_t/a_w)(1 - de/dalpha))
        # - the divisor is the tail's share of the lift slope, which the
        # usual textbook approximation (no divisor) leaves out; fus is the
        # destabilising fuselage term.
        self.a_tail = lift_curve_slope(build.tail_aspect_ratio)
        self.v_h = (self.s_tail_h * self.tail_arm
                    / max(self.S * wing.mac_m, 1e-12))
        tail_term = (self.eta_tail * self.v_h
                     * (self.a_tail / max(self.a_wing, 1e-9))
                     * (1.0 - self.downwash))
        fus_term = (build.fuselage_np_shift_k * body.fuse_width_m ** 2
                    * body.fuse_length_m
                    / max(self.S * wing.mac_m * self.a_wing, 1e-12))
        lift_share = (self.eta_tail * self.s_tail_h / max(self.S, 1e-12)
                      * (self.a_tail / max(self.a_wing, 1e-9))
                      * (1.0 - self.downwash))
        self.x_np_over_c = 0.25 + (tail_term - fus_term) / (1.0 + lift_share)
        # The CG is placed, not assumed: the layout puts it the target
        # static margin ahead of this neutral point (layout.py), and
        # Aircraft calls set_cg() with where the parts actually put it.
        self.x_cg_over_c = self.x_np_over_c - build.target_static_margin
        self.static_margin = build.target_static_margin

    def set_cg(self, x_cg_over_c: float) -> None:
        """Move the CG to where the layout's parts put it (MAC fraction)."""
        self.x_cg_over_c = x_cg_over_c
        self.static_margin = self.x_np_over_c - x_cg_over_c

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
    def _section_cd(self, v: float) -> float:
        """Built-wing section drag per unit planform (both surfaces)."""
        wing = self.wing
        re_w = max(v, 3.0) * wing.chord_m / self.env.nu
        re_scale = (wing.airfoil.re_ref / max(re_w, 3.0e4)) ** 0.15
        return wing.airfoil.cd_min * re_scale * self.build.profile_drag_factor

    def _turbulent_section_cd(self, v: float) -> float:
        """The same section with no laminar run: fully turbulent flat-plate
        friction on both surfaces times Raymer's thickness form factor,
        with the same built-wing allowance."""
        wing = self.wing
        re_w = max(v, 3.0) * wing.chord_m / self.env.nu
        return (2.0 * skin_friction_cf(re_w, 0.0)
                * tail_form_factor(wing.airfoil.t_over_c)
                * self.build.profile_drag_factor)

    def _washed_section_cd(self, v: float) -> float:
        """Section drag inside the slipstream: the smooth-air value moved
        `slipstream_transition_frac` of the way to fully turbulent."""
        cd = self._section_cd(v)
        return cd + self.transition_frac * max(
            self._turbulent_section_cd(v) - cd, 0.0)

    def _washed_strip_m2(self, d: float) -> float:
        """Wetted-equivalent wing planform inside the slipstreams (each
        surface counted as half), on the wing's own area scale.  A nose
        tractor washes a strip one diameter wide, less the surface the
        fuselage covers; wing propellers a strip each, at most their share
        of the span."""
        if d <= 0.0:
            return 0.0
        n_cov = min(max(self.build.wing_covered_surfaces, 0), 2)
        scale = self.S / max(self.wing.planform_within_m2(self.wing.span_m),
                             1e-12)
        a = 0.0
        if self._n_nose:
            a += (self.wing.planform_within_m2(d)
                  - 0.5 * n_cov * self.wing.planform_within_m2(
                      min(self.body.fuse_width_m, d)))
        a += sum(min(d, self._share) * self.wing.chord_at(y)
                 for y in self._wing_ys)
        return max(a, 0.0) * scale

    def _fuselage_drag_area(self, v: float) -> float:
        """Skin-friction drag area of the fuselage's exposed skin, m^2.
        Behind a nose propeller its laminar run is lost in proportion too."""
        body = self.body
        re_f = max(v, 3.0) * body.fuse_length_m / self.env.nu
        lam = 0.10 * ((1.0 - self.transition_frac) if self._fuse_washed
                      else 1.0)
        cf_f = skin_friction_cf(re_f, lam)
        f = body.fineness
        ff_f = 1.0 + 60.0 / max(f, 2.0) ** 3 + max(f, 2.0) / 400.0
        return cf_f * ff_f * max(body.fuse_wetted_m2 - self.fuse_covered_m2,
                                 0.0)

    def drag_parts(self, v: float) -> dict:
        """Parasite drag coefficients by component, on the wing area:
        wing, fuselage, tails, excrescence, junctions, discrete."""
        body, b = self.body, self.build
        v = max(v, 3.0)
        cd_wing = self._section_cd(v) * self.wing_wetted_frac
        cd_fuse = self._fuselage_drag_area(v) / self.S

        s_tail = body.s_horiz_m2 + body.s_vert_m2
        c_tail = math.sqrt(max(s_tail, 1e-4) / b.tail_aspect_ratio)
        re_t = v * c_tail / self.env.nu
        cf_t = skin_friction_cf(re_t, 0.20)
        tc_t = b.tail_t_over_c
        # Raymer: wetted area of a thin surface 1.977 + 0.52 t/c times its
        # planform, with his thickness form factor
        cd_tail = (cf_t * tail_form_factor(tc_t) * (1.977 + 0.52 * tc_t)
                   * s_tail / self.S)

        # winglets (m110): two thin surfaces, root chord the tip chord
        cd_wl = 0.0
        if self.wing.winglet_h_m > 0.0:
            c_wl = 0.75 * self.wing.tip_chord_m
            cf_wl = skin_friction_cf(v * c_wl / self.env.nu, 0.20)
            tw = b.winglet_t_over_c
            cd_wl = (cf_wl * tail_form_factor(tw) * (1.977 + 0.52 * tw)
                     * 2.0 * self.wing.winglet_area_m2 / self.S)
        area_extra = (b.gear_drag_area_m2
                      + self.nacelle_drag_area_m2 * self.n_wing_motors
                      + self.extra_drag_area_m2)
        # the wing behind the propeller, with its laminar run cut short
        cd_trans = (max(self._washed_section_cd(v) - self._section_cd(v),
                        0.0) * self.washed_strip_m2 / self.S)
        skin = cd_wing + cd_fuse + cd_tail + cd_trans + cd_wl
        return {"wing": cd_wing, "fuselage": cd_fuse, "tails": cd_tail,
                "winglets": cd_wl,
                "slipstream transition": cd_trans,
                "excrescence": skin * b.excrescence_frac,
                "junctions": self.junction_drag_area_m2 / self.S,
                "discrete": area_extra / self.S}

    def cd0(self, v: float) -> float:
        """Parasite drag: the airframe's skin friction plus its leakage and
        protuberance allowance, the computed junction interference, and
        the discrete items (gear, nacelles, the deployed sensor)."""
        return sum(self.drag_parts(v).values())

    # ------------------------------------------------------------------
    # Installed thrust: the slipstream scrubbing the airframe
    # ------------------------------------------------------------------
    def slipstream_drag_area_m2(self, prop_diameter_m: float,
                                v: float) -> float:
        """Drag area (CD*A, m^2) of the airframe inside the propeller
        slipstreams, at airspeed v (its skin friction falls with Reynolds
        number).

        A nose tractor washes the whole fuselage and the wing across a
        strip one propeller diameter wide (a top wing sits within half a
        fuselage height of the thrust line); wing-mounted propellers each
        wash a strip of wing one diameter wide (at most their share of the
        span), and their nacelles.  An odd motor count has one of each
        kind of propeller (layout.motor_stations).  The tails are left out:
        downwash and height carry most of them clear.
        """
        b = self.build
        d = max(prop_diameter_m, 0.0)
        # the washed wing, at its turbulent-slipstream section drag
        f = (self._washed_section_cd(v) * self._washed_strip_m2(d)
             * (1.0 + b.excrescence_frac))
        if self._n_nose:
            f += (self._fuselage_drag_area(v) * (1.0 + b.excrescence_frac)
                  + self.wing_junction_m2)
        f += (self.nacelle_drag_area_m2 * len(self._wing_ys)
              + self.nacelle_junction_m2)
        return f

    def install_factor_at(self, v: float) -> float:
        """installed_thrust_factor for this model's propellers at v."""
        if self.prop_diameter_m <= 0.0:
            return 1.0
        return self.installed_thrust_factor(self.prop_diameter_m, v)

    def installed_thrust_factor(self, prop_diameter_m: float,
                                v: float) -> float:
        """Installed thrust / free-air thrust, from slipstream scrubbing.

        Actuator-disc momentum: far downstream the slipstream's dynamic
        pressure exceeds the free stream's by T/A (from V_w^2 - V^2 =
        2T/(rho A)).  Everything of drag area f inside it therefore drags
        (T/A) f more than in clean air - a fixed fraction f/A of the
        thrust, whatever the speed.  So the installed thrust is

            T_installed = T (1 - f / A_disc)

        This replaces a flat 0.95.  Taking the far-wake increment all
        along the fuselage slightly overstates it near the nose, where the
        slipstream has not finished accelerating.
        """
        d = max(prop_diameter_m, 1e-3)
        a_disc = max(self.n_motors, 1) * math.pi * 0.25 * d * d
        f = self.slipstream_drag_area_m2(d, v)
        return max(0.80, 1.0 - f / a_disc)

    # ------------------------------------------------------------------
    # Slipstream lift
    # ------------------------------------------------------------------
    _GAIN_STEP = 0.05          # far-wake speed increment grid for the gains
    _GAIN_MAX = 3.0            # ... and its ceiling (near-static thrust)

    def slipstream_gains(self, v: float, thrust_n: float
                         ) -> "tuple[float, float]":
        """(lift, CLmax) of the wing with the propellers making `thrust_n`
        (installed, all of them) at airspeed v, each as a ratio to the
        power-off wing.

        Momentum theory gives the slipstream: far behind the disc its speed
        is V sqrt(1 + T/(q A)); at the wing it has built up a fraction
        (1 + x/sqrt(R^2 + x^2))/2 of that increment, x the axial distance
        from the disc, and contracted by continuity.  The wing's response
        is the lifting line with the local airspeed in it
        (liftingline.blown_wing), which is what keeps a jet only a
        propeller wide from being credited with its full dynamic pressure.
        Tabulated on the far-wake increment and interpolated.

        The slipstream is not uniform: the thrust, and with it the speed
        increment, follows the propeller's radial loading - none inside
        the spinner, a peak near 70% radius, none at the tip
        (_thrust_shape).  And the body behind the propeller takes its
        middle: the strip the fuselage (nose) or nacelle (wing) covers
        blows no wing.  Each exposed strip goes to the lifting line as one
        band at the speed that keeps its dynamic-pressure lift.
        """
        if thrust_n <= 0.0 or not self._jets:
            return 1.0, 1.0
        r = 0.5 * self.prop_diameter_m
        a_disc = math.pi * r * r
        tp = (thrust_n / max(self.install_factor_at(v), 1e-6)
              / max(self.n_motors, 1))
        ratio = tp / max(self.q(max(v, 0.5)) * a_disc, 1e-9)
        # sqrt(1 + x) - 1 written as x / (sqrt(1 + x) + 1): the same
        # number, without the cancellation at a light thrust load
        du = min(ratio / (math.sqrt(1.0 + ratio) + 1.0), self._GAIN_MAX)
        z = du / self._GAIN_STEP
        i = int(math.floor(z))
        t = z - i
        g0 = self._gains_at(i)
        if t < 1e-12:
            return g0
        g1 = self._gains_at(i + 1)
        return (g0[0] + t * (g1[0] - g0[0]), g0[1] + t * (g1[1] - g0[1]))

    def _gains_at(self, i: int) -> "tuple[float, float]":
        hit = self._gain_cache.get(i)
        if hit is not None:
            return hit
        du = i * self._GAIN_STEP             # far-wake speed increment, 2a
        jets = self._jet_bands(du)
        eta_b = min(self.body.fuse_width_m / self.wing.span_m, 0.5)
        cl0, f0 = blown_wing(self.AR, self.wing.taper, (), self._re_exp,
                             eta_body=eta_b)
        cl, fs = blown_wing(self.AR, self.wing.taper, tuple(jets),
                            self._re_exp, eta_body=eta_b)
        out = (cl / cl0, fs / f0)
        self._gain_cache[i] = out
        return out

    def set_power_table(self, vs, thrust_n, torque_nm) -> None:
        """Full-throttle total thrust and per-propeller torque against
        airspeed, for the swirl and the roll trim (Aircraft sets it)."""
        self._pt_v, self._pt_t, self._pt_q = (list(vs), list(thrust_n),
                                              list(torque_nm))

    def _torque(self, v: float, thrust_n: float) -> float:
        """Per-propeller torque making `thrust_n` (all propellers) at v:
        the full-throttle torque scaled by the share of full thrust."""
        vs = self._pt_v
        if not vs or thrust_n <= 0.0:
            return 0.0
        x = min(max(v, vs[0]), vs[-1])
        i = min(int((x - vs[0]) / (vs[1] - vs[0])), len(vs) - 2)
        t = (x - vs[i]) / (vs[i + 1] - vs[i])
        tf = self._pt_t[i] + t * (self._pt_t[i + 1] - self._pt_t[i])
        qf = self._pt_q[i] + t * (self._pt_q[i + 1] - self._pt_q[i])
        return qf * min(thrust_n / max(tf, 1e-6), 1.5)

    def swirl_angle(self, v: float, thrust_n: float) -> float:
        """The swirl angle in a propeller's slipstream, rad: the
        angular-momentum flux Q = rho u v_t r A at the torque radius
        (0.75 R), over the axial speed V + w (momentum theory)."""
        q_p = self._torque(v, thrust_n)
        if q_p <= 0.0 or self.prop_diameter_m <= 0.0:
            return 0.0
        r = 0.5 * self.prop_diameter_m
        a = math.pi * r * r
        rho = self.env.rho
        t_p = thrust_n / max(self.n_motors, 1)
        w = 0.5 * (-v + math.sqrt(v * v + 2.0 * t_p / (rho * a)))
        u = max(v + w, 1.0)
        return q_p / (rho * a * 0.75 * r * u * u)

    def _terms(self, i: int):
        hit = self._terms_cache.get(i)
        if hit is None:
            hit = blown_wing_terms(self.AR, self.wing.taper,
                                   tuple(self._jet_bands(i * self._GAIN_STEP)))
            self._terms_cache[i] = hit
        return hit

    def induced_factor(self, v: float, thrust_n: float, cl: float) -> float:
        """Induced drag with the propellers making `thrust_n`, over the
        planform's at the same lift: 1, less what swirl recovery gives
        back (blown_wing_terms).

        Only the swirl's share is credited.  The lifting line with the
        jet's faster air in it would also lower the induced drag of the
        axial slipstream alone - lift made on less circulation - but that
        leaves the jet's own wake out of the account: an axial slipstream
        by itself does not lower the span's induced drag, swirl recovery
        does (Kroo 1986, propeller-wing integration).  So the axial part
        stays the planform's L^2 / (q pi b^2 e), as it always was, and
        the factor carries the swirl's change to it.  A nose propeller's
        swirl is antisymmetric and not credited (_jet_bands)."""
        if thrust_n <= 0.0 or not self._jets or abs(cl) < 1e-6:
            return 1.0
        r = 0.5 * self.prop_diameter_m
        tp = thrust_n / max(self.install_factor_at(v), 1e-6) / max(self.n_motors, 1)
        ratio = tp / max(self.q(max(v, 0.5)) * math.pi * r * r, 1e-9)
        du = min(ratio / (math.sqrt(1.0 + ratio) + 1.0), self._GAIN_MAX)
        z = du / self._GAIN_STEP
        i = int(math.floor(z))
        t = z - i
        a0, a1 = self._terms(i), self._terms(i + 1)
        cla, clp, saa, sap, spp, wa, wp = (x + t * (y - x)
                                           for x, y in zip(a0, a1))
        z0 = self._terms(0)
        eps = self.swirl_angle(v, thrust_n) if self._wing_ys else 0.0
        if eps <= 0.0:
            return 1.0
        alpha = (cl - eps * clp) / max(cla, 1e-9)
        cdi = (4.0 * math.pi * self.AR
               * (alpha * alpha * saa + 2.0 * alpha * eps * sap
                  + eps * eps * spp)
               - eps * (alpha * wa + eps * wp))
        alpha0 = cl / max(cla, 1e-9)                    # the same, no swirl
        cdi_ax = 4.0 * math.pi * self.AR * alpha0 * alpha0 * saa
        cdi0 = 4.0 * math.pi * self.AR * z0[2] * (cl / max(z0[0], 1e-9)) ** 2
        return max(1.0 + (cdi - cdi_ax) / max(cdi0, 1e-12), 0.05)

    def roll_trim_drag(self, v: float, thrust_n: float) -> float:
        """Drag of holding the wings level against an unpaired motor's
        torque - the nose motor; wing pairs counter-rotate.  The ailerons
        balance it with an antisymmetric lift, whose induced drag is
        32 Cl^2 / (pi AR) on q S, Cl = Q / (q S b)."""
        if not self._n_nose:
            return 0.0
        q_t = self._n_nose * self._torque(v, thrust_n)
        if q_t <= 0.0:
            return 0.0
        q = self.q(max(v, 0.5))
        c_roll = q_t / (q * self.S * self.wing.span_m)
        return q * self.S * 32.0 * c_roll * c_roll / (math.pi * self.AR)

    def _thrust_shape(self, x: float) -> float:
        """Thrust per unit disc area at radius fraction x, relative to the
        disc's mean: the Hough-Ordway distribution, r* sqrt(1 - r*) with
        r* = (x - hub) / (1 - hub) - the standard body-force propeller of
        actuator-disc CFD.  Zero inside the spinner and at the tip."""
        xh = self.build.prop_hub_frac
        if x <= xh or x >= 1.0:
            return 0.0
        rs = (x - xh) / (1.0 - xh)
        return rs * math.sqrt(1.0 - rs) / self._ho_mean

    def _jet_bands(self, du: float) -> list:
        """The blown bands on the half span, (eta centre, eta half-width,
        speed ratio), for a mean far-wake increment du."""
        r = 0.5 * self.prop_diameter_m
        half_b = 0.5 * self.wing.span_m
        ratio = (1.0 + du) ** 2 - 1.0        # the disc's mean T / (q A)
        jets = []
        for y, h, x, cap, block in self._jets:
            f = 0.5 * (1.0 + x / math.sqrt(r * r + x * x))
            u = 1.0 + f * du
            r_s = r * math.sqrt((1.0 + 0.5 * du) / u)   # continuity
            # (r - h)(r + h), not r^2 - h^2: exact where the jet grazes
            half = min(math.sqrt(max((r_s - h) * (r_s + h), 0.0)), cap)
            if half <= block or u <= 1.0:
                continue
            # the exposed strip each side of the body: its speed is the
            # rms of the local speed across it (lift goes as u^2)
            n = 16
            acc = 0.0
            for k in range(n):
                s = block + (k + 0.5) * (half - block) / n
                loc = ratio * self._thrust_shape(min(math.hypot(s, h) / r_s,
                                                     1.0))
                acc += (1.0 + f * loc / (math.sqrt(1.0 + loc) + 1.0)) ** 2
            u_eq = math.sqrt(acc / n)
            mid, hw = 0.5 * (block + half), 0.5 * (half - block)
            edge = 0.25 * half / half_b     # the whole slipstream's layer
            # Swirl: a wing propeller turns inboard-up (the pair counter-
            # rotate), so its slipstream rises on the inboard side and
            # falls outboard.  A nose propeller's is antisymmetric across
            # the span: it cancels in the symmetric lifting line to first
            # order and is left out.
            if y == 0.0:
                jets.append((mid / half_b, hw / half_b, u_eq, edge, 0.0))
            else:
                jets.append((abs(y - mid) / half_b, hw / half_b, u_eq,
                             edge, 1.0))
                jets.append(((y + mid) / half_b, hw / half_b, u_eq, edge,
                             -1.0))
        return jets

    # ------------------------------------------------------------------
    # Trim
    # ------------------------------------------------------------------
    def tail_load(self, v: float, lift_n: float,
                  thrust_n: float = 0.0) -> float:
        """Horizontal tail lift required to trim, in newtons.

        Moments about the CG (nose-up positive, x aft), with the tail arm
        l_t measured from the wing's aerodynamic centre to the tail's - as
        build_body and the layout measure it - so the tail sits
        l_t - (x_cg - x_ac) behind the CG:

            M_ac + L_w (x_cg - x_ac) - L_t (l_t - (x_cg - x_ac)) = 0

        Writing k = (x_cg - x_ac)/c (c the mean aerodynamic chord) and
        using L_w + L_t = L (the tail is part of the lift balance), the
        L_t k c terms cancel:

            L_t = c (Cm_ac q S + L k) / l_t

        (M_ac = Cm_ac q S c on the MAC is the wing's pitching moment about
        its AC exactly, for a constant section Cm_ac.)

        Two terms, and both matter:

        * the **camber** term, Cm_ac q S c, is negative for a cambered
          section and scales with dynamic pressure only;
        * the **CG** term, L k c, scales with the lift being carried, so in
          an n-g turn the tail load grows with n.  Modelling only the camber
          term - as this did previously - holds the tail load constant
          through a 3 g turn, which is wrong.

        Negative means download.  With the propellers making `thrust_n`,
        the blown wing's camber moment grows with its lift, and a thrust
        line dz above the CG (thrust_line_dz_m) pitches the nose down, a
        third term (m116):

            L_t = [c (Cm_ac q S + L k) - T dz] / l_t
        """
        q = self.q(v)
        c = self.wing.mac_m
        k = self.x_cg_over_c - 0.25          # CG aft of the wing AC, MACs
        g = self.slipstream_gains(v, thrust_n)[0] if thrust_n > 0.0 else 1.0
        return ((c * (self.wing.airfoil.cm_ac * q * self.S * g + lift_n * k)
                 - thrust_n * self.thrust_line_dz_m) / self.tail_arm)

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

    def cl_required(self, v: float, lift_n: float,
                    thrust_n: float = 0.0) -> float:
        """Wing CL needed (on the free stream).

        The tail carries part of the load, so the wing carries L - L_t.
        With a download (L_t negative) the wing carries more than the
        aircraft weight.  With the propellers making `thrust_n` the blown
        wing makes the same lift at a lower CL (slipstream_gains).
        """
        g = self.slipstream_gains(v, thrust_n)[0] if thrust_n > 0.0 else 1.0
        qs = self.q(v) * self.S * g
        wing_lift = lift_n - self.tail_load(v, lift_n, thrust_n)
        return wing_lift / max(qs, 1e-9)

    def cd_from_cl(self, cl: float, v: float) -> float:
        """Wing drag coefficient: profile bucket + stall rise + induced."""
        af = self.wing.airfoil
        cd = self.cd0(v)
        # the section's own polar when the catalogue has one, else the
        # parametric bucket + stall rise - see Airfoil.cd_excess
        cd += af.cd_excess(cl)
        cd += self.k_induced * cl * cl
        return cd

    def drag(self, v: float, lift_n: float, thrust_n: float = 0.0) -> float:
        """Total drag at this speed carrying this much lift, with the
        propellers making `thrust_n` (0: power off).

        Power on, the blown wing works at a lower CL, which lowers the
        section's drag above its bucket.  The induced drag stays that of
        the whole lift on the span, L^2/(q pi b^2 e): the slipstream moves
        where the lift is made, not how much the span has to carry, and a
        loading peaked behind the propeller is if anything less elliptic.
        """
        q = self.q(v)
        cl = self.cl_required(v, lift_n, thrust_n)
        if thrust_n > 0.0:
            g = self.slipstream_gains(v, thrust_n)[0]
            cl_free = cl * g                  # the same lift on q S
            af = self.wing.airfoil
            d = q * self.S * (self.cd0(v) + af.cd_excess(cl)
                              + self.k_induced * cl_free * cl_free
                              * self.induced_factor(v, thrust_n, cl_free))
            d += self.roll_trim_drag(v, thrust_n)
        else:
            d = q * self.S * self.cd_from_cl(cl, v)
        # induced drag of the trimming tail load, L_t^2 / (q pi b_t^2 e_t)
        lt = self.tail_load(v, lift_n, thrust_n)
        if abs(lt) > 1e-9 and self.b_tail > 1e-6:
            d += lt * lt / (self.eta_tail * q * math.pi * self.b_tail ** 2
                            * self.build.tail_efficiency)
        return d

    @property
    def alpha_zero_lift(self) -> float:
        """The wing's zero-lift angle, rad: the section's own, from its
        XFOIL polar (Airfoil.alpha_l0_deg, m109).  A hand-entered section
        without one falls back on the camber-moment estimate of thin-airfoil
        theory (parabolic camber), 4 cm_ac / pi - far off for aft-loaded
        high-lift sections (-19.5 deg for the S1223 against -10.0)."""
        af = self.wing.airfoil
        if af.alpha_l0_deg is not None:
            return math.radians(af.alpha_l0_deg)
        return 4.0 * af.cm_ac / math.pi

    def cl_at_ground_alpha(self) -> float:
        """The most wing lift coefficient the aeroplane can reach with its
        wheels on the runway: the 3-D lift slope times
        BuildStandard.max_ground_alpha_deg, measured from zero lift (m104).
        A low-aspect-ratio wing makes its lift slowly with angle, so this -
        not its CLmax - can set its liftoff and touchdown speeds."""
        a = math.radians(self.build.max_ground_alpha_deg) - self.alpha_zero_lift
        return max(self.a_wing * a, 0.0)

    def max_lift(self, v: float, thrust_n: float = 0.0,
                 ground: bool = False) -> float:
        """Maximum supportable weight at this speed, after trim, with the
        propellers making `thrust_n` (0: power off).  `ground`: with the
        wheels on the runway, where the wing can reach no more than
        cl_at_ground_alpha (m104).

        The wing stalls at L_w = q S CLmax.  Substituting that and the trim
        relation L_t = (Cm_ac q S c + L k c - T dz) / l_t (tail_load: AC-to-AC
        arm, c the MAC, dz the thrust line over the CG) into L = L_w + L_t
        and solving for L:

            L = [ q S CLmax l_t + Cm_ac q S c - T dz ] / (l_t - k c)

        Power on, CLmax and the camber moment carry the slipstream's gains.
        """
        q = self.q(v)
        c = self.wing.mac_m
        k = self.x_cg_over_c - 0.25          # CG aft of the wing AC, MACs
        g_l, g_max = (self.slipstream_gains(v, thrust_n) if thrust_n > 0.0
                      else (1.0, 1.0))
        qs = q * self.S
        lw_max = qs * self.cl_max_at_v(v) * g_max
        if ground:
            lw_max = min(lw_max, qs * self.cl_at_ground_alpha() * g_l)
        # (a CG at the tail's AC could not be trimmed at all; the static-
        # margin checks reject anything near it long before)
        return ((lw_max * self.tail_arm
                 + self.wing.airfoil.cm_ac * qs * g_l * c
                 - thrust_n * self.thrust_line_dz_m)
                / max(self.tail_arm - k * c, 0.05 * self.tail_arm))

    def v_stall_powered(self, weight_n: float, thrust_at,
                        load_factor: float = 1.0,
                        ground: bool = False) -> float:
        """Stall speed with the propellers at the thrust `thrust_at(v)` -
        full throttle, for a takeoff or a hard turn.  Lower than the
        power-off stall, by the slipstream's lift.  `ground`: the slowest
        the runway attitude allows (max_lift)."""
        target = load_factor * weight_n

        def residual(v: float) -> float:
            return self.max_lift(v, thrust_at(v), ground) - target

        lo, hi = 1.0, 120.0
        if residual(lo) >= 0.0:
            return lo
        if residual(hi) <= 0.0:
            return hi
        return brent(residual, lo, hi, xtol=1e-9)

    def v_stall(self, weight_n: float, load_factor: float = 1.0,
                ground: bool = False) -> float:
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
            return self.max_lift(v, 0.0, ground) - target

        lo, hi = 1.0, 120.0
        if residual(lo) >= 0.0:
            return lo
        if residual(hi) <= 0.0:
            return hi
        return brent(residual, lo, hi, xtol=1e-12)

    def v_stall_in_air(self, weight_n: float, rho: float, nu: float,
                       load_factor: float = 1.0,
                       ground: bool = False) -> float:
        """v_stall for this wing, trim and CG in other air, power off (m118):
        the density sets the dynamic pressure and the viscosity the
        Reynolds number, which is all of the air the power-off stall sees."""
        keep = self.env
        self.env = replace(keep, rho=rho, nu=nu)
        try:
            return self.v_stall(weight_n, load_factor, ground)
        finally:
            self.env = keep

    def ld_ratio(self, v: float, weight_n: float) -> float:
        return weight_n / max(self.drag(v, weight_n), 1e-9)

    def best_ld_speed(self, weight_n: float) -> float:
        """Speed for best L/D: coarse scan to bracket, then golden section,
        up to 3x the minimum-drag speed of the parabolic polar,
        sqrt(2W/(rho S)) (k/CD0)^0.25 - well past where L/D peaks."""
        lo = self.v_stall(weight_n) * 1.02
        v_md = (math.sqrt(2.0 * weight_n / (self.env.rho * self.S))
                * (self.k_induced / max(self.cd0(lo), 1e-6)) ** 0.25)
        hi = max(3.0 * v_md, 1.5 * lo)
        v, _ = scan_then_refine(lambda x: self.ld_ratio(x, weight_n),
                                lo, hi, n_scan=20, tol=1e-6)
        return v


def build_body(payload_bay_length_m: float, bay_width_m: float,
               bay_height_m: float, wing: WingGeometry,
               build: BuildStandard, nose_length_m: float = 0.25,
               cg_minus_ac_m: float = 0.0, vol_h: Optional[float] = None,
               vol_v: Optional[float] = None) -> BodyGeometry:
    """Fuselage and tails from the payload bay, the nose and the tail arm.

    Nothing is a fixed proportion.  The cross-section wraps the payload
    bay.  The tail arm is the lightest one (layout.tail_arm_optimum_m) but
    never shorter than the tailcone needs to taper from the bay to the tail
    post at a non-separating angle.  The tails follow from the arm and the
    volume coefficients.  The payload bay is centred on the CG, which sits
    `cg_minus_ac_m` behind the wing's aerodynamic centre, so the fuselage
    runs: nose + bay + tailcone to the tail's trailing edge.  The nose
    length comes from the balance (layout.solve_layout); the default is
    only a first guess for that iteration.
    """
    from .layout import (section_perimeters, tail_arm_optimum_m,
                         tail_chord_m, tailcone_min_length_m)
    width = max(bay_width_m + 0.035, 0.10)
    height = max(bay_height_m + 0.045, 0.11)
    s = wing.area_m2

    # the tail volumes: the build standard's fixed ones, or what the
    # requirements call for (Aircraft._tail_volumes_required, m111)
    vh = build.vol_coeff_h if vol_h is None else vol_h
    vv = build.vol_coeff_v if vol_v is None else vol_v

    def c_h_at(arm_m: float) -> float:
        return tail_chord_m(vh * s * wing.mac_m / arm_m, build)

    arm_opt = tail_arm_optimum_m(wing, width, height, build, vh, vv)
    # tail trailing edge at least the minimum tailcone behind the bay end,
    # measured from the wing AC (which sits cg_minus_ac ahead of the bay
    # centre): arm + 0.75 c_h >= bay/2 + cg_minus_ac + tailcone_min
    need = (0.5 * payload_bay_length_m + cg_minus_ac_m
            + tailcone_min_length_m(width, height, build))
    arm_min = max(need, 0.10)
    for _ in range(40):                       # c_h shrinks as the arm grows
        arm_min = max(need - 0.75 * c_h_at(max(arm_min, 0.05)), 0.10)
    arm = max(arm_opt, arm_min)
    s_h = vh * s * wing.mac_m / arm
    s_v = vv * s * wing.span_m / arm
    c_h = tail_chord_m(s_h, build)
    nose = max(nose_length_m, 0.0)
    length = (nose + 0.5 * payload_bay_length_m - cg_minus_ac_m
              + arm + 0.75 * c_h)
    tailcone = max(length - nose - payload_bay_length_m, 0.0)
    p0, pe = section_perimeters(width, height, build)
    return BodyGeometry(length, width, height, arm, s_h, s_v, nose,
                        tailcone, pe, arm_opt, arm_min, p0)
