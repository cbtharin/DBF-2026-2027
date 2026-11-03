"""Where every part sits - and so where the centre of gravity is.

Nothing here is a fixed proportion of the fuselage.  The layout follows
from the parts and three design rules:

1.  **The payload bay is centred on the CG.**  Mission 1 flies empty,
    Mission 2 full, Mission 3 with the sensor towed from the bay; with the
    payload at the CG, loading or unloading it does not move the CG, so
    one airframe is balanced for all three.  (The sensor's tow point is at
    the bay.)

2.  **The wing is placed so that CG gives the target static margin.**  The
    neutral point is computed (aerodynamics.py); the wing's mean chord is
    positioned so the CG sits `target_static_margin` of a chord ahead of it.

3.  **The battery is the trim weight.**  The nose holds the motor (one
    nose tractor, for an odd motor count - see `motor_stations`), then a
    battery-and-avionics bay.  The battery slides in that
    bay to balance the aeroplane about the payload.  If it cannot go far
    enough forward, the nose is lengthened until it can.  If the aeroplane
    is nose-heavy even with the battery against the payload bay - a long
    bay puts the nose far ahead of the CG - the battery moves behind the
    bay, into the front of the tailcone, as far as the tapering tailcone
    still has room for it.  If even that cannot balance it, the CG ends up
    forward of target and the static margin is reported as it is -
    larger, never faked.

The tail arm is chosen too, not assumed: it is the arm that makes the
tails plus the tailcone lightest (a longer arm needs smaller tails for
the same volume coefficient, but more tailcone), and never shorter than
the tailcone needs to taper from the payload bay to the tail post without
the flow separating off it.

Part sizes come from their masses: a battery pack's volume from LiPo pack
density (makers' dimensions give 2.1 g/cm^3), its shape from typical pack
proportions; a motor's length from outrunner density.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

#: Typical LiPo brick proportions, length : width : height.
PACK_PROPORTIONS = (4.0, 1.3, 1.0)


@dataclass
class Station:
    name: str
    mass_kg: float
    x_m: float                 # centroid, metres aft of the nose
    x0_m: float = 0.0          # extent along the fuselage, for drawing
    x1_m: float = 0.0


@dataclass
class Layout:
    stations: List[Station]
    nose_length_m: float       # nose tip to the front of the payload bay
    motor_section_m: float     # motor + mount, or the nose cap for twins
    battery_length_m: float
    avionics_length_m: float
    bay_start_m: float
    bay_end_m: float
    tailcone_length_m: float
    fuse_length_m: float
    x_wing_le_mac_m: float
    x_wing_le_root_m: float
    x_wing_ac_m: float
    x_tail_ac_m: float
    tail_arm_m: float
    x_cg_target_m: float
    x_cg_m: float              # Mission 2, loaded
    x_cg_empty_m: float        # Mission 1
    x_cg_m3_m: float           # Mission 3, sensor towed from the bay
    x_battery_m: float
    x_main_gear_m: float
    x_nose_gear_m: float
    static_margin: float       # actual, Mission 2
    static_margin_empty: float
    balanced: bool
    note: str = ""
    moments: Dict[str, float] = field(default_factory=dict)


# ==========================================================================
# Part sizes from their masses
# ==========================================================================
def motor_length_m(mass_kg: float, build) -> float:
    """Outrunner can length, taking length ~ diameter."""
    m = max(mass_kg, 0.0)
    return (4.0 * m / (math.pi * build.motor_density_kg_m3)) ** (1.0 / 3.0)


def pack_box_m(mass_kg: float, build) -> Tuple[float, float, float]:
    """(length, width, height) of one pack from its mass."""
    vol = max(mass_kg, 1e-6) / build.pack_density_kg_m3
    a, b, c = PACK_PROPORTIONS
    h = (vol / (a * b * c)) ** (1.0 / 3.0)
    return a * h, b * h, c * h


def tail_chord_m(s_h: float, build) -> float:
    """Horizontal tail chord from its area and aspect ratio."""
    return math.sqrt(max(s_h, 1e-9) / build.tail_aspect_ratio)


#: Where the wing propellers go (Design.wing_prop_layout):
#:   "mid"      each in the middle of an equal share of the exposed semi-span
#:   "tip"      the outermost each side on the wingtip - its axis on the
#:              tip, the disc overhanging, since rule 3.1 limits the wing's
#:              span and says nothing of propellers (with
#:              rules.wingspan_counts_props, its disc's edge on the tip
#:              instead) - and the rest evenly inboard of it
#:   "stagger"  neighbours alternately ahead of and behind each other, so
#:              their discs may overlap in front view, up to
#:              STAGGER_MAX_OVERLAP of a diameter: the innermost 3 cm clear
#:              of the fuselage, the outermost at the tip, the rest between
WING_PROP_LAYOUTS = ("mid", "tip", "stagger")
STAGGER_MAX_OVERLAP = 0.30
PROP_GAP_M = 0.03


def disc_overlap_fraction(step_m: float, prop_d_m: float) -> float:
    """The share of a propeller disc's area its neighbour's disc covers in
    front view, the centres `step_m` apart: the lens of two equal circles
    over one circle's area."""
    r = 0.5 * prop_d_m
    if prop_d_m <= 0.0 or step_m >= prop_d_m:
        return 0.0
    s = max(step_m, 0.0)
    lens = (2.0 * r * r * math.acos(s / (2.0 * r))
            - 0.5 * s * math.sqrt(max(4.0 * r * r - s * s, 0.0)))
    return lens / (math.pi * r * r)


def motor_stations(n_motors: int, span_m: float, fuse_width_m: float,
                   layout: str = "mid", prop_d: float = 0.0,
                   props_in_span: bool = False
                   ) -> Tuple[int, List[float], float]:
    """Where the motors go: (motors on the nose, spanwise stations of the
    wing motors, the span share each wing motor has).

    An odd count puts one motor on the nose - a propeller cannot turn on
    the centreline anywhere else - and every other motor in a wing
    nacelle, the same number each side.  Each wing motor sits in the
    middle of an equal share of the exposed semi-span, from the fuselage
    side to the tip.  So 1 is a nose tractor, 2 a twin on the wings with a
    nose cap, 3 a nose motor and one each side, 4 two each side."""
    n = max(int(n_motors), 1)
    n_nose = n % 2
    per_side = (n - n_nose) // 2
    if per_side == 0:
        return n_nose, [], 0.0
    share = 0.5 * max(span_m - fuse_width_m, 0.0) / per_side
    ys = [0.5 * fuse_width_m + share * (j + 0.5) for j in range(per_side)]
    r = 0.5 * prop_d
    if prop_d > 0.0 and layout in ("tip", "stagger"):
        # (the "share" returned is what the overlap check and the blown
        # bands use: the centre spacing, or twice the inner disc's
        # distance to the fuselage side if that is tighter)
        # the outermost on the wingtip: its axis there, or (if the discs
        # count in the span) its disc's edge
        y_out = 0.5 * span_m - (r if props_in_span else 0.0)
        y_in = (0.5 * fuse_width_m + PROP_GAP_M + r
                if layout == "stagger" else ys[0])
        if per_side == 1:
            ys = [y_out]
            share = 2.0 * (y_out - 0.5 * fuse_width_m)
        else:
            step = (y_out - y_in) / (per_side - 1)
            ys = [y_in + step * j for j in range(per_side)]
            share = (step if layout == "stagger"
                     else min(step, 2.0 * (y_in - 0.5 * fuse_width_m)))
    return n_nose, [-y for y in reversed(ys)] + ys, share


def motor_diameter_m(mass_kg: float) -> float:
    """An outrunner's can diameter from its mass, for the width of the
    nacelle around it: a power law through KDE Direct's published sizes,
    28.3 mm at 42 g (2306XF) and 77.2 mm at 630 g (6815XF) - d grows
    about as m^0.37, a little faster than the m^(1/3) of a scaled copy
    because big motors are flatter."""
    return 0.0283 * (max(mass_kg, 1e-3) / 0.042) ** 0.371


# Copper, for the wing-motor power leads
COPPER_RESISTIVITY_OHM_M = 1.72e-8     # at 20 C
COPPER_DENSITY_KG_M3 = 8960.0
LEAD_JACKET_FACTOR = 1.25              # silicone jacket on the copper's mass


def wing_power_leads(n_motors: int, span_m: float, fuse_width_m: float,
                     i_rated_a: float, build, layout: str = "mid",
                     prop_d: float = 0.0, props_in_span: bool = False
                     ) -> Tuple[float, float, List[float]]:
    """Wing motors' power leads: resistance, mass and lengths.

    Returns (mean lead resistance per motor, ohm; total lead mass, kg;
    each wing motor's lead length, m).

    A pair runs the motor's spanwise station plus
    BuildStandard.lead_route_m, in copper sized to the motor's rated
    current at BuildStandard.lead_current_density_a_mm2: resistance
    2 L rho / A, mass 2 L A density x the jacket.  For a given length the
    two multiply to a constant, so the current density is the trade
    between them.  A nose motor's ESC sits by the battery and its leads
    are the fuse harness.  The resistance is averaged over every motor,
    nose included, because one operating point describes them all
    (PropulsionSystem); the total loss is right to first order.
    """
    n = max(int(n_motors), 1)
    _n_nose, ys, _share = motor_stations(n, span_m, fuse_width_m, layout,
                                         prop_d, props_in_span)
    area = (max(i_rated_a, 1.0)
            / max(build.lead_current_density_a_mm2, 0.1) * 1e-6)   # m^2
    lengths = [abs(y) + max(build.lead_route_m, 0.0) for y in ys]
    r_sum = sum(2.0 * ln * COPPER_RESISTIVITY_OHM_M / area for ln in lengths)
    mass = sum(2.0 * ln * area * COPPER_DENSITY_KG_M3 * LEAD_JACKET_FACTOR
               for ln in lengths)
    return r_sum / n, mass, lengths


def section_perimeters(width: float, height: float, build) -> Tuple[float, float]:
    """(full-section perimeter, tail-post perimeter) of the fuselage.

    The full section is a rectangle with rounded corners, radius
    `fuse_corner_radius_frac` of the smaller side: each corner swaps two
    straight runs of r for a quarter circle, 2(w + h) - (8 - 2 pi) r.  (A
    flat 0.92 on the sharp-cornered perimeter stood in for this.)"""
    r = (min(max(build.fuse_corner_radius_frac, 0.0), 0.5)
         * min(width, height))
    return (2.0 * (width + height) - (8.0 - 2.0 * math.pi) * r,
            math.pi * build.tail_post_diameter_m)


def tailcone_min_length_m(width: float, height: float, build) -> float:
    """Shortest tailcone that tapers from the full section to the tail post
    without exceeding the boat-tail angle at which the flow separates."""
    d_eq = math.sqrt(4.0 * width * height / math.pi)
    taper = max(d_eq - build.tail_post_diameter_m, 0.0)
    return 0.5 * taper / math.tan(math.radians(build.tailcone_max_angle_deg))


def tail_arm_optimum_m(wing, width: float, height: float, build,
                       vol_h=None, vol_v=None) -> float:
    """The tail arm that makes horizontal tail + fin + tailcone lightest.

    With the volume coefficients held, tail areas go as 1/l and the
    tailcone grows one-for-one with l:
        m(l) = (rho_h Vh S c + rho_v Vv S b) / l + lambda l   (c = MAC)
        dm/dl = 0  ->  l* = sqrt((rho_h Vh S c + rho_v Vv S b) / lambda)
    lambda is the tailcone's mass per metre: its mean perimeter times the
    fuselage areal density, plus the tail boom and control links.
    """
    s, c, b = wing.area_m2, wing.mac_m, wing.span_m
    vh = build.vol_coeff_h if vol_h is None else vol_h
    vv = build.vol_coeff_v if vol_v is None else vol_v
    a = (build.h_tail_areal_density * vh * s * c
         + build.v_tail_areal_density * vv * s * b)
    p0, pe = section_perimeters(width, height, build)
    lam = (build.fuse_areal_density * 0.5 * (p0 + pe)
           + build.tail_boom_kg_per_m)
    return math.sqrt(a / max(lam, 1e-9))


def wing_le_root_from_mac(x_le_mac: float, wing) -> float:
    """Root leading-edge station for a straight-tapered wing whose
    quarter-chord line is unswept (as the three-view draws it)."""
    lam = wing.taper
    c_r = wing.root_chord_m
    c_t = c_r * lam
    y_mac = wing.span_m / 6.0 * (1.0 + 2.0 * lam) / (1.0 + lam)
    return x_le_mac - 0.25 * (c_r - c_t) * (2.0 * y_mac / wing.span_m)


# ==========================================================================
# The balance
# ==========================================================================
def solve_layout(ac, cg_minus_ac_m: float) -> Layout:
    """Place every part and balance the aeroplane (rules in the module
    docstring).

    `ac` must already carry its body, wing, mass breakdown, propulsion and
    payload geometry.  `cg_minus_ac_m` is how far the target CG sits aft
    of the wing's aerodynamic centre: (x_np/c - SM_target - 0.25) c.
    """
    b = ac.cfg.build
    body, wing = ac.body, ac.wing
    c = wing.mac_m                  # the reference chord for pitch
    bay_l = ac.bay_length_m
    bay_w, bay_h = ac.bay_width_m, ac.bay_height_m
    nm = max(ac.design.n_motors, 1)
    n_nose, _ys, _share = motor_stations(nm, wing.span_m, body.fuse_width_m)
    f_nose = n_nose / float(nm)          # share of the propulsion on the nose
    f_wing = 1.0 - f_nose
    items = dict(ac.mass_breakdown.items)

    # ---- contingency is spread over the airframe it is contingency for --
    contingency = items.pop("contingency", 0.0)
    base = sum(items.values())
    scale = 1.0 + contingency / max(base, 1e-9)

    def m(key_prefix: str) -> float:
        return scale * sum(v for k, v in items.items() if k.startswith(key_prefix))

    payload_kg = ac.m2_payload_kg
    m_motors = m("motors")
    m_props = m("propellers")
    m_escs = m("ESCs")
    m_mounts = m("motor mounts")
    m_nacelles = m("wing nacelles")
    m_leads = m("wing power leads")
    m_packs = m("propulsion packs")
    m_harness = m("fuse harness")
    m_rxb = m("Rx/servo battery")
    m_avion = m("avionics")
    m_servos = m("servos")
    m_deploy = m("sensor deploy")
    m_restraint = m("payload restraint")
    m_wing_skin = m("wing skin")
    m_spar = m("wing spar")
    m_fuse = m("fuselage")
    m_ht = m("horizontal tail")
    m_vt = m("vertical tail")
    m_boom = m("tail boom")
    m_gear = m("landing gear") + m("tall gear")

    # ---- nose: motor section, then the battery/avionics bay -------------
    d_eq = math.sqrt(4.0 * body.fuse_width_m * body.fuse_height_m / math.pi)
    l_mot = motor_length_m(ac.motor.mass_kg, b)
    if n_nose:
        motor_section = l_mot + b.motor_mount_length_m
    else:
        motor_section = 0.5 * d_eq             # a nose cap, no motor
    pl, pw, _ph = pack_box_m(ac.pack.mass_kg, b)
    across = max(1, int(bay_w // max(pw, 1e-6)))
    rows = int(math.ceil(ac.n_systems / float(across)))
    l_bat = rows * pl
    m_av = m_harness + m_rxb + m_avion + m_escs * f_nose
    l_av = max(m_av / b.avionics_density_kg_m3
               / max(0.5 * bay_w * bay_h, 1e-6), 0.03)
    nose_min = motor_section + l_bat + l_av

    # ---- the parts that move with the payload bay -----------------------
    # Positions relative to the CG target x_p (the bay centroid).
    rel_ac = -cg_minus_ac_m                    # wing AC, relative to x_p
    rel_le_mac = rel_ac - 0.25 * c
    rel_le_root = wing_le_root_from_mac(rel_le_mac, wing)
    c_h = tail_chord_m(body.s_horiz_m2, b)
    rel_tail = rel_ac + body.tail_arm_m
    rel_end = rel_tail + 0.75 * c_h            # fuselage ends at the tail TE
    rel_bay_end = 0.5 * bay_l
    l_tc = max(rel_end - rel_bay_end, 0.0)
    main_gear_rel = (math.tan(math.radians(b.gear_tipback_deg))
                     * (ac.gear_height_m + 0.5 * body.fuse_height_m))
    p0, pe = section_perimeters(body.fuse_width_m, body.fuse_height_m, b)

    def stations(nose: float, x_b: float) -> List[Station]:
        x_p = nose + 0.5 * bay_l
        aft = x_b > nose                       # battery behind the bay
        st: List[Station] = []

        def add(name, mass, x, x0=None, x1=None):
            if mass > 1e-12:
                st.append(Station(name, mass, x,
                                  x if x0 is None else x0,
                                  x if x1 is None else x1))
        # propulsion
        if n_nose:
            add("motor", m_motors * f_nose, 0.5 * l_mot, 0.0, l_mot)
            add("propeller + spinner", (m_props + m_mounts) * f_nose, 0.0)
        if f_wing > 0.0:
            x_nac = x_p + rel_le_root - (0.5 * l_mot + b.motor_mount_length_m)
            # staggered props (4+ motors): every other nacelle sits 0.6 D
            # forward on a boom (Aircraft: 'stagger extensions')
            stag = (ac.design.wing_prop_layout == "stagger"
                    and len(_ys) >= 4)
            l_ext = 0.6 * ac.prop.diameter_m if stag else 0.0
            for frac, dx_ in (((0.5, 0.0), (0.5, -l_ext)) if stag
                              else ((1.0, 0.0),)):
                x_n = x_nac + dx_
                add("motors (on the wing)",
                    (m_motors + m_mounts) * f_wing * frac, x_n)
                add("propellers", m_props * f_wing * frac, x_n - 0.5 * l_mot)
                add("ESCs (in the nacelles)", m_escs * f_wing * frac, x_n)
                add("wing nacelles", m_nacelles * frac, x_n)
            if stag:
                add("stagger extensions", m("stagger extensions"),
                    x_nac - 0.5 * l_ext)
            # the leads run from the battery to the nacelles (with
            # stagger, half of them to the forward ones)
            if stag:
                add("wing power leads", 0.5 * m_leads, 0.5 * (x_b + x_nac))
                add("wing power leads", 0.5 * m_leads,
                    0.5 * (x_b + x_nac - l_ext))
            else:
                add("wing power leads", m_leads, 0.5 * (x_b + x_nac))
        x_nose_gear = motor_section
        # battery and avionics share the bay between motor and payload
        add("propulsion battery", m_packs, x_b,
            x_b - 0.5 * l_bat, x_b + 0.5 * l_bat)
        if aft:                                 # nose holds avionics only
            x_av = 0.5 * (motor_section + nose)
        else:
            g0 = (motor_section, x_b - 0.5 * l_bat)      # gap ahead
            g1 = (x_b + 0.5 * l_bat, nose)               # gap behind
            la, lb = max(g0[1] - g0[0], 0.0), max(g1[1] - g1[0], 0.0)
            x_av = ((0.5 * (g0[0] + g0[1]) * la
                     + 0.5 * (g1[0] + g1[1]) * lb) / max(la + lb, 1e-9))
        add("ESC, Rx, Rx battery, fuse harness", m_av, x_av)
        # payload bay
        add("payload (containers + simulators)", payload_kg, x_p,
            nose, nose + bay_l)
        add("sensor deploy mechanism + restraint", m_deploy + m_restraint,
            x_p)
        add("bay servos (door, sensor)", m_servos * 2.0 / 6.0, x_p)
        # wing
        x_le = x_p + rel_le_mac
        add("wing skin/ribs", m_wing_skin, x_le + 0.42 * c, x_le, x_le + c)
        add("wing spar", m_spar, x_le + 0.28 * c)
        add("aileron servos", m_servos * 2.0 / 6.0, x_le + 0.70 * c)
        # tail
        x_t = x_p + rel_tail
        add("horizontal tail", m_ht, x_t + 0.15 * c_h)
        add("vertical tail", m_vt, x_t + 0.15 * c_h)
        x_tc0 = nose + bay_l
        add("tail boom/links", m_boom, x_tc0 + 0.5 * l_tc)
        add("tail servos", m_servos * 2.0 / 6.0, x_tc0 + 0.5 * l_tc)
        # fuselage shell: constant section to the bay end, then a linear
        # taper to the tail post - the area-weighted centroid of that skin
        a1 = p0 * x_tc0
        a2 = 0.5 * (p0 + pe) * l_tc
        x2 = x_tc0 + l_tc * (p0 + 2.0 * pe) / (3.0 * max(p0 + pe, 1e-9))
        x_f = (a1 * 0.5 * x_tc0 + a2 * x2) / max(a1 + a2, 1e-12)
        add("fuselage shell", m_fuse, x_f, 0.0, x_tc0 + l_tc)
        # gear: tricycle, mains behind the CG at the tip-back angle
        add("main gear", 0.8 * m_gear, x_p + main_gear_rel)
        add("nose gear", 0.2 * m_gear, x_nose_gear)
        return st

    def moment(nose: float, x_b: float) -> float:
        x_p = nose + 0.5 * bay_l
        return sum(s.mass_kg * (s.x_m - x_p) for s in stations(nose, x_b))

    # Where a battery behind the bay still fits: the tailcone tapers at the
    # boat-tail angle, and the pack's cross-section diagonal has to fit.
    _pl, pw_, ph_ = pack_box_m(ac.pack.mass_kg, b)
    diag = math.hypot(pw_, ph_)
    s_fit = ((d_eq - diag)
             / (2.0 * math.tan(math.radians(b.tailcone_max_angle_deg))))
    aft_ok = s_fit >= l_bat and l_tc >= l_bat

    def xb_fwd(nose: float) -> Tuple[float, float]:
        """Battery in the nose: right behind the motor .. against the bay."""
        return motor_section + 0.5 * l_bat, nose - 0.5 * l_bat

    def xb_aft(nose: float) -> Tuple[float, float]:
        """Battery behind the bay, where the tailcone still holds it."""
        x0 = nose + bay_l + 0.5 * l_bat
        return x0, nose + bay_l + min(s_fit, l_tc) - 0.5 * l_bat

    def lengthen(xb_of) -> Tuple[float, float, bool]:
        """Tail-heavy with the battery at its forward limit: lengthen the
        nose until it balances there."""
        n0 = motor_section + l_av + (0.0 if xb_of is xb_aft else l_bat)
        n_lo, n_hi = n0, n0 + 5.0
        if moment(n_hi, xb_of(n_hi)[0]) > 0.0:
            return n_hi, xb_of(n_hi)[0], False
        for _ in range(60):
            mid = 0.5 * (n_lo + n_hi)
            if moment(mid, xb_of(mid)[0]) > 0.0:
                n_lo = mid
            else:
                n_hi = mid
        return n_hi, xb_of(n_hi)[0], True

    def slide(nose: float, lo: float, hi: float) -> float:
        a_, b_ = lo, hi
        for _ in range(60):
            mid = 0.5 * (a_ + b_)
            if moment(nose, mid) > 0.0:
                b_ = mid
            else:
                a_ = mid
        return 0.5 * (a_ + b_)

    def margin(nose: float, x_b: float) -> float:
        """Static margin this arrangement actually gives (Mission 2)."""
        st_ = stations(nose, x_b)
        x_p_ = nose + 0.5 * bay_l
        x_cg_ = (sum(q.mass_kg * q.x_m for q in st_)
                 / sum(q.mass_kg for q in st_))
        return b.target_static_margin + (x_p_ - x_cg_) / c

    # Candidate arrangements.  moment > 0 means the CG is aft of target
    # (tail-heavy).  Each is (nose, battery x, note).
    cands = []
    nose = nose_min
    lo, hi = xb_fwd(nose)
    if moment(nose, lo) > 0.0:
        # tail-heavy even with the battery right behind the motor
        n_, xb_, ok_ = lengthen(xb_fwd)
        cands.append((n_, xb_, "nose lengthened to %.0f mm so the battery "
                      "can balance the tail" % (1000 * n_) if ok_ else
                      "tail-heavy even with a 5 m nose - check the tails"))
    elif moment(nose, hi) > 0.0:
        cands.append((nose, slide(nose, lo, hi),
                      "balanced by the battery position in the nose"))
    else:
        # nose-heavy with the battery against the bay: keep the shortest
        # nose as one option, and try the battery behind the bay
        cands.append((nose, hi, "battery against the payload bay, nose as "
                      "short as its parts allow; CG slightly forward of "
                      "target"))
        if aft_ok:
            n_b = motor_section + l_av
            lo_b, hi_b = xb_aft(n_b)
            if moment(n_b, lo_b) > 0.0:
                n_, xb_, ok_ = lengthen(xb_aft)
                cands.append((n_, xb_, "battery behind the payload bay, "
                              "nose lengthened to %.0f mm to balance"
                              % (1000 * n_)))
            elif moment(n_b, hi_b) < 0.0:
                cands.append((n_b, hi_b, "battery as far aft as the "
                              "tailcone holds it; CG forward of target"))
            else:
                cands.append((n_b, slide(n_b, lo_b, hi_b), "balanced by "
                              "the battery, placed behind the payload bay"))

    # Choose: the shortest fuselage whose margin is inside the band
    # [target, target + static_margin_band]; if none is, the one closest
    # to the band.  A little extra margin is better than a longer, heavier
    # fuselage built only to hit the target exactly.
    t_lo = b.target_static_margin - 1e-6
    t_hi = b.target_static_margin + b.static_margin_band + 1e-6
    scored = []
    for n_, xb_, note_ in cands:
        sm_ = margin(n_, xb_)
        miss = max(t_lo - sm_, sm_ - t_hi, 0.0)
        scored.append((miss, n_, xb_, note_, sm_))
    scored.sort(key=lambda r: (r[0], r[1]))
    _miss, nose, x_b, note, _sm = scored[0]
    balanced = _miss <= 0.0

    st = stations(nose, x_b)
    x_p = nose + 0.5 * bay_l
    mtot = sum(s.mass_kg for s in st)
    x_cg = sum(s.mass_kg * s.x_m for s in st) / mtot
    empty = [s for s in st if not s.name.startswith("payload")]
    m_e = sum(s.mass_kg for s in empty)
    x_cg_e = sum(s.mass_kg * s.x_m for s in empty) / max(m_e, 1e-9)
    x_cg_3 = ((m_e * x_cg_e + ac.m3_payload_kg * x_p)
              / max(m_e + ac.m3_payload_kg, 1e-9))
    sm_t = b.target_static_margin
    x_le_mac = x_p + rel_le_mac
    x_ac = x_p + rel_ac
    return Layout(
        stations=st, nose_length_m=nose, motor_section_m=motor_section,
        battery_length_m=l_bat, avionics_length_m=l_av,
        bay_start_m=nose, bay_end_m=nose + bay_l, tailcone_length_m=l_tc,
        fuse_length_m=nose + bay_l + l_tc,
        x_wing_le_mac_m=x_le_mac,
        x_wing_le_root_m=x_p + rel_le_root, x_wing_ac_m=x_ac,
        x_tail_ac_m=x_p + rel_tail, tail_arm_m=body.tail_arm_m,
        x_cg_target_m=x_p, x_cg_m=x_cg, x_cg_empty_m=x_cg_e, x_cg_m3_m=x_cg_3,
        x_battery_m=x_b, x_main_gear_m=x_p + main_gear_rel,
        x_nose_gear_m=motor_section,
        static_margin=sm_t + (x_p - x_cg) / c,
        static_margin_empty=sm_t + (x_p - x_cg_e) / c,
        balanced=balanced, note=note,
        moments={s.name: s.mass_kg * (s.x_m - x_p) for s in st})
