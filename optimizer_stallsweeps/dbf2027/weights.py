"""Mass buildup.

Structure mass is where student optimisers usually cheat, and it is exactly
the term that decides both missions.  So the spar is sized station by
station along the span for every wing load case - flight at the ultimate
load factor, the tip-lift test, shear and torsion, each with the same 1.5
factor of safety - not from a wave of the hand,
and everything else is an areal density you can calibrate against last year's
airframe by weighing the parts.

Calibration hook: `BuildStandard.wing_areal_density` etc. in config.py.
Weigh your existing wing, divide by planform area, put that number in.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

from .aerodynamics import BodyGeometry, WingGeometry
from .config import BuildStandard


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


#: Stations per half-wing at which the spar is sized, root to tip.  The
#: mass converges to well under a gram by 24 (checked in validate.py).
SPAR_STATIONS = 24


def _tip_angle(eta: float) -> float:
    """theta, with eta = cos(theta), taken from the station's distance to
    the tip: 2 asin(sqrt((1 - eta) / 2)) keeps full relative accuracy as
    eta -> 1, where acos(eta) loses it (its slope is infinite there).
    For eta >= 1/2, 1 - eta is itself exact in floating point."""
    e = min(max(eta, 0.0), 1.0)
    return 2.0 * math.asin(math.sqrt(0.5 * (1.0 - e)))


def _elliptic_moment(eta: float) -> float:
    """Bending moment of an elliptic half-wing lift distribution at
    eta = y/(b/2), as a fraction of n*W*b/pi.  4/(3 pi) of (W/2)(b/2) at
    the root, falling to zero at the tip.

    With eta = cos(theta) it is f/6, f = (9/4) sin(theta) +
    (1/4) sin(3 theta) - 3 theta cos(theta) - the same quantity the
    textbook form (1 - eta^2)^1.5 / 3 - eta (pi/4 - (eta sqrt(1 - eta^2)
    + asin eta) / 2) gives, but that form takes the difference of two
    terms that both tend to theta^3 / 3 at the tip while their difference
    goes as theta^5 / 15, and loses its digits there.  Near the tip f is
    summed from its series, whose theta and theta^3 terms vanish exactly:
        f = sum_{k>=2} (-1)^k [(9 + 3^(2k+1))/4 - 3(2k+1)]
                                theta^(2k+1) / (2k+1)!
    """
    th = _tip_angle(eta)
    if th >= 0.5:
        return (2.25 * math.sin(th) + 0.25 * math.sin(3.0 * th)
                - 3.0 * th * math.cos(th)) / 6.0
    t2 = th * th
    power = th * t2 * t2                      # theta^5
    fact = 120.0                              # 5!
    p3 = 243.0                                # 3^5
    total, sign = 0.0, 1.0
    for k in range(2, 30):
        term = sign * ((9.0 + p3) / 4.0 - 3.0 * (2 * k + 1)) * power / fact
        total += term
        if abs(term) <= 1e-17 * abs(total):
            break
        power *= t2
        fact *= (2 * k + 2) * (2 * k + 3)
        p3 *= 9.0
        sign = -sign
    return total / 6.0


def _x_minus_sin(x: float) -> float:
    """x - sin(x), summed from its series below x = 1, where it is the
    difference of two nearly equal numbers."""
    if x >= 1.0:
        return x - math.sin(x)
    x2 = x * x
    term = x * x2 / 6.0
    total = 0.0
    k = 1
    while abs(term) > 1e-17 * max(abs(total), 1e-300):
        total += term
        term *= -x2 / ((2 * k + 2) * (2 * k + 3))
        k += 1
    return total


def _elliptic_shear(eta: float) -> float:
    """Shear of the same distribution, as a fraction of its root value:
    (2 theta - sin 2 theta) / pi with eta = cos(theta).  (The textbook
    1 - (2/pi)(eta sqrt(1 - eta^2) + asin eta) is the same number, taken
    as 1 minus something that tends to 1 at the tip.)"""
    return _x_minus_sin(2.0 * _tip_angle(eta)) / math.pi


def _equivalent_span_ratio(wing: WingGeometry, build: BuildStandard) -> float:
    """The span of the wing whose induced drag the end-plated one has,
    over the real span: sqrt(AR_eff / AR) = sqrt(1 + gain h/b) (m110),
    1 with no winglets."""
    h = max(wing.winglet_h_m, 0.0)
    return math.sqrt(1.0 + build.winglet_ar_gain * h / wing.span_m) if h > 0.0 else 1.0


def _ellipse_area(x: float) -> float:
    """integral_0^x sqrt(1 - u^2) du."""
    x = min(max(x, 0.0), 1.0)
    return 0.5 * (x * math.sqrt(max(1.0 - x * x, 0.0)) + math.asin(x))


def _end_plated_moment(eta: float, r: float) -> float:
    """Bending moment at eta = y/(b/2) of an end-plated wing's lift, as a
    fraction of n W b / pi, as _elliptic_moment (m117): the inboard part of
    an ellipse over the equivalent span r b, truncated at the real tip and
    carrying the whole lift.  r = 1 is the elliptic wing."""
    if r <= 1.0 + 1e-12:
        return _elliptic_moment(eta)
    u0, u1 = min(max(eta, 0.0), 1.0) / r, 1.0 / r
    b0 = (max(1.0 - u0 * u0, 0.0)) ** 1.5 / 3.0     # -B(u0)
    b1 = (max(1.0 - u1 * u1, 0.0)) ** 1.5 / 3.0     # -B(u1)
    a0, a1 = _ellipse_area(u0), _ellipse_area(u1)
    return max(math.pi * r * (b0 - b1 - u0 * (a1 - a0)) / (4.0 * a1), 0.0)


def _end_plated_shear(eta: float, r: float) -> float:
    """Shear of the same distribution as a fraction of its root value, as
    _elliptic_shear (m117)."""
    if r <= 1.0 + 1e-12:
        return _elliptic_shear(eta)
    a1 = _ellipse_area(1.0 / r)
    return max((a1 - _ellipse_area(min(max(eta, 0.0), 1.0) / r)) / a1, 0.0)


@dataclass
class SparDesign:
    """A carbon tube spar sized station by station along the half-span.

    Every wing load case carries the same factor of safety
    (`BuildStandard.load_safety_factor`, 1.5):

      flight bending  n_ult = FoS x n_limit, elliptic lift - the V-n case;
                      with winglets the inboard part of an ellipse over the
                      equivalent span, which loads the tips (m117)
      wing-tip test   RULE 4.2.4: the aeroplane lifted by its tips at
                      MGTOW, a tip reaction of W/2 each side.  Its moment
                      falls off linearly, not elliptically, so it governs
                      the OUTBOARD spar even where flight governs the root
      shear+torsion   flight or tip-test shear, plus the section pitching
                      moment at the dive speed, all carried by the tube

    `fos` is the factor actually achieved against each case's UNFACTORED
    load, at the station where it is smallest.  It is never below the
    requirement when the spar closes, and above it wherever minimum gauge
    rather than load sets the wall.
    """
    mass_kg: float
    t_req_m: float              # thickest wall any load case asks for
    t_used_m: float             # thickest wall used (>= minimum gauge)
    gauge_limited: bool         # minimum gauge governs at every station
    closes: bool                # the tube fits its wall at every station
    r_root_m: float
    r_tip_m: float
    fos_required: float
    fos: Dict[str, float] = field(default_factory=dict)
    fos_at_eta: Dict[str, float] = field(default_factory=dict)
    worst_eta: float = 0.0      # where the tube comes closest to not closing


def spar_design(wing: WingGeometry, mtow_n: float, build: BuildStandard,
                n_ultimate: float, q_dive_pa: float = 0.0,
                point_masses: Sequence[Tuple[float, float]] = ()
                ) -> SparDesign:
    """Size the spar for every wing load case, at every station.

    `n_ultimate` comes from the V-n diagram in loads.py, not a constant.
    `mtow_n` is the heaviest configuration (Mission 2), which is also the
    MGTOW the tip test is run at.  `q_dive_pa` is the dynamic pressure at
    the design dive speed, for the torsion case; 0 leaves torsion out.

    Thin-wall tube: section modulus pi r^2 t, shear stress V/(pi r t) at
    the neutral axis, torsional shear T/(2 pi r^2 t).  The tube radius
    follows the aerofoil depth (0.40 x local thickness, never under 6 mm),
    so on a tapered wing the tube tapers too; the wall at each station is
    the thickest any case needs, and never below minimum gauge.
    Conservative throughout: no inertia relief from the wing's own weight.
    The wing motors' is counted: `point_masses` is (y, kg) for each pod on
    one half-wing (motor, propeller, ESC, mount and nacelle), and its
    weight - times the load factor in flight, at 1 g in the tip-lift test -
    takes off the shear and bending inboard of it.  Leaving it out charged
    a wing-motor layout for spar it does not need.
    """
    fos = max(build.load_safety_factor, 1.0)
    w = mtow_n
    half = 0.5 * wing.span_m
    lam = min(max(wing.taper, 0.05), 1.0)
    c_root = wing.root_chord_m
    tc = wing.airfoil.t_over_c
    cm = abs(wing.airfoil.cm_ac)
    sig, tau = build.spar_sigma_allow_pa, build.spar_tau_allow_pa
    t_min = build.spar_min_wall_m
    n_lim = n_ultimate / fos

    def chord(e: float) -> float:
        return c_root * (1.0 - (1.0 - lam) * e)

    def radius(e: float) -> float:
        return max(0.40 * tc * chord(e), 0.006)

    def torque(e: float) -> float:
        # section pitching moment outboard of e, unfactored: q c^2 cm per
        # unit span, with the chord linear in e
        c0, c1 = chord(e), chord(1.0)
        integ = (1.0 - e) * (c0 * c0 + c0 * c1 + c1 * c1) / 3.0
        return q_dive_pa * cm * integ * half

    pods = [(abs(y), m) for y, m in point_masses]
    r_eq = _equivalent_span_ratio(wing, build)
    n = SPAR_STATIONS
    r_list: List[float] = []
    t_list: List[float] = []
    t_req_max, t_used_max = 0.0, 0.0
    closes, worst_eta, worst_ratio = True, 0.0, 0.0
    fos_got = {"flight bending": float("inf"), "wing-tip test": float("inf"),
               "shear + torsion": float("inf")}
    fos_eta = {k: 0.0 for k in fos_got}
    for i in range(n + 1):
        e = i / float(n)
        r = radius(e)
        # unfactored loads at this station, less the weight of the wing
        # pods outboard of it (their shear and their moment about it)
        y = e * half
        g_out = 9.80665 * sum(m for yi, m in pods if yi > y)
        gm_out = 9.80665 * sum(m * (yi - y) for yi, m in pods if yi > y)
        m_flight = abs(n_lim * w * wing.span_m / math.pi
                       * _end_plated_moment(e, r_eq) - n_lim * gm_out)
        m_tip = abs(0.5 * w * half * (1.0 - e) - gm_out)     # RULE 4.2.4
        v_flight = abs(n_lim * 0.5 * w * _end_plated_shear(e, r_eq)
                       - n_lim * g_out)
        v_tip = abs(0.5 * w - g_out) if e < 1.0 else 0.0
        # shear flow per unit wall, i.e. stress x wall
        shear = (max(v_flight, v_tip) / (math.pi * r)
                 + torque(e) / (2.0 * math.pi * r * r))
        t_req = fos * max(max(m_flight, m_tip) / (sig * math.pi * r * r),
                          shear / tau)
        t = max(t_req, t_min)
        r_list.append(r)
        t_list.append(t)
        t_req_max = max(t_req_max, t_req)
        t_used_max = max(t_used_max, t)
        if t / r > worst_ratio:
            worst_ratio, worst_eta = t / r, e
        if t >= r:
            closes = False
        # the factor of safety this wall actually gives each case
        cap_m = sig * math.pi * r * r * t
        for key, got in (("flight bending", cap_m / m_flight
                          if m_flight > 1e-9 else float("inf")),
                         ("wing-tip test", cap_m / m_tip
                          if m_tip > 1e-9 else float("inf")),
                         ("shear + torsion", tau * t / shear
                          if shear > 1e-12 else float("inf"))):
            if got < fos_got[key]:
                fos_got[key], fos_eta[key] = got, e

    # two half-wings, trapezoidal rule along the span
    dy = half / n
    per_half = sum(0.5 * dy * 2.0 * math.pi
                   * (r_list[i] * t_list[i] + r_list[i + 1] * t_list[i + 1])
                   for i in range(n))
    mass = 2.0 * per_half * build.spar_rho * build.spar_fitting_factor
    return SparDesign(mass, t_req_max, t_used_max, t_req_max < t_min, closes,
                      r_list[0], r_list[-1], fos, fos_got, fos_eta,
                      worst_eta)


def spar_sizing(wing: WingGeometry, mtow_n: float, build: BuildStandard,
                n_ultimate: float, q_dive_pa: float = 0.0,
                point_masses: Sequence[Tuple[float, float]] = ()):
    """Returns (mass_kg, t_required_m, t_used_m, gauge_limited).

    `gauge_limited` is True when every load case, everywhere along the
    span, asks for less wall than you can actually lay up or buy.  That is
    common on this class of aeroplane and worth knowing: it means the spar
    is sized by buildability, not by load, so raising the load factor
    costs nothing and lowering it saves nothing.
    """
    sd = spar_design(wing, mtow_n, build, n_ultimate, q_dive_pa,
                     point_masses)
    return sd.mass_kg, sd.t_req_m, sd.t_used_m, sd.gauge_limited


def spar_mass(wing: WingGeometry, mtow_n: float, build: BuildStandard,
              n_ultimate: float, q_dive_pa: float = 0.0,
              point_masses: Sequence[Tuple[float, float]] = ()) -> float:
    return spar_design(wing, mtow_n, build, n_ultimate, q_dive_pa,
                       point_masses).mass_kg


def airframe_masses(wing: WingGeometry, body: BodyGeometry, mtow_n: float,
                    build: BuildStandard, n_ultimate: float,
                    q_dive_pa: float = 0.0,
                    point_masses: Sequence[Tuple[float, float]] = ()
                    ) -> MassBreakdown:
    mb = MassBreakdown()
    mb.add("wing skin/ribs", build.wing_areal_density * wing.area_m2)
    mb.add("wing spar", spar_mass(wing, mtow_n, build, n_ultimate,
                                  q_dive_pa, point_masses))
    if wing.winglet_h_m > 0.0:              # m110 (placed with the wing skin)
        mb.add("wing skin (winglets)",
               2.0 * (build.winglet_areal_density * wing.winglet_area_m2
                      + build.winglet_fixed_kg))
    mb.add("fuselage", build.fuse_areal_density * body.fuse_wetted_m2)
    # Split, so the elevator and the rudder can be calibrated separately.
    mb.add("horizontal tail (elevator)",
           build.h_tail_areal_density * body.s_horiz_m2)
    mb.add("vertical tail (rudder)",
           build.v_tail_areal_density * body.s_vert_m2)
    mb.add("tail boom/links", build.tail_boom_fixed_kg
           + build.tail_boom_kg_per_m * body.tail_arm_m)
    mb.add("landing gear",
           build.gear_mass_fixed_kg + build.gear_mass_frac * (mtow_n / 9.80665))
    mb.add("servos", build.n_servos_base * build.servo_mass_kg)
    mb.add("avionics/wiring", build.avionics_kg)     # RULE 3.2.3b
    # RULE 3.1.2b - permanently installed for ALL missions, so its mass
    # is carried even in Mission 1, which has no payload at all.
    mb.add("sensor deploy mechanism", build.deploy_mech_kg)
    return mb


def converge_empty_mass(wing: WingGeometry, body: BodyGeometry,
                        propulsion_items: dict, payload_mass_kg: float,
                        build: BuildStandard, n_ultimate: float,
                        tol: float = 1e-10, maxiter: int = 60,
                        q_dive_pa: float = 0.0,
                        point_masses: Sequence[Tuple[float, float]] = ()):
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
                             n_ultimate, q_dive_pa, point_masses)
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
