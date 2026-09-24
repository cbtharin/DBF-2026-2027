"""Mass buildup.

Structure mass is where student optimisers usually cheat, and it is exactly
the term that decides both missions.  So the spar is sized from the actual
root bending moment at the ultimate load factor, not from a wave of the hand,
and everything else is an areal density you can calibrate against last year's
airframe by weighing the parts.

Calibration hook: `BuildStandard.wing_areal_density` etc. in config.py.
Weigh your existing wing, divide by planform area, put that number in.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict

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
    # Split, so the elevator and the rudder can be calibrated separately.
    mb.add("horizontal tail (elevator)",
           build.h_tail_areal_density * body.s_horiz_m2)
    mb.add("vertical tail (rudder)",
           build.v_tail_areal_density * body.s_vert_m2)
    mb.add("tail boom/links", 0.030 + 0.045 * body.tail_arm_m)
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
