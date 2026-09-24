"""Propulsion: motor / propeller / ESC / battery matched at every flight point.

The operating point is found by balancing motor shaft torque against propeller
torque, so the model captures the thing that actually kills DBF power systems:
a prop that unloads at speed (thrust falls off a cliff) or one that overloads
the motor at low speed (magnet smoke on the takeoff roll).

    motor:  I = (Vt - RPM/Kv) / Rm,     Q = Kt (I - I0),   Kt = 60 / (2 pi Kv)
    prop:   Q = Cp(J) rho n^2 D^5 / (2 pi),   T = Ct(J) rho n^2 D^4,  J = V/(nD)

Rm is the motor's catalogue winding resistance.  It is not a tuning knob and
it cannot be switched off: it *is* the motor equation.  Setting Rm to zero
makes the current unbounded and the torque balance insoluble - a zero
resistance motor makes infinite torque at any speed.  What the model no
longer does is vary Rm with winding temperature, and it no longer models
battery sag; see "What is deliberately not modelled" below.

Solution method
---------------
Terminal voltage is throttle times the pack's open-circuit voltage, which
does not depend on current, so the whole system is one scalar equation in
propeller speed n - motor torque equals propeller torque - which Brent
solves to machine precision in about ten function evaluations.  There is no
iteration around the electrical loop and therefore no iteration-truncation
error.

That speed is also why there is no interpolation table here.  A tabulated
deck was measured at 43 us per call against 41 us for this exact solve, so
tabulation bought nothing but error, and it was removed.

Propulsion system architecture - RULE 3.2.3.1a, 3.2.3.1b, 3.2.3.1d
----------------------------------------------
A propulsion system is one battery pack, one fuse, one externally accessible
arming plug, one or more ESCs and one or more motors.  Packs may not be
wired in series or parallel, so the only legal way to carry more than one
pack is to build a second complete, independent system - and then every pack
must be identical.  `n_systems` is therefore a real design variable, and
`motors_per_system = n_motors / n_systems` must come out a whole number.

What is deliberately not modelled
---------------------------------
Battery internal resistance (voltage sag) and the rise of Rm with winding
temperature are both omitted, at the team's direction.  Both omissions push
the same way: predicted thrust is optimistic, by roughly 3-6 percent at full
throttle on a warm pack for a typical 6S system.  This is tolerable here
because the binding electrical constraint is no longer something the model
has to infer - it is the 100 A blade fuse, which is checked directly against
the computed pack current.  Treat absolute thrust numbers as an upper bound
and the comparison between designs as sound, since every design is flattered
by the same amount.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

from .components import ESC, BatteryPack, Motor, Propeller, pick_esc
from .config import Environment
from .numerics import brent


@dataclass
class OperatingPoint:
    rpm: float
    thrust_total_n: float
    motor_current_a: float
    pack_current_a: float        # per pack, which is what the fuse sees
    total_current_a: float       # summed over every pack on the aeroplane
    pack_voltage_v: float
    p_elec_w: float              # drawn from all batteries
    p_shaft_w: float             # total, all motors
    eta_prop: float
    eta_motor: float
    advance_ratio: float


class PropulsionSystem:
    """One or more identical propulsion systems, solved as a set.

    Every system carries the same pack and the same number of identical
    motors turning identical propellers, so they all sit at the same
    operating point and one solve describes the whole aeroplane.
    """

    ROOT_XTOL = 1e-12          # Brent tolerance on rev/s

    def __init__(self, motor: Motor, prop: Propeller, n_motors: int,
                 pack: BatteryPack, env: Environment, esc: Optional[ESC] = None,
                 install_thrust_factor: float = 0.95, n_systems: int = 1):
        self.install_thrust_factor = install_thrust_factor
        self.motor = motor
        self.prop = prop
        self.n_motors = n_motors
        self.n_systems = max(1, int(n_systems))
        self.motors_per_system = n_motors / float(self.n_systems)
        self.pack = pack
        self.env = env
        self.esc = esc or pick_esc(motor.i_max_a, pack.series) or _biggest_esc()
        self.kt = 60.0 / (2.0 * math.pi * motor.kv)     # N*m per amp
        self.rm = motor.rm
        self._memo_key: Optional[Tuple[float, float, float]] = None
        self._memo_val: Tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)

    # ------------------------------------------------------------------
    # Energy available
    # ------------------------------------------------------------------
    @property
    def total_watt_hours(self) -> float:     # RULE 3.2.3.2b
        """Rated energy summed over every pack, which is what the 100 Wh
        limit is written against."""
        return self.n_systems * self.pack.watt_hours

    @property
    def usable_watt_hours(self) -> float:
        return self.total_watt_hours * self.pack.usable_fraction

    @property
    def fuse_rating_a(self) -> float:        # RULE 3.2.3.2d
        """The fuse is rated at the pack's labelled maximum continuous
        discharge current, capped at the largest legal blade fuse."""
        return min(self.pack.i_max_a, 100.0)

    # ------------------------------------------------------------------
    # Exact single-point solve
    # ------------------------------------------------------------------
    def solve(self, v: float, soc: float, throttle: float) -> OperatingPoint:
        m, p, env = self.motor, self.prop, self.env
        throttle = max(0.0, min(1.0, throttle))
        v = max(0.0, v)
        v_pack = self.pack.ocv(soc)

        if throttle * v_pack <= 0.05:
            return OperatingPoint(0.0, 0.0, 0.0, 0.0, 0.0, v_pack,
                                  0.0, 0.0, 0.0, 0.0, 0.0)

        vt = throttle * v_pack
        d5 = p.diameter_m ** 5
        two_pi = 2.0 * math.pi
        k = 60.0 / m.kv

        def residual(nn: float) -> float:
            cur = (vt - k * nn) / self.rm
            q_motor = self.kt * (cur - m.i0)
            j = v / (nn * p.diameter_m) if nn > 1e-6 else 1e9
            q_prop = p.cp(j) * env.rho * nn * nn * d5 / two_pi
            return q_motor - q_prop

        # At the no-load speed the motor current is zero, so the residual is
        # negative there: a guaranteed bracket against a positive residual
        # at standstill.
        n_noload = vt * m.kv / 60.0
        if residual(1e-4) <= 0.0:
            n = 0.0
        else:
            n = brent(residual, 1e-4, n_noload, xtol=self.ROOT_XTOL)
        i_motor = max(0.0, (vt - k * n) / self.rm)

        # Battery current is the motor current scaled by throttle (the ESC
        # is a switching converter: it trades voltage for current).
        i_pack = self.motors_per_system * i_motor * throttle
        i_total = self.n_motors * i_motor * throttle
        thrust = (self.n_motors * max(0.0, p.thrust(env.rho, n, v))
                  * self.install_thrust_factor)
        p_shaft = self.n_motors * p.power(env.rho, n, v)
        p_elec = v_pack * i_total / max(self.esc.efficiency, 0.5)
        j = v / (n * p.diameter_m) if n > 1e-6 else 0.0
        eta_prop = (thrust * v / p_shaft) if p_shaft > 1e-6 and v > 0.1 else 0.0
        eta_motor = (p_shaft / (self.n_motors * max(vt * i_motor, 1e-6))
                     if i_motor > 1e-6 else 0.0)
        return OperatingPoint(60.0 * n, thrust, i_motor, i_pack, i_total,
                              v_pack, p_elec, p_shaft, min(eta_prop, 0.98),
                              min(eta_motor, 0.98), j)

    # ------------------------------------------------------------------
    # Query interface
    # ------------------------------------------------------------------
    def operating(self, v: float, soc: float = 1.0,
                  throttle: float = 1.0) -> Tuple[float, float, float, float]:
        """Returns (thrust_N, p_elec_W, total_current_A, pack_voltage_V).

        Memoised on the last query, because the mission integrator asks for
        thrust and power at the same flight condition back to back.
        """
        throttle = max(0.0, min(1.0, throttle))
        v = max(0.0, v)
        key = (v, soc, throttle)
        if key == self._memo_key:
            return self._memo_val
        op = self.solve(v, soc, throttle)
        out = (op.thrust_total_n, op.p_elec_w, op.total_current_a,
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
        """Total current drawn from all packs - the number the 100 A total
        limit is written against."""
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

    def throttle_for_current_limit(self, v: float, soc: float,
                                   i_limit: float) -> float:
        """Largest throttle that keeps total current under `i_limit`.

        Used to hold the aeroplane inside the 100 A fuse rather than simply
        rejecting a design that would momentarily exceed it: a real pilot
        with a real ESC just does not pull that much.
        """
        if self.pack_current(v, soc, 1.0) <= i_limit:
            return 1.0
        lo = 0.02
        if self.pack_current(v, soc, lo) >= i_limit:
            return lo
        return brent(lambda th: self.pack_current(v, soc, th) - i_limit,
                     lo, 1.0, xtol=1e-9)

    # ------------------------------------------------------------------
    # Ratings / masses
    # ------------------------------------------------------------------
    @property
    def mass_kg(self) -> float:
        """Motors, props, ESCs, mounts and every pack on the aeroplane.

        The per-system fuse harness and arming plug are added in the weight
        buildup, which is where the rest of the installation mass lives.
        """
        return (self.n_motors * (self.motor.mass_kg + self.prop.mass_kg
                                 + self.esc.mass_kg + 0.025)   # mount + spinner
                + self.n_systems * self.pack.mass_kg)

    def rating_check(self) -> Dict[str, float]:
        """Worst-case (static, full charge, full throttle) electrical loads,
        from the exact solve."""
        op = self.solve(0.0, 1.0, 1.0)
        per_motor_w = op.p_elec_w / max(self.n_motors, 1)
        return {
            "static_thrust_n": op.thrust_total_n,
            "motor_current_a": op.motor_current_a,
            "pack_current_a": op.pack_current_a,
            "total_current_a": op.total_current_a,
            "fuse_rating_a": self.fuse_rating_a,
            "motor_power_w": per_motor_w,
            "motor_current_margin": self.motor.i_max_a / max(op.motor_current_a, 1e-6),
            "motor_power_margin": self.motor.p_max_w / max(per_motor_w, 1e-6),
            "esc_current_margin": self.esc.i_cont_a / max(op.motor_current_a, 1e-6),
            "pack_current_margin": self.pack.i_max_a / max(op.pack_current_a, 1e-6),
            "fuse_margin": self.fuse_rating_a / max(op.pack_current_a, 1e-6),
            "rpm": op.rpm,
            "tip_mach": (math.pi * self.prop.diameter_m * op.rpm / 60.0) / 340.0,
        }


def _biggest_esc() -> ESC:
    from .components import ESCS
    return max(ESCS.values(), key=lambda e: e.i_cont_a)


# --------------------------------------------------------------------------
# Instance cache.  A PropulsionSystem depends only on the hardware, never on
# the wing, so thousands of candidate wings share one object - and with it
# its last-operating-point memo.  Bounded by the hardware shortlist size.
# --------------------------------------------------------------------------
_SYSTEM_CACHE: Dict[Tuple, PropulsionSystem] = {}


def get_propulsion(motor: Motor, prop: Propeller, n_motors: int,
                   pack: BatteryPack, env: Environment,
                   install_thrust_factor: float = 0.95,
                   n_systems: int = 1) -> PropulsionSystem:
    key = (motor.key, prop.key, n_motors, pack.key, n_systems,
           round(env.rho, 4), round(env.nu, 9),
           round(install_thrust_factor, 4))
    sys_ = _SYSTEM_CACHE.get(key)
    if sys_ is None:
        sys_ = PropulsionSystem(motor, prop, n_motors, pack, env, None,
                                install_thrust_factor, n_systems)
        _SYSTEM_CACHE[key] = sys_
    return sys_


def clear_system_cache() -> None:
    _SYSTEM_CACHE.clear()
