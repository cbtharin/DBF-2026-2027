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
resistance motor makes infinite torque at any speed.  Rm is taken at the
winding's temperature (rm_scale: Aircraft.prop_sys_at re-solves the laps
in 2 C bands as the windings heat), and the pack's internal resistance, the
harness and each wing motor's leads are series resistance in the same
solve - battery sag, solved exactly (PropulsionSystem.solve).

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

Propulsion system architecture - RULE 3.2.3.1a, 3.2.3.1b, 3.2.3.1d, 3.2.3.1e
----------------------------------------------
A propulsion system is one battery pack, one fuse, one externally accessible
arming plug, one or more ESCs and one or more motors.  Packs may not be
wired in series or parallel, so the only legal way to carry more than one
pack is to build a second complete, independent system - and then every pack
must be identical.  `n_systems` is therefore a real design variable, and
`motors_per_system = n_motors / n_systems` must come out a whole number.
A pack with more than one output connector is still one pack: its leads are
combined before its one fuse, and it may not feed a second system (3.2.3.1e).
So every pack's whole current passes one fuse - there is no way to split a
pack across two systems for two fuses' worth of current.

Sag and heating
---------------
Both are modelled (they were once left out, which made full-throttle thrust
optimistic by roughly 3-6 percent on a warm 6S pack): the pack's terminal
voltage falls with its current through its internal resistance and the
harness, and Rm rises 0.393 %/C with the winding temperature the mission
tracks.  The 100 A blade fuse is still checked directly against the
computed pack current.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace
import bisect
from typing import Dict, Optional, Tuple

from .components import ESC, BatteryPack, Motor, Propeller, pick_esc
from .config import BuildStandard, Environment
from .numerics import brent, scan_then_refine


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
    p_motor_loss_w: float = 0.0  # all motors: copper + iron, heats windings


class PropulsionSystem:
    """One or more identical propulsion systems, solved as a set.

    Every system carries the same pack and the same number of identical
    motors turning identical propellers, so they all sit at the same
    operating point and one solve describes the whole aeroplane.
    """

    ROOT_XTOL = 1e-12          # Brent tolerance on rev/s

    def __init__(self, motor: Motor, prop: Propeller, n_motors: int,
                 pack: BatteryPack, env: Environment, esc: Optional[ESC] = None,
                 install_thrust_factor: float = 1.0, n_systems: int = 1,
                 r_harness_ohm: float = 0.003, rm_scale: float = 1.0,
                 max_fuse_a: float = 100.0, esc_airflow_mps: float = 10.0,
                 mission_s: float = 300.0, r_lead_ohm: float = 0.0,
                 i0_speed_share: float = BuildStandard.i0_speed_share):
        self.install_thrust_factor = install_thrust_factor
        # A number, or a table ((v, k), ...) of installed / free-air thrust
        # against airspeed (Aircraft builds it from the slipstream drag at
        # each speed), interpolated in installed_factor().
        self._k_v = ([v_ for v_, _k in install_thrust_factor]
                     if isinstance(install_thrust_factor, tuple) else None)
        # the ESC's working conditions: how long it runs hard and how much
        # air reaches it (ESC.sustained_a)
        self.esc_airflow_mps = esc_airflow_mps
        self.mission_s = mission_s
        # RULE 3.2.3.2d - the largest fuse the rules allow (Rules.
        # max_fuse_current_a)
        self.max_fuse_a = max_fuse_a
        self.motor = motor
        self.prop = prop
        self.n_motors = n_motors
        self.n_systems = max(1, int(n_systems))
        self.motors_per_system = n_motors / float(self.n_systems)
        self.pack = pack
        self.env = env
        self.kt = 60.0 / (2.0 * math.pi * motor.kv)     # N*m per amp
        # Winding resistance at the winding's temperature: the catalogue
        # value is a cold one, and copper gains 0.393 %/C (motor_thermal).
        self.rm = motor.rm * max(rm_scale, 0.1)
        # Battery sag: every system's current flows through its pack's
        # internal resistance and its fuse/plug/lead harness.
        self.r_sys = pack.r_internal_ohm + max(r_harness_ohm, 0.0)
        # ...and each wing motor's own DC leads out to its nacelle: the mean
        # per motor (a nose motor's are in the harness), carrying only that
        # motor's current (layout.wing_power_leads).
        self.r_lead = max(r_lead_ohm, 0.0)
        # The no-load current rises with speed (Motor.i0_at).
        self.i0_speed_share = i0_speed_share
        self._memo_key: Optional[Tuple[float, float, float]] = None
        self._memo_val: Tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)
        self._rating: Optional[Dict[str, float]] = None
        # The ESC is sized on the current the motor REALLY draws - its peak
        # over every airspeed, not the motor's datasheet maximum - in the
        # densest and hottest air it will fly in (the named Tucson
        # conditions, or this air if it is harsher), so every weather
        # preset carries the same, hardware-fixed ESC.  The ESC's own
        # efficiency only scales the battery power, not the motor current,
        # so a provisional one does for the sizing solve.
        if esc is not None:
            self.esc = esc
        else:
            self.esc = _PROVISIONAL_ESC
            rho_s, t_s = _esc_sizing_air(env)
            self.env = replace(env, rho=rho_s)
            i_pk = self.solve(self.peak_current_speed(), 1.0,
                              1.0).motor_current_a
            self.env = env
            self.esc = (pick_esc(i_pk, pack.series, mission_s,
                                 esc_airflow_mps, t_s)
                        or _strongest_esc(pack.series, mission_s,
                                          esc_airflow_mps, t_s))

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
    def fuse_rating_a(self) -> float:        # RULE 3.2.3.2d, 3.2.3.1e
        """Each propulsion system's fuse: the pack's listed C rating times
        its capacity, up to the largest legal fuse - 100 A
        (Rules.max_fuse_current_a).  The rule also accepts the label's own
        maximum discharge current; the team sizes by C x Ah, and every
        catalogue pack that states both gives the same figure."""
        return min(self.pack.capacity_ah * self.pack.c_rating,
                   self.max_fuse_a)

    # ------------------------------------------------------------------
    # Exact single-point solve
    # ------------------------------------------------------------------
    def solve(self, v: float, soc: float, throttle: float) -> OperatingPoint:
        m, p, env = self.motor, self.prop, self.env
        throttle = max(0.0, min(1.0, throttle))
        v = max(0.0, v)
        v_ocv = self.pack.ocv(soc)

        if throttle * v_ocv <= 0.05:
            return OperatingPoint(0.0, 0.0, 0.0, 0.0, 0.0, v_ocv,
                                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        # Battery sag, solved exactly.  The pack's terminal voltage is
        #   V = Voc - R I_pack,   I_pack = (motors per system) thr I_motor
        # (the ESC trades voltage for current), less the drop in the
        # motor's own leads, R_l thr I, and the motor sees thr of that:
        #   thr (Voc - R m thr I - R_l thr I) = E + Rm I
        #   =>  I = (thr Voc - E) / (Rm + m thr^2 R + thr^2 R_l)
        # so sag is an extra series resistance on the motor, and the whole
        # system is still one scalar root-find in speed - the no-load
        # current, which depends on that speed, included.
        voc_t = throttle * v_ocv
        r_eff = (self.rm + self.motors_per_system * throttle * throttle
                 * self.r_sys + throttle * throttle * self.r_lead)
        share = self.i0_speed_share
        d5 = p.diameter_m ** 5
        two_pi = 2.0 * math.pi
        k = 60.0 / m.kv

        def residual(nn: float) -> float:
            cur = (voc_t - k * nn) / r_eff
            q_motor = self.kt * (cur - m.i0_at(nn, share))
            j = v / (nn * p.diameter_m) if nn > 1e-6 else 1e9
            # the blade read at this trial rpm (Propeller.rpm_tables)
            q_prop = p.cp(j, 60.0 * nn) * env.rho * nn * nn * d5 / two_pi
            return q_motor - q_prop

        # At the no-load speed the motor current is zero - and so is the sag
        # - so the residual is negative there: a guaranteed bracket against
        # a positive residual at standstill.
        n_noload = voc_t * m.kv / 60.0
        if residual(1e-4) <= 0.0:
            n = 0.0
        else:
            n = brent(residual, 1e-4, n_noload, xtol=self.ROOT_XTOL)
        i_motor = max(0.0, (voc_t - k * n) / r_eff)

        # Battery current is the motor current scaled by throttle (the ESC
        # is a switching converter: it trades voltage for current).
        i_pack = self.motors_per_system * i_motor * throttle
        i_total = self.n_motors * i_motor * throttle
        v_pack = v_ocv - self.r_sys * i_pack          # loaded terminal volts
        # at the motor: after its own leads (their loss is not motor heat)
        vt = throttle * (v_pack - self.r_lead * throttle * i_motor)
        thrust = (self.n_motors * max(0.0, p.thrust(env.rho, n, v))
                  * self.installed_factor(v))
        p_shaft = self.n_motors * p.power(env.rho, n, v)
        # Energy leaves the pack's chemistry at the OPEN-CIRCUIT voltage: the
        # I^2 R lost inside the pack and harness is energy spent too.  That
        # is how sag costs endurance, not just thrust.
        p_elec = v_ocv * i_total / max(self.esc.efficiency, 0.5)
        j = v / (n * p.diameter_m) if n > 1e-6 else 0.0
        eta_prop = (thrust * v / p_shaft) if p_shaft > 1e-6 and v > 0.1 else 0.0
        eta_motor = (p_shaft / (self.n_motors * max(vt * i_motor, 1e-6))
                     if i_motor > 1e-6 else 0.0)
        # what the motors do not turn into shaft power heats them
        p_loss = max(self.n_motors * vt * i_motor - p_shaft, 0.0)
        return OperatingPoint(60.0 * n, thrust, i_motor, i_pack, i_total,
                              v_pack, p_elec, p_shaft, min(eta_prop, 0.98),
                              min(eta_motor, 0.98), j, p_loss)

    # ------------------------------------------------------------------
    # Query interface
    # ------------------------------------------------------------------
    def operating(self, v: float, soc: float = 1.0,
                  throttle: float = 1.0) -> Tuple[float, float, float, float, float]:
        """Returns (thrust_N, p_elec_W, total_current_A, pack_voltage_V,
        motor_loss_W).

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
               op.pack_voltage_v, op.p_motor_loss_w)
        self._memo_key, self._memo_val = key, out
        return out

    def thrust(self, v: float, soc: float = 1.0, throttle: float = 1.0) -> float:
        return self.operating(v, soc, throttle)[0]

    def power_elec(self, v: float, soc: float = 1.0,
                   throttle: float = 1.0) -> float:
        return self.operating(v, soc, throttle)[1]

    def pack_heat(self, v: float, soc: float = 1.0,
                  throttle: float = 1.0) -> float:
        """Joule heat in ONE pack, W: its own share of the current through
        its own DC resistance (the harness heats the harness, not the
        pack)."""
        i_pack = self.operating(v, soc, throttle)[2] / self.n_systems
        return i_pack * i_pack * self.pack.r_internal_ohm

    def motor_loss(self, v: float, soc: float = 1.0,
                   throttle: float = 1.0) -> float:
        """Heat into the windings, all motors, W."""
        return self.operating(v, soc, throttle)[4]

    def pack_current(self, v: float, soc: float = 1.0,
                     throttle: float = 1.0) -> float:
        """Total current drawn from all packs.  (No rule limits it: each
        pack's own fuse is the limit - pack_current_a in rating_check.)"""
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

    def installed_factor(self, v: float) -> float:
        """Installed thrust / free-air thrust at airspeed v: the number
        given, or the table interpolated linearly (held at its ends)."""
        k = self.install_thrust_factor
        if self._k_v is None:
            return k
        if v <= k[0][0]:
            return k[0][1]
        if v >= k[-1][0]:
            return k[-1][1]
        i = bisect.bisect_right(self._k_v, v) - 1
        (v0, k0), (v1, k1) = k[i], k[i + 1]
        return k0 + (k1 - k0) * (v - v0) / max(v1 - v0, 1e-12)

    def speed_ceiling(self, soc: float = 1.0) -> float:
        """The airspeed above which the propellers make no thrust, m/s:
        their zero-thrust advance ratio (the largest over the rpm tables)
        at the no-load rpm the pack can give at this state of charge.  No
        level flight is possible above it."""
        n_free = self.pack.ocv(soc) * self.motor.kv / 60.0
        return max(self.prop.j_zero_max * n_free * self.prop.diameter_m, 1.0)

    def peak_current_speed(self) -> float:
        """Airspeed at which full charge and full throttle draw the most
        current, m/s.

        Searched from standstill to the propeller's zero-thrust speed at
        the pack's no-load rpm - every speed the aeroplane could be at,
        dives included.  Current rises with the propeller's power
        coefficient, so the peak sits where CP(J) peaks: above J = 0 for 48
        of the 59 measured propellers.
        """
        v_top = self.speed_ceiling(1.0)
        v_pk, i_pk = scan_then_refine(
            lambda v: self.solve(v, 1.0, 1.0).total_current_a,
            0.0, v_top, n_scan=24, tol=1e-3)
        return v_pk if i_pk > self.solve(0.0, 1.0, 1.0).total_current_a             else 0.0

    def rating_check(self) -> Dict[str, float]:
        """Worst-case electrical loads - full charge, full throttle, at the
        airspeed that draws the most current - and the static thrust.

        The worst case is NOT the static point.  A measured propeller's
        power coefficient usually peaks above J = 0, so the current climbs
        as the aeroplane gathers speed before it falls away: 67 A static
        and 74 A at 22 m/s on the recommended design.  This checked only
        the static point, which was the worst case for the parametric
        propeller it was written for (whose CP falls from J = 0) and is
        not for measured ones.  Every current and power margin, and the
        100 A fuse and total checks, now use the peak.  Cached: it depends
        only on the hardware and the air.
        """
        if self._rating is not None:
            return dict(self._rating)
        static = self.solve(0.0, 1.0, 1.0)
        v_pk = self.peak_current_speed()
        op = self.solve(v_pk, 1.0, 1.0) if v_pk > 0.0 else static
        per_motor_w = op.p_elec_w / max(self.n_motors, 1)
        # The fastest the propeller can be driven: Kv x the full-charge
        # open-circuit voltage, the motor's speed with no load and no
        # losses.  The propeller approaches it as it unloads (a
        # full-throttle dive), and for the maker's structural rpm limit,
        # where there is one, the bound must not read low.
        rpm_free = self.motor.kv * self.pack.ocv(1.0)
        self._rating = {
            "static_thrust_n": static.thrust_total_n,
            "static_current_a": static.total_current_a,
            "peak_current_speed_mps": v_pk,
            "motor_current_a": op.motor_current_a,
            "pack_current_a": op.pack_current_a,
            "total_current_a": op.total_current_a,
            "fuse_rating_a": self.fuse_rating_a,
            "motor_power_w": per_motor_w,
            "motor_current_margin": self.motor.i_max_a / max(op.motor_current_a, 1e-6),
            "motor_power_margin": self.motor.p_max_w / max(per_motor_w, 1e-6),
            # the ESC holds its SUSTAINED rating (for the mission, at its
            # cooling airflow, in this air), not its published one
            "esc_sustained_a": self.esc.sustained_a(
                self.mission_s, self.esc_airflow_mps, self.env.temp_c),
            "esc_current_margin": self.esc.sustained_a(
                self.mission_s, self.esc_airflow_mps, self.env.temp_c)
                / max(op.motor_current_a, 1e-6),
            "pack_current_margin": self.pack.i_max_a / max(op.pack_current_a, 1e-6),
            "fuse_margin": self.fuse_rating_a / max(op.pack_current_a, 1e-6),
            "rpm": static.rpm,
            "rpm_no_load": rpm_free,
            "prop_rpm_margin": (self.prop.max_rpm / max(rpm_free, 1.0)
                                if self.prop.max_rpm > 0.0 else 99.0),
            "tip_mach": (math.pi * self.prop.diameter_m * static.rpm / 60.0)
                        / 340.0,
        }
        return dict(self._rating)


def _biggest_esc() -> ESC:
    from .components import ESCS
    return max(ESCS.values(), key=lambda e: e.i_cont_a)


def _strongest_esc(cells: int, mission_s: float, airflow_mps: float,
                   ambient_c: float) -> ESC:
    """No ESC holds the current: the one that holds the most at this cell
    count (the rating check then fails the design, with the numbers)."""
    from .components import ESCS
    fits = [e for e in ESCS.values() if e.min_cells <= cells <= e.max_cells]
    return max(fits or list(ESCS.values()),
               key=lambda e: e.sustained_a(mission_s, airflow_mps, ambient_c))


# the sizing solve only needs an efficiency; the real ESC replaces it
_PROVISIONAL_ESC = ESC("provisional", 1e9, 0.0, 0.96, 99)


def _esc_sizing_air(env: Environment) -> Tuple[float, float]:
    """(density, temperature) to size the ESC in: the densest and hottest
    of the named Tucson conditions, or this air if it is harsher."""
    from .config import CONDITIONS
    rho, temp = env.rho, env.temp_c
    for name, (alt, t_c, wind) in CONDITIONS.items():
        if name.startswith("tucson"):
            e = Environment.at_density_altitude(alt, t_c, wind)
            rho, temp = max(rho, e.rho), max(temp, t_c)
    return rho, temp


# --------------------------------------------------------------------------
# Instance cache.  A PropulsionSystem depends only on the hardware, never on
# the wing, so thousands of candidate wings share one object - and with it
# its last-operating-point memo.  Bounded by the hardware shortlist size.
# --------------------------------------------------------------------------
_SYSTEM_CACHE: Dict[Tuple, PropulsionSystem] = {}


def get_propulsion(motor: Motor, prop: Propeller, n_motors: int,
                   pack: BatteryPack, env: Environment,
                   install_thrust_factor: float = 1.0,
                   n_systems: int = 1,
                   r_harness_ohm: float = 0.003,
                   rm_scale: float = 1.0,
                   max_fuse_a: float = 100.0,
                   esc: Optional[ESC] = None,
                   esc_airflow_mps: float = 10.0,
                   mission_s: float = 300.0,
                   r_lead_ohm: float = 0.0,
                   i0_speed_share: float = BuildStandard.i0_speed_share
                   ) -> PropulsionSystem:
    # The whole motor, propeller and pack records are the key, not their
    # names: an entry edited under the same key in the Catalogues tab (a
    # measured Rm or I0, a new mass) must not be answered by a system
    # cached before the edit.  Keying on the names alone did exactly that
    # for motors and propellers.
    key = (motor, prop, n_motors, pack, n_systems,
           round(env.rho, 4), round(env.nu, 9),
           (tuple((round(a, 3), round(k, 5)) for a, k in install_thrust_factor)
            if isinstance(install_thrust_factor, tuple)
            else round(install_thrust_factor, 4)), round(r_harness_ohm, 7),
           round(rm_scale, 4), round(max_fuse_a, 6),
           esc, round(esc_airflow_mps, 4), round(mission_s, 3),
           round(r_lead_ohm, 7), round(i0_speed_share, 4))
    sys_ = _SYSTEM_CACHE.get(key)
    if sys_ is None:
        sys_ = PropulsionSystem(motor, prop, n_motors, pack, env, esc,
                                install_thrust_factor, n_systems,
                                r_harness_ohm, rm_scale, max_fuse_a,
                                esc_airflow_mps, mission_s, r_lead_ohm,
                                i0_speed_share)
        _SYSTEM_CACHE[key] = sys_
    return sys_


def clear_system_cache() -> None:
    _SYSTEM_CACHE.clear()


# ==========================================================================
# Motor heating
# ==========================================================================
def motor_thermal(motor: Motor, build) -> Tuple[float, float]:
    """(thermal resistance K/W, heat capacity J/K) of one motor.

    A lumped model: the winding and can heat with the motor's own losses
    and shed heat to the prop wash,

        C dT/dt = P_loss - (T - T_air) / R

    C is the motor's mass times a copper/steel/magnet specific heat.  R
    scales with surface area, k / m^(2/3), with k calibrated on the
    catalogue motors whose makers rate them continuously - at that rated
    current they sit at the temperature limit.  That puts every motor on
    one footing whatever the maker's rating basis (15 s, 30 s, 180 s or
    continuous); the short ratings remain as hard current and power limits.
    """
    m = max(motor.mass_kg, 0.02)
    return (build.motor_thermal_k / m ** (2.0 / 3.0),
            m * build.motor_specific_heat_j_kgk)


def winding_step(t_w: float, t_air: float, p_loss_w: float, dt: float,
                 r_th: float, c_th: float) -> float:
    """Exact step of the lumped model over dt at constant loss power."""
    if dt <= 0.0:
        return t_w
    x = dt / max(r_th * c_th, 1e-6)
    # 1 - exp(-x) as -expm1(-x): a short step keeps its full precision
    return (t_air + (t_w - t_air) * math.exp(-x)
            - p_loss_w * r_th * math.expm1(-x))
