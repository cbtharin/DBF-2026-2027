"""The Design vector and the Aircraft it produces.

`Design` is exactly the list of things you can change:
    airfoil, span, chord, motor type, propeller, motor count, battery pack
    (chosen whole off the shelf), the number of independent propulsion
    systems, sensor mass, and the cargo you choose to fly in Mission 2.

Battery packs are COTS items picked from a catalogue, not built from cells:
the rules forbid altered packs and forbid wiring packs in series or
parallel, so "cells in series" and "cells in parallel" are not design
variables any more - the pack you buy is the pack you fly.

`Aircraft` closes the loop: the payload count sets the bay (from one fixed
container size), the bay sets the fuselage, the fuselage sets the tail arm,
the tails and wing set the drag, a V-n diagram sets the ultimate load
factor, that sizes the spar, the spar sets the weight, and the weight feeds
back into both the load factor and the structure until it converges.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

from .aerodynamics import AeroModel, WingGeometry, build_body
from .components import AIRFOILS, MOTORS, PACKS, PROPS
from .config import Config
from .loads import LoadCase, design_load_factors
from .numerics import brent
from .payload import PayloadGeometry, best_bay_arrangement, size_payload
from .propulsion import PropulsionSystem, get_propulsion
from .units import G0, KG2LB
from .weights import MassBreakdown, converge_empty_mass, spar_sizing


@dataclass(frozen=True)
class Design:
    # --- geometry -----------------------------------------------------
    span_m: float
    chord_m: float
    airfoil: str
    # --- propulsion ---------------------------------------------------
    motor: str
    prop: str
    n_motors: int
    # --- energy -------------------------------------------------------
    # A catalogue key: one unaltered, commercially procured pack.
    pack: str
    # --- payload ------------------------------------------------------
    # The maximum sensor weight declared at Tech Inspection.  Mission 2
    # must fly this weight, every simulator must match it, and the Ground
    # Mission is scored on it.
    sensor_mass_kg: float
    n_simulators: int
    taper: float = 1.0
    # Mission 3 may fly a LIGHTER sensor than the declared maximum (rule
    # 3.1.1).  0 means "fly the declared maximum".
    sensor_m3_kg: float = 0.0
    # The declared Ground Mission drop height, in whole inches, up to 60.
    drop_height_in: float = 60.0
    # Independent propulsion systems, each with its own pack, fuse and
    # arming plug.  Packs may not be paralleled, so building a second
    # system is the only way to carry a second pack.  Must divide
    # n_motors evenly.  Last in the list only because it has a default.
    n_systems: int = 1

    def key(self) -> Tuple:
        return (round(self.span_m, 4), round(self.chord_m, 4), self.airfoil,
                self.motor, self.prop, self.n_motors, self.pack,
                self.n_systems, round(self.sensor_mass_kg, 4),
                self.n_simulators, round(self.sensor_m3_kg, 4),
                round(self.drop_height_in, 2), round(self.taper, 3))

    def discrete_key(self) -> Tuple:
        return (self.airfoil, self.motor, self.prop, self.n_motors,
                self.pack, self.n_systems)


class Aircraft:
    def __init__(self, design: Design, cfg: Config):
        self.design = design
        self.cfg = cfg
        r, b, env = cfg.rules, cfg.build, cfg.env

        self.airfoil = AIRFOILS[design.airfoil]
        self.motor = MOTORS[design.motor]
        self.prop = PROPS[design.prop]
        # RULE 3.2.3a, 3.2.3.1c - an atomic COTS catalogue entry.  There
        # is no way to express cells in series or parallel here, because
        # that would be a custom pack; and one key for the whole
        # aeroplane means multiple packs cannot differ.
        self.pack = PACKS[design.pack]
        # One pack per propulsion system, every pack identical.
        self.n_systems = max(1, int(design.n_systems))
        self.total_watt_hours = self.n_systems * self.pack.watt_hours

        # ---- payload: the container is sized by the sensor -----------
        self.payload_geom: PayloadGeometry = size_payload(
            design.sensor_mass_kg, r, design.drop_height_in)
        container_tare = self.payload_geom.container_tare_kg
        # RULE 3.1.3a: "The shipping container simulator must be the
        # same weight as the sensor shipping container plus the maximum
        # weight sensor declared in Tech Inspection within +/- 1 ounce."
        #
        # So simulator mass is NOT a design variable.  Every simulator
        # weighs exactly what the loaded container weighs, and Mission 2's
        # scored weight is simply (1 + n) loaded containers.  Treating it
        # as free lets the optimiser fly a stack of light boxes, which is
        # not a legal aeroplane.
        self.sensor_max_kg = design.sensor_mass_kg
        self.loaded_container_kg = self.sensor_max_kg + container_tare
        self.sim_mass_kg = self.loaded_container_kg
        self.payload_sim_total_kg = design.n_simulators * self.sim_mass_kg
        self.m2_scored_kg = self.loaded_container_kg + self.payload_sim_total_kg
        self.m2_payload_kg = self.m2_scored_kg
        # RULE 3.1.1c: the sensor weight may be reduced for a Mission 3
        # attempt, but never raised above the declared maximum.  A lighter
        # sensor flies more laps; M3 scores laps x weight, so there is a
        # genuine optimum in between and it is a design variable.
        # 0 means "fly the declared maximum".  An explicit value is taken
        # as given and NOT clamped: a request to fly more than the declared
        # maximum is illegal (rule 3.1.1) and static_checks() rejects it,
        # rather than being silently corrected into a legal design that is
        # not the one that was asked for.
        self.sensor_m3_kg = (design.sensor_m3_kg if design.sensor_m3_kg > 0
                             else self.sensor_max_kg)
        self.m3_payload_kg = self.sensor_m3_kg
        self.sensor_density_kg_m3 = self.payload_geom.required_density(
            design.sensor_mass_kg, r.container_fill_fraction)

        # ---- geometry ------------------------------------------------
        self.wing = WingGeometry(design.span_m, design.chord_m,
                                 self.airfoil, design.taper)
        pg = self.payload_geom
        # RULE 3.1.1d - the sensor is carried internally.  The bay is
        # always at least one container envelope, and rule 3.1.3c makes
        # the sensor fit inside a container, so the Mission 3 stowed
        # sensor fits by construction: it is strictly smaller than the
        # loaded box Mission 2 already has to carry.
        n_boxes = 1 + design.n_simulators
        rows, cols, layers, bay_l, bay_w, bay_h = best_bay_arrangement(
            n_boxes, pg.box_l, pg.box_w, pg.box_h)
        self.bay_rows, self.bay_cols, self.bay_layers = rows, cols, layers
        self.body = build_body(bay_l, bay_w, bay_h, self.wing, b)

        # ---- landing gear height is set by propeller clearance -------
        if design.n_motors == 1:
            need = (0.5 * self.prop.diameter_m + b.prop_clearance_margin_m
                    - 0.5 * self.body.fuse_height_m)
            self.gear_height_m = max(b.gear_height_m, need)
        else:
            self.gear_height_m = b.gear_height_m
        gear_scale = self.gear_height_m / b.gear_height_m
        gear_extra_drag = b.gear_drag_area_m2 * (gear_scale - 1.0)
        self.gear_extra_mass_kg = 0.55 * (self.gear_height_m - b.gear_height_m)

        # ---- aero ----------------------------------------------------
        self.aero_clean = AeroModel(self.wing, self.body, env, b,
                                    design.n_motors,
                                    extra_drag_area_m2=gear_extra_drag)
        # RULE 3.1.2, 3.1.1a
        # Mission 3 tows the sensor on a line at least 1.5 x wingspan
        # long, so the deployed drag is the sensor body PLUS
        # that line.  A 2.7 m length of 1.5 mm line at Cd 1.1 is about
        # 0.0045 m2 of CD*A - the same order as the sensor itself, and it
        # scales with span, so a bigger wing is penalised twice in M3.
        self.tow_line_length_m = r.tow_line_span_multiple * self.wing.span_m
        self.tow_line_drag_area_m2 = (r.tow_line_cd
                                      * r.tow_line_diameter_m
                                      * self.tow_line_length_m)
        self.deployed_drag_area_m2 = (r.sensor_deployed_drag_area_m2
                                      + self.tow_line_drag_area_m2)
        self.aero_deployed = AeroModel(self.wing, self.body, env, b,
                                       design.n_motors,
                                       extra_drag_area_m2=gear_extra_drag
                                       + self.deployed_drag_area_m2)

        # ---- propulsion ----------------------------------------------
        self.prop_sys: PropulsionSystem = get_propulsion(
            self.motor, self.prop, design.n_motors, self.pack, env,
            b.install_thrust_factor, self.n_systems)

        # ---- mass ----------------------------------------------------
        nm = design.n_motors
        esc = self.prop_sys.esc
        propulsion_items = {
            "motors (%d)" % nm: nm * self.motor.mass_kg,
            "propellers (%d)" % nm: nm * self.prop.mass_kg,
            "ESCs (%d)" % nm: nm * esc.mass_kg,
            "motor mounts/spinners": nm * 0.025,
            "propulsion packs (%d)" % self.n_systems:
                self.n_systems * self.pack.mass_kg,
            # 3.2.3.1: one blade-fuse harness and one externally accessible
            # arming plug per propulsion system.
            "fuse harness + arming plug (%d)" % self.n_systems:
                self.n_systems * (b.fuse_harness_mass_kg
                                  + b.arming_plug_mass_kg),
            # 3.2.3: a separate Rx/servo battery is mandatory.  It carries
            # no propulsion energy but it does carry mass.
            "Rx/servo battery": b.rx_battery_mass_kg,
        }
        if self.gear_extra_mass_kg > 1e-9:
            propulsion_items["tall gear penalty"] = self.gear_extra_mass_kg
        self.propulsion_items = propulsion_items

        # The spar is sized by the ultimate load factor, which comes from a
        # V-n diagram, which needs the weight and the speed - which need the
        # spar:
        #   spar mass  <- ultimate load factor <- weight and Vmax
        #   weight     <- spar mass
        # Resolve by alternating: converge the mass, rebuild the load case,
        # converge again.  A few passes settle it to well under a gram.
        n_override = b.ultimate_load_factor_override
        n_ult = n_override if n_override else 4.0
        mb = None
        for _ in range(6):
            mb, empty, mtow_m2 = converge_empty_mass(
                self.wing, self.body, propulsion_items, self.m2_payload_kg,
                b, n_ult)
            w = mtow_m2 * G0
            self.load_case: LoadCase = design_load_factors(
                self.aero_clean, w, self._max_level_speed(w), b)
            n_new = (n_override if n_override
                     else self.load_case.n_ultimate)
            if abs(n_new - n_ult) < 1e-6:
                n_ult = n_new
                break
            n_ult = n_new
        mb, empty, mtow_m2 = converge_empty_mass(
            self.wing, self.body, propulsion_items, self.m2_payload_kg,
            b, n_ult)
        self.n_ultimate = n_ult
        self.n_limit = (n_ult / b.load_safety_factor if n_override
                        else self.load_case.n_limit)

        (_sm, self.spar_t_req_m, self.spar_t_used_m,
         self.spar_gauge_limited) = spar_sizing(
            self.wing, mtow_m2 * G0, b, n_ult)

        self.mass_breakdown: MassBreakdown = mb
        self.empty_mass_kg = empty
        self.mtow_m2_kg = mtow_m2
        self.mtow_m3_kg = empty + self.m3_payload_kg
        self.mtow_m1_kg = empty            # Mission 1 flies empty

        self.w_m2_n = self.mtow_m2_kg * G0
        self.w_m3_n = self.mtow_m3_kg * G0
        self.w_m1_n = self.mtow_m1_kg * G0

    def _max_level_speed(self, weight_n: float) -> float:
        """Maximum level speed, solved locally.

        mission.py imports this module, so it cannot be imported back; this
        is the same Brent solve on thrust minus drag.
        """
        aero = self.aero_clean

        def excess(v: float) -> float:
            return self.prop_sys.thrust(v, 0.85, 1.0) - aero.drag(v, weight_n)

        lo = aero.v_stall(weight_n) * 1.02
        hi = 65.0
        if excess(lo) <= 0.0:
            return lo
        if excess(hi) > 0.0:
            return hi
        return brent(excess, lo, hi, xtol=1e-6)

    # ------------------------------------------------------------------
    @property
    def wing_mass_kg(self) -> float:
        return (self.mass_breakdown.items.get("wing skin/ribs", 0.0)
                + self.mass_breakdown.items.get("wing spar", 0.0))

    @property
    def motor_mass_kg(self) -> float:
        return self.design.n_motors * self.motor.mass_kg

    @property
    def wing_area_m2(self) -> float:
        return self.wing.area_m2

    @property
    def aspect_ratio(self) -> float:
        return self.wing.aspect_ratio

    def wing_loading_m2(self) -> float:
        return self.w_m2_n / self.wing_area_m2

    def static_thrust_n(self) -> float:
        return self.prop_sys.rating_check()["static_thrust_n"]

    def usable_energy_j(self) -> float:
        """Energy you can actually get out of every pack, chemistry
        dependent.  Li-ion gives up more of its rated capacity than LiPo
        because it is not being asked for anything like its peak current."""
        return (self.total_watt_hours * 3600.0 * self.pack.usable_fraction)

    # ------------------------------------------------------------------
    # Hard constraints that do not need a mission simulation
    # ------------------------------------------------------------------
    def static_checks(self) -> Dict[str, Any]:
        r, b = self.cfg.rules, self.cfg.build
        rc = self.prop_sys.rating_check()
        d = self.design

        clearance = (self.gear_height_m + 0.5 * self.body.fuse_height_m
                     - 0.5 * self.prop.diameter_m) if d.n_motors == 1 else 1.0
        if d.n_motors >= 2:
            lateral = self.wing.span_m / d.n_motors
            prop_gap = lateral - self.prop.diameter_m
        else:
            prop_gap = 1.0

        fails: List[str] = []
        # --- rules limits, checked first: these are not trade-offs -------
        # RULE 3.2.3.2b
        if self.total_watt_hours > r.max_watt_hours + 1e-9:
            fails.append("propulsion energy %.1f Wh over the %.0f Wh limit"
                         % (self.total_watt_hours, r.max_watt_hours))
        # RULE 3.2.3.2c
        if self.pack.watt_hours > r.max_pack_watt_hours + 1e-9:
            fails.append("single pack %.1f Wh over the %.0f Wh FAA limit"
                         % (self.pack.watt_hours, r.max_pack_watt_hours))
        # RULE 3.2.3.2a
        if self.pack.chemistry not in r.battery_chemistry_allowed:
            fails.append("%s packs are not in the searched set"
                         % self.pack.chemistry)
        if d.n_motors > r.max_motors:
            fails.append("%d motors, limit is %d" % (d.n_motors, r.max_motors))
        # RULE 3.2.3.1a, 3.2.3.1b - one pack per system, and packs may
        # not be paralleled, so a second pack means a second full system.
        if self.n_systems > d.n_motors:
            fails.append("more propulsion systems than motors")
        if abs(d.n_motors / float(self.n_systems)
               - round(d.n_motors / float(self.n_systems))) > 1e-9:
            fails.append("%d motors do not divide evenly across %d "
                         "propulsion systems" % (d.n_motors, self.n_systems))
        # RULE 3.2.1a - checked at Mission 2 weight, the heaviest case.
        if self.mtow_m2_kg > r.max_gross_mass_kg + 1e-9:
            fails.append("gross weight %.1f lb over the %.0f lb limit"
                         % (self.mtow_m2_kg * KG2LB,
                            r.max_gross_mass_kg * KG2LB))
        # RULE 3.1
        # The tolerance is not cosmetic.  6 ft is 1.8288000000000002 m in
        # binary, and a limit that has been through any text field comes
        # back as 1.8288 - so a design sitting exactly ON the limit, which
        # is where the optimiser puts it, was being rejected by 2e-16 m.
        # Every other limit here already carried this tolerance.
        if self.wing.span_m > r.max_wingspan_m + 1e-9:
            fails.append("span %.4f m over the %.4f m limit"
                         % (self.wing.span_m, r.max_wingspan_m))
        # --- current --- RULE 3.2.3.2d ---
        # The 100 A blade fuse is the hard electrical stop.
        if rc["total_current_a"] > r.max_total_current_a + 1e-9:
            fails.append("draws %.0f A total, limit is %.0f A"
                         % (rc["total_current_a"], r.max_total_current_a))
        if rc["pack_current_a"] > self.prop_sys.fuse_rating_a + 1e-9:
            fails.append("blows its %.0f A fuse (%.0f A per pack)"
                         % (self.prop_sys.fuse_rating_a, rc["pack_current_a"]))
        # --- hardware ratings --------------------------------------------
        if rc["motor_current_margin"] < 1.0:
            fails.append("motor over current (%.0f A draw)" % rc["motor_current_a"])
        if rc["motor_power_margin"] < 1.0:
            fails.append("motor over power")
        if rc["esc_current_margin"] < 1.0:
            fails.append("ESC over current")
        if rc["pack_current_margin"] < 1.0:
            fails.append("pack over C-rating (%.0f A)" % rc["pack_current_a"])
        if rc["tip_mach"] > 0.82:
            fails.append("prop tip Mach %.2f" % rc["tip_mach"])
        if self.gear_height_m > 0.30:
            fails.append("landing gear would need %.0f mm legs for prop clearance"
                         % (self.gear_height_m * 1000))
        if prop_gap < 0.03:
            fails.append("propellers overlap")
        # the tail has to be able to trim without stalling itself
        v_trim = self.aero_clean.v_stall(self.w_m2_n) * 1.3
        ct = abs(self.aero_clean.tail_cl(v_trim, self.w_m2_n))
        ct_turn = abs(self.aero_clean.tail_cl(
            v_trim, self.n_limit * self.w_m2_n))
        if max(ct, ct_turn) > b.tail_cl_max:
            fails.append("tail needs CL %.2f to trim (max %.2f) - bigger "
                         "tail or longer arm" % (max(ct, ct_turn), b.tail_cl_max))
        sm = self.aero_clean.static_margin
        if sm < b.min_static_margin:
            fails.append("static margin %.3f c - %s" % (
                sm, "unstable" if sm < 0 else "too little for a model"))
        if sm > b.max_static_margin:
            fails.append("static margin %.3f c - excessive, all trim drag" % sm)
        if self.aspect_ratio < 3.0 or self.aspect_ratio > 16.0:
            fails.append("aspect ratio %.1f outside modelled range" % self.aspect_ratio)
        if self.wing.chord_m < 0.12:
            fails.append("chord too small (low Re)")
        # the payload has to physically fit in the fuselage
        if self.payload_geom.box_w > self.body.fuse_width_m + 1e-9:
            fails.append("container too wide for the fuselage")
        # There is NO rules cap on sensor mass, so the only ceiling is
        # physics: the sensor has to be made of something.  A heavy sensor
        # in a fixed box is simply a dense one, and the report says what
        # density it would need.  Only an impossible density is rejected.
        # RULE 3.1.1b: the sensor must be at least six inches long.  The
        # envelope is derived from its mass and packing density, so a very
        # light sensor can come out too short to be legal.
        if self.payload_geom.sensor_l < r.sensor_min_length_in * 0.0254:
            fails.append("sensor is %.1f in long, minimum is %.0f in"
                         % (self.payload_geom.sensor_l / 0.0254,
                            r.sensor_min_length_in))
        # RULE 3.1.1c: Mission 3 may fly a lighter sensor, never a
        # heavier one than the weight declared at Tech Inspection.
        if self.sensor_m3_kg > self.sensor_max_kg + 1e-9:
            fails.append("Mission 3 sensor is heavier than the declared "
                         "maximum")
        if self.sensor_density_kg_m3 > r.max_sensor_density_kg_m3:
            fails.append("sensor would need %.0f kg/m3 to fit the container "
                         "(denser than lead)" % self.sensor_density_kg_m3)

        return {
            "ok": not fails,
            "failures": fails,
            "ratings": rc,
            "prop_ground_clearance_m": clearance,
            "prop_lateral_gap_m": prop_gap,
            "structure_ok": self.structure_ok,
        }

    # ------------------------------------------------------------------
    @property
    def structure_ok(self) -> bool:
        """Will the wing hold together?

        The spar is sized to carry the ultimate load with the working
        stress already inside the allowable, so the answer is yes by
        construction - unless the sizing cannot close, which happens in
        exactly one way: the wall it asks for is thicker than the tube's
        own radius.  At that point it is not a tube any more, and no
        amount of layup fixes it; the wing is too long, too thin or too
        heavily loaded.  That is the whole question the team asked for:
        a verdict, not a stress report.
        """
        r_spar = max(0.40 * self.wing.thickness_root_m, 0.006)
        return self.spar_t_used_m < r_spar

    # ------------------------------------------------------------------
    def summary(self) -> Dict[str, Any]:
        rc = self.prop_sys.rating_check()
        return {
            "span_m": self.wing.span_m,
            "chord_m": self.wing.chord_m,
            "area_m2": self.wing_area_m2,
            "AR": self.aspect_ratio,
            "airfoil": self.airfoil.name,
            "CLmax_3d": self.aero_clean.cl_max_at_v(
                self.aero_clean.v_stall(self.w_m2_n)),
            "oswald_e": self.aero_clean.e,
            "CD0_at_20ms": self.aero_clean.cd0(20.0),
            "motor": self.motor.name,
            "n_motors": self.design.n_motors,
            "prop": self.prop.key,
            "pack": self.pack.key,
            "pack_name": self.pack.name,
            "n_systems": self.n_systems,
            "pack_V_nom": self.pack.v_nominal,
            "pack_Wh": self.pack.watt_hours,
            "total_Wh": self.total_watt_hours,
            "chemistry": self.pack.chemistry,
            "pack_c_rating": self.pack.c_rating,
            "fuse_a": self.prop_sys.fuse_rating_a,
            "esc": self.prop_sys.esc.key,
            "wing_mass_kg": self.wing_mass_kg,
            "motor_mass_kg": self.motor_mass_kg,
            "battery_mass_kg": self.n_systems * self.pack.mass_kg,
            "empty_kg": self.empty_mass_kg,
            "mtow_m2_kg": self.mtow_m2_kg,
            "mtow_m3_kg": self.mtow_m3_kg,
            "m2_scored_lb": self.m2_scored_kg * KG2LB,
            "sensor_lb": self.sensor_max_kg * KG2LB,
            "sensor_m3_lb": self.sensor_m3_kg * KG2LB,
            "sim_mass_lb": self.sim_mass_kg * KG2LB,
            "drop_height_in": self.design.drop_height_in,
            "tow_line_m": self.tow_line_length_m,
            "deployed_drag_area_m2": self.deployed_drag_area_m2,
            "container_lb": self.payload_geom.container_tare_kg * KG2LB,
            "sensor_density_kg_m3": self.sensor_density_kg_m3,
            "static_margin": self.aero_clean.static_margin,
            "neutral_point_c": self.aero_clean.x_np_over_c,
            "taper": self.design.taper,
            "tail_cl_cruise": self.aero_clean.tail_cl(
                self.aero_clean.v_stall(self.w_m2_n) * 1.3, self.w_m2_n),
            "static_thrust_n": rc["static_thrust_n"],
            "static_current_a": rc["total_current_a"],
            "static_T_over_W_m2": rc["static_thrust_n"] / self.w_m2_n,
            "wing_loading_Nm2": self.wing_loading_m2(),
            "n_limit": self.n_limit,
            "n_ultimate": self.n_ultimate,
            "load_driver": getattr(self.load_case, "driver", "override"),
            "spar_gauge_limited": self.spar_gauge_limited,
            "structure_ok": self.structure_ok,
        }
