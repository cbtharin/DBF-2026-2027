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

import json
import math
import os
from dataclasses import dataclass, fields, replace
from typing import Any, Dict, List, Optional, Tuple

from .aerodynamics import AeroModel, WingGeometry, build_body
from .components import AIRFOILS, MOTORS, PACKS, PROPS
from .config import Config, Environment
from .loads import LoadCase, design_load_factors
from .numerics import brent
from .payload import PayloadGeometry, best_bay_arrangement, size_payload
from .propulsion import PropulsionSystem, get_propulsion
from .units import G0, KG2LB, MS2MPH
from .layout import Layout, solve_layout
from .liftingline import clear_caches as _clear_lifting_line
from .propulsion import clear_system_cache as _clear_systems
from .weights import (MassBreakdown, SparDesign, converge_empty_mass,
                      spar_design)


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
    # where the wing propellers go: "mid", "tip" or "stagger" (layout.motor_stations)
    wing_prop_layout: str = "mid"
    # height of the vertical winglet at each tip, m (0 = none, m110)
    winglet_h_m: float = 0.0
    # horizontal tail volume chosen, S_h l_t / (S c): never below what the
    # requirements call for (m111); 0 = just the requirements (m112)
    tail_vol_h: float = 0.0

    def key(self) -> Tuple:
        return (round(self.span_m, 4), round(self.chord_m, 4), self.airfoil,
                self.motor, self.prop, self.n_motors, self.pack,
                self.n_systems, round(self.sensor_mass_kg, 4),
                self.n_simulators, round(self.sensor_m3_kg, 4),
                round(self.drop_height_in, 2), round(self.taper, 3),
                round(self.winglet_h_m, 4), round(self.tail_vol_h, 4))

    def discrete_key(self) -> Tuple:
        return (self.airfoil, self.motor, self.prop, self.n_motors,
                self.pack, self.n_systems, self.wing_prop_layout)


#: The design both front ends open on, relative to the project root (the
#: folder above this package) so it resolves from any working directory.
OPENING_DESIGN_FILE = "examples/optimised_v1.json"


def opening_design() -> Tuple[Design, str]:
    """The design the GUI and the web page open on, and where it came from.

    Read from examples/optimised_v1.json when it is there, so both open on
    whatever the last search actually recommended rather than a frozen
    claim about what a good aeroplane looks like.

    Every Design field in the file is passed through.  The GUI used to
    rebuild the Design one named field at a time and left `taper` out, so
    it opened on a rectangular wing instead of the optimised 0.45 taper -
    and scored a different aeroplane, which at the edge of the envelope
    failed Mission 2.  The web page did not read the file at all and
    opened on a hard-coded design from several revisions ago.  Both now
    come through here.
    """
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(here, *OPENING_DESIGN_FILE.split("/"))
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        known = {f.name for f in fields(Design)}
        d = Design(**{k: v for k, v in data.items() if k in known})
        if d.sensor_m3_kg <= 0.0:
            d = replace(d, sensor_m3_kg=d.sensor_mass_kg)
        return d, OPENING_DESIGN_FILE
    except (OSError, ValueError, TypeError):
        return (Design(span_m=1.55, chord_m=1.55 / 5.74, airfoil="SD7062",
                       motor="SCORPION_SII_4025_440", prop="APC_12x8", n_motors=1,
                       pack="GA_6S3300_60C", sensor_mass_kg=0.9,
                       n_simulators=3, taper=1.0, sensor_m3_kg=0.9,
                       drop_height_in=60.0, n_systems=1),
                "built-in fallback")


def gear_height_for(prop_d: float, prop_axis_m: float, build) -> float:
    """Landing-gear height for a propeller of diameter prop_d whose axis
    sits prop_axis_m above the fuselage bottom: the standard leg, or
    taller if the disc would not clear the ground by the build margin."""
    need = 0.5 * prop_d + build.prop_clearance_margin_m - prop_axis_m
    return max(build.gear_height_m, need)


def gear_extra_mass_kg(gear_height_m: float, build) -> float:
    """Mass of the gear's extra height over the standard leg, kg."""
    return 0.55 * (gear_height_m - build.gear_height_m)


def propulsion_mass_items(motor, prop, esc, pack, n_motors: int,
                          n_systems: int, n_wing_motors: int,
                          lead_mass_kg: float, gear_extra_kg: float,
                          n_stagger_ext: int, build) -> Dict[str, float]:
    """Every mass the propulsion brings, itemised as the weight buildup
    lists it, kg.  Aircraft builds its propulsion items from this, and the
    stage-1 screen (optimizer._installed_mass) weighs each candidate set
    with it - one formula, so the two cannot drift apart.  The weight
    buildup adds contingency (BuildStandard.contingency_frac) on top."""
    from .layout import motor_diameter_m
    b, nm = build, n_motors
    items = {
        "motors (%d)" % nm: nm * motor.mass_kg,
        "propellers (%d)" % nm: nm * prop.mass_kg,
        "ESCs (%d)" % nm: nm * esc.mass_kg,
        # scaled with the can's face: 25 g for a 62 mm can
        "motor mounts/spinners": nm * b.motor_mount_kg * (
            motor_diameter_m(motor.mass_kg) / 0.0621) ** 2,
        "propulsion packs (%d)" % n_systems: n_systems * pack.mass_kg,
        # 3.2.3.1: one blade-fuse harness and one externally accessible
        # arming plug per propulsion system.
        "fuse harness + arming plug (%d)" % n_systems:
            n_systems * (b.fuse_harness_mass_kg + b.arming_plug_mass_kg),
        # 3.2.3: a separate Rx/servo battery is mandatory.  It carries
        # no propulsion energy but it does carry mass.
        "Rx/servo battery": b.rx_battery_mass_kg,
    }
    if n_wing_motors:
        # each wing motor's pod, and its power leads
        items["wing nacelles (%d)" % n_wing_motors] = (
            n_wing_motors * b.nacelle_mass_frac
            * (motor.mass_kg + esc.mass_kg))
        items["wing power leads (%d)" % n_wing_motors] = lead_mass_kg
    if gear_extra_kg > 1e-9:
        items["tall gear penalty"] = gear_extra_kg
    if n_stagger_ext:
        # the forward nacelles' booms: 0.6 D each, 100 g/m
        items["stagger extensions (%d)" % n_stagger_ext] = (
            n_stagger_ext * 0.6 * prop.diameter_m * 0.10)
    return items


class Aircraft:
    def __init__(self, design: Design, cfg: Config):
        # Every aircraft starts from empty shared tables (lifting line,
        # propulsion systems).  They are keyed on rounded inputs but hold
        # the first caller's exact answer, so without this a design scored
        # after a near neighbour got the neighbour's numbers, and a search
        # depended on which worker scored what (m99).  Costs nothing
        # measurable: each design refills them in its first flight.
        _clear_lifting_line()
        _clear_systems()
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
        # what the packs really hold (measured capacity when known) - the
        # missions draw on this; the rules check total_watt_hours
        self.total_energy_wh = self.n_systems * self.pack.energy_wh

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
                                 self.airfoil, design.taper,
                                 winglet_h_m=max(design.winglet_h_m, 0.0))
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
        self.bay_length_m, self.bay_width_m, self.bay_height_m = (
            bay_l, bay_w, bay_h)
        # The cross-section wraps the bay whatever the length turns out to
        # be, so the gear can be sized from it before the layout is known.
        body0 = build_body(bay_l, bay_w, bay_h, self.wing, b)

        # ---- landing gear height is set by propeller clearance -------
        # A nose propeller turns on the fuselage centreline; wing
        # propellers on the wing's thrust line - the fuselage top for a
        # top wing, its middle for a mid wing.  Both must clear the ground.
        from .layout import motor_stations
        n_nose, ys, _share = motor_stations(design.n_motors, self.wing.span_m,
                                            body0.fuse_width_m)
        axis = [0.5 * body0.fuse_height_m] if n_nose else []
        if ys:
            axis.append(body0.fuse_height_m
                        if b.wing_covered_surfaces == 1
                        else 0.5 * body0.fuse_height_m)
        # lowest propeller axis above the fuselage bottom (the screen
        # sizes other propellers' gear from it: optimizer._installed_mass)
        self.prop_axis_m = min(axis)
        self.gear_height_m = gear_height_for(self.prop.diameter_m,
                                             self.prop_axis_m, b)
        # the extra height over the standard leg, as a share of it - a
        # difference, rather than (ratio - 1)
        gear_extra_drag = (b.gear_drag_area_m2
                           * (self.gear_height_m - b.gear_height_m)
                           / b.gear_height_m)
        self.gear_extra_mass_kg = gear_extra_mass_kg(self.gear_height_m, b)

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
        # ---- wing-motor power leads -----------------------------------
        # A wing motor's ESC rides in its nacelle, so its DC leads run out
        # from the battery bay: their resistance goes into the power solve
        # and their mass into the weight (layout.wing_power_leads).
        from .layout import wing_power_leads
        self.lead_r_ohm, self.lead_mass_kg, self.lead_lengths_m = (
            wing_power_leads(design.n_motors, self.wing.span_m,
                             body0.fuse_width_m, self.motor.i_max_a, b,
                             design.wing_prop_layout, self.prop.diameter_m,
                             r.wingspan_counts_props))
        self.n_wing_motors = len(self.lead_lengths_m)
        from .layout import motor_diameter_m
        # a wing nacelle is the motor's can and a cowl round it
        # (layout.motor_diameter_m); its width hides part of the slipstream
        # from the wing, as the fuselage does a nose propeller's
        from .layout import motor_diameter_m
        self.nacelle_width_m = (1.25 * motor_diameter_m(self.motor.mass_kg)
                                if self.n_wing_motors else 0.0)
        # ---- propulsion ----------------------------------------------
        # Free-air thrust until the airframe exists; the installed factor
        # (slipstream scrubbing) is set from it in the loop below.
        self.install_thrust_factor = 1.0
        self.prop_sys: PropulsionSystem = get_propulsion(
            self.motor, self.prop, design.n_motors, self.pack, env,
            self.install_thrust_factor, self.n_systems,
            b.harness_resistance_ohm, max_fuse_a=r.max_fuse_current_a,
            esc_airflow_mps=b.esc_airflow_mps,
            mission_s=r.mission_window_s, r_lead_ohm=self.lead_r_ohm,
            i0_speed_share=b.i0_speed_share)
        self._ps_free = self.prop_sys       # free-air thrust, for the tables
        # the ESC is hardware: every later variant (the installed-thrust
        # rebuild, hot windings) carries the one sized here
        self.esc = self.prop_sys.esc

        # ---- mass ----------------------------------------------------
        nm = design.n_motors
        esc = self.prop_sys.esc
        from .layout import motor_stations as _ms, disc_overlap_fraction
        _n0, _ys0, _share0 = _ms(nm, self.wing.span_m, body0.fuse_width_m,
                                 design.wing_prop_layout, self.prop.diameter_m,
                                 r.wingspan_counts_props)
        stag = design.wing_prop_layout == "stagger" and len(_ys0) >= 4
        n_ext = len(_ys0) // 2 if stag else 0
        propulsion_items = propulsion_mass_items(
            self.motor, self.prop, esc, self.pack, nm, self.n_systems,
            self.n_wing_motors, self.lead_mass_kg, self.gear_extra_mass_kg,
            n_ext, b)
        # Staggered wing props: every other nacelle runs 0.6 D further
        # forward (a carbon boom and its fairing, 100 g/m: in the items
        # above), and the aft props work partly in their neighbours'
        # slipstreams.
        self.stagger_extra_drag_m2 = 0.0
        self.stagger_thrust_factor = 1.0
        self._stagger = None                # (overlap share, jet development)
        if stag:
            d_p = self.prop.diameter_m
            l_ext = 0.6 * d_p
            # their skin: Cf 0.005 on the extra wetted length of the pod
            self.stagger_extra_drag_m2 = (n_ext * 0.005 * math.pi
                                          * self.nacelle_width_m * l_ext)
            # Each overlapped pair puts part of the aft disc in the forward
            # one's slipstream: that part flies faster air, so it makes
            # less thrust - linearly towards the zero-thrust speed V0 - at
            # the momentum-theory speed-up of the forward jet, 0.6 D
            # downstream - at every airspeed, full throttle
            # (_stagger_factor, folded into the installed-thrust table).
            phi = disc_overlap_fraction(_share0, d_p)
            per_side = len(_ys0) // 2
            phi_sys = 2.0 * (per_side - 1) * phi / max(nm, 1)
            r_p = 0.5 * d_p
            f_dev = 0.5 * (1.0 + l_ext / math.sqrt(r_p * r_p + l_ext * l_ext))
            self._stagger = (phi_sys, f_dev)
        self.propulsion_items = propulsion_items
        # Each wing pod hangs on the spar outboard of the root - motor,
        # propeller, ESC, mount and nacelle - and relieves the bending
        # inboard of it (weights.spar_design), in flight and in the tip
        # test.  The power leads run the span and are left out of it.
        from .layout import motor_stations
        _nn, _ys, _sh = motor_stations(nm, self.wing.span_m,
                                       body0.fuse_width_m,
                                       design.wing_prop_layout,
                                       self.prop.diameter_m,
                                       r.wingspan_counts_props)
        pod_kg = (self.motor.mass_kg + self.prop.mass_kg + esc.mass_kg
                  + b.motor_mount_kg * (motor_diameter_m(self.motor.mass_kg)
                                        / 0.0621) ** 2
                  + b.nacelle_mass_frac
                  * (self.motor.mass_kg + esc.mass_kg))
        self.wing_point_masses = [(y, pod_kg) for y in _ys if y > 0.0]

        # ---- layout, aero and mass, solved together --------------------
        # Nothing about where the parts sit is assumed.  The fuselage length
        # depends on the nose, the nose on the balance, the balance on every
        # mass and on the neutral point - and the masses and the neutral
        # point depend on the fuselage.  Alternate until the nose length and
        # the CG offset from the wing stop moving (layout.py explains the
        # rules the balance follows).
        nose, cg_used = 0.25, 0.0
        # tail volumes: fixed, or re-sized from the requirements each pass
        # until they settle with everything else (m111)
        vol_h, vol_v = (design.tail_vol_h if design.tail_vol_h > 0.0
                        else b.vol_coeff_h), b.vol_coeff_v
        self.tail_sizing: Dict[str, Any] = {"mode": b.tail_sizing}
        for _it in range(30):
            self.body = build_body(bay_l, bay_w, bay_h, self.wing, b,
                                   nose, cg_used, vol_h, vol_v)
            self.aero_clean = AeroModel(self.wing, self.body, env, b,
                                        design.n_motors,
                                        extra_drag_area_m2=gear_extra_drag
                                        + self.stagger_extra_drag_m2,
                                        prop_diameter_m=self.prop.diameter_m,
                                        nacelle_width_m=self.nacelle_width_m,
                                        wing_prop_layout=design.wing_prop_layout,
                                        props_in_span=r.wingspan_counts_props)
            self.aero_deployed = AeroModel(self.wing, self.body, env, b,
                                           design.n_motors,
                                           extra_drag_area_m2=gear_extra_drag
                                           + self.deployed_drag_area_m2
                                           + self.stagger_extra_drag_m2,
                                           prop_diameter_m=self.prop.diameter_m,
                                           nacelle_width_m=self.nacelle_width_m,
                                           wing_prop_layout=design.wing_prop_layout,
                                           props_in_span=r.wingspan_counts_props)
            # Installed thrust: what the slipstream scrubbing this airframe
            # costs (AeroModel.installed_thrust_factor), and the staggered
            # props' overlap loss, at every airspeed (_install_table).
            self.install_thrust_factor = self._install_table(self.aero_clean)
            self.prop_sys = get_propulsion(
                self.motor, self.prop, design.n_motors, self.pack, env,
                self.install_thrust_factor, self.n_systems,
                b.harness_resistance_ohm, max_fuse_a=r.max_fuse_current_a,
                esc=self.esc, esc_airflow_mps=b.esc_airflow_mps,
                mission_s=r.mission_window_s, r_lead_ohm=self.lead_r_ohm,
                i0_speed_share=b.i0_speed_share)
            n_override = b.ultimate_load_factor_override
            n_ult = n_override if n_override else 4.0
            q_dive = 0.0
            mb = None
            for _ in range(6):
                mb, empty, mtow_m2 = converge_empty_mass(
                    self.wing, self.body, propulsion_items, self.m2_payload_kg,
                    b, n_ult, q_dive_pa=q_dive,
                    point_masses=self.wing_point_masses)
                w = mtow_m2 * G0
                self.load_case: LoadCase = design_load_factors(
                    self.aero_clean, w, self._max_level_speed(w), b)
                n_new = (n_override if n_override
                         else self.load_case.n_ultimate)
                q_new = 0.5 * env.rho * self.load_case.v_dive ** 2
                if (abs(n_new - n_ult) < 1e-6
                        and abs(q_new - q_dive) < 1e-6 * max(q_new, 1.0)):
                    n_ult, q_dive = n_new, q_new
                    break
                n_ult, q_dive = n_new, q_new
            mb, empty, mtow_m2 = converge_empty_mass(
                self.wing, self.body, propulsion_items, self.m2_payload_kg,
                b, n_ult, q_dive_pa=q_dive,
                point_masses=self.wing_point_masses)
            self.n_ultimate = n_ult
            self.n_limit = (n_ult / b.load_safety_factor if n_override
                            else self.load_case.n_limit)
            self.q_dive_pa = q_dive

            self.mass_breakdown = mb
            cg_new = ((self.aero_clean.x_np_over_c - b.target_static_margin
                       - 0.25) * self.wing.mac_m)
            self.layout: Layout = solve_layout(self, cg_new)
            if b.tail_sizing == "requirements":
                vh_new, vv_new, self.tail_sizing = self._tail_volumes_required(
                    mtow_m2 * G0, (empty + self.m3_payload_kg) * G0,
                    empty * G0)
            else:
                vh_new, vv_new = vol_h, vol_v
            if (abs(self.layout.nose_length_m - nose) < 1e-6
                    and abs(cg_new - cg_used) < 1e-7
                    and abs(vh_new - vol_h) < 1e-6
                    and abs(vv_new - vol_v) < 1e-7):
                break
            nose, cg_used = self.layout.nose_length_m, cg_new
            vol_h, vol_v = vh_new, vv_new
        # The CG is where the parts put it - on target when the battery
        # can balance the aeroplane, forward of it when it cannot.
        lay = self.layout
        c_mac = self.wing.mac_m
        self.aero_clean.set_cg((lay.x_cg_m - lay.x_wing_le_mac_m) / c_mac)
        self.aero_deployed.set_cg((lay.x_cg_m3_m - lay.x_wing_le_mac_m)
                                  / c_mac)

        # Each propeller's full-throttle thrust and torque against airspeed
        # (cold windings, 70% charge): the swirl in a wing propeller's
        # slipstream and the roll trim against an unpaired motor are both
        # set by the torque (AeroModel.set_power_table).
        # from standstill to the propellers' speed ceiling
        _vc = self.prop_sys.speed_ceiling(0.7)
        _vs = [_vc * i / 30.0 for i in range(31)]
        _tt, _qq = [], []
        for _v in _vs:
            _op = self.prop_sys.solve(_v, 0.7, 1.0)
            _om = 2.0 * math.pi * _op.rpm / 60.0
            _tt.append(_op.thrust_total_n)
            _qq.append(_op.p_shaft_w / max(design.n_motors, 1) / _om
                       if _om > 1e-6 else 0.0)
        for _a in (self.aero_clean, self.aero_deployed):
            _a.set_power_table(_vs, _tt, _qq)

        # RULE 4.2.4 - the spar, sized for flight, the tip-lift test, shear
        # and torsion, each at the factor of safety, at MGTOW.
        self.spar: SparDesign = spar_design(
            self.wing, mtow_m2 * G0, b, n_ult, q_dive,
            point_masses=self.wing_point_masses)
        self.spar_t_req_m = self.spar.t_req_m
        self.spar_t_used_m = self.spar.t_used_m
        self.spar_gauge_limited = self.spar.gauge_limited

        self.mass_breakdown: MassBreakdown = mb
        self.empty_mass_kg = empty
        self.mtow_m2_kg = mtow_m2
        self.mtow_m3_kg = empty + self.m3_payload_kg
        self.mtow_m1_kg = empty            # Mission 1 flies empty

        self.w_m2_n = self.mtow_m2_kg * G0
        self.w_m3_n = self.mtow_m3_kg * G0
        self.w_m1_n = self.mtow_m1_kg * G0
        # for reports: the stagger loss at the Mission 2 top speed
        if self._stagger:
            self.stagger_thrust_factor = self._stagger_factor(
                self._max_level_speed(self.w_m2_n))

    def _stagger_factor(self, v: float) -> float:
        """Thrust kept at airspeed v by staggered wing props whose discs
        overlap (1 for every other layout).  The overlapped share of the
        aft disc works in the forward jet, faster by the momentum-theory
        increment du developed 0.6 D downstream, and its thrust falls
        linearly towards the zero-thrust speed V0:
            k = max(0.5, 1 - phi_sys du / (V0 - V)),
            du = f_dev V (sqrt(1 + T_p/(q A)) - 1)
        at full throttle, 70% charge, free-air thrust."""
        if not self._stagger:
            return 1.0
        phi_sys, f_dev = self._stagger
        d_p = self.prop.diameter_m
        nm = max(self.design.n_motors, 1)
        v_e = max(v, 0.5)
        op = self._ps_free.solve(v_e, 0.7, 1.0)
        r_p = 0.5 * d_p
        a_disc = math.pi * r_p * r_p
        q = 0.5 * self.cfg.env.rho * v_e * v_e
        du = f_dev * v_e * (math.sqrt(1.0 + op.thrust_total_n / nm
                                      / (q * a_disc)) - 1.0)
        v0 = self.prop.j_zero * (op.rpm / 60.0) * d_p
        return max(0.5, 1.0 - phi_sys * du / max(v0 - v_e, 1.0))

    def _install_table(self, aero) -> Tuple[Tuple[float, float], ...]:
        """Installed / free-air thrust against airspeed, 25 points from
        standstill to the propellers' speed ceiling: the slipstream
        scrubbing (AeroModel.installed_thrust_factor, whose skin friction
        falls with speed) times the stagger loss."""
        v_hi = self._ps_free.speed_ceiling(1.0)
        out = []
        for i in range(25):
            v = v_hi * i / 24.0
            k = (aero.installed_thrust_factor(self.prop.diameter_m, v)
                 * self._stagger_factor(v))
            out.append((round(v, 3), round(k, 5)))
        return tuple(out)

    def install_thrust_factor_at(self, v: float) -> float:
        """Installed / free-air thrust at airspeed v."""
        return self.prop_sys.installed_factor(v)

    def _max_level_speed(self, weight_n: float) -> float:
        """Maximum level speed, solved locally.

        mission.py imports this module, so it cannot be imported back; this
        is the same Brent solve on thrust minus drag.
        """
        aero = self.aero_clean

        def excess(v: float) -> float:
            t = self.prop_sys.thrust(v, 0.85, 1.0)
            return t - aero.drag(v, weight_n, t)      # power-on drag

        lo = aero.v_stall(weight_n) * 1.02
        # no thrust, so no level flight, above the propellers' ceiling
        hi = max(self.prop_sys.speed_ceiling(0.85), lo * 1.001)
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

    @property
    def depth_of_discharge(self) -> float:
        """Fraction of rated energy one mission may use: the team's setting,
        capped by the chemistry's ceiling."""
        return max(0.05, min(self.cfg.build.depth_of_discharge,
                             self.pack.usable_fraction))

    def fuse_peak_a(self) -> float:
        """The most current any one pack draws in this air - what each
        system's fuse must hold (RULE 3.2.3.2d).

        Full charge, full throttle, at the airspeed that draws the most
        (PropulsionSystem.rating_check), with the windings no warmer than
        the air: every flight starts with them at air temperature, and on
        a morning colder than the catalogue's 20 C their resistance is
        lower and the current higher than the cold-catalogue figure.  The
        pack only discharges and the windings only warm from there, so
        nothing later in a flight draws more."""
        i_pk = self.prop_sys.rating_check()["pack_current_a"]
        t_air = self.cfg.env.temp_c
        if t_air < 20.0:
            i_pk = max(i_pk, self.prop_sys_at(t_air).rating_check()
                       ["pack_current_a"])
        return i_pk

    def prop_sys_at(self, winding_c: float) -> PropulsionSystem:
        """The propulsion system with its windings at this temperature:
        copper resistance scaled from the catalogue's cold (20 C) value.
        Quantised to 2 C so a mission reuses a handful of solved systems."""
        b = self.cfg.build
        t = 2.0 * round(winding_c / 2.0)
        return get_propulsion(
            self.motor, self.prop, self.design.n_motors, self.pack,
            self.cfg.env, self.install_thrust_factor, self.n_systems,
            b.harness_resistance_ohm,
            1.0 + b.copper_temp_coeff * (t - 20.0),
            self.cfg.rules.max_fuse_current_a, self.esc,
            b.esc_airflow_mps, self.cfg.rules.mission_window_s,
            self.lead_r_ohm, b.i0_speed_share)

    def usable_energy_j(self) -> float:
        """Energy one mission may take out of every pack - a fresh, fully
        charged pack each mission - at the depth of discharge."""
        return self.total_energy_wh * 3600.0 * self.depth_of_discharge

    # ------------------------------------------------------------------
    # Hard constraints that do not need a mission simulation
    # ------------------------------------------------------------------
    def _tail_trim_need(self, w2: float, w3: float, w1: float
                        ) -> Tuple[float, str]:
        """The largest tail load per unit dynamic pressure (m2 of tail at
        CL 1, eta 1) the flight envelope asks for to trim, and where: each
        loading at its own CG (Mission 2, Mission 3 towing, Mission 1
        empty), at the stall at 1 g - in the linear trim model the same
        tail CL as the manoeuvre corner, n g at sqrt(n) x stall - and at
        the dive speed at 1 g, where the camber moment's download is
        largest (m111)."""
        a, lay, c = self.aero_clean, self.layout, self.wing.mac_m
        keep = a.x_cg_over_c
        worst, where = 0.0, ""
        try:
            for name, w, x in (("Mission 2", w2, lay.x_cg_m),
                               ("Mission 3", w3, lay.x_cg_m3_m),
                               ("Mission 1", w1, lay.x_cg_empty_m)):
                a.set_cg((x - lay.x_wing_le_mac_m) / c)
                for case, v in (("stall", a.v_stall(w)),
                                ("dive", self.load_case.v_dive)):
                    r = abs(a.tail_load(v, w)) / a.q(v)
                    if r > worst:
                        worst, where = r, "%s %s" % (name, case)
        finally:
            a.set_cg(keep)
        return worst, where

    def tail_trim_cl(self) -> Tuple[float, str]:
        """The horizontal tail's largest CL to trim across the flight
        envelope (_tail_trim_need), and where."""
        need, where = self._tail_trim_need(self.w_m2_n, self.w_m3_n,
                                           self.w_m1_n)
        return (need / (self.cfg.build.tail_eta
                        * max(self.body.s_horiz_m2, 1e-9)), where)

    def _tail_volumes_required(self, w2: float, w3: float, w1: float
                               ) -> Tuple[float, float, Dict[str, Any]]:
        """The horizontal and vertical tail volumes the requirements call
        for (m111; BuildStandard.tail_sizing), on the current CGs, neutral
        point, weights (Missions 2, 3 and 1) and propulsion - the m111
        notes in ASSUMPTIONS.md.  Returns (Vh, Vv, what each asked for)."""
        from .liftingline import lift_curve_slope
        b = self.cfg.build
        a, body, wing, lay = self.aero_clean, self.body, self.wing, self.layout
        s, c, span = wing.area_m2, wing.mac_m, wing.span_m
        lt = max(body.tail_arm_m, 0.05)
        eta = b.tail_eta
        # -- horizontal: trim across the envelope
        need, trim_at = self._tail_trim_need(w2, w3, w1)
        sh_trim = need / (eta * b.tail_trim_cl_frac * b.tail_cl_max)
        # -- horizontal: rotation about the main wheels at the liftoff
        # speed (mission.simulate_takeoff's), Mission 2, full throttle
        thr = lambda v: self.prop_sys.thrust(v, 1.0, 1.0)
        keep = a.x_cg_over_c
        a.set_cg((lay.x_cg_m - lay.x_wing_le_mac_m) / c)
        try:
            v_r = max(b.rotate_margin * a.v_stall_powered(w2, thr),
                      a.v_stall_powered(w2, thr, ground=True))
        finally:
            a.set_cg(keep)
        q_r = a.q(v_r)
        t_r = thr(v_r)
        g_l = a.slipstream_gains(v_r, t_r)[0] if t_r > 0.0 else 1.0
        d_mg = lay.x_main_gear_m - lay.x_cg_m           # CG ahead of wheels
        d_ac = lay.x_main_gear_m - lay.x_wing_ac_m      # wing AC ahead
        l_mt = max(lay.x_tail_ac_m - lay.x_main_gear_m, 0.05)
        h_cg = self.gear_height_m + 0.5 * body.fuse_height_m
        h_t = self.gear_height_m + self.prop_axis_m     # thrust line height
        l_w = q_r * s * b.takeoff_cl_ground * g_l
        m_ac = self.airfoil.cm_ac * q_r * s * c * g_l
        fric = b.rolling_friction * max(w2 - l_w, 0.0)
        # the roll's acceleration acts at the CG (d'Alembert): the thrust
        # pitches the nose down only by its line's height above the CG,
        # the friction that slows the roll by the CG's height; the drag, a
        # little above the CG, would help and is left out
        moment = (w2 * d_mg + t_r * (h_t - h_cg) + fric * h_cg
                  - l_w * d_ac - m_ac)
        sh_rot = max(moment, 0.0) / (eta * q_r * b.tail_cl_max * l_mt)
        sh_min = b.vol_coeff_h_min * s * c / lt
        # the team's choice may be bigger: less trim download (m112)
        sh_pick = max(self.design.tail_vol_h, 0.0) * s * c / lt
        sh = max(sh_trim, sh_rot, sh_min, sh_pick)
        # -- vertical: directional stability against the fuselage
        l_tc = min(body.tailcone_length_m, body.fuse_length_m)
        fus_vol = (body.fuse_width_m * body.fuse_height_m
                   * (body.fuse_length_m - l_tc + l_tc / 3.0))
        cnb_fus = (-1.3 * fus_vol / (s * span)
                   * body.fuse_height_m / max(body.fuse_width_m, 1e-6))
        a_v = lift_curve_slope(1.55 * 0.5 * b.tail_aspect_ratio)
        sv_stab = (max(b.cn_beta_target - cnb_fus, 0.0) * s * span
                   / (eta * a_v * lt))
        # -- vertical: one wing motor out
        sv_oei = 0.0
        ys = [abs(y) for y, _m in self.wing_point_masses]
        if b.fin_motor_out and self.n_wing_motors >= 2 and ys:
            v_mc = b.motor_out_speed_factor * a.v_stall(w2)
            q_mc = a.q(v_mc)
            t_one = thr(v_mc) / max(self.design.n_motors, 1)
            d_wm = (b.windmill_drag_coeff * q_mc * math.pi
                    * (0.5 * self.prop.diameter_m) ** 2)
            sv_oei = ((t_one + d_wm) * max(ys)
                      / (eta * q_mc * lt * b.fin_cl_max))
        sv_min = b.vol_coeff_v_min * s * span / lt
        sv = max(sv_stab, sv_oei, sv_min)
        vh, vv = sh * lt / (s * c), sv * lt / (s * span)
        h_by = max((sh_trim, "trim (%s)" % trim_at), (sh_rot, "rotation"),
                   (sh_min, "floor"), (sh_pick, "choice"))[1]
        v_by = max((sv_stab, "stability"), (sv_oei, "motor out"),
                   (sv_min, "floor"))[1]
        info = {"mode": "requirements", "vol_h": vh, "vol_v": vv,
                "s_h_trim_m2": sh_trim, "s_h_rotation_m2": sh_rot,
                "s_h_floor_m2": sh_min, "s_h_choice_m2": sh_pick,
                "s_v_stability_m2": sv_stab,
                "s_v_motor_out_m2": sv_oei, "s_v_floor_m2": sv_min,
                "h_sized_by": h_by, "v_sized_by": v_by,
                "v_rotate_mps": v_r, "rotation_moment_nm": moment,
                "cn_beta_fuselage": cnb_fus}
        return vh, vv, info

    def stall_air(self) -> Tuple[float, float, str]:
        """(density, kinematic viscosity, name) of the air the stall limit
        is checked in: BuildStandard.stall_check_conditions, or the
        mission's own when that is blank (m118)."""
        name = self.cfg.build.stall_check_conditions
        if not name:
            return self.cfg.env.rho, self.cfg.env.nu, "the mission's"
        env = Environment.at_site(name)
        return env.rho, env.nu, name

    def stall_speed_factored(self) -> float:
        """The slowest the aeroplane can hold its Mission 2 (heaviest)
        weight power off, with the lift it can use divided by BuildStandard.
        stall_lift_safety_factor, m/s - what max_stall_mps limits - in the
        air stall_check_conditions names (the hottest weather preset, m118).
        The lift it can use is its CLmax, or what its wing makes at the
        runway attitude (max_ground_alpha_deg) if that is less (m104)."""
        rho, nu, _name = self.stall_air()
        return self.aero_clean.v_stall_in_air(
            self.w_m2_n, rho, nu,
            load_factor=self.cfg.build.stall_lift_safety_factor, ground=True)

    def static_checks(self) -> Dict[str, Any]:
        r, b = self.cfg.rules, self.cfg.build
        rc = self.prop_sys.rating_check()
        d = self.design

        from .layout import (motor_stations, STAGGER_MAX_OVERLAP,
                             PROP_GAP_M)
        n_nose, ys, share = motor_stations(d.n_motors, self.wing.span_m,
                                           self.body.fuse_width_m,
                                           d.wing_prop_layout,
                                           self.prop.diameter_m,
                                           r.wingspan_counts_props)
        axis = [0.5 * self.body.fuse_height_m] if n_nose else []
        if ys:
            axis.append(self.body.fuse_height_m
                        if b.wing_covered_surfaces == 1
                        else 0.5 * self.body.fuse_height_m)
        clearance = self.gear_height_m + min(axis) - 0.5 * self.prop.diameter_m
        # Wing propellers: the tightest of the gap between neighbours on a
        # side (share - D) and the inner disc's clearance from the fuselage
        # side ((share - D) / 2).
        prop_gap = 0.5 * (share - self.prop.diameter_m) if ys else 1.0
        stagger_overlap = 0.0
        if ys and d.wing_prop_layout == "stagger" and len(ys) >= 4:
            # staggered: neighbours may overlap (fore and aft of each
            # other); only the inner disc must clear the fuselage
            prop_gap = (min(abs(y) for y in ys) - 0.5 * self.prop.diameter_m
                        - 0.5 * self.body.fuse_width_m)
            stagger_overlap = max(0.0, 1.0 - share / self.prop.diameter_m)

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
        # RULE 3.2.1g - the rules allow multiple motors and propellers
        # with no number; max_motors (4) is the team's cap.  One propeller
        # per motor, so it caps the propellers too.
        if d.n_motors > r.max_motors:
            fails.append("%d motors, limit is %d" % (d.n_motors, r.max_motors))
        # RULE 3.2.3.1a, 3.2.3.1b - one pack per system, and packs may
        # not be paralleled, so a second pack means a second full system.
        # RULE 3.2.3.1e - and never the other way round: one pack (even one
        # with two output leads) is never split across two systems, so
        # n_systems is always the number of separate packs carried.
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
        span_counted = self.wing.span_m + (
            2.0 * self.wing.winglet_h_m if b.winglets_count_in_span else 0.0)
        if span_counted > r.max_wingspan_m + 1e-9:
            fails.append("span %.4f m over the %.4f m limit"
                         % (span_counted, r.max_wingspan_m))
        # --- current --- RULE 3.2.3.2d, 3.2.3.1d ---
        # Each system's blade fuse is the hard electrical stop.  There is
        # no rule on the total across systems: each has its own fuse.
        i_fuse = self.fuse_peak_a()
        if i_fuse > self.prop_sys.fuse_rating_a + 1e-9:
            fails.append("blows its %.0f A fuse (%.0f A per pack)"
                         % (self.prop_sys.fuse_rating_a, i_fuse))
        # --- hardware ratings --------------------------------------------
        # The motor's rated cell count.  Its current and power limits are
        # checked below, but they say nothing about the winding insulation
        # or the rotor's speed limit, which is what a voltage rating is.
        if self.pack.series > self.motor.max_cells:
            fails.append("%dS pack on a motor rated for %dS at most"
                         % (self.pack.series, self.motor.max_cells))
        if self.pack.series > self.prop_sys.esc.max_cells:
            fails.append("%dS pack over the ESC's %dS rating"
                         % (self.pack.series, self.prop_sys.esc.max_cells))
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
        if rc["prop_rpm_margin"] < 1.0:
            fails.append("propeller over its maker's %.0f rpm limit (%.0f "
                         "rpm unloaded)" % (self.prop.max_rpm,
                                            rc["rpm_no_load"]))
        # --- structure --- RULE 4.2.4 ---
        # Every wing load case at the factor of safety, everywhere on the
        # span.  The spar is sized to meet it, so this fails only when the
        # tube it needs cannot exist: a wall thicker than its own radius.
        if not self.spar.closes:
            fails.append("wing spar cannot carry the %.2fx factor of safety "
                         "at %.0f%% span - the wall it needs is thicker "
                         "than the section is deep"
                         % (self.spar.fos_required, 100 * self.spar.worst_eta))
        if self.gear_height_m > 0.30:
            fails.append("landing gear would need %.0f mm legs for prop clearance"
                         % (self.gear_height_m * 1000))
        if prop_gap < PROP_GAP_M - 1e-9:
            fails.append("propellers overlap each other or the fuselage")
        if stagger_overlap > STAGGER_MAX_OVERLAP + 1e-9:
            fails.append("staggered propellers overlap %.0f%% of a diameter "
                         "(at most %.0f%%)" % (100 * stagger_overlap,
                                              100 * STAGGER_MAX_OVERLAP))
        # The team's stall limit (BuildStandard.max_stall_mps): power off
        # at the Mission 2 weight, on the lift over its factor of safety.
        if b.max_stall_mps < float("inf"):
            v_sf = self.stall_speed_factored()
            if v_sf > b.max_stall_mps + 1e-9:
                rho_s, nu_s, air = self.stall_air()
                fails.append("stall %.1f mph on lift / %.2f in %s air, over "
                             "the %.1f mph limit (%.1f mph true)"
                             % (v_sf * MS2MPH, b.stall_lift_safety_factor, air,
                                b.max_stall_mps * MS2MPH,
                                self.aero_clean.v_stall_in_air(
                                    self.w_m2_n, rho_s, nu_s) * MS2MPH))
        # the tail has to be able to trim without stalling itself, across
        # the flight envelope (tail_trim_cl, m111)
        ct, ct_at = self.tail_trim_cl()
        if ct > b.tail_cl_max:
            fails.append("tail needs CL %.2f to trim at the %s (max %.2f) - "
                         "bigger tail or longer arm" % (ct, ct_at, b.tail_cl_max))
        sm = self.aero_clean.static_margin
        if sm < b.min_static_margin:
            fails.append("static margin %.3f c - %s" % (
                sm, "unstable" if sm < 0 else "too little for a model"))
        if sm > b.max_static_margin:
            fails.append("static margin %.3f c - excessive, all trim drag" % sm)
        # Mission 1 flies empty.  With the payload on the CG the margin is
        # the same; when the layout could not balance exactly it is not.
        sm1 = self.layout.static_margin_empty
        if not b.min_static_margin <= sm1 <= b.max_static_margin:
            fails.append("static margin empty (Mission 1) %.3f c, outside "
                         "%.2f-%.2f" % (sm1, b.min_static_margin,
                                        b.max_static_margin))
        if self.aspect_ratio < 1.0 or self.aspect_ratio > 16.0:
            fails.append("aspect ratio %.1f outside modelled range" % self.aspect_ratio)
        if self.wing.chord_m < 0.12:
            fails.append("chord too small (low Re)")
        # the team's chord limit (BuildStandard.max_chord_m, m106): the
        # root chord, the wing's largest
        if self.wing.root_chord_m > b.max_chord_m + 1e-9:
            fails.append("chord %.1f in over the %.1f in limit"
                         % (self.wing.root_chord_m / 0.0254,
                            b.max_chord_m / 0.0254))
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
        """Will the wing hold together, with the factor of safety?

        The spar is sized so that every wing load case - flight, the
        tip-lift test, shear and torsion - sits at or above the factor of
        safety everywhere on the span, so the answer is yes by
        construction - unless the sizing cannot close, which happens in
        exactly one way: somewhere along the span the wall it asks for is
        thicker than the tube's own radius.  At that point it is not a
        tube any more, and no amount of layup fixes it; the wing is too
        long, too thin or too heavily loaded.
        """
        return self.spar.closes

    # ------------------------------------------------------------------
    def tail_dimensions(self) -> Dict[str, float]:
        """What is built (m115), in metres: the stabiliser at
        tail_aspect_ratio and the fin at half of it, as the aerodynamics
        and the weights take them, both at tail_t_over_c; stations from
        the nose tip (both tails' quarter chords at the tail AC)."""
        b, lay = self.cfg.build, self.layout
        s_h, s_v = self.body.s_horiz_m2, self.body.s_vert_m2
        c_h = math.sqrt(s_h / b.tail_aspect_ratio)
        c_v = math.sqrt(s_v / (0.5 * b.tail_aspect_ratio))
        return {"tail_area_m2": s_h, "tail_span_m": s_h / c_h,
                "tail_chord_m": c_h, "tail_thickness_m": b.tail_t_over_c * c_h,
                "fin_area_m2": s_v, "fin_height_m": s_v / c_v,
                "fin_chord_m": c_v, "fin_thickness_m": b.tail_t_over_c * c_v,
                "tail_arm_m": self.body.tail_arm_m,
                "fuse_length_m": self.body.fuse_length_m,
                "fuse_width_m": self.body.fuse_width_m,
                "fuse_height_m": self.body.fuse_height_m,
                "x_wing_le_m": lay.x_wing_le_root_m,
                "x_tail_le_m": lay.x_tail_ac_m - 0.25 * c_h,
                "x_fin_le_m": lay.x_tail_ac_m - 0.25 * c_v,
                "x_cg_m": lay.x_cg_m}

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
            "winglet_h_m": self.design.winglet_h_m,
            "tail_vol_h_chosen": self.design.tail_vol_h,
            "vol_coeff_h": self.body.s_horiz_m2 * self.body.tail_arm_m
            / (self.wing.area_m2 * self.wing.mac_m),
            "vol_coeff_v": self.body.s_vert_m2 * self.body.tail_arm_m
            / (self.wing.area_m2 * self.wing.span_m),
            **self.tail_dimensions(),
            "tail_trim_cl": self.tail_trim_cl()[0],
            "tail_h_sized_by": self.tail_sizing.get("h_sized_by", "volume"),
            "tail_v_sized_by": self.tail_sizing.get("v_sized_by", "volume"),
            "AR_eff": self.aero_clean.AR_eff,
            "stall_m2_mps": self.aero_clean.v_stall(self.w_m2_n),
            "stall_factored_mps": self.stall_speed_factored(),
            "cl_at_ground_alpha": self.aero_clean.cl_at_ground_alpha(),
            "tail_cl_cruise": self.aero_clean.tail_cl(
                self.aero_clean.v_stall(self.w_m2_n) * 1.3, self.w_m2_n),
            "static_thrust_n": rc["static_thrust_n"],
            "static_current_a": rc["static_current_a"],
            "peak_current_a": rc["total_current_a"],
            "static_T_over_W_m2": rc["static_thrust_n"] / self.w_m2_n,
            "wing_loading_Nm2": self.wing_loading_m2(),
            "n_limit": self.n_limit,
            "n_ultimate": self.n_ultimate,
            "load_driver": getattr(self.load_case, "driver", "override"),
            "spar_gauge_limited": self.spar_gauge_limited,
            "structure_ok": self.structure_ok,
        }
