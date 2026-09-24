"""Environment, course, rules and build-standard configuration.

Rule values here were supplied by the team for the 2026-27 season.  Items
still tagged `# VERIFY` are ones nobody has confirmed against the rulebook.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from .units import FT2M, IN2M, LB2KG


# --------------------------------------------------------------------------
# Atmosphere / field
# --------------------------------------------------------------------------
@dataclass
class Environment:
    rho: float = 1.225            # kg/m^3 air density
    nu: float = 1.46e-5           # m^2/s kinematic viscosity
    g: float = 9.80665            # m/s^2
    wind_mps: float = 4.5         # steady wind along the runway
    field_elev_m: float = 0.0

    @property
    def density_altitude_ft(self) -> float:
        """The ISA altitude at which the air is this thin, in feet.

        The number a pilot actually talks about.  Inverting the standard
        atmosphere for density: rho/rho0 = (1 - 6.87535e-6 h)^4.2559,
        with h in feet, so h = (1 - (rho/rho0)^(1/4.2559)) / 6.87535e-6.
        """
        ratio = max(self.rho, 1e-6) / 1.225
        return (1.0 - ratio ** (1.0 / 4.2559)) / 6.87535e-6

    @staticmethod
    def at_density_altitude(alt_m: float, temp_c: float = 15.0,
                            wind_mps: float = 4.5) -> "Environment":
        """ISA density for a field elevation + temperature, which is what
        actually bites you at the competition site."""
        t_k = temp_c + 273.15
        p = 101325.0 * (1.0 - 2.25577e-5 * alt_m) ** 5.25588
        rho = p / (287.058 * t_k)
        mu = 1.458e-6 * t_k ** 1.5 / (t_k + 110.4)
        return Environment(rho=rho, nu=mu / rho, wind_mps=wind_mps,
                           field_elev_m=alt_m)


# --------------------------------------------------------------------------
# Course
# --------------------------------------------------------------------------
Segment = Tuple[str, float]       # ("straight", metres) | ("turn", degrees)


@dataclass
class Course:
    """2026-27 course as described by the team:

        two 1000 ft straightaways, a 180 degree turn at each end, and a
        horizontal 360 degree loop in the middle of the second straightaway.

    So the second straight is flown as 500 ft, the 360, then another 500 ft.
    Total straight distance 2000 ft; total turning 720 degrees.

    Turn *time* is `angle / turn_rate`, and turn rate is `g sqrt(n^2-1)/V`.
    That is a function of how hard each aeroplane can turn, so the lap time
    genuinely depends on each design's own best turning method - exactly as
    it should.  The ground track length of a turn never enters the time.
    """
    # RULE 3.3.1: "Upwind and downwind markers will be 500 ft from the
    # starting line."  The start/finish line is in the MIDDLE, so a lap is
    # 500 ft out to the upwind marker, a 180, 1000 ft back past the line to
    # the downwind marker with the 360 at midfield, a 180, then 500 ft home.
    #
    # Straights total 2000 ft and turning totals 720 degrees either way, but
    # the split matters: the aeroplane accelerates out of four turns per lap
    # rather than three, and it never gets one long 1000 ft run to build up
    # speed on.
    segments: List[Segment] = field(default_factory=lambda: [
        ("straight", 500.0 * FT2M),    # start/finish line to upwind marker
        ("turn", 180.0),               # upwind turn
        ("straight", 500.0 * FT2M),    # back down to the start/finish line
        ("turn", 360.0),               # horizontal 360 at midfield
        ("straight", 500.0 * FT2M),    # on to the downwind marker
        ("turn", 180.0),               # downwind turn
        ("straight", 500.0 * FT2M),    # home to the start/finish line
    ])
    headings: List[int] = field(
        default_factory=lambda: [+1, 0, -1, 0, -1, 0, +1])
    pattern_altitude_m: float = 30.0

    @property
    def straight_length_m(self) -> float:
        return sum(L for kind, L in self.segments if kind == "straight")

    @property
    def total_turn_deg(self) -> float:
        return sum(L for kind, L in self.segments if kind == "turn")


# --------------------------------------------------------------------------
# Rules / scoring constants
# --------------------------------------------------------------------------
@dataclass
class Rules:
    # --- flight windows -----------------------------------------------------
    # RULE 3.3.3.1, 3.3.3.2, 3.3.3.3 - a 5 minute window for every
    # flight mission; M1 needs 3 laps and M2 needs 5.
    mission_window_s: float = 300.0
    m1_required_laps: int = 3
    m2_required_laps: int = 5

    # --- field constraints --------------------------------------------------
    # No takeoff distance limit in 2026-27.  The ground roll is still computed
    # and reported, because a 300 ft roll is a real operational risk even when
    # it is not a rules violation.
    takeoff_distance_m: float = float("inf")   # RULE 3.2.1b
    # RULE 3.1
    # Maximum allowable wingspan: 6 ft.  This is a real rule, and it is the
    # constraint that stops the wing growing - without it, area scales as
    # span squared while the spar only scales linearly, and the optimiser
    # walks off to a 10 ft aeroplane.
    max_wingspan_m: float = 6.0 * FT2M

    # --- mass ---------------------------------------------------------------
    # RULE 3.2.1a, 4.2.4
    # 55 lb is the hard gross-weight ceiling.  Anything heavier is simply
    # not a legal aeroplane, so it is rejected outright rather than scored.
    max_gross_mass_kg: float = 55.0 * LB2KG

    # --- energy -------------------------------------------------------------
    # RULE 3.2.3.2b - total propulsion energy summed over every pack.
    max_watt_hours: float = 100.0
    # RULE 3.2.3.2c - and each individual pack must also clear the FAA
    # hand-carry limit: the same number, but applied per pack, not to the
    # sum.  Two legal packs can still break the airframe total.
    max_pack_watt_hours: float = 100.0
    # RULE 3.2.3.2a - NiCd and NiMH are legal but are not searched: at
    # DBF power levels they lose on both mass and current.
    battery_chemistry_allowed: Tuple[str, ...] = ("LiPo", "LiIon")

    # --- propulsion system --------------------------------------------------
    # RULE 3.2.3.1a, 3.2.3.1b, 3.2.3.1c
    # A propulsion system is one battery, one fuse, one arming plug,
    # one or more ESCs and one or more motors.  At most one battery pack per
    # system, and packs may not be wired in series or parallel - so the only
    # way to carry more than one pack is to build more than one complete,
    # independent system, and every pack must then be identical.
    max_motors: int = 2
    # RULE 3.2.3.2d
    # The fuse is a 100 A blade fuse at most, and it is rated at the pack's
    # maximum continuous discharge current.  100 A is therefore the hard
    # ceiling on current drawn through any one propulsion system.
    max_fuse_current_a: float = 100.0
    # Total draw across all systems, which is what the team specified.
    max_total_current_a: float = 100.0

    # --- payload ------------------------------------------------------------
    # No rules cap on sensor or simulator mass, so these are search bounds
    # wide enough that the physics, not the bound, decides.  If a result pins
    # to one of them, the report says so.
    # No rules cap on sensor mass.  This is only a search box, set wide
    # enough that the physics decides; if a result pins here, the report
    # says so and the bound should be widened again.
    sensor_mass_min_kg: float = 0.15
    sensor_mass_max_kg: float = 20.00
    sensor_deployed_drag_area_m2: float = 0.0045   # CD*A while deployed

    # Shipping container: ONE fixed assumed box for every design, rather
    # than resizing it per sensor.  The dimensions come from the drop-test
    # calculation for a roughly 1 kg sensor: a 5 ft drop reaches 5.47 m/s,
    # holding the sensor under 75 g needs a 20 mm crush stroke, and foam
    # crushes usefully through about 70 percent of its thickness, so 29 mm
    # of foam per face around a 165 x 83 x 83 mm sensor envelope.  Rounded
    # to a buildable box.  Change these four numbers if the real container
    # differs; nothing else has to change.
    container_length_m: float = 0.230
    container_width_m: float = 0.150
    container_height_m: float = 0.150
    container_tare_kg: float = 0.5 * LB2KG
    # Nominal packing density of a sensor built out of electronics, battery
    # and structure.  Used to REPORT the density a given sensor mass would
    # need, not to cap it: a heavy sensor is simply a dense one.
    sensor_packing_density_kg_m3: float = 800.0
    container_fill_fraction: float = 0.55   # usable internal volume fraction
    # The only real ceiling is physics.  A sensor cannot be denser than the
    # dense metal you would ballast it with; lead is 11340 kg/m3.  This
    # should essentially never bind - it exists so the optimiser cannot
    # propose a 50 lb brick in a 9 inch box.
    max_sensor_density_kg_m3: float = 11340.0

    # Simulators have the container's external dimensions (that is what they
    # simulate) and whatever mass you load them to.  No rules cap on either
    # the count or the mass, so both are physics-limited.
    max_simulators: int = 40
    sim_mass_max_kg: float = 3.00

    # --- Ground Mission --- RULE 3.3.4 -------------------------------------
    # GM = 0.5 + (sensor weight x drop height) / Max(...).  Six drops, two
    # on each of three faces, onto concrete, at a height the team declares
    # in whole inches up to 60.  The sensor must be at its declared maximum
    # weight, and any physical damage fails the mission.
    #
    # This is a scored mission and it is part of the total, so it belongs
    # in the model: it makes sensor weight pay in THREE places, not two.
    max_drop_height_in: float = 60.0
    # RULE 3.1.3d - rigging ring to the bottom of the container must be
    # under 15 in or the fixture cannot reach the full drop height.
    max_rigging_plus_container_in: float = 15.0
    # Deceleration the sensor is designed to survive, in g.  The container
    # foam is sized to hold the sensor under this on the declared drop.
    sensor_survivable_g: float = 75.0

    # --- sensor --- RULE 3.1.1a --------------------------------------------
    # The sensor is a towed body, deployed on a tow line and recovered
    # before landing.  It is not a free-form lump of mass.
    sensor_min_length_in: float = 6.0          # RULE 3.1.1b
    # RULE 3.1.2 - the tow line must deploy to at least 1.5 x wingspan,
    # measured from the exit to the sensor's forward tip.  That is a long
    # piece of line in the airstream, and it is charged as drag.
    tow_line_span_multiple: float = 1.5
    tow_line_diameter_m: float = 0.0015    # ~1.5 mm braided line
    tow_line_cd: float = 1.10              # circular cylinder, subcritical Re

    # --- normalisers --------------------------------------------------------
    # "self" normalises each mission against the best design found in this
    # simulation, which is what the team asked for.  "field" uses the two
    # estimates below instead, for when you know what the competition can do.
    normalisation: str = "self"
    field_best_m2_lb_per_s: float = 0.060
    field_best_m3_lap_lb: float = 22.0
    # Ground Mission: sensor weight (lb) x declared drop height (in).
    # A 2 lb sensor dropped the full 60 in is 120 lb*in, so this is a
    # deliberately modest guess at what the field will manage.
    field_best_gm_lb_in: float = 120.0


# --------------------------------------------------------------------------
# Build standard (what your shop actually produces)
# --------------------------------------------------------------------------
@dataclass
class BuildStandard:
    # Areal densities: the single most valuable thing to calibrate against
    # your own shop.  Weigh last year's parts and divide by area.
    wing_areal_density: float = 1.10      # kg/m^2 of planform: ribs+skin+covering
    # Horizontal and vertical tails are separate numbers because they are
    # built differently: the elevator surface is usually a flat plate with
    # one hinge line, while the fin carries the rudder, the tailwheel
    # steering and often the antenna, so it comes out heavier per unit
    # area.  Both include hinges, horns and linkages.
    #
    # NOTE these are EFFECTIVE densities, and at 1.80 they are well above
    # the wing's 1.10 kg/m^2 - which is backwards for most builds, where a
    # tail is lighter per unit area than a wing.  The value is inherited
    # from an earlier model that applied a blanket x2 factor here; it is
    # preserved so results do not move silently.  Weigh your own tails and
    # replace it.
    h_tail_areal_density: float = 1.80    # kg/m^2 of planform (elevator)
    v_tail_areal_density: float = 1.80    # kg/m^2 of planform (rudder)
    fuse_areal_density: float = 1.40      # kg/m^2 of wetted area
    spar_sigma_allow_pa: float = 600e6    # carbon tube working stress
    spar_rho: float = 1600.0              # kg/m^3
    spar_min_wall_m: float = 0.0006       # minimum practical layup/tube wall
    spar_fitting_factor: float = 1.8      # webs, joiners, root fittings

    # The ultimate load factor is DERIVED per design from a V-n diagram
    # (see loads.py): the worse of the manoeuvre case and the Pratt gust
    # case, capped by what the wing can aerodynamically generate at the
    # design dive speed, times a 1.5 safety factor.  Set the override to a
    # number if tech inspection specifies a wing test instead.
    ultimate_load_factor_override: Optional[float] = None
    load_safety_factor: float = 1.50       # ultimate / limit, standard
    load_factor_floor: float = 2.50        # sanity floor on limit load
    pilot_overshoot_factor: float = 1.15   # a pilot flying for time overshoots
    dive_speed_factor: float = 1.25        # V_D / V_cruise
    gust_cruise_mps: float = 25.0 * FT2M   # rough-air discrete gust at Vc
    gust_dive_mps: float = 12.5 * FT2M     # half of it at V_D, per FAR 23

    gear_mass_frac: float = 0.030         # of MTOW
    gear_mass_fixed_kg: float = 0.080
    servo_mass_kg: float = 0.020
    n_servos_base: int = 6                # 2 ail, 1 elev, 1 rud, 1 door, 1 sensor
    avionics_kg: float = 0.100            # rx, wiring, telemetry
    # 3.2.3: a SEPARATE receiver/servo battery is required, and it does not
    # count against the 100 Wh propulsion budget - but it does have mass.
    # A 2S 1000 mAh LiFe/LiPo Rx pack plus its switch harness.
    rx_battery_mass_kg: float = 0.075
    # Per propulsion system: one blade fuse in its own harness, plus one
    # externally accessible arming plug.
    fuse_harness_mass_kg: float = 0.045
    arming_plug_mass_kg: float = 0.030
    deploy_mech_kg: float = 0.150         # winch/door for the M3 sensor
    payload_structure_frac: float = 0.080 # restraint structure per kg of payload
    contingency_frac: float = 0.080       # honest margin on the whole airframe

    # drag areas (CD*A, m^2)
    gear_drag_area_m2: float = 0.0040
    nacelle_drag_area_m2: float = 0.0015  # per motor beyond the first
    misc_interference_frac: float = 0.10

    # fuselage proportions
    fuse_nose_m: float = 0.22
    fuse_tailcone_m: float = 0.55
    fuse_min_length_m: float = 0.90
    tail_arm_frac: float = 0.45           # of fuselage length
    vol_coeff_h: float = 0.55
    vol_coeff_v: float = 0.040
    tail_aspect_ratio: float = 4.0
    # The tail sits in the wing's wake and behind the fuselage, so it sees
    # less dynamic pressure than the wing.  0.9 is conventional.
    tail_eta: float = 0.90
    # Stability.  The neutral point is now COMPUTED, so static margin is an
    # output and these are real limits on it.
    min_static_margin: float = 0.05     # below this it is not flyable
    max_static_margin: float = 0.35     # above this it is a dart, all trim drag
    # The fuselage is destabilising: it moves the neutral point forward.
    # Raymer-style term, calibrated so a typical DBF pod shifts it ~0.05c.
    fuselage_np_shift_k: float = 0.60
    # CG position aft of the wing aerodynamic centre, in chords.  Positive
    # means the CG is behind the wing AC, which unloads the tail.  A
    # conventional stable layout sits a little aft of the wing AC and well
    # forward of the neutral point; 0.05 c is typical.  Set 0.0 to model
    # trim from wing camber alone.
    cg_aft_of_ac_frac: float = 0.05
    tail_cl_max: float = 1.10        # symmetric section plus elevator
    tail_efficiency: float = 0.85    # tail span efficiency for induced drag

    # --- realism factors: these exist to stop the model flattering itself ---
    cl_max_3d_factor: float = 0.82        # built wing vs 2-D section data
    # A freewheeling propeller is a disc of drag once the power is cut.
    # Coefficient on disc area, per propeller.  With measured propeller
    # data (which extends past J_zero into negative thrust) this could be
    # computed instead of assumed.
    windmill_drag_coeff: float = 0.08
    # Nobody flies the perfect corner speed on every lap.  Lap time is
    # divided by this.  0.93 means 7.5 percent off a flawless lap; set it
    # from your own flight logs, and set it to 1.0 to model a perfect pilot.
    pilot_efficiency: float = 0.93
    install_thrust_factor: float = 0.95   # blockage, non-uniform inflow

    # operating
    rolling_friction: float = 0.035       # paved runway, pneumatic wheels
    braking_friction: float = 0.30        # wheel brakes / drag on rollout
    takeoff_cl_ground: float = 0.35
    rotate_margin: float = 1.15           # V_rotate / V_stall
    approach_margin: float = 1.30         # V_approach / V_stall
    descent_rate_mps: float = 3.0
    flare_time_s: float = 3.0
    max_bank_deg: float = 70.0
    # Fallback only: the mission uses the derived LIMIT load factor, because
    # you should not exceed limit load in normal flight.
    max_turn_load_factor: float = 3.5
    prop_clearance_margin_m: float = 0.05
    gear_height_m: float = 0.130          # main gear leg length


@dataclass
class Config:
    env: Environment = field(default_factory=Environment)
    course: Course = field(default_factory=Course)
    rules: Rules = field(default_factory=Rules)
    build: BuildStandard = field(default_factory=BuildStandard)
