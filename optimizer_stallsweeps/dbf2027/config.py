"""Environment, course, rules and build-standard configuration.

Rule values here were supplied by the team for the 2026-27 season.  Items
still tagged `# VERIFY` are ones nobody has confirmed against the rulebook.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .units import FT2M, IN2M, LB2KG


# --------------------------------------------------------------------------
# Atmosphere / field
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Contest site
# --------------------------------------------------------------------------
# Rules 2.4 and the schedule table in section 2: the 2026-27 fly-off is at
# the Tucson International Modelplex Park Association (TIMPA), Tucson AZ,
# 15-18 April 2027, hosted by RTX.
#
# These are CLIMATE NORMALS for Tucson in mid-April, not a forecast - no
# forecast exists for 2027, and the rules themselves point teams at
# historical weather rather than a prediction.  They are a reasonable
# design point, not a promise about the day.
#
# ELEVATION - measured, not assumed.  TIMPA sits west of the city in the
# Avra Valley, well below Tucson International (2,643 ft).  The field is
# at 32.2647 N, 111.2740 W; the USGS 3DEP point-elevation service returns
# 667.06 m there (2,188 ft), on a 1 m raster.  The parcel is 160 acres and
# flat, so the runway is within a few metres of that.
#
#   https://epqs.nationalmap.gov/v1/json?x=-111.273994&y=32.264718&units=Meters
#
# TEMPERATURE - NWS Tucson / Davis-Monthan climate normals, 1992-2021:
# April mean daily high 83 F (28.3 C), mean daily low 53 F (11.7 C).  The
# three conditions below are built from those two numbers plus a hot day
# that is warm for April without being a record (April records reach about
# 100 F / 38 C).
#
# WIND - April is Tucson's windiest month, mean hourly 8.6 mph (3.8 m/s).
# The nominal 5 m/s is above that mean because the wind peaks in the
# afternoon, when most rounds are flown; mornings are calmer.
SITE_NAME = "TIMPA, Tucson AZ (15-18 April 2027)"
SITE_ELEVATION_M = 667.0          # 2,188 ft - USGS 3DEP at 32.2647N 111.2740W

#: Named conditions: (field elevation m, temperature C, wind m/s).
#: Tucson in April is dry and, by local standards, windy; the spread
#: through a flying day is wide enough that a design which only works at
#: one of these is not really finished.  Density altitudes, at 667 m:
#:
#:   morning   12 C   ~2,350 ft    rho 1.143   (April normal low, 53 F)
#:   nominal   28 C   ~4,170 ft    rho 1.082   (April normal high, 83 F)  <- default
#:   hot       36 C   ~5,040 ft    rho 1.054   (warm April day, not a record)
#:   ISA sea level    0 ft         rho 1.225   (for comparison only)
#:
#: So flying at Tucson costs roughly 12 percent of sea-level density on a
#: nominal afternoon.  That is not a detail: it is less lift at a given
#: speed, less thrust from the same propeller, and a higher stall speed.
CONDITIONS: Dict[str, Tuple[float, float, float]] = {
    "tucson-morning": (SITE_ELEVATION_M, 12.0, 3.0),
    "tucson": (SITE_ELEVATION_M, 28.0, 5.0),
    "tucson-hot": (SITE_ELEVATION_M, 36.0, 7.0),
    "isa": (0.0, 15.0, 4.5),
}
DEFAULT_CONDITIONS = "tucson"


@dataclass
class Environment:
    # Defaults are the Tucson nominal above, NOT sea-level ISA: the
    # aeroplane is going to be flown in Arizona in April, and designing
    # it at sea level then discovering it 4,000 ft of density altitude
    # later is the wrong order to do things in.
    # These two must equal at_site("tucson"); selftest asserts it, because
    # a dataclass default cannot call the function that derives it and
    # hand-copied constants drift.
    rho: float = 1.082329         # kg/m^3 air density  (667 m, 28 C)
    nu: float = 1.710606e-5       # m^2/s kinematic viscosity
    g: float = 9.80665            # m/s^2
    wind_mps: float = 5.0         # steady wind along the runway
    field_elev_m: float = 667.0
    # Air temperature: the motors start the mission at it and cool to it.
    temp_c: float = 28.0

    @property
    def density_altitude_ft(self) -> float:
        """The ISA altitude at which the air is this thin, in feet.

        The number a pilot actually talks about.  Inverting the standard
        atmosphere for density: rho/rho0 = (1 - 6.87535e-6 h)^4.2559,
        with h in feet, so h = (1 - (rho/rho0)^(1/4.2559)) / 6.87535e-6.
        """
        ratio = max(self.rho, 1e-6) / 1.225
        # 1 - r^(1/4.2559) as -expm1(ln r / 4.2559): no cancellation near
        # sea level, where r -> 1
        return -math.expm1(math.log(ratio) / 4.2559) / 6.87535e-6

    @staticmethod
    def at_site(name: str = DEFAULT_CONDITIONS) -> "Environment":
        """One of the named conditions in `CONDITIONS`.

        Unknown names fall back to the default rather than raising, so a
        typo in a saved file costs you the right weather, not the run.
        """
        alt, temp, wind = CONDITIONS.get(name,
                                         CONDITIONS[DEFAULT_CONDITIONS])
        return Environment.at_density_altitude(alt, temp, wind)

    @staticmethod
    def at_density_altitude(alt_m: float, temp_c: float = 27.0,
                            wind_mps: float = 5.0) -> "Environment":
        """ISA density for a field elevation + temperature, which is what
        actually bites you at the competition site."""
        t_k = temp_c + 273.15
        p = 101325.0 * (1.0 - 2.25577e-5 * alt_m) ** 5.25588
        rho = p / (287.058 * t_k)
        mu = 1.458e-6 * t_k ** 1.5 / (t_k + 110.4)
        return Environment(rho=rho, nu=mu / rho, wind_mps=wind_mps,
                           field_elev_m=alt_m, temp_c=temp_c)


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
    # RULE 3.1 says only "Maximum allowable wingspan is 6 feet": nothing on
    # how it is measured, nor that propellers or nacelles count.  Taken
    # literally it is the wing's own tip-to-tip span, so a wingtip
    # propeller may sit with its axis on the tip and its disc overhanging
    # (layout.motor_stations).  True counts the discs in the span instead:
    # their edge then stops at the tip.
    wingspan_counts_props: bool = False

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
    # RULE 3.2.1g - "The airplane may use multiple motors and/or
    # propellers", with no number.  4 is the TEAM's cap (one propeller
    # per motor, so 4 propellers too).
    max_motors: int = 4
    # RULE 3.2.3.2d
    # RULE 3.2.3.2d - "The maximum current rating for the fuse is the
    # maximum continuous discharge current rating of the battery pack
    # (battery capacity * C-rating or maximum labeled discharge current)
    # up to 100 amps."  So each system's fuse is min(pack rating, this),
    # and that is the ceiling on the current through any one system.
    max_fuse_current_a: float = 100.0
    # There is NO rules limit on the current summed over every system.
    # Each system has its own fuse (3.2.3.1), and 3.2.3.1e only stops one
    # pack being split across two systems - so two systems, each with its
    # own identical pack, may draw up to 100 A each.  (A 100 A total that
    # was here was a team limit, not a rule; it has been removed.)
    # There is NO rules limit on voltage.  The only ceilings in the model
    # are hardware ratings - each motor's and ESC's cell count
    # (Motor.max_cells, ESC.max_cells) - and the catalogue, which stops at
    # 12S because nothing COTS above that fits under 100 Wh.

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
    # Tube allowables at ULTIMATE load.  The loads already carry the 1.5
    # factor of safety, so these are failure stresses, not knocked down
    # again: 600 MPa is the compression side of a roll-wrapped carbon
    # tube in bending (the tension side is stronger), and 50 MPa is a
    # conservative in-plane shear strength for the same tube.  CALIBRATE
    # against the maker's figures for the tube you buy.
    spar_sigma_allow_pa: float = 600e6
    spar_tau_allow_pa: float = 50e6
    spar_rho: float = 1600.0              # kg/m^3
    spar_min_wall_m: float = 0.0006       # minimum practical layup/tube wall
    spar_fitting_factor: float = 1.8      # webs, joiners, root fittings

    # The ultimate load factor is DERIVED per design from a V-n diagram
    # (see loads.py): the worse of the manoeuvre case and the Pratt gust
    # case, capped by what the wing can aerodynamically generate at the
    # design dive speed, times the factor of safety.  Set the override to
    # a number if tech inspection specifies a wing test instead.
    ultimate_load_factor_override: Optional[float] = None
    # Factor of safety on EVERY wing load case - flight (ultimate = FoS x
    # limit), the wing-tip lift test (rule 4.2.4), shear and torsion.  The
    # spar is sized so that none of them is below it anywhere on the span.
    load_safety_factor: float = 1.50
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
    # The pod that carries a wing motor and its ESC ahead of the wing -
    # firewall, pylon, cowl and connectors - as a share of the motor + ESC
    # mass it holds.  A nose motor mounts on the fuselage's own firewall,
    # so it has none.
    nacelle_mass_frac: float = 0.15
    deploy_mech_kg: float = 0.150         # winch/door for the M3 sensor
    payload_structure_frac: float = 0.080 # restraint structure per kg of payload
    contingency_frac: float = 0.080       # honest margin on the whole airframe

    # drag areas (CD*A, m^2)
    gear_drag_area_m2: float = 0.0040
    # A wing nacelle's drag coefficient on its frontal area - the motor
    # can and a cowl round it (Aircraft.nacelle_width_m).  0.20 is a cowled
    # pod with a cooling inlet: a faired body of revolution is 0.05-0.1, a
    # bare motor can 0.5 or more.  It was a flat 0.0015 m2 per nacelle
    # whatever the motor - a coefficient of 0.44 on a 66 mm pod, and
    # twelve of those on twelve small ones, half the aeroplane's drag.
    nacelle_cd: float = 0.20
    # A measured drag area per nacelle, m2, overrides nacelle_cd.
    nacelle_drag_area_m2: Optional[float] = None
    # Each motor's mount and spinner, for a 62 mm can (a ~350 g motor);
    # scaled with the can's face area, so a 28 mm racing motor's is 5 g.
    motor_mount_kg: float = 0.025
    # The spinner: no thrust inside it, so no slipstream increment behind
    # it.  A propeller's thrust rises from zero there to a peak near 70%
    # of the radius and falls to zero at the tip (the Hough-Ordway
    # distribution, AeroModel._thrust_shape); 0.15 of the radius is a
    # typical spinner on these props.
    prop_hub_frac: float = 0.15
    # Leakage and protuberance drag on the airframe's skin-friction drag:
    # hinge gaps, control horns and linkages, fasteners, hatch and payload
    # door seams, the antenna.  None of it is computable without the built
    # parts.  Junction interference (wing-fuselage, tail-fuselage) is NOT
    # in here any more - it is computed (AeroModel.junction_drag_area_m2).
    # This was a blanket 10% interference on the whole drag, gear and
    # sensor included; it stays 10% on the airframe until a glide test
    # measures the built aeroplane.
    excrescence_frac: float = 0.10

    # Wing on the fuselage.  The fuselage covers part of the planform
    # (fuselage width x root chord), and three things follow from that,
    # all computed in AeroModel:
    #   * span efficiency drops by 1 - 2 (d/b)^2 (Kroo): the lift dips
    #     where the fuselage interrupts the span;
    #   * the covered wing surfaces and the fuselage top under the wing
    #     are not wetted, so they carry no skin friction;
    #   * CLmax: this fraction of the covered planform's lift survives at
    #     the stall.  At small d/b the fuselage carries nearly all of it
    #     (1.0) - right for a top wing, whose upper surface over the
    #     fuselage is still clean, and for a filleted junction.  An
    #     unfaired junction that separates early is nearer 0.5; 0 counts
    #     only the exposed wing, the pessimistic bound.
    fuselage_lift_carryover: float = 1.0
    # Wing surfaces the fuselage covers: 1 for a top (or bottom) wing, 2
    # for a mid wing passing through the fuselage.
    wing_covered_surfaces: int = 1
    # A propeller's slipstream is turbulent, so the wing behind it loses
    # part of the laminar run XFOIL credits it with.  The blade wakes make
    # the transition intermittent rather than immediate, so the washed
    # section's drag is taken this far from the smooth-air value towards
    # fully turbulent (AeroModel._washed_section_cd).  0.5 is the middle of
    # the two bounds - no effect (0) and fully turbulent (1) - until a
    # measurement narrows it.
    slipstream_transition_frac: float = 0.5
    # Tail section thickness - sets the tails' form factor (Raymer) and
    # their fuselage junction drag.
    tail_t_over_c: float = 0.09
    # Fuselage section corners, as a radius over the smaller side.  The
    # section perimeter is 2(w + h) - (8 - 2 pi) r.  0.19 reproduces the
    # 0.92 perimeter factor this replaced.
    fuse_corner_radius_frac: float = 0.19

    # Fuselage layout.  Nothing here is a length or a position: the nose,
    # tailcone, tail arm, wing station and CG are all CALCULATED
    # (layout.py) from the parts, these physical properties, and the
    # target static margin below.
    # Pack density from makers' dimensions (MaxAmps 6S 4500: 165 x 46 x
    # 36 mm, 574 g; Thunder Power 8S 3300: 141 x 53 x 43 mm, 687 g).
    pack_density_kg_m3: float = 2100.0
    # Outrunner can, taking length ~ diameter (Scorpion SII-4025: about
    # 50 mm x 49 mm, 353 g).  Sets the nose's motor section.
    motor_density_kg_m3: float = 3700.0
    # ESC, receiver, Rx battery and fuse harness, packed with wiring.
    avionics_density_kg_m3: float = 1000.0
    motor_mount_length_m: float = 0.025   # firewall + mount behind the can
    # The tailcone tapers to a tail post of this diameter, no steeper than
    # this boat-tail half-angle (beyond ~12-15 deg the flow separates),
    # which sets the shortest tailcone for a given cross-section.
    tail_post_diameter_m: float = 0.030
    tailcone_max_angle_deg: float = 12.0
    tail_boom_fixed_kg: float = 0.030     # tail post fittings
    tail_boom_kg_per_m: float = 0.045     # boom and control links
    # Tricycle main gear sits this far behind the CG, as a tip-back angle
    # from the CG height: the conventional 15 deg.
    gear_tipback_deg: float = 15.0
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
    # The static margin the layout is balanced to - a design requirement,
    # like the factor of safety.  The wing is placed so the CG (on the
    # payload bay) sits this far ahead of the computed neutral point, in
    # mean chords.  10% is a conventional stable, responsive margin.  If
    # the parts cannot balance there (nose-heavy), the real, larger
    # margin is reported instead.
    target_static_margin: float = 0.10
    # How far above the target the layout may leave the margin rather than
    # lengthen the fuselage (or move the battery behind the bay) to hit it
    # exactly: 0.10-0.15 is all comfortably stable and responsive.
    static_margin_band: float = 0.05
    tail_cl_max: float = 1.10        # symmetric section plus elevator
    tail_efficiency: float = 0.85    # tail span efficiency for induced drag
    # Winglets (Design.winglet_h_m, m110): the effective aspect ratio is
    # AR (1 + winglet_ar_gain h/b) - Hoerner's end-plate relation (Fluid-
    # Dynamic Lift, ch. 3).  Each winglet's root chord is the wing's tip
    # chord, its taper 0.5.
    winglet_ar_gain: float = 1.9
    # Tails sized from requirements (m111) - Aircraft._tail_volumes_required.
    # "volume" keeps the fixed vol_coeff_h / vol_coeff_v.
    tail_sizing: str = "requirements"
    tail_trim_cl_frac: float = 0.80      # share of tail CLmax used to trim
    vol_coeff_h_min: float = 0.30        # pitch damping/gusts not modelled
    vol_coeff_v_min: float = 0.020
    cn_beta_target: float = 0.06         # /rad, directional stability
    fin_cl_max: float = 0.90             # fin + rudder side-force CLmax
    fin_motor_out: bool = True           # fin holds one wing motor out
    motor_out_speed_factor: float = 1.20  # x stall, one wing motor out
    winglet_areal_density: float = 1.2   # kg/m^2 of planform
    winglet_fixed_kg: float = 0.010      # each: its fitting
    winglet_t_over_c: float = 0.08
    # A vertical winglet adds no projected span.  True if the organisers
    # count winglet height in the 6 ft limit (rule 3.1 does not say).
    winglets_count_in_span: bool = False

    # --- realism factors: these exist to stop the model flattering itself ---
    # Built wing vs 2-D section data: gaps, finish, washout.  The planform
    # part - where the wing stalls first - is computed separately
    # (liftingline.stall_factor: about 0.92 at taper 0.45, AR 6, since the
    # 2026-10-03 lifting-line fix; 0.93 before).  0.88 x 0.93 was the 0.82
    # blanket 3-D factor this replaced; with the corrected lifting line a
    # taper-0.45 wing gets 0.88 x 0.92 = 0.81, and other tapers are priced
    # by the physics.
    cl_max_3d_factor: float = 0.88
    # Built wing vs smooth-section minimum drag.  The airfoil data is
    # XFOIL, which assumes a perfect surface; film covering, sag between
    # ribs and an imperfect leading edge all add drag it cannot see.  This
    # was a hidden literal 1.05; it is the knob to set once a built wing
    # has been measured (a glide or tow test gives it directly).
    profile_drag_factor: float = 1.05
    # Battery sag.  Each propulsion system's current also flows through its
    # blade fuse, arming plug, connectors and leads - about 3 mOhm for a
    # 100 A blade fuse, two high-current connectors and ~0.3 m of 10 AWG
    # round trip.  Added to the pack's own internal resistance.
    harness_resistance_ohm: float = 0.003
    # Wing-motor power leads.  A nose motor's ESC sits by the battery, so
    # its leads are the harness above.  A wing motor's ESC rides in its
    # nacelle, and its DC leads run from the battery bay out along the
    # wing: the motor's spanwise station plus this much routing (bay to
    # wing root, up into the wing, into the nacelle).  The copper is sized
    # to the motor's rated current at this current density; its resistance
    # and its mass follow from copper's (layout.wing_power_leads), so a
    # lead's loss and its weight trade against each other here.
    lead_route_m: float = 0.20
    lead_current_density_a_mm2: float = 12.0
    # Air reaching the ESCs, m/s: a cooling inlet in the fuselage (or the
    # prop wash over a nacelle) gives them a fraction of the flight speed.
    # The ESC is chosen to hold its motor's peak current, with 20%
    # headroom, for the mission window at this airflow - each ESC's
    # published rating moved to those conditions along Scorpion's derating
    # curve (ESC.sustained_a).  Measure it, or duct the ESCs.
    esc_airflow_mps: float = 10.0
    # Motor heating (propulsion.motor_thermal).  The winding heats with its
    # own losses (I^2 R + iron) and sheds heat to the prop wash; copper
    # resistance rises 0.393 %/C, and a mission whose winding would pass
    # the limit is a failed mission (Mission 3 stops a lap short).
    # 110 C: below N-grade magnet and enamel limits, the usual ceiling for
    # a hobby outrunner flown hard.
    motor_temp_limit_c: float = 110.0
    # Battery heating: each pack heats with I^2 R through its DC resistance
    # (BatteryPack.r_cell_ohm) and, over a five-minute flight inside the
    # fuselage, keeps it.  60 C is the usual LiPo discharge ceiling in the
    # makers' datasheets; a mission whose pack would pass it fails.  This,
    # not the label's C rating (which independent tests put at about half
    # the real sustained figure), is what limits a small pack.
    pack_temp_limit_c: float = 60.0
    pack_specific_heat_j_kgk: float = 1000.0   # LiPo pouch cells
    motor_specific_heat_j_kgk: float = 500.0   # copper / steel / magnet mix
    copper_temp_coeff: float = 0.00393
    # No-load current against speed.  Makers quote I0 at one voltage (10 V;
    # Hacker 8.4 V - Motor.i0_volts); flown on 6S a motor turns about twice
    # as fast.  Its no-load loss is part hysteresis and bearing friction,
    # which draw a fixed current, and part eddy current, whose current
    # grows in proportion to speed (a loss rising with speed squared, over
    # a back-EMF proportional to speed).  This is the eddy share at the
    # maker's test speed:  I0(n) = I0 (1 - s + s n / n_test).  0.5 puts I0
    # at about 1.6x the quoted figure at 6S.  Measure I0 at two voltages
    # to calibrate it; 0 is the old constant-I0 model.
    i0_speed_share: float = 0.5
    # Thermal resistance to the air = k / mass^(2/3) (surface area).  k is
    # the median of the catalogue motors rated CONTINUOUSLY (Scorpion,
    # Neu): at their rated current they reach the limit from 25 C.  The
    # 15-30 s ratings (Hacker, SunnySky) are burst limits, not thermal
    # data, so every motor is held to this common, size-based cooling.
    motor_thermal_k: float = 0.168
    # Depth of discharge: the most of each pack's energy one mission may
    # ever use - takeoff, laps and landing together - the floor nothing
    # passes.  Every mission flies a fresh, fully charged pack, so this is
    # set by what the team will do to a pack: 0.95 is about 3.5 V per cell
    # at rest.  Capped by each chemistry's ceiling
    # (BatteryPack.usable_fraction, 0.95).  Where the laps END is the
    # landing point (land_at_discharge, 88%), so a flight normally lands
    # with about 11% left; this floor only binds if the landing point is
    # raised past it.  The ESC low-voltage cutoff is checked on top.
    depth_of_discharge: float = 0.95
    # Extra energy held back beyond the depth of discharge.  This used to
    # be a hidden 8%, when the pack allowance was also standing in for the
    # sag that was not modelled; with sag modelled it is 0.  Raise it only
    # for a margin you can name.
    energy_reserve: float = 0.0
    # Turn for landing once any pack in any propulsion system has used
    # this share of its energy: no lap is started that would take a pack
    # past it, and the landing is flown out of what is left.  The team's
    # call, 2026-10-01: 88%, so every mission crosses the start/finish line
    # for the last time with at least 12% in every pack - room for the
    # landing circuit and for model and pack error.  depth_of_discharge
    # stays the floor that nothing, landing included, may pass.
    land_at_discharge: float = 0.88
    # Mission 2 used to plan on finishing with 2% of the pack still inside
    # the depth-of-discharge budget, a margin for model error.  The landing
    # point above replaced it on 2026-10-01, so it is 0; it stays a knob
    # for extra Mission 2 margin.
    m2_energy_margin: float = 0.0
    # Every flight planned at full throttle for its whole length - the
    # team's factor of safety, 2026-10-02.  There is no cruise speed for
    # the pilot to hold, so a pilot who flies harder than the plan cannot
    # run a pack past the landing point, and the motors, ESCs and packs
    # are heated by the hardest flight there is.  Every second from
    # throttle-up to the last crossing of the line is charged at
    # full-throttle power and heat: the straights, the turns, and the
    # glide into each turn (flown as a glide, so the lap time is the real
    # one, but charged as if the throttle never came back).  The landing
    # circuit, outside the scored window, keeps its own power.  False
    # brings back the planned cruise: Mission 2 at the fastest speed that
    # still finishes, Mission 3 at the speed that gives the most laps.
    full_throttle_flights: bool = True
    # A freewheeling propeller is a disc of drag once the power is cut.
    # For a propeller with a measured table this is COMPUTED
    # (Propeller.windmill_drag_coeff: the table's last points extended to
    # zero torque); this coefficient on disc area is only the fallback
    # for a propeller without one.
    windmill_drag_coeff: float = 0.08
    # Nobody flies the perfect corner speed on every lap.  Lap time is
    # divided by this.  0.95 (the team's call, 2026-10-02; 0.93 before)
    # means 5.3 percent off a flawless lap; set it from your own flight
    # logs, and set it to 1.0 to model a perfect pilot.
    pilot_efficiency: float = 0.95
    # (The installed-thrust factor that was here, a flat 0.95, is now
    # computed per design from the slipstream scrubbing the airframe
    # behind each propeller: AeroModel.installed_thrust_factor.)

    # operating
    rolling_friction: float = 0.035       # paved runway, pneumatic wheels
    braking_friction: float = 0.30        # wheel brakes / drag on rollout
    takeoff_cl_ground: float = 0.35
    rotate_margin: float = 1.15           # V_rotate / V_stall
    approach_margin: float = 1.30         # V_approach / V_stall
    # The fastest stall the team will fly, m/s: power off, at the Mission 2
    # (heaviest) weight, in the mission's air - checked on the wing's lift
    # divided by stall_lift_safety_factor, the way every wing load carries
    # load_safety_factor.  The aeroplane must hold its weight at this speed
    # on 1/1.5 of its maximum lift, so its true stall is about
    # 1/sqrt(1.5) = 0.82 of the limit.  inf = no limit (--max-stall-mph).
    max_stall_mps: float = float("inf")
    stall_lift_safety_factor: float = 1.50
    # The air the stall limit is checked in, a named condition (m118): the
    # hottest weather preset, so a design sized on the nominal day still
    # meets the limit on a hot afternoon, where the same wing stalls
    # faster in thinner air.  "" = the mission's own air.
    stall_check_conditions: str = "tucson-hot"
    # The most wing angle of attack the aeroplane can reach with its wheels
    # on the runway, deg: its tail-strike pitch plus the wing's incidence -
    # set it from the gear and tail geometry.  A low-aspect-ratio wing makes
    # its lift slowly with angle (at AR 2 it needs about 17 deg to lift
    # off), so this, not its CLmax, can set its liftoff and touchdown
    # speeds and the stall limit (AeroModel.cl_at_ground_alpha, m104).
    max_ground_alpha_deg: float = 15.0
    # The longest wing chord the team will build, m - the root chord, the
    # wing's largest (transport, the build table, the hangar).  inf = no
    # limit (--max-chord-in, m106).
    max_chord_m: float = float("inf")
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


# --------------------------------------------------------------------------
# What each build-standard number IS
# --------------------------------------------------------------------------
#: Every BuildStandard field, sorted by what kind of number it is, so the
#: GUI and the report can say which numbers are safety margins (the
#: team's choice of how much to hold back) and which are estimates
#: waiting for a measurement.  Validation checks every field is here.
#:
#:   margin     deliberate: held back against a failure that costs the
#:              mission, the aeroplane or a rules violation
#:   calibrate  an estimate of a real effect or property the model cannot
#:              compute; replace it with a measurement
#:   design     a choice about how the aeroplane is built or laid out
#:   operating  how it is flown
#:   physics    a physical constant
ROLE_TITLES: Dict[str, str] = {
    "margin": "Safety margins (deliberate - held back on purpose)",
    "calibrate": "Estimates to replace with measurements",
    "design": "Design and build choices",
    "operating": "How it is flown",
    "physics": "Physical constants",
}

FIELD_ROLES: Dict[str, Tuple[str, str]] = {
    # ---- margins ---------------------------------------------------------
    "load_safety_factor": ("margin", "on every wing load case"),
    "load_factor_floor": ("margin", "lowest limit load the spar is sized to"),
    "pilot_overshoot_factor": ("margin", "turn g a pilot overshoots by"),
    "dive_speed_factor": ("margin", "design dive speed / cruise"),
    "gust_cruise_mps": ("margin", "design gust at cruise"),
    "gust_dive_mps": ("margin", "design gust at dive speed"),
    "contingency_frac": ("margin", "unallocated mass on the airframe"),
    "min_static_margin": ("margin", "least static margin that is flyable"),
    "max_static_margin": ("margin", "most before it is all trim drag"),
    "target_static_margin": ("margin", "static margin the layout balances to"),
    "motor_temp_limit_c": ("margin", "winding temperature ceiling"),
    "pack_temp_limit_c": ("margin", "battery temperature ceiling"),
    "depth_of_discharge": ("margin", "share of each pack a mission may use"),
    "energy_reserve": ("margin", "extra energy held back, every mission"),
    "land_at_discharge": ("margin", "turn for landing at this share used"),
    "m2_energy_margin": ("margin", "Mission 2 energy left for model error"),
    "full_throttle_flights": ("margin", "every flight charged at full "
                                        "throttle, end to end"),
    "rotate_margin": ("margin", "rotation speed / stall speed"),
    "approach_margin": ("margin", "approach speed / stall speed"),
    "max_stall_mps": ("margin", "fastest stall flown, on the factored lift"),
    "stall_lift_safety_factor": ("margin", "lift divided by this for the "
                                           "stall limit"),
    "stall_check_conditions": ("margin", "air the stall limit is checked in "
                                         "('' = the mission's)"),
    "max_ground_alpha_deg": ("design", "most wing angle of attack on the "
                                       "runway: tail strike + incidence"),
    "max_chord_m": ("design", "longest wing chord built (the root)"),
    "prop_clearance_margin_m": ("margin", "propeller tip to the ground"),
    # ---- estimates to calibrate -------------------------------------------
    "wing_areal_density": ("calibrate", "weigh built wing panels"),
    "h_tail_areal_density": ("calibrate", "weigh a built stabiliser"),
    "v_tail_areal_density": ("calibrate", "weigh a built fin"),
    "fuse_areal_density": ("calibrate", "weigh a built fuselage"),
    "spar_sigma_allow_pa": ("calibrate", "tube maker's compression strength"),
    "spar_tau_allow_pa": ("calibrate", "tube maker's shear strength"),
    "spar_rho": ("calibrate", "tube density"),
    "spar_fitting_factor": ("calibrate", "joiners, webs, root fittings"),
    "gear_mass_frac": ("calibrate", "weigh the gear"),
    "gear_mass_fixed_kg": ("calibrate", "weigh the gear"),
    "servo_mass_kg": ("calibrate", "servo data sheet"),
    "avionics_kg": ("calibrate", "weigh receiver and wiring"),
    "rx_battery_mass_kg": ("calibrate", "weigh the Rx pack"),
    "fuse_harness_mass_kg": ("calibrate", "weigh the fuse harness"),
    "arming_plug_mass_kg": ("calibrate", "weigh the arming plug"),
    "nacelle_mass_frac": ("calibrate", "weigh a built wing nacelle"),
    "deploy_mech_kg": ("calibrate", "weigh the deploy mechanism"),
    "payload_structure_frac": ("calibrate", "weigh the payload restraint"),
    "gear_drag_area_m2": ("calibrate", "tow or glide test"),
    "nacelle_drag_area_m2": ("calibrate", "tow or glide test"),
    "nacelle_cd": ("calibrate", "tow or glide test, on frontal area"),
    "motor_mount_kg": ("calibrate", "weigh a mount and spinner"),
    "prop_hub_frac": ("design", "spinner radius / propeller radius"),
    "excrescence_frac": ("calibrate", "hinges, horns, seams: glide test"),
    "fuselage_lift_carryover": ("calibrate", "lift over the fuselage at stall"),
    "slipstream_transition_frac": ("calibrate",
                                   "laminar run lost behind the prop, 0-1"),
    "pack_density_kg_m3": ("calibrate", "pack mass / box volume"),
    "motor_density_kg_m3": ("calibrate", "motor mass / can volume"),
    "avionics_density_kg_m3": ("calibrate", "measure the avionics bay"),
    "tail_boom_fixed_kg": ("calibrate", "weigh the tail post fittings"),
    "tail_boom_kg_per_m": ("calibrate", "weigh boom and links per metre"),
    "tail_eta": ("calibrate", "tail dynamic pressure ratio"),
    "winglet_ar_gain": ("calibrate", "effective AR gain per h/b: tow or glide test"),
    "tail_sizing": ("design", "tails from requirements or fixed volumes"),
    "tail_trim_cl_frac": ("margin", "share of tail CLmax used to trim"),
    "vol_coeff_h_min": ("margin", "least horizontal tail volume"),
    "vol_coeff_v_min": ("margin", "least vertical tail volume"),
    "cn_beta_target": ("margin", "directional stability the fin must give"),
    "fin_cl_max": ("calibrate", "fin + rudder side-force coefficient"),
    "fin_motor_out": ("design", "fin sized to hold one wing motor out"),
    "motor_out_speed_factor": ("margin", "x stall where one wing motor out is held"),
    "winglet_areal_density": ("calibrate", "weigh a built winglet"),
    "winglet_fixed_kg": ("calibrate", "weigh a winglet fitting"),
    "winglet_t_over_c": ("design", "winglet section thickness"),
    "winglets_count_in_span": ("design", "rule 3.1 reading: winglets in the span"),
    "fuselage_np_shift_k": ("calibrate", "fuselage pitch term: AVL or flight"),
    "tail_cl_max": ("calibrate", "tail section + elevator"),
    "tail_efficiency": ("calibrate", "tail span efficiency"),
    "cl_max_3d_factor": ("calibrate", "built wing vs section: stall test"),
    "profile_drag_factor": ("calibrate", "built wing vs XFOIL: glide test"),
    "harness_resistance_ohm": ("calibrate", "measure fuse + plug + leads"),
    "i0_speed_share": ("calibrate", "no-load current at two voltages"),
    "esc_airflow_mps": ("design", "cooling air reaching the ESCs"),
    "motor_specific_heat_j_kgk": ("calibrate", "motor heat capacity"),
    "pack_specific_heat_j_kgk": ("calibrate", "pack heat capacity"),
    "motor_thermal_k": ("calibrate", "winding temperature after a run"),
    "windmill_drag_coeff": ("calibrate", "fallback when no prop table"),
    "pilot_efficiency": ("calibrate", "your flight logs"),
    "rolling_friction": ("calibrate", "runway surface and wheels"),
    "braking_friction": ("calibrate", "brakes and surface"),
    "takeoff_cl_ground": ("calibrate", "ground-roll attitude"),
    # ---- design and build choices -----------------------------------------
    "spar_min_wall_m": ("design", "thinnest tube wall you can buy"),
    "n_servos_base": ("design", "servo count"),
    "wing_covered_surfaces": ("design", "1 top wing, 2 mid wing"),
    "tail_t_over_c": ("design", "tail section thickness"),
    "fuse_corner_radius_frac": ("design", "fuselage corner radius / side"),
    "motor_mount_length_m": ("design", "firewall and mount"),
    "lead_route_m": ("design", "battery bay to a wing nacelle, extra run"),
    "lead_current_density_a_mm2": ("design", "wing lead copper per amp"),
    "tail_post_diameter_m": ("design", "tailcone end"),
    "tailcone_max_angle_deg": ("design", "steepest non-separating taper"),
    "gear_tipback_deg": ("design", "main gear behind the CG"),
    "vol_coeff_h": ("design", "horizontal tail volume"),
    "vol_coeff_v": ("design", "vertical tail volume"),
    "tail_aspect_ratio": ("design", "tail aspect ratio"),
    "static_margin_band": ("design", "tolerance above the target"),
    "gear_height_m": ("design", "main gear leg length"),
    # ---- operating --------------------------------------------------------
    "descent_rate_mps": ("operating", "approach descent rate"),
    "flare_time_s": ("operating", "flare duration"),
    "max_bank_deg": ("operating", "steepest bank flown"),
    "max_turn_load_factor": ("operating", "fallback turn load factor"),
    # ---- physics ----------------------------------------------------------
    "copper_temp_coeff": ("physics", "copper resistance per C"),
    # not numeric: the GUI does not show it, the report does when it is set
    "ultimate_load_factor_override": ("margin", "tech-inspection wing test"),
}


def field_role(name: str) -> Tuple[str, str]:
    """(role, one-line description) of a BuildStandard field."""
    return FIELD_ROLES.get(name, ("design", ""))
