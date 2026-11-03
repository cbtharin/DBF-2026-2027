# Assumptions

Every number this model produces rests on one of four kinds of assumption.
They are not equally trustworthy, and they are tagged accordingly:

| Tag | Meaning | How to fix it |
|---|---|---|
| **RULE** | A competition rule | Confirmed by the team unless marked otherwise |
| **CALIBRATE** | A shop-specific number | Weigh or measure your own hardware |
| **MODEL** | A physics simplification | Understand the limit; upgrade if it binds |
| **NUMERICAL** | Discretisation / solver error | Measured; see the bottom section |

Honest summary: **numerical error is negligible** (parts in 10^7 on lap
time). **Model error is a few percent.** The one open rule is what the sensor
is allowed to be: with no cap on its mass, the optimiser makes it as heavy as
the 100 Wh and the 300 s window allow (the 20 kg search bound is set wide
enough not to bind), because sensor weight scores in both missions. See
section 8.

---

## 1. RULE — supplied by the team for 2026-27

| Assumption | Value | Status |
|---|---|---|
| Fly-off site | TIMPA, Tucson AZ, 15-18 April 2027 (rules 2.4) | given |
| Air density | 1.082 kg/m3 - 667 m, 28 C, ~4,170 ft density altitude | climate normals, NOT a forecast |
| Site elevation | 2,188 ft (667 m) - USGS 3DEP at 32.2647N 111.2740W | measured |
| Temperature | 28 C nominal (April normal high 83 F); 12 C morning, 36 C hot | NWS normals 1992-2021 |
| Wind | 5.0 m/s nominal; April is Tucson's windiest month | climate normals |
| Humidity | not modelled; dry air is denser, so this is conservative by ~0.3% | our choice |
| Lap course | 2 × 1000 ft straights, 180° turn at each end, 360° loop at the midpoint of the second straight | given |
| Total straights / turning | 2000 ft / 720° | derived |
| Takeoff distance limit | **none** | given |
| Battery energy cap | 100 Wh total, **and** 100 Wh per individual pack | given |
| Battery chemistry | LiPo, Li-ion, NiMH, NiCd all legal; **only LiPo and Li-ion searched** | given / our choice |
| Motors | COTS, unmodified, brushed or brushless; propellers commercially produced | given |
| Battery packs | unaltered, commercially procured **COTS** only; no custom packs; labels must state chemistry, mAh, V and C rating | given |
| Pack wiring | one pack per propulsion system; packs may **not** be wired in series or parallel; identical packs if more than one system | given |
| Fuse | blade fuse per system, rated at the pack's C x Ah, **100 A maximum**; each system's peak current (full charge, full throttle, worst airspeed, windings at the air temperature, every weather day's propeller) must stay inside it | given |
| Total current | **no rules limit**: each system's own fuse is the limit (a 100 A total was a team limit and has been removed) | given |
| Rx/servo battery | separate pack required, any chemistry, outside the 100 Wh budget; ESC BEC disabled | given |
| Maximum motors | **no number in the rules**; 4 is the team's cap (`Rules.max_motors`), raised by `--motor-counts` to the largest count asked for | given / team cap |
| Maximum gross weight | **55 lb** | given |
| Mission window | 300 s, M2 and M3 | given |
| Ground Mission | six drops on concrete, declared height up to 60 in; GM = 0.5 + (weight x height)/Max | given |
| Total Mission Score | **M1 + M2 + M3 + GM**, ceiling 7.5 | given |
| Simulator weight | **fixed** at sensor + container tare, +/- 1 oz (3.1.3) | given |
| Mission 3 sensor | may be lighter than the declared max, never heavier (3.1.1) | given |
| Sensor | towed, min 6 in long, own battery and lights, carried internally | given |
| Tow line | at least 1.5 x wingspan when deployed (3.1.2) | given |
| Course | markers 500 ft each side of the start/finish line | given |
| Wing tip load test | at MGTOW, roughly a 2.5 g case, at Tech Inspection | given |
| M2 required laps | 5 | given |
| Max simulators | **no rules cap** — drag and weight limit it | given |
| Max sensor / simulator mass | **no rules cap of any kind** | given |
| Shipping container tare | 0.5 lb floor | given |
| Landing | required to score; **not** in the scored time window | from the rules text |
| Maximum wingspan | **6 ft (72 in)**, the wing's own span tip to tip: a wingtip propeller's disc may overhang the tip (`Rules.wingspan_counts_props`, off) | given, read literally |
| Ultimate load factor | **derived per design** from a V-n diagram | see section 2a |

**Turn distance.** Turn *time* is `angle / turn_rate` with
`turn_rate = g·√(n²−1)/V`, so it depends on how hard each aeroplane can turn —
exactly as you said it should. The ground-track length of a turn never enters
the time calculation, so nothing has to be assumed about turn radius.

**Not modelled:** Mission 1's exact requirements (treated as 3 laps,
pass/fail), the design report score, the total competition score, and tech
inspection beyond the payload count.

---

## 2a. DERIVED - design load factors (no longer assumed)

The spar used to be sized to a hand-picked 4.5 g.  It is now sized from a
V-n diagram built for each aeroplane, and the limit load factor is whichever
of three cases governs:

1. **Manoeuvre** - the sustained turn the mission demands, `n = 1/cos(bank)`
   at the 70 deg bank limit, times a 1.15 pilot overshoot factor -> 3.36 g.
2. **Gust** - the Pratt discrete-gust formula, at cruise with a 25 ft/s
   rough-air gust and at the design dive speed (1.25 Vc) with half of it:

   ```
   mu = 2(W/S)/(rho c a g)      Kg = 0.88 mu/(5.3+mu)
   n  = 1 + Kg rho V a U / (2 W/S)
   ```

3. **Aerodynamic ceiling** - you cannot pull more g than the wing can make.
   At the dive speed that is `n = q_D S CLmax / W`, and the limit is capped
   there, because designing for a load the wing physically cannot generate
   is wasted weight.

`ultimate = 1.5 x limit`, the factor of safety. The same 1.5 is applied
to every wing load case: flight bending, the wing-tip lift test at MGTOW
(rule 4.2.4), and shear plus torsion, each checked at 24 stations along
the span (see section 5).

**This changed the answer.** Light wing loadings are gust-critical, so a
lightly loaded DBF wing needs *more* than 4.5 g, not less:

| design | governed by | ultimate |
|---|---|---|
| 1.6 m span, AR 6, light | gust at cruise | 5.98 g |
| 1.8 m span, AR 12 | gust at cruise | 5.54 g |
| 1.0 m span, AR 6 | manoeuvre | 5.04 g |
| 1.8 m span, heavy payload | wing CLmax ceiling | 4.74 g |

The heavily loaded aeroplane needs the *least*, because it cannot generate
much g before stalling.  The old flat 4.5 g was too low for the light cases
and too high for the heavy ones.

Normal flight is held inside **limit** load, not ultimate, so the mission's
turn model uses the derived limit rather than a fixed 3.5 g cap.

**Caveat worth knowing:** on this class of aeroplane the spar is often
**gauge limited** - the strength requirement falls below the thinnest wall
you can actually lay up or buy, so the load factor stops mattering.  The
report says when that is the case.  Set
`BuildStandard.ultimate_load_factor_override` if tech inspection specifies
a wing test instead.

---

## 2. MODEL — payload and container geometry

The container is **one fixed assumed box**, the same for every design:
**230 x 150 x 150 mm at 0.5 lb tare.**

Those dimensions come from the drop-test calculation for a roughly 1 kg
sensor — a 5 ft drop reaches 5.47 m/s, holding the sensor under 75 g needs a
20 mm crush stroke, and foam crushes usefully through about 70% of its
thickness, so 29 mm of foam per face around a 165 x 83 x 83 mm sensor —
rounded to a buildable box. It is **not** recomputed per design. Change the
four numbers in `Rules.container_*` if the real container differs.

- **Simulators share that external envelope**, because that is what they
  simulate, and carry whatever mass you load them to.
- The sensor's own envelope is estimated at a nominal 800 kg/m3, but that
  is NOT a cap. There is no sensor mass ceiling: a heavier sensor in a fixed
  box is simply a denser one, and the report states the density it would
  need (water 1000, aluminium 2700, steel 7850, lead 11340 kg/m3) so you can
  judge whether it is buildable. The only rejection is a density greater
  than lead, which is physically impossible rather than merely awkward.

So each extra simulator adds a container envelope of bay volume, which
lengthens or widens the fuselage, which adds wetted area and drag. That is
the coupling you asked to be sure was there.

**Bay packing is chosen, not assumed.** The model tries 1–4 boxes across and
1–3 layers high and picks the arrangement with the least wetted area, with a
penalty (×1.6) for a bay wider or taller than 0.55 of its length.

### Fuselage length and part placement are calculated

Nothing about where the parts sit, or how long the fuselage is, is a fixed
proportion (`dbf2027/layout.py`; the report's LAYOUT section lists every
station).

- **Payload bay on the CG**, so the empty, loaded and sensor-towing
  aeroplanes balance alike.
- **Wing placed for the target static margin**: the CG 0.10 of the mean
  aerodynamic chord (MAC, the reference chord for pitch) ahead of the
  computed neutral point, with up to 0.05 more accepted rather than a
  longer fuselage. It is a design requirement, like the factor of safety.
  The neutral point counts the tail's own share of the lift slope,
  `x_np/c = 0.25 + [η V_h (a_t/a_w)(1 − dε/dα) − fus] / (1 + η (S_h/S)(a_t/a_w)(1 − dε/dα))`
  on the MAC, fus the fuselage's destabilising term. Without the divisor
  (the usual textbook approximation, used until 2026-10-03) it sat about
  0.02 MAC too far aft, so a layout balanced to the 10% target really had
  about 8%.
- **Battery as the trim weight.** It sits between the motor and the bay,
  or behind the bay when a long bay makes the aeroplane nose-heavy. The
  nose is lengthened only if the battery cannot balance the tail.
- **Part sizes from their masses.** Battery length comes from pack volume
  at 2.1 g/cm³, which is what the makers' dimensions give (MaxAmps 6S 4500:
  165 x 46 x 36 mm, 574 g). Motor length comes from outrunner density,
  about 3.7 g/cm³ (Scorpion SII-4025: about 50 x 49 mm, 353 g).
- **Tail arm**: the lightest tails-plus-tailcone,
  `l* = sqrt((ρ_h V_h S c + ρ_v V_v S b)/λ)` (c the MAC), measured from
  the wing's aerodynamic centre to the tail's. It is no shorter than a
  tailcone tapering from the bay to a 30 mm tail post at 12° (steeper and
  the flow separates).
- **Tails are sized from their requirements (m111), and the horizontal
  one may be bigger by choice (m112).** The tail volumes are no longer
  fixed coefficients. They were V_h 0.55 and V_v 0.040, and
  `tail_sizing = "volume"` still uses those. The requirements are
  re-solved inside the layout loop, on the current CGs, neutral point,
  weights and propulsion.
  - **Horizontal tail: the largest of three requirements.**
    - **Trim across the flight envelope,** at 0.8 of the tail's CLmax
      (`tail_trim_cl_frac`). Each loading (Mission 2, Mission 3 towing,
      Mission 1 empty) is checked at its own CG, at two points:
      - the stall at 1 g. In the linear trim model the manoeuvre corner,
        n g at √n × the stall, asks the tail for the same CL;
      - the dive speed, where the camber moment's download is largest.
    - **Rotation.** At the liftoff speed the elevator, at the tail's full
      CLmax, must pitch the aeroplane about its main wheels:
      `D_t l_mt ≥ W d_mg + T (h_T − h_cg) + μ N h_cg − L_w d_ac − M_ac`.
      Every arm comes from the layout's own stations. The thrust line is
      counted from the CG because the roll's acceleration acts there.
    - **A floor,** V_h 0.30 (`vol_coeff_h_min`), since pitch damping and
      gusts are not modelled.
  - **The horizontal tail volume is also a search variable** (`tail_vol_h`,
    0.30 to 0.90, `--tail-vh`), because the smallest tail that meets the
    requirements is not the best one.
    - **Why bigger can pay:** a bigger tail moves the neutral point aft.
      The CG follows it, and the camber moment is trimmed by the wing's
      own lift instead of a tail download.
    - **Example:** on the v1.16 heavy lifter (DAE21, 11.2 kg sensor),
      going from the requirement's V_h 0.30 to 0.80 cut the cruise
      download from 20 N to 1 N and Mission 2 from 215 s to 198 s. It
      cost 0.32 kg more tail and tailcone.
    - **The search** finds where that trade stops paying.
    - **Saved designs:** the designs in `examples/` were optimised with
      the old fixed tail, so they record `tail_vol_h` 0.55.
  - **Fin: the largest of three requirements.**
    - **Directional stability:** Cn_β ≥ 0.06 /rad (`cn_beta_target`)
      against the fuselage's destabilising `−1.3 Vol/(S b) (D/W)` (Raymer
      eq. 16.47). The fin's slope is taken at 1.55 × its geometric aspect
      ratio (end-plating, Raymer).
    - **One wing motor out,** when there are two or more wing motors
      (`fin_motor_out`). At 1.2 × the stall (`motor_out_speed_factor`)
      the fin must hold the outermost motor's thrust plus the dead
      propeller's windmill drag, at a fin + rudder CL of 0.9
      (`fin_cl_max`).
    - **A floor,** V_v 0.020.
    - **A deep fuselage needs a big fin.** The four-simulator baseline
      stacks its bay two high (345 × 185 mm) and needs about three times
      the old fixed fin.
    - **The fin stays at its requirement:** it does not trim, so a bigger
      one only weighs and drags more.
  - **The tail check** uses the same envelope at the tail's full CLmax. It
    used to ask for the limit load at 1.3 × the stall, a lift the wing
    cannot make there: 1.3 × the stall holds 1.69 g at most.
  - **The thrust line is trimmed in flight (m116).** Wing propellers on a
    top wing sit half a fuselage depth above the CG, so full throttle
    pitches the nose down and the tail pushes down more:
    `L_t = [c (Cm_ac q S + L k) − T dz] / l_t`, with the same term in the
    lift at the stall.
    - **The height:** `dz` is the thrust line over the CG, weighted over
      the motors. The nose propeller is on the centreline and adds none; a
      mid wing's propellers add none.
    - **The CG's height is assumed,** at mid-fuselage, as in the rotation
      case. A top wing with its motors raises it, so `dz` is an upper
      bound.
    - **What it costs a two-motor design:** about 1% of the weight more
      download at full throttle, and about 1% more cruise drag. On the v2
      40–50 mph winners that was enough to lose a Mission 3 lap they held
      by seconds (0.17–0.22 of score before re-optimising).
  - **Not modelled:** crosswind landing (rudder authority), the
    propeller's own destabilising side force, and dynamic modes.
  - **CALIBRATE** `fin_cl_max`. The requirements' margins are design
    choices.
- **Fuselage = nose + bay + tailcone.** The wetted area includes the taper.
- **Gear**: tricycle, mains 15° tip-back behind the CG.

Still assumed, and small: the motor mount length (25 mm), the tail-post
diameter, the electronics' packing density, typical pack proportions
(4 : 1.3 : 1), and the structural centroids of the wing (0.42 c skin,
0.28 c spar) and tails.

---

## 3. MODEL — aerodynamics, with every airfoil judged fairly

Component parasite-drag buildup + section polar + lifting-line induced drag
+ trim. Two effects exist specifically so high-lift sections are judged on
what they deliver rather than on their catalogue numbers:

1. **Reynolds-dependent CLmax.** Example (examples/optimised_v1.json): the
   wing runs at Re 413k–654k on its mean chord, from its Mission 2 stall to
   its top speed (250k at the tip at stall), against section data quoted at
   500k, and high-camber sections lose more at low Re. The penalty exponent scales with the section's own camber, so S1223 is
   docked harder than MH32. Quoting S1223 at CLmax 2.2 on a 9 inch chord is
   exactly the optimism this is meant to prevent.
2. **Trim drag.** A cambered section has a nose-down moment the tail must
   balance with a download, so the wing carries *more* than the aircraft
   weight and the tail makes induced drag of its own. S1223 has Cm = −0.27
   against MH32's −0.06. Without this, high-lift sections win on lift and pay
   nothing. The moments are taken about the CG with the tail arm measured
   from the wing's aerodynamic centre to the tail's, as the layout measures
   it, and the MAC as the pitch reference: `L_t = c (Cm_ac q S + L k) / l_t`,
   k the CG's distance aft of the wing's AC in MACs. (Until 2026-10-03 the
   balance was written for an arm measured from the CG, on the mean chord.)

Other choices, all deliberately sober:

- Section data: **all 32 airfoils from airfoiltools.com** — the XFOIL
  polar at Re 500,000, Ncrit 9. The model reads the (cl, cd) table itself;
  a fitted parabola was off by more than 50% for the high-camber sections.
  Airfoils added by hand fall back to the fitted parameters. **XFOIL is a
  smooth-surface prediction** — its minimum drag is 22–53% below the
  previous catalogue's values.
- 3-D CLmax = 2-D section CLmax × **0.88** (built wing: hinge gaps,
  washout, surface finish) × the **planform stall factor** × the
  **fuselage factor**. The stall factor is computed: the lifting-line
  loading against each section's own Reynolds number, so the wing stalls
  where its first section does. At AR 6.5 it is about 0.88 for a
  rectangular wing, 0.92 at taper 0.45 and 0.78–0.81 at taper 0.20, where
  the small tips stall first. (0.88 was chosen when the factor came out
  0.93 at taper 0.45, so that 0.88 × 0.93 kept the old blanket 0.82; since
  the 2026-10-03 lifting-line fix it gives 0.81 there.) Example
  (examples/optimised_v1.json): 0.913 at taper 0.43 and AR 5.1, stalling
  at 66% of the semi-span. That tip-stall penalty is what lets taper run
  down to 0.20. It is the only finite-span factor: a second one,
  `min(1, 0.80 + 0.04 AR)`, counted the same effect twice below AR 5 and has been removed. The
  fuselage factor is `1 − (1 − fuselage_lift_carryover) × S_cov/S`. At
  the default carry-over of 1.0 (right for a top wing), it is 1.
- **The wing is rectangular for now (m101, 2026-10-08).** The search holds
  taper at 1.0 (`DesignSpace.taper_bounds = (1, 1)`), and Nelder–Mead
  leaves it out of the simplex. The model still prices any taper: the
  lifting line, and the tip stall above. `--taper 0.2 1.0`, or the GUI's
  taper fields, search it again. Seeds saved with a taper fly
  rectangular.
- **Aspect ratio 1.5 to 14 (m104; it was 3.5 to 14).** The bottom is
  a 1.22 m chord on the 6 ft span. `--ar` changes the range; the model
  itself goes down to 1.0 (m106), but below about 2 the lifting line is
  only approximate. `--max-chord-in` caps the root chord instead
  (`max_chord_m`), which floors the aspect ratio at span ÷ chord.
  - **The catch:** a stubby wing makes its lift slowly with angle. The
    finite-wing slope is 2.6 /rad at AR 2 against 4.5 at AR 6 (for a 2π
    section), so lifting off at 1.15 × its stall would take about 17° at
    AR 2 and 22° at AR 1.5.
  - **The runway attitude limit:** `max_ground_alpha_deg` (15°: the
    tail-strike pitch plus the wing incidence, `--max-ground-alpha`) caps
    the wing's angle of attack with the wheels down. The lift the wing
    makes there is the 3-D slope × that angle, measured from zero lift.
  - **That lift sets three things:** the liftoff speed (the larger of
    1.15 × the power-on stall and where the attitude holds the weight),
    the approach (1.3 × the larger of the stall and that speed), and the
    stall limit.
  - **The zero-lift angle and the section lift slope come from the
    polars (m109).** Each is a straight-line fit of Cl on α over −2° to
    +6° on the same airfoiltools XFOIL polars (Re 500,000, Ncrit 9) the
    drag tables come from, with R² 0.990–0.9999
    (results/experiment_2026-10-02/data/airfoil_lift_xfoil.json).
    - **Why:** the angle had been estimated from the camber moment,
      4 cm_ac/π (thin-airfoil theory). That is far off for aft-loaded
      high-lift sections: −19.5° for the S1223 against −10.0° from its
      polar, and −9.8° for the DAE21 against −6.0°.
    - **What it changed:** the old estimate overstated the lift a
      high-camber wing makes at the runway attitude, and with it the
      stall-limited results.
    - **Hand-entered sections:** with the catalogue columns left blank,
      the model falls back on the estimate and 2π.
  - **The wing's 3-D slope is DATCOM's (Polhamus) form:**
    `2π AR / (2 + √((AR/k)² + 4))`, with k = a0/2π the section's slope
    over 2π. It tends to the section's own slope at long aspect ratio and
    to πAR/2 at short. The Helmbold form used before,
    `a0 AR / (2 + √(AR² + 4))`, gives the same at a0 = 2π but is wrong at
    long aspect ratio for any other a0. The tails keep 2π.
  - **Normal wings are unaffected.** At AR 4–6 the attitude isn't
    limiting, and nothing changes.
  - **Not modelled:** low-aspect-ratio vortex lift (below about AR 2), so
    CLmax there is conservative.
- **Root chord = 2c/(1 + λ)**, c the mean chord S/b, so the chords that
  get built have the area that gets simulated. (The mean-aerodynamic-chord
  formula used before drew 5% less wing at taper 0.43, 13% at 0.20.)
- Wing profile drag is the section `cd_min`, Re-corrected, ×
  `profile_drag_factor` (**1.05**, the built-wing allowance over a smooth
  XFOIL section). **CALIBRATE** — a glide or tow test on a built wing
  gives it directly.
  The wing surface that sits on the fuselage is not wetted: one surface
  of the covered strip for a top wing (`wing_covered_surfaces` 1), both
  for a mid wing.
- Fuselage: flat-plate `Cf` with a fineness form factor, 10% laminar run.
  Rounded-corner section (radius 0.19 of the smaller side, replacing a
  flat 0.92). The top the wing sits on is not wetted, but it is weighed.
- Tails: Raymer's thickness form factor at t/c 0.09 (1.19, was a flat
  1.25).
- **Junction interference is computed**: Hoerner, `ΔD/(q t²) = 0.75 t/c
  − 0.0003/(t/c)²` per root. Two wing roots (half each for a top wing),
  two stabiliser roots, one fin. It was part of a flat +10% on all drag.
- **Excrescences +10%** (`excrescence_frac`), on the skin-friction terms
  only: hinges, horns, linkages, seams, the antenna. **CALIBRATE** with a
  glide test.
- Discrete drag areas (CD·A): gear 0.0040 m² (more, in proportion, for
  legs taller than the standard 130 mm), deployed sensor 0.0045 m² plus its
  tow line (Cd 1.1 on a 1.5 mm line 1.5 spans long), and each wing
  nacelle: `nacelle_cd` **0.20** on its frontal area, a cowled pod 1.25 ×
  the motor's can diameter, the can's diameter from its mass
  (`layout.motor_diameter_m`, a power law through KDE Direct's sizes). A
  measured `nacelle_drag_area_m2` overrides it. (It was a flat 0.0015 m²
  per nacelle whatever the motor.) Example
  (results/exp8_C_cap2_2026-10-03/best_design.json, a 275 g Neu 6510M):
  0.0008 m² per nacelle. **CALIBRATE.**
- Span efficiency from the lifting line, **× (1 − 2(d/b)²) for the
  fuselage** (Kroo). Example (examples/optimised_v1.json): 0.992 × 0.980 =
  0.972. The lifting line's collocation rows are `[2 n μ + sin θ]`, as its
  `Γ = 4 b V Σ Aₙ sin nθ` normalisation requires; until 2026-10-03 they
  were `[n μ + sin θ]`, which solved a wing of twice the aspect ratio (a
  rectangular AR 6 came out e = 0.907 against the textbook 0.954). The span
  efficiency, the stall factor and the blown-wing ratios all come from
  that solve.
- **Winglets are a search variable (m110).** Each is a vertical
  end-plate at the tip; its height is a share of the span (`--winglets`,
  0 to 0.15 by default; 0 is none).
  - **The benefit:** Hoerner's end-plate relation gives the effective
    aspect ratio, `AR_eff = AR (1 + 1.9 h/b)` (`winglet_ar_gain`; it
    holds to h/b of about 0.2). AR_eff sets the induced drag,
    `k = 1/(π AR_eff e)`, the wing's lift slope and the downwash at the
    tail. CLmax is unchanged.
  - **The cost:** what they weigh (`winglet_areal_density` 1.2 kg/m² on
    h × 0.75 of the tip chord each, plus 10 g of fixing), their skin drag
    at t/c 0.08 (a "winglets" item in the drag breakdown) and a junction
    each.
  - **The span rule:** they do not count against the span limit unless
    `winglets_count_in_span` is set. Check the rules' wording.
  - **The spar carries their tip loading (m117).** An end plate loads the
    tips, which is how it cuts the induced drag.
    - **The lift distribution:** the flight bending and shear take the
      inboard part of an ellipse over the equivalent span,
      `b √(AR_eff / AR)`, cut off at the real tips and carrying the whole
      lift. For h/b = 0.15 that is 6.6% more root bending, and more
      outboard.
    - **Where it matters:** most of this class's spars are at minimum
      gauge, which covers it; only a spar sized by its load gets heavier.
    - **Not sized:** the winglets' own side load on the tip joint.
  - **CALIBRATE** `winglet_ar_gain` with a tow or glide test.

**Net effect of the computed corrections**, when they replaced the flat
factors (on the design then recommended): CD0 at 20 m/s fell from 0.0289
to 0.0271 (−6%). The biggest changes were the gear and sensor no longer
carrying 10% interference, and the skin the wing and fuselage cover.
Installed thrust rose from the flat 0.95 to 0.963, and induced drag rose
2%. The old flat factors were partly acting as a hidden margin; that
margin now has to come from measurement. `excrescence_frac` and
`profile_drag_factor` are the numbers a glide test sets. Example with
today's code (examples/optimised_v1.json): CD0 0.0282 at 20 m/s and
0.0266 at its 31.4 m/s top speed.

**The propeller slipstream is modelled five ways.** The drag with the
propellers making thrust (power-on drag) is used wherever the motors run:
the climb, the turns, the powered straights and both top-speed solves; the
takeoff roll and its rotation speed carry the slipstream's lift. Only the
glide into a turn and the landing use the power-off drag.

- **Scrubbing:** the airframe inside the slipstream drags more, so
  installed thrust is `T (1 − f/A)`, worked out at every airspeed up to
  the propellers' speed ceiling (section 4).
- **Turbulence:** the wing strip behind the propeller loses half its
  laminar-run advantage (`slipstream_transition_frac` 0.5, between none
  and fully turbulent), and behind a nose propeller the fuselage loses
  half its 10% laminar run. **CALIBRATE** if you can measure it. Example
  (examples/optimised_v1.json): the wing strip's share is 2.6% of CD0 at
  20 m/s.
- **Lift:** that strip lifts more with the propellers making thrust. This
  comes from the lifting line with the jet's airspeed in it. The jet's
  speed follows the propeller's radial loading (the Hough–Ordway shape:
  none inside a spinner of 0.15 of the radius, `prop_hub_frac`, a peak
  near 70% of the radius, none at the tip), and the body behind the
  propeller - the fuselage behind a nose propeller, the nacelle behind a
  wing one - takes the middle of the slipstream, so only the wing either
  side of it is blown. The stall search skips the wing over the fuselage.
  Example (examples/optimised_v1.json): ×1.04 lift and ×1.03 CLmax at
  28 m/s, full throttle, 70% charge, and a Mission 2 stall of 19.0 m/s at
  full throttle against 19.8 m/s power off.
- **Swirl recovery:** wing propellers turn inboard-up, each pair
  counter-rotating, and the wing behind them takes back part of their
  swirl as a cut in induced drag (the swirl angle from the propeller's
  torque by momentum theory, put through the lifting line). Only the
  swirl's share is credited: an axial slipstream by itself does not lower
  the span's induced drag. Example
  (results/exp8_C_cap2_2026-10-03/best_design.json, two wingtip 16x12s):
  induced drag ×0.88 at 25 m/s full throttle, ×0.93 at its 33.5 m/s top
  speed.
- **Roll trim:** an unpaired nose motor's torque is held by the ailerons,
  whose antisymmetric lift costs induced drag `32 Cl²/(π AR)` on `q S`,
  `Cl = Q/(q S b)`. Wing pairs counter-rotate and cost none.

**Known limits:** a nose propeller's swirl (antisymmetric across the span,
it cancels to first order in the symmetric lifting line), swirl's effect
on where the wing stalls, its sideways push on the fin, and the tails in
the slipstream are not modelled; no flutter or aeroelastic divergence,
which matters because nothing else in the model punishes very high aspect
ratio.

---

## 4. MODEL — propulsion

Motor torque balanced against propeller torque, solved exactly:

```
motor:  I = (Vt − RPM/Kv)/Rm      Q = Kt(I − I0(n))   Kt = 60/(2π·Kv)
        I0(n) = I0 · (1 − s + s · n/n_test),   n_test = Kv · V_test
prop:   Q = Cp(J, rpm)·ρ·n²·D⁵/(2π)    T = Ct(J, rpm)·ρ·n²·D⁴   J = V/(nD)
```

The propeller is read at the rpm it is turning: the torque balance reads
CP at each trial rpm, so every thrust and power in the mission comes from
the blade at its own rpm (see "Propellers" below).

**The no-load current rises with speed.** Makers quote I0 at one voltage
(`Motor.i0_volts`: 10 V for most, Hacker 8.4 V, T-Motor's AT5330-A 18 V,
Cobra's C-5330s 20 V), and on 6S a motor turns about
twice as fast. Its no-load loss is part hysteresis and bearing friction,
which draw a fixed current, and part eddy current, whose current grows in
proportion to speed. `s` is the eddy share at the maker's test speed
(`i0_speed_share`, 0.5: I0 comes out about 1.6x the quoted figure at 6S).
Measuring I0 at two voltages calibrates it; 0 is the old constant I0.

The ESC gives the motor `throttle` times the pack's loaded voltage, and the
sag is an extra series resistance on the motor ("Battery sag" below), so
the electrical loop does not close on itself and the whole system is one
scalar root-find in `n`. There is no iteration error anywhere in the power
system, and the operating point is solved, never interpolated: the only
tables are the propeller's own (in J and rpm) and the installed-thrust
factor against airspeed.

### Motor data

Every motor comes from the makers' own pages. 129 were entered (99 of
them on 2026-09-29 for 1-4 motor designs, each checked for efficiency,
copper loss and Km-vs-mass before entry). The same day's pruning study
(`results/catalogue_prune_2026-09-29/`) kept the 38 that earn a place at
some motor count: Scorpion, T-Motor, Cobra, KDE Direct, Hacker, SunnySky,
AXI and Neu. Seven more are retired: loadable for saved designs, never
searched. The motor-count experiment (`results/experiment_2026-10-02/`)
then added 14 KDE Direct motors (2026-10-02: the 1806 to the 7215, for
small props on many motors and props over 21 in; the 26 g 1806 was removed
again on 2026-10-04) and 59 low-Kv motors for
the big high-pitch props (2026-10-03: 19 from T-Motor, Scorpion, Hacker,
AXI, Cobra, SunnySky and Dualsky, and 40 Neu winds), each passing the
catalogue's efficiency, copper-loss and Km-vs-mass checks: the catalogue
now holds 117, 110 of them searched. The packs were pruned the same way,
86 to 38 searched (66 since 28 8S and 9S packs were added on 2026-10-05 for the battery study), and the ESCs 55 to 21 (22 with the KDE UAS20LV the
experiment added). Each entry's `source` names the page and the rating
basis. The previous catalogue was unsourced and wrong where it mattered:
a "Cobra C-4120/12 465KV" that does not exist (the real one is 850 Kv), a Scorpion SII-4025-330 at 500 g against a real 353 g, a T-Motor
AT4125-540 at 0.038 Ω against a real 0.014 Ω.

**The makers do not rate on a common basis**, and the model treats their
limits as hard at full throttle while a Mission 2 runs ~250–300 s near full
power:

| maker | current / power limit means |
|---|---|
| Scorpion, Cobra | continuous |
| KDE Direct, T-Motor | 180 s |
| AXI | 60 s |
| SunnySky | 10–60 s (30 s for most) |
| Hacker | 15 s maximum *power* only; current derived as P / (rated cells × 3.7 V), P / 22.2 V for the first five |
| Dualsky | 15 s burst only |
| Neu | "continuous watts" as the power limit; current = min(saturation, W / 22.2 V) for the first five, the saturation "Max Amps" for the 40 low-Kv winds |

So a SunnySky, Hacker or Dualsky motor is held to a looser limit than a
Scorpion.
**CALIBRATE:** a motor that looks marginal on a short-duration rating should
be bench-run for five minutes at the design current before it is trusted.
Neu also rounds `Rm` to 0.01 Ω, one significant figure for its
low-resistance windings.

The 2026-09-29 motors were **screened in the densest expected air** (the
Tucson morning preset), because a propeller absorbs more power in dense air
and the motor limit binds there. The search itself screens and shortlists
at the nominal condition (each ESC sized in the densest and hottest Tucson
air), then flies every finalist on all three presets with a propeller swap
where needed (rule 3.2.1e), gives each a weather-robust version, its
sensor trimmed 1% at a time until all three fly, and recommends the
finalist whose robust version scores most on its worst preset — see
README.md, section 15.

### Battery sag is modelled (restored in v1.16)

Every pack has a DC internal resistance, and each propulsion system's
current also flows through its fuse, arming plug and leads:

    V_loaded = V_oc(SoC) - (R_pack + R_harness) x I_pack

The ESC feeds the motor throttle x V_loaded, which makes sag an extra
series resistance of (motors per pack) x throttle^2 x R on the motor. A
wing motor's leads add throttle^2 x R_lead of their own (below).
The power system is still one exact scalar root-find; it does not iterate,
and sag is not tabulated. Energy is charged at the **open-circuit**
voltage, so the I^2R heat inside the pack and wiring is spent energy too.
That is how sag costs endurance as well as thrust.

| input | value | where |
|---|---|---|
| LiPo cell resistance (DC) | 19 × (50/C)^0.65 mOhm*Ah, fitted to six independently measured packs — 3.7 mOhm for a 4500 mAh 60C cell. The DC figure under a sustained load is ~2× the 1 kHz AC one, falling as the cells warm, so a published AC figure is entered ×1.7 | `BatteryPack.r_cell_ohm` |
| Li-ion cell resistance (DC) | 60 mOhm*Ah — 15 mOhm for a 4000 mAh pack (Molicel P45B 12 DC at 4.5 Ah) | same |
| Measured value | overrides the estimate — `r_cell_mohm` column in the pack catalogue (enter DC, or AC × 1.7) | Catalogues tab |
| Pack heating | ∫I²R per pack over each mission, 1,000 J/kg·K, no cooling; a pack past 60 °C ends the mission (`pack_temp_limit_c`) | `fly_mission` |
| True capacity | a pack's measured capacity, where a test gives one, sets its energy; the 100 Wh rule keeps the label | `BatteryPack.energy_wh` |
| Harness | 3 mOhm per system: 100 A blade fuse, two connectors, ~0.3 m of 10 AWG round trip | `BuildStandard.harness_resistance_ohm` |
| Wing-motor leads | each wing motor's ESC rides in its nacelle, so its DC leads run its spanwise station + 0.20 m from the battery bay, in copper sized at 12 A/mm2 of the motor's rated current: R = 2 L rho/A, mass 2 L A x copper density x 1.25 (jacket). The mean per motor goes into the power solve; a nose motor's are the harness | `layout.wing_power_leads` |
| Wing nacelles | 15% of the motor + ESC mass each pod carries (none for a nose motor) | `BuildStandard.nacelle_mass_frac` |
| Nacelle junctions | Hoerner junction drag on each side of a wing nacelle at the local section thickness, inside the slipstream's scrubbed drag | `AeroModel.nacelle_junction_m2` |
| Nacelle drag | `nacelle_cd` 0.20 on the frontal area of a pod 1.25 × the motor's can diameter (section 3); a measured drag area overrides it | `AeroModel.nacelle_drag_area_m2` |
| Motor mounts | each motor's mount and spinner: 25 g for a 62 mm can (about a 350 g motor), scaled with the can's face area | `BuildStandard.motor_mount_kg` |

What it does (6S, 60% charge, full throttle, static; the current code and
catalogue):

| system | R total | current, no sag -> sag | thrust |
|---|---|---|---|
| SunnySky X4120 + 13x10, GA 6S 4500 60C | 25.5 mOhm | 57.9 -> 52.0 A | **-10.4%** |
| 2 x Scorpion 4025-440 + 13x8, 2 x TP Elite 6S 2250 55C | 50.6 mOhm per system | 72.1 -> 63.2 A | **-12.8%** |

Sag also *lowers* the current, so fewer designs hit the 100 A fuse.
Two packs each carry half the current and sag half as far, which the model
now rewards.

### What differs between 1, 2, 3 and 4 motors

Modelled, for every motor count:

- each motor's own torque balance, current and power limits, winding
  temperature and ESC; each pack's own sag, heating and fuse;
- the no-load current at each motor's actual speed (above);
- a nose motor for an odd count and the rest in wing nacelles
  (`layout.motor_stations`), with the gear sized to the lowest propeller.
  The wing propellers sit in the middle of equal shares of the exposed
  semi-span, on the wingtips (the outermost with its axis on the tip, its
  disc overhanging) or, with two or more a side, staggered fore and aft
  (every other nacelle 0.6 D ahead on a boom; neighbours may overlap up
  to 30% of a diameter, the overlapped share of an aft disc losing thrust
  in the forward jet at every airspeed); the search tries every layout a
  count allows;
- per wing motor: its nacelle's drag (Cd 0.20 on the pod's frontal area),
  its two Hoerner junctions, its pod's mass (`nacelle_mass_frac`), the
  pod's weight relieving the spar inboard of it, and its power leads'
  length, resistance and mass from its station (`layout.wing_power_leads`);
- each motor's mount and spinner, scaled with its can;
- slipstream scrubbing per propeller (a nose propeller washes the whole
  fuselage; a wing propeller its strip of wing and its nacelle), the
  slipstream's lift with the fuselage or nacelle taking its middle, its
  turbulent transition, swirl recovery behind the counter-rotating,
  inboard-up wing pairs, the roll trim against an unpaired nose motor's
  torque, windmilling per propeller;
- per system: pack, fuse harness and arming plug mass; energy in whole
  packs.

Not modelled, with why:

- compressibility at the propeller tips beyond what the tables carry
  (each propeller's tested rpm levels, and APC's computed rpm trend above
  them): a tip Mach number past 0.82 is rejected, and the leading designs
  turn their tips at Mach 0.35–0.53 static (examples/optimised_v1.json
  0.53, results/exp8_C_cap2_2026-10-03 0.39);
- a nose propeller's swirl (antisymmetric across the span, it cancels to
  first order in the symmetric lifting line) and swirl's effect on the
  stall; swirl recovery behind wing propellers and the roll trim against
  a nose motor's torque are modelled (above);
- ESC efficiency against current: 96% flat, since no maker publishes it;
- a folding propeller folding back when the power is cut: there is no
  folded state, so the CAM folding props are charged a freewheeling disc's
  drag in the glide like any other;
- an rpm limit for the APC props: the catalogue has none for them (APC
  publishes its own limits, which were not entered), so only the tip
  Mach number is checked for them; the CAM props' limits are (above);
- a blown fuse on one wing system leaves asymmetric thrust: nothing sizes
  the fin for it - check the rudder can hold it on approach.

On the leading designs of 2026-09-30, each re-fitted (propeller and
payload), the new terms cost the 3-motor design 1.2% of its score and the
1-motor design 0.6% (results/thorough_2026-09-30/motor_count_model_check).

**The resistance is estimated, and it is the number to measure first.**
Most chargers read cell IR. A load-step test (volts at rest, then volts at
a known current) gives the DC value directly. Enter it per cell in the
`r_cell_mohm` column.

**The ESC low-voltage cutoff is a hard limit.** If the loaded voltage at
full throttle in the turn falls below 3.2 V/cell (LiPo) or 2.8 V/cell
(Li-ion), the ESC would cut the power, so the flight ends there. A pack
that sags below cutoff on a full charge cannot fly the aircraft at all.

**Motor heating is modelled** (propulsion.motor_thermal):
- **Model:** one lumped thermal mass per motor, mass × 500 J/kg·K. It
  heats with the motor's losses (electrical input minus shaft power)
  and sheds heat through `R = k / m^(2/3)`.
- **k = 0.168 K/W·kg^(2/3):** the median of the continuously-rated
  catalogue motors (Scorpion, Neu), taken as reaching 110 °C at rated
  current from 25 °C.
- **Temperature effects:** windings start at air temperature (the presets
  carry 12 / 28 / 36 °C), and `Rm` rises 0.393%/°C.
- **Failure:** a mission that would take the windings past 110 °C fails.
- **Short makers' ratings** (15 s Hacker and Dualsky, 10–60 s SunnySky,
  60 s AXI, 180 s KDE and T-Motor) are burst limits, not thermal data.
  They remain hard limits; cooling is the common, size-based estimate.
  CALIBRATE `motor_thermal_k` from a winding temperature measured after a
  full-length run.

### Battery usage: 95% per mission

Every mission flies a fresh, fully charged pack. So the share of a pack
one mission may use is set by what the team is willing to do to a pack,
not by what a later flight needs.

- `BuildStandard.depth_of_discharge` = **0.95** (Advanced panel), the
  team's choice: takeoff, laps and landing together may use 95% of rated
  energy, landing with 5%, about 3.5 V/cell at rest.
- Each chemistry has a ceiling, `BatteryPack.usable_fraction` = 0.95, and
  the setting cannot go past it.
- **The cost:** 5% of a 100 Wh budget is about 5 Wh, **about a third of
  one Mission 2 lap** (about 14 Wh). There is no energy left for a
  go-around, and deep discharges age packs, so use fresh packs for scored
  flights.
- The **ESC low-voltage cutoff** (3.2 V/cell loaded) is still enforced,
  so a sagging pack cannot use the last few percent anyway.

What it buys, on the v1.15 recommended airframe with nothing else changed:
the heaviest sensor that still finishes Mission 2 goes from 22.3 lb at 80%
to 25.1 lb at 90% (gross 35.9 to 39.1 lb), and 26.3 lb at 95%, where the
300 s window binds instead. At the design's own 21 lb sensor, Mission 3
goes from 4 laps to 5.

**The landing point: 88%** (`land_at_discharge`, the team's call on
2026-10-01). Once any pack in any propulsion system has used 88% of its
energy, the airplane turns for landing. No lap is started that would take
a pack past it, and the landing is flown out of the 12% that is left. All
packs are identical and carry the same current, so they reach it
together. Every mission therefore crosses the line for the last time with
at least 12 Wh of 100. That covers the landing circuit (about half a lap,
roughly 8 Wh) and some model and pack error. The 95% depth of discharge
stays the floor nothing may pass.

It replaced Mission 2's 2% planning margin (`m2_energy_margin`, now 0,
still a knob).

**Every flight charged at full throttle** (`full_throttle_flights`, the
team's call on 2026-10-02). On top of the landing point, every second of
every flight is charged at full-throttle power and heat, so the pilot has
no cruise speed to hold and cannot run a pack short by flying harder than
the plan. Against planning each flight's speed, it cost the two leading
designs about 8% of score before the search re-optimised them: neither
could fly Mission 2 at any payload on its own propeller, and the best on
smaller ones (a 16x8 on the Neu 6521, 12x8s on three KDE 4014XFs) fly
Mission 2 in 260-300 s, against the 300 s window. `BuildStandard.energy_reserve` is 0 for every mission.
Previously 8% was withheld on top of an 80% limit, leaving a flight 73.6%
of the pack.

The landing's energy is budgeted **before** the last lap is started. A
flight never flies a lap that would leave it unable to land within the
budget, so a battery-limited Mission 3 ends within one lap of the limit.

### Remaining deliberate pessimism

- **Installed thrust is computed** from slipstream scrubbing, `T (1 −
  f/A_disc)`: the airframe of drag area `f` inside the slipstream drags
  `T/A` of extra dynamic pressure. It is worked out at every airspeed from
  standstill to the propellers' speed ceiling (the skin friction of `f`
  falls with speed), floored at 0.80, and interpolated in the power solve.
  Example (examples/optimised_v1.json): 0.937 static, 0.959 at its
  31.4 m/s top speed; examples/baseline.json: 0.884 and 0.928. It was a
  flat 95%.
- **Windmilling drag is computed** from each propeller's table, extended
  to zero torque (0.014–0.11 on disc area across the catalogue; the props
  under 8 in, whose short tables extrapolated to 0.54, were removed on
  2026-10-04; the 0.08 in the build standard is only the fallback for a
  propeller without a table).
- **Usable battery energy is the depth of discharge**: 95%, landing
  included (see "Battery usage" above).
- ESC efficiency 96% (not published by any maker). The ESC is the
  lightest that holds its motor's peak current with 20% headroom for the
  300 s mission window at 10 m/s of cooling air (`esc_airflow_mps`), in
  the densest and hottest Tucson air. Its published rating is moved to
  those conditions along Scorpion's derating table (the only one
  published): a 60 s / 32 m/s rating keeps 56% (58% at the nominal
  28 °C); Castle's 5-minute rating keeps 97% (100%).

**Propellers are UIUC wind-tunnel data where it exists, and APC's
computed data corrected to it elsewhere.** 81 props, each a `(J, CT, CP)`
table (static point plus the J sweep), none under 8 in. 59 are measured by
UIUC: 27 of the 31 APC Thin Electric sizes from 10 to 21 in, the 8–9 in
Thin Electrics, five 8 in APC Sport props, and 20 Aero-naut CAM Carbon
folding props of 12–16 in (UIUC volume 3). (The props under 8 in - GWS,
Micro Invent and the smallest APC Sport and Free Flight - were removed on
2026-10-04: none can carry this aeroplane, none reached a finalist, and
their short tables gave unreliable windmill drag.) 22 are
APC's own computed tables corrected by the ratio their measured neighbours
needed: the four Thin Electric sizes UIUC did not test (15x8, 15x10,
17x10, 20x13), the eight over 21 in (to 27x13) and ten high-pitch props of
17–21 in (P/D 0.65–1.1). The APC 15x13 was removed: APC makes no such Thin
Electric prop, and its table had been extrapolated from the 16x12. The old
parametric fit overstated static CT by 10–50% and static CP by 30–80%.
Each prop carries its maker's listed weight (a CAM prop: its blade pair
plus a 40 g yoke and spinner). The CAM props also carry Aero-naut's rpm
limits (9,000–13,000 rpm): a set whose motor could drive its propeller
past the limit unloaded - the no-load speed at full charge, which a
full-throttle dive approaches - is rejected, in the screen and in the
aircraft's checks, as an over-current motor is.

**Each propeller is read at the rpm it is turning**
(`dbf2027/prop_rpm_tables.json`): linear between its two nearest rpm
levels, held at the nearest level outside them. The measured 12–21 in
Thin Electrics carry every rpm UIUC tested, each sweep carried past its
last point on the next lower rpm's measured shape; above the fastest rpm
tested, APC's computed rpm trend for the same prop is applied to UIUC's
level, which errs low (measured efficiency rises 3.2% per ~1000 rpm, APC's
code says 2.0%). The CAM props have no computed data and are held at
their fastest tested rpm (5,000–7,000), which also errs low. The
APC-computed props carry every 1000-rpm block from 2000 rpm. The 8–11 in
Thin Electrics and the 8 in APC Sport props have one table, the
highest-rpm sweep. Thrust-stand your
shortlisted prop and load it with `--props-csv`; `--props` restricts a
search to named propellers.

**The propellers' speed ceiling.** A propeller makes no thrust past its
zero-thrust advance ratio (the largest over its rpm tables), and it turns
no faster than the motor's no-load rpm on the pack's open-circuit voltage,
so there is no level flight above `J_zero · n_no-load · D`. That ceiling
bounds the top-speed searches, the installed-thrust and full-throttle
power tables and the quick reject (section 7), in place of the fixed
speeds used before 2026-10-03. Example
(examples/optimised_v1.json, an APC 16x8 on a 362 Kv Neu, 6S): 38.6 m/s on
a full pack.

### Battery packs are catalogue items

Custom packs are banned and packs may not be paralleled, so "cell × series ×
parallel" is not a design variable any more. A pack is an atomic entry: 80
real, purchasable packs from 22 to 99.9 Wh (72 searched, 8 retired by the
2026-09-29 pruning; 28 of them the 8S and 9S packs added on 2026-10-05, 6
more - 7S, 13S and 8S/9S Li-ion - from a second search that day), each
carrying the four numbers its label must state plus a mass.

**LiPo, 3S to 13S (75 packs, six of them LiHV).** The first were the
most energy-dense COTS packs at each cell count under the 100 Wh per-pack
limit, from about 1,100 maker listings (MaxAmps, Thunder Power, CNHL,
retrieved 2026-09-27); more came from the makers' pages and independent
tests on 2026-09-29. They run 116–179 Wh/kg. Six older 6S entries are
kept, retired, because saved designs name them.

- **Why higher voltage helps.** The cap is in Wh, so every cell count
  carries the same energy. But twice the voltage is half the current for
  the same power: half the load on the 100 A fuse and a quarter of the
  I²R sag loss.
- **The limit.** At 9S–13S the only packs under 100 Wh are MaxAmps'
  2,000–2,800 mAh line (67–96 Wh, 25–27C); a second search of about 20
  makers on 2026-10-05 found no other (several makers' sites could not be
  reached). Every high-C pack at those cell
  counts starts above 100 Wh (Gaoneng 12S 2800, HRB 10S 3000, Thunder
  Power 10S 3300).
- **7S.** The one-piece pack nearest the cap is MaxAmps' Graphene 3750
  (97.1 Wh, 142.5 A). Thunder Power's 7S packs are sold split, two
  sub-packs joined by a series interconnect, which rule 3.2.3.1b may
  forbid - confirm before using one.

**A high-voltage pack needs hardware rated for it.**
- **Motors:** each carries `max_cells`, its maker's rated cell count, and a
  pack with more cells is rejected. They run from 4S to 14S: the large
  KDE motors, the Hacker A60-6XS, the Scorpion SII-4035-250 and -380,
  most T-Motor AT motors and most of the Neu winds take 12S or more; the
  Neu motors run 7–12S (from their maximum voltage), the small KDE
  multirotor motors 4–6S, most SunnySky motors only 6S. The Hacker
  A40-14L and A50-16S pages give no rating, so they are held to 6S until
  verified.
- **ESCs:** 12S jobs get an HV unit: Castle's Phoenix Edge HV or Edge
  Lite HV (40–160 A), or KDE's UAS55HVC or UAS95HVC.

**C rating.** MaxAmps publishes a continuous current rather than a C
number, so C is that current divided by capacity (25–32C).

**Li-ion (5 packs, 4 searched)** runs 193–217 Wh/kg at 10–11C. That is a
genuine trade against LiPo's current, and the optimiser makes it rather
than being told the answer. The two above 6S (Nexus Battery Systems 8S
and 9S 1P Molicel P30B, 86 and 97 Wh, 30 A) are fixed catalogue products
but built to order - confirm against rule 3.2.3a.

**Every entry is marked *VERIFY*.** Masses and C ratings drift between
production runs — confirm against the current datasheet before ordering, and
load your own shortlist with `--packs-csv`.

NiCd and NiMH are legal under the rules but are not in the catalogue: at
these power levels they lose on both mass and current, and the team asked
for lithium only.

---

## 5. CALIBRATE — structures and weight

The mass breakdown is **itemised** — wing skin, wing spar, motors,
propellers, ESCs, battery and fuselage all appear as separate lines rather
than a lumped "propulsion group", so you can check each against reality.

- **Stall speed: a limit on the factored lift (m101).**
  - **The check:** `max_stall_mps` is the fastest stall the team will
    fly. It's taken power off, at the Mission 2 weight (the heaviest),
    with the wing's maximum lift divided by `stall_lift_safety_factor`:
    the same 1.5 the wing loads carry. So the aeroplane must hold its
    weight at the limit speed on 2/3 of its lift, and its true stall is
    about 1/√1.5 = 0.82 of the limit.
  - **In the hottest weather preset's air (m118).** The limit is on the
    true airspeed, so it is checked where the stall is fastest:
    `stall_check_conditions` = `tucson-hot` (36 °C) by default,
    `--stall-air` on the command line, `mission` for the mission's own
    air.
    - **Why:** the v3 sweep sized every design on the nominal day (28 °C)
      and put it exactly on the limit there. The hot day's air is 2.6%
      thinner, so the same wing stalls about 1.5% faster, and every
      stall-limited finalist broke the limit by 0.5–0.9 mph on the hot
      day. It took 6–11% less sensor to fly all three presets.
    - **How:** the same wing, trim and CG at that air's density and
      viscosity (`AeroModel.v_stall_in_air`), so the Reynolds number
      moves with it. That is exactly what the weather check's hot preset
      finds, so a design on the limit flies every preset.
    - **What it costs:** about 1.7% of score at 40–50 mph with only the
      sensor trimmed (v3 winners), before the search re-optimises around
      it.
  - **How it's solved:** as `AeroModel.v_stall(W, load_factor=1.5)`, so
    CLmax's Reynolds dependence and the trim download are included.
  - **Default:** off (inf). `--max-stall-mph` sets it, and the GUI's
    Advanced tab shows both numbers.
  - **Trade study:** `results/stall_limit_2026-10-08`.
  - **Searching under a tight limit (m107, m108).** With a chord cap or
    a stall limit, the coarse sweep keeps each sample's aspect ratio at or
    above span ÷ max chord. A sample that fails only on stall has its
    sensor brought down to the heaviest that meets the limit, and is
    scored there. Otherwise the flyable region (full span, AR at the cap,
    a light sensor) is too small for random samples to hit. Without
    either limit, nothing changes.
  - **The fit is a secant (m120).** The stall goes very nearly as the
    square root of the Mission 2 weight, so v² is close to linear in the
    sensor: a first step from that law and a bracketed secant land within
    10⁻⁶ of the limit in 2-4 aircraft builds, where the bisection it
    replaced took 15 and stopped within 3×10⁻⁵. The aeroplane built for
    the stall test is the one flown. Under a stall limit 81% of coarse
    samples need the fit, and it was most of the coarse sweep's time: the
    sweep's worker is about 1.7× faster, with the same samples flying.
  - **Nothing in the coarse sweep flies (m121).** The search goes on
    from its seed designs - scored and made flyable as always, and the
    run normalised against them - and stops only if none of them flies
    either. v6's 25 mph search sampled nothing that flew and stopped
    before trying v5's 25 mph winner, which flies it.
- **Spar: a 1.5 factor of safety on every wing load case.** A carbon
  tube, sized at 24 stations along the half-span for the worst of:
  - flight bending at `n_ult = 1.5 × n_limit`, with elliptic lift (with
    winglets, the end-plated distribution: m117);
  - the wing-tip lift test at MGTOW (rule 4.2.4), × 1.5. Its moment,
    `(W/2)(b/2)(1 − η)`, falls off linearly, so it governs the outboard
    spar even where flight governs the root;
  - shear plus torsion (section `c_m` at the dive speed), × 1.5, all
    carried by the tube.

  Each wing motor's pod (motor, propeller, ESC, mount and nacelle) hangs
  on the spar at its station, and its weight - times the load factor in
  flight, at 1 g in the tip test - takes off the shear and bending inboard
  of it.

  Thin-wall tube: bending `t = M/(σπr²)`, shear `V/(πr)/τ`, torsion
  `T/(2πr²)/τ`. The radius follows the local airfoil depth
  (0.40 × thickness), so the tube tapers with the wing. Minimum wall
  0.6 mm, ×1.8 for webs and fittings, mass integrated along the span
  (this replaced a flat ×0.6 taper factor).
  σ 600 MPa and τ 50 MPa are **failure** stresses, because the 1.5 is in
  the loads (CALIBRATE both to the tube you buy); ρ 1600 kg/m³. No
  inertia relief from the wing's weight. A wing whose tube cannot close
  somewhere on the span is rejected.
- **Areal densities (CALIBRATE):** wing 1.10 kg/m², tails 1.80 kg/m²
  (stabiliser and fin set separately; an inherited effective value, above
  the wing's, so weigh your own), fuselage 1.40 kg/m² of wetted area.
- **Landing gear** 0.080 kg + 3% MTOW, plus 0.55 kg/m of extra leg when a big
  propeller forces taller gear.
- Fixed: 6 servos at 20 g, 0.10 kg avionics, 0.15 kg sensor deploy mechanism,
  8% of payload for restraint structure, **8% contingency**.
- Per motor: its mount and spinner, 25 g for a 62 mm can scaled with the
  can's face area (`motor_mount_kg`), and its propeller at the maker's
  listed weight (section 4).

**Not modelled:** aeroelastic deflection, fatigue, landing loads, build
variability. (The CG is placed, not assumed: section 2.)

---

## 6. MODEL — mission and operations

- **Every second charged at full throttle** (`full_throttle_flights`,
  2026-10-02): the energy and heat of every second from throttle-up to the
  last crossing of the line are those of full throttle at that airspeed,
  the glides and turns included. A factor of safety: there is no cruise
  speed for the pilot to hold. The landing circuit keeps its own power.
- **Straights:** full throttle on the power-on drag, then power off at the
  last moment that still allows arrival at the turn at corner speed (flown
  as a glide, charged as full throttle). The glide is planned as it is
  flown: over the ground, with the wind, and with the freewheeling
  propellers' windmilling drag on top of the airframe's.
- **Turns:** sustained level turns at the speed that maximises turn rate,
  at full throttle with the slipstream's lift, limited by whichever of
  CLmax, thrust, structure (the derived limit load factor, section 2a) or
  bank (≤ 70°) binds first. No energy trade through the turn.
- **Takeoff:** μ = 0.035, ground CL = 0.35 (raised by the slipstream's
  lift), rotate at 1.15 × the full-throttle stall speed, induced drag
  reduced by the standard ground-effect factor. The acceleration is
  checked at 40 speeds between brake release and rotation first: a roll
  whose acceleration reaches zero below rotation speed fails at once
  instead of being integrated. Ground roll is reported even though there
  is no rules limit, because a 300 ft roll is an operational risk
  regardless.
- **Climb** to a 30 m pattern at full throttle, by energy height: first
  the speed change from liftoff to the best-climb speed, flown level, then
  the climb at that speed's specific excess power `V (T − D)/W`. Its
  energy is always charged. A design with no excess power at its
  best-climb speed, or anywhere between liftoff and it, cannot climb, and
  the mission fails there. The ground distance covered is credited against
  the first lap, at the lap's mean ground speed and never more than the
  climb's own time.
- **Landing is simulated** — descent at 3 m/s, approach at 1.3 V_stall, a 3 s
  flare, then a braked rollout (μ = 0.30) integrated with RK4. Its energy
  (the descent, on the throttle for 35% of the approach drag; the flare and
  rollout at idle) is charged to the battery and a design that cannot land
  on its remaining charge is infeasible. Its duration is reported and
  included in `total_flight_time_s`.
  **It is deliberately NOT added to the scored mission time**, because the
  rules you quoted stop the clock as the aircraft crosses the line at the end
  of the final lap: "landing is not part of the 5 minute time window". If
  that reading is wrong, one line in `fly_mission` moves it into the total.
- **Wind** is a steady component along the runway; it changes ground speed on
  the straights and correctly does not change turn rate. No gusts, no
  gradient, no crosswind.
- State of charge is held constant within a lap and re-simulated when it
  drifts by more than 0.08, or when the windings warm into another 2 °C
  band. The sustained turn and the top speed are worked out once per
  flight, at 70% charge.
- **Pilot efficiency 0.95** (`pilot_efficiency`, the team's call on
  2026-10-02; 0.93 before): each lap's time and energy are divided by it,
  5.3% off a flawless lap. One number, the same on every lap - calibrate
  it from flight logs.
- **Not modelled:** pilot variability, dropped laps, the sensor light-mode
  sequencing (assumed always achievable), the risk of failing to stow.

---

## 7. NUMERICAL — measured, not asserted

Run `python grid_convergence.py` to reproduce.

| Source | Method | Measured error |
|---|---|---|
| Mission ODE (lap time, energy) | RK4 + bisection event location | **1×10⁻⁵ s** across dx from 24 m to 0.75 m |
| Propulsion operating point | Brent on the torque balance | **5×10⁻¹⁴ N·m** residual; **3×10⁻¹³ rev/s** vs independent bisection |
| Mass fixed point | undamped iteration, guard | **0.000 kg** residual |
| Glide-to-corner distance, climb speed-up | 5-point Gauss–Legendre | exact to polynomial degree 9 |
| Corner speed, best climb, best L/D | scan + golden section | bracketed to 10⁻⁴–10⁻⁶ |

**What was wrong and got fixed:**

1. The straight-leg integrator was **explicit Euler**, first-order, and at
   coarse steps had an outright bug (when the time step hit its cap the
   distance advanced anyway) giving a **42 s / 20% error** at dx = 16 m.
   Since screening ran there, it was biasing which designs survived. Now RK4
   with both events located by bisection: step-size independent and faster.
2. Battery sag was iterated, then tabulated, then interpolated — two stacked
   error sources, 7–20% RMS at part throttle. It was reduced to a closed
   form, making the power system one scalar root-find, exact and 4× faster;
   the interpolation table measured *slower* than the exact solve, so it was
   deleted. Sag was later removed at the team's direction and then restored
   in v1.16, still as the same closed form: an exact, single scalar
   root-find (see "Battery sag is modelled").
3. Corner speed and best climb speed were bare grid scans returning the best
   *sample*. Now bracketed and refined.
4. The mass loop was damped 50/50, raising its contraction factor from 0.04
   to 0.52. Undamped with a divergence guard: four iterations instead of 20+.
5. **Adding the landing fields to `MissionResult` shifted its positional
   arguments**, so four failure paths wrote their reason string into
   `landing_time_s` and left `reason` empty. Fixed with keyword arguments.
6. The search **sampled span and chord independently** (aspect ratios from
   1.6 to 20) and **simulator count uniformly over 0–40** (expected payload
   ~32 kg). 98% of samples were unflyable. Now aspect ratio is sampled
   directly and the chord derived, and payload is sampled against the
   lifting capacity the power system actually has: the least of a static
   thrust/weight of 0.30, the energy to fly the shortest Mission 2 course
   (every turn at its tightest radius, at the slowest speed that still
   finishes in the window) and the 55 lb limit, less a rough empty mass
   and the loaded container.
7. A **cheap pre-screen** rejects hopeless designs in 0.02 ms instead of
   145 ms. It uses genuine physical *lower* bounds — minimum drag from the
   parabolic polar, `D_min = 2W√(CD0·k)` at the lowest CD0 the aeroplane
   can reach (at its speed ceiling), and the propellers' speed ceiling
   (section 4) against both the stall speed and five laps' straights — so
   it cannot reject something the simulator would accept. The drag floor
   is taken at **0.80** of that minimum: with the propellers making thrust,
   the slipstream's lift and swirl recovery can take the drag below the
   power-off polar (to 0.92 of it, measured on 1–2 motor designs), and
   0.80 keeps the floor a necessary condition. Audited against the full
   simulator when it was introduced: **zero false negatives.**

**Remaining numerical approximations**, accepted deliberately: SoC frozen
within a lap; the first-lap climb credit is a distance correction rather than
an integrated climbing lap; turn ground-track wind drift ignored (turn *rate*
is exact); Nelder-Mead is local, and the staged search with a diverse
shortlist is the defence, not a guarantee. The propulsion screen rates
each motor, propeller and pack set on a Mission 2 flown - takeoff, climb,
top speed, fastest sustained turn, laps with the wind - on the best known
airframe of its motor count (the best seed; a count with none borrows the
best seed's airframe), its weight changed by the difference in the
propulsion installation as the full model weighs it - motors, propellers,
ESCs, mounts scaled to the motor, packs, fuse harnesses, wing nacelles and
power leads, the gear a bigger nose propeller needs, all with their
contingency - and its power solve carrying its own wing leads' resistance
(since 2026-10-04: before, every set carried the reference set's nacelles,
leads and contingency). The airframe itself (wing, spar, drag) stays the
reference's, so a set that would suit a very different airframe can rank
low.
With no seeds it falls back to the set's thrust power at its own best
speed above the slowest airspeed that could finish in the window. Since
2026-10-01 every search also refines the known good designs (`examples/`,
`--seeds`) ahead of the sampled ones - at most six per motor count, each
made flyable by trimming its sensor 3% at a time, and in a restricted
search (`--props`, `--motors`, `--airfoils`) only those whose parts it
allows - so it never finishes below a seed it refined.

**Conditioning (audited 2026-10-01).** Every formula was checked for
rounding error that the arithmetic itself amplifies:

- **Linear solve:** the lifting-line matrices have condition number 1,193
  at worst (24 terms, AR 1.5, taper 0.05; infinity norm) and 401 with the
  scaled partial pivoting `solve_linear` now uses (806 and 284 before the
  2026-10-03 lifting-line fix). That bounds the relative error near
  1e-13. Scaled pivoting also keeps a merely rescaled row from taking a
  near-zero pivot (a test case on which plain partial pivoting was 100%
  wrong).
- **Cancellation, rewritten:**
  - slipstream increment: sqrt(1 + x) - 1 as x / (sqrt(1 + x) + 1);
  - jet half-width: (r - h)(r + h) instead of r^2 - h^2;
  - turn rate and bank: sqrt((n - 1)(n + 1)) and atan instead of
    sqrt(n^2 - 1) and acos(1/n);
  - winding step: -expm1(-x) for 1 - exp(-x);
  - density altitude and Richardson / GCI: `expm1`;
  - tall-gear drag: a difference instead of a ratio minus 1;
  - elliptic spar moment and shear: their tip angle taken as
    2 asin(sqrt((1 - eta)/2)), and series near the tip. The old moment
    was 11x too large at eta = 1 - 1e-6, where both of its terms tend to
    theta^3/3 and their difference goes as theta^5/15.
- **Results:** none moved. The spar stations agree with the old formulas
  to 2e-16, and each fix is checked against a 60-digit reference in
  `selftest.py`.
- **Physical, not numerical:** some quantities are ill-conditioned in
  themselves. Climb rate and acceleration near a ceiling or top speed are
  small differences of thrust and drag. Bank angle as the load factor
  nears 1 has an infinite slope. A few percent of model error in thrust
  or drag moves those far more than rounding ever could.

---

## 8. The thing to decide next

With no cap on size or payload mass, the optimiser takes all the span the
rules allow and as heavy a sensor as the energy and the window allow, and
two results need your judgement rather than more physics:

**Span is settled: 6 ft, enforced as a hard rule.** The search never samples
beyond it, and an over-span wing is rejected. This is what stops the wing
growing - area scales as span² while the spar only scales linearly, so
without the rule the optimum ran away to 10 ft. It is the wing's own span:
a wingtip propeller may overhang it (section 1).

**Sensor mass runs as high as the energy and the window allow, and it
beats carrying simulators.** The sensor counts toward the Mission 2 scored
weight *and* multiplies the Mission 3 score, while a simulator only does the former. So the model's
advice is "make the sensor as heavy as you are allowed to", and it prefers
one huge sensor to a fuselage full of simulators. A 25 lb "sensor"
(examples/optimised_v1.json; 31 lb in results/exp8_C_cap2_2026-10-03) is a
brick with lights on it. **The binding constraint here is the rules'
definition of what the sensor has to be and do — deployable, stowable,
remotely commanded lights — not physics.** Find that definition and set
`--sensor-max` from it; it will change the answer more than anything else in
this document.

Then, in order: thrust-stand data for your shortlisted prop at flight rpm
(above the fastest rpm UIUC tested, the catalogue carries APC's computed
rpm trend, which errs low), structural areal densities
(weigh last year's wing), and the sensor
density assumption (800 kg/m³), which sets the sensor's estimated length -
checked against the 6 in minimum - and the density the report quotes; the
container itself is one fixed box (section 2).
