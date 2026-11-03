# Stall-limit sweep (v6): the pick at each limit

The design files here are the best design at each true stall limit from the
v6 sweep (run on Sapelo2, job 49309634, October 2026). Each one is the
highest-scoring design that flies all three weather presets as built, the
same pick the results workbook and the proposal figure use. Each file has
been re-scored from the JSON and reproduces the sweep's numbers.

The search found nothing at 20 or 25 mph, so there are no files for those.

This folder is not read by the search: only the designs at the top level of
`examples/` are used as seeds.

## The designs

Every pick has 1 motor and 1 propulsion system (one pack, one fuse, one
arming plug). All six use a 6 ft span, a rectangular wing and a mid-wing
motor.

| Limit | File | Airfoil | Chord | Motor | Propeller | Battery | Sensor, M2 (M3) | M2 gross | Stall, hot day | Score of 7.5 |
|---|---|---|---|---|---|---|---|---|---|---|
| 30 mph | `stall30mph.json` | S1223 | 22.8 in | Neu 6510/11/242 | APC 20x13 | Tattu R-Line 6S 4500 95C | 18.51 (18.51) lb | 37.6 lb | 30.0 mph | 6.10 |
| 35 mph | `stall35mph.json` | S1223 | 18.5 in | Neu 6514/9.5/200 | APC 25x12.5 | Tattu R-Line 6S 4500 95C | 23.71 (23.65) lb | 41.6 lb | 35.0 mph | 6.61 |
| 40 mph | `stall40mph.json` | S1223 | 15.1 in | Scorpion S-5028-220 | APC 22x12W | Tattu R-Line 6S 4500 95C | 27.51 (27.44) lb | 44.4 lb | 40.0 mph | 6.96 |
| 45 mph | `stall45mph.json` | S1223 | 12.4 in | Neu 6514/8/237 | APC 20.5x14EPN (20x13 on the cool morning) | Ovonic 6S 4500 100C | 30.45 (30.37) lb | 45.9 lb | 45.0 mph | 7.19 |
| 50 mph | `stall50mph.json` | S1223 | 10.9 in | Neu 6514/8/237 | APC 20.5x14EPN | Ovonic 6S 4500 100C | 31.56 (31.56) lb | 46.5 lb | 48.3 mph | 7.26 |
| none | `nolimit.json` | SD7032 | 12.1 in | Neu 6514/6.5/292 | APC 18x12 | Tattu R-Line 6S 4500 95C | 34.09 (30.86) lb | 49.6 lb | 55.0 mph | 7.34 |

The catalogue keys in the files are `NEU_6510M_11_242`, `NEU_6514M_9P5_200`,
`SCORPION_S5028_220`, `NEU_6514M_8_237`, `NEU_6514M_6P5_292`,
`TATTU_RL_6S4500_95C` and `OVONIC_6S4500_100C`. Every part is defined in
`dbf2027/components.py`.

The stall is the true power-off stall at the Mission 2 weight, in the 36 C
afternoon air (lift / 1.0). The 50 mph pick doesn't reach its limit.

On the results workbook's scale (exp8), the scores are 6.341, 6.927, 7.325,
7.584, 7.673 and 7.767.

## Propulsion, for checking against

All at full charge and full throttle, in the nominal Tucson air (28 C, 667 m:
1.082 kg/m3, a density altitude of about 4,170 ft).

| Limit | Static thrust | Static current | Pack voltage under load | Static rpm | Peak current, and the airspeed it peaks at |
|---|---|---|---|---|---|
| 30 mph | 13.2 lbf | 43.2 A | 24.35 V | 5,233 | 47.6 A at 33 mph |
| 35 mph | 18.9 lbf | 56.9 A | 24.08 V | 4,133 | 60.2 A at 23 mph |
| 40 mph | 18.6 lbf | 60.6 A | 24.01 V | 4,855 | 62.7 A at 23 mph |
| 45 and 50 mph | 14.9 lbf | 48.9 A | 24.26 V | 5,252 | 54.0 A at 36 mph |
| none | 14.4 lbf | 51.8 A | 24.18 V | 6,621 | 57.4 A at 39 mph |

The current peaks above zero airspeed because a measured propeller's power
coefficient peaks there. In the denser 12 C morning air every current is a
little higher: the 30 mph pick peaks at 50 A.

## Battery voltage

There is no fixed voltage. Each pack follows an open-circuit voltage curve
against its state of charge, and sags under load:

- **Open circuit, per cell (LiPo):** 4.20 V full, falling to 3.35 V empty.
  For these 6S packs that's 25.2 V full, 24.3 V at 70% charge and 21.8 V at
  the landing point (12% left).
- **Under load:** the pack voltage is the open-circuit voltage minus the
  current times the pack's internal resistance plus the harness (about
  0.02 ohm on these packs). That's solved together with the motor and
  propeller, so it's about 24.0-24.4 V at full throttle on a full pack.
- **Which charge each calculation uses:**
  - the current, fuse, ESC and motor-rating checks use full charge, full
    throttle, at the airspeed that draws the most current;
  - takeoff uses full charge;
  - the laps track the charge lap by lap from the energy used, so the
    voltage and thrust fall as the flight goes on;
  - the turn and top-speed capability use 70% charge;
  - the ESC low-voltage cutoff is checked on every lap.
- High-voltage (LiHV) packs in the catalogue use 4.35 V full per cell. None
  of these picks use one.

## Scoring a design

From the top folder of the project:

```bash
python dbf.py --evaluate examples/stall_limits_v6/stall30mph.json --max-stall-mph 30 --stall-lift-fos 1.0 --stall-air tucson-hot
```

Change the file and the `--max-stall-mph` number for each limit. For
`nolimit.json`, leave off the three stall flags. It takes a few seconds and
prints the weights, propulsion, each mission, the weather plan and the score.

**Leave `--stall-lift-fos 1.0` in.** The default is 1.5, which checks the
limit against 2/3 of the wing's lift. Without it, a 30 mph limit means a
true stall of about 24.5 mph, and these designs fail it.

**Don't compare the TOTAL MISSION SCORE line.** A single design is scored
against built-in guesses of the field's best (0.060 lb/s, 22 lap·lb,
120 lb·in), which every one of these beats, so they all print 7.500.
Compare the raw values in the SCORE section instead, or turn them into the
proposal figure's score of 7.5:

```
score = 1 + (1 + M2raw / 0.13311) + (2 + M3raw / 157.81) + (0.5 + GMraw / 2132.76)
```

For `stall30mph.json`, the raw values are 0.0661 lb/s, 92.54 lap·lb and
1110.4 lb·in, which gives 6.10.

## How the sweep was run

One Slurm job per limit on Sapelo2 (`cluster/slurm/jobs_nochordcap_v6.txt`):

```bash
sbatch --array=1-8 --mem=48G --time=24:00:00 slurm/array.sbatch slurm/jobs_nochordcap_v6.txt
```

Each job runs the search with these arguments (here, the 30 mph line):

```
--thorough --top 25 --finalists 32 --motor-counts 1 2 --keep-propulsion 600 --refine 100
--swap-parents 12 --seeds seeds --ar 1 14 --max-stall-mph 30 --stall-lift-fos 1.0
--fill-cores --stall-air tucson-hot
```

The jobs took 1 to 10 hours each on a full cluster node, so don't rerun the
search on a PC to check these; score the design files instead.

### Settings changed from the defaults

These change which designs are allowed:

| Flag | Used | Default | What it does |
|---|---|---|---|
| `--max-stall-mph` | 30-50 | no limit | the stall limit |
| `--stall-lift-fos` | 1.0 | 1.5 | true stall instead of stall on 2/3 of the lift |
| `--stall-air` | tucson-hot | tucson-hot | the 36 C air the limit is checked in (stated, not changed) |
| `--motor-counts` | 1 2 | 1-4 | only 1- and 2-motor designs searched |
| `--ar` | 1 14 | 1.5 14 | aspect ratios searched |
| `--max-chord-in` | not used | no cap | v6 has no chord cap (the earlier v5 sweep used 24 in) |

The rest only set how wide the search is, not how a design is scored:
`--thorough`, `--top`, `--finalists`, `--keep-propulsion`, `--refine`,
`--swap-parents`, `--fill-cores`, and `--seeds` (the cluster bundle's folder
of past winners).

### Defaults worth knowing

These weren't changed, but they set the numbers above:

- **Air:** Tucson at 667 m. The nominal day is 28 C with a 5 m/s wind. Every
  pick also has to fly the 12 C morning and the 36 C afternoon (the weather
  presets), and may change its propeller between them (rule 3.2.1e).
- **Battery:** each flight turns for landing at 88% of the pack used, and
  the energy limit is 100 Wh with a 100 A fuse.
- **Laps:** lap times are divided by 0.95 for the pilot, and turns are
  banked no more than 70 degrees.
- **Motors:** the winding resistance rises with temperature, which is
  tracked through each flight.
