# DBF 2026-27 optimiser — the 20-file core

These twenty files are the whole program. No `pip install`, no launcher
script, no test files, no example designs, no other folders.

## Which one starts it

**`gui.py`.** You never run it directly, though — you run the package that
contains it:

```bash
python -m dbf2027.gui
```

That opens the desktop application: search space, component catalogues,
charts, three-view, dimensions, score breakdown, rules traceability, and it
runs its own optimiser searches. Every other file in the list is a module it
imports. There is nothing else to open.

`gui.py` cannot be run on its own. It uses relative imports (`from . import
plots`), so `python gui.py` fails with:

```
ImportError: attempted relative import with no known parent package
```

That is expected, not a broken file.

## The folder layout is not optional

The twenty files must sit inside a folder named **exactly** `dbf2027`, and
you run the command from the folder *above* it:

```
dbf-optimizer/            <-- run "python -m dbf2027.gui" from HERE
└── dbf2027/              <-- this name must be exactly "dbf2027"
    ├── __init__.py
    ├── gui.py
    └── ... the other 18
```

Two ways this goes wrong, both of which produce the same unhelpful error:

| what you did | what happens |
|---|---|
| all 20 files loose in one folder | `ModuleNotFoundError: No module named 'dbf2027'` |
| ran the command from *inside* `dbf2027/` | `ModuleNotFoundError: No module named 'dbf2027'` |

If you see that message, the files are almost certainly fine and the layout
is wrong. Check that `dbf2027` is the folder name and that you are one level
above it when you run the command.

## Requirements

- **Python 3.8+** (verified on 3.13)
- **tkinter**, which ships with the standard CPython installer on Windows and
  macOS. On Linux it is a separate package — `sudo apt install python3-tk`.
- Nothing else. Pure standard library: no numpy, no matplotlib, no scipy.

## The 20 files

### Entry point

| file | lines | what it is |
|---|---|---|
| `gui.py` | 1537 | Desktop front end: change a number, see what it does to the aeroplane. **This is what starts the program.** |
| `__init__.py` | 45 | Package exports. Imports `plots`, which is why `plots.py` is never optional. |

### Configuration and units

| file | lines | what it is |
|---|---|---|
| `config.py` | 434 | Environment, course, rules and build-standard configuration. Every number you might want to change lives here. |
| `units.py` | 46 | Unit conversions. All internal computation is strictly SI. |

### Physics

| file | lines | what it is |
|---|---|---|
| `aerodynamics.py` | 371 | Aerodynamic model: lift, drag and the drag polar. |
| `liftingline.py` | 104 | Prandtl lifting-line: span efficiency from the actual planform. |
| `propulsion.py` | 301 | Motor, propeller, ESC and battery matched at every flight point. |
| `mission.py` | 740 | Mission simulation: takeoff roll, climb, laps, energy budget. |
| `loads.py` | 110 | Design load factors, derived rather than assumed. |
| `weights.py` | 140 | Mass buildup. |
| `payload.py` | 145 | Payload geometry: sensor envelope, shipping container, bay packing. |

### The aeroplane, and choosing one

| file | lines | what it is |
|---|---|---|
| `aircraft.py` | 499 | The `Design` vector and the `Aircraft` it produces. |
| `components.py` | 761 | Component catalogues: airfoils, motors, propellers, cells, ESCs. |
| `scoring.py` | 481 | Mission scoring per the 2026-27 rules. |
| `optimizer.py` | 741 | Staged search over the design space. |

### Numerics

| file | lines | what it is |
|---|---|---|
| `numerics.py` | 435 | Dependency-free numerical methods: RK4, Brent, Gauss–Legendre, PCHIP, golden section, Latin hypercube. |

### Output

| file | lines | what it is |
|---|---|---|
| `report.py` | 656 | Human-readable output: rankings, a full design dossier, a shopping list. |
| `plots.py` | 754 | Chart data and the arithmetic that turns it into pixels. |
| `catalogue.py` | 273 | Function catalogue: what every callable in the package is, and where. |
| `rulebook.py` | 569 | The 2026-27 rules, and where each one is enforced in this code. |

**9,142 lines total.**

## Why three of these are in the set

Not every file earns its place for an obvious reason:

- **`plots.py`** cannot be removed even though nothing about the physics needs
  charts — `__init__.py` imports it, so importing the package pulls it in.
- **`catalogue.py`** and **`rulebook.py`** are here only because `gui.py`
  builds its Catalogues and Rules tabs from them.

Drop the GUI and all three become removable, leaving a 17-module headless
solver — but you would then need a command-line entry point, which is not in
this set.

## Default flight conditions

The model defaults to the contest site, not sea-level ISA, because designing
at sea level and discovering 4,000 ft of density altitude later is the wrong
order to do things in.

| preset | elevation | temp | wind | density | density altitude |
|---|---|---|---|---|---|
| `tucson-morning` | 667 m | 12 °C | 3.0 m/s | 1.143 | 2,347 ft |
| **`tucson`** (default) | 667 m | 28 °C | 5.0 m/s | **1.082** | **4,171 ft** |
| `tucson-hot` | 667 m | 36 °C | 7.0 m/s | 1.054 | 5,038 ft |
| `isa` | 0 m | 15 °C | 4.5 m/s | 1.225 | 0 ft |

Elevation is the USGS 3DEP value at 32.2647 N, 111.2740 W — 667.06 m,
2,188 ft. Temperatures are NWS Tucson climate normals for April 1992–2021:
mean daily high 83 °F, mean daily low 53 °F. These are normals, not a
forecast. All of it is editable in the GUI, and all of it lives in
`config.py`.

## What this set does not include

Everything below is absent, and the GUI works without all of it. Each entry
is the file you would need to add.

| absent | what you lose |
|---|---|
| `run_optimizer.py` | the command-line search — `--quick`, `--thorough`, `--evaluate` |
| `dbf.py` | the dispatcher that fronts every mode with one command |
| `webapp.py`, `webpage.py` | the browser front end |
| `selftest.py`, `validate.py`, `guitest.py`, `webtest.py` | the four test suites, 1,329 checks |
| `grid_convergence.py`, `analyse_results.py` | numerical error study, assumptions audit |
| `examples/*.json` | the saved designs the GUI loads at startup |

Two of those are worth knowing about specifically:

- **You have not lost the ability to search.** `gui.py` calls `optimizer.py`
  directly and never touches `run_optimizer.py`. The Search tab works.
- **Without `examples/`**, the GUI opens on a built-in baseline design
  instead of a saved one, and starts normally. Add
  `examples/optimised_v1.json` next to the `dbf2027` folder if you want it to
  open on the optimised design instead.
