"""DBF 2026-27 aircraft sizing, mission simulation and design optimiser.

Quick start
-----------
    from dbf2027 import Config, Design, Aircraft, evaluate_scores

    cfg = Config()
    d = Design(span_m=1.55, chord_m=0.27, airfoil="SD7062",
               motor="COBRA_C4120_12", prop="APC_12x8", n_motors=1,
               pack="GA_6S3300_60C", n_systems=1,
               sensor_mass_kg=0.9, n_simulators=3)
    print(evaluate_scores(Aircraft(d, cfg), cfg).total_mission)
"""
from .aerodynamics import AeroModel, BodyGeometry, WingGeometry
from .aircraft import Aircraft, Design
from .components import (AIRFOILS, ESCS, MOTORS, PACKS, PROPS, BatteryPack,
                         legal_packs, load_motors_csv, load_packs_csv,
                         load_props_csv)
from .config import BuildStandard, Config, Course, Environment, Rules
from .mission import (MissionResult, mission_1, mission_2, mission_3,
                      turn_capability)
from .optimizer import (DesignSpace, SearchSettings, evaluate_design,
                        optimise, sensitivity)
from .plots import Chart, Mapper, Series, dimensions, three_view
from .scoring import ScoreCard, evaluate_scores, rescore

__version__ = "1.10.0"

__all__ = [
    "AeroModel", "BodyGeometry", "WingGeometry",
    "Aircraft", "Design",
    "AIRFOILS", "ESCS", "MOTORS", "PACKS", "PROPS", "BatteryPack",
    "legal_packs", "load_motors_csv", "load_packs_csv", "load_props_csv",
    "BuildStandard", "Config", "Course", "Environment", "Rules",
    "MissionResult", "mission_1", "mission_2", "mission_3", "turn_capability",
    "DesignSpace", "SearchSettings", "evaluate_design", "optimise",
    "sensitivity",
    "ScoreCard", "evaluate_scores", "rescore",
    "Chart", "Mapper", "Series", "dimensions", "three_view",
]

# `catalogue` and `gui` are deliberately NOT imported here: catalogue
# walks this package (importing it from here would be circular) and gui
# pulls in tkinter, which not every environment has.  Import them
# directly - `from dbf2027 import gui` - when you want them.
