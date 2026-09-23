"""The 2026-27 rules, and where each one is enforced in this code.

Why this exists
---------------
A comment saying "rule 3.1.3 says simulators must match the container" is
worth something the day it is written and nothing two refactors later,
because nothing checks it.  So the rules live here as data, the code is
annotated with `RULE <id>` markers, and `check_annotations()` reads the
source back and verifies the two agree:

  * every rule that this model is supposed to enforce has at least one
    marker in the source;
  * every marker in the source names a rule that actually exists;
  * every rule cites the section of the rules package it came from.

`python dbf.py --rules` prints the traceability matrix.  `--validate`
fails the build if it has drifted.  That turns "we documented the rules"
into something that stays true.

How to add a rule
-----------------
1. Add a `Rule(...)` below with its section number and a short quote.
2. Put `# RULE <id>` on the line or block that enforces it.
3. Run `python dbf.py --rules` and confirm it is no longer listed as
   unenforced.

Kinds
-----
``model``    changes a computed number (mass, drag, time, energy)
``limit``    a hard constraint; a design that breaks it is rejected
``score``    part of a scoring formula
``build``    a physical or procedural requirement that no simulator can
             check - it belongs on a build and tech-inspection checklist,
             and it is listed here so it is not silently forgotten

A ``build`` rule needs no source marker.  Everything else does.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# Ids carry an optional letter suffix (3.1.3a) because one numbered
# section often states several separable requirements.  The marker regex
# has to accept that suffix, or "RULE 3.1.3a" silently registers as a
# marker for the non-existent rule "3.1.3".
_ID = r"(?:\d+\.)*\d+[a-z]?"
RULE_MARKER = re.compile(r"RULE\s+(" + _ID + r"(?:\s*,\s*" + _ID + r")*)")


@dataclass(frozen=True)
class Rule:
    """One requirement from the rules package."""
    id: str                      # section number, e.g. "3.2.3.1"
    title: str
    text: str                    # short quote or faithful paraphrase
    kind: str                    # model | limit | score | build
    note: str = ""               # how this model handles it, or why not


RULES: Tuple[Rule, ...] = (
    # ==================================================================
    # 3.1  Mission and payload requirements
    # ==================================================================
    Rule("3.1", "Maximum wingspan",
         "Maximum allowable wingspan is 6 feet.",
         "limit",
         "Clamped in the search bounds AND rejected in static_checks, so "
         "it binds whether a design is sampled or hand-entered."),
    Rule("3.1.1a", "Sensor is towed and recoverable",
         "The sensor is a towed device that must be deployable from and "
         "recoverable by the airplane in flight.",
         "model",
         "Mission 3 flies the deployed drag of the sensor plus its tow "
         "line; the sensor is stowed for takeoff and landing."),
    Rule("3.1.1b", "Minimum sensor length",
         "The sensor must have a minimum overall length of six inches.",
         "limit",
         "The envelope is derived from sensor mass and packing density, "
         "so a very light sensor comes out too short and is rejected."),
    Rule("3.1.1c", "Mission 3 sensor weight",
         "The sensor weight may be changed for any Mission 3 attempt; the "
         "weight may not exceed the maximum declared weight in Tech "
         "Inspection.",
         "model",
         "`sensor_m3_kg` is a design variable bounded by the declared "
         "maximum. Lighter buys laps; M3 scores laps x weight, so there "
         "is a real optimum in between."),
    Rule("3.1.1d", "Sensor carried internally",
         "The sensor must be carried internally to the airplane.",
         "model",
         "The sensor volume is packed into the fuselage bay alongside "
         "the containers, which is what sets fuselage length."),
    Rule("3.1.1e", "Sensor lights and battery",
         "Three lights with OFF / SOLID ON / FLASHING ON modes, commanded "
         "over the tow line; control electronics and a battery internal "
         "to the sensor.",
         "build",
         "Mass is inside the sensor mass the team declares. The light "
         "controller and the three modes are named in the parts list; "
         "nothing here can verify they work."),
    Rule("3.1.2", "Tow line length",
         "The deployment and recovery mechanism must deploy the sensor "
         "with a tow line to a distance equal to or greater than 1.5 "
         "times the airplane wingspan.",
         "model",
         "Charged as drag: line length scales with span, so a bigger "
         "wing is penalised twice in Mission 3."),
    Rule("3.1.2b", "Deployment mechanism internal and permanent",
         "The deployment and recovery mechanism must be internal to the "
         "airplane and permanently installed in the Mission 3 location "
         "for all missions.",
         "model",
         "Its mass is carried in every mission, including Mission 1."),
    Rule("3.1.3a", "Simulator weight must match the loaded container",
         "The shipping container simulator must be the same weight as the "
         "sensor shipping container plus the maximum weight sensor "
         "declared in Tech Inspection within +/- 1 ounce.",
         "model",
         "NOT a design variable. M2 scored weight is (1 + n) loaded "
         "containers, so each extra kg of sensor is carried (1 + n) "
         "times. Treating it as free let the optimiser fly a heavy "
         "sensor beside a stack of light boxes, which is illegal."),
    Rule("3.1.3b", "Simulator size must match",
         "The shipping container simulators must be the same shape and "
         "size as the sensor shipping container within +/- 1/8 inch.",
         "model",
         "Every box uses one container envelope for bay packing and drag."),
    Rule("3.1.3c", "Container is a rectangular prism that encloses the sensor",
         "The sensor shipping container must fully enclose the sensor and "
         "must be a rectangular prism.",
         "limit",
         "A sensor that will not fit the box at a physical density is "
         "rejected."),
    Rule("3.1.3d", "Rigging and container depth",
         "The total distance from the rigging ring to the bottom surface "
         "of any container side should be less than 15 inches to assure "
         "teams will have the ability to drop up to the maximum drop "
         "height.",
         "limit",
         "Checked in the Ground Mission: too deep a container cannot "
         "reach the declared drop height."),
    Rule("3.1.3e", "Container markings",
         "TOP marked on the top surface; each pair of opposing sides "
         "marked O or E.",
         "build",
         "Named in the parts list. Nothing here can verify a marking."),

    # ==================================================================
    # 3.2  Airplane requirements
    # ==================================================================
    Rule("3.2.1a", "Gross weight limit",
         "The airplane TOGW (take-off gross weight with payload) must be "
         "less than 55-lb.",
         "limit",
         "Checked at Mission 2 weight, which is the heaviest case, and "
         "it also bounds how much payload the search bothers sampling."),
    Rule("3.2.1b", "Ground rolling takeoff, no distance limit",
         "The airplane will use ground rolling takeoff and landing. There "
         "is no limit on take-off distance this year.",
         "model",
         "The roll is integrated and reported because a 300 ft roll is a "
         "real operational risk, but it is not a constraint."),
    Rule("3.2.1c", "No externally assisted takeoff",
         "All energy for take-off must come from the on-board propulsion "
         "battery pack(s).",
         "model",
         "Takeoff energy is drawn from the same pack budget as the laps."),
    Rule("3.2.1d", "COTS electric motor, commercial propeller",
         "Unmodified over-the-counter brushed or brushless electric "
         "motor; commercially produced propeller, which may be clipped.",
         "build",
         "The catalogues hold real parts; every entry carries a `source` "
         "field marked VERIFY."),
    Rule("3.2.1e", "Propeller may change between flights",
         "The propeller (diameter/pitch) may be changed for each flight "
         "attempt.",
         "build",
         "NOT exploited: this model flies one propeller for all missions, "
         "which is conservative. Per-mission propellers would score "
         "higher and are legal."),
    Rule("3.2.1f", "No autopilot, GPS or cameras",
         "No autopilots/flight controllers, no onboard GPS, no cameras.",
         "build",
         "Nothing in the model implies any of them."),
    Rule("3.2.2", "Arming plug",
         "An external motor arming plug for each propulsion system, "
         "separate from the Rx switch, on the top surface, at least 6 "
         "inches from any propeller plane, between the fuse and the ESC.",
         "build",
         "One per propulsion system is costed in the mass buildup and "
         "named in the parts list; placement is a build item."),
    Rule("3.2.3a", "COTS battery packs only",
         "All battery packs must be un-altered and commercially procured "
         "as COTS battery packs. Custom battery packs are not allowed.",
         "model",
         "Packs are atomic catalogue entries. There is no way to express "
         "cells in series or parallel, because that would be a custom "
         "pack."),
    Rule("3.2.3b", "Separate receiver battery, BEC disabled",
         "A separate battery is required for the Rx/Servos. If the ESC "
         "has a BEC it must be disabled.",
         "model",
         "Its mass is carried; it does not count against the 100 Wh "
         "propulsion budget."),
    Rule("3.2.3c", "Charging sack and insulated connectors",
         "Fully insulated connectors; batteries stored and charged in a "
         "commercial labelled charging sack.",
         "build",
         "Named in the parts list with the rule attached."),
    Rule("3.2.3.1a", "One pack per propulsion system",
         "There shall be a maximum of one battery pack connected to a "
         "propulsion system.",
         "model",
         "`n_systems` is a design variable; each system carries exactly "
         "one pack."),
    Rule("3.2.3.1b", "No series or parallel pack wiring",
         "Each battery pack must be independently connected to its own "
         "propulsion system. Batteries may not be connected in series or "
         "parallel.",
         "model",
         "The only way to carry a second pack is a second complete "
         "system, which also costs a fuse and an arming plug."),
    Rule("3.2.3.1c", "Identical packs across systems",
         "All propulsion battery packs must be identical (same "
         "manufacturer, part number, size, voltage, power, rating).",
         "model",
         "One pack key for the whole aeroplane, so they cannot differ."),
    Rule("3.2.3.1d", "Fuse per system, blade style",
         "Each battery/propulsion system is required to have its own fuse "
         "and arming plug. The fuse must be a blade style fuse in a "
         "separate harness between the battery and arming plug.",
         "model",
         "Mass costed per system; the fuse rating is checked against the "
         "current the system actually draws."),
    Rule("3.2.3.2a", "Battery chemistry",
         "Propulsion batteries may be NiCAD, NiMH or lithium-based "
         "chemistries.",
         "limit",
         "Only LiPo and Li-ion are catalogued, at the team's direction: "
         "nickel loses on both mass and current at these power levels."),
    Rule("3.2.3.2b", "Total propulsion energy",
         "Total Propulsion Energy stored on airplane (sum of all "
         "propulsion batteries) shall not exceed 100 Watt-hours.",
         "limit",
         "Summed over every pack, checked before anything is flown."),
    Rule("3.2.3.2c", "Per-pack energy limit",
         "Individual battery packs may not exceed the FAA limits for hand "
         "carry of 100 Watt-hours (rated capacity * rated voltage) per "
         "battery pack.",
         "limit",
         "Checked per pack as well as on the sum - a second system of "
         "legal packs can still break the airframe total."),
    Rule("3.2.3.2d", "Fuse current rating",
         "The maximum current rating for the fuse is the maximum "
         "continuous discharge current rating of the battery pack "
         "(capacity * C-rating) up to 100 amps.",
         "limit",
         "The binding electrical constraint in practice. Checked at "
         "static full throttle, which is where it bites."),

    # ==================================================================
    # 3.3  Mission operations
    # ==================================================================
    Rule("3.3.1", "Flight course",
         "Upwind and downwind markers will be 500 ft from the starting "
         "line.",
         "model",
         "The start/finish line is in the MIDDLE: 500 out, 180, 1000 "
         "back with a 360 at midfield, 180, 500 home. Four turns to "
         "accelerate out of per lap, not three."),
    Rule("3.3.2", "Staging",
         "Five minute assembly limit in the staging box; batteries and "
         "payloads installed there; no work afterwards.",
         "build",
         "An operational constraint on the team, not on the aeroplane. "
         "It does argue for a payload bay you can load fast."),
    Rule("3.3.3a", "Missions flown in order",
         "The Flight Missions must be flown in order. A new mission may "
         "not be flown until the team has obtained a successful score for "
         "the preceding mission.",
         "score",
         "Scored as all-or-nothing per mission; a failed mission scores "
         "0, not a small number."),
    Rule("3.3.3b", "Mission 1 may be waived",
         "Teams may elect to waive a Mission 1 attempt and proceed "
         "directly to Mission 2. A successful Mission 2 then also results "
         "in a Mission 1 score.",
         "build",
         "A contest-day tactic, not a design variable. The model always "
         "flies M1, which is the conservative reading."),
    Rule("3.3.3c", "Same configuration for all missions",
         "The airplane must be flown in the same configuration for all "
         "three missions unless specifically cited herein.",
         "model",
         "One airframe, one propeller, one motor count across M1, M2 and "
         "M3; payload and the M3 sensor weight change, as the rules "
         "allow. NOT modelled: the FAQ also permits swapping the "
         "BATTERY between missions ('unless specifically allowed ... "
         "i.e., batteries, payloads, payload support'), and this model "
         "flies one pack throughout. Measured cost of that on the "
         "current optimum: nil - both scored missions want the full "
         "100 Wh, and the only mission that could use a lighter pack is "
         "M1, which scores 1.0 either way."),
    Rule("3.3.3d", "Landing required, and outside the window",
         "A successful landing must be completed to get a score. Landing "
         "is not part of the 5 minute time window.",
         "model",
         "Landing is simulated, its energy is charged to the pack, and "
         "it decides feasibility - but its time is excluded from the "
         "scored window."),
    Rule("3.3.3.1", "Mission 1 - Staging Flight",
         "No payload. 3 laps within a 5-minute window. M1 = 1.0 for a "
         "successful mission.",
         "score",
         "Pass/fail, so it cannot be optimised - only achieved."),
    Rule("3.3.3.2", "Mission 2 - Delivery Flight",
         "Payload is the sensor in shipping container plus optional "
         "simulators. 5 laps in a 5-minute window. "
         "M2 = 1 + N(weight/time) / Max(weight/time).",
         "score",
         "Weight is the scored payload only, not gross weight. Time "
         "stops crossing the line on the fifth lap."),
    Rule("3.3.3.3", "Mission 3 - Sensor Flight",
         "Sensor deployed, light modes commanded through the flight, "
         "stowed before landing. 5-minute window. "
         "M3 = 2 + N(#laps x sensor weight) / Max(...).",
         "score",
         "As many laps as fit the window, at the M3 sensor weight."),
    Rule("3.3.4", "Ground Mission",
         "Six drops onto concrete, two on each of three faces, from a "
         "height declared in whole inches up to 60. The sensor must be at "
         "its declared maximum weight. Any physical damage is a failure. "
         "GM = 0.5 + N(weight x height) / Max(weight x height).",
         "score",
         "Not a flight, so there is nothing to integrate: it turns on "
         "whether the container foam can hold the sensor under the "
         "deceleration it survives. This is the THIRD place sensor "
         "weight pays."),

    # ==================================================================
    # 4  Tech inspection
    # ==================================================================
    Rule("4.2.4", "Wing tip load test",
         "The airplane will be lifted at each wing tip at MGTOW to verify "
         "adequate wing strength (roughly equivalent to a 2.5 g load "
         "case). Verify the MGTOW does not exceed 55 lbs.",
         "model",
         "The spar is sized from a derived V-n ultimate near 5 g, which "
         "is well above the 2.5 g the test demands. Deliberately "
         "conservative: the test is a floor, not the design case."),

    # ==================================================================
    # 6  Scoring
    # ==================================================================
    Rule("6.3", "Total score",
         "Competition Score = Total Report Score * Total Mission Score + "
         "Participation Score. Total Mission Score = M1 + M2 + M3 + GM.",
         "score",
         "Mission score is a MULTIPLIER on the report score, so a point "
         "of mission score is worth proportionally more the better the "
         "report is. Ceiling is 7.5."),
    Rule("6.2", "Units and rounding",
         "Scoring uses US English units. Time to 2 dp (s), length to 2 dp "
         "(in), weight to 2 dp (oz and lb).",
         "build",
         "The model works in SI internally and reports both. It does not "
         "round to the scoring precision, because rounding a comparison "
         "between candidates would only lose information."),
)


RULES_BY_ID: Dict[str, Rule] = {r.id: r for r in RULES}

# A rule of these kinds must be traceable to a line of source.
ENFORCED_KINDS = ("model", "limit", "score")


# ==========================================================================
# Reading the annotations back out of the source
# ==========================================================================
def _source_files(root: Optional[str] = None) -> List[str]:
    here = os.path.dirname(os.path.abspath(__file__))
    root = root or os.path.dirname(here)
    out = []
    for folder in (here, root):
        for name in sorted(os.listdir(folder)):
            if name.endswith(".py"):
                out.append(os.path.join(folder, name))
    return out


def scan_annotations(root: Optional[str] = None) -> Dict[str, List[str]]:
    """Map rule id -> ["file.py:line", ...] for every `RULE <id>` marker.

    This module is skipped: it is the registry, not an enforcement site,
    and every id appears in it by definition.
    """
    found: Dict[str, List[str]] = {}
    me = os.path.abspath(__file__)
    for path in _source_files(root):
        if os.path.abspath(path) == me:
            continue
        try:
            with open(path, encoding="utf-8") as fh:
                lines = fh.readlines()
        except OSError:
            continue
        short = os.path.basename(path)
        for n, line in enumerate(lines, 1):
            for m in RULE_MARKER.finditer(line):
                for rid in re.split(r"\s*,\s*", m.group(1)):
                    found.setdefault(rid, []).append("%s:%d" % (short, n))
    return found


def check_annotations(root: Optional[str] = None) -> Tuple[List[str], List[str]]:
    """Return (unenforced rule ids, marker ids that match no rule).

    Both lists empty means the code and the rulebook agree.
    """
    found = scan_annotations(root)
    unenforced = [r.id for r in RULES
                  if r.kind in ENFORCED_KINDS and r.id not in found]
    unknown = sorted(set(found) - set(RULES_BY_ID))
    return unenforced, unknown


def coverage(root: Optional[str] = None) -> Dict[str, Dict]:
    """Every rule with its enforcement sites, for a report or a UI."""
    found = scan_annotations(root)
    return {r.id: {"rule": r, "sites": found.get(r.id, [])} for r in RULES}


def _wrap(text: str, width: int, indent: str) -> List[str]:
    words, lines, cur = text.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > width:
            lines.append(indent + cur)
            cur = w
        else:
            cur = (cur + " " + w) if cur else w
    if cur:
        lines.append(indent + cur)
    return lines


def as_text(width: int = 78, verbose: bool = True,
            root: Optional[str] = None) -> str:
    """The traceability matrix, for a terminal."""
    cov = coverage(root)
    unenforced, unknown = check_annotations(root)
    out = ["=" * width,
           "DBF 2026-27 RULES - WHAT THE MODEL ENFORCES, AND WHERE",
           "=" * width, "",
           "kind:  model = changes a number   limit = rejects a design",
           "       score = a scoring formula  build = checklist item, not",
           "                                          checkable here", ""]
    by_kind: Dict[str, List[Rule]] = {}
    for r in RULES:
        by_kind.setdefault(r.kind, []).append(r)

    for kind in ("limit", "model", "score", "build"):
        group = by_kind.get(kind, [])
        if not group:
            continue
        out.append("-" * width)
        out.append("%s  (%d)" % (kind.upper(), len(group)))
        out.append("-" * width)
        for r in group:
            sites = cov[r.id]["sites"]
            flag = "" if (sites or kind == "build") else "   << NOT ENFORCED"
            out.append("  %-10s %s%s" % (r.id, r.title, flag))
            if verbose:
                out += _wrap('"' + r.text + '"', width - 6, "      ")
                if r.note:
                    out += _wrap("-> " + r.note, width - 6, "      ")
            if sites:
                shown = ", ".join(sites[:4])
                if len(sites) > 4:
                    shown += ", +%d more" % (len(sites) - 4)
                out.append("      enforced at: %s" % shown)
            elif kind == "build":
                out.append("      (build / tech-inspection checklist)")
            out.append("")

    out.append("=" * width)
    n_sites = sum(len(v["sites"]) for v in cov.values())
    out.append("%d rules registered, %d annotation sites in the source"
               % (len(RULES), n_sites))
    if unenforced:
        out.append("UNENFORCED (no `RULE <id>` marker found): "
                   + ", ".join(unenforced))
    if unknown:
        out.append("UNKNOWN markers (no such rule): " + ", ".join(unknown))
    if not unenforced and not unknown:
        out.append("Annotations and rulebook agree.")
    out.append("=" * width)
    return "\n".join(out)


def as_markdown(root: Optional[str] = None) -> str:
    """The same matrix, for committing next to the README."""
    cov = coverage(root)
    unenforced, unknown = check_annotations(root)
    out = ["# Rules compliance matrix", "",
           "Every requirement from the 2026-27 rules package that bears on "
           "this model, and the line of source that enforces it. Generated "
           "by `python dbf.py --rules --markdown`; `--validate` fails if "
           "an annotation drifts from the registry in "
           "[`dbf2027/rulebook.py`](dbf2027/rulebook.py).", "",
           "| kind | meaning |", "|---|---|",
           "| `model` | changes a computed number |",
           "| `limit` | a hard constraint; a design that breaks it is rejected |",
           "| `score` | part of a scoring formula |",
           "| `build` | physical or procedural; no simulator can check it, "
           "listed so it is not forgotten |", ""]
    for kind, label in (("limit", "Hard limits"), ("model", "Modelled"),
                        ("score", "Scoring"),
                        ("build", "Build and tech inspection")):
        group = [r for r in RULES if r.kind == kind]
        if not group:
            continue
        out += ["## %s" % label, "",
                "| Rule | Requirement | How this model handles it | Enforced at |",
                "|---|---|---|---|"]
        for r in group:
            sites = cov[r.id]["sites"]
            where = ("`" + "`, `".join(sites[:3]) + "`") if sites else "—"
            if len(sites) > 3:
                where += " +%d" % (len(sites) - 3)
            out.append("| **%s**<br>%s | %s | %s | %s |"
                       % (r.id, r.title, r.text.replace("|", "\\|"),
                          r.note.replace("|", "\\|"), where))
        out.append("")
    out.append("---")
    out.append("")
    if unenforced or unknown:
        if unenforced:
            out.append("**UNENFORCED:** " + ", ".join(unenforced))
        if unknown:
            out.append("**UNKNOWN MARKERS:** " + ", ".join(unknown))
    else:
        out.append("All enforceable rules are annotated in the source, and "
                   "every annotation names a registered rule.")
    return "\n".join(out) + "\n"


def main(argv: Optional[List[str]] = None) -> int:
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--markdown" in argv:
        with open("RULES.md", "w", encoding="utf-8") as fh:
            fh.write(as_markdown())
        unenforced, unknown = check_annotations()
        print("wrote RULES.md  (%d rules)" % len(RULES))
        if unenforced or unknown:
            print("  WARNING: unenforced %s, unknown %s"
                  % (unenforced, unknown))
        return 0
    print(as_text(verbose="--brief" not in argv))
    unenforced, unknown = check_annotations()
    return 1 if (unenforced or unknown) else 0


if __name__ == "__main__":
    raise SystemExit(main())
