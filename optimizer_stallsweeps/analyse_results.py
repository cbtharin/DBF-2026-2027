#!/usr/bin/env python3
"""Audit a finished search: did every option get a fair hearing?

A catalogue is only useful if the search actually explores it.  This reads
`ranked_designs.csv` and reports which airfoils, motors, propellers, motor
counts and chemistries survived to the finals, how they rank, and how much
of the spread is real rather than noise.

    python analyse_results.py results_v4/ranked_designs.csv
"""
from __future__ import annotations

import collections
import csv
import sys

from dbf2027.components import AIRFOILS, MOTORS, PACKS


def load(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def f(row, key, default=0.0):
    try:
        return float(row.get(key, default) or default)
    except ValueError:
        return default


def tally(rows, key, label, universe=None):
    print()
    print("%s" % label)
    print("-" * len(label))
    counts = collections.Counter(r[key] for r in rows)
    best = {}
    for r in rows:
        v = f(r, "objective")
        if r[key] not in best or v > best[r[key]]:
            best[r[key]] = v
    order = sorted(counts, key=lambda k: -best.get(k, -1e9))
    print("  %-26s %7s %10s" % ("option", "finals", "best score"))
    for k in order:
        print("  %-26s %7d %10.3f" % (k[:26], counts[k], best[k]))
    if universe:
        missing = [k for k in universe if k not in counts]
        if missing:
            print("  did not reach the finals: %s" % ", ".join(
                sorted(m[:20] for m in missing)))


def main(argv):
    path = argv[1] if len(argv) > 1 else "results_v4/ranked_designs.csv"
    rows = load(path)
    if not rows:
        print("no rows in %s" % path)
        return 1

    print("=" * 70)
    print("SEARCH AUDIT: %s  (%d finalists)" % (path, len(rows)))
    print("=" * 70)

    objs = sorted((f(r, "objective") for r in rows), reverse=True)
    print()
    print("Score spread")
    print("------------")
    print("  best            %.4f" % objs[0])
    print("  10th            %.4f" % objs[min(9, len(objs) - 1)])
    print("  worst finalist  %.4f" % objs[-1])
    within = sum(1 for o in objs if o >= objs[0] * 0.98)
    print("  within 2%% of the best: %d of %d designs" % (within, len(objs)))
    if within > len(objs) // 3:
        print("  -> the objective is FLAT. Pick on buildability, not score.")

    tally(rows, "d_airfoil", "Airfoils in the finals", AIRFOILS.keys())
    tally(rows, "motor", "Motors in the finals")
    tally(rows, "d_prop", "Propellers in the finals")
    tally(rows, "d_n_motors", "Motor count")
    tally(rows, "chemistry", "Battery chemistry")

    print()
    print("Design-variable ranges across the finals")
    print("----------------------------------------")
    for key, label, scale, unit in (
            ("d_span_m", "span", 39.3701, "in"),
            ("AR", "aspect ratio", 1.0, ""),
            ("d_sensor_mass_kg", "sensor", 2.20462, "lb"),
            ("d_n_simulators", "simulators", 1.0, ""),
            ("d_drop_height_in", "drop height", 1.0, "in"),
            ("total_Wh", "battery", 1.0, "Wh"),
            ("mtow_m2_kg", "M2 gross", 2.20462, "lb"),
            ("n_ultimate", "ultimate load", 1.0, "g")):
        vals = [f(r, key) for r in rows if r.get(key)]
        if not vals:
            continue
        print("  %-16s %8.2f to %8.2f %s"
              % (label, min(vals) * scale, max(vals) * scale, unit))

    print()
    print("What limited the turn, and what sized the spar")
    print("----------------------------------------------")
    for key, label in (("m2_turn_limit", "Mission 2 turn limited by"),
                       ("load_driver", "ultimate load governed by")):
        c = collections.Counter(r.get(key, "?") for r in rows)
        print("  %s:" % label)
        for k, n in c.most_common():
            print("      %-26s %d" % (k, n))
    gauge = collections.Counter(r.get("spar_gauge_limited", "?") for r in rows)
    print("  spar gauge limited: %s" % dict(gauge))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
