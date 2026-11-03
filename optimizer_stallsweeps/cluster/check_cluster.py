#!/usr/bin/env python3
"""Check a cluster install against the PC, and estimate search times.

    python check_cluster.py            # on the cluster (in a job: see slurm/check.sbatch)

1. Same answers: scores a fixed set of designs (the full search's 1- and
   2-motor winners and perturbations of them, some of which fail) and runs
   two short refinements, then compares every number with the PC's
   (reference/check_reference.json, written when the bundle was built).
   "identical" means bit for bit.  Differences in the last digits come from
   the platform's maths library or a Python older than 3.12 (which sums
   floats differently) - harmless for any one score, but a refinement can
   then take a slightly different path, so searches are close rather than
   identical.
2. Speed: times the designs with every allotted core busy and scales the
   PC's full-search stage times (exp8: 9.9 h on 15 workers) to this node.

    python check_cluster.py --make-reference   # on the PC, by make_cluster_bundle.py
"""
from __future__ import annotations

import json
import math
import os
import platform
import sys
import time
from dataclasses import asdict, replace
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import hpc  # noqa: E402
from dbf2027.aircraft import Design  # noqa: E402
from dbf2027.config import Config  # noqa: E402
from dbf2027.optimizer import (DesignSpace, SearchSettings, _refine_one,  # noqa: E402
                               evaluate_design)

REF_FILE = os.path.join(HERE, "reference", "check_reference.json")
DX = 12.0
#: the full 1-2 motor search (exp8, 2026-10-03) on the PC: 15 workers,
#: i5-12600K.  Stage wall times, s, and the tasks behind them.
EXP8 = {"workers": 15, "screen_s": 2314.0, "coarse_s": 17319.0,
        "refine_s": 10258.0, "refine_tasks": 112, "swap_s": 5621.0,
        "swap_refine_tasks": 18, "final_weather_s": 73.0, "total_s": 35586.0}


def _designs(winners):
    """The winners and fixed perturbations of them (lighter and heavier
    sensors, shorter spans, other taper) - some fly, some do not."""
    out = []
    for d in winners:
        out.append(d)
        for k_span, k_sens, taper in ((0.80, 0.85, 0.6), (0.95, 1.10, 0.9),
                                      (0.70, 0.60, 0.35), (1.00, 1.25, None),
                                      (0.90, 0.95, 0.75)):
            out.append(replace(d, span_m=d.span_m * k_span,
                               chord_m=d.chord_m * k_span,
                               sensor_mass_kg=d.sensor_mass_kg * k_sens,
                               sensor_m3_kg=d.sensor_m3_kg * k_sens,
                               taper=d.taper if taper is None else taper))
    return out


def _row(r):
    return {k: (asdict(v) if k == "design" else v) for k, v in r.items()
            if isinstance(v, (int, float, str, bool, type(None))) or k == "design"}


def _eval(task):
    d, refs = task
    return _row(evaluate_design(d, Config(), DX, *refs))


def _refine(task):
    d, refs = task
    cfg = Config()
    row = evaluate_design(d, cfg, DX, *refs)
    sp = DesignSpace()
    sp.motor_counts = [1, 2]
    st = SearchSettings.thorough()
    st.nm_iters = 4
    return _row(_refine_one((row, cfg, sp, st, refs[0], refs[1], refs[2])))


def _load_test(designs, refs, workers):
    """Seconds of one worker per design with every worker busy."""
    reps = max(1, int(math.ceil(3.0 * workers / len(designs))))
    tasks = [(d, refs) for d in designs] * reps
    t0 = time.time()
    with Pool(workers) as pool:
        pool.map(_eval, tasks, chunksize=1)
    return (time.time() - t0) * workers / len(tasks)


def make_reference(winners, refs, workers):
    designs = _designs(winners)
    t0 = time.time()
    rows = [_eval((d, refs)) for d in designs]
    serial = (time.time() - t0) / len(designs)
    with Pool(min(2, workers)) as pool:
        refined = pool.map(_refine, [(d, refs) for d in winners])
    loaded = _load_test(designs, refs, workers)
    os.makedirs(os.path.dirname(REF_FILE), exist_ok=True)
    json.dump({"made": time.strftime("%Y-%m-%d %H:%M"),
               "platform": platform.platform(),
               "python": "%s %s" % (platform.python_implementation(),
                                    platform.python_version()),
               "refs": refs, "dx": DX, "workers": workers,
               "serial_s_per_design": serial,
               "loaded_s_per_design": loaded,
               "designs": [asdict(d) for d in designs],
               "rows": rows, "refined": refined},
              open(REF_FILE, "w"), indent=1)
    print("reference written: %d designs, %d refinements, %.2f s/design "
          "serial, %.2f s/design with %d workers busy"
          % (len(designs), len(refined), serial, loaded, workers))


def _compare(a, b, path=""):
    """(number of values compared, number identical, max relative
    difference, first mismatch)."""
    n = same = 0
    worst, first = 0.0, None
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                n += 1
                first = first or "%s.%s missing" % (path, k)
                worst = max(worst, float("inf"))
                continue
            n2, s2, w2, f2 = _compare(a[k], b[k], "%s.%s" % (path, k))
            n, same, worst = n + n2, same + s2, max(worst, w2)
            first = first or f2
        return n, same, worst, first
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        if len(a) != len(b):
            return 1, 0, float("inf"), "%s: %d vs %d items" % (path, len(a), len(b))
        for i, (x, y) in enumerate(zip(a, b)):
            n2, s2, w2, f2 = _compare(x, y, "%s[%d]" % (path, i))
            n, same, worst = n + n2, same + s2, max(worst, w2)
            first = first or f2
        return n, same, worst, first
    if isinstance(a, bool) or isinstance(b, bool) or isinstance(a, str) \
            or a is None or b is None:
        ok = a == b
        return 1, int(ok), (0.0 if ok else float("inf")), \
            (None if ok else "%s: %r vs %r" % (path, a, b))
    a, b = float(a), float(b)
    if a == b:
        return 1, 1, 0.0, None
    rel = abs(a - b) / max(abs(a), abs(b), 1e-300)
    return 1, 0, rel, "%s: %.17g vs %.17g" % (path, a, b)


def estimate(k, workers):
    """Wall time of an exp8-sized search on `workers` workers, each `k`
    times the PC's per-worker time."""
    e = EXP8
    w = e["workers"]
    t_refine = e["refine_s"] / math.ceil(e["refine_tasks"] / float(w))
    swap_refine = math.ceil(e["swap_refine_tasks"] / float(w)) * t_refine
    spread = (e["screen_s"] + e["coarse_s"]) * w \
        + max(e["swap_s"] - swap_refine, 0.0) * w
    rounds = (math.ceil(e["refine_tasks"] / float(workers))
              + math.ceil(e["swap_refine_tasks"] / float(workers)))
    return k * (spread / workers + rounds * t_refine + e["final_weather_s"])


def check(workers):
    ref = json.load(open(REF_FILE))
    refs = ref["refs"]
    designs = [Design(**d) for d in ref["designs"]]
    print("reference: %s, %s, made %s" % (ref["python"], ref["platform"], ref["made"]))
    t0 = time.time()
    with Pool(workers) as pool:
        rows = pool.map(_eval, [(d, refs) for d in designs], chunksize=1)
        refined = pool.map(_refine, [(designs[i], refs) for i in
                                     range(0, len(designs), len(designs) // 2)])
    print("scored %d designs and ran %d short refinements in %.0f s"
          % (len(rows), len(refined), time.time() - t0))
    ok = True
    for label, mine, theirs in (("design scores", rows, ref["rows"]),
                                ("refinements", refined, ref["refined"])):
        n, same, worst, first = _compare(mine, theirs)
        if same == n:
            verdict = "identical (%d values, bit for bit)" % n
        elif worst < 1e-9:
            verdict = ("same to %.0e (%d of %d values bit for bit) - last-digit "
                       "differences from this platform's maths library" % (worst, same, n))
        else:
            ok = False
            verdict = "DIFFERENT: %d of %d values match; largest relative " \
                      "difference %.3g; first: %s" % (same, n, worst, first)
        print("  %-14s %s" % (label, verdict))
    loaded = _load_test(designs, refs, workers)
    k = loaded / ref["loaded_s_per_design"]
    print("speed: %.2f s per design per worker with all %d busy (PC: %.2f s "
          "with %d) - each worker runs at %.2fx the PC's"
          % (loaded, workers, ref["loaded_s_per_design"], ref["workers"],
             1.0 / k))
    print("estimated wall time of a full exp8-size search (PC: %.1f h):"
          % (EXP8["total_s"] / 3600.0))
    for n in sorted({workers, 32, 64, 128, 256, 512}):
        print("  %4d workers  %5.1f h%s" % (n, estimate(k, n) / 3600.0,
                                            "   <- this job" if n == workers else ""))
    print("  (refinement runs one Nelder-Mead per configuration, so past about "
          "%d workers only the screen, coarse sweep and swap scoring speed up)"
          % EXP8["refine_tasks"])
    return 0 if ok else 1


def main():
    hpc.use_fork()
    workers = hpc.allotted_cpus()
    if "--workers" in sys.argv:
        workers = int(sys.argv[sys.argv.index("--workers") + 1])
    if "--make-reference" in sys.argv:
        src = sys.argv[sys.argv.index("--make-reference") + 1]
        rows = [r for r in json.load(open(src)) if r.get("feasible")]
        winners = [Design(**{k[2:]: v for k, v in
                             next(x for x in rows if x["d_n_motors"] == n).items()
                             if k.startswith("d_")}) for n in (1, 2)]
        refs = [rows[0]["m2_raw_lb_per_s"] / rows[0]["m2_ratio"],
                rows[0]["m3_raw_lap_lb"] / rows[0]["m3_ratio"],
                rows[0]["gm_raw_lb_in"] / rows[0]["gm_ratio"]]
        make_reference(winners, refs, workers)
        return 0
    hpc.banner(workers, "check")
    return check(workers)


if __name__ == "__main__":
    sys.exit(main())
