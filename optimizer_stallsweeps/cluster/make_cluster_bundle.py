#!/usr/bin/env python3
"""Build the cluster bundle from this project.

    python cluster/make_cluster_bundle.py [--out DIR] [--no-reference]

Writes DIR (default: dbf-optimizer-cluster next to this project) and
DIR.tar.gz to upload.  The bundle is:

  dbf2027/        the modules a search imports, copied byte for byte (the
                  GUI, web app and rulebook reader are left out) - so the
                  physics is the PC's
  run_optimizer.py, catalogues/, examples/   as in the project
  seeds/          the saved designs the PC's full searches were seeded with
  hpc.py, run_cluster.py, check_cluster.py, sapelo2_env.sh, install_pypy.sh,
  slurm/          the cluster layer (from cluster/ in this project)
  reference/      the PC's own answers for check_cluster.py
  BUNDLE_INFO.json  what was copied, from where, with checksums

Rebuild it whenever the main code changes: nothing in the bundle is edited
by hand, so it never drifts from the project.
"""
from __future__ import annotations

import ast
import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
RESULTS = os.path.join(SRC, "results")
#: seed designs: every search's best design, plus the saved top designs of
#: the full searches (the seeds of exp8 and the battery study)
SEED_SOURCES = ["*/best_design.json", "thorough_2026-09-30/retuned",
                "thorough_2026-10-01_seeded/top25",
                "exp3_A_bigprops_1to4_2026-10-02/top25",
                "exp3_C_cap2_2026-10-02/top25",
                "exp3_C_cap2_nacCd010_2026-10-03/top25",
                "exp3_C_cap2_nacCd030_2026-10-03/top25",
                "exp6_C_cap2_2026-10-03/top25",
                "exp8_C_cap2_2026-10-03/top25",
                "exp8_C_cap2_seed7_2026-10-04/top25"]
REFERENCE_RUN = os.path.join(RESULTS, "exp8_C_cap2_2026-10-03", "ranked_designs.json")
ADDONS = ["hpc.py", "run_cluster.py", "check_cluster.py", "sapelo2_env.sh",
          "install_pypy.sh"]
EXECUTABLE = (".sh", "run_cluster.py", "check_cluster.py")
#: the full-search settings (exp8) the generated jobs files use
FULL = ("--thorough --top 25 --finalists 32 --motor-counts 1 2 "
        "--keep-propulsion 600 --refine 100 --swap-parents 12 --seeds seeds")


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


def search_modules():
    """The dbf2027 modules a search imports (run_optimizer and the report
    writers it calls), found by importing them in a fresh interpreter."""
    code = ("import sys; sys.argv = ['x']; import run_optimizer, dbf2027.report; "
            "print(' '.join(sorted(m.split('.')[-1] for m in sys.modules "
            "if m.startswith('dbf2027.'))))")
    out = subprocess.run([sys.executable, "-B", "-c", code], cwd=SRC,
                         capture_output=True, text=True,
                         env=dict(os.environ, PYTHONPATH=SRC))
    if out.returncode != 0:
        raise SystemExit(out.stderr)
    return out.stdout.split()


def copy(src, dst, files):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)
    files[os.path.relpath(dst, OUT).replace(os.sep, "/")] = sha(dst)


def collect_seeds(files):
    sys.path.insert(0, SRC)
    from dbf2027.aircraft import Design
    seen, n = set(), 0
    for pat in SEED_SOURCES:
        for p in sorted(glob.glob(os.path.join(RESULTS, pat))):
            for f in ([p] if os.path.isfile(p) else sorted(glob.glob(os.path.join(p, "*.json")))):
                try:
                    data = json.load(open(f, encoding="utf-8"))
                    Design(**data)
                except (OSError, ValueError, TypeError):
                    continue
                key = json.dumps(data, sort_keys=True)
                if key in seen:
                    continue
                seen.add(key)
                rel = os.path.relpath(f, RESULTS).replace(os.sep, "__")
                copy(f, os.path.join(OUT, "seeds", rel), files)
                n += 1
    return n


def battery_jobs():
    """One full-scale search per cell count, every pack of that count."""
    sys.path.insert(0, SRC)
    from dbf2027.components import PACKS, RETIRED
    by_s = {}
    for k, p in PACKS.items():
        if k not in RETIRED:
            by_s.setdefault(p.series, []).append(k)
    lines = ["# name | run_cluster.py arguments: one full-scale search per cell count, every",
             "# searchable pack of that count (LiPo, LiHV and Li-ion).  Generated %s."
             % time.strftime("%Y-%m-%d"),
             "#     sbatch --array=1-%d slurm/array.sbatch slurm/jobs_battery.txt" % len(by_s)]
    for s in sorted(by_s):
        lines.append("battery_%dS | %s --no-seed-refine --packs %s"
                     % (s, FULL, ",".join(sorted(by_s[s]))))
    return "\n".join(lines) + "\n"


def oldest_python(paths):
    """The oldest Python 3 grammar every bundled module parses under."""
    for minor in range(8, 14):
        try:
            for p in paths:
                ast.parse(open(p, encoding="utf-8").read(), p,
                          feature_version=(3, minor))
            return "3.%d" % minor
        except SyntaxError:
            continue
    return "3.13"


def main():
    global OUT
    args = sys.argv[1:]
    OUT = os.path.abspath(args[args.index("--out") + 1] if "--out" in args
                          else os.path.join(os.path.dirname(SRC), "dbf-optimizer-cluster"))
    # generated parts are replaced; results/ and logs/ of earlier runs stay
    for sub in ("dbf2027", "seeds", "slurm", "catalogues", "examples", "reference"):
        shutil.rmtree(os.path.join(OUT, sub), ignore_errors=True)
    files = {}
    mods = search_modules()
    for m in ["__init__"] + mods:
        copy(os.path.join(SRC, "dbf2027", m + ".py"), os.path.join(OUT, "dbf2027", m + ".py"), files)
    for f in sorted(os.listdir(os.path.join(SRC, "dbf2027"))):
        if not f.endswith((".py", ".pyc")) and os.path.isfile(os.path.join(SRC, "dbf2027", f)):
            copy(os.path.join(SRC, "dbf2027", f), os.path.join(OUT, "dbf2027", f), files)
    copy(os.path.join(SRC, "run_optimizer.py"), os.path.join(OUT, "run_optimizer.py"), files)
    for f in glob.glob(os.path.join(SRC, "catalogues", "*")) + glob.glob(os.path.join(SRC, "examples", "*.json")):
        copy(f, os.path.join(OUT, os.path.relpath(f, SRC)), files)
    for f in ADDONS:
        copy(os.path.join(HERE, f), os.path.join(OUT, f), files)
    copy(os.path.join(HERE, "README-cluster.md"), os.path.join(OUT, "README.md"), files)
    for f in glob.glob(os.path.join(HERE, "slurm", "*")):
        copy(f, os.path.join(OUT, "slurm", os.path.basename(f)), files)
    with open(os.path.join(OUT, "slurm", "jobs_battery.txt"), "w", newline="\n") as fh:
        fh.write(battery_jobs())
    files["slurm/jobs_battery.txt"] = sha(os.path.join(OUT, "slurm", "jobs_battery.txt"))
    n_seeds = collect_seeds(files)
    for d in ("logs", "results"):
        os.makedirs(os.path.join(OUT, d), exist_ok=True)
        open(os.path.join(OUT, d, ".keep"), "a").close()
    py = [os.path.join(OUT, f) for f in files if f.endswith(".py")]
    min_py = oldest_python(py)

    if "--no-reference" not in args:
        print("writing the PC's reference answers (about 3 minutes)...")
        rc = subprocess.run([sys.executable, "-B", "check_cluster.py",
                             "--make-reference", REFERENCE_RUN], cwd=OUT).returncode
        if rc != 0:
            raise SystemExit("reference run failed")
        files["reference/check_reference.json"] = sha(os.path.join(OUT, "reference", "check_reference.json"))

    info = {"built": time.strftime("%Y-%m-%d %H:%M"), "source": SRC,
            "python_used": sys.version.split()[0],
            "grammar_parses_from": min_py,
            "modules": mods, "seed_designs": n_seeds, "files": files}
    json.dump(info, open(os.path.join(OUT, "BUNDLE_INFO.json"), "w"), indent=1)

    tgz = OUT + ".tar.gz"
    with tarfile.open(tgz, "w:gz") as tar:
        def keep(ti):
            if "__pycache__" in ti.name or ti.name.endswith(".pyc"):
                return None
            base = os.path.basename(ti.name)
            if ti.isfile() and (base.endswith(EXECUTABLE)):
                ti.mode = 0o755
            return ti
        # results/ and logs/ go in empty: earlier runs are not re-uploaded
        for item in sorted(set(f.split("/")[0] for f in files) | {"BUNDLE_INFO.json"}) \
                + ["logs/.keep", "results/.keep"]:
            tar.add(os.path.join(OUT, item), "dbf-optimizer-cluster/" + item, filter=keep)
    print("bundle: %s (%d files, %d seed designs, modules: %s)"
          % (OUT, len(files), n_seeds, " ".join(mods)))
    print("upload: %s (%.1f MB); grammar parses on Python %s and later"
          % (tgz, os.path.getsize(tgz) / 1e6, min_py))


if __name__ == "__main__":
    main()
