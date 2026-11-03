# DBF optimiser: cluster bundle (UGA Sapelo2)

This bundle is the main project's search code, copied file for file, with a thin layer for running it on a Slurm cluster. A design scores the same here as on the PC. What changes is how many cores do the work:

| Layout | How | Best for |
|---|---|---|
| One node, all its cores | `slurm/search.sbatch` | one full search |
| Many searches at once | `slurm/array.sbatch` + a jobs file | seeds, motor counts, battery classes |
| One search over several nodes | `slurm/mpi_search.sbatch` (MPI) | the biggest single searches |

**What was left out:** the GUI, web app, rulebook reader, test suites, docs and old results. None of them runs during a search.

**Same answer at any worker count:** since m99 (2026-10-07), every design is scored from empty shared tables. Before that, a design scored right after a near neighbour could pick up the neighbour's cached lifting-line and propulsion numbers. So a search on 64 workers could rank designs differently, in the last digits, from the same search on the PC's 15.

## Time estimates

A full search at the exp8 settings took 9.9 h on the PC (15 workers):
- 49% coarse sweep;
- 29% refinement;
- 16% swap test;
- 6% propulsion screen.

The screen, coarse sweep and swap scoring split perfectly across cores. Refinement runs one Nelder–Mead per configuration (about 100 of them). Each is a chain of evaluations, one after another, so past about 100 workers refinement stops getting faster.

`check_cluster.py` measures this node's speed and prints the estimate for 32 to 512 workers.

## 1. Build and upload (on the PC)

```bash
python cluster/make_cluster_bundle.py
scp ../dbf-optimizer-cluster.tar.gz MyID@xfer.gacrc.uga.edu:/scratch/MyID/
```

Then on Sapelo2:

```bash
ssh MyID@sapelo2.gacrc.uga.edu
cd /scratch/MyID && tar xzf dbf-optimizer-cluster.tar.gz && cd dbf-optimizer-cluster
```

Rebuild and re-upload whenever the main code or catalogue changes. Nothing in the bundle is edited by hand. `BUNDLE_INFO.json` records the source, date and checksums.

## 2. Check it (about 5 minutes)

```bash
sbatch slurm/check.sbatch
cat logs/dbf-check_*.out
```

It checks two things:

- **Same answers as the PC?** It scores 12 fixed designs and runs two short refinements, then compares every number with the PC's own answers in `reference/`.
  - "identical" means bit for bit.
  - "same to 1e-12" means the last digits differ, from the platform's maths library. That's harmless for any one score, but a refinement can take a slightly different path, so a whole search lands close to the PC's rather than identical.
- **How fast?** It times the node with every allotted core busy, then estimates a full search's wall time.

## 3. Run

```bash
sbatch slurm/search.sbatch                 # one full search (exp8 settings)
sbatch slurm/search.sbatch --seed 7        # extra arguments pass through to run_optimizer.py
sbatch --array=1-4 slurm/array.sbatch slurm/jobs_example.txt    # four searches in parallel
sbatch --array=1-8%4 slurm/array.sbatch slurm/jobs_battery.txt  # battery classes, 4 at a time
sbatch --array=1-6 --mem=48G --time=24:00:00 slurm/array.sbatch slurm/jobs_chord24_v2.txt  # 24 in chord, 30-50 mph stall, new physics
sbatch --array=1-6 --mem=48G --time=24:00:00 slurm/array.sbatch slurm/jobs_chord24_v3.txt  # the same, filling all 64 cores
sbatch --array=1-5 --mem=48G --time=24:00:00 slurm/array.sbatch slurm/jobs_chord24_v4.txt  # stall limit in the hot day's air
sbatch --array=1-7 --mem=48G --time=24:00:00 slurm/array.sbatch slurm/jobs_chord24_v5.txt  # the limit on the true stall (lift factor 1.0), 20-50 mph
sbatch --array=1-8 --mem=48G --time=24:00:00 slurm/array.sbatch slurm/jobs_nochordcap_v6.txt  # v5 with no chord cap, and a fresh no-limit search
squeue --me                                # what's queued or running
seff <jobid>                               # CPU and memory efficiency after it ends
```

- **Results** land in `results/<name>_<jobid>/`, in the same files as on the PC: `ranked_designs.json`, `top25/`, `search_log.txt` and the others.
- **Copy them home before the scratch purge:**

  ```bash
  scp -r MyID@xfer.gacrc.uga.edu:/scratch/MyID/dbf-optimizer-cluster/results/<run> .
  ```

- **Jobs files** have one search per line, in the form `name | run_cluster.py arguments`. `jobs_battery.txt` is generated when the bundle is built: one full-scale search per cell count, over every pack of that count. `jobs_chord24_v2.txt` reruns the 24 in chord / stall-limit sweep on the 2026-10-14 physics (section lift from the polars, winglets, tails from their requirements), with a no-limit reference search as its sixth line. `jobs_chord24_v3.txt` is the same sweep with `--fill-cores`. `jobs_chord24_v4.txt` reruns its five stall limits with the limit checked in the hottest weather preset's air (m118); v3's no-limit search still stands, as the change only touches the stall check. `jobs_chord24_v5.txt` is v4 with `--stall-lift-fos 1.0`, so each limit is the true stall, power off at the Mission 2 weight on the wing's full lift (v4's limits, on lift / 1.5, were true stalls of 25-42 mph), from 20 to 50 mph. 20 mph is below what v4's lightest airframes reach with the smallest legal sensor, so that line may stop after its coarse sweep with "no feasible design". `jobs_nochordcap_v6.txt` is v5 with no chord cap (the chord bounded only by the 6 ft span and aspect ratio 1-14), plus a fresh no-limit search on the current code as its eighth line.
- **Faster under a stall limit (m120, m121).** The coarse sweep is 50-75% of a stall-limited run, and most of it went on fitting each sample's sensor to the limit. The fit is now a secant (2-4 aircraft builds, not 15), and the sweep's worker about 1.7× faster; expect a 50 mph search to drop from about 10 h to about 7 h. A search whose coarse sweep flies nothing now goes on from its seed designs instead of stopping - v6's 25 mph line stopped that way.
- **Fill the cores** with `--fill-cores` (m114). The v2 30 mph search used 68% of its 64 cores (seff, job 49175948_1):
  - its refinement had only 26 configurations with a flying sample, so 38 cores sat idle for 25 minutes;
  - the swap test refined 18 of the 1,746 swaps that beat their parent;
  - the weather check had 32 finalists.

  A refinement runs many minutes on one core, so a stage lasts as long as its last round. `--fill-cores` rounds the refinements, the swaps' refinements and the finalists up to whole rounds of the workers (100 → 128, 18 → 64, 32 → 64). With too few configurations, the best are refined from more than one start. The wall time stays about the same and more of the promising designs are refined. `--swap-refine N` sets the swaps' refinements directly.
- **Cores, memory and time** are set at the top of each `.sbatch` file: 64 cores, 32 GB and 12 h by default. A full search took 2.0 h on 64 cores (job 48862616, 2026-10-08). seff reported 1.7 GB, but that was one process: with 8 GB, two of five searches lost their workers (2026-10-14). Use `sinfo -p batch -o "%c %m %D"` to see node sizes.
  - A 128-core node roughly halves the coarse sweep.
  - Asking for fewer cores gets through the queue sooner.

## 4. Faster still: PyPy (optional)

The code is pure standard-library Python with no numpy, so PyPy can run it unchanged, often several times faster.

```bash
./install_pypy.sh pypy3.11-v7.3.20-linux64      # current name at https://pypy.org/download.html
DBF_PYPY=1 sbatch slurm/check.sbatch            # does PyPy match the PC, and how fast is it?
DBF_PYPY=1 sbatch slurm/search.sbatch
```

Use it only if the check says "identical" or "same to 1e-…". PyPy may add floats in a different order from CPython 3.12+, which would show up there.

## 5. Several nodes (MPI, optional)

Edit `--nodes` in `slurm/mpi_search.sbatch`, then:

```bash
sbatch slurm/mpi_search.sbatch
```

It loads the cluster's `mpi4py` module. `run_cluster.py --mpi` then puts an MPI pool in place of the process pool.

**Catalogue CSVs:** with MPI, the workers import the package fresh. So files you add with `--motors-csv` or `--packs-csv` at run time don't reach them. Put those files in `catalogues/` instead, which loads at import, as on the PC.

## Python

- The default is the newest CPython 3.12+ module on the cluster (`sapelo2_env.sh`).
- **Why 3.12 or newer:** 3.12 changed how `sum()` adds floats, and the PC runs 3.13. On 3.11 and older, results differ in the last digits.
- To use another interpreter, set `DBF_PYTHON=/path/to/python`.
