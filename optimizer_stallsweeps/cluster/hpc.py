"""Cluster execution layer for the optimiser (UGA Sapelo2 / any Slurm cluster).

The physics is not touched: the bundle's dbf2027 package is a verbatim copy
of the main project's, so a design scores the same here as on the PC.  This
module only decides WHERE the evaluations run:

  * one node (default): a fork-started process pool sized to the cores
    Slurm gave the job - not the node's whole core count, which on a shared
    node would oversubscribe the cores other jobs own;
  * several nodes (--mpi): an mpi4py.futures pool over every MPI rank, put
    in place of optimizer._pool_map.  The search code looks the function up
    by name each time it maps a stage, so nothing else changes.

Every stage after the propulsion screen is a map over independent designs,
so both layouts return exactly what the PC's pool returns, in the same
order.
"""
from __future__ import annotations

import multiprocessing as mp
import os
import platform
import socket
import sys
import time
from typing import Any, Callable, List, Optional


def allotted_cpus() -> int:
    """Cores this job may use: the CPU affinity Slurm's cgroup sets,
    capped by SLURM_CPUS_PER_TASK when that is set."""
    n = 0
    if hasattr(os, "sched_getaffinity"):
        try:
            n = len(os.sched_getaffinity(0))
        except OSError:
            n = 0
    if n <= 0:
        n = os.cpu_count() or 1
    s = os.environ.get("SLURM_CPUS_PER_TASK", "")
    if s.isdigit() and int(s) > 0:
        n = min(n, int(s))
    return max(1, n)


def use_fork() -> None:
    """Fork-start the workers on Linux: they inherit the loaded catalogues
    and start in milliseconds (Python 3.14 makes forkserver the default)."""
    if sys.platform.startswith("linux"):
        try:
            mp.set_start_method("fork", force=True)
        except (RuntimeError, ValueError):
            pass


def start_mpi():
    """The MPI worker pool and its size (every rank but this one).  Launch
    with:  srun python -m mpi4py.futures run_cluster.py --mpi ..."""
    try:
        from mpi4py import MPI
        from mpi4py.futures import MPIPoolExecutor
    except ImportError as exc:
        raise SystemExit("--mpi needs mpi4py (module spider mpi4py): %s" % exc)
    size = MPI.COMM_WORLD.Get_size()
    if size < 2:
        raise SystemExit("--mpi needs at least 2 MPI ranks; launch with "
                         "srun python -m mpi4py.futures run_cluster.py --mpi")
    return MPIPoolExecutor(), size - 1


def install_executor(executor, workers: int) -> None:
    """Put an executor (concurrent.futures interface) in place of the
    search's process pool.  Same contract as optimizer._pool_map: results
    in task order; a handful of tasks run here, serially."""
    from dbf2027 import optimizer as opt

    def _pool_map(func: Callable, tasks: List[Any], processes: int,
                  on_done: Optional[Callable[[int, int], None]] = None,
                  min_parallel: int = 8) -> List[Any]:
        n = len(tasks)
        if processes == 1 or n < max(min_parallel, 2):
            out = []
            for i, t in enumerate(tasks):
                out.append(func(t))
                if on_done is not None:
                    on_done(i + 1, n)
            return out
        chunk = max(1, n // (workers * 8))
        out = list(executor.map(func, tasks, chunksize=chunk))
        if on_done is not None:
            on_done(n, n)
        return out

    opt._pool_map = _pool_map
    opt._n_processes = lambda settings: workers


def banner(workers: int, mode: str) -> None:
    impl = platform.python_implementation()
    print("cluster run  host %s, %s %s, %d workers (%s), Slurm job %s"
          % (socket.gethostname(), impl, platform.python_version(), workers,
             mode, os.environ.get("SLURM_JOB_ID", "-")))
    if impl == "CPython" and sys.version_info < (3, 12):
        print("  NOTE: CPython %s sums floats differently from 3.12+ (the PC "
              "runs 3.13), so results can differ in the last digits and a "
              "refinement can take a different path.  Load Python >= 3.12."
              % platform.python_version())
    sys.stdout.flush()


class Timer:
    def __init__(self):
        self.t0 = time.time()

    def __str__(self):
        s = time.time() - self.t0
        return "%.0f s (%.2f h)" % (s, s / 3600.0)
