#!/usr/bin/env python3
"""Run the optimiser on a cluster node (or several).

Takes every run_optimizer.py argument (python run_cluster.py --help lists
them), plus:

  --mpi    spread the workers over every MPI rank, across nodes.  Launch as
           srun python -m mpi4py.futures run_cluster.py --mpi ...

Without --mpi it uses every core Slurm gave the job on this node; pass
--processes N to choose yourself.
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import hpc  # noqa: E402


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    mpi = "--mpi" in argv
    argv = [a for a in argv if a != "--mpi"]
    if "-h" in argv or "--help" in argv:
        print(__doc__)
    hpc.use_fork()
    if mpi:
        executor, workers = hpc.start_mpi()
        hpc.install_executor(executor, workers)
        mode = "MPI ranks"
    else:
        if "--processes" in argv:
            workers = int(argv[argv.index("--processes") + 1])
        else:
            workers = hpc.allotted_cpus()
            argv += ["--processes", str(workers)]
        mode = "process pool"
    hpc.banner(workers, mode)
    import run_optimizer
    timer = hpc.Timer()
    rc = run_optimizer.main(argv)
    print("cluster run finished in %s" % timer)
    if mpi:
        executor.shutdown()
    return rc or 0


if __name__ == "__main__":
    sys.exit(main())
