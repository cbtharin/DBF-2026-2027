# Interpreter for the optimiser on Sapelo2.  Source it from a job script:
#     source ./sapelo2_env.sh
# It sets $PY.  The optimiser is pure standard-library Python: nothing to
# pip install.
#
#   default        the newest CPython >= 3.12 among the cluster's modules
#                  (3.12 changed how sum() adds floats; the PC runs 3.13, so
#                  an older Python gives results that differ in the last
#                  digits and refinements that wander slightly)
#   DBF_PYPY=1     PyPy from ~/pypy-dbf (install_pypy.sh) - usually several
#                  times faster on this code.  Run check_cluster.py with it
#                  first: it says whether PyPy reproduces the PC's numbers.
#   DBF_PYTHON=... any interpreter you name
module purge >/dev/null 2>&1
if [ -n "$DBF_PYTHON" ]; then
    PY="$DBF_PYTHON"
elif [ -n "$DBF_PYPY" ]; then
    PY="$HOME/pypy-dbf/bin/pypy3"
    [ -x "$PY" ] || { echo "no PyPy at $PY - run ./install_pypy.sh first" >&2; return 1 2>/dev/null || exit 1; }
else
    mod=$(module -t avail Python 2>&1 | grep -E '^Python/3\.(1[2-9]|[2-9][0-9])\.' | sort -V | tail -1)
    if [ -n "$mod" ]; then
        module load "$mod"
    else
        echo "WARNING: no Python >= 3.12 module found (module spider Python); using $(command -v python3)" >&2
    fi
    PY=python3
fi
export PY
"$PY" -c 'import platform, sys; print("python:", platform.python_implementation(), platform.python_version(), sys.executable)'
