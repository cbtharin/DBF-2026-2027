#!/usr/bin/env python3
"""Launcher for the DBF 2026-27 optimiser's desktop front end.

    python dbf_gui.py

Equivalent to `python -m dbf2027.gui`; this file exists so the GUI can be
started by double-clicking it, and so the `__main__` guard below is in a
module that is safe for Windows to re-import.  It has to be: the search
uses `multiprocessing`, Windows starts child processes by re-importing
`__main__`, and without the guard each worker would open its own window.
"""
from multiprocessing import freeze_support

from dbf2027.gui import main

if __name__ == "__main__":
    freeze_support()
    raise SystemExit(main())
