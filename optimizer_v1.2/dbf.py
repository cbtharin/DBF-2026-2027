#!/usr/bin/env python3
"""DBF 2026-27 mission simulator, design optimiser and validation suite.

One entry point for everything.  The code itself lives in the `dbf2027/`
package next to this file; this is only a dispatcher, so that every
command is one obvious line no matter which module implements it.

    python dbf.py --gui                   the desktop front end
    python dbf.py --web                   the same thing in a browser
    python dbf.py --quick                 fast search
    python dbf.py                         default search
    python dbf.py --thorough              wide search
    python dbf.py --evaluate design.json  score one aeroplane, print its parts
    python dbf.py --catalogue             list every function in the package
    python dbf.py --rules                 every rule, and where it is enforced
    python dbf.py --selftest              physics and scoring invariants
    python dbf.py --validate              catalogue, bounds, determinism
    python dbf.py --guitest               the GUI, driven headlessly
    python dbf.py --webtest               the web front end, over real HTTP
    python dbf.py --convergence           numerical error study
    python dbf.py --audit ranked.csv      which options reached the finals
    python dbf.py --version

Pure Python standard library: no pip install, no numpy, no matplotlib.
Python 3.8+.  The physics, every derivation and every assumption are
documented in README.md; FUNCTIONS.md indexes the code.

This file used to be a generated single-file build of the whole package.
It is not any more - the package ships as modules, which is easier to
read, to test and to change.
"""
from __future__ import annotations

import os
import sys

__version__ = "1.11.0"

# Run correctly when invoked from another directory.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    if "--version" in argv:
        print("DBF 2026-27 optimiser v%s" % __version__)
        return 0

    if "--gui" in argv:
        from multiprocessing import freeze_support
        from dbf2027.gui import main as gui_main
        freeze_support()
        return gui_main()

    if "--web" in argv:
        from multiprocessing import freeze_support
        from dbf2027.webapp import main as web_main
        freeze_support()
        return web_main(argv)

    if "--rules" in argv:
        from dbf2027.rulebook import main as rules_main
        i = argv.index("--rules")
        return rules_main(argv[i + 1:])

    if "--catalogue" in argv:
        from dbf2027.catalogue import main as catalogue_main
        i = argv.index("--catalogue")
        return catalogue_main(argv[i + 1:])

    if "--selftest" in argv:
        from selftest import main as selftest_main
        return selftest_main()

    if "--validate" in argv:
        from validate import main as validate_main
        return validate_main()

    if "--guitest" in argv:
        from guitest import main as guitest_main
        return guitest_main()

    if "--webtest" in argv:
        from webtest import main as webtest_main
        return webtest_main()

    if "--convergence" in argv:
        from grid_convergence import main as convergence_main
        return convergence_main()

    if "--audit" in argv:
        from analyse_results import main as audit_main
        i = argv.index("--audit")
        path = (argv[i + 1] if i + 1 < len(argv)
                else "results/ranked_designs.csv")
        return audit_main(["audit", path])

    from run_optimizer import main as cli_main
    return cli_main(argv)


if __name__ == "__main__":
    from multiprocessing import freeze_support
    freeze_support()
    sys.exit(main())
