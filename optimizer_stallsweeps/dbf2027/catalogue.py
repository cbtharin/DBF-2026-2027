"""Function catalogue: what every callable in the package is, and where.

Built by introspection rather than maintained by hand, so it cannot drift
out of date the way a hand-written index does.  It walks the package's
modules, reads signatures with `inspect`, takes the first line of each
docstring as the summary, and emits either plain text, Markdown, or a
list of records for the GUI to display and filter.

    python -m dbf2027.catalogue                 print it
    python -m dbf2027.catalogue --markdown      write FUNCTIONS.md
    python -m dbf2027.catalogue --search trim   only entries matching "trim"

Anything whose name starts with an underscore is internal and left out
unless you ask for it with `include_private`.
"""
from __future__ import annotations

import importlib
import inspect
import sys
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

# Order is the dependency order of the package: units first, report last.
# It is also a reasonable reading order for someone new to the code.
MODULE_ORDER = [
    "units", "numerics", "config", "components", "payload", "liftingline",
    "loads", "aerodynamics", "weights", "layout", "propulsion", "aircraft",
    "mission",
    "scoring", "optimizer", "report", "plots", "gui", "webapp",
]

MODULE_BLURB = {
    "units": "Unit conversions. Everything internal to the model is SI.",
    "numerics": "Root finding, quadrature, integration, interpolation and "
                "optimisation, all written out rather than imported.",
    "config": "Environment, course, rules and build standard - every number "
              "you might want to change lives in one of these four.",
    "components": "Catalogues: airfoils, motors, propellers, COTS battery "
                  "packs and ESCs.",
    "payload": "The fixed shipping container, the sensor envelope it has to "
               "hold, and how boxes pack into a fuselage bay.",
    "liftingline": "Prandtl lifting line: span efficiency and lift-curve "
                   "slope from the planform.",
    "loads": "The V-n diagram - manoeuvre, gust and CLmax cases - which "
             "sets the load factor the spar is sized for.",
    "aerodynamics": "Lift, drag, Reynolds-dependent CLmax, longitudinal "
                    "trim and the neutral point.",
    "weights": "Itemised mass buildup and spar sizing.",
    "layout": "Where every part sits: nose, bay, wing and tail stations, battery trim, CG.",
    "propulsion": "Motor/propeller torque balance and the propulsion-system "
                  "architecture the rules impose.",
    "aircraft": "The design vector, and the sized aeroplane it produces.",
    "mission": "Takeoff, climb, turns, laps and landing, flown as an ODE.",
    "scoring": "The three mission scoring formulas.",
    "optimizer": "The four-stage search over the design space.",
    "report": "Dossier, ranking table, Pareto front and bill of materials.",
    "plots": "Chart geometry: series in data units, mapped to pixels.",
    "gui": "The desktop front end.",
    "webapp": "The web front end: a standard-library HTTP server and a "
              "JSON API over the same model.",
}


@dataclass
class Entry:
    """One catalogue record: a function, class or method."""
    module: str
    qualname: str
    kind: str                 # "function", "class" or "method"
    signature: str
    summary: str
    doc: str
    lineno: int

    @property
    def label(self) -> str:
        return "%s.%s" % (self.module, self.qualname)

    def matches(self, query: str) -> bool:
        q = query.lower().strip()
        if not q:
            return True
        return (q in self.qualname.lower() or q in self.module.lower()
                or q in self.summary.lower() or q in self.doc.lower())


def _first_line(obj: Any) -> str:
    doc = inspect.getdoc(obj) or ""
    for line in doc.splitlines():
        if line.strip():
            return line.strip()
    return ""


def _signature(obj: Any) -> str:
    try:
        sig = inspect.signature(obj)
    except (TypeError, ValueError):
        return "(...)"
    # Annotations are long and mostly noise in an index; drop them but keep
    # the defaults, which are the part you actually want to see.
    parts = []
    for name, par in sig.parameters.items():
        if name in ("self", "cls"):
            continue
        text = name
        if par.kind is par.VAR_POSITIONAL:
            text = "*" + name
        elif par.kind is par.VAR_KEYWORD:
            text = "**" + name
        if par.default is not par.empty:
            text += "=%r" % (par.default,)
        parts.append(text)
    return "(%s)" % ", ".join(parts)


def _lineno(obj: Any) -> int:
    try:
        return inspect.getsourcelines(obj)[1]
    except (OSError, TypeError):
        return 0


def collect(modules: Optional[List[str]] = None,
            include_private: bool = False,
            include_methods: bool = True) -> List[Entry]:
    """Walk the package and return one Entry per public callable.

    Only objects actually *defined* in the module are listed, so a name
    imported from elsewhere does not appear twice.
    """
    names = modules if modules is not None else MODULE_ORDER
    out: List[Entry] = []
    for short in names:
        try:
            mod = importlib.import_module("." + short, package=__package__)
        except Exception:                     # gui needs tkinter; skip if absent
            continue
        mod_name = mod.__name__.rsplit(".", 1)[-1]
        for name, obj in sorted(vars(mod).items()):
            if name.startswith("_") and not include_private:
                continue
            if getattr(obj, "__module__", None) != mod.__name__:
                continue
            if inspect.isfunction(obj):
                out.append(Entry(mod_name, name, "function", _signature(obj),
                                 _first_line(obj), inspect.getdoc(obj) or "",
                                 _lineno(obj)))
            elif inspect.isclass(obj):
                out.append(Entry(mod_name, name, "class", _signature(obj),
                                 _first_line(obj), inspect.getdoc(obj) or "",
                                 _lineno(obj)))
                if not include_methods:
                    continue
                for mname, meth in sorted(vars(obj).items()):
                    if mname.startswith("_") and not include_private:
                        continue
                    if isinstance(meth, property):
                        fn = meth.fget
                        if fn is None:
                            continue
                        out.append(Entry(
                            mod_name, "%s.%s" % (name, mname), "property",
                            "", _first_line(fn), inspect.getdoc(fn) or "",
                            _lineno(fn)))
                    elif inspect.isfunction(meth):
                        out.append(Entry(
                            mod_name, "%s.%s" % (name, mname), "method",
                            _signature(meth), _first_line(meth),
                            inspect.getdoc(meth) or "", _lineno(meth)))
                    elif isinstance(meth, staticmethod):
                        fn = meth.__func__
                        out.append(Entry(
                            mod_name, "%s.%s" % (name, mname), "method",
                            _signature(fn), _first_line(fn),
                            inspect.getdoc(fn) or "", _lineno(fn)))
    return out


def counts(entries: List[Entry]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for e in entries:
        out[e.kind] = out.get(e.kind, 0) + 1
    return out


def by_module(entries: List[Entry]) -> Dict[str, List[Entry]]:
    out: Dict[str, List[Entry]] = {}
    for e in entries:
        out.setdefault(e.module, []).append(e)
    return out


def as_text(entries: List[Entry], width: int = 78) -> str:
    """Plain text, grouped by module, for a terminal."""
    lines = ["=" * width, "DBF 2026-27 FUNCTION CATALOGUE", "=" * width]
    grouped = by_module(entries)
    for short in MODULE_ORDER:
        group = grouped.get(short)
        if not group:
            continue
        lines.append("")
        lines.append("%s.py  -  %s" % (short, MODULE_BLURB.get(short, "")))
        lines.append("-" * width)
        for e in group:
            head = "  %s%s" % (e.qualname, e.signature)
            if e.kind == "class":
                head = "  class %s" % e.qualname
            elif e.kind == "property":
                head = "  %s  [property]" % e.qualname
            lines.append(head)
            if e.summary:
                lines.append("      %s" % e.summary)
    c = counts(entries)
    lines.append("")
    lines.append("-" * width)
    lines.append("  %d entries: %s"
                 % (len(entries),
                    ", ".join("%d %s" % (v, k) for k, v in sorted(c.items()))))
    return "\n".join(lines)


def as_markdown(entries: List[Entry]) -> str:
    """Markdown, for committing next to the README."""
    c = counts(entries)
    out = ["# Function catalogue", "",
           "Every public callable in `dbf2027/`, generated by "
           "`python -m dbf2027.catalogue --markdown`. "
           "%d entries: %s."
           % (len(entries),
              ", ".join("%d %s" % (v, k) for k, v in sorted(c.items()))),
           "",
           "| Module | What it is for |", "|---|---|"]
    grouped = by_module(entries)
    for short in MODULE_ORDER:
        if short in grouped:
            out.append("| [`%s`](#%s) | %s |"
                       % (short, short, MODULE_BLURB.get(short, "")))
    for short in MODULE_ORDER:
        group = grouped.get(short)
        if not group:
            continue
        out += ["", "## %s" % short, "", MODULE_BLURB.get(short, ""), "",
                "| Name | Kind | Summary |", "|---|---|---|"]
        for e in group:
            name = "`%s%s`" % (e.qualname, e.signature)
            if e.kind in ("class", "property"):
                name = "`%s`" % e.qualname
            out.append("| %s | %s | %s |"
                       % (name, e.kind, e.summary.replace("|", "\\|")))
    return "\n".join(out) + "\n"


def main(argv: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    query = ""
    if "--search" in argv:
        i = argv.index("--search")
        query = argv[i + 1] if i + 1 < len(argv) else ""
    entries = collect(include_private="--private" in argv)
    if query:
        entries = [e for e in entries if e.matches(query)]
    if "--markdown" in argv:
        path = "FUNCTIONS.md"
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(as_markdown(entries))
        print("wrote %s  (%d entries)" % (path, len(entries)))
        return 0
    print(as_text(entries))
    return 0


if __name__ == "__main__":
    sys.exit(main())
