"""The optimiser as a web page, served from a standard-library HTTP server.

    python dbf.py --web                  http://127.0.0.1:8000
    python dbf.py --web --host 0.0.0.0   reachable from your team's laptops
    python dbf.py --web --port 9000

Same model, same charts, different front end.  The page is one HTML
document with a few hundred lines of vanilla JavaScript; the physics stays
in Python and the browser only renders what the server sends back.

Why server-rendered SVG
-----------------------
Charts come back as SVG built by `plots.chart_svg`, which goes through the
exact same `Mapper` the desktop GUI uses to stroke a tkinter Canvas.  So a
chart is identical in both front ends, there is no JavaScript charting
library to add as a dependency, and there is only one place a chart can be
wrong.  The browser does `innerHTML = svg` and nothing else.

Scope, honestly
---------------
This serves on your own machine.  `--host 0.0.0.0` puts it on your local
network, which is enough for a design review round a table.  It is **not**
hardened for the public internet: there is no authentication, no TLS and
no rate limiting, and a search will happily use every core on the machine.
See the README for what deploying it publicly would actually take.
"""
from __future__ import annotations

import json
import threading
import traceback
from dataclasses import fields
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any, Dict, List, Optional, Tuple

from . import plots
from .aircraft import Aircraft, Design, opening_design
from .catalogue import as_text as catalogue_text, collect as catalogue_collect
from .components import AIRFOILS, MOTORS, PACKS, PROPS
from .config import (CONDITIONS, DEFAULT_CONDITIONS, SITE_NAME, Config,
                     Environment)
from .optimizer import DesignSpace, SearchSettings, optimise
from .report import print_bom, print_design_report, print_pareto, print_ranking
from .scoring import evaluate_scores, score_explanation
from .units import KG2LB

MAX_BODY = 1 << 20          # a design request is a few hundred bytes


def json_safe(obj: Any) -> Any:
    """Replace values that are valid Python floats but not valid JSON.

    Infinity becomes the string "inf", which `float()` reads straight back
    on the way in, so an unlimited field round-trips through the form
    unchanged.  NaN becomes null, because nothing sensible reads it back.
    """
    if isinstance(obj, float):
        if obj != obj:
            return None
        if obj == float("inf"):
            return "inf"
        if obj == float("-inf"):
            return "-inf"
        return obj
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_safe(v) for v in obj]
    return obj


# ==========================================================================
# Turning a request into a Config and a Design
# ==========================================================================
def config_from(payload: Dict[str, Any]) -> Config:
    """Build a Config from the form, ignoring anything we do not know.

    Unknown keys are skipped rather than trusted: this is the boundary
    where browser input meets the model, and `setattr` on an arbitrary
    name would let a request invent attributes the model never declared.
    """
    cfg = Config()
    env = payload.get("env", {})
    built = Environment.at_density_altitude(
        float(env.get("altitude_m", 0.0)),
        float(env.get("temp_c", 15.0)),
        float(env.get("wind_mps", 4.5)))
    # An explicit density overrides the one derived from elevation and
    # temperature; viscosity still follows the temperature.
    try:
        rho = float(env.get("rho") or 0.0)
    except (TypeError, ValueError):
        rho = 0.0
    cfg.env = (Environment(rho=rho, nu=built.nu, wind_mps=built.wind_mps,
                           field_elev_m=built.field_elev_m,
                           temp_c=built.temp_c)
               if rho > 0 else built)
    for owner in ("rules", "course", "build"):
        target = getattr(cfg, owner)
        allowed = {f.name: f for f in fields(target)}
        for name, value in (payload.get(owner) or {}).items():
            f = allowed.get(name)
            if f is None:
                continue
            current = getattr(target, name)
            if isinstance(current, bool) or not isinstance(current,
                                                           (int, float)):
                continue
            try:
                setattr(target, name,
                        int(float(value)) if isinstance(current, int)
                        else float(value))
            except (TypeError, ValueError):
                continue
    return cfg


def space_from(payload: Dict[str, Any]) -> DesignSpace:
    """The bounds a search may explore, from the form.

    These are NOT the design: `design` is one aeroplane that /evaluate
    scores exactly, while `space` is the box /search roams inside.  The
    two were conflated before - a search silently used a default
    DesignSpace and ignored everything the user had typed.
    """
    sp = DesignSpace()
    src = payload.get("space") or {}

    def pair(name, lo_default, hi_default):
        try:
            lo = float(src.get(name + "_lo", lo_default))
            hi = float(src.get(name + "_hi", hi_default))
        except (TypeError, ValueError):
            return (lo_default, hi_default)
        return (min(lo, hi), max(lo, hi))

    sp.span_bounds = pair("span", *sp.span_bounds)
    sp.ar_bounds = pair("ar", *sp.ar_bounds)
    sp.taper_bounds = pair("taper", *sp.taper_bounds)
    sp.m3_frac_bounds = pair("m3_frac", *sp.m3_frac_bounds)
    sp.drop_height_bounds = pair("drop", *sp.drop_height_bounds)
    r = Config().rules
    sp.sensor_bounds = pair("sensor", r.sensor_mass_min_kg,
                            r.sensor_mass_max_kg)

    def ints(name, fallback):
        out = []
        for part in str(src.get(name, "")).split(","):
            part = part.strip()
            if part:
                try:
                    out.append(int(float(part)))
                except ValueError:
                    pass
        return out or fallback

    sp.motor_counts = ints("motor_counts", sp.motor_counts)
    sp.system_counts = ints("system_counts", sp.system_counts)
    return sp


def design_from(payload: Dict[str, Any]) -> Design:
    d = payload.get("design", {})
    span = float(d.get("span_m", 1.55))
    ar = max(float(d.get("aspect_ratio", 6.0)), 0.5)
    return Design(
        span_m=span, chord_m=span / ar,
        airfoil=str(d.get("airfoil", "SD7062")),
        motor=str(d.get("motor", "SCORPION_SII_4025_440")),
        prop=str(d.get("prop", "APC_12x8")),
        n_motors=int(d.get("n_motors", 1)),
        pack=str(d.get("pack", "GA_6S3300_60C")),
        sensor_mass_kg=float(d.get("sensor_mass_kg", 0.9)),
        n_simulators=int(d.get("n_simulators", 0)),
        sensor_m3_kg=float(d.get("sensor_m3_kg", 0.0) or 0.0),
        drop_height_in=float(d.get("drop_height_in", 60.0)),
        taper=float(d.get("taper", 1.0)),
        n_systems=int(d.get("n_systems", 1)),
        wing_prop_layout=str(d.get("wing_prop_layout", "mid") or "mid"),
        winglet_h_m=float(d.get("winglet_h_m", 0.0) or 0.0),
        tail_vol_h=float(d.get("tail_vol_h", 0.0) or 0.0))


def capture(fn, *args, **kw) -> str:
    import io
    from contextlib import redirect_stdout
    buf = io.StringIO()
    with redirect_stdout(buf):
        fn(*args, **kw)
    return buf.getvalue()


# ==========================================================================
# The API
# ==========================================================================
def _opening_payload() -> Dict[str, Any]:
    """The opening design as form values: every Design field, plus the
    aspect ratio the form edits instead of chord, plus its source."""
    d, src = opening_design()
    out: Dict[str, Any] = {f.name: getattr(d, f.name) for f in fields(Design)}
    out["aspect_ratio"] = d.span_m / d.chord_m
    out["source"] = src
    return out


def api_meta() -> Dict[str, Any]:
    """Everything the page needs to build its form: catalogues and every
    editable field on Rules, Course and BuildStandard."""
    cfg = Config()
    sp0 = DesignSpace()
    groups = []
    for owner, title in (("rules", "Rules and limits"),
                         ("course", "Course geometry"),
                         ("build", "Build standard and modelling factors")):
        obj = getattr(cfg, owner)
        items = []
        for f in fields(obj):
            v = getattr(obj, f.name)
            if isinstance(v, bool) or not isinstance(v, (int, float)):
                continue
            items.append({"name": f.name, "value": v,
                          "int": isinstance(v, int)})
        groups.append({"owner": owner, "title": title, "fields": items})
    return {
        "airfoils": sorted(AIRFOILS),
        "motors": sorted(MOTORS),
        "props": sorted(PROPS),
        "packs": [{"key": k, "label": "%s  -  %.1f Wh, %.0f A, %.0f g"
                   % (p.name, p.watt_hours, p.i_max_a, p.mass_kg * 1000)}
                  for k, p in sorted(PACKS.items())],
        "space": {"span": list(sp0.span_bounds), "ar": list(sp0.ar_bounds),
                  "taper": list(sp0.taper_bounds),
                  "m3_frac": list(sp0.m3_frac_bounds),
                  "drop": list(sp0.drop_height_bounds),
                  "sensor": [cfg.rules.sensor_mass_min_kg,
                             cfg.rules.sensor_mass_max_kg],
                  "motor_counts": ",".join(str(v) for v in sp0.motor_counts),
                  "system_counts": ",".join(str(v)
                                            for v in sp0.system_counts),
                  "catalogue": [len(sp0.airfoils), len(sp0.motors),
                                len(sp0.props), len(sp0.packs)]},
        "site": SITE_NAME,
        # The same opening design the desktop GUI uses, so the two front
        # ends start on the same aeroplane.  The form used to hard-code a
        # design from several revisions ago.
        "opening": _opening_payload(),
        "conditions": {k: list(v) for k, v in sorted(CONDITIONS.items())},
        "default_conditions": DEFAULT_CONDITIONS,
        "advanced": groups,
        "charts": [n for n, _ in plots.CHART_BUILDERS] + ["Speeds vs weight"],
    }


def api_evaluate(payload: Dict[str, Any]) -> Dict[str, Any]:
    cfg = config_from(payload)
    design = design_from(payload)
    # A single design has nothing to normalise against, so it is always
    # scored on the field estimates - same rule as the desktop GUI.
    cfg.rules.normalisation = "field"
    ac = Aircraft(design, cfg)
    checks = ac.static_checks()
    card = evaluate_scores(ac, cfg, 1.5)

    imperial = bool(payload.get("imperial", True))
    name = payload.get("chart") or plots.CHART_BUILDERS[0][0]
    try:
        if name == "Speeds vs weight":
            chart = plots.weight_sweep(ac, cfg, imperial=imperial)
        else:
            chart = dict(plots.CHART_BUILDERS)[name](ac, imperial=imperial)
        chart_svg = plots.chart_svg(chart)
        series = [{"label": s.label, "xs": s.xs, "ys": s.ys}
                  for s in chart.series]
        axes = {"x": chart.x_label, "y": chart.y_label, "title": chart.title}
    except Exception as exc:
        chart_svg = '<p class="err">chart failed: %s</p>' % exc
        series, axes = [], {}

    return {
        "ok": bool(checks["ok"]),
        "failures": checks["failures"],
        "structure_ok": bool(ac.structure_ok),
        "report": capture(print_design_report, design, cfg, 1.5,
                          "DESIGN EVALUATION"),
        "bom": capture(print_bom, design, cfg),
        "score_blocks": score_explanation(ac, cfg, card),
        "dimensions": plots.dimensions(ac),
        "chart_svg": chart_svg,
        "chart_series": series,
        "chart_axes": axes,
        "three_view_svg": plots.three_view_svg(ac),
        "env": {"rho": cfg.env.rho,
                "density_altitude_ft": cfg.env.density_altitude_ft,
                "isa_pct": 100.0 * cfg.env.rho / 1.225},
        "summary": {
            "m1": card.m1, "m2": card.m2, "m3": card.m3,
            "total": card.total_mission,
            "mtow_lb": ac.mtow_m2_kg * KG2LB,
            "watt_hours": ac.total_watt_hours,
            "span_in": ac.wing.span_m * 39.3701,
        },
    }


class SearchJob:
    """One background search, polled by the page."""

    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.thread: Optional[threading.Thread] = None
        self.state = "idle"
        self.log = ""
        self.best: Optional[Dict[str, Any]] = None

    def running(self) -> bool:
        return self.thread is not None and self.thread.is_alive()

    def start(self, cfg: Config, mode: str,
              space: Optional[DesignSpace] = None) -> bool:
        if self.running():
            return False
        settings = (SearchSettings.quick() if mode == "quick" else
                    SearchSettings.thorough() if mode == "thorough" else
                    SearchSettings())
        with self.lock:
            self.state = "running"
            self.log = ("running the %s search - this takes minutes and "
                        "uses every core\n" % mode)
            self.best = None

        def work() -> None:
            try:
                # Capture the progress log AND keep the result from the same
                # run - calling optimise twice would double a 20 minute
                # search for nothing.
                import io
                from contextlib import redirect_stdout
                buf = io.StringIO()
                with redirect_stdout(buf):
                    res = optimise(cfg, space or DesignSpace(), settings,
                                   verbose=True)
                text = buf.getvalue()
                text += capture(print_ranking, res.ranked, 12,
                                "TOP DESIGNS (%s-normalised)" % res.mode)
                text += capture(print_pareto, res.pareto)
                # best design flyable on every weather preset
                best = ((res.best_row.get("design_robust")
                         or res.best_row["design"]) if res.best_row else None)
                with self.lock:
                    self.log = text
                    self.state = "done"
                    self.best = ({f.name: getattr(best, f.name)
                                  for f in fields(Design)} if best else None)
            except Exception:
                with self.lock:
                    self.log = traceback.format_exc()
                    self.state = "error"

        self.thread = threading.Thread(target=work, daemon=True)
        self.thread.start()
        return True

    def status(self) -> Dict[str, Any]:
        with self.lock:
            return {"state": self.state, "log": self.log, "best": self.best}


JOB = SearchJob()
_CATALOGUE = None


def api_catalogue(query: str = "") -> Dict[str, Any]:
    global _CATALOGUE
    if _CATALOGUE is None:
        _CATALOGUE = catalogue_collect()
    entries = [e for e in _CATALOGUE if e.matches(query)]
    return {"text": catalogue_text(entries), "count": len(entries),
            "total": len(_CATALOGUE)}


# ==========================================================================
# HTTP plumbing
# ==========================================================================
class Handler(BaseHTTPRequestHandler):
    server_version = "dbf2027"

    def log_message(self, fmt, *args):       # quieter than the default
        if self.path.startswith("/api/search"):
            return
        print("  %s %s" % (self.command, self.path))

    # ---- helpers ----
    def _send(self, code: int, body: bytes, ctype: str) -> None:
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        try:
            self.wfile.write(body)
        except (BrokenPipeError, ConnectionAbortedError):
            pass                              # the tab was closed mid-reply

    def _json(self, obj: Any, code: int = 200) -> None:
        # allow_nan=False makes this fail loudly rather than emitting bare
        # Infinity/NaN, which Python's json writes happily and every
        # browser rejects as invalid JSON.  `takeoff_distance_m` is
        # genuinely infinite - there is no takeoff limit this year - so
        # this is not hypothetical; it broke the page on first load.
        body = json.dumps(json_safe(obj), allow_nan=False)
        self._send(code, body.encode("utf-8"),
                   "application/json; charset=utf-8")

    def _body(self) -> Dict[str, Any]:
        n = int(self.headers.get("Content-Length") or 0)
        if n <= 0 or n > MAX_BODY:
            return {}
        try:
            return json.loads(self.rfile.read(n).decode("utf-8"))
        except (ValueError, UnicodeDecodeError):
            return {}

    # ---- routes ----
    def do_GET(self) -> None:
        path = self.path.split("?", 1)[0]
        if path in ("/", "/index.html"):
            self._send(200, PAGE.encode("utf-8"), "text/html; charset=utf-8")
        elif path == "/api/meta":
            self._json(api_meta())
        elif path == "/api/search":
            self._json(JOB.status())
        elif path == "/api/rules":
            from .rulebook import as_text as rules_text
            self._json({"text": rules_text()})
        elif path == "/api/catalogue":
            q = ""
            if "?" in self.path:
                from urllib.parse import parse_qs
                q = (parse_qs(self.path.split("?", 1)[1]).get("q")
                     or [""])[0]
            self._json(api_catalogue(q))
        else:
            self._json({"error": "not found"}, 404)

    def do_POST(self) -> None:
        path = self.path.split("?", 1)[0]
        try:
            if path == "/api/evaluate":
                self._json(api_evaluate(self._body()))
            elif path == "/api/search":
                body = self._body()
                started = JOB.start(config_from(body),
                                    body.get("mode", "quick"),
                                    space_from(body))
                self._json({"started": started, **JOB.status()})
            else:
                self._json({"error": "not found"}, 404)
        except Exception as exc:
            # A bad design is a normal outcome, not a server fault: report
            # it as data so the page can show it next to the form.
            self._json({"ok": False, "failures": ["%s: %s"
                                                  % (type(exc).__name__, exc)],
                        "report": traceback.format_exc()})


def serve(host: str = "127.0.0.1", port: int = 8000,
          open_browser: bool = True) -> int:
    httpd = ThreadingHTTPServer((host, port), Handler)
    shown = "127.0.0.1" if host in ("0.0.0.0", "") else host
    url = "http://%s:%d" % (shown, port)
    print("=" * 70)
    print("DBF 2026-27 optimiser - web front end")
    print("=" * 70)
    print("  %s" % url)
    if host == "0.0.0.0":
        print("  also reachable from other machines on this network")
        print("  NOTE: no authentication and no TLS - keep it on a network")
        print("        you trust.")
    print("  Ctrl-C to stop.")
    print()
    if open_browser:
        import webbrowser
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")
    finally:
        httpd.server_close()
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    host, port, open_browser = "127.0.0.1", 8000, True
    if "--host" in argv:
        host = argv[argv.index("--host") + 1]
    if "--port" in argv:
        port = int(argv[argv.index("--port") + 1])
    if "--no-browser" in argv:
        open_browser = False
    return serve(host, port, open_browser)


# The page itself lives in its own module so this one stays readable.
from .webpage import PAGE           # noqa: E402  (circular-free, data only)


if __name__ == "__main__":
    raise SystemExit(main())
