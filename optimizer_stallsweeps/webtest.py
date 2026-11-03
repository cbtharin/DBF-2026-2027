#!/usr/bin/env python3
"""Drive the web front end over real HTTP, with no browser.

    python webtest.py      (or: python dbf.py --webtest)

Starts the actual server on a spare port, then talks to it the way the
page does: fetch the metadata, post a design, read back the report, the
score blocks, the dimensions and the SVG.  Checks the numbers that come
out of the socket match the model, not just that a request succeeded.
"""
from __future__ import annotations

import json
import socket
import sys
import threading
import traceback
import urllib.error
import urllib.request

PASS, FAIL = 0, 0


def check(name: str, ok: bool, detail: str = "") -> None:
    global PASS, FAIL
    if ok:
        PASS += 1
        print("  PASS  %s" % name)
    else:
        FAIL += 1
        print("  FAIL  %s%s" % (name, ("   [%s]" % detail) if detail else ""))


def section(name: str) -> None:
    print("\n%s\n%s" % (name, "-" * len(name)))


def free_port() -> int:
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    return port


def main() -> int:
    from http.server import ThreadingHTTPServer
    from dbf2027 import webapp
    from dbf2027.aircraft import Aircraft, Design
    from dbf2027.scoring import evaluate_scores

    port = free_port()
    httpd = ThreadingHTTPServer(("127.0.0.1", port), webapp.Handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = "http://127.0.0.1:%d" % port

    def get(path: str):
        with urllib.request.urlopen(base + path, timeout=60) as r:
            return r.status, r.read().decode("utf-8")

    def post(path: str, obj):
        req = urllib.request.Request(
            base + path, data=json.dumps(obj).encode("utf-8"),
            headers={"Content-Type": "application/json"}, method="POST")
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.status, json.loads(r.read().decode("utf-8"))

    try:
        # --------------------------------------------------------------
        section("A. The page itself")
        status, html = get("/")
        check("GET / returns 200", status == 200, str(status))
        check("it is a complete HTML document",
              html.startswith("<!DOCTYPE html>") and "</html>" in html)
        for needed in ("id=\"go\"", "id=\"chart\"", "id=\"advanced\"",
                       "/api/evaluate", "/api/meta", "/api/search",
                       "/api/catalogue"):
            check("the page wires up %s" % needed, needed in html)
        # A container id that collides with a field id makes $("#x") return
        # the container; .value is undefined, the key vanishes from the
        # POST and the server quietly uses its own default.  That shipped
        # once: the Propulsion panel was id="prop", so the propeller on
        # screen was not the propeller being scored.
        import re as _re
        ids = [i for i in _re.findall(r'id="([^"]+)"', html)
               if "${" not in i]          # skip JS template placeholders
        dupes = sorted({i for i in ids if ids.count(i) > 1})
        check("no duplicate element ids in the page", not dupes, str(dupes))
        field_ids = ["airfoil", "span_m", "aspect_ratio", "taper", "motor",
                     "prop", "n_motors", "pack", "n_systems",
                     "sensor_mass_kg", "n_simulators", "sensor_m3_kg",
                     "drop_height_in"]
        check("no container shadows a design field id",
              not [f for f in field_ids if ('id="%s"' % f) in html],
              str([f for f in field_ids if ('id="%s"' % f) in html]))
        check("the form reports a dropped field instead of silently "
              "sending a partial design",
              "is missing or shadowed" in html)

        check("no external resources are loaded",
              "http://" not in html.replace("http://www.w3.org", "")
              and "cdn" not in html.lower())

        # --------------------------------------------------------------
        section("B. Metadata drives the form")
        status, body = get("/api/meta")
        meta = json.loads(body)
        check("GET /api/meta returns 200", status == 200)
        from dbf2027 import components as _C
        check("every catalogue is offered",
              len(meta["airfoils"]) == len(_C.AIRFOILS) == 32
              and len(meta["motors"]) == len(_C.MOTORS)
              and len(meta["props"]) == len(_C.PROPS)
              and len(meta["packs"]) == len(_C.PACKS),
              "%d/%d/%d/%d" % (len(meta["airfoils"]), len(meta["motors"]),
                               len(meta["props"]), len(meta["packs"])))
        owners = {g["owner"] for g in meta["advanced"]}
        check("rules, course and build are all editable",
              owners == {"rules", "course", "build"}, str(owners))
        names = {f["name"] for g in meta["advanced"] for f in g["fields"]}
        for need in ("pattern_altitude_m", "wing_areal_density",
                     "fuse_areal_density", "h_tail_areal_density",
                     "v_tail_areal_density", "max_watt_hours",
                     "max_fuse_current_a", "max_gross_mass_kg"):
            check("%s is exposed" % need, need in names)
        check("five charts are offered", len(meta["charts"]) == 5)
        # Python writes bare Infinity/NaN, which browsers reject outright.
        # takeoff_distance_m really is infinite, so this is a live risk.
        check("the JSON is strictly valid - no Infinity or NaN",
              "Infinity" not in body and "NaN" not in body)
        json.loads(body, parse_constant=lambda c: (_ for _ in ()).throw(
            ValueError("bare %s in JSON" % c)))
        check("a strict parser accepts it", True)
        inf_fields = [f for g in meta["advanced"] for f in g["fields"]
                      if f["value"] == "inf"]
        check("an unlimited field is served as \"inf\"", bool(inf_fields),
              str([f["name"] for f in inf_fields]))
        check("pack labels carry the Wh figure",
              all("Wh" in p["label"] for p in meta["packs"]))

        # The form must open on the same design as the desktop GUI, not a
        # design typed into the page that goes stale as the optimum moves.
        with open("examples/optimised_v1.json", encoding="utf-8") as fh:
            rec = json.load(fh)
        op = meta.get("opening", {})
        lost = [k for k, v in rec.items()
                if k not in op or (op[k] != v if isinstance(v, str)
                                   else abs(float(op[k]) - float(v)) > 1e-9)]
        check("the page opens on the recommended design, every field",
              not lost, "lost or changed: %s" % lost)
        check("including the optimised taper",
              abs(float(op.get("taper", 1.0)) - float(rec["taper"])) < 1e-9,
              str(op.get("taper")))
        check("and it names where that design came from",
              "optimised_v1.json" in str(op.get("source", "")))
        _, page = get("/")
        check("the form reads the opening design, not typed-in values",
              "m.opening" in page
              and 'row("airfoil","airfoil","MH32"' not in page)

        # --------------------------------------------------------------
        section("C. Evaluating a design over HTTP")
        # The design the page opens on, so this section exercises what a
        # user sees first.  The old hand-typed fixture (MH32, Hacker A50)
        # stopped completing Mission 2 on the airfoiltools polars.
        design = {k: str(op[k]) for k in (
            "airfoil", "span_m", "aspect_ratio", "taper", "motor", "prop",
            "n_motors", "pack", "n_systems", "sensor_mass_kg",
            "n_simulators", "sensor_m3_kg", "drop_height_in")}
        # The contest site, which is also what the form defaults to.
        payload = {"design": design,
                   "env": {"wind_mps": "5.0", "altitude_m": "667",
                           "temp_c": "28"},
                   "rules": {}, "course": {}, "build": {},
                   "imperial": True, "chart": "Drag and thrust vs speed"}
        status, j = post("/api/evaluate", payload)
        check("POST /api/evaluate returns 200", status == 200)
        check("the design is feasible", j["ok"] is True, str(j["failures"]))
        check("the report came back", "DESIGN EVALUATION" in j["report"])
        check("the parts list came back", "BILL OF MATERIALS" in j["bom"])
        check("the parts list has no prices", "$" not in j["bom"])
        check("dimensions came back", len(j["dimensions"]) > 20,
              "%d rows" % len(j["dimensions"]))
        check("the chart is SVG",
              j["chart_svg"].startswith("<svg") and "<polyline" in j["chart_svg"])
        check("the three-view is SVG",
              j["three_view_svg"].startswith("<svg")
              and "<polygon" in j["three_view_svg"])
        # The axis label and the footnote both sit below the plot frame;
        # at the old default height they were drawn on top of each other.
        import re as _re2
        ys = [float(v) for v in
              _re2.findall(r'<text[^>]* y="([0-9.]+)"', j["chart_svg"])]
        check("nothing is drawn outside the SVG canvas",
              ys and max(ys) <= 560, "lowest text at y=%.0f" % max(ys))
        check("the axis label and footnote do not collide",
              len({round(y) for y in ys if y > 520}) == len(
                  [y for y in ys if y > 520]),
              str(sorted(round(y) for y in ys if y > 520)))

        check("chart series came back as data too",
              len(j["chart_series"]) >= 3
              and len(j["chart_series"][0]["xs"]) > 20)

        # the numbers on the wire must match the model, not be re-derived
        from dbf2027.config import Config
        cfg = webapp.config_from(payload)
        cfg.rules.normalisation = "field"
        ac = Aircraft(webapp.design_from(payload), cfg)
        card = evaluate_scores(ac, cfg, 1.5)
        s = j["summary"]
        check("the served total matches the model",
              abs(s["total"] - card.total_mission) < 1e-9,
              "%.4f vs %.4f" % (s["total"], card.total_mission))
        check("the served gross weight matches the model",
              abs(s["mtow_lb"] - ac.mtow_m2_kg * 2.2046226) < 1e-4)
        from dbf2027.config import (CONDITIONS as _COND,
                                    DEFAULT_CONDITIONS as _DC)
        check("meta names the contest site",
              "Tucson" in meta.get("site", ""), meta.get("site"))
        check("meta offers the named conditions",
              set(meta["conditions"]) == set(_COND)
              and meta["default_conditions"] == _DC)
        check("the response reports the air density it flew",
              "env" in j and 1.0 < j["env"]["rho"] < 1.15,
              str(j.get("env")))
        check("and the site costs thousands of feet of density altitude",
              j["env"]["density_altitude_ft"] > 3000.0,
              "%.0f ft" % j["env"]["density_altitude_ft"])
        thin = json.loads(json.dumps(payload))
        thin["env"]["rho"] = "0.95"
        _, jt = post("/api/evaluate", thin)
        check("a posted air density overrides elevation and temperature",
              abs(jt["env"]["rho"] - 0.95) < 1e-9, str(jt["env"]["rho"]))
        check("thinner air is a higher density altitude",
              jt["env"]["density_altitude_ft"]
              > j["env"]["density_altitude_ft"] + 1000)
        # The form is built by JS at runtime, so the served document
        # carries the row() call rather than the rendered input.
        check("the page builds an air-density box",
              'row("rho"' in html and "rhonote" in html)

        check("the served energy matches the model",
              abs(s["watt_hours"] - ac.total_watt_hours) < 1e-9)

        # --------------------------------------------------------------
        section("D. The scoring panel is served as structured data")
        blocks = j["score_blocks"]
        check("five blocks: M1, M2, M3, Ground Mission, total",
              len(blocks) == 5, "%d" % len(blocks))
        titles = " ".join(b["title"] for b in blocks)
        check("each mission is named",
              "Mission 1" in titles and "Mission 2" in titles
              and "Mission 3" in titles)
        formulas = " ".join(b["formula"] for b in blocks)
        check("M2's formula is shown",
              "M2 = 1 + N(weight/time) / Max(weight/time)" in formulas)
        check("M3's formula is shown",
              "M3 = 2 + N(#laps x sensor weight)" in formulas)
        check("every block has rows and a result",
              all(b["rows"] and b["result"] for b in blocks))
        check("the Ground Mission is scored and shown",
              any("Ground Mission" in b["title"] for b in blocks))
        check("GM's formula is shown",
              any("GM = 0.5 + N(sensor weight x drop height)" in b["formula"]
                  for b in blocks))
        check("the last block totals all four events",
              blocks[-1]["result"].startswith("TOTAL =")
              and "M1 + M2 + M3 + GM" in blocks[-1]["formula"])

        # --------------------------------------------------------------
        section("E. Advanced values sent from the form reach the model")
        p2 = json.loads(json.dumps(payload))
        p2["rules"]["max_watt_hours"] = "40"
        _, j2 = post("/api/evaluate", p2)
        check("a 40 Wh cap makes this 99.9 Wh design infeasible",
              j2["ok"] is False
              and any("Wh" in f for f in j2["failures"]),
              str(j2["failures"])[:70])

        # The recommended design sits on the energy limit and cannot
        # afford a 45 m climb, so compare two altitudes at a payload with
        # margin.  Comparing failed flights compared the times the battery
        # ran out, which is how this once "passed" backwards.
        p_lo = json.loads(json.dumps(payload))
        for k in ("sensor_mass_kg", "sensor_m3_kg"):
            p_lo["design"][k] = "7.5"
        p3 = json.loads(json.dumps(p_lo))
        # 90 m: the faster climb covers ground towards the first lap, so
        # the time a higher pattern costs is small - 0.07 s for 15 m more,
        # under the 0.1 s the page shows.
        p3["course"]["pattern_altitude_m"] = "90"
        _, j_lo = post("/api/evaluate", p_lo)
        _, j3 = post("/api/evaluate", p3)
        t_low = [b for b in j_lo["score_blocks"] if "Mission 2" in b["title"]][0]
        t_high = [b for b in j3["score_blocks"] if "Mission 2" in b["title"]][0]

        def lap_time(block):
            for label, value, _ in block["rows"]:
                if label.startswith("time for"):
                    return float(value.split()[0])
            return -1.0
        check("both altitude runs complete Mission 2, so there is a time "
              "to compare", lap_time(t_low) > 0 and lap_time(t_high) > 0,
              "%.2f / %.2f s" % (lap_time(t_low), lap_time(t_high)))
        check("raising the pattern altitude costs time",
              lap_time(t_high) > lap_time(t_low),
              "%.2f vs %.2f s" % (lap_time(t_high), lap_time(t_low)))

        p4 = json.loads(json.dumps(payload))
        p4["build"]["h_tail_areal_density"] = "3.60"
        _, j4 = post("/api/evaluate", p4)

        def dim(rows, key):
            for r in rows:
                if r[0] == key:
                    return r[1]
            return ""
        check("doubling elevator density raises the empty mass",
              dim(j4["dimensions"], "Empty mass")
              != dim(j["dimensions"], "Empty mass"),
              "%s -> %s" % (dim(j["dimensions"], "Empty mass"),
                            dim(j4["dimensions"], "Empty mass")))

        # a value the model never declared must be ignored, not set
        p5 = json.loads(json.dumps(payload))
        p5["rules"]["not_a_real_field"] = "1"
        p5["build"]["__class__"] = "evil"
        status5, j5 = post("/api/evaluate", p5)
        check("unknown config keys are ignored, not trusted",
              status5 == 200 and j5["ok"] is True)

        # --------------------------------------------------------------
        section("F. Charts render for every option and both unit systems")
        for name in meta["charts"]:
            for imp in (True, False):
                pc = json.loads(json.dumps(payload))
                pc["chart"], pc["imperial"] = name, imp
                _, jc = post("/api/evaluate", pc)
                svg = jc["chart_svg"]
                check("%s renders (%s)" % (name, "imp" if imp else "SI"),
                      svg.startswith("<svg") and svg.count("<polyline") >= 1,
                      svg[:60])

        # --------------------------------------------------------------
        section("G. Design vs search space are distinct over HTTP")
        check("meta carries the search-space defaults",
              set(meta["space"]) >= {"span", "ar", "taper", "sensor",
                                     "m3_frac", "drop", "motor_counts",
                                     "system_counts"})
        base_sp = webapp.space_from({})
        narrowed = webapp.space_from({"space": {
            "span_lo": "1.20", "span_hi": "1.35",
            "motor_counts": "1", "system_counts": "1",
            "drop_lo": "20", "drop_hi": "30"}})
        from dbf2027.config import Config as _C
        b = narrowed.continuous_bounds(_C())
        check("a posted span range reaches the optimiser",
              abs(b[0][0] - 1.20) < 1e-9 and abs(b[0][1] - 1.35) < 1e-9,
              str(b[0]))
        check("posted motor counts reach the optimiser",
              narrowed.motor_counts == [1] and narrowed.system_counts == [1])
        check("posted drop-height bounds reach the optimiser",
              abs(b[5][0] - 20) < 1e-9 and abs(b[5][1] - 30) < 1e-9)
        check("an inverted posted range is repaired",
              webapp.space_from({"space": {"span_lo": "2", "span_hi": "1"}}
                                ).span_bounds == (1.0, 2.0))
        check("garbage in a bound falls back to the default",
              webapp.space_from({"space": {"span_lo": "abc"}}).span_bounds
              == base_sp.span_bounds)
        check("the design does NOT come from the space",
              webapp.design_from(payload).span_m
              != narrowed.span_bounds[0])
        check("the page has a search-space card and posts it",
              'id="space"' in html and "out.space" in html)
        check("the page says the search ignores the design fields",
              "a search ignores them" in html)

        # --------------------------------------------------------------
        section("H. Rules traceability over HTTP")
        status, body = get("/api/rules")
        rules_txt = json.loads(body)["text"]
        check("the rules matrix is served", status == 200
              and len(rules_txt) > 3000, "%d chars" % len(rules_txt))
        check("annotations and rulebook agree",
              "Annotations and rulebook agree" in rules_txt)
        for rid in ("3.1.3a", "3.3.4", "3.2.3.2d", "6.3"):
            check("rule %s appears" % rid, rid in rules_txt)
        check("it cites source locations",
              ".py:" in rules_txt and "enforced at" in rules_txt)
        check("the page has a Rules tab", "Rules" in html
              and "/api/rules" in html)

        # --------------------------------------------------------------
        section("I. Function catalogue over HTTP")
        status, body = get("/api/catalogue")
        cat = json.loads(body)
        check("the whole catalogue is served", cat["count"] > 150,
              "%d" % cat["count"])
        status, body = get("/api/catalogue?q=trim")
        cat2 = json.loads(body)
        check("filtering narrows it", 0 < cat2["count"] < cat["count"],
              "%d of %d" % (cat2["count"], cat["count"]))
        check("the filtered text mentions the trim solver",
              "tail_load" in cat2["text"])

        # --------------------------------------------------------------
        section("J. Errors and edges")
        bad = json.loads(json.dumps(payload))
        bad["design"]["span_m"] = "not a number"
        status, jb = post("/api/evaluate", bad)
        check("garbage input returns a message, not a 500",
              status == 200 and jb["ok"] is False, str(status))
        check("the failure is explained", bool(jb["failures"]))
        # A failed response must not carry a design's worth of panels, and
        # the page must CLEAR the ones it does not get - otherwise the
        # previous aeroplane stays on screen with nothing to say so.
        check("a failed response carries no stale dimensions",
              not jb.get("dimensions"))
        check("a failed response carries no stale chart",
              not jb.get("chart_svg"))
        check("a failed response carries no stale score blocks",
              not jb.get("score_blocks"))
        check("the page clears panels the response omits",
              'j.chart_svg ||' in html and 'else $("#dims")' in html)
        try:
            get("/api/nope")
            check("unknown routes 404", False)
        except urllib.error.HTTPError as e:
            check("unknown routes 404", e.code == 404, str(e.code))
        status, js = get("/api/search")
        check("the search endpoint reports idle before any run",
              json.loads(js)["state"] == "idle")

        print("\n" + "=" * 60)
        print("%d passed, %d failed" % (PASS, FAIL))
        print("=" * 60)
        return 1 if FAIL else 0
    finally:
        httpd.shutdown()
        httpd.server_close()


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc()
        sys.exit(2)
