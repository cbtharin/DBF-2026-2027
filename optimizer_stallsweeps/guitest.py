#!/usr/bin/env python3
"""Exercise the GUI without anyone sitting in front of it.

    python guitest.py

Builds the real window, drives it the way a user would - type in a field,
press Evaluate, switch charts, filter the catalogue, resize - and checks
what came out.  Everything runs through `update()` rather than
`mainloop()`, so the window opens, is driven, and closes on its own.

This exists because a GUI is the easiest place in a project for something
to rot silently: a renamed model attribute breaks a panel that nobody
notices until a demo.  Every check here is one of those.
"""
from __future__ import annotations

import json
import sys
import traceback

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


def text_of(widget) -> str:
    return widget.get("1.0", "end")


def main() -> int:
    try:
        import tkinter as tk
    except ImportError:
        print("tkinter is not available - GUI tests skipped")
        return 0

    from dataclasses import replace

    from dbf2027 import plots
    from dbf2027.aircraft import Design
    from dbf2027.components import PACKS
    from dbf2027.config import Config
    from dbf2027.gui import App

    root = tk.Tk()
    root.geometry("1400x900")
    app = App(root)
    root.update()

    # ------------------------------------------------------------------
    section("A. The window builds and evaluates its default design")
    check("an aircraft was built", app.aircraft is not None)
    check("the report panel has content", len(text_of(app.txt_report)) > 2000,
          "%d chars" % len(text_of(app.txt_report)))
    check("the report is the design dossier",
          "DESIGN EVALUATION" in text_of(app.txt_report))
    check("the parts list has content",
          "BILL OF MATERIALS" in text_of(app.txt_bom))
    check("the parts list carries no prices",
          "$" not in text_of(app.txt_bom))
    check("the dimensions table is populated",
          len(app.tv_dims.get_children()) > 20,
          "%d rows" % len(app.tv_dims.get_children()))
    check("the status bar reports a score", "total" in app.status.get(),
          app.status.get())
    check("a design is loaded at startup", app.design is not None)
    check("and its source is named",
          "optimised_v1.json" in app.design_source.cget("text")
          or "fallback" in app.design_source.cget("text"),
          app.design_source.cget("text")[:50])

    # The design values belong to the Report, not to the left panel: they
    # used to appear in both, with nothing saying which was authoritative.
    check("the left panel has no design widgets",
          not any(hasattr(app, n) for n in
                  ("f_span", "f_airfoil", "f_motor", "f_prop", "f_pack",
                   "f_sensor", "f_nsim", "f_ar", "f_taper", "f_nmotors",
                   "f_nsys", "f_m3sensor", "f_drop")),
          str([n for n in ("f_span", "f_airfoil", "f_motor", "f_sensor")
               if hasattr(app, n)]))
    check("it keeps only a one-line caption naming the design",
          app.design.airfoil in app.design_label.cget("text")
          and app.design.motor in app.design_label.cget("text"))
    check("the caption is a caption, not a readout",
          len(app.design_label.cget("text")) < 200,
          "%d chars" % len(app.design_label.cget("text")))
    check("the Report still states the design in full",
          all(t in text_of(app.txt_report)
              for t in ("airfoil", "span", "chord", "aspect ratio")))
    # Air density is what reaches the model; elevation and temperature
    # are only a way of arriving at it, so it is shown and editable.
    from dbf2027.config import Environment as _Env
    from dbf2027.config import CONDITIONS as _COND, DEFAULT_CONDITIONS as _DC
    check("the conditions preset defaults to the contest site",
          app.f_cond.get() == _DC, app.f_cond.get())
    check("and the boxes match that preset",
          (app.f_alt.get(), app.f_temp.get(), app.f_wind.get())
          == _COND[_DC])
    app.f_cond.set("tucson-hot"); app._apply_conditions()
    check("picking a preset loads all three boxes",
          (app.f_alt.get(), app.f_temp.get(), app.f_wind.get())
          == _COND["tucson-hot"])
    hot_rho = app.f_rho.get()
    app.f_cond.set("tucson-morning"); app._apply_conditions()
    check("a cool morning is denser than a hot afternoon",
          app.f_rho.get() > hot_rho,
          "%.4f vs %.4f" % (app.f_rho.get(), hot_rho))
    app.f_cond.set(_DC); app._apply_conditions()

    check("air density is on the conditions panel", hasattr(app, "f_rho"))
    check("it starts at the contest site density, not sea level",
          abs(app.f_rho.get() - _Env.at_site(_DC).rho) < 1e-3,
          "%.4f" % app.f_rho.get())
    check("and says what that means in density altitude",
          "density altitude" in app.rho_note.cget("text"),
          app.rho_note.cget("text")[:50])
    app.f_alt.set(1500.0); app.f_temp.set(35.0); app._sync_density()
    hot = app.f_rho.get()
    check("changing elevation and temperature recomputes it",
          hot < 1.0, "%.4f kg/m3" % hot)
    check("the derived value reaches the model",
          abs(app.read_config().env.rho - hot) < 1e-9)
    check("the note follows it", "8235 ft" in app.rho_note.cget("text")
          or "ft density altitude" in app.rho_note.cget("text"))
    app.f_rho.set(1.050)
    check("a typed density overrides elevation and temperature",
          abs(app.read_config().env.rho - 1.050) < 1e-9)
    check("but viscosity still follows the temperature",
          abs(app.read_config().env.nu
              - _Env.at_density_altitude(1500.0, 35.0).nu) < 1e-12)
    app.f_cond.set(_DC); app._apply_conditions()
    check("the preset restores the site conditions",
          abs(app.f_rho.get() - _Env.at_site(_DC).rho) < 1e-3)

    check("conditions and scoring stay on the left",
          all(hasattr(app, n) for n in ("f_wind", "f_alt", "f_temp",
                                        "f_norm", "f_fm2", "f_fm3",
                                        "f_fmg")))

    # ------------------------------------------------------------------
    section("B. Every advanced field maps to a real model attribute")
    cfg = Config()
    missing = [(o, n) for (o, n) in app.adv
               if not hasattr(getattr(cfg, o), n)]
    check("no advanced field is orphaned", not missing, str(missing[:3]))
    n_rules = sum(1 for o, _ in app.adv if o == "rules")
    n_build = sum(1 for o, _ in app.adv if o == "build")
    check("rules fields are exposed", n_rules >= 15, "%d" % n_rules)
    check("build fields are exposed", n_build >= 30, "%d" % n_build)
    # the build standard is grouped by kind, each field with its note
    check("the margins get their own group on the Advanced tab",
          ("build", "load_safety_factor") in app.adv
          and ("build", "m2_energy_margin") in app.adv
          and ("build", "excrescence_frac") in app.adv
          and ("build", "fuselage_lift_carryover") in app.adv)
    tops = [app.tv_margins.item(i, "text")
            for i in app.tv_margins.get_children()]
    check("the Margins tab shows margins, computed corrections, estimates",
          any("Safety margins" in t for t in tops)
          and any("Computed for this design" in t for t in tops)
          and any("measurements" in t for t in tops), str(tops[:3]))
    _comp = [c for c in app.tv_margins.get_children()
             if "Computed" in app.tv_margins.item(c, "text")][0]
    check("the computed corrections are listed with what they replaced",
          len(app.tv_margins.get_children(_comp)) >= 8
          and all(str(app.tv_margins.item(k, "values")[2]).startswith("was")
                  for k in app.tv_margins.get_children(_comp)))
    check("advanced values round-trip into a Config",
          app.read_config().rules.max_watt_hours == 100.0)

    # changing an advanced number must actually reach the model
    app.adv[("rules", "max_watt_hours")].set(40.0)
    cfg2 = app.read_config()
    check("an edited advanced field reaches the Config",
          cfg2.rules.max_watt_hours == 40.0)
    app.evaluate()
    check("a 40 Wh cap makes the 100 Wh design infeasible",
          "CONSTRAINTS VIOLATED" in text_of(app.txt_report))
    app.adv[("rules", "max_watt_hours")].set(100.0)
    app.evaluate()
    check("restoring the cap makes it legal again",
          "all satisfied" in text_of(app.txt_report))

    # ------------------------------------------------------------------
    section("C. The design is state, driven by search or file")
    base = app.design
    app.write_design(replace(base, span_m=1.40, chord_m=1.40 / 5.55),
                     "test")
    app.evaluate()
    check("setting the design changes the aeroplane",
          abs(app.aircraft.wing.span_m - 1.40) < 1e-9,
          "%.3f" % app.aircraft.wing.span_m)
    check("the caption follows the design",
          "1.400 m span" in app.design_label.cget("text"),
          app.design_label.cget("text")[:60])
    check("and names where it came from",
          "test" in app.design_source.cget("text"))

    app.write_design(replace(base, n_simulators=1), "test")
    app.evaluate()
    one = app.aircraft.m2_scored_kg
    loaded = app.aircraft.loaded_container_kg
    app.write_design(replace(base, n_simulators=2), "test")
    app.evaluate()
    two = app.aircraft.m2_scored_kg
    check("an extra simulator adds exactly one loaded container",
          abs((two - one) - loaded) < 1e-9,
          "%.4f vs %.4f kg" % (two - one, loaded))

    app.write_design(replace(base, n_motors=2, n_systems=2), "test")
    app.evaluate()
    two_sys = app.aircraft.total_watt_hours
    app.write_design(replace(base, n_motors=2, n_systems=1), "test")
    app.evaluate()
    one_sys = app.aircraft.total_watt_hours
    check("two propulsion systems carry twice the energy",
          abs(two_sys - 2 * one_sys) < 1e-6,
          "%.1f vs %.1f Wh" % (two_sys, one_sys))

    app.write_design(base, "restored")
    app.evaluate()
    check("the original design is restored", app.aircraft is not None)

    # ------------------------------------------------------------------
    section("D. Charts draw, for every chart and both unit systems")
    app.out_nb.select(2)              # the charts tab must be laid out
    root.update()
    names = [n for n, _ in plots.CHART_BUILDERS] + ["Speeds vs weight"]
    for name in names:
        for imp in (True, False):
            app.chart_name.set(name)
            app.imperial.set(imp)
            app._draw_chart()
            root.update()
            items = len(app.chart_canvas.find_all())
            check("%s draws (%s)" % (name, "imperial" if imp else "SI"),
                  items > 20 and app.chart_canvas.chart is not None,
                  "%d canvas items" % items)
    app.imperial.set(True)

    app.out_nb.select(3)
    root.update()
    app.three_view.redraw()
    root.update()
    check("the three-view drew something",
          len(app.three_view.find_all()) > 8,
          "%d items" % len(app.three_view.find_all()))

    # resize must not throw
    root.geometry("1100x700")
    root.update()
    app.chart_canvas.redraw()
    app.three_view.redraw()
    root.update()
    check("charts survive a resize", len(app.chart_canvas.find_all()) > 20)

    # ------------------------------------------------------------------
    section("E. Chart geometry is correct, not just non-empty")
    chart = plots.forces_vs_speed(app.aircraft)
    m = plots.Mapper(chart, 800, 500)
    x0, x1, y0, y1 = chart.bounds()
    check("the mapper puts the x range inside the plot frame",
          abs(m.px(x0) - m.left) < 1e-6 and abs(m.px(x1) - m.right) < 1e-6)
    check("the y axis is flipped for a canvas",
          m.py(y0) > m.py(y1))
    check("ticks land inside their range",
          all(x0 - 1e-9 <= t <= x1 + 1e-9
              for t in plots.nice_ticks(x0, x1)))
    check("ticks are evenly spaced",
          len({round(b - a, 9) for a, b in
               zip(plots.nice_ticks(0, 97), plots.nice_ticks(0, 97)[1:])}) == 1)
    check("a degenerate range does not divide by zero",
          plots.Mapper(plots.Chart("t", "x", "y"), 400, 300).px(0.5) > 0)
    check("thrust and drag are both plotted",
          any("thrust" in s.label for s in chart.series)
          and any("drag" in s.label for s in chart.series))
    check("drag curves are drawn for three weights",
          sum(1 for s in chart.series if "drag" in s.label) == 3)

    polar = plots.lift_vs_drag(app.aircraft)
    check("the polar has varying lift, not a flat line",
          max(polar.series[0].ys) > 3 * min(polar.series[0].ys))
    check("the polar marks both mission weights",
          len(polar.h_lines) == 2)

    ld = plots.lift_to_drag(app.aircraft)
    peak = max(ld.series[-1].ys)
    check("peak L/D is physically plausible", 4.0 < peak < 30.0,
          "%.1f" % peak)

    # heavier aeroplanes must stall faster - a real invariant, via the chart
    sweep = plots.weight_sweep(app.aircraft, app.cfg)
    stall = [s for s in sweep.series if "stall" in s.label][0]
    check("stall speed rises with weight",
          all(b >= a - 1e-9 for a, b in zip(stall.ys, stall.ys[1:])))

    # ------------------------------------------------------------------
    section("F. Three-view geometry is to scale")
    plan, side, ext = plots.three_view(app.aircraft)
    ys = [p[1] for s in plan for p in s.points]
    span_drawn = max(ys) - min(ys)
    check("the plan view spans the real wingspan",
          abs(span_drawn - app.aircraft.wing.span_m) < 0.05,
          "%.3f vs %.3f m" % (span_drawn, app.aircraft.wing.span_m))
    xs = [p[0] for s in side for p in s.points]
    check("the side view is about a fuselage long",
          max(xs) >= app.aircraft.body.fuse_length_m - 0.3,
          "%.2f vs %.2f m" % (max(xs), app.aircraft.body.fuse_length_m))
    check("both views have shapes", len(plan) >= 4 and len(side) >= 5)

    # ------------------------------------------------------------------
    section("G. The scoring panel shows the formulas and the arithmetic")
    txt = text_of(app.txt_score)
    check("the scoring panel rendered", len(txt) > 800, "%d chars" % len(txt))
    for needed in ("M2 = 1 + N(weight/time) / Max(weight/time)",
                   "M3 = 2 + N(#laps x sensor weight)",
                   "Total Mission Score = M1 + M2 + M3 + GM"):
        check("it states the formula: %s" % needed[:34], needed in txt)
    check("it shows the payload weight M2 actually scores",
          "payload weight" in txt and "what M2 scores" in txt)
    check("it shows the time the ratio divides by", "time for 5 laps" in txt)
    check("it names the reference being used",
          "reference" in txt and ("field" in txt or "this run" in txt))
    check("it reports the total", "TOTAL =" in txt)
    # the arithmetic on screen must match the ScoreCard, not be re-derived
    from dbf2027.scoring import evaluate_scores as _ev
    card = _ev(app.aircraft, app.cfg, 1.5)
    check("the panel's total matches the scorecard",
          ("TOTAL = %.3f" % card.total_mission) in txt,
          "%.3f" % card.total_mission)
    check("a saturated ratio is flagged as such",
          ("SATURATED" in txt) == (card.m2_saturated or card.m3_saturated))

    # ------------------------------------------------------------------
    # The weather plan (rule 3.2.1e): a button on the Report tab, and the
    # plan it runs covers every preset for the design on screen.
    check("the Report tab has a Weather plan button",
          hasattr(app, "weather_btn") and callable(app.show_weather_plan))
    import io as _io
    from contextlib import redirect_stdout as _rs
    from dbf2027.report import print_weather_plan as _pwp
    with _rs(_io.StringIO()) as _buf:
        _plan = _pwp(app.design, app.cfg, 12.0)
    check("and the plan covers all three weather presets",
          set(_plan["presets"]) == {"tucson-morning", "tucson", "tucson-hot"}
          and "WEATHER PLAN" in _buf.getvalue(),
          str(sorted(_plan["presets"])))

    section("H. Course altitude and the areal densities are editable")
    check("pattern altitude is exposed", ("course", "pattern_altitude_m")
          in app.adv)
    for name in ("wing_areal_density", "fuse_areal_density",
                 "h_tail_areal_density", "v_tail_areal_density"):
        check("%s is exposed" % name, ("build", name) in app.adv)
    # This is a PLUMBING check: does the course field reach the model?
    # Whether a taller pattern costs time is physics, and selftest owns
    # it there with a design that has energy margin to spare - the
    # recommended design sits on the energy limit by construction, so a
    # taller climb can legitimately make it infeasible rather than slow.
    base_alt = app.adv[("course", "pattern_altitude_m")].get()
    app.evaluate()
    low = _ev(app.aircraft, app.read_config(), 1.5)
    t_low = low.r2.time_s if low.r2 else 0.0
    ok_low = bool(low.r2 and low.r2.ok)
    app.adv[("course", "pattern_altitude_m")].set(base_alt + 8.0)
    app.evaluate()
    high = _ev(app.aircraft, app.read_config(), 1.5)
    t_high = high.r2.time_s if high.r2 else 0.0
    ok_high = bool(high.r2 and high.r2.ok)
    check("the course altitude reaches the mission model",
          (t_high != t_low) or (ok_high != ok_low),
          "%.2f s (%s) vs %.2f s (%s)" % (t_low, ok_low, t_high, ok_high))
    check("a taller pattern costs time or costs feasibility",
          (not ok_high) or t_high >= t_low - 1e-9,
          "%.2f -> %.2f s" % (t_low, t_high))
    app.adv[("course", "pattern_altitude_m")].set(base_alt)
    app.evaluate()
    check("restoring the altitude restores the design",
          app.aircraft is not None)

    h_before = app.aircraft.mass_breakdown.items["horizontal tail (elevator)"]
    v_before = app.aircraft.mass_breakdown.items["vertical tail (rudder)"]
    app.adv[("build", "h_tail_areal_density")].set(3.60)
    app.evaluate()
    h_after = app.aircraft.mass_breakdown.items["horizontal tail (elevator)"]
    v_after = app.aircraft.mass_breakdown.items["vertical tail (rudder)"]
    # The elevator mass is its density times its area, both halves of
    # which the model now owns: a heavier tail also moves the lightest tail
    # arm (layout.tail_arm_optimum_m) aft, which shrinks both tails.
    _bd = app.aircraft.body
    check("the elevator density reaches the elevator mass",
          abs(h_after - 3.60 * _bd.s_horiz_m2) < 1e-9 and h_after > h_before,
          "%.4f vs %.4f kg" % (h_after, h_before))
    check("and the rudder mass is still its own density times its area",
          abs(v_after - app.aircraft.cfg.build.v_tail_areal_density
              * _bd.s_vert_m2) < 1e-9)
    app.adv[("build", "h_tail_areal_density")].set(1.80)
    app.evaluate()

    # ------------------------------------------------------------------
    section("I. Normal tab is a design; Search space tab is the bounds")
    # These are different things and the GUI used to conflate them: the
    # Normal tab looked like it steered the search, and did not - the
    # search built a default DesignSpace in code and ignored the form.
    from dbf2027.optimizer import DesignSpace as _DS
    from dbf2027.config import Config as _Cfg
    sp = app.read_space()
    check("the search space is readable from the form",
          isinstance(sp, _DS))
    b0 = sp.continuous_bounds(_Cfg())
    app.sp_span_lo.set(1.40)
    app.sp_span_hi.set(1.55)
    b1 = app.read_space().continuous_bounds(_Cfg())
    check("editing the span bounds reaches the optimiser",
          abs(b1[0][0] - 1.40) < 1e-9 and abs(b1[0][1] - 1.55) < 1e-9,
          str(b1[0]))
    check("a user lower bound is not silently dragged down",
          b1[0][0] > b0[0][0], "%.3f vs %.3f" % (b1[0][0], b0[0][0]))
    d_before = app.design
    app.write_design(replace(d_before, span_m=9.99), "test")
    b2 = app.read_space().continuous_bounds(_Cfg())
    check("the loaded design does NOT change the search bounds",
          b2[0] == b1[0], str(b2[0]))
    app.write_design(d_before, "restored")

    app.sp_motors.set("1")
    check("motor counts are editable",
          app.read_space().motor_counts == [1])
    app.sp_motors.set("1,2")
    app.sp_span_lo.set(2.5)
    app.sp_span_hi.set(1.5)
    inv = app.read_space().continuous_bounds(_Cfg())
    check("an inverted range is repaired, not left inverted",
          inv[0][0] <= inv[0][1], str(inv[0]))
    app.reset_search_space()
    check("reset restores the defaults",
          app.read_space().continuous_bounds(_Cfg())[0] == b0[0],
          str(app.read_space().continuous_bounds(_Cfg())[0]))
    check("the 6 ft rule still caps the upper bound whatever is typed",
          _DS(span_bounds=(0.8, 99.0)).continuous_bounds(_Cfg())[0][1]
          <= _Cfg().rules.max_wingspan_m + 1e-12)
    check("the search space is the only design input on the left",
          hasattr(app, "sp_span_lo") and not hasattr(app, "f_span"))

    # ------------------------------------------------------------------
    section("J. The GUI holds no state that outlived the form")
    # A stale panel that looks right is worse than an empty one, because
    # nothing tells you it is stale.  Bad input used to leave the
    # Dimensions, Score, Charts, Three-view and Parts list panels showing
    # a complete, plausible aeroplane from the previous evaluation - and
    # `self.aircraft` still pointed at it, so changing chart type redrew
    # the OLD design.
    app.evaluate()
    dims_live = [app.tv_dims.item(i)["values"] for i in
                 app.tv_dims.get_children()]
    check("a good design populates the dimensions table",
          len(dims_live) > 20 and app.aircraft is not None)
    app.f_wind.set("not a number")
    app.evaluate()
    check("bad input is reported", "bad input" in app.status.get())
    check("and the aircraft handle is dropped", app.aircraft is None)
    check("but the design itself is kept, since it is input",
          app.design is not None)
    check("and the dimensions table is cleared",
          len(app.tv_dims.get_children()) <= 1)
    check("and the score panel no longer shows a total",
          "TOTAL =" not in text_of(app.txt_score))
    check("and the parts list is cleared",
          "BILL OF MATERIALS" not in text_of(app.txt_bom))
    check("and the chart is cleared", app.chart_canvas.chart is None)
    check("and the three-view is cleared", app.three_view.data is None)
    check("redrawing a chart with no aircraft does not resurrect one",
          (app._draw_chart() or app.chart_canvas.chart) is None)
    app.f_wind.set(4.5)
    app.evaluate()
    check("a good design restores every panel",
          app.aircraft is not None
          and len(app.tv_dims.get_children()) > 20
          and "TOTAL =" in text_of(app.txt_score)
          and app.chart_canvas.chart is not None)

    # An unbuildable design must clear the panels too, not just print a
    # traceback into the report.
    keep = app.design
    app.write_design(replace(keep, motor="NO_SUCH_MOTOR"), "test")
    app.evaluate()
    check("a design that cannot be built clears the panels",
          app.aircraft is None and "could not be built"
          in text_of(app.txt_report))
    check("and the charts go with it", app.chart_canvas.chart is None)
    app.write_design(keep, "restored")
    app.evaluate()
    check("recovering from that works", app.aircraft is not None)

    # The opening design must not be a frozen claim about what is good.
    import json as _json
    with open("examples/optimised_v1.json", encoding="utf-8") as fh:
        rec = _json.load(fh)
    d0, src0 = App._opening_design_and_source()
    check("the window opens on the recommended design, not a hard-coded one",
          d0.airfoil == rec["airfoil"] and d0.motor == rec["motor"]
          and d0.prop == rec["prop"],
          "%s/%s vs %s/%s" % (d0.airfoil, d0.motor,
                              rec["airfoil"], rec["motor"]))
    check("and it says the file is where it came from",
          "optimised_v1.json" in src0, src0)
    check("its chord is derived from the file's aspect ratio",
          abs(d0.span_m / d0.chord_m
              - rec["span_m"] / rec["chord_m"]) < 1e-9)
    # Checking a handful of fields let `taper` go missing: the window
    # opened on a rectangular wing and scored a different aeroplane.
    # Every field the file carries must survive the load.
    from dataclasses import fields as _fields
    from dbf2027 import Design as _Design
    _lost = []
    for _f in _fields(_Design):
        if _f.name not in rec:
            continue
        _a, _b = getattr(d0, _f.name), rec[_f.name]
        _same = (_a == _b if isinstance(_b, str)
                 else abs(float(_a) - float(_b)) < 1e-9)
        if not _same:
            _lost.append(_f.name)
    check("every design field survives the load, taper included",
          not _lost, "lost or changed: %s" % _lost)

    # The search's normalisers must be visible and applied, not held in a
    # private attribute nobody reads.
    check("no dead reference attributes remain",
          not hasattr(app, "ref_m2") and not hasattr(app, "ref_m3"))
    for name in ("field_best_m2_lb_per_s", "field_best_m3_lap_lb",
                 "field_best_gm_lb_in"):
        check("the %s normaliser is editable" % name.replace("field_best_", ""),
              ("rules", name) in app.adv)

    # Timers must not outlive the window.
    check("pending timers are tracked so they can be cancelled",
          isinstance(app._after_ids, list) and app._after_ids)

    # ------------------------------------------------------------------
    section("K. Catalogues can be browsed and extended")
    from dbf2027 import components as _C
    # Every catalogue write in this section goes to a scratch folder, never
    # the real catalogues/ - a test must not rewrite the motors someone has
    # entered.  The real file is snapshotted and compared at the end.
    import os as _os, shutil as _sh, tempfile as _tf
    _real_dir = _C.user_dir()
    _real_csv = _os.path.join(_real_dir, "motors.csv")
    _snapshot = (open(_real_csv, "rb").read()
                 if _os.path.exists(_real_csv) else None)
    _real_user_dir = _C.user_dir
    _sandbox = _tf.mkdtemp(prefix="dbf_guitest_")
    _C.user_dir = lambda: _sandbox
    _C.reload_catalogues()
    check("the catalogue tab exists", hasattr(app, "cat_tree"))
    for label, key, n in (("Airfoils", "airfoils", 20),
                          ("Motors", "motors", 20),
                          ("Propellers", "props", 28),
                          ("Battery packs", "packs", 17),
                          ("ESCs", "escs", 5)):
        app._select_catalogue(label)
        rows = len(app.cat_tree.get_children())
        cols = list(app.cat_tree.cget("columns"))
        check("%s list every entry" % label, rows >= n, "%d rows" % rows)
        check("%s columns come from the spec" % label,
              cols == _C.CATALOGUES_BY_KEY[key].fieldnames(), str(cols[:3]))
        check("%s form has a box per column" % label,
              set(app.cat_fields) == set(cols))
        check("%s dropdown agrees with the table" % label,
              app.cat_combo.get() == label, app.cat_combo.get())

    # Adding one must persist, apply immediately, and be usable.
    app._select_catalogue("Motors")
    before = len(_C.MOTORS)
    for name, val in (("key", "GUITEST_M"), ("name", "GUI test motor"),
                      ("kv", 620), ("rm", 0.045), ("i0", 1.2),
                      ("mass_kg", 0.230), ("i_max_a", 48),
                      ("p_max_w", 1050), ("shaft_mm", 5)):
        app.cat_fields[name].set(val)
    app._add_catalogue_entry()
    check("adding an entry puts it in the live catalogue",
          "GUITEST_M" in _C.MOTORS and len(_C.MOTORS) == before + 1)
    check("and the values round-trip through the constructor",
          abs(_C.MOTORS["GUITEST_M"].kv - 620) < 1e-9
          and abs(_C.MOTORS["GUITEST_M"].rm - 0.045) < 1e-12)
    check("it is marked as yours, not built in",
          _C.is_user_entry("motors", "GUITEST_M")
          and not _C.is_user_entry("motors", "SCORPION_SII_4025_440"))
    check("it is written to the user CSV",
          any(r.get("key") == "GUITEST_M" for r in _C.user_rows("motors")))
    check("it appears in the table", "GUITEST_M" in
          app.cat_tree.get_children())
    check("and is tinted as a user row",
          "user" in app.cat_tree.item("GUITEST_M", "tags"))
    check("the source field says it needs verifying",
          "VERIFY" in _C.MOTORS["GUITEST_M"].source)

    # It must be a real, usable component, not just a table row.
    d_keep = app.design
    app.write_design(replace(d_keep, motor="GUITEST_M"), "catalogue test")
    app.evaluate()
    check("a design can be built on the added motor",
          app.aircraft is not None
          and app.aircraft.motor.key == "GUITEST_M")
    check("and the search space offers it",
          "GUITEST_M" in _DS().motors)
    app.write_design(d_keep, "restored")
    app.evaluate()

    # Reload from disk must keep it; delete must remove it.
    _C.reload_catalogues()
    check("it survives a reload from disk", "GUITEST_M" in _C.MOTORS)
    app._select_catalogue("Motors")
    app.cat_tree.selection_set("GUITEST_M")
    _C.remove_user_entry("motors", "GUITEST_M")
    app._refresh_catalogue()
    check("deleting removes it from the catalogue and the CSV",
          "GUITEST_M" not in _C.MOTORS
          and not any(r.get("key") == "GUITEST_M"
                      for r in _C.user_rows("motors")))
    check("built-in entries cannot be deleted",
          _C.remove_user_entry("motors", "SCORPION_SII_4025_440") is False
          and "SCORPION_SII_4025_440" in _C.MOTORS)
    _C.user_dir = _real_user_dir
    _C.reload_catalogues()
    app._refresh_catalogue()
    _sh.rmtree(_sandbox, ignore_errors=True)
    check("and your real catalogues/ folder was never touched",
          (open(_real_csv, "rb").read() if _os.path.exists(_real_csv)
           else None) == _snapshot)

    # ------------------------------------------------------------------
    section("L. The rules panel is traceable, not decorative")
    rtxt = text_of(app.txt_rules)
    check("the rules panel rendered", len(rtxt) > 3000, "%d chars" % len(rtxt))
    from dbf2027.rulebook import RULES, check_annotations
    check("annotations and rulebook agree",
          check_annotations() == ([], []), str(check_annotations()))
    check("it says so on screen", "Annotations and rulebook agree" in rtxt)
    for rid in ("3.1.3a", "3.3.4", "3.2.3.2d", "6.3"):
        check("rule %s appears" % rid, rid in rtxt)
    check("it cites source locations", ".py:" in rtxt and "enforced at" in rtxt)
    check("build-only rules are marked as such",
          "build / tech-inspection checklist" in rtxt)
    check("every registered rule is listed",
          all(r.id in rtxt for r in RULES),
          str([r.id for r in RULES if r.id not in rtxt])[:60])

    # ------------------------------------------------------------------
    section("M. Function catalogue")
    check("the catalogue found the package",
          len(app._catalogue) > 100, "%d entries" % len(app._catalogue))
    check("the catalogue panel rendered",
          "FUNCTION CATALOGUE" in text_of(app.txt_cat))
    app.cat_query.set("trim")
    app._draw_catalogue()
    filtered = text_of(app.txt_cat)
    check("filtering narrows the catalogue",
          len(filtered) < len(catalogue_full(app)),
          "%d vs %d chars" % (len(filtered), len(catalogue_full(app))))
    check("the filter finds the trim solver", "tail_load" in filtered)
    app.cat_query.set("")
    app._draw_catalogue()
    kinds = {e.kind for e in app._catalogue}
    check("classes, functions, methods and properties are all indexed",
          kinds >= {"class", "function", "method", "property"}, str(kinds))
    check("no private names leak in",
          not any(e.qualname.startswith("_") for e in app._catalogue))
    check("every entry names a real module",
          all(e.module for e in app._catalogue))

    # ------------------------------------------------------------------
    section("N. Save / load round trip")
    import tempfile
    path = tempfile.mktemp(suffix=".json")
    d_before = app.read_design()
    with open(path, "w", encoding="utf-8") as fh:
        from dataclasses import fields as _fields
        json.dump({f.name: getattr(d_before, f.name)
                   for f in _fields(Design)}, fh)
    with open(path, encoding="utf-8") as fh:
        d_after = Design(**json.load(fh))
    app.write_design(d_after, "round trip")
    check("a design survives the JSON round trip exactly",
          app.read_design().key() == d_before.key())
    check("with no precision lost to a text box",
          app.read_design().span_m == d_before.span_m,
          "%.17g vs %.17g" % (app.read_design().span_m, d_before.span_m))

    # ------------------------------------------------------------------
    section("O. Background work shows progress")
    import time as _time
    check("the status bar has a progress bar",
          hasattr(app, "progress")
          and str(app.progress.cget("mode")) == "determinate")
    app._progress_begin("starting the quick search")
    app.search_btn.configure(state="disabled")
    app.search_q.put(("progress", (0.10, "Stage 2 of 5, coarse sweep: "
                                         "10 of 100 designs")))
    app.search_q.put(("progress", (0.40, "Stage 3 of 5, local refinement: "
                                         "3 of 12 configurations")))
    for _ in range(40):                 # let the real poll loop drain it
        root.update()
        _time.sleep(0.02)
    _pt = app.progress_text.get()
    check("a search's progress reaches the bar",
          abs(float(app.progress.cget("value")) - 400.0) < 1e-6,
          str(app.progress.cget("value")))
    check("with the stage, the percentage and the time left",
          "Stage 3 of 5" in _pt and "40%" in _pt and "left" in _pt, _pt)
    check("a progress message does not re-enable Run search",
          str(app.search_btn.cget("state")) == "disabled")
    app.search_btn.configure(state="normal")
    app._progress_busy("weather plan: flying three presets")
    check("the weather plan shows a busy bar",
          str(app.progress.cget("mode")) == "indeterminate")
    app._progress_end("weather plan done in 0:05")
    check("which stops, full, when the work is done",
          str(app.progress.cget("mode")) == "determinate"
          and float(app.progress.cget("value")) == 1000.0
          and "done" in app.progress_text.get())
    app._progress_end()
    check("and clears back to idle", app.progress_text.get() == ""
          and float(app.progress.cget("value")) == 0.0)

    root.update()
    root.destroy()

    print("\n" + "=" * 60)
    print("%d passed, %d failed" % (PASS, FAIL))
    print("=" * 60)
    return 1 if FAIL else 0


def catalogue_full(app) -> str:
    from dbf2027.catalogue import as_text
    return as_text(app._catalogue)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc()
        sys.exit(2)
