"""Chart data and the arithmetic that turns it into pixels.

Deliberately free of tkinter.  This module computes *what* to draw - the
series, the axis ranges, the tick positions, the mapping from data units to
canvas coordinates - and `gui.py` does nothing but stroke lines along the
numbers it produces.  Two reasons: the geometry is then testable without a
display, and the same series can be written to CSV or drawn by something
else later.

There is no matplotlib here because there is no matplotlib anywhere in this
project: standard library only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Sequence, Tuple

from .units import KG2LB, MS2MPH, N2LBF


# ==========================================================================
# Series and axes
# ==========================================================================
@dataclass
class Series:
    """One named curve in data units, plus how to draw it."""
    label: str
    xs: List[float]
    ys: List[float]
    colour: str = "#2f6fdb"
    dash: Tuple[int, ...] = ()
    width: int = 2
    # A marker drawn at one point, e.g. the design cruise condition.
    marker: Optional[Tuple[float, float]] = None


@dataclass
class Chart:
    """A complete figure: title, axis labels, series and axis ranges."""
    title: str
    x_label: str
    y_label: str
    series: List[Series] = field(default_factory=list)
    x_range: Optional[Tuple[float, float]] = None
    y_range: Optional[Tuple[float, float]] = None
    # Horizontal rules, e.g. "this is the weight you have to lift".
    h_lines: List[Tuple[float, str, str]] = field(default_factory=list)
    footnote: str = ""

    def bounds(self) -> Tuple[float, float, float, float]:
        """(xlo, xhi, ylo, yhi), from the data unless overridden."""
        xs = [x for s in self.series for x in s.xs]
        ys = [y for s in self.series for y in s.ys]
        ys += [v for v, _, _ in self.h_lines]
        if not xs or not ys:
            return (0.0, 1.0, 0.0, 1.0)
        xlo, xhi = (self.x_range if self.x_range else (min(xs), max(xs)))
        ylo, yhi = (self.y_range if self.y_range else (min(ys), max(ys)))
        if xhi - xlo < 1e-12:
            xhi = xlo + 1.0
        if yhi - ylo < 1e-12:
            yhi = ylo + 1.0
        if self.y_range is None:            # a little headroom, and include 0
            ylo = min(ylo, 0.0)
            yhi += 0.08 * (yhi - ylo)
        return (xlo, xhi, ylo, yhi)


def nice_ticks(lo: float, hi: float, target: int = 6) -> List[float]:
    """Tick positions at 1, 2, 2.5 or 5 times a power of ten.

    The usual trick: take the raw spacing the target count implies, then
    round it up to the next "human" number so the labels read 0, 5, 10
    rather than 0, 4.7, 9.4.
    """
    span = hi - lo
    if span <= 0 or not math.isfinite(span):
        return [lo]
    raw = span / max(target, 1)
    mag = 10.0 ** math.floor(math.log10(raw))
    for mult in (1.0, 2.0, 2.5, 5.0, 10.0):
        if raw <= mult * mag:
            step = mult * mag
            break
    else:
        step = 10.0 * mag
    first = math.ceil(lo / step) * step
    out, t = [], first
    while t <= hi + 1e-9 * step:
        out.append(0.0 if abs(t) < 1e-12 else t)
        t += step
    return out


def tick_label(v: float, step: float) -> str:
    """Format a tick so it has exactly as many decimals as it needs."""
    if step >= 1.0:
        return "%.0f" % v
    if step >= 0.1:
        return "%.1f" % v
    if step >= 0.01:
        return "%.2f" % v
    return "%.3f" % v


class Mapper:
    """Data units to pixels, for a plot area inset in a canvas.

    y is flipped, because canvases count downward and graphs count upward.
    """

    def __init__(self, chart: Chart, width: int, height: int,
                 pad_left: int = 72, pad_right: int = 18,
                 pad_top: int = 38, pad_bottom: int = 52):
        self.chart = chart
        self.x0, self.x1, self.y0, self.y1 = chart.bounds()
        self.left, self.right = pad_left, width - pad_right
        self.top, self.bottom = pad_top, height - pad_bottom
        if self.right <= self.left:
            self.right = self.left + 1
        if self.bottom <= self.top:
            self.bottom = self.top + 1

    @property
    def plot_width(self) -> int:
        return self.right - self.left

    @property
    def plot_height(self) -> int:
        return self.bottom - self.top

    def px(self, x: float) -> float:
        f = (x - self.x0) / (self.x1 - self.x0)
        return self.left + f * self.plot_width

    def py(self, y: float) -> float:
        f = (y - self.y0) / (self.y1 - self.y0)
        return self.bottom - f * self.plot_height

    def points(self, xs: Sequence[float],
               ys: Sequence[float]) -> List[float]:
        """Flat [x0, y0, x1, y1, ...] ready for a canvas line."""
        out: List[float] = []
        for x, y in zip(xs, ys):
            if not (math.isfinite(x) and math.isfinite(y)):
                continue
            out.extend((self.px(x), self.py(y)))
        return out


# ==========================================================================
# The charts themselves
# ==========================================================================
PALETTE = ["#2f6fdb", "#d9534f", "#2e9e5b", "#b8860b", "#7b52c7",
           "#0d9aa8", "#c2571a"]


def _speed_axis(ac, weight_n: float, n: int = 70) -> List[float]:
    """Speeds from just below stall to comfortably past the top speed.

    The upper end is taken from the aeroplane's own maximum level speed so
    the interesting part of the curve fills the axis instead of being
    squashed into the left third.
    """
    try:
        v_lo = max(4.0, 0.80 * ac.aero_clean.v_stall(weight_n))
    except Exception:
        v_lo = 8.0
    try:
        v_hi = 1.30 * ac._max_level_speed(weight_n)
    except Exception:
        v_hi = 45.0
    v_hi = max(v_hi, v_lo + 10.0)
    return [v_lo + (v_hi - v_lo) * i / (n - 1) for i in range(n)]


def forces_vs_speed(ac, weights_kg: Optional[Sequence[float]] = None,
                    imperial: bool = True) -> Chart:
    """Drag against airspeed at several weights, with thrust over the top.

    This is the chart that answers "how fast can it go and how hard is it
    working to get there".  Where the thrust curve crosses a drag curve is
    the maximum level speed at that weight; the minimum of each drag curve
    is its best-endurance speed, and the point where a line from the origin
    is tangent to it is best range.

    Drag is the *trimmed* value, so it includes the induced drag of the tail
    load needed to balance the aeroplane at that speed - which is why the
    curves do not collapse onto one another when you scale weight.
    """
    if weights_kg is None:
        weights_kg = [ac.mtow_m1_kg, ac.mtow_m3_kg, ac.mtow_m2_kg]
    fy = N2LBF if imperial else 1.0
    fx = MS2MPH if imperial else 1.0
    chart = Chart(
        title="Drag and thrust against airspeed",
        x_label="airspeed  [%s]" % ("mph" if imperial else "m/s"),
        y_label="force  [%s]" % ("lbf" if imperial else "N"))

    vs = _speed_axis(ac, ac.mtow_m2_kg * 9.80665)
    for i, m in enumerate(sorted(set(round(w, 4) for w in weights_kg))):
        w_n = m * 9.80665
        xs, ys = [], []
        for v in vs:
            try:
                d = ac.aero_clean.drag(v, w_n)
            except Exception:
                continue
            if math.isfinite(d) and d > 0:
                xs.append(v * fx)
                ys.append(d * fy)
        if xs:
            chart.series.append(Series(
                "drag at %.1f lb" % (m * KG2LB) if imperial
                else "drag at %.1f kg" % m,
                xs, ys, PALETTE[i % len(PALETTE)]))

    xs, ys = [], []
    for v in vs:
        try:
            t = ac.prop_sys.thrust(v, 0.75, 1.0)
        except Exception:
            continue
        xs.append(v * fx)
        ys.append(t * fy)
    if xs:
        chart.series.append(Series("thrust available (full throttle)",
                                   xs, ys, "#333333", dash=(6, 4)))
    chart.footnote = ("Thrust is at 75% state of charge. A drag curve "
                      "crossing it is that weight's top speed.")
    return chart


def lift_vs_drag(ac, v_mps: Optional[float] = None,
                 imperial: bool = True) -> Chart:
    """Lift against drag at a fixed airspeed - the drag polar in force units.

    Sweeping *speed* at a fixed weight would be pointless here, because in
    level flight lift always equals weight and the curve would be a flat
    line.  So this sweeps lift coefficient instead, at one airspeed, which
    is what a polar actually is: how much drag you pay for each newton of
    lift you ask the wing for.

    Horizontal rules mark the two mission weights.  Where a rule crosses
    the curve is the drag you will be pulling all lap at that weight, and
    if a rule sits above the top of the curve the wing cannot hold the
    aeroplane up at this speed at all.
    """
    a = ac.aero_clean
    if v_mps is None:
        v_mps = 1.35 * a.v_stall(ac.mtow_m2_kg * 9.80665)
    fy = N2LBF if imperial else 1.0
    q_s = a.q(v_mps) * ac.wing.area_m2
    cl_max = a.cl_max_at_v(v_mps)
    chart = Chart(
        title="Lift against drag at %.0f %s"
              % (v_mps * (MS2MPH if imperial else 1.0),
                 "mph" if imperial else "m/s"),
        x_label="drag  [%s]" % ("lbf" if imperial else "N"),
        y_label="lift  [%s]" % ("lbf" if imperial else "N"))

    xs, ys = [], []
    for i in range(80):
        cl = cl_max * (i + 1) / 80.0
        lift = q_s * cl
        try:
            drag = a.drag(v_mps, lift)
        except Exception:
            continue
        if math.isfinite(drag) and drag > 0:
            xs.append(drag * fy)
            ys.append(lift * fy)
    chart.series.append(Series("trimmed polar", xs, ys, PALETTE[0]))

    for i, (label, w_kg) in enumerate((("Mission 3 weight", ac.mtow_m3_kg),
                                       ("Mission 2 weight", ac.mtow_m2_kg))):
        chart.h_lines.append((w_kg * 9.80665 * fy,
                              "%s  %.1f lb" % (label, w_kg * KG2LB),
                              PALETTE[(i + 1) % len(PALETTE)]))
    chart.footnote = ("Each point is one angle of attack. A weight line "
                      "above the curve means the wing cannot lift it at "
                      "this speed.")
    return chart


def lift_to_drag(ac, imperial: bool = True) -> Chart:
    """L/D against airspeed at each mission weight.

    The peak is the speed to fly for range; it moves right as the aeroplane
    gets heavier, which is exactly why Mission 2 and Mission 3 want
    different cruise speeds.
    """
    fx = MS2MPH if imperial else 1.0
    chart = Chart(title="Lift-to-drag ratio against airspeed",
                  x_label="airspeed  [%s]" % ("mph" if imperial else "m/s"),
                  y_label="L / D")
    vs = _speed_axis(ac, ac.mtow_m2_kg * 9.80665, n=90)
    for i, (label, w_kg) in enumerate((("empty", ac.mtow_m1_kg),
                                       ("Mission 3", ac.mtow_m3_kg),
                                       ("Mission 2", ac.mtow_m2_kg))):
        w_n = w_kg * 9.80665
        xs, ys = [], []
        for v in vs:
            try:
                d = ac.aero_clean.drag(v, w_n)
            except Exception:
                continue
            if math.isfinite(d) and d > 1e-6:
                xs.append(v * fx)
                ys.append(w_n / d)
        if xs:
            chart.series.append(Series(
                "%s  (%.1f lb)" % (label, w_kg * KG2LB), xs, ys,
                PALETTE[i % len(PALETTE)]))
    chart.footnote = "Peak L/D is the speed to fly for range."
    return chart


def power_vs_speed(ac, imperial: bool = True) -> Chart:
    """Electrical power required against airspeed, at Mission 2 weight.

    Power, not drag, is what empties the battery, and the 100 Wh cap is
    usually what stops a design - so this is the curve to look at when a
    design runs out of energy before it runs out of time.
    """
    fx = MS2MPH if imperial else 1.0
    chart = Chart(title="Power required against airspeed",
                  x_label="airspeed  [%s]" % ("mph" if imperial else "m/s"),
                  y_label="power  [W]")
    vs = _speed_axis(ac, ac.mtow_m2_kg * 9.80665, n=70)
    for i, (label, w_kg) in enumerate((("Mission 3", ac.mtow_m3_kg),
                                       ("Mission 2", ac.mtow_m2_kg))):
        w_n = w_kg * 9.80665
        xs, ys = [], []
        for v in vs:
            try:
                d = ac.aero_clean.drag(v, w_n)
                thr = ac.prop_sys.throttle_for_thrust(v, 0.75, d)
                p = ac.prop_sys.power_elec(v, 0.75, thr)
            except Exception:
                continue
            if math.isfinite(p) and p > 0:
                xs.append(v * fx)
                ys.append(p)
        if xs:
            chart.series.append(Series(
                "%s  (%.1f lb)" % (label, w_kg * KG2LB), xs, ys,
                PALETTE[i % len(PALETTE)]))
    usable_w = ac.usable_energy_j() / max(ac.cfg.rules.mission_window_s, 1.0)
    chart.h_lines.append((usable_w,
                          "average power that empties the pack in %.0f s"
                          % ac.cfg.rules.mission_window_s, "#d9534f"))
    chart.footnote = ("Below the red line the aeroplane finishes the window "
                      "on one charge; above it, energy runs out first.")
    return chart


def weight_sweep(ac, cfg, lo_kg: Optional[float] = None,
                 hi_kg: Optional[float] = None, n: int = 26,
                 imperial: bool = True) -> Chart:
    """Stall, best-L/D and maximum level speed as gross weight varies.

    This is the "how much can I actually carry" chart: the stall and top
    speed curves converge as weight rises, and where they meet is the point
    the aeroplane can no longer fly at all.
    """
    fx = MS2MPH if imperial else 1.0
    lo = lo_kg if lo_kg else ac.empty_mass_kg * 1.02
    hi = hi_kg if hi_kg else min(cfg.rules.max_gross_mass_kg,
                                 max(ac.mtow_m2_kg * 1.6, lo * 1.6))
    chart = Chart(title="Speeds against gross weight",
                  x_label="gross weight  [%s]" % ("lb" if imperial else "kg"),
                  y_label="airspeed  [%s]" % ("mph" if imperial else "m/s"))
    fw = KG2LB if imperial else 1.0
    xs_s, ys_s, xs_t, ys_t, xs_l, ys_l = [], [], [], [], [], []
    for i in range(n):
        m = lo + (hi - lo) * i / (n - 1)
        w_n = m * 9.80665
        try:
            v_s = ac.aero_clean.v_stall(w_n)
        except Exception:
            continue
        xs_s.append(m * fw)
        ys_s.append(v_s * fx)
        # best L/D speed, by scan
        best_v, best_ld = 0.0, 0.0
        for k in range(40):
            v = v_s * (1.0 + 2.2 * k / 39.0)
            try:
                d = ac.aero_clean.drag(v, w_n)
            except Exception:
                continue
            if d > 1e-6 and w_n / d > best_ld:
                best_ld, best_v = w_n / d, v
        if best_v > 0:
            xs_l.append(m * fw)
            ys_l.append(best_v * fx)
        try:
            v_max = ac._max_level_speed(w_n)
            xs_t.append(m * fw)
            ys_t.append(v_max * fx)
        except Exception:
            pass
    chart.series.append(Series("stall speed", xs_s, ys_s, PALETTE[1]))
    chart.series.append(Series("best L/D speed", xs_l, ys_l, PALETTE[2],
                               dash=(5, 3)))
    chart.series.append(Series("maximum level speed", xs_t, ys_t, PALETTE[0]))
    chart.footnote = ("Where stall meets top speed the aeroplane cannot fly. "
                      "Gross weight limit is %.0f lb."
                      % (cfg.rules.max_gross_mass_kg * KG2LB))
    return chart


CHART_BUILDERS = [
    ("Drag and thrust vs speed", forces_vs_speed),
    ("Lift vs drag", lift_vs_drag),
    ("L/D vs speed", lift_to_drag),
    ("Power required vs speed", power_vs_speed),
]


# ==========================================================================
# Three-view geometry
# ==========================================================================
@dataclass
class Shape:
    """A closed polygon in metres, in the plane of one view."""
    points: List[Tuple[float, float]]
    fill: str = ""
    outline: str = "#333333"
    width: int = 1


def three_view(ac) -> Tuple[List[Shape], List[Shape], Tuple[float, float]]:
    """Plan and side outlines of the aeroplane, in metres.

    Returns (plan shapes, side shapes, (extent_x, extent_y)) with the origin
    at the nose on the centreline.  Crude but to scale: enough to see that
    the tail arm is sane, that the props clear each other and the ground,
    and that the fuselage is the length the bay forced it to be.
    """
    w, b = ac.wing, ac.body
    span, c_root = w.span_m, w.root_chord_m
    c_tip = c_root * w.taper
    L = b.fuse_length_m
    # Wing quarter-chord sits a little ahead of mid-fuselage.
    x_le = 0.30 * L
    plan: List[Shape] = []

    # fuselage in plan
    hw = 0.5 * b.fuse_width_m
    plan.append(Shape([(0.0, 0.0), (0.16 * L, -hw), (0.86 * L, -hw),
                       (L, -0.30 * hw), (L, 0.30 * hw), (0.86 * L, hw),
                       (0.16 * L, hw), (0.0, 0.0)], "#eef2f8", "#334", 1))

    # wing: a straight-tapered half panel, mirrored
    sweep = 0.25 * (c_root - c_tip)       # quarter-chord kept straight
    for sgn in (-1.0, 1.0):
        plan.append(Shape([
            (x_le, 0.0),
            (x_le + sweep, sgn * 0.5 * span),
            (x_le + sweep + c_tip, sgn * 0.5 * span),
            (x_le + c_root, 0.0)], "#cfe0f7", "#2f6fdb", 1))

    # horizontal tail at the tail arm
    x_h = x_le + 0.25 * c_root + b.tail_arm_m
    ar_t = 4.0
    b_h = math.sqrt(max(b.s_horiz_m2, 1e-6) * ar_t)
    c_h = max(b.s_horiz_m2, 1e-6) / max(b_h, 1e-6)
    plan.append(Shape([(x_h - 0.5 * c_h, -0.5 * b_h),
                       (x_h + 0.5 * c_h, -0.5 * b_h),
                       (x_h + 0.5 * c_h, 0.5 * b_h),
                       (x_h - 0.5 * c_h, 0.5 * b_h)], "#dfe8d8", "#2e9e5b", 1))

    # propeller discs
    n_m = ac.design.n_motors
    r_p = 0.5 * ac.prop.diameter_m
    if n_m == 1:
        centres = [(0.0, 0.0)]                       # nose-mounted tractor
    else:
        # Wing-mounted tractors sit ahead of the leading edge, far enough
        # that the disc clears it; spaced so the discs cannot overlap.
        y_m = max(0.25 * span, 0.55 * ac.prop.diameter_m)
        centres = [(x_le - 0.55 * r_p, -y_m), (x_le - 0.55 * r_p, y_m)]
    for cx, cy in centres:
        plan.append(Shape(_circle(cx, cy, r_p), "", "#b03a2e", 1))

    # ---- side view -------------------------------------------------
    side: List[Shape] = []
    hh = 0.5 * b.fuse_height_m
    side.append(Shape([(0.0, 0.0), (0.12 * L, -hh), (0.80 * L, -hh),
                       (L, -0.25 * hh), (L, 0.35 * hh), (0.80 * L, hh),
                       (0.12 * L, hh), (0.0, 0.0)], "#eef2f8", "#334", 1))
    t = w.airfoil.t_over_c * c_root
    side.append(Shape([(x_le, -0.15 * t), (x_le + c_root, 0.0),
                       (x_le + 0.35 * c_root, -0.9 * t)],
                      "#cfe0f7", "#2f6fdb", 1))
    # vertical tail
    ar_v = 1.6
    h_v = math.sqrt(max(b.s_vert_m2, 1e-6) * ar_v)
    c_v = max(b.s_vert_m2, 1e-6) / max(h_v, 1e-6)
    side.append(Shape([(x_h - 0.5 * c_v, -hh), (x_h + 0.6 * c_v, -hh),
                       (x_h + 0.35 * c_v, -hh - h_v),
                       (x_h - 0.2 * c_v, -hh - h_v)],
                      "#e8dfef", "#7b52c7", 1))
    # horizontal tail edge-on
    side.append(Shape([(x_h - 0.5 * c_h, -0.1 * hh), (x_h + 0.5 * c_h, 0.0),
                       (x_h - 0.5 * c_h, 0.1 * hh)], "#dfe8d8", "#2e9e5b", 1))
    # gear and ground line
    g = ac.gear_height_m
    side.append(Shape([(x_le + 0.3 * c_root, hh),
                       (x_le + 0.3 * c_root, hh + g)], "", "#555", 2))
    side.append(Shape(_circle(x_le + 0.3 * c_root, hh + g, 0.05), "", "#555", 2))
    side.append(Shape([(0.10 * L, hh), (0.10 * L, hh + 0.75 * g)], "", "#555", 2))
    side.append(Shape(_circle(0.10 * L, hh + 0.75 * g, 0.035), "", "#555", 2))
    if n_m == 1:
        side.append(Shape([(0.0, -r_p), (0.0, r_p)], "", "#b03a2e", 2))

    return plan, side, (max(L, span), max(span, b.fuse_height_m + g + 0.2))


def _circle(cx: float, cy: float, r: float, n: int = 36):
    return [(cx + r * math.cos(2 * math.pi * i / n),
             cy + r * math.sin(2 * math.pi * i / n)) for i in range(n + 1)]


def dimensions(ac) -> List[Tuple[str, str, str]]:
    """(label, SI, imperial) rows for the dimensions panel."""
    w, b, d = ac.wing, ac.body, ac.design
    IN = 39.37007874
    FT2 = 10.76391
    rows = [
        ("Wing span", "%.3f m" % w.span_m, "%.1f in  (%.2f ft)"
         % (w.span_m * IN, w.span_m * IN / 12.0)),
        ("Wing area", "%.4f m2" % w.area_m2, "%.0f in2  (%.2f ft2)"
         % (w.area_m2 * IN * IN, w.area_m2 * FT2)),
        ("Mean chord", "%.3f m" % w.chord_m, "%.2f in" % (w.chord_m * IN)),
        ("Root chord", "%.3f m" % w.root_chord_m,
         "%.2f in" % (w.root_chord_m * IN)),
        ("Tip chord", "%.3f m" % (w.root_chord_m * w.taper),
         "%.2f in" % (w.root_chord_m * w.taper * IN)),
        ("Taper ratio", "%.3f" % w.taper, ""),
        ("Aspect ratio", "%.2f" % w.aspect_ratio, ""),
        ("Root thickness", "%.1f mm" % (1000 * w.thickness_root_m),
         "%.2f in" % (w.thickness_root_m * IN)),
        ("", "", ""),
        ("Fuselage length", "%.3f m" % b.fuse_length_m,
         "%.1f in" % (b.fuse_length_m * IN)),
        ("Fuselage width", "%.3f m" % b.fuse_width_m,
         "%.1f in" % (b.fuse_width_m * IN)),
        ("Fuselage height", "%.3f m" % b.fuse_height_m,
         "%.1f in" % (b.fuse_height_m * IN)),
        ("Wetted area", "%.3f m2" % b.fuse_wetted_m2,
         "%.2f ft2" % (b.fuse_wetted_m2 * FT2)),
        ("Fineness ratio", "%.2f" % b.fineness, ""),
        ("", "", ""),
        ("Tail arm", "%.3f m" % b.tail_arm_m,
         "%.1f in" % (b.tail_arm_m * IN)),
        ("Horizontal tail area", "%.4f m2" % b.s_horiz_m2,
         "%.0f in2" % (b.s_horiz_m2 * IN * IN)),
        ("Vertical tail area", "%.4f m2" % b.s_vert_m2,
         "%.0f in2" % (b.s_vert_m2 * IN * IN)),
        ("", "", ""),
        ("Propeller diameter", "%.3f m" % ac.prop.diameter_m,
         "%.1f x %.1f in" % (ac.prop.diameter_m * IN, ac.prop.pitch_m * IN)),
        ("Landing gear height", "%.0f mm" % (1000 * ac.gear_height_m),
         "%.2f in" % (ac.gear_height_m * IN)),
        ("Payload bay", "%d x %d x %d boxes"
         % (ac.bay_rows, ac.bay_cols, ac.bay_layers), ""),
        ("", "", ""),
        ("Empty mass", "%.3f kg" % ac.empty_mass_kg,
         "%.2f lb" % (ac.empty_mass_kg * KG2LB)),
        ("Mission 2 gross", "%.3f kg" % ac.mtow_m2_kg,
         "%.2f lb" % (ac.mtow_m2_kg * KG2LB)),
        ("Mission 3 gross", "%.3f kg" % ac.mtow_m3_kg,
         "%.2f lb" % (ac.mtow_m3_kg * KG2LB)),
        ("Wing loading (M2)", "%.1f N/m2" % ac.wing_loading_m2(),
         "%.2f oz/ft2" % (ac.wing_loading_m2() / 9.80665 * 35.274 / FT2)),
    ]
    return rows


# ==========================================================================
# SVG rendering
# ==========================================================================
# The desktop GUI strokes these same charts onto a tkinter Canvas.  The web
# front end cannot, so the identical geometry is emitted as SVG instead.
# Both go through `Mapper`, so a chart looks the same either way and there
# is only one place where a chart can be wrong.
def _esc(text: str) -> str:
    return (str(text).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def chart_svg(chart: Chart, width: int = 960, height: int = 560) -> str:
    """Render a `Chart` as a standalone SVG string.

    A deeper bottom margin than the canvas uses: an SVG carries both the
    axis label and the footnote below the frame, and at the default height
    the two collided.
    """
    m = Mapper(chart, width, height, pad_bottom=76)
    p: List[str] = [
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" '
        'width="100%%" height="%d" font-family="system-ui,Segoe UI,sans-serif">'
        % (width, height, height),
        '<rect width="%d" height="%d" fill="#fff"/>' % (width, height),
        '<text x="%d" y="22" text-anchor="middle" font-size="15" '
        'font-weight="600" fill="#1c2430">%s</text>'
        % (width // 2, _esc(chart.title)),
    ]

    xt = nice_ticks(m.x0, m.x1)
    yt = nice_ticks(m.y0, m.y1)
    xstep = (xt[1] - xt[0]) if len(xt) > 1 else 1.0
    ystep = (yt[1] - yt[0]) if len(yt) > 1 else 1.0
    for t in xt:
        x = m.px(t)
        p.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                 'stroke="#eceff3"/>' % (x, m.top, x, m.bottom))
        p.append('<text x="%.1f" y="%.1f" text-anchor="middle" '
                 'font-size="11" fill="#5b6672">%s</text>'
                 % (x, m.bottom + 16, tick_label(t, xstep)))
    for t in yt:
        y = m.py(t)
        p.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                 'stroke="#eceff3"/>' % (m.left, y, m.right, y))
        p.append('<text x="%.1f" y="%.1f" text-anchor="end" font-size="11" '
                 'fill="#5b6672">%s</text>'
                 % (m.left - 8, y + 4, tick_label(t, ystep)))
    p.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" '
             'stroke="#c3ccd6"/>'
             % (m.left, m.top, m.plot_width, m.plot_height))
    p.append('<text x="%.1f" y="%.1f" text-anchor="middle" font-size="12" '
             'fill="#333">%s</text>'
             % ((m.left + m.right) / 2, m.bottom + 40, _esc(chart.x_label)))
    p.append('<text x="16" y="%.1f" text-anchor="middle" font-size="12" '
             'fill="#333" transform="rotate(-90 16 %.1f)">%s</text>'
             % ((m.top + m.bottom) / 2, (m.top + m.bottom) / 2,
                _esc(chart.y_label)))

    for value, label, colour in chart.h_lines:
        if not (m.y0 <= value <= m.y1):
            continue
        y = m.py(value)
        p.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                 'stroke-dasharray="4 4"/>' % (m.left, y, m.right, y, colour))
        p.append('<text x="%.1f" y="%.1f" text-anchor="end" font-size="11" '
                 'fill="%s">%s</text>'
                 % (m.right - 6, y - 6, colour, _esc(label)))

    for s in chart.series:
        pts = m.points(s.xs, s.ys)
        if len(pts) < 4:
            continue
        d = " ".join("%.1f,%.1f" % (pts[i], pts[i + 1])
                     for i in range(0, len(pts), 2))
        dash = (' stroke-dasharray="%s"' % " ".join(str(v) for v in s.dash)
                if s.dash else "")
        p.append('<polyline points="%s" fill="none" stroke="%s" '
                 'stroke-width="%d"%s/>' % (d, s.colour, s.width, dash))

    ly = m.top + 14
    for s in chart.series:
        p.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                 'stroke-width="%d"/>'
                 % (m.left + 12, ly, m.left + 38, ly, s.colour, s.width))
        p.append('<text x="%.1f" y="%.1f" font-size="11" fill="#333">%s</text>'
                 % (m.left + 44, ly + 4, _esc(s.label)))
        ly += 16
    if chart.footnote:
        p.append('<text x="%.1f" y="%d" font-size="11" fill="#7a838d">%s</text>'
                 % (m.left, height - 10, _esc(chart.footnote)))
    p.append("</svg>")
    return "".join(p)


def three_view_svg(ac, width: int = 900, height: int = 620) -> str:
    """Plan and side views as SVG, both on one scale, with a metre rule."""
    plan, side, _ext = three_view(ac)
    pad = 30
    panel_h = (height - 3 * pad) // 2

    def extent(shapes):
        xs = [q[0] for sh in shapes for q in sh.points]
        ys = [q[1] for sh in shapes for q in sh.points]
        return min(xs), max(xs), min(ys), max(ys)

    scales = []
    for shapes in (plan, side):
        x0, x1, y0, y1 = extent(shapes)
        scales.append(min((width - 2 * pad) / max(x1 - x0, 1e-6),
                          panel_h / max(y1 - y0, 1e-6)))
    scale = min(scales) * 0.92

    p = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" '
         'width="100%%" height="%d" font-family="system-ui,sans-serif">'
         % (width, height, height),
         '<rect width="%d" height="%d" fill="#fff"/>' % (width, height)]

    for idx, (shapes, title) in enumerate(((plan, "PLAN"), (side, "SIDE"))):
        x0, x1, y0, y1 = extent(shapes)
        oy = pad + idx * (panel_h + pad)
        cx = pad + ((width - 2 * pad) - (x1 - x0) * scale) / 2
        cy = oy + (panel_h - (y1 - y0) * scale) / 2
        p.append('<text x="%d" y="%.1f" font-size="11" font-weight="600" '
                 'fill="#6b7580">%s</text>' % (pad, oy - 8, title))
        for sh in shapes:
            pts = " ".join("%.1f,%.1f" % (cx + (q[0] - x0) * scale,
                                          cy + (q[1] - y0) * scale)
                           for q in sh.points)
            if sh.fill and len(sh.points) >= 3:
                p.append('<polygon points="%s" fill="%s" stroke="%s" '
                         'stroke-width="%d"/>'
                         % (pts, sh.fill, sh.outline, sh.width))
            else:
                p.append('<polyline points="%s" fill="none" stroke="%s" '
                         'stroke-width="%d"/>' % (pts, sh.outline, sh.width))
        if idx == 0:
            span_px = (y1 - y0) * scale
            len_px = (x1 - x0) * scale
            p.append(_dim_v_svg(cx - 14, cy, cy + span_px,
                                "span %.0f in" % (ac.wing.span_m * 39.3701)))
            p.append(_dim_h_svg(cy + span_px + 16, cx, cx + len_px,
                                "length %.0f in"
                                % (ac.body.fuse_length_m * 39.3701)))

    bar_m = 0.5
    while bar_m * scale > (width - 2 * pad) * 0.35:
        bar_m *= 0.5
    bx, by = pad, height - 12
    p.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#333" '
             'stroke-width="2"/>' % (bx, by, bx + bar_m * scale, by))
    p.append('<text x="%.1f" y="%.1f" font-size="11" fill="#333">'
             '%.2f m  (%.1f in)</text>'
             % (bx + bar_m * scale + 8, by + 4, bar_m, bar_m * 39.3701))
    p.append("</svg>")
    return "".join(p)


def _dim_h_svg(y, x0, x1, text):
    return ('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#98a2ad"/>'
            '<text x="%.1f" y="%.1f" text-anchor="middle" font-size="11" '
            'fill="#6b7580">%s</text>'
            % (x0, y, x1, y, (x0 + x1) / 2, y + 14, _esc(text)))


def _dim_v_svg(x, y0, y1, text):
    return ('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#98a2ad"/>'
            '<text x="%.1f" y="%.1f" text-anchor="middle" font-size="11" '
            'fill="#6b7580" transform="rotate(-90 %.1f %.1f)">%s</text>'
            % (x, y0, x, y1, x - 6, (y0 + y1) / 2, x - 6, (y0 + y1) / 2,
               _esc(text)))
