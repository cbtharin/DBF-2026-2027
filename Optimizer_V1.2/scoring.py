"""Mission scoring per the 2026-27 rules.

    M1 = 1.0                      (successful flight, no payload)
    M2 = 1 + N(weight/time) / Max(weight/time)
    M3 = 2 + N(laps * sensor weight) / Max(laps * sensor weight)
    GM = 0.5 + N(sensor weight * drop height) / Max(...)

    Total Mission Score = M1 + M2 + M3 + GM, so the ceiling is 7.5, not
    6.0.  The Ground Mission is not a flight but it is scored the same way
    and it counts the same, and it is the third place sensor weight pays -
    which is why the optimiser drives the sensor as heavy as the aeroplane
    can carry.

Normalisation
-------------
`Max(...)` is the best score in the field, which nobody knows in advance.
Two modes:

  "self"  (default) - normalise against the best design found in this
          simulation.  The best design scores exactly 1.0 on that ratio by
          construction, and every other design is reported as a fraction of
          it.  This answers "which of our candidates is best, and by how
          much", which is the question a design team actually has.

  "field" - normalise against an estimate of the competition
          (`rules.field_best_*`).  Use this when you know what the field can
          do and want an absolute score.  The ratio is capped at 1.0, because
          Max(...) includes you: beating your estimate just makes you the
          new Max.

Two things about the formulas drive the whole design either way:

1.  A completed Mission 3 is worth at least 2.0 points before the ratio term
    even starts.  A failed mission scores 0, not a small number.  Reliability
    beats cleverness, so the optimiser scores failures as 0.

2.  Raw units are lbf per second for M2 and lap*lbf for M3.  Both are ratios,
    so any consistent unit works; lb/s is what the scoring sheet shows.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

from .aircraft import Aircraft
from .config import Config
from .mission import MissionResult, mission_1, mission_2, mission_3
from .units import KG2LB

# Weight on performance beyond the reference.  In "field" mode the score
# saturates at the reference, but margin past it is insurance against having
# underestimated the competition, so the optimiser still values it weakly.
INSURANCE_WEIGHT = 0.05


@dataclass
class ScoreCard:
    feasible: bool
    m1: float
    m2: float
    m3: float
    total_mission: float
    m2_raw_lb_per_s: float
    m3_raw_lap_lb: float
    gm: float = 0.0
    gm_raw_lb_in: float = 0.0
    gm_ratio: float = 0.0
    gm_saturated: bool = False
    r1: Optional[MissionResult] = None
    r2: Optional[MissionResult] = None
    r3: Optional[MissionResult] = None
    reason: str = ""
    m2_ratio: float = 0.0          # raw / reference, uncapped
    m3_ratio: float = 0.0
    m2_saturated: bool = False     # only meaningful in "field" mode
    m3_saturated: bool = False
    mode: str = "self"

    @property
    def objective(self) -> float:
        """What the optimiser maximises (M1 is a constant for any flyer).

        The Ground Mission is included: it is part of the total score and
        it pulls the design the same way Mission 3 does, towards a heavy
        sensor in a container that survives a 60 inch drop.
        """
        excess = (max(0.0, self.m2_ratio - 1.0)
                  + max(0.0, self.m3_ratio - 1.0)
                  + max(0.0, self.gm_ratio - 1.0))
        return self.m2 + self.m3 + self.gm + INSURANCE_WEIGHT * excess


def _score_from_ratios(m1: float, r2: float, r3: float, mode: str,
                       rg: float = 0.0):
    """Apply the formulas.  In field mode every ratio caps at 1.0.

    RULE 6.3 - Total Mission Score = M1 + M2 + M3 + GM, ceiling 7.5.
    RULE 3.3.3a - a failed mission scores 0, not a small number.
    """
    if mode == "field":
        gm = (0.5 + min(rg, 1.0)) if rg > 0 else 0.0
        m2 = (1.0 + min(r2, 1.0)) if r2 > 0 else 0.0
        m3 = (2.0 + min(r3, 1.0)) if r3 > 0 else 0.0
    else:
        gm = (0.5 + rg) if rg > 0 else 0.0
        m2 = (1.0 + r2) if r2 > 0 else 0.0
        m3 = (2.0 + r3) if r3 > 0 else 0.0
    return m2, m3, gm, m1 + m2 + m3 + gm


def ground_mission(ac, cfg: Config):
    """RULE 3.3.4 - six drops of the loaded container onto concrete.

    Returns (ok, raw_lb_in, reason).

    It is a drop test, not a flight, so there is no trajectory to
    integrate.  What decides it is whether the container's foam can hold
    the sensor under the deceleration it survives, from the height the
    team declared.  Both of those are already in the payload model, so
    this is a comparison, not a new simulation.

    The score is sensor weight x drop height, and the sensor must be at
    its declared maximum, so this pulls in exactly the same direction as
    Mission 3: build the heaviest sensor the aeroplane can still fly.
    """
    from .payload import crush_stroke_m
    r = cfg.rules
    h_in = min(max(ac.design.drop_height_in, 0.0), r.max_drop_height_in)
    if h_in <= 0.0:
        return False, 0.0, "no drop height declared"

    # Foam crushes usefully through about 70 percent of its thickness.
    stroke_needed = crush_stroke_m(h_in, r.sensor_survivable_g)
    stroke_available = 0.70 * ac.payload_geom.pad_m
    if stroke_available < stroke_needed:
        return (False, 0.0,
                "container foam is %.0f mm, needs %.0f mm to hold the "
                "sensor under %.0f g from %.0f in"
                % (1000 * ac.payload_geom.pad_m,
                   1000 * stroke_needed / 0.70, r.sensor_survivable_g, h_in))

    # RULE 3.1.3d - the rigging ring to the bottom of the container must
    # clear 15 in, or the fixture cannot reach the declared height.
    container_in = ac.payload_geom.box_h / 0.0254
    if container_in > r.max_rigging_plus_container_in:
        return (False, 0.0,
                "container is %.1f in deep; with rigging it exceeds the "
                "%.0f in the drop fixture allows"
                % (container_in, r.max_rigging_plus_container_in))

    return True, ac.sensor_max_kg * KG2LB * h_in, ""


def evaluate_scores(ac: Aircraft, cfg: Config, dx: float = 6.0,
                    ref_m2: Optional[float] = None,
                    ref_m3: Optional[float] = None,
                    ref_gm: Optional[float] = None) -> ScoreCard:
    """Fly all three missions, run the Ground Mission, and score them.

    `ref_m2` / `ref_m3` / `ref_gm` are the normalisers.  In "self" mode the
    optimiser passes the best raw values found so far; they are only a
    scale for the objective, and the final report renormalises against the
    finished population.
    """
    r = cfg.rules
    mode = r.normalisation
    fb2 = ref_m2 if ref_m2 is not None else r.field_best_m2_lb_per_s
    fb3 = ref_m3 if ref_m3 is not None else r.field_best_m3_lap_lb
    fbg = ref_gm if ref_gm is not None else r.field_best_gm_lb_in
    fb2 = max(fb2, 1e-9)
    fb3 = max(fb3, 1e-9)
    fbg = max(fbg, 1e-9)

    checks = ac.static_checks()
    if not checks["ok"]:
        return ScoreCard(False, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         reason="; ".join(checks["failures"]), mode=mode)

    # RULE 3.3.3c - one Aircraft object flies all three missions, so the
    # airframe, power system and propeller cannot differ between them.
    # Only the payload and the Mission 3 sensor weight change, which is
    # exactly what the rules allow.
    r1 = mission_1(ac, cfg, dx)
    r2 = mission_2(ac, cfg, dx)
    r3 = mission_3(ac, cfg, dx)

    m1 = 1.0 if r1.ok else 0.0                       # RULE 3.3.3.1

    if r2.ok and r2.time_s > 0.0:
        # RULE 3.3.3.2 - scored payload weight over flight time.
        raw2 = (ac.m2_scored_kg * KG2LB) / r2.time_s
        ratio2 = raw2 / fb2
    else:
        raw2, ratio2 = 0.0, 0.0

    if r3.ok and r3.laps > 0:
        # RULE 3.3.3.3, 3.1.1c - Mission 3 flies the sensor weight
        # actually chosen for it, which may be below the declared maximum.
        raw3 = r3.laps * (ac.sensor_m3_kg * KG2LB)
        ratio3 = raw3 / fb3
    else:
        raw3, ratio3 = 0.0, 0.0

    gm_ok, raw_g, gm_reason = ground_mission(ac, cfg)
    ratio_g = (raw_g / fbg) if gm_ok else 0.0

    m2, m3, gm, total = _score_from_ratios(m1, ratio2, ratio3, mode, ratio_g)
    reasons = [x.reason for x in (r1, r2, r3) if not x.ok and x.reason]
    if not gm_ok and gm_reason:
        reasons.append("ground mission: " + gm_reason)
    return ScoreCard(r1.ok and r2.ok and r3.ok and gm_ok, m1, m2, m3, total,
                     raw2, raw3, gm, raw_g, ratio_g,
                     ratio_g > 1.0, r1, r2, r3, "; ".join(reasons),
                     ratio2, ratio3, ratio2 > 1.0, ratio3 > 1.0, mode)


def rescore(card: ScoreCard, ref_m2: float, ref_m3: float,
            mode: Optional[str] = None,
            ref_gm: Optional[float] = None) -> ScoreCard:
    """Re-apply the ratio terms with different normalisers, without
    re-running the simulation."""
    mode = mode or card.mode
    r2 = card.m2_raw_lb_per_s / max(ref_m2, 1e-9) if card.m2_raw_lb_per_s > 0 else 0.0
    r3 = card.m3_raw_lap_lb / max(ref_m3, 1e-9) if card.m3_raw_lap_lb > 0 else 0.0
    rg = (card.gm_raw_lb_in / max(ref_gm, 1e-9)
          if (ref_gm and card.gm_raw_lb_in > 0) else card.gm_ratio)
    m2, m3, gm, total = _score_from_ratios(card.m1, r2, r3, mode, rg)
    return ScoreCard(card.feasible, card.m1, m2, m3, total,
                     card.m2_raw_lb_per_s, card.m3_raw_lap_lb,
                     gm, card.gm_raw_lb_in, rg, rg > 1.0,
                     card.r1, card.r2, card.r3, card.reason,
                     r2, r3, r2 > 1.0, r3 > 1.0, mode)


def population_references(rows: Sequence[Dict]) -> tuple:
    """Best raw M2, M3 and Ground Mission across a set of designs."""
    best2 = max([r.get("m2_raw_lb_per_s", 0.0) for r in rows] or [0.0])
    best3 = max([r.get("m3_raw_lap_lb", 0.0) for r in rows] or [0.0])
    bestg = max([r.get("gm_raw_lb_in", 0.0) for r in rows] or [0.0])
    return max(best2, 1e-9), max(best3, 1e-9), max(bestg, 1e-9)


def apply_references(rows: List[Dict], ref_m2: float, ref_m3: float,
                     ref_gm: float, mode: str = "self") -> None:
    """Rescore every row against the given references, in place.

    Rows carry raw scores AND an objective computed from whatever
    references were current when they were evaluated.  Mixing rows scored
    against different references in one sorted list compares numbers that
    do not share a scale - the smaller the references, the larger every
    ratio, so a row scored against generous references outranks a better
    design scored against strict ones.  Anything that combines two
    batches of rows must put them on the same scale first.
    """
    for row in rows:
        r2 = row.get("m2_raw_lb_per_s", 0.0) / max(ref_m2, 1e-9)
        r3 = row.get("m3_raw_lap_lb", 0.0) / max(ref_m3, 1e-9)
        rg = row.get("gm_raw_lb_in", 0.0) / max(ref_gm, 1e-9)
        m1 = row.get("m1", 0.0)
        m2, m3, gm, total = _score_from_ratios(m1, r2, r3, mode, rg)
        row["m2"], row["m3"], row["gm"] = m2, m3, gm
        row["total_mission"] = total
        row["m2_ratio"], row["m3_ratio"], row["gm_ratio"] = r2, r3, rg
        excess = (max(0.0, r2 - 1.0) + max(0.0, r3 - 1.0)
                  + max(0.0, rg - 1.0))
        row["objective"] = ((m2 + m3 + gm + INSURANCE_WEIGHT * excess)
                            if row.get("feasible") else -1e9)


def normalise_population(rows: List[Dict], mode: str = "self") -> tuple:
    """Rescore every row against the best in the population, in place.

    This is what makes "self" mode honest: the normaliser is the finished
    population, not a moving target from mid-search.
    """
    ref2, ref3, refg = population_references(rows)
    apply_references(rows, ref2, ref3, refg, mode)
    return ref2, ref3, refg


# ==========================================================================
# Explaining a score, step by step
# ==========================================================================
def score_explanation(ac, cfg, card: ScoreCard,
                      ref_m2: Optional[float] = None,
                      ref_m3: Optional[float] = None,
                      ref_gm: Optional[float] = None) -> List[Dict]:
    """The scoring arithmetic as a list of blocks, for a UI to render.

    Every front end wants the same thing here - the formula as the rules
    write it, then the same formula with this aeroplane's numbers in it,
    then the result - so it is built once, as data, and the tkinter and
    web front ends both format it.

    Each block is {"title", "formula", "rows": [(label, value, note)],
    "result", "note"}.
    """
    r = cfg.rules
    mode = card.mode
    if ref_m2 is None:
        ref_m2 = r.field_best_m2_lb_per_s
    if ref_m3 is None:
        ref_m3 = r.field_best_m3_lap_lb
    if ref_gm is None:
        ref_gm = r.field_best_gm_lb_in
    ref_label = ("the best design in this run" if mode == "self"
                 else "your estimate of the field's best")

    blocks: List[Dict] = []

    # ---------------- Mission 1 ----------------
    r1 = card.r1
    blocks.append({
        "title": "Mission 1 - flight",
        "formula": "M1 = 1.0 for a successful flight, 0 otherwise",
        "rows": [
            ("laps flown", "%d" % (r1.laps if r1 else 0), "3 required"),
            ("completed", "yes" if (r1 and r1.ok) else "no",
             "" if (r1 and r1.ok) else (r1.reason if r1 else "not flown")),
        ],
        "result": "M1 = %.3f" % card.m1,
        "note": "No payload. It is pass/fail, so it cannot be optimised - "
                "it only has to be achieved.",
    })

    # ---------------- Mission 2 ----------------
    r2 = card.r2
    t2 = r2.time_s if r2 else 0.0
    w2 = ac.m2_scored_kg * KG2LB
    raw2 = card.m2_raw_lb_per_s
    rows2 = [
        ("sensor (declared max)", "%.3f lb" % (ac.sensor_max_kg * KG2LB),
         "M2 must fly the max"),
        ("container tare", "%.3f lb"
         % (ac.payload_geom.container_tare_kg * KG2LB), "fixed 0.5 lb"),
        ("loaded container", "%.3f lb" % (ac.loaded_container_kg * KG2LB),
         "sensor + container"),
        ("%d simulators x %.3f lb" % (ac.design.n_simulators,
                                      ac.sim_mass_kg * KG2LB),
         "%.3f lb" % (ac.design.n_simulators * ac.sim_mass_kg * KG2LB),
         "each must match the loaded container (3.1.3)"),
        ("payload weight  N", "%.3f lb" % w2, "what M2 scores"),
        ("time for %d laps  t" % r.m2_required_laps, "%.1f s" % t2,
         "landing is NOT in this"),
        ("raw  N / t", "%.5f lb/s" % raw2, ""),
        ("reference  Max(N/t)", "%.5f lb/s" % ref_m2, ref_label),
        ("ratio", "%.4f" % card.m2_ratio,
         "capped at 1.0" if mode == "field" else "uncapped in self mode"),
    ]
    blocks.append({
        "title": "Mission 2 - delivery flight",
        "formula": "M2 = 1 + N(weight/time) / Max(weight/time)",
        "rows": rows2,
        "result": "M2 = 1 + %s = %.3f"
                  % (("min(%.4f, 1)" % card.m2_ratio) if mode == "field"
                     else "%.4f" % card.m2_ratio, card.m2),
        "note": ("SATURATED - this design meets or beats the reference, so "
                 "the ratio is held at 1.0 and the score cannot rise "
                 "further." if card.m2_saturated and mode == "field" else
                 "Five laps inside the 300 s window, carrying the sensor "
                 "container and every simulator."),
    })

    # ---------------- Mission 3 ----------------
    r3 = card.r3
    laps3 = r3.laps if r3 else 0
    sens_lb = ac.sensor_m3_kg * KG2LB
    rows3 = [
        ("laps flown  #", "%d" % laps3,
         "as many as fit in %.0f s" % r.mission_window_s),
        ("sensor weight", "%.3f lb" % sens_lb,
         "may be under the declared max (3.1.1)"),
        ("declared max", "%.3f lb" % (ac.sensor_max_kg * KG2LB),
         "%.0f%% of it is flown" % (100.0 * sens_lb
                                    / max(ac.sensor_max_kg * KG2LB, 1e-9))),
        ("tow line", "%.2f m" % ac.tow_line_length_m,
         "min 1.5 x span (3.1.2), and it is drag"),
        ("raw  # x weight", "%.3f lap*lb" % card.m3_raw_lap_lb, ""),
        ("reference  Max(...)", "%.3f lap*lb" % ref_m3, ref_label),
        ("ratio", "%.4f" % card.m3_ratio,
         "capped at 1.0" if mode == "field" else "uncapped in self mode"),
    ]
    blocks.append({
        "title": "Mission 3 - sensor flight",
        "formula": "M3 = 2 + N(#laps x sensor weight) / Max(#laps x sensor weight)",
        "rows": rows3,
        "result": "M3 = 2 + %s = %.3f"
                  % (("min(%.4f, 1)" % card.m3_ratio) if mode == "field"
                     else "%.4f" % card.m3_ratio, card.m3),
        "note": ("SATURATED - already at the ceiling against this reference."
                 if card.m3_saturated and mode == "field" else
                 "The sensor is deployed and the score multiplies laps by "
                 "its weight, so a heavy sensor and a fast aeroplane both "
                 "pay - which is the whole tension in this design."),
    })

    # ---------------- Ground Mission ----------------
    gm_ok, _raw, gm_reason = ground_mission(ac, cfg)
    h_in = min(ac.design.drop_height_in, r.max_drop_height_in)
    rows_g = [
        ("declared drop height", "%.0f in" % h_in,
         "whole inches, max %.0f" % r.max_drop_height_in),
        ("sensor weight", "%.3f lb" % (ac.sensor_max_kg * KG2LB),
         "declared max, same as M2"),
        ("container foam", "%.0f mm per face" % (1000 * ac.payload_geom.pad_m),
         "holds %.0f g through 70%% crush" % r.sensor_survivable_g),
        ("drops", "6", "two on each of three faces, onto concrete"),
        ("survives", "yes" if gm_ok else "NO", gm_reason),
        ("raw  weight x height", "%.1f lb*in" % card.gm_raw_lb_in, ""),
        ("reference  Max(...)", "%.1f lb*in" % ref_gm, ref_label),
        ("ratio", "%.4f" % card.gm_ratio,
         "capped at 1.0" if mode == "field" else "uncapped in self mode"),
    ]
    blocks.append({
        "title": "Ground Mission - drop test",
        "formula": "GM = 0.5 + N(sensor weight x drop height) / "
                   "Max(weight x height)",
        "rows": rows_g,
        "result": "GM = 0.5 + %s = %.3f"
                  % (("min(%.4f, 1)" % card.gm_ratio) if mode == "field"
                     else "%.4f" % card.gm_ratio, card.gm),
        "note": ("SATURATED - already at the ceiling against this reference."
                 if card.gm_saturated and mode == "field" else
                 "Not a flight, but it is scored and it counts. This is the "
                 "THIRD place sensor weight pays, which is why the optimum "
                 "sensor is as heavy as the aeroplane can still fly."),
    })

    # ---------------- total ----------------
    blocks.append({
        "title": "Total",
        "formula": "Total Mission Score = M1 + M2 + M3 + GM",
        "rows": [("M1", "%.3f" % card.m1, "max 1.0"),
                 ("M2", "%.3f" % card.m2, "max 2.0"),
                 ("M3", "%.3f" % card.m3, "max 3.0"),
                 ("GM", "%.3f" % card.gm, "max 1.5")],
        "result": "TOTAL = %.3f" % card.total_mission,
        "note": "Ceiling is 7.500: 1 + 2 + 3 + 1.5, reached only by a design "
                "that is the best in the field at all three scored events. "
                "The competition score is Total Report Score x Total Mission "
                "Score + Participation, so this number is a multiplier.",
    })
    return blocks


def score_explanation_text(ac, cfg, card: ScoreCard,
                           ref_m2: Optional[float] = None,
                           ref_m3: Optional[float] = None,
                           ref_gm: Optional[float] = None,
                           width: int = 74) -> str:
    """The same explanation, formatted for a terminal or a Text widget."""
    out = []
    for b in score_explanation(ac, cfg, card, ref_m2, ref_m3, ref_gm):
        out.append("=" * width)
        out.append(b["title"].upper())
        out.append("=" * width)
        out.append("  " + b["formula"])
        out.append("")
        for label, value, note in b["rows"]:
            out.append("    %-34s %14s   %s" % (label, value, note))
        out.append("")
        out.append("  " + b["result"])
        if b["note"]:
            out.append("")
            for line in _wrap(b["note"], width - 4):
                out.append("  " + line)
        out.append("")
    return "\n".join(out)


def _wrap(text: str, width: int) -> List[str]:
    words, lines, cur = text.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > width:
            lines.append(cur)
            cur = w
        else:
            cur = (cur + " " + w) if cur else w
    if cur:
        lines.append(cur)
    return lines
