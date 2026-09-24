"""Component catalogues: airfoils, motors, propellers, cells, ESCs.

IMPORTANT
---------
The numbers below are *representative* values for real, commonly-available
DBF-class hardware.  They are good enough to rank configurations, and they are
NOT a substitute for the manufacturer data sheet or your own dyno / thrust-stand
run.  Before you cut aluminium, replace the entries you actually intend to buy
with measured values (see `load_motors_csv` / `load_props_csv`).

Every record carries `source` so you can see what still needs verifying.
"""
from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

from .units import IN2M


# ==========================================================================
# Airfoils
# ==========================================================================
@dataclass(frozen=True)
class Airfoil:
    name: str
    cl_max_2d: float        # section clmax at Re ~ 300k
    cd_min: float           # section minimum profile drag at re_ref
    cl_at_cd_min: float
    k_drag: float           # cd = cd_min + k*(cl-cl_opt)^2
    cl_break: float         # cl beyond which drag rises sharply
    t_over_c: float
    cm_ac: float
    re_ref: float = 300_000.0
    note: str = ""


AIRFOILS: Dict[str, Airfoil] = {
    "S1223": Airfoil("S1223", 2.20, 0.0215, 1.05, 0.0130, 1.85, 0.121, -0.28,
                     note="Selig high-lift. Huge CLmax, draggy, sharp stall. Heavy-haul choice."),
    "E423": Airfoil("E423", 1.95, 0.0165, 0.90, 0.0105, 1.65, 0.125, -0.24,
                    note="Eppler high-lift. Friendlier than S1223, still payload-oriented."),
    "S1210": Airfoil("S1210", 2.05, 0.0195, 1.00, 0.0120, 1.75, 0.120, -0.26,
                     note="High lift, slightly cleaner than S1223."),
    "SD7062": Airfoil("SD7062", 1.55, 0.0122, 0.60, 0.0082, 1.30, 0.140, -0.14,
                      note="Thick, docile, good spar depth. Common DBF compromise."),
    "SD7037": Airfoil("SD7037", 1.35, 0.0102, 0.45, 0.0072, 1.15, 0.092, -0.10,
                      note="Low drag, sport. Favours lap-count missions."),
    "CLARKY": Airfoil("CLARKY", 1.42, 0.0118, 0.40, 0.0085, 1.20, 0.117, -0.08,
                      note="Flat bottom, trivially buildable, forgiving."),
    "MH32": Airfoil("MH32", 1.28, 0.0092, 0.35, 0.0068, 1.05, 0.087, -0.06,
                    note="Fast, thin, low drag. Needs a real spar."),
    "NACA4412": Airfoil("NACA4412", 1.48, 0.0112, 0.45, 0.0080, 1.25, 0.120, -0.09,
                        note="Classic, well documented, easy to validate in XFOIL."),
    # --- added to broaden the trade study to 20 sections ----------------
    "FX63137": Airfoil("FX63137", 1.92, 0.0185, 0.95, 0.0112, 1.62, 0.137, -0.25,
                       note="Wortmann FX 63-137. High lift, thick, popular for heavy-lift."),
    "GOE398": Airfoil("GOE398", 1.78, 0.0170, 0.85, 0.0108, 1.50, 0.110, -0.20,
                      note="Gottingen 398. High lift, classic, easy to build."),
    "S1218": Airfoil("S1218", 1.88, 0.0180, 0.92, 0.0115, 1.58, 0.113, -0.22,
                     note="Selig high-lift, a little cleaner than S1223."),
    "NACA6412": Airfoil("NACA6412", 1.62, 0.0135, 0.62, 0.0092, 1.38, 0.120, -0.14,
                        note="High camber NACA. Good lift, heavy pitching moment."),
    "E214": Airfoil("E214", 1.38, 0.0108, 0.48, 0.0076, 1.18, 0.111, -0.09,
                    note="Eppler general purpose. Balanced, forgiving."),
    "E193": Airfoil("E193", 1.32, 0.0100, 0.42, 0.0072, 1.12, 0.104, -0.07,
                    note="Eppler sport. Low drag, mild stall, easy to build."),
    "E205": Airfoil("E205", 1.26, 0.0094, 0.38, 0.0068, 1.06, 0.106, -0.06,
                    note="Eppler low drag. Fast, needs speed to work."),
    "S3021": Airfoil("S3021", 1.30, 0.0098, 0.40, 0.0070, 1.10, 0.092, -0.07,
                     note="Selig sport. Good low-Re behaviour for its drag."),
    "SA7038": Airfoil("SA7038", 1.36, 0.0104, 0.45, 0.0074, 1.16, 0.093, -0.09,
                      note="Somers/Selig. Similar to SD7037, slightly more lift."),
    "RG15": Airfoil("RG15", 1.18, 0.0088, 0.32, 0.0064, 0.98, 0.089, -0.05,
                    note="Rolf Girsberger. Fast sailplane section, very low drag."),
    "SD8020": Airfoil("SD8020", 1.06, 0.0080, 0.22, 0.0060, 0.88, 0.080, -0.02,
                      note="Thin, near-symmetric, lowest drag here. Speed only."),
    "AG35": Airfoil("AG35", 1.04, 0.0078, 0.25, 0.0062, 0.86, 0.073, -0.04,
                    note="Drela thin section. Minimum drag, minimum structure depth."),
}


# ==========================================================================
# Motors (brushless outrunners)
# ==========================================================================
@dataclass(frozen=True)
class Motor:
    key: str
    name: str
    kv: float               # rpm per volt, no load
    rm: float               # ohm, effective line resistance used by the model
    i0: float               # A, no-load current
    mass_kg: float
    i_max_a: float          # manufacturer continuous current
    p_max_w: float          # manufacturer continuous power
    shaft_mm: float
    source: str = "catalogue estimate - VERIFY"


MOTORS: Dict[str, Motor] = {
    "SCORPION_SII_4020_420": Motor("SCORPION_SII_4020_420", "Scorpion SII-4020-420KV",
                                   420, 0.033, 1.6, 0.405, 70, 1900, 6.0),
    "SCORPION_SII_3020_890": Motor("SCORPION_SII_3020_890", "Scorpion SII-3020-890KV",
                                   890, 0.052, 1.4, 0.205, 50, 1000, 5.0),
    "HACKER_A50_14L": Motor("HACKER_A50_14L", "Hacker A50-14L V4",
                            355, 0.041, 1.5, 0.425, 60, 1700, 8.0),
    "HACKER_A40_12L": Motor("HACKER_A40_12L", "Hacker A40-12L V4",
                            610, 0.055, 1.3, 0.222, 45, 1000, 5.0),
    "COBRA_C4120_12": Motor("COBRA_C4120_12", "Cobra C-4120/12 465KV",
                            465, 0.036, 1.5, 0.305, 58, 1450, 6.0),
    "COBRA_C3520_14": Motor("COBRA_C3520_14", "Cobra C-3520/14 780KV",
                            780, 0.048, 1.2, 0.190, 44, 900, 5.0),
    "TMOTOR_AT4130_300": Motor("TMOTOR_AT4130_300", "T-Motor AT4130 300KV",
                               300, 0.030, 1.8, 0.408, 75, 2000, 8.0),
    "TMOTOR_AT4125_540": Motor("TMOTOR_AT4125_540", "T-Motor AT4125 540KV",
                               540, 0.038, 1.6, 0.345, 65, 1550, 6.0),
    "TMOTOR_AT2814_900": Motor("TMOTOR_AT2814_900", "T-Motor AT2814 900KV",
                               900, 0.065, 1.1, 0.152, 36, 700, 4.0),
    "SUNNYSKY_X4112S_400": Motor("SUNNYSKY_X4112S_400", "SunnySky X4112S 400KV",
                                 400, 0.035, 1.5, 0.300, 60, 1400, 6.0),
    "NEU_1512_1Y": Motor("NEU_1512_1Y", "Neu 1512/1.5Y (direct drive wind)",
                         420, 0.025, 2.0, 0.390, 80, 2200, 6.0),
    "EFLITE_POWER_32": Motor("EFLITE_POWER_32", "E-flite Power 32 770KV",
                             770, 0.058, 1.3, 0.200, 42, 850, 5.0),
    # --- added to broaden the trade study to 20 motors ------------------
    # Smaller motors matter now that 3 and 4 motor layouts are searched.
    "SCORPION_SII_4025_330": Motor("SCORPION_SII_4025_330", "Scorpion SII-4025-330KV",
                                   330, 0.028, 1.9, 0.500, 85, 2300, 6.0),
    "HACKER_A60_7XS": Motor("HACKER_A60_7XS", "Hacker A60-7XS V4",
                            365, 0.026, 2.2, 0.640, 95, 2800, 8.0),
    "TMOTOR_AT3520_880": Motor("TMOTOR_AT3520_880", "T-Motor AT3520 880KV",
                               880, 0.058, 1.3, 0.195, 42, 850, 5.0),
    "TMOTOR_AT2820_1050": Motor("TMOTOR_AT2820_1050", "T-Motor AT2820 1050KV",
                                1050, 0.072, 1.0, 0.128, 32, 600, 4.0),
    "COBRA_C2826_10": Motor("COBRA_C2826_10", "Cobra C-2826/10 1130KV",
                            1130, 0.078, 1.0, 0.118, 30, 560, 4.0),
    "COBRA_C3530_14": Motor("COBRA_C3530_14", "Cobra C-3530/14 750KV",
                            750, 0.045, 1.3, 0.225, 48, 1000, 5.0),
    "SUNNYSKY_X3520_520": Motor("SUNNYSKY_X3520_520", "SunnySky X3520 520KV",
                                520, 0.042, 1.4, 0.238, 50, 1100, 5.0),
    "AXI_2826_10": Motor("AXI_2826_10", "AXi 2826/10 920KV",
                         920, 0.068, 1.1, 0.177, 35, 700, 5.0),
}


# ==========================================================================
# Propellers
# ==========================================================================
@dataclass(frozen=True)
class Propeller:
    key: str
    diameter_m: float
    pitch_m: float
    mass_kg: float
    ct0: float              # static thrust coefficient
    cp0: float              # static power coefficient
    j_zero: float           # advance ratio at zero thrust
    kp: float = 0.75        # Cp roll-off: Cp(j_zero) = (1-kp) * cp0
    source: str = "parametric p/D fit - REPLACE with UIUC/APC data"

    @property
    def p_over_d(self) -> float:
        return self.pitch_m / self.diameter_m

    def ct(self, j: float) -> float:
        return max(0.0, self.ct0 * (1.0 - j / self.j_zero))

    def cp(self, j: float) -> float:
        x = min(j / self.j_zero, 1.25)
        return max(0.20 * self.cp0, self.cp0 * (1.0 - self.kp * x))

    def thrust(self, rho: float, n_rps: float, v: float) -> float:
        if n_rps <= 1e-6:
            return 0.0
        j = v / (n_rps * self.diameter_m)
        return self.ct(j) * rho * n_rps ** 2 * self.diameter_m ** 4

    def power(self, rho: float, n_rps: float, v: float) -> float:
        if n_rps <= 1e-6:
            return 0.0
        j = v / (n_rps * self.diameter_m)
        return self.cp(j) * rho * n_rps ** 3 * self.diameter_m ** 5


def _make_prop(d_in: float, p_in: float) -> Propeller:
    """Parametric APC-thin-electric-like coefficients from pitch/diameter.

    Anchored on three things a real propeller has to satisfy at once:

      * static figure of merit  FOM = ct0^1.5 / (1.2533 cp0)  lands at ~0.71
        for a low-pitch prop falling to ~0.52 for a high-pitch one;
      * peak propulsive efficiency  eta_max = 0.444 * j_zero * ct0/cp0  runs
        ~0.70 (p/D 0.5) to ~0.78 (p/D 0.83);
      * peak efficiency occurs near 0.8x geometric pitch speed, i.e. at
        J/j_zero = 2/3 with the linear Cp roll-off above.

    Check: 13x8 at 8000 rpm gives 7.8 lbf static on 700 W (FOM 0.63), peak
    eta 0.74 at 22 m/s.  Fine for ranking; measure the two or three props you
    shortlist before you trust an absolute number.
    """
    pod = p_in / d_in
    ct0 = 0.115 + 0.030 * pod
    ratio = min(2.70, max(1.40, 2.48 - 2.27 * (pod - 0.50)))
    cp0 = ct0 / ratio
    jz = 1.15 * pod + 0.06
    mass = 0.0009 * (d_in ** 2.2) / 10.0 + 0.012
    return Propeller(key="APC_%gx%g" % (d_in, p_in), diameter_m=d_in * IN2M,
                     pitch_m=p_in * IN2M, mass_kg=mass,
                     ct0=ct0, cp0=cp0, j_zero=jz)


_PROP_SIZES = [
    (10, 5), (10, 7), (11, 5.5), (11, 7), (11, 8),
    (12, 6), (12, 8), (12, 10), (13, 6.5), (13, 8), (13, 10),
    (14, 7), (14, 8.5), (14, 10), (14, 12),
    (15, 8), (15, 10), (15, 13),
    (16, 8), (16, 10), (16, 12),
    (17, 10), (17, 12), (18, 8), (18, 10), (18, 12),
    (20, 10), (20, 13),
]

PROPS: Dict[str, Propeller] = {}
for _d, _p in _PROP_SIZES:
    _pr = _make_prop(_d, _p)
    PROPS[_pr.key] = _pr


# ==========================================================================
# Batteries - commercial off-the-shelf packs only
# ==========================================================================
@dataclass(frozen=True)
class BatteryPack:
    """An unaltered, commercially procured pack, as its label reads.

    RULE 3.2.3a, 3.2.3.1b.  The rules forbid custom packs and forbid
    wiring packs in series or parallel, so a pack is an atomic item: you buy this exact
    thing or you do not use it.  Everything here is what the manufacturer's
    label has to state - chemistry, capacity, voltage and C rating.

    VERIFY every entry against the current manufacturer datasheet before
    ordering.  Masses and C ratings drift between production runs.
    """
    key: str
    name: str                  # what the label says
    chemistry: str             # "LiPo" or "LiIon" - lithium only
    series: int                # cells in series, as labelled
    capacity_ah: float
    c_rating: float            # labelled continuous discharge C
    mass_kg: float             # manufacturer figure, pack only
    v_nom_cell: float
    v_full_cell: float
    v_min_cell: float
    usable_fraction: float
    source: str = "manufacturer catalogue - VERIFY"

    @property
    def v_nominal(self) -> float:
        return self.series * self.v_nom_cell

    @property
    def v_full(self) -> float:
        return self.series * self.v_full_cell

    @property
    def v_min(self) -> float:
        return self.series * self.v_min_cell

    @property
    def watt_hours(self) -> float:
        """Rated capacity x rated voltage - RULE 3.2.3.2b, 3.2.3.2c define
        both the airframe limit and the FAA hand-carry limit this way."""
        return self.v_nominal * self.capacity_ah

    @property
    def i_max_a(self) -> float:
        """Labelled maximum continuous discharge, capacity x C rating.

        RULE 3.2.3.2d sizes the fuse from exactly this number."""
        return self.capacity_ah * self.c_rating

    def ocv(self, soc: float) -> float:
        """Open-circuit pack voltage against state of charge."""
        s_ = max(0.0, min(1.0, soc))
        per_cell = (self.v_min_cell
                    + (self.v_full_cell - self.v_min_cell)
                    * (0.18 * s_ + 0.82 * (s_ ** 0.45)))
        return per_cell * self.series


def _lipo(key, name, series, mah, c_rating, grams):
    return BatteryPack(key, name, "LiPo", series, mah / 1000.0, c_rating,
                       grams / 1000.0, 3.70, 4.20, 3.35, 0.80)


def _liion(key, name, series, mah, c_rating, grams):
    return BatteryPack(key, name, "LiIon", series, mah / 1000.0, c_rating,
                       grams / 1000.0, 3.60, 4.20, 3.00, 0.85)


# Representative packs that are actually purchasable and sit under the
# 100 Wh per-pack limit.  Wh shown is rated capacity x rated voltage.
_PACK_LIST = (
    # --- 6S LiPo, the workhorse range -------------------------------------
    _lipo("TP_6S2200_25C", "Thunder Power TP2200-6SPX25", 6, 2200, 25, 352),
    _lipo("TP_6S3300_70C", "Thunder Power G8 6S 3300 70C", 6, 3300, 70, 542),
    _lipo("TP_6S4000_70C", "Thunder Power G8 6S 4000 70C", 6, 4000, 70, 655),
    _lipo("GA_6S3300_60C", "Gens Ace 6S 3300 60C", 6, 3300, 60, 556),
    _lipo("GA_6S4000_60C", "Gens Ace 6S 4000 60C", 6, 4000, 60, 668),
    _lipo("GA_6S4500_60C", "Gens Ace 6S 4500 60C", 6, 4500, 60, 748),
    _lipo("TUR_6S3000_65C", "Turnigy Graphene 6S 3000 65C", 6, 3000, 65, 500),
    _lipo("PU_6S2250_45C", "Pulse Ultra 6S 2250 45C", 6, 2250, 45, 368),
    # --- 5S and 4S, for lower-Kv or lower-voltage systems -----------------
    _lipo("GA_5S4000_60C", "Gens Ace 5S 4000 60C", 5, 4000, 60, 556),
    _lipo("TP_5S3300_70C", "Thunder Power G8 5S 3300 70C", 5, 3300, 70, 455),
    _lipo("GA_4S5000_60C", "Gens Ace 4S 5000 60C", 4, 5000, 60, 548),
    _lipo("GA_4S4000_60C", "Gens Ace 4S 4000 60C", 4, 4000, 60, 442),
    _lipo("TUR_4S5000_65C", "Turnigy Graphene 4S 5000 65C", 4, 5000, 65, 545),
    _lipo("GA_3S5000_60C", "Gens Ace 3S 5000 60C", 3, 5000, 60, 412),
    # --- Li-ion: far better Wh/kg, far worse current --------------------
    _liion("GA_6S3500_LI", "Gens Ace Li-ion 6S 3500 (21700)", 6, 3500, 10, 392),
    _liion("GA_6S4000_LI", "Gens Ace Li-ion 6S 4000 (21700)", 6, 4000, 10, 447),
    _liion("GA_4S6000_LI", "Gens Ace Li-ion 4S 6000 (21700)", 4, 6000, 10, 448),
)

PACKS = {p.key: p for p in _PACK_LIST}


def legal_packs(max_wh_total: float, n_systems: int = 1):   # RULE 3.2.3.2b
    """Packs that keep the whole aeroplane inside the energy limit.

    Total propulsion energy is the sum over packs, and every pack must also
    be under the FAA 100 Wh hand-carry limit on its own.
    """
    out = []
    for pk in PACKS.values():
        if pk.watt_hours > 100.0 + 1e-9:
            continue
        if pk.watt_hours * n_systems > max_wh_total + 1e-9:
            continue
        out.append(pk)
    return out


def load_packs_csv(path: str) -> int:
    """CSV columns: key,name,chemistry,series,capacity_mah,c_rating,mass_g

    Use this to put your own shortlist in, straight off the labels.
    """
    if not os.path.exists(path):
        return 0
    n = 0
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            chem = row.get("chemistry", "LiPo").strip()
            mk = _liion if chem.lower().startswith("lii") else _lipo
            pk = mk(row["key"], row.get("name", row["key"]),
                    int(row["series"]), float(row["capacity_mah"]),
                    float(row["c_rating"]), float(row["mass_g"]))
            PACKS[row["key"]] = pk
            n += 1
    return n


# ==========================================================================
# ESC
# ==========================================================================
@dataclass(frozen=True)
class ESC:
    key: str
    i_cont_a: float
    mass_kg: float
    efficiency: float
    max_cells: int


ESCS: Dict[str, ESC] = {
    "ESC_45A": ESC("ESC_45A", 45, 0.045, 0.960, 6),
    "ESC_60A": ESC("ESC_60A", 60, 0.063, 0.960, 6),
    "ESC_80A": ESC("ESC_80A", 80, 0.082, 0.955, 8),
    "ESC_100A": ESC("ESC_100A", 100, 0.105, 0.950, 8),
    "ESC_120A": ESC("ESC_120A", 120, 0.125, 0.950, 12),
}


def pick_esc(current_a: float, cells: int) -> Optional[ESC]:
    """Smallest ESC with 20 percent headroom that supports the cell count."""
    need = current_a * 1.20
    best = None
    for e in ESCS.values():
        if e.i_cont_a >= need and e.max_cells >= cells:
            if best is None or e.i_cont_a < best.i_cont_a:
                best = e
    return best


# ==========================================================================
# CSV overrides - measured data wins over the estimates above
# ==========================================================================
def load_motors_csv(path: str) -> int:
    """CSV columns: key,name,kv,rm,i0,mass_kg,i_max_a,p_max_w,shaft_mm"""
    if not os.path.exists(path):
        return 0
    n = 0
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            MOTORS[row["key"]] = Motor(
                row["key"], row.get("name", row["key"]), float(row["kv"]),
                float(row["rm"]), float(row["i0"]), float(row["mass_kg"]),
                float(row["i_max_a"]), float(row["p_max_w"]),
                float(row.get("shaft_mm", 5)),
                source="measured/" + os.path.basename(path))
            n += 1
    return n


def load_props_csv(path: str) -> int:
    """CSV columns: key,diameter_in,pitch_in,mass_kg,ct0,cp0,j_zero[,kp]

    Fit ct0/cp0/j_zero from your own thrust-stand sweep or from the UIUC
    propeller database (Brandt and Selig).  That is the single highest-value
    data upgrade you can make to this model.
    """
    if not os.path.exists(path):
        return 0
    n = 0
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            PROPS[row["key"]] = Propeller(
                row["key"], float(row["diameter_in"]) * IN2M,
                float(row["pitch_in"]) * IN2M, float(row["mass_kg"]),
                float(row["ct0"]), float(row["cp0"]), float(row["j_zero"]),
                float(row.get("kp", 0.80)),
                source="measured/" + os.path.basename(path))
            n += 1
    return n


# ==========================================================================
# Catalogue registry and user additions
# ==========================================================================
# Every catalogue is described here as data - its columns, how to build an
# entry from a row of text, and how to turn an entry back into one.  The
# GUI generates its table and its add/edit form from this, so a column
# added here appears in the interface without touching the interface, and
# a CSV written by hand and a row typed into the GUI go through exactly
# the same constructor.
#
# User additions live in `catalogues/*.csv` next to dbf.py and are loaded
# at package import, which matters: a search runs in worker processes that
# re-import the package, and a motor that existed only in the parent
# process would vanish halfway through a sweep.
@dataclass(frozen=True)
class Column:
    """One editable field of a catalogue entry."""
    name: str
    kind: type
    default: Any = ""
    help: str = ""


@dataclass(frozen=True)
class CatalogueSpec:
    key: str                       # "motors"
    label: str                     # "Motors"
    singular: str
    columns: Tuple[Column, ...]
    note: str = ""

    @property
    def store(self) -> Dict[str, Any]:
        return _STORES[self.key]()

    def build(self, row: Dict[str, str]):
        return _BUILDERS[self.key](row)

    def to_row(self, obj) -> Dict[str, str]:
        return _TO_ROW[self.key](obj)

    def fieldnames(self) -> List[str]:
        return [c.name for c in self.columns]


def _f(row, name, default=0.0):
    try:
        return float(row.get(name, default) or default)
    except (TypeError, ValueError):
        return float(default)


def _i(row, name, default=0):
    try:
        return int(float(row.get(name, default) or default))
    except (TypeError, ValueError):
        return int(default)


def _s(row, name, default=""):
    return str(row.get(name, default) or default)


_STORES = {"airfoils": lambda: AIRFOILS, "motors": lambda: MOTORS,
           "props": lambda: PROPS, "packs": lambda: PACKS,
           "escs": lambda: ESCS}

_BUILDERS = {
    "airfoils": lambda r: Airfoil(
        _s(r, "name"), _f(r, "cl_max_2d", 1.4), _f(r, "cd_min", 0.009),
        _f(r, "cl_at_cd_min", 0.3), _f(r, "k_drag", 0.02),
        _f(r, "cl_break", 1.0), _f(r, "t_over_c", 0.12),
        _f(r, "cm_ac", -0.08), _f(r, "re_ref", 300000.0),
        _s(r, "note", "user entry - VERIFY")),
    "motors": lambda r: Motor(
        _s(r, "key"), _s(r, "name") or _s(r, "key"), _f(r, "kv", 500),
        _f(r, "rm", 0.04), _f(r, "i0", 1.4), _f(r, "mass_kg", 0.25),
        _f(r, "i_max_a", 50), _f(r, "p_max_w", 1000),
        _f(r, "shaft_mm", 5), source="user entry - VERIFY"),
    "props": lambda r: Propeller(
        _s(r, "key"), _f(r, "diameter_in", 12) * IN2M,
        _f(r, "pitch_in", 8) * IN2M, _f(r, "mass_kg", 0.035),
        _f(r, "ct0", 0.11), _f(r, "cp0", 0.045), _f(r, "j_zero", 0.75),
        _f(r, "kp", 0.80), source="user entry - VERIFY"),
    "packs": lambda r: (
        _liion if _s(r, "chemistry", "LiPo").lower().startswith("lii")
        else _lipo)(
            _s(r, "key"), _s(r, "name") or _s(r, "key"), _i(r, "series", 6),
            _f(r, "capacity_mah", 3000), _f(r, "c_rating", 60),
            _f(r, "mass_g", 500)),
    "escs": lambda r: ESC(
        _s(r, "key"), _f(r, "i_cont_a", 60), _f(r, "mass_kg", 0.06),
        _f(r, "efficiency", 0.96), _i(r, "max_cells", 6)),
}

_TO_ROW = {
    "airfoils": lambda a: {
        "name": a.name, "cl_max_2d": a.cl_max_2d, "cd_min": a.cd_min,
        "cl_at_cd_min": a.cl_at_cd_min, "k_drag": a.k_drag,
        "cl_break": a.cl_break, "t_over_c": a.t_over_c, "cm_ac": a.cm_ac,
        "re_ref": a.re_ref, "note": a.note},
    "motors": lambda m: {
        "key": m.key, "name": m.name, "kv": m.kv, "rm": m.rm, "i0": m.i0,
        "mass_kg": m.mass_kg, "i_max_a": m.i_max_a, "p_max_w": m.p_max_w,
        "shaft_mm": m.shaft_mm},
    "props": lambda p: {
        "key": p.key, "diameter_in": p.diameter_m / IN2M,
        "pitch_in": p.pitch_m / IN2M, "mass_kg": p.mass_kg, "ct0": p.ct0,
        "cp0": p.cp0, "j_zero": p.j_zero, "kp": p.kp},
    "packs": lambda b: {
        "key": b.key, "name": b.name, "chemistry": b.chemistry,
        "series": b.series, "capacity_mah": b.capacity_ah * 1000,
        "c_rating": b.c_rating, "mass_g": b.mass_kg * 1000},
    "escs": lambda e: {
        "key": e.key, "i_cont_a": e.i_cont_a, "mass_kg": e.mass_kg,
        "efficiency": e.efficiency, "max_cells": e.max_cells},
}

CATALOGUES: Tuple[CatalogueSpec, ...] = (
    CatalogueSpec("airfoils", "Airfoils", "airfoil", (
        Column("name", str, "", "the dictionary key as well"),
        Column("cl_max_2d", float, 1.40, "section CLmax at Re 300k"),
        Column("cd_min", float, 0.0090, "minimum profile drag"),
        Column("cl_at_cd_min", float, 0.30, "CL at that minimum"),
        Column("k_drag", float, 0.020, "cd = cd_min + k(cl-cl_opt)^2"),
        Column("cl_break", float, 1.00, "CL where drag rises sharply"),
        Column("t_over_c", float, 0.12, "thickness, sets spar depth"),
        Column("cm_ac", float, -0.080, "pitching moment; drives trim"),
        Column("re_ref", float, 300000.0, "Reynolds the data is at"),
        Column("note", str, "", "free text shown in the report")),
        "From XFOIL or wind tunnel data at a Reynolds number near "
        "300k. t_over_c matters twice: drag and spar depth."),
    CatalogueSpec("motors", "Motors", "motor", (
        Column("key", str, "", "unique id"),
        Column("name", str, "", "as the manufacturer lists it"),
        Column("kv", float, 500.0, "rpm per volt, no load"),
        Column("rm", float, 0.040, "ohm; NOT optional, see propulsion.py"),
        Column("i0", float, 1.4, "A, no-load current"),
        Column("mass_kg", float, 0.250, ""),
        Column("i_max_a", float, 50.0, "manufacturer continuous"),
        Column("p_max_w", float, 1000.0, "manufacturer continuous"),
        Column("shaft_mm", float, 5.0, "")),
        "Rm is the winding resistance and the model cannot run without "
        "it: with Rm = 0 the torque balance has no solution."),
    CatalogueSpec("props", "Propellers", "propeller", (
        Column("key", str, "", "unique id, e.g. APC_13x8"),
        Column("diameter_in", float, 12.0, ""),
        Column("pitch_in", float, 8.0, ""),
        Column("mass_kg", float, 0.035, ""),
        Column("ct0", float, 0.110, "static thrust coefficient"),
        Column("cp0", float, 0.045, "static power coefficient"),
        Column("j_zero", float, 0.75, "advance ratio at zero thrust"),
        Column("kp", float, 0.80, "Cp roll-off with J")),
        "The propeller fit is the single largest source of model error. "
        "Measured ct0/cp0/j_zero from a thrust stand or the UIUC "
        "database are worth more here than anywhere else in the model."),
    CatalogueSpec("packs", "Battery packs", "pack", (
        Column("key", str, "", "unique id"),
        Column("name", str, "", "exactly as the label reads"),
        Column("chemistry", str, "LiPo", "LiPo or LiIon"),
        Column("series", int, 6, "cells in series, as labelled"),
        Column("capacity_mah", float, 3000.0, "rated capacity"),
        Column("c_rating", float, 60.0, "labelled continuous discharge"),
        Column("mass_g", float, 500.0, "manufacturer figure")),
        "COTS packs only (rule 3.2.3a). Copy the four numbers off the "
        "label: chemistry, mAh, voltage (via series) and C rating. Wh "
        "and the fuse rating are derived, so they cannot disagree."),
    CatalogueSpec("escs", "ESCs", "ESC", (
        Column("key", str, "", "unique id"),
        Column("i_cont_a", float, 60.0, "continuous current"),
        Column("mass_kg", float, 0.060, ""),
        Column("efficiency", float, 0.960, ""),
        Column("max_cells", int, 6, "maximum pack series count")),
        "Chosen automatically for a design: the smallest ESC with 20 "
        "percent headroom that supports the cell count."),
)

CATALOGUES_BY_KEY: Dict[str, CatalogueSpec] = {c.key: c for c in CATALOGUES}


def user_dir() -> str:
    """Where user-added catalogue entries live: `catalogues/` by dbf.py."""
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(os.path.dirname(here), "catalogues")


def user_csv_path(key: str) -> str:
    return os.path.join(user_dir(), "%s.csv" % key)


def user_rows(key: str) -> List[Dict[str, str]]:
    """Rows currently in the user CSV for one catalogue."""
    path = user_csv_path(key)
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return [dict(r) for r in csv.DictReader(fh)]


def write_user_rows(key: str, rows: List[Dict[str, Any]]) -> str:
    """Replace the user CSV for one catalogue.  Returns the path."""
    spec = CATALOGUES_BY_KEY[key]
    os.makedirs(user_dir(), exist_ok=True)
    path = user_csv_path(key)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=spec.fieldnames())
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in spec.fieldnames()})
    return path


def add_user_entry(key: str, row: Dict[str, Any]) -> str:
    """Add or replace one entry, persist it, and apply it now.

    The same row goes through the same constructor whether it arrives
    from the GUI or from a CSV someone wrote by hand, so the two can
    never drift.
    """
    spec = CATALOGUES_BY_KEY[key]
    obj = spec.build({k: str(v) for k, v in row.items()})
    ident = getattr(obj, "key", None) or getattr(obj, "name")
    rows = [r for r in user_rows(key)
            if (r.get("key") or r.get("name")) != ident]
    rows.append({k: row.get(k, "") for k in spec.fieldnames()})
    path = write_user_rows(key, rows)
    spec.store[ident] = obj
    return path


def remove_user_entry(key: str, ident: str) -> bool:
    """Remove a USER entry.  Built-in entries are left alone.

    Returns True if a row was removed.  Deleting a built-in would make
    the shipped catalogue depend on local state, so it is refused; the
    entry simply is not in the user CSV to remove.
    """
    rows = user_rows(key)
    kept = [r for r in rows if (r.get("key") or r.get("name")) != ident]
    if len(kept) == len(rows):
        return False
    write_user_rows(key, kept)
    reload_catalogues()
    return True


_BUILTIN: Dict[str, Dict[str, Any]] = {}


def is_user_entry(key: str, ident: str) -> bool:
    """True when `ident` came from a user CSV rather than the package."""
    return ident not in _BUILTIN.get(key, {})


def load_user_catalogues() -> Dict[str, int]:
    """Apply every `catalogues/*.csv`.  Returns counts per catalogue.

    Called at package import so that worker processes see the same
    catalogue as the parent: a search re-imports the package in every
    worker, and an entry that existed only in the parent would disappear
    partway through a sweep.
    """
    if not _BUILTIN:
        for spec in CATALOGUES:
            _BUILTIN[spec.key] = dict(spec.store)
    counts = {}
    for spec in CATALOGUES:
        n = 0
        for row in user_rows(spec.key):
            try:
                obj = spec.build(row)
            except Exception:
                continue
            ident = getattr(obj, "key", None) or getattr(obj, "name")
            if ident:
                spec.store[ident] = obj
                n += 1
        counts[spec.key] = n
    return counts


def reload_catalogues() -> Dict[str, int]:
    """Reset to the shipped catalogues, then re-apply the user CSVs."""
    if not _BUILTIN:
        return load_user_catalogues()
    for spec in CATALOGUES:
        spec.store.clear()
        spec.store.update(_BUILTIN[spec.key])
    return load_user_catalogues()


def load_airfoils_csv(path: str) -> int:
    """CSV columns: name,cl_max_2d,cd_min,cl_at_cd_min,k_drag,cl_break,
    t_over_c,cm_ac[,re_ref,note]"""
    if not os.path.exists(path):
        return 0
    n = 0
    spec = CATALOGUES_BY_KEY["airfoils"]
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            obj = spec.build(row)
            AIRFOILS[obj.name] = obj
            n += 1
    return n


def load_escs_csv(path: str) -> int:
    """CSV columns: key,i_cont_a,mass_kg,efficiency,max_cells"""
    if not os.path.exists(path):
        return 0
    n = 0
    spec = CATALOGUES_BY_KEY["escs"]
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            obj = spec.build(row)
            ESCS[obj.key] = obj
            n += 1
    return n


# Apply user additions as soon as the catalogues exist.
load_user_catalogues()
