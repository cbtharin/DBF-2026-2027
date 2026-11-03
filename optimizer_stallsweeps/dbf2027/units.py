"""Unit conversions.  All internal computation in this package is strictly SI
(metre, kilogram, second, newton).  Conversions live here so that the physics
modules never carry mixed units.
"""
from __future__ import annotations

# length
IN2M = 0.0254
M2IN = 1.0 / IN2M
FT2M = 0.3048
M2FT = 1.0 / FT2M

# mass / force
LB2KG = 0.45359237
KG2LB = 1.0 / LB2KG
OZ2KG = 0.028349523125
KG2OZ = 1.0 / OZ2KG
G0 = 9.80665                      # m/s^2
LBF2N = LB2KG * G0
N2LBF = 1.0 / LBF2N

# speed
MPH2MS = 0.44704
MS2MPH = 1.0 / MPH2MS
KT2MS = 0.514444

# area
IN2_2_M2 = IN2M ** 2
FT2_2_M2 = FT2M ** 2

# energy
WH2J = 3600.0
J2WH = 1.0 / WH2J


def kg2lbf(m_kg: float) -> float:
    """Mass in kg -> weight in pound-force at standard gravity."""
    return m_kg * KG2LB


def lbf2kg(w_lb: float) -> float:
    return w_lb * LB2KG


def n2lb(force_n: float) -> float:
    return force_n * N2LBF
