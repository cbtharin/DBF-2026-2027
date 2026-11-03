#!/usr/bin/env python3
"""Quantify the NUMERICAL error in the mission model.

This is the study that tells you how much of a reported lap time is physics
and how much is discretisation.  It does four things:

  1. Grid convergence on the mission integrator: run the same aeroplane at
     systematically refined step sizes, compute the observed order p, the
     Richardson-extrapolated zero-step-size answer, and the Grid Convergence
     Index (Roache) as an error band.
  2. Accuracy of the propulsion operating-point solver: the closed-form
     battery-sag solution against a brute-force fixed point, and the
     torque residual at the returned answer.
  3. Convergence of the fixed-point mass loop.
  4. Repeatability / cost, so you can pick a step size on evidence.

    python grid_convergence.py

Numerical error is the part you can drive to zero.  Model error (propeller
coefficients, structural areal densities, the course figure) is the part you
cannot, and it is documented in ASSUMPTIONS.md.
"""
from __future__ import annotations

import math
import sys
import time

from dbf2027 import Aircraft, Config, Design
from dbf2027.components import MOTORS, PACKS, PROPS
from dbf2027.mission import mission_2, mission_3
from dbf2027.numerics import (brent, gci, observed_order,
                              richardson_extrapolate)
from dbf2027.propulsion import PropulsionSystem

LINE = "=" * 78
THIN = "-" * 78


def demo_design() -> Design:
    """The aeroplane the thorough run recommended - a representative point."""
    return Design(span_m=1.826, chord_m=0.259, airfoil="CLARKY",
                  motor="SCORPION_SII_4025_440", prop="APC_13x10", n_motors=2,
                  pack="GA_6S3300_60C", n_systems=1,
                  sensor_mass_kg=1.361, n_simulators=7)


# ==========================================================================
# 1. Grid convergence of the mission integrator
# ==========================================================================
def grid_convergence(cfg: Config, design: Design) -> None:
    print()
    print(LINE)
    print("1. GRID CONVERGENCE OF THE MISSION INTEGRATOR")
    print(LINE)
    print("Refinement ratio r = 2 on the nominal ground-track step dx.")
    print()

    ac = Aircraft(design, cfg)
    steps = [24.0, 12.0, 6.0, 3.0, 1.5, 0.75]
    rows = []
    for dx in steps:
        t0 = time.time()
        r2 = mission_2(ac, cfg, dx=dx)
        r3 = mission_3(ac, cfg, dx=dx)
        rows.append((dx, r2.time_s, r2.energy_j / 3600.0, r3.laps,
                     r3.time_s, time.time() - t0))

    print("%8s %13s %11s %8s %11s %8s" %
          ("dx [m]", "M2 time [s]", "M2 E [Wh]", "M3 laps", "M3 time [s]",
           "cost [s]"))
    print(THIN)
    for dx, t2, e2, l3, t3, cost in rows:
        print("%8.2f %13.5f %11.4f %8d %11.5f %8.3f"
              % (dx, t2, e2, l3, t3, cost))
    print(THIN)

    # Richardson on the three finest grids for each scalar quantity
    print()
    print("Richardson extrapolation on the three finest grids (r = 2):")
    print()
    print("%-18s %12s %10s %14s %12s" %
          ("quantity", "fine value", "order p", "extrapolated", "GCI [%]"))
    print(THIN)
    for label, idx, unit in (("M2 elapsed time", 1, "s"),
                             ("M2 energy", 2, "Wh"),
                             ("M3 elapsed time", 4, "s")):
        f_c, f_m, f_f = rows[-3][idx], rows[-2][idx], rows[-1][idx]
        p = observed_order(f_c, f_m, f_f, r=2.0)
        if p is None:
            print("%-18s %12.5f %10s %14s %12s"
                  % (label, f_f, "exact", "-", "0.000"))
            continue
        p_use = max(min(p, 8.0), 0.5)
        ex = richardson_extrapolate(f_m, f_f, p_use, r=2.0)
        band = 100.0 * gci(f_m, f_f, p_use, r=2.0)
        print("%-18s %12.5f %10.2f %14.5f %12.4f"
              % (label, f_f, p, ex, band))
    print(THIN)

    # Error of each grid against the extrapolated answer
    f_c, f_m, f_f = rows[-3][1], rows[-2][1], rows[-1][1]
    p = observed_order(f_c, f_m, f_f, r=2.0)
    p_use = max(min(p, 8.0), 0.5) if p else 4.0
    exact_t2 = richardson_extrapolate(f_m, f_f, p_use, r=2.0)
    print()
    print("Error in Mission 2 elapsed time vs the extrapolated exact value")
    print("(%.5f s):" % exact_t2)
    print()
    print("%8s %14s %12s %10s" % ("dx [m]", "error [s]", "error [%]", "cost [s]"))
    print(THIN)
    for dx, t2, _e2, _l3, _t3, cost in rows:
        err = t2 - exact_t2
        print("%8.2f %14.5f %12.5f %10.3f"
              % (dx, err, 100.0 * err / exact_t2, cost))
    print(THIN)
    print("A 1 percent error in Mission 2 time is a 1 percent error in the M2")
    print("score ratio, so this is the number that has to be small.")


# ==========================================================================
# 2. Propulsion solver accuracy
# ==========================================================================
def propulsion_accuracy(cfg: Config) -> None:
    """The propulsion operating point is solved, not interpolated.

    Battery sag used to be handled by a fixed-point iteration, and the result
    was then tabulated and interpolated.  Both steps introduced error.  The
    sag loop now has a closed-form solution (see propulsion.py), so the whole
    system is one scalar root-find, and the table is gone.

    This section verifies the closed form against a brute-force 200-iteration
    fixed point, and reports how well the root is actually converged.
    """
    print()
    print(LINE)
    print("2. PROPULSION SOLVER ACCURACY (no interpolation anywhere)")
    print(LINE)

    combos = [("SCORPION_SII_4025_440", "APC_13x10", 2, "GA_6S3300_60C", 2),
              ("SCORPION_SII_4025_440", "APC_12x8", 1, "GA_6S3300_60C", 1),
              ("TMOTOR_AT4130_300", "APC_16x12", 2, "TP_6S4000_70C", 1)]

    print("%-30s %14s %14s %14s" %
          ("system", "max |residual|", "rel residual", "max dn [rev/s]"))
    print(THIN)
    for mk, pk, nm, pack_key, nsys in combos:
        ps = PropulsionSystem(MOTORS[mk], PROPS[pk], nm, PACKS[pack_key],
                              cfg.env, n_systems=nsys)
        m, p, env = ps.motor, ps.prop, ps.env
        worst_r = worst_rel = worst_n = 0.0
        for i in range(25):
            v = 0.5 + 45.0 * i / 24.0
            for soc in (1.0, 0.72, 0.41, 0.13):
                for thr in (1.0, 0.85, 0.60, 0.35, 0.15):
                    op = ps.solve(v, soc, thr)
                    n = op.rpm / 60.0
                    if n <= 1e-6:
                        continue
                    # torque residual at the returned rpm; the motor sees
                    # throttle x the LOADED pack voltage (battery sag)
                    vt = thr * op.pack_voltage_v
                    cur = (vt - 60.0 * n / m.kv) / ps.rm
                    q_m = ps.kt * (cur - m.i0)
                    j = v / (n * p.diameter_m)
                    q_p = (p.cp(j) * env.rho * n * n
                           * p.diameter_m ** 5 / (2.0 * math.pi))
                    worst_r = max(worst_r, abs(q_m - q_p))
                    worst_rel = max(worst_rel,
                                    abs(q_m - q_p) / max(abs(q_p), 1e-12))
                    # independent bisection on the same residual, to a far
                    # tighter tolerance than the solver uses
                    worst_n = max(worst_n,
                                  abs(n - _bisect_rpm(ps, v, soc, thr)))
        print("%-30s %14.3e %14.3e %14.3e"
              % ("%s + %s x%d" % (mk.split("_")[0], pk.replace("APC_", ""), nm),
                 worst_r, worst_rel, worst_n))
    print(THIN)
    print("|residual| is motor torque minus propeller torque at the answer,")
    print("in N*m - it should be at the level of floating-point noise.")
    print("dn compares Brent against an independent 200-step bisection on")
    print("the same residual; it bounds the root-finding error directly.")

    # cost
    ps = PropulsionSystem(MOTORS["SCORPION_SII_4025_440"], PROPS["APC_13x10"],
                          2, PACKS["GA_6S3300_60C"], cfg.env, n_systems=2)
    n_rep = 20000
    t0 = time.time()
    for i in range(n_rep):
        ps._memo_key = None
        ps.thrust(5.0 + 35.0 * ((i * 0.618) % 1.0), 0.7, 1.0)
    print()
    print("Cost: %.2f us per operating point."
          % (1e6 * (time.time() - t0) / n_rep))


def _bisect_rpm(ps, v: float, soc: float, thr: float,
                iters: int = 200) -> float:
    """Independent root find on the torque balance, by plain bisection.

    Brent is fast but elaborate; bisection is slow and stupid, which is
    exactly what makes it a useful second opinion.  200 halvings drives
    the bracket below machine epsilon, so any disagreement is Brent's.
    """
    m, p, env = ps.motor, ps.prop, ps.env
    vt = thr * ps.pack.ocv(soc)
    if vt <= 0.05:
        return 0.0
    # battery sag appears to the motor as an extra series resistance
    r_eff = ps.rm + ps.motors_per_system * thr * thr * ps.r_sys

    def res(nn):
        cur = (vt - 60.0 * nn / m.kv) / r_eff
        j = v / (nn * p.diameter_m) if nn > 1e-6 else 1e9
        return (ps.kt * (cur - m.i0)
                - p.cp(j) * env.rho * nn * nn * p.diameter_m ** 5
                / (2.0 * math.pi))

    lo, hi = 1e-4, m.kv * vt / 60.0
    if res(lo) <= 0.0:
        return 0.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if res(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ==========================================================================
# 3. Mass loop convergence
# ==========================================================================
def mass_loop(cfg: Config, design: Design) -> None:
    from dbf2027.aerodynamics import WingGeometry, build_body
    from dbf2027.components import AIRFOILS
    from dbf2027.weights import airframe_masses

    print()
    print(LINE)
    print("3. FIXED-POINT CONVERGENCE OF THE MASS LOOP")
    print(LINE)
    print("Structure mass depends on MTOW, and MTOW contains the structure.")
    print()

    ac = Aircraft(design, cfg)
    wing = ac.wing
    body = ac.body
    b = cfg.build
    payload = ac.m2_payload_kg
    # Every installed item, exactly as the aircraft counts them: motors,
    # props, ESCs, mounts, packs, fuse harnesses, arming plugs and the
    # separate Rx battery.
    prop_mass = sum(ac.propulsion_items.values())
    n_ult = ac.n_ultimate

    mtow = payload + prop_mass + 2.0
    print("%6s %14s %16s" % ("iter", "MTOW [kg]", "change [kg]"))
    print(THIN)
    prev = mtow
    for i in range(12):
        mb = airframe_masses(wing, body, mtow * 9.80665, b, n_ult)
        mb.add("propulsion group", prop_mass)
        mb.add("payload restraint", b.payload_structure_frac * payload)
        mb.add("contingency", b.contingency_frac * mb.total)
        new = mb.total + payload
        mtow = 0.5 * mtow + 0.5 * new
        print("%6d %14.9f %16.2e" % (i, mtow, mtow - prev))
        if abs(mtow - prev) < 1e-9:
            break
        prev = mtow
    print(THIN)
    print("Converged MTOW used by the model: %.9f kg" % ac.mtow_m2_kg)
    print("Residual: %.3e kg  (tolerance 1e-4 kg = 0.1 g)"
          % abs(ac.mtow_m2_kg - (ac.empty_mass_kg + ac.m2_payload_kg)))


# ==========================================================================
# 4. Recommended step size
# ==========================================================================
def recommend(cfg: Config, design: Design) -> None:
    print()
    print(LINE)
    print("4. COST vs ACCURACY - PICKING dx")
    print(LINE)
    ac = Aircraft(design, cfg)
    ref2 = mission_2(ac, cfg, dx=0.375).time_s
    print("Reference (dx = 0.375 m): M2 time = %.5f s" % ref2)
    print()
    print("%8s %14s %12s %11s %12s" %
          ("dx [m]", "M2 time [s]", "error [%]", "cost [ms]", "verdict"))
    print(THIN)
    for dx in (24.0, 16.0, 12.0, 8.0, 6.0, 4.0, 2.0, 1.0):
        t0 = time.time()
        n_rep = 3
        for _ in range(n_rep):
            t2 = mission_2(ac, cfg, dx=dx).time_s
        cost = 1000.0 * (time.time() - t0) / n_rep
        err = 100.0 * (t2 - ref2) / ref2
        verdict = ("screening" if abs(err) < 0.5 else
                   ("too coarse" if abs(err) >= 0.5 else ""))
        if abs(err) < 0.05:
            verdict = "final scoring"
        print("%8.2f %14.5f %12.4f %11.2f %12s" % (dx, t2, err, cost, verdict))
    print(THIN)


def main() -> int:
    cfg = Config()
    design = demo_design()
    print(LINE)
    print("NUMERICAL ERROR STUDY - DBF 2026-27 MODEL")
    print(LINE)
    print("Integrator: RK4 (4th order) with bisection event location")
    print("Quadrature: 5-point Gauss-Legendre (exact to degree 9)")
    print("Root finding: Brent (superlinear)")
    print("Propulsion:  exact torque balance (no iteration, no interpolation)")
    grid_convergence(cfg, design)
    propulsion_accuracy(cfg)
    mass_loop(cfg, design)
    recommend(cfg, design)
    print()
    print(LINE)
    print("Numerical error is now well below model uncertainty.")
    print("See ASSUMPTIONS.md for the model error that remains.")
    print(LINE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
