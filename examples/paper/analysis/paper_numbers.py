#!/usr/bin/env python3
"""Numbers quoted in the text of Section 3 of the paper, computed from the stored case records.

usage: paper_numbers.py [--runs DIR] [--sections NAME,...] [matrix [aux]]

  --runs      directory holding the study directories                                     default: runs
  --sections  comma-separated subset of fields, energy, conservation, segments, gas_spring, cd_ratio, cd_details,
              cost, studies, axial                                                       default: all
  matrix aux  study directories of the parametric matrix and of the sensitivity studies   default: matrixC8ax
              studiesC8ax (reference mesh); "matrixC8 studiesC8" gives the values with 1 mm axial cells,
              "matrix studies" those of the medium mesh
"""

import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figdata import GAMMA, R_GAS, REF1, cd_orifice, configure, read_ts, resolved, summary  # noqa: E402

RUNS = "runs"               # directory holding the study directories (set in main)
MAT = "matrixC8ax"          # parametric matrix (set in main)
AUX = "studiesC8ax"         # sensitivity studies (set in main)
NS, VS = (2, 4, 6, 8, 10), ("0.025", "0.1", "0.4", "1.0", "2.0", "4.0")


def case(n, v, root=None):
    """Case directory of the matrix case with N orifices at the mean speed v (string, as in the case name)."""
    return os.path.join(RUNS, MAT if root is None else root, "N%d_v%s" % (n, v))


def at(ts, t):
    """Value of every column at the last sample not after t (end of the segment that ends at t)."""
    i = int(np.searchsorted(ts["time_s"], t + 1e-12, side="right")) - 1
    return {k: v[i] for k, v in ts.items()}


def u_ideal(cfg):
    """Ideal orifice velocity on the nominal orifice area."""
    d = cfg["_derived"]
    return d["mean_speed_mps"] * d["orifice_area_ratio"]


def fields_2ms():
    """Orifice velocities and chamber pressures at 2 m/s at the instants of the field figures."""
    print("== fields at 2 m/s (compression t/T = 0.25 0.5 0.75 1; rebound 1.25 1.5 1.75 2)")
    for n in (2, 6, 10):
        c = case(n, "2.0")
        cfg, ts = resolved(c), read_ts(c)
        T, uo = cfg["_derived"]["stroke_time_s"], u_ideal(cfg)
        for lab, fr in (("comp", (0.25, 0.5, 0.75, 1.0)), ("reb", (1.25, 1.5, 1.75, 2.0))):
            r = [at(ts, f * T) for f in fr]
            print("N%-2d %s Umax %s | Umean %s | Umean/uo %s | dp kPa %s | p_avg kPa %s | Ma_mean %s" % (
                n, lab, " ".join("%.0f" % x["Umax_orifice"] for x in r),
                " ".join("%.1f" % x["U_orifice_mean_mps"] for x in r),
                " ".join("%.3f" % (x["U_orifice_mean_mps"] / uo) for x in r),
                " ".join("%.2f" % (x["dp_12_Pa"] / 1e3) for x in r),
                " ".join("%.1f" % (x["p_avg"] / 1e3) for x in r),
                " ".join("%.3f" % x["Ma_orifice_mean"] for x in r)))
        sel = ts["time_s"] <= T
        dp = ts["p1"] - ts["p2"]
        print("     max(p1-p2) compression %.2f kPa; rebound p2-p1 at 1.25T..2T: %s kPa; p_avg(2T)/p0-1 = %.3f%%; "
              "p_avg(T) %.1f kPa" % (np.nanmax(dp[sel]) / 1e3,
                                    " ".join("%.1f" % (-at(ts, f * T)["p1"] / 1e3 + at(ts, f * T)["p2"] / 1e3)
                                             for f in (1.25, 1.5, 1.75, 2.0)),
                                    100 * (ts["p_avg"][-1] / cfg["gas"]["initial_pressure_Pa"] - 1),
                                    at(ts, T)["p_avg"] / 1e3))


def energy():
    """Energy balance of one cycle: internal-energy rise at the initial volume against the dissipated energy."""
    print("== energy balance")
    rows = []
    for n in NS:
        for v in VS:
            c = case(n, v)
            cfg, s, ts = resolved(c), summary(c), read_ts(c)
            dU = (ts["p_avg"][-1] - cfg["gas"]["initial_pressure_Pa"]) * ts["V_analytic_m3"][0] / (GAMMA - 1)
            rows.append((n, v, dU, s["E_diss_J"]))
    for n, v, dU, E in rows:
        if v == "2.0":
            print("  N%d v2: dU %.4f E_d %.4f diff %.4f J" % (n, dU, E, dU - E))
    big = [(dU / E - 1) * 100 for n, v, dU, E in rows if E > 0.01]
    print("  E_d>0.01 J: %d cases, dU/E_d-1 = %.2f .. %.2f %%" % (len(big), min(big), max(big)))
    Es = [E for *_, E in rows]
    print("  E_d range %.3g .. %.3g J" % (min(Es), max(Es)))


def conservation():
    """Mass changes at the mesh switches and over the cycle, volume error, switch transient of the force."""
    print("== conservation")
    jumps, drifts, verr = [], [], []
    for n in NS:
        for v in VS:
            s = summary(case(n, v))
            jumps.append((max(map(abs, s["mass_jumps_rel"])), n, v))
            drifts.append((s["mass_drift_max_rel"], n, v))
            verr.append(s["volume_error_max_rel"])
    print("  max jump %.2g (N%d v%s); max drift %.2g (N%d v%s); median drift %.2g; V error %.4f..%.4f" % (
        *max(jumps), *max(drifts), float(np.median([d[0] for d in drifts])), min(verr), max(verr)))
    # switch transient of the piston force: largest deviation from a straight line fitted to the
    # preceding window, within five acoustic transit times after each switch
    for n, v in ((2, "1.0"), (2, "4.0")):
        c = case(n, v)
        cfg, s, ts = resolved(c), summary(c), read_ts(c)
        tau = s["acoustic_time_s"]
        amp = 0.0
        for seg in cfg["_derived"]["schedule"]["segments"][1:]:
            t0 = seg["start"]
            pre = (ts["time_s"] > t0 - 5 * tau) & (ts["time_s"] <= t0)
            post = (ts["time_s"] > t0) & (ts["time_s"] <= t0 + 5 * tau)
            if pre.sum() < 3 or post.sum() < 1:
                continue
            k = np.polyfit(ts["time_s"][pre], ts["Fx"][pre], 1)
            amp = max(amp, float(np.max(np.abs(ts["Fx"][post] - np.polyval(k, ts["time_s"][post])))))
        print("  N%d v%s: switch transient %.3f N, F_damp %.4f N, L/c0 %.3f ms" % (
            n, v, amp, s["F_damp_mean_N"], 1e3 * tau))


def segments_study():
    """Number of sub-meshes per stroke, K = 3 and 6 against the matrix case (K = 4)."""
    print("== sub-meshes per stroke (N2, 1 m/s)")
    ref = summary(case(2, "1.0"))
    for name in ("seg_K3", "seg_K6"):
        p = os.path.join(RUNS, AUX, name)
        if not os.path.exists(os.path.join(p, "results", "summary.json")):
            continue
        s = summary(p)
        print("  %s F %.4f (%+.2f%%) E %.4f (%+.2f%%) jump %.2g drift %.2g switches %d" % (
            name, s["F_damp_mean_N"], 100 * (s["F_damp_mean_N"] / ref["F_damp_mean_N"] - 1), s["E_diss_J"],
            100 * (s["E_diss_J"] / ref["E_diss_J"] - 1), max(map(abs, s["mass_jumps_rel"])),
            s["mass_drift_max_rel"], len(s["mass_jumps_rel"])))
    print("  K=4 F %.4f E %.4f jump %.2g drift %.2g switches %d" % (
        ref["F_damp_mean_N"], ref["E_diss_J"], max(map(abs, ref["mass_jumps_rel"])), ref["mass_drift_max_rel"],
        len(ref["mass_jumps_rel"])))


def gas_spring():
    """Polytropic exponent of the low-speed cases."""
    ns = [summary(case(n, v))["n_polytropic"] for n in NS for v in ("0.025", "0.1")]
    print("== gas spring: n = %.4f .. %.4f (v <= 0.1, %d cases)" % (min(ns), max(ns), len(ns)))


def cd_and_ratio():
    """Effective discharge coefficient, damping-to-spring ratio against Eq. (ratio), pre-charge, turbulence model,
    sinusoidal motion and time-step studies."""
    print("== discharge coefficient and damping ratio")
    rows = []
    for n in NS:
        for v in VS:
            c = case(n, v)
            cfg, s = resolved(c), summary(c)
            d, gas = cfg["_derived"], cfg["gas"]
            T0, p0 = gas["initial_temperature_K"], gas["initial_pressure_Pa"]
            c0 = math.sqrt(GAMMA * R_GAS * T0)
            mu0 = gas["sutherland_As"] * T0 ** 1.5 / (T0 + gas["sutherland_Ts_K"])
            uo = u_ideal(cfg)
            re = p0 / (R_GAS * T0) * uo * 2e-3 * cfg["geometry"]["orifice"]["radius_mm"] / mu0
            afam = d["areas_mm2"]["piston_front"] / d["areas_mm2"]["rod"]
            rows.append(dict(n=n, v=v, re=re, ma=uo / c0, cd=s["Cd_effective"], ratio=s["damping_to_spring_ratio"],
                             afam=afam, F=s["F_damp_mean_N"], E=s["E_diss_J"]))
    mid = [r for r in rows if 0.02 < r["ma"] < 0.3]
    cd = float(np.mean([r["cd"] for r in mid]))
    print("  Cd mean %.4f +- %.4f (%d cases, 0.02<Ma<0.3); Re range of these %.0f..%.0f" % (
        cd, float(np.std([r["cd"] for r in mid])), len(mid), min(r["re"] for r in mid), max(r["re"] for r in mid)))
    dev = [100 * (r["ratio"] / (0.5 * GAMMA * r["afam"] * (r["ma"] / cd) ** 2) - 1) for r in mid]
    print("  ratio vs Eq.(ratio) with mean Cd: max |dev| %.1f%%, mean %.1f%%" % (max(map(abs, dev)), np.mean(dev)))
    for r in sorted(rows, key=lambda r: r["re"]):
        print("   N%-2d v%-5s Re %8.1f Ma %.4f Cd %.4f F %.4g E %.4g" % (r["n"], r["v"], r["re"], r["ma"], r["cd"],
                                                                     r["F"], r["E"]))
    by = {(r["n"], r["v"]): r for r in rows}
    a, b = by[(2, "1.0")], by[(4, "2.0")]
    print("  pair N2 1 m/s / N4 2 m/s: Ma %.4f %.4f, ratio diff %.1f%%" % (a["ma"], b["ma"],
                                                                         100 * (b["ratio"] / a["ratio"] - 1)))
    print("  max Ma %.3f; Cd at max Ma %.4f" % (max(r["ma"] for r in rows), max(rows, key=lambda r: r["ma"])["cd"]))
    # aux: 20 bar, laminar, sinusoidal
    for name, n, v in (("p20_N2_v1.0", 2, "1.0"), ("p20_N6_v1.0", 6, "1.0"), ("p20_N10_v1.0", 10, "1.0")):
        p = os.path.join(RUNS, AUX, name)
        if not os.path.exists(os.path.join(p, "results", "summary.json")):
            continue
        s, r = summary(p), by[(n, v)]
        d_ratio = s["damping_to_spring_ratio"] / r["ratio"] - 1
        d_cd = (r["cd"] / s["Cd_effective"]) ** 2 - 1
        print("  %s ratio %+.2f%%, (Cd1/Cd20)^2-1 %+.2f%%, diff %.2f%%; Cd20 %.4f" % (
            name, 100 * d_ratio, 100 * d_cd, 100 * (d_ratio - d_cd), s["Cd_effective"]))
    for name, n, v in (("lam_N2_v0.1", 2, "0.1"), ("lam_N2_v1.0", 2, "1.0")):
        p = os.path.join(RUNS, AUX, name)
        if not os.path.exists(os.path.join(p, "results", "summary.json")):
            continue
        s, r = summary(p), by[(n, v)]
        print("  %s F %+.2f%%, Cd %.4f vs %.4f, Re %.0f" % (name, 100 * (s["F_damp_mean_N"] / r["F"] - 1),
                                                         s["Cd_effective"], r["cd"], r["re"]))
    p = os.path.join(RUNS, AUX, "sin_N4_v1.0")
    if os.path.exists(os.path.join(p, "results", "summary.json")):
        print("  sinusoidal N4 1 m/s: F ratio %.3f" % (summary(p)["F_damp_mean_N"] / by[(4, "1.0")]["F"]))
    for name, n, v in (("co_0.5", 2, "1.0"), ("co05_N2_v1.0", 2, "1.0"), ("co05_N2_v4.0", 2, "4.0")):
        p = os.path.join(RUNS, AUX, name)
        if not os.path.exists(os.path.join(p, "results", "summary.json")):
            continue
        s, r = summary(p), by[(n, v)]
        print("  %s: F %+.2f%%, Cd %+.2f%%, E %+.2f%%" % (
            name, 100 * (s["F_damp_mean_N"] / r["F"] - 1), 100 * (s["Cd_effective"] / r["cd"] - 1),
            100 * (s["E_diss_J"] / r["E"] - 1)))


def cost():
    """Cells, time steps, wall time, cores and core-hours of the matrix; core-hours of the other studies."""
    print("== cost")
    for n in NS:
        cells, steps, wall, cores = [], [], [], set()
        for v in VS:
            c = case(n, v)
            mr = json.load(open(os.path.join(c, "mesh", "mesh_report.json")))
            cells += [m["hexahedra"] for m in mr]
            s = summary(c)
            steps.append(s["solver_steps"])
            wall.append(s["solver_wall_s"] / 3600)
            cores.add(resolved(c)["run"]["n_procs"])
        print("  N%-2d cells %.0f..%.0f k, steps %.1f..%.1f k, wall %.1f..%.1f h, cores %s" % (
            n, min(cells) / 1e3, max(cells) / 1e3, min(steps) / 1e3, max(steps) / 1e3, min(wall), max(wall),
            sorted(cores)))
    ms, ch = [], 0.0
    for n in NS:
        for v in VS:
            c = case(n, v)
            st = json.load(open(os.path.join(c, "state.json")))["stages"]
            ms.append(st["mesh"]["wall_s"])
            ch += summary(c)["solver_wall_s"] * resolved(c)["run"]["n_procs"] / 3600
    print("  mesh stage %.0f..%.0f s; solver core-hours %.0f" % (min(ms), max(ms), ch))
    for root in (REF1, AUX):
        cs = [os.path.join(RUNS, root, c) for c in sorted(os.listdir(os.path.join(RUNS, root)))
              if os.path.exists(os.path.join(RUNS, root, c, "results", "summary.json"))]
        print("  %s: %d cases, solver core-hours %.0f" % (root, len(cs), sum(
            summary(c)["solver_wall_s"] * resolved(c)["run"]["n_procs"] / 3600 for c in cs)))


def m0_cells(c):
    """Cells of the first sub-mesh M0."""
    return json.load(open(os.path.join(c, "mesh", "mesh_report.json")))[0]["hexahedra"]


def studies_ref():
    """Verification and sensitivity studies on the reference mesh (15_studies_ref) against their references."""
    print("== studies (%s) against their references" % AUX)
    aux = lambda name: os.path.join(RUNS, AUX, name)  # noqa: E731
    pairs = [("co05_N2_v1.0", case(2, "1.0")), ("co05_N2_v4.0", case(2, "4.0")), ("rep_N2_v1.0", case(2, "1.0"))]
    pairs += [("ax%s_N%d_v1.0" % (h, n), case(n, "1.0")) for n in (2, 10) for h in ("02", "005")]
    pairs += [("orif_C%d" % c, case(2, "1.0")) for c in (6, 10, 12)]
    pairs += [("orifN10_C%d" % c, case(10, "1.0")) for c in (6, 10, 12)]
    pairs += [("p20_N2_C12", aux("p20_N2_v1.0")), ("p20_N10_C12", aux("p20_N10_v1.0"))]
    for name, ref in pairs:
        if not all(os.path.exists(os.path.join(c, "results", "summary.json")) for c in (aux(name), ref)):
            print("  %-14s not available" % name)
            continue
        a, b = summary(aux(name)), summary(ref)
        print("  %-14s F %.5g vs %.5g (%+.2f%%), Cd %.4f vs %.4f (%+.2f%%), Cd_o %+.2f%%, E %+.2f%%, steps %d vs %d,"
              " cells %.1fk vs %.1fk" % (
                  name, a["F_damp_mean_N"], b["F_damp_mean_N"], 100 * (a["F_damp_mean_N"] / b["F_damp_mean_N"] - 1),
                  a["Cd_effective"], b["Cd_effective"], 100 * (a["Cd_effective"] / b["Cd_effective"] - 1),
                  100 * (cd_orifice(aux(name)) / cd_orifice(ref) - 1),
                  100 * (a["E_diss_J"] / b["E_diss_J"] - 1), a["solver_steps"], b["solver_steps"],
                  m0_cells(aux(name)) / 1e3, m0_cells(ref) / 1e3))
    # relative damping at 20 bar against 1 bar on the reference level and on the finest in-plane level
    for n, fine1 in ((2, aux("orif_C12")), (10, aux("orifN10_C12"))):
        f20, c20 = aux("p20_N%d_C12" % n), aux("p20_N%d_v1.0" % n)
        if not all(os.path.exists(os.path.join(c, "results", "summary.json")) for c in (f20, fine1, c20)):
            continue
        c8 = summary(c20)["F_damp_mean_N"] / 20 / summary(case(n, "1.0"))["F_damp_mean_N"]
        c12 = summary(f20)["F_damp_mean_N"] / 20 / summary(fine1)["F_damp_mean_N"]
        print("  N%-2d F(20 bar)/(20 F(1 bar)) - 1: C8 %+.1f%%, C12 %+.1f%%" % (n, 100 * (c8 - 1), 100 * (c12 - 1)))


def axial_effect():
    """Reference matrix against the matrix with 1 mm axial cells on the same cross-sections, case by case."""
    print("== axial refinement: %s against %s (same cross-sections)" % (MAT, REF1))
    rows = []
    for n in NS:
        for v in VS:
            a, b = case(n, v), case(n, v, REF1)
            if not all(os.path.exists(os.path.join(c, "results", "summary.json")) for c in (a, b)):
                continue
            sa, sb = summary(a), summary(b)
            ch = [summary(c)["solver_wall_s"] * resolved(c)["run"]["n_procs"] / 3600 for c in (a, b)]
            r = dict(n=n, v=v, dF=sa["F_damp_mean_N"] / sb["F_damp_mean_N"] - 1,
                     dCd=sa["Cd_effective"] / sb["Cd_effective"] - 1, dCdo=cd_orifice(a) / cd_orifice(b) - 1,
                     dE=sa["E_diss_J"] / sb["E_diss_J"] - 1, cells=m0_cells(a) / m0_cells(b),
                     steps=sa["solver_steps"] / sb["solver_steps"], coreh=ch[0] / ch[1], ch=ch)
            rows.append(r)
            print("   N%-2d v%-5s F %+.1f%% Cd %+.1f%% Cd_o %+.1f%% E %+.1f%% cells x%.2f steps x%.2f core-h x%.2f"
                  " (%.1f vs %.1f)" % (n, v, 100 * r["dF"], 100 * r["dCd"], 100 * r["dCdo"], 100 * r["dE"],
                                       r["cells"], r["steps"], r["coreh"], ch[0], ch[1]))
    if not rows:
        return
    for v in VS:
        sel = [r for r in rows if r["v"] == v]
        if sel:
            print("  v %-5s dF %+.1f..%+.1f%%, dCd %+.1f..%+.1f%%" % (
                v, 100 * min(r["dF"] for r in sel), 100 * max(r["dF"] for r in sel),
                100 * min(r["dCd"] for r in sel), 100 * max(r["dCd"] for r in sel)))
    print("  all: dF %+.1f..%+.1f%%, dE %+.1f..%+.1f%%, cells x%.2f..%.2f, core-h x%.2f..%.2f, total %.0f vs %.0f"
          % (100 * min(r["dF"] for r in rows), 100 * max(r["dF"] for r in rows),
             100 * min(r["dE"] for r in rows), 100 * max(r["dE"] for r in rows),
             min(r["cells"] for r in rows), max(r["cells"] for r in rows),
             min(r["coreh"] for r in rows), max(r["coreh"] for r in rows),
             sum(r["ch"][0] for r in rows), sum(r["ch"][1] for r in rows)))


def cd_isentropic(c):
    """Orifice coefficient with the mass flux of an isentropic nozzle at the chamber pressures (0.4-0.6 T)."""
    ts, cfg = read_ts(c), resolved(c)
    T = cfg["_derived"]["stroke_time_s"]
    sel = (ts["time_s"] > 0.4 * T) & (ts["time_s"] < 0.6 * T)
    p1, p2, T1, md = (float(np.mean(ts[k][sel])) for k in ("p1", "p2", "T1", "mdot"))
    a_o = cfg["geometry"]["orifice"]["count"] * math.pi * (cfg["geometry"]["orifice"]["radius_mm"] * 1e-3) ** 2
    r = max(p2 / p1, (2.0 / (GAMMA + 1.0)) ** (GAMMA / (GAMMA - 1.0)))
    g = p1 / math.sqrt(R_GAS * T1) * math.sqrt(2.0 * GAMMA / (GAMMA - 1.0) * (r ** (2.0 / GAMMA) - r ** ((GAMMA + 1.0) / GAMMA)))
    return abs(md) / (a_o * g), p2 / p1


def cd_details():
    """Orifice coefficient C_d,o, model with the coefficient of each case, isentropic-nozzle coefficient."""
    print("== discharge coefficient details")
    rows = []
    for n in NS:
        for v in VS:
            c = case(n, v)
            cfg, s = resolved(c), summary(c)
            d, gas = cfg["_derived"], cfg["gas"]
            T0 = gas["initial_temperature_K"]
            uo = u_ideal(cfg)
            ma = uo / math.sqrt(GAMMA * R_GAS * T0)
            afam = d["areas_mm2"]["piston_front"] / d["areas_mm2"]["rod"]
            rows.append(dict(n=n, v=v, ma=ma, cd=s["Cd_effective"], cdo=cd_orifice(c), ratio=s["damping_to_spring_ratio"],
                             afam=afam, c=c))
    mid = [r for r in rows if 0.02 < r["ma"] < 0.3]
    cdo = [r["cdo"] for r in mid]
    cd = float(np.mean([r["cd"] for r in mid]))
    print("  Cd_o mean %.4f +- %.4f; Cd_o/Cd mean %.3f" % (np.mean(cdo), np.std(cdo),
                                                           np.mean([r["cdo"] / r["cd"] for r in mid])))
    print("  at 1 m/s Cd: " + " ".join("N%d %.4f" % (r["n"], r["cd"]) for r in rows if r["v"] == "1.0"))
    for n in NS:
        own = [100 * (0.5 * GAMMA * r["afam"] * (r["ma"] / r["cd"]) ** 2 / r["ratio"] - 1) for r in mid if r["n"] == n]
        if own:
            print("  N%-2d model with own Cd vs CFD ratio: %+.1f .. %+.1f %%" % (n, min(own), max(own)))
    for r in rows:
        if r["n"] == 2 and r["v"] in ("0.4", "1.0", "2.0", "4.0"):
            ci, pr = cd_isentropic(r["c"])
            print("  N2 v%s isentropic Cd %.3f at p2/p1 %.3f" % (r["v"], ci, pr))
    r4 = [r for r in rows if r["n"] == 2 and r["v"] == "4.0"][0]
    line = 0.5 * GAMMA * r4["afam"] * (r4["ma"] / cd) ** 2
    own = 0.5 * GAMMA * r4["afam"] * (r4["ma"] / r4["cd"]) ** 2
    print("  N2 v4: CFD/line-1 %+.1f%%, own-Cd model/CFD-1 %+.1f%%" % (100 * (r4["ratio"] / line - 1), 100 * (own / r4["ratio"] - 1)))
    by = {(r["n"], r["v"]): r for r in rows}
    print("  pair N8 4 m/s / N2 1 m/s ratio diff %+.1f%%" % (100 * (by[(8, "4.0")]["ratio"] / by[(2, "1.0")]["ratio"] - 1)))
    nd = {k: r["ratio"] / (0.5 * GAMMA * r["afam"] * (r["ma"] / cd) ** 2) for k, r in by.items()}
    print("  ratio / Eq.(ratio), against N2 1 m/s: N4 2 m/s %+.1f%%, N8 4 m/s %+.1f%%" % (
        100 * (nd[(4, "2.0")] / nd[(2, "1.0")] - 1), 100 * (nd[(8, "4.0")] / nd[(2, "1.0")] - 1)))
    low = [100 * (r["ratio"] / (0.5 * GAMMA * r["afam"] * (r["ma"] / cd) ** 2) - 1) for r in rows if r["ma"] < 0.02]
    print("  Ma<0.02: CFD/line-1 %+.0f .. %+.0f %%" % (min(low), max(low)))
    for n in (2, 6, 10):
        p = os.path.join(RUNS, AUX, "p20_N%d_v1.0" % n)
        s, cfg = summary(p), resolved(p)
        d = cfg["_derived"]
        ma = u_ideal(cfg) / math.sqrt(GAMMA * R_GAS * cfg["gas"]["initial_temperature_K"])
        line = 0.5 * GAMMA * d["areas_mm2"]["piston_front"] / d["areas_mm2"]["rod"] * (ma / cd) ** 2
        print("  p20 N%-2d CFD/line-1 %+.1f%%, Cd %.4f" % (n, 100 * (s["damping_to_spring_ratio"] / line - 1), s["Cd_effective"]))


SECTIONS = {"fields": fields_2ms, "energy": energy, "conservation": conservation, "segments": segments_study,
            "gas_spring": gas_spring, "cd_ratio": cd_and_ratio, "cd_details": cd_details, "cost": cost,
            "studies": studies_ref, "axial": axial_effect}


def main(argv=None):
    global RUNS, MAT, AUX
    ap = argparse.ArgumentParser(description="Numbers quoted in the text of Section 3 of the paper.")
    ap.add_argument("--runs", default="runs", help="directory holding the study directories (default: runs)")
    ap.add_argument("--sections", default="", help="comma-separated subset of %s (default: all)"
                    % ",".join(SECTIONS))
    ap.add_argument("matrix", nargs="?", default="matrixC8ax", help="matrix study (default: matrixC8ax)")
    ap.add_argument("aux", nargs="?", default="studiesC8ax", help="sensitivity studies (default: studiesC8ax)")
    args = ap.parse_args(argv)
    sel = [x for x in args.sections.split(",") if x]
    bad = [x for x in sel if x not in SECTIONS]
    if bad:
        ap.error("unknown section(s): %s" % " ".join(bad))
    RUNS, MAT, AUX = args.runs, args.matrix, args.aux
    configure(runs=RUNS)
    print("matrix:", MAT, " aux:", AUX)
    for name, fn in SECTIONS.items():
        if not sel or name in sel:
            fn()


if __name__ == "__main__":
    sys.exit(main())
