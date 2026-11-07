#!/usr/bin/env python3
"""Verification numbers on the reference mesh (matrixC8ax, studiesC8ax, furtherC8ax): discrete orifice area,
in-plane discretisation uncertainty, mesh-study tables, checks of the reduced model and of the force decomposition,
and the further cases of 17_further_ref.yaml.

usage: verification_numbers.py [--runs DIR] [--tables DIR] [section ...]

  --runs    directory holding the study directories                                        default: runs
  --tables  directory into which the section inplane writes inplane_bias.csv, inplane_cdo.csv and inplane_est.csv
            (read by figdata.py)                                                           default: tables

sections (default: all except levels and uncertainty, which can be run on their own):
  area         discrete orifice area of every mesh level (polygon of n boundary segments per orifice) and the
               coefficients on that area
  kinematic    low-Mach model with the uniform-pressure flow rate Q = |u_p| (A_s - A_m V1/V), the force area
               A_f - A_m V1/V and the gas density V0/V: factor kappa, C_d,o / C_d, and the prediction of the damping
               ratio with an independent orifice coefficient (Lichtarowicz et al.)
  branches     spring part of the two branches: bias of the branch average 1/2 (F_c - F_r) against the force
               minus the uniform-pressure part, F - p_avg A_m
  window       mean damping force over central windows of 90 ... 50 % of the stroke, decay of the start transient
  hyst         thermal-hysteresis loss of the gas spring with isothermal walls (linear theory) against E_d
  gci          sensitivity of the in-plane extrapolation to the representative size and to the area correction
  meshstudy    rows of the mesh-study table (coefficients on the discrete orifice area)
  inplane      in-plane bias and uncertainty of the reference mesh from the levels C = 8, 12, 16; writes the tables
  levels       in-plane levels of every condition (supplementary table)
  further      further cases: preserved orifice area, conservative field transfer, laminar flow on C = 16,
               orifice count at equal total orifice area
  guidance     accuracy and cost of the mesh levels at every condition
  uncertainty  LaTeX rows of the in-plane uncertainty table
"""

import argparse
import cmath
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figdata import GAMMA, R_GAS, cd_long_orifice, cd_orifice, configure, gci, read_ts, resolved, summary  # noqa: E402

RUNS = "runs"               # directory holding the study directories (set in main)
TABLES = "tables"           # directory of the derived CSV tables (set in main)
MAT, AUX = "matrixC8ax", "studiesC8ax"
FURTHER = "furtherC8ax"     # further cases on the reference mesh (17_further_ref.yaml)
FURTHER1 = "furtherC8"      # further cases with 1 mm axial cells (13_further_C8.yaml)
NS, VS = (2, 4, 6, 8, 10), ("0.025", "0.1", "0.4", "1.0", "2.0", "4.0")


def case(n, v, root=MAT):
    """Case directory of the matrix case with N orifices at the mean speed v (string, as in the case name)."""
    return os.path.join(RUNS, root, "N%d_v%s" % (n, v))


def area_ratio(c):
    """A_o,h / (N pi r_o^2) of the cross-section of case c: regular polygon with the mean number of boundary
    segments per orifice, n = orifice-wall faces / (N * axial segments in the piston)."""
    cfg = resolved(c)
    N = cfg["geometry"]["orifice"]["count"]
    ns = max(2, int(round(cfg["geometry"]["piston_length_mm"] / cfg["mesh"]["piston_cell_size_mm"])))
    rep = json.load(open(os.path.join(c, "mesh", "mesh_report.json")))[0]
    if "orifice_area_ratio" in rep:
        return rep["orifice_area_ratio"], rep["faces"]["orificeWall"] / (N * ns)
    n = rep["faces"]["orificeWall"] / (N * ns)
    return n / (2.0 * math.pi) * math.sin(2.0 * math.pi / n), n


def state(c):
    """Ideal orifice velocity, Mach and Reynolds numbers of case c on the discrete orifice area."""
    cfg = resolved(c)
    d, g, gas = cfg["_derived"], cfg["geometry"], cfg["gas"]
    T0, p0 = gas["initial_temperature_K"], gas["initial_pressure_Pa"]
    c0 = math.sqrt(GAMMA * R_GAS * T0)
    mu0 = gas["sutherland_As"] * T0 ** 1.5 / (T0 + gas["sutherland_Ts_K"])
    a, _ = area_ratio(c)
    uo = d["mean_speed_mps"] * d["orifice_area_ratio"] / a          # ideal orifice velocity on the discrete area
    re = p0 / (R_GAS * T0) * uo * 2e-3 * g["orifice"]["radius_mm"] / mu0
    return dict(cfg=cfg, d=d, a=a, uo=uo, ma=uo / c0, re=re, p0=p0, c0=c0)


# ---------------------------------------------------------------------------------------------- area
def sec_area():
    print("== discrete orifice area A_o,h/(N pi r_o^2), n boundary segments per orifice")
    for lab, c in (("medium N2", "matrix/N2_v1.0"), ("medium N10", "matrix/N10_v1.0"),
                   ("C6", AUX + "/orif_C6"), ("C8 (ref) N2", MAT + "/N2_v1.0"), ("C8 N10", MAT + "/N10_v1.0"),
                   ("C10", AUX + "/orif_C10"), ("C12", AUX + "/orif_C12"), ("C12 N10", AUX + "/orifN10_C12"),
                   ("C16", AUX + "/orif_C16")):
        a, n = area_ratio(os.path.join(RUNS, c))
        print("  %-12s n %5.2f  a %.4f  (F ~ 1/a^2: %+.1f%%)" % (lab, n, a, 100 * (a ** -2 - 1)))
    am = area_ratio(os.path.join(RUNS, "matrixC8/N2_v1.0"))[0]
    a1 = area_ratio(os.path.join(RUNS, "matrix/N2_v1.0"))[0]
    fm, f8 = summary(os.path.join(RUNS, "matrix/N2_v1.0"))["F_damp_mean_N"], \
        summary(os.path.join(RUNS, "matrixC8/N2_v1.0"))["F_damp_mean_N"]
    print("  medium vs C8 (1 mm axial), N2 1 m/s: F %+.1f%%, area alone (a8/am)^2-1 %+.1f%%" % (
        100 * (fm / f8 - 1), 100 * ((am / a1) ** 2 - 1)))
    for n in (6, 10):
        fm = summary(os.path.join(RUNS, "matrix/N%d_v1.0" % n))["F_damp_mean_N"]
        f8 = summary(os.path.join(RUNS, "matrixC8/N%d_v1.0" % n))["F_damp_mean_N"]
        a1 = area_ratio(os.path.join(RUNS, "matrix/N%d_v1.0" % n))[0]
        print("  medium vs C8 (1 mm axial), N%d 1 m/s: F %+.1f%%, area alone %+.1f%%" % (
            n, 100 * (fm / f8 - 1), 100 * ((am / a1) ** 2 - 1)))
    f16 = summary(os.path.join(RUNS, AUX, "orif_C16"))["F_damp_mean_N"]
    f8 = summary(case(2, "1.0"))["F_damp_mean_N"]
    a16 = area_ratio(os.path.join(RUNS, AUX, "orif_C16"))[0]
    print("  C8 -> C16 (0.1 mm axial), N2 1 m/s: F %+.1f%%, area alone %+.1f%%" % (
        100 * (f16 / f8 - 1), 100 * ((am / a16) ** 2 - 1)))
    cds = [summary(case(n, v))["Cd_effective"] for n in NS for v in VS if 0.02 < state(case(n, v))["ma"] < 0.3]
    print("  matrix (16 cases): C_d nominal %.4f, on the discrete area %.4f" % (np.mean(cds), np.mean(cds) / am))


# ----------------------------------------------------------------------------------------- kinematic
def kappa(cfg, a=1.0, npts=200):
    """Stroke average (central 90 %) of (V0/V) q^2 (1 - A_m V1 / (A_f V)), q = (A_s - A_m V1/V)/(A_s - A_m)."""
    d, g = cfg["_derived"], cfg["geometry"]
    A = d["areas_mm2"]
    As, Am = A["cylinder"], A["rod"]
    Ao = A["orifices"] * a
    Af = As - Ao
    L, lp, S, x0 = g["cylinder_length_mm"], g["piston_length_mm"], cfg["motion"]["stroke_mm"], d["x0_mm"]
    V = lambda x: As * x + Ao * lp + (As - Am) * (L - x - lp)          # noqa: E731
    V0 = V(x0)
    xs = np.linspace(x0 - 0.95 * S, x0 - 0.05 * S, npts)
    k = []
    for x in xs:
        v1 = As * x + 0.5 * Ao * lp          # chamber 1 plus half of the orifice volume
        q = (As - Am * v1 / V(x)) / (As - Am)
        k.append(V0 / V(x) * q * q * (1.0 - Am * v1 / (Af * V(x))))
    return float(np.mean(k)), Af / Am


def sec_kinematic():
    print("== low-Mach model with the uniform-pressure flow rate (kinematic factor kappa)")
    rows = []
    for n in NS:
        for v in VS:
            c = case(n, v)
            st, s = state(c), summary(c)
            k, afam = kappa(st["cfg"], st["a"])
            cdo = cd_orifice(c) / st["a"]
            cd = s["Cd_effective"] / st["a"]
            ratio = s["damping_to_spring_ratio"]
            re_l = st["re"] / cdo if cdo == cdo else float("nan")
            # long-orifice coefficient at the Reynolds number of the ideal jet velocity, Re_L = Re_o / C_d,o
            ld = st["cfg"]["geometry"]["piston_length_mm"] / (2 * st["cfg"]["geometry"]["orifice"]["radius_mm"])
            cl = cd_long_orifice(st["re"], ld)
            line = lambda cc: 0.5 * GAMMA * afam * k * (st["ma"] / cc) ** 2      # noqa: E731
            rows.append(dict(n=n, v=v, ma=st["ma"], re=st["re"], k=k, cdo=cdo, cd=cd, ratio=ratio, cl=cl,
                             own=line(cdo) / ratio - 1, corr=line(cl) / ratio - 1, re_l=re_l, afam=afam))
    print("  kappa = %.3f .. %.3f, sqrt(kappa) = %.3f .. %.3f" % (
        min(r["k"] for r in rows), max(r["k"] for r in rows), min(r["k"] for r in rows) ** 0.5,
        max(r["k"] for r in rows) ** 0.5))
    mid = [r for r in rows if 0.02 < r["ma"] < 0.3]
    print("  16 cases: C_d,o/C_d %.3f +- %.3f, C_d,o %.4f +- %.4f, C_d %.4f (discrete area)" % (
        np.mean([r["cdo"] / r["cd"] for r in mid]), np.std([r["cdo"] / r["cd"] for r in mid]),
        np.mean([r["cdo"] for r in mid]), np.std([r["cdo"] for r in mid]), np.mean([r["cd"] for r in mid])))
    print("  model with own C_d,o and kappa vs CFD: %+.1f .. %+.1f %% (mean %+.1f)" % (
        100 * min(r["own"] for r in mid), 100 * max(r["own"] for r in mid), 100 * np.mean([r["own"] for r in mid])))
    print("  model with Lichtarowicz C_d,o and kappa vs CFD: %+.1f .. %+.1f %% (mean %+.1f)" % (
        100 * min(r["corr"] for r in mid), 100 * max(r["corr"] for r in mid),
        100 * np.mean([r["corr"] for r in mid])))
    print("  equivalent C_d of Eq.(ratio) from the correlation, C_d,L/sqrt(kappa): %.3f .. %.3f" % (
        min(r["cl"] / r["k"] ** 0.5 for r in mid), max(r["cl"] / r["k"] ** 0.5 for r in mid)))
    for r in sorted(rows, key=lambda r: r["ma"]):
        print("   N%-2d v%-5s Ma %.4f Re %7.0f kappa %.3f Cd %.3f Cdo %.3f CdL %.3f own %+6.1f%% corr %+6.1f%%" % (
            r["n"], r["v"], r["ma"], r["re"], r["k"], r["cd"], r["cdo"], r["cl"], 100 * r["own"], 100 * r["corr"]))
    return rows


# ------------------------------------------------------------------------------------------ branches
def sec_branches():
    print("== branch average 1/2 (F_c - F_r) against F - p_avg A_m (uniform-pressure spring part removed)")
    for n in NS:
        out = []
        for v in VS:
            c = case(n, v)
            cfg, ts, s = resolved(c), read_ts(c), summary(c)
            d = cfg["_derived"]
            T, S, x0 = d["stroke_time_s"], cfg["motion"]["stroke_mm"], d["x0_mm"]
            Am = d["areas_mm2"]["rod"] * 1e-6
            skip = s["switch_skip_s"]
            sw = [g["start"] for g in d["schedule"]["segments"][1:]]
            ok = np.ones(len(ts["time_s"]), bool)
            for t in sw:
                ok &= ~((ts["time_s"] >= t) & (ts["time_s"] <= t + skip))
            xs = np.linspace(x0 - 0.95 * S, x0 - 0.05 * S, 200)
            br = {}
            for lab, sel in (("c", (ts["time_s"] < T) & ok), ("r", (ts["time_s"] > T) & ok)):
                o = np.argsort(ts["piston_x_mm"][sel])
                x = ts["piston_x_mm"][sel][o]
                br[lab] = (np.interp(xs, x, ts["Fx"][sel][o]), np.interp(xs, x, ts["p_avg"][sel][o] * Am))
            fd_branch = np.mean(0.5 * (br["c"][0] - br["r"][0]))
            fd_direct = np.mean(0.5 * ((br["c"][0] - br["c"][1]) - (br["r"][0] - br["r"][1])))
            out.append("v%s %+.1f%%" % (v, 100 * (fd_branch / fd_direct - 1)))
        print("  N%-2d branch/direct-1: %s" % (n, ", ".join(out)))


# -------------------------------------------------------------------------------------------- window
def sec_window():
    print("== mean damping force over central windows (relative to the 90 % window); start transient")
    for v in ("1.0", "2.0", "4.0"):
        for n in NS:
            c = case(n, v)
            cfg, ts, s = resolved(c), read_ts(c), summary(c)
            d = cfg["_derived"]
            T, S, x0 = d["stroke_time_s"], cfg["motion"]["stroke_mm"], d["x0_mm"]
            skip = s["switch_skip_s"]
            sw = [g["start"] for g in d["schedule"]["segments"][1:]]
            ok = np.ones(len(ts["time_s"]), bool)
            for t in sw:
                ok &= ~((ts["time_s"] >= t) & (ts["time_s"] <= t + skip))
            br = {}
            for lab, sel in (("c", (ts["time_s"] < T) & ok), ("r", (ts["time_s"] > T) & ok)):
                o = np.argsort(ts["piston_x_mm"][sel])
                br[lab] = (ts["piston_x_mm"][sel][o], ts["Fx"][sel][o])
            res = []
            for f in (0.9, 0.8, 0.7, 0.6, 0.5):
                xs = np.linspace(x0 - S + 0.5 * (1 - f) * S, x0 - 0.5 * (1 - f) * S, 200)
                fd = 0.5 * (np.interp(xs, *br["c"]) - np.interp(xs, *br["r"]))
                res.append(np.mean(fd))
            # start transient: deviation of F_c from a quadratic fit over the central 60 %, near the start (x0)
            xc, fc = br["c"]
            m = (xc > x0 - 0.8 * S) & (xc < x0 - 0.2 * S)
            pf = np.polyfit(xc[m], fc[m], 2)
            dev = np.abs(fc - np.polyval(pf, xc)) / s["F_damp_mean_N"]
            near = xc > x0 - 0.5 * S
            big = xc[near & (dev > 0.05)]
            L5 = x0 - big.min() if len(big) else 0.0
            # reversal transient on the rebound branch near x0 - S
            xr, fr = br["r"]
            m = (xr > x0 - 0.8 * S) & (xr < x0 - 0.2 * S)
            pr = np.polyfit(xr[m], fr[m], 2)
            devr = np.abs(fr - np.polyval(pr, xr)) / s["F_damp_mean_N"]
            bigr = xr[(xr < x0 - 0.5 * S) & (devr > 0.05)]
            R5 = bigr.max() - (x0 - S) if len(bigr) else 0.0
            print("  N%-2d v%s F90 %.4g | 80 %+.2f%% 70 %+.2f%% 60 %+.2f%% 50 %+.2f%% | >5%% of F_damp up to %.1f mm "
                  "after start, %.1f mm after reversal" % (n, v, res[0], *[100 * (r / res[0] - 1) for r in res[1:]],
                                                            L5, R5))


# ---------------------------------------------------------------------------------------------- hyst
def sec_hyst():
    """Linear theory of a gas spring with isothermal walls (thermal penetration depth delta = sqrt(2 alpha/omega)):
    V'/V0 = -(p'/(gamma p0)) [1 + (gamma-1) f], f = tanh(z)/z, z = (1+i) r_h/delta (plane layer of the hydraulic
    radius r_h = V/A_w); the loss per cycle is pi gamma p0 V_a^2/V0 |Im 1/(1+(gamma-1) f)| with the volume
    amplitude V_a = A_m S/2 and omega = pi/T for a stroke of duration T."""
    print("== thermal hysteresis with isothermal walls (estimate) against the orifice dissipation E_d")
    k_N2 = 0.0259          # W/(m K), N2 at 300 K
    rows = []
    for n in NS:
        for v in VS:
            c = case(n, v)
            cfg, s = resolved(c), summary(c)
            d, g = cfg["_derived"], cfg["geometry"]
            A = {k: x * 1e-6 for k, x in d["areas_mm2"].items()}
            L, lp, S = g["cylinder_length_mm"] * 1e-3, g["piston_length_mm"] * 1e-3, cfg["motion"]["stroke_mm"] * 1e-3
            x = (d["x0_mm"] - 0.5 * cfg["motion"]["stroke_mm"]) * 1e-3     # mid-stroke
            rs, rm, ro = g["cylinder_radius_mm"] * 1e-3, g["rod_radius_mm"] * 1e-3, g["orifice"]["radius_mm"] * 1e-3
            V = A["cylinder"] * x + A["orifices"] * lp + A["annulus"] * (L - x - lp)
            Aw = (2 * math.pi * rs * L + A["cylinder"] + A["annulus"] + A["piston_front"] + A["piston_back"]
                  + 2 * math.pi * rm * (L - x - lp) + n * 2 * math.pi * ro * lp)
            rh = V / Aw
            p0 = cfg["gas"]["initial_pressure_Pa"]
            rho0 = p0 / (R_GAS * cfg["gas"]["initial_temperature_K"])
            alpha = k_N2 / (rho0 * cfg["gas"]["cp_J_per_kgK"])
            om = math.pi / d["stroke_time_s"]
            delta = math.sqrt(2 * alpha / om)
            z = (1 + 1j) * rh / delta
            f = cmath.tanh(z) / z
            G = 1.0 / (1.0 + (GAMMA - 1.0) * f)
            Va = A["rod"] * S / 2
            E = math.pi * GAMMA * p0 * Va ** 2 / V * abs(G.imag)
            rows.append((n, v, E, s["E_diss_J"], rh / delta))
    for v in VS:
        r = [x for x in rows if x[1] == v]
        print("  v%-5s r_h/delta %.2f  E_hyst %.2e J  E_d %s" % (v, r[0][4], r[0][2], " ".join(
            "N%d %.1e" % (x[0], x[3]) for x in r)))
    dom = [x for x in rows if x[3] > 10 * x[2]]
    print("  E_d > 10 E_hyst in %d of 30 cases: %s" % (len(dom), ", ".join("N%d v%s" % x[:2] for x in dom)))
    dom = [x for x in rows if x[3] > x[2]]
    print("  E_d > E_hyst in %d of 30 cases" % len(dom))


# ----------------------------------------------------------------------------------------------- gci
def sec_gci():
    print("== in-plane extrapolation, N = 2, 1 m/s: sensitivity")
    names = [("C6", AUX + "/orif_C6", 6), ("C8", MAT + "/N2_v1.0", 8), ("C10", AUX + "/orif_C10", 10),
             ("C12", AUX + "/orif_C12", 12), ("C16", AUX + "/orif_C16", 16)]
    F, Cdo, CdoA, nq, a, ns = [], [], [], [], [], []
    for lab, c, C in names:
        p = os.path.join(RUNS, c)
        s = summary(p)
        aa, nseg = area_ratio(p)
        F.append(s["F_damp_mean_N"])
        Cdo.append(cd_orifice(p))
        CdoA.append(cd_orifice(p) / aa)
        a.append(aa)
        ns.append(nseg)
        print("  %-4s n_seg %4.1f a %.4f F %.4f Cd_o %.4f Cd_o(area) %.4f" % (lab, nseg, aa, F[-1], Cdo[-1], CdoA[-1]))
    Cs = np.array([6, 8, 10, 12, 16], float)
    nqs = np.array([48, 132, 161, 224, 429], float)
    for lab, h in (("h ~ 1/C", 1 / Cs), ("h ~ n_q^-1/2", nqs ** -0.5)):
        for qlab, f in (("F", F), ("Cd_o nominal", Cdo), ("Cd_o discrete area", CdoA)):
            for tri in ((4, 3, 2), (4, 3, 1), (3, 2, 1)):
                try:
                    p, fe, g = gci([f[i] for i in tri], [h[i] for i in tri])
                    print("  %-12s %-18s C%d/%d/%d: p %.2f ext %.4f GCI %.1f%% C8-ext %+.1f%%" % (
                        lab, qlab, *[int(Cs[i]) for i in tri], p, fe, 100 * g, 100 * (f[1] / fe - 1)))
                except ValueError as e:
                    print("  %-12s %-18s C%d/%d/%d: %s" % (lab, qlab, *[int(Cs[i]) for i in tri], e))
    # least squares f = f0 + c h^p over the four finest levels (Eca and Hoekstra 2014, basic form)
    for lab, h in (("h ~ 1/C", 1 / Cs), ("h ~ n_q^-1/2", nqs ** -0.5), ("h ~ 1/n_seg", 1 / np.array(ns))):
        for qlab, f in (("F", F), ("Cd_o discrete area", CdoA)):
            hh, ff = h[1:], np.array(f[1:])
            best = None
            for p in np.linspace(0.5, 4.0, 701):
                X = np.vstack([np.ones_like(hh), hh ** p]).T
                coef, res, *_ = np.linalg.lstsq(X, ff, rcond=None)
                r = float(np.sum((X @ coef - ff) ** 2))
                if best is None or r < best[0]:
                    best = (r, p, coef[0])
            sd = math.sqrt(best[0] / (len(ff) - 3))
            print("  LSQ %-12s %-18s C8..C16: p %.2f ext %.4f (s %.4f) C8-ext %+.1f%%" % (
                lab, qlab, best[1], best[2], sd, 100 * (f[1] / best[2] - 1)))
            # five levels C6..C16
            hh, ff = h, np.array(f)
            best = None
            for p in np.linspace(0.3, 4.0, 741):
                X = np.vstack([np.ones_like(hh), hh ** p]).T
                coef, *_ = np.linalg.lstsq(X, ff, rcond=None)
                r = float(np.sum((X @ coef - ff) ** 2))
                if best is None or r < best[0]:
                    best = (r, p, coef[0])
            print("  LSQ %-12s %-18s C6..C16: p %.2f ext %.4f (s %.4f)" % (
                lab, qlab, best[1], best[2], math.sqrt(best[0] / (len(ff) - 3))))



# ----------------------------------------------------------------------------------------- meshstudy
def sec_meshstudy():
    """Rows of the mesh-study table with the discrete orifice area: a, C_d and C_d,o on a A_o."""
    print("== mesh-study table (coefficients on the discrete orifice area)")
    rows = [("Coarse", "studies/mesh_coarse"), ("Medium", "studies/mesh_medium"), ("Fine", "studies/mesh_fine"),
            ("Extra fine", "studies/mesh_xfine"), ("Axial fine", "studies/mesh_axialFine"),
            ("Section fine", "studies/mesh_sectionFine"),
            ("C6 1mm", "studies/orif_C6"), ("C8 1mm", "studies/orif_C8"), ("C10 1mm", FURTHER1 + "/orif_C10"),
            ("C12 1mm", "studies/orif_C12"),
            ("C6", AUX + "/orif_C6"), ("C8 ref", MAT + "/N2_v1.0"), ("C10", AUX + "/orif_C10"),
            ("C12", AUX + "/orif_C12"), ("C16", AUX + "/orif_C16"),
            ("ax0.2", AUX + "/ax02_N2_v1.0"), ("ax0.05", AUX + "/ax005_N2_v1.0"), ("Co0.5", AUX + "/co05_N2_v1.0"),
            ("N10 C6", AUX + "/orifN10_C6"), ("N10 C8 ref", MAT + "/N10_v1.0"), ("N10 C10", AUX + "/orifN10_C10"),
            ("N10 C12", AUX + "/orifN10_C12"), ("N10 ax0.2", AUX + "/ax02_N10_v1.0"),
            ("N10 ax0.05", AUX + "/ax005_N10_v1.0")]
    ref = {2: summary(case(2, "1.0"))["F_damp_mean_N"], 10: summary(case(10, "1.0"))["F_damp_mean_N"]}
    vals = {}
    for lab, c in rows:
        p = os.path.join(RUNS, c)
        s, cfg = summary(p), resolved(p)
        a, n = area_ratio(p)
        rep = json.load(open(os.path.join(p, "mesh", "mesh_report.json")))[0]
        N = s["orifice_count"]
        hf = cfg["mesh"].get("face_cell_size_mm") or cfg["mesh"]["piston_cell_size_mm"]
        v = dict(a=a, n=n, nq=rep["faces"]["orificeInlet"] // N, cells=rep["hexahedra"] / 1e3,
                 steps=s["solver_steps"] / 1e3, F=s["F_damp_mean_N"], cd=s["Cd_effective"] / a,
                 cdo=cd_orifice(p) / a, hf=hf, d=100 * (s["F_damp_mean_N"] / ref[N] - 1), N=N)
        vals[lab] = v
        print("  %-12s & $%.2g$ & $%d$ & $%.3f$ & $%.1f$ & $%.1f$ & $%.4g$ & $%.3f$ & $%.3f$ & $%+.1f$ \\\\" % (
            lab, hf, v["nq"], a, v["cells"], v["steps"], v["F"], v["cd"], v["cdo"], v["d"]))
    for N, labs in ((2, ("C16", "C12", "C10")), (10, ("N10 C12", "N10 C10", "N10 C8 ref"))):
        Cs = [16, 12, 10] if N == 2 else [12, 10, 8]
        for key in ("F", "cd", "cdo"):
            p, fe, g = gci([vals[l][key] for l in labs], [1.0 / C for C in Cs])
            print("  N%d %s: p %.2f ext %.4f GCI %.1f%% C8-ext %+.1f%%" % (
                N, key, p, fe, 100 * g, 100 * (vals["C8 ref" if N == 2 else "N10 C8 ref"][key] / fe - 1)))



# ------------------------------------------------------------------------------- in-plane uncertainty
TRIPLES = [  # (label, C = 8, 12, 16 case directories, relative to the runs directory)
    ("N2 0.1", MAT + "/N2_v0.1", AUX + "/orif_C12_v0.1", FURTHER + "/orif_C16_v0.1"),
    ("N2 0.4", MAT + "/N2_v0.4", FURTHER + "/orif_C12_v0.4", FURTHER + "/orif_C16_v0.4"),
    ("N2 1", MAT + "/N2_v1.0", AUX + "/orif_C12", AUX + "/orif_C16"),
    ("N2 4", MAT + "/N2_v4.0", AUX + "/orif_C12_v4.0", FURTHER + "/orif_C16_v4.0"),
    ("N10 1", MAT + "/N10_v1.0", AUX + "/orifN10_C10", AUX + "/orifN10_C12"),   # levels C = 8, 10, 12
    ("N10 2", MAT + "/N10_v2.0", FURTHER + "/orifN10_C12_v2.0", FURTHER + "/orifN10_C16_v2.0"),
    ("N10 4", MAT + "/N10_v4.0", AUX + "/orifN10_C12_v4.0", FURTHER + "/orifN10_C16_v4.0"),
    ("N2 1 20bar", AUX + "/p20_N2_v1.0", AUX + "/p20_N2_C12", FURTHER + "/p20_N2_C16"),
    ("N10 1 20bar", AUX + "/p20_N10_v1.0", AUX + "/p20_N10_C12", FURTHER + "/p20_N10_C16"),
]


def done(c):
    """True if the case c (relative to the runs directory) has a case summary."""
    return c is not None and os.path.exists(os.path.join(RUNS, c, "results", "summary.json"))


EXT = {}     # label -> extrapolation of the damping force (filled by sec_inplane)
MEDIUM = {"N2 0.1": "matrix/N2_v0.1", "N2 0.4": "matrix/N2_v0.4", "N2 1": "matrix/N2_v1.0", "N2 4": "matrix/N2_v4.0",
          "N10 1": "matrix/N10_v1.0", "N10 2": "matrix/N10_v2.0", "N10 4": "matrix/N10_v4.0",
          "N2 1 20bar": "studies/p20_N2_v1.0", "N10 1 20bar": "studies/p20_N10_v1.0"}


def extrap(f, h):
    """Extrapolation of three levels (f, h ordered fine -> coarse) to h = 0: Richardson with the apparent order
    when it is at least one; without an apparent order (no power law fits the three levels) or with an order below
    one, for which the power law extrapolates without bound, the first-order least-squares expansion f0 + c h, one
    of the expansions of Eca and Hoekstra (2014). Returns (f_ext, order or None)."""
    try:
        p, fe, _ = gci(f, h)
        if p >= 1.0:
            return fe, p
    except ValueError:
        pass
    return np.polyfit(h, f, 1)[1], None


def bounded_gci(f, h):
    """Grid convergence index of the finest level (Celik et al. 2008, safety factor 1.25) with the apparent order
    bounded to 1 <= p <= 2 (2: formal order of the schemes; without an apparent order, or below one, first order),
    absolute. f, h ordered fine -> coarse."""
    try:
        p, _, _ = gci(f, h)
    except ValueError:
        p = None
    pb = min(max(p, 1.0), 2.0) if p else 1.0
    return 1.25 * abs(f[1] - f[0]) / ((h[1] / h[0]) ** pb - 1.0)


HCHOICES = ("1/C", "nq", "nseg")


def sec_inplane():
    """In-plane bias and uncertainty of the reference mesh at several orifice Reynolds numbers from the levels
    C = 8, 12, 16 (N = 10 at 1 m/s: 8, 10, 12), with three choices of the representative size: h ~ 1/C, h ~ n_q^-1/2
    (quadrilaterals per orifice) and h ~ 1/n_seg (boundary segments per orifice). Per choice: extrapolation (extrap)
    and the bounded GCI of the finest level. Estimate: midpoint of the extrapolations; uncertainty U: the largest
    GCI. Writes inplane_bias.csv (bias of C = 8 against the estimate, for interpolation), inplane_cdo.csv (C_d,o
    estimate and U) and inplane_est.csv (estimates against Eq. (ratio) on the full orifice area) into the tables
    directory."""
    print("== in-plane bias and uncertainty of the reference mesh (C = 8/12/16)")
    out = []
    for lab, c8, c12, c16 in TRIPLES:
        if not (done(c8) and done(c12) and done(c16)):
            print("  %-12s incomplete" % lab)
            continue
        cs = [os.path.join(RUNS, c) for c in (c16, c12, c8)]          # fine -> coarse
        st = state(cs[2])
        F = [summary(c)["F_damp_mean_N"] for c in cs]
        cdo = [cd_orifice(c) / area_ratio(c)[0] for c in cs]
        nq = [json.load(open(os.path.join(c, "mesh", "mesh_report.json")))[0]["faces"]["orificeInlet"] /
              resolved(c)["geometry"]["orifice"]["count"] for c in cs]
        nseg = [area_ratio(c)[1] for c in cs]
        Cs = [12, 10, 8] if lab == "N10 1" else [16, 12, 8]
        H = {"1/C": [1.0 / C for C in Cs], "nq": [x ** -0.5 for x in nq], "nseg": [1.0 / x for x in nseg]}
        ests, cests, Us, Ucs, ps, ls1 = [], [], [], [], {}, []
        for hl in HCHOICES:
            h = H[hl]
            fe, p = extrap(F, h)
            ce, pc = extrap(cdo, h)
            ests.append(fe)
            cests.append(ce)
            Us.append(bounded_gci(F, h))
            Ucs.append(bounded_gci(cdo, h))
            ps[hl] = (p, pc)
            ls1 += ["F " + hl] * (p is None) + ["C " + hl] * (pc is None)
        fm, cm = 0.5 * (max(ests) + min(ests)), 0.5 * (max(cests) + min(cests))
        U, Uc = max(Us), max(Ucs)
        # two finest levels equal to 0.1 %: the index carries no information; take the change from the coarsest level
        if abs(F[0] - F[1]) < 1e-3 * abs(F[0]):
            U = max(U, abs(F[0] - F[2]))
        if abs(cdo[0] - cdo[1]) < 1e-3 * abs(cdo[0]):
            Uc = max(Uc, abs(cdo[0] - cdo[2]))
        covers = all(fm - U <= x <= fm + U for x in ests + [F[0]]) and all(cm - Uc <= x <= cm + Uc for x in cests + [cdo[0]])
        # against Eq. (ratio) on the full orifice area (a = 1) with the long-orifice coefficient
        d = st["d"]
        uo1 = d["mean_speed_mps"] * d["orifice_area_ratio"]
        ma1, re1 = uo1 / st["c0"], st["re"] * uo1 / st["uo"]
        k1, afam1 = kappa(st["cfg"], 1.0)
        model = 0.5 * GAMMA * afam1 * k1 * (ma1 / cd_long_orifice(re1)) ** 2
        ratio8 = summary(cs[2])["damping_to_spring_ratio"]
        dev = ratio8 * fm / F[2] / model - 1.0
        print("  %-12s Re %7.0f F16/12/8 %.4g %.4g %.4g | n_q %s n_seg %s" % (
            lab, st["re"], *F, "/".join("%.0f" % x for x in nq), "/".join("%.1f" % x for x in nseg)))
        for hl, fe, ce, u, uc in zip(HCHOICES, ests, cests, Us, Ucs):
            p, pc = ps[hl]
            print("      %-4s p %5s ext %.4g GCI %4.1f%% | Cdo p %5s ext %.4f GCI %4.1f%%" % (
                hl, "%.2f" % p if p else "LS1", fe, 100 * u / F[0], "%.2f" % pc if pc else "LS1", ce, 100 * uc / cdo[0]))
        print("      F_est %.4g +- %.3g (%.1f%%), C8 %+.1f%% (C8/C16 %+.1f%%), C12 %+.1f%%, C16 %+.1f%%; Cdo_est %.4f +- %.4f "
              "(%.1f%%), C8 %+.1f%%; covers %s; model (a=1) dev %+.1f%% (%+.1f..%+.1f)" % (
                  fm, U, 100 * U / fm, 100 * (F[2] / fm - 1), 100 * (F[2] / F[0] - 1), 100 * (F[1] / fm - 1),
                  100 * (F[0] / fm - 1), cm, Uc, 100 * Uc / cm, 100 * (cdo[2] / cm - 1), covers,
                  100 * dev, 100 * (dev - (1 + dev) * U / fm), 100 * (dev + (1 + dev) * U / fm)))
        out.append((st["re"], F[2] / fm - 1, F[2] / (fm + U) - 1, F[2] / (fm - U) - 1, lab))
        EXT[lab] = dict(re=st["re"], fm=fm, fh=0.5 * (max(ests) - min(ests)), lo=min(ests), hi=max(ests), U=U,
                        cm=cm, ch=0.5 * (max(cests) - min(cests)), Uc=Uc, cases=(c8, c12, c16), ls1=ls1, F=F, cdo=cdo,
                        nq=nq, nseg=nseg, ps=ps, ests=ests, cests=cests, Us=Us, Ucs=Ucs, dev=dev, ma1=ma1,
                        p20="20bar" in lab)
    if out:
        os.makedirs(TABLES, exist_ok=True)
        path_bias, path_cdo, path_est = (os.path.join(TABLES, x)
                                         for x in ("inplane_bias.csv", "inplane_cdo.csv", "inplane_est.csv"))
        with open(path_bias, "w") as f:
            f.write("# Re, bias F(C=8)/F_est-1 against the estimate (midpoint of the extrapolations with h~1/C, "
                    "h~nq^-1/2, h~1/nseg), against est+U, against est-U, condition\n")
            for r in sorted(out):
                f.write("%.6g,%.5f,%.5f,%.5f,%s\n" % r)
        print("  wrote %s (%d conditions)" % (path_bias, len(out)))
        with open(path_cdo, "w") as f:
            f.write("# Re, C_d,o estimate (midpoint of the extrapolations), uncertainty U (largest bounded GCI), condition\n")
            for lab, e in sorted(EXT.items(), key=lambda x: x[1]["re"]):
                f.write("%.6g,%.5f,%.5f,%s\n" % (e["re"], e["cm"], e["Uc"], lab))
        with open(path_est, "w") as f:
            f.write("# Ma (full area), Re, estimate over Eq. (ratio) on the full orifice area, uncertainty U of that "
                    "ratio, p0 > 1 bar, condition\n")
            for lab, e in sorted(EXT.items(), key=lambda x: x[1]["re"]):
                f.write("%.6g,%.6g,%.5f,%.5f,%d,%s\n" % (e["ma1"], e["re"], 1 + e["dev"], (1 + e["dev"]) * e["U"] / e["fm"],
                                                        int(e["p20"]), lab))
        print("  wrote %s, %s" % (path_cdo, path_est))


def sec_levels():
    """Supplementary table: the in-plane levels of every condition with four significant digits."""
    if not EXT:
        sec_inplane()
    print("== in-plane levels (supplementary table)")
    for lab, e in sorted(EXT.items(), key=lambda x: (x[1]["p20"], x[1]["re"])):
        Cs = [12, 10, 8] if lab == "N10 1" else [16, 12, 8]
        for i in (2, 1, 0):
            print("  %-12s & $%d$ & $%.0f$ & $%.1f$ & $%.4g$ & $%.4f$ \\\\" % (
                lab, Cs[i], e["nq"][i], e["nseg"][i], e["F"][i], e["cdo"][i]))
        print("      orders F: %s | Cdo: %s" % (
            ", ".join("%s %s" % (h, "%.2f" % e["ps"][h][0] if e["ps"][h][0] else "LS1") for h in HCHOICES),
            ", ".join("%s %s" % (h, "%.2f" % e["ps"][h][1] if e["ps"][h][1] else "LS1") for h in HCHOICES)))


def sec_further():
    """Further cases of 17_further_ref.yaml: orifice discs scaled to the nominal area, conservative field transfer at
    the mesh switches, laminar flow on the level C = 16, orifice counts at equal total orifice area."""
    print("== further cases")
    def F(c):
        return summary(os.path.join(RUNS, c))["F_damp_mean_N"]
    # orifice area preserved
    for lab, a, b in (("area C8ax N2 v1", FURTHER + "/area_N2_v1.0", MAT + "/N2_v1.0"),
                      ("medium, area preserved", FURTHER + "/medarea_N2_v1.0", FURTHER + "/med_N2_v1.0")):
        if done(a) and done(b):
            pa, pb = os.path.join(RUNS, a), os.path.join(RUNS, b)
            aa, ab = area_ratio(pa)[0], area_ratio(pb)[0]
            print("  %-24s F %+.2f%% (area alone (a_b/a_a)^2-1 %+.2f%%), Cd_o on discrete area %.4f vs %.4f" % (
                lab, 100 * (F(a) / F(b) - 1), 100 * ((ab / aa) ** 2 - 1), cd_orifice(pa) / aa, cd_orifice(pb) / ab))
    for ref in ("matrixC8/N2_v1.0", "studies/orif_C8", "matrix/N2_v1.0"):
        for c in (FURTHER + "/med_N2_v1.0", FURTHER + "/medarea_N2_v1.0"):
            if done(c) and done(ref):
                print("  %s against %s: F %+.2f%%" % (c, ref, 100 * (F(c) / F(ref) - 1)))
    # conservative transfer
    a, b = FURTHER + "/cvw_N2_v4.0", MAT + "/N2_v4.0"
    if done(a) and done(b):
        sa, sb = summary(os.path.join(RUNS, a)), summary(os.path.join(RUNS, b))
        print("  cellVolumeWeight N2 v4: jumps %s (interp %s), drift %.2e (interp %.2e), F %+.3f%%, E %+.3f%%" % (
            ["%.1e" % x for x in sa["mass_jumps_rel"]], ["%.1e" % x for x in sb["mass_jumps_rel"]],
            sa["mass_drift_max_rel"], sb["mass_drift_max_rel"], 100 * (sa["F_damp_mean_N"] / sb["F_damp_mean_N"] - 1),
            100 * (sa["E_diss_J"] / sb["E_diss_J"] - 1)))
    # laminar on C = 16 (C_d,o on the discrete orifice area of C = 16, a = 0.9966)
    a, b = FURTHER + "/lam_N2_C16", AUX + "/orif_C16"
    if done(a) and done(b):
        print("  laminar vs SST on C = 16, N2 1 m/s: F %+.1f%%, Cd_o %.4f vs %.4f" % (
            100 * (F(a) / F(b) - 1), cd_orifice(os.path.join(RUNS, a)) / 0.9966, cd_orifice(os.path.join(RUNS, b)) / 0.9966))
    # equal total orifice area
    for v in ("1.0", "4.0"):
        cs = [(2, FURTHER + "/eqA_N2_v%s" % v), (6, MAT + "/N6_v%s" % v), (10, FURTHER + "/eqA_N10_v%s" % v)]
        cs = [(n, c) for n, c in cs if done(c)]
        if len(cs) < 2 or cs[0][0] == 10 and len(cs) < 2 or not any(n == 6 for n, _ in cs):
            continue
        cs.sort(key=lambda x: x[0] != 6)
        cs = cs[:1] + sorted(cs[1:])
        base = summary(os.path.join(RUNS, cs[0][1]))["damping_to_spring_ratio"]
        for n, c in cs:
            p = os.path.join(RUNS, c)
            st, s = state(p), summary(p)
            k, afam = kappa(st["cfg"], st["a"])
            ld = st["cfg"]["geometry"]["piston_length_mm"] / (2 * st["cfg"]["geometry"]["orifice"]["radius_mm"])
            cl = cd_long_orifice(st["re"], ld)
            model = 0.5 * GAMMA * afam * k * (st["ma"] / cl) ** 2
            print("  equal area v%s N%-2d r_o %.3f Ma %.4f Re %6.0f ratio %.4g (%+.1f%% vs N6) Cd_o %.3f "
                  "model %.4g CFD/model %.3f l/d %.2f" % (
                      v, n, st["cfg"]["geometry"]["orifice"]["radius_mm"], st["ma"], st["re"],
                      s["damping_to_spring_ratio"], 100 * (s["damping_to_spring_ratio"] / base - 1),
                      cd_orifice(p) / st["a"], model, s["damping_to_spring_ratio"] / model,
                      st["cfg"]["geometry"]["piston_length_mm"] / (2 * st["cfg"]["geometry"]["orifice"]["radius_mm"])))


def sec_guidance():
    """Accuracy and cost of the mesh levels: bias of the damping force against the extrapolated estimate at each
    condition, and the work (cells of M0 x solver steps) relative to the reference mesh of the same condition."""
    if not EXT:
        sec_inplane()
    print("== accuracy and cost guidance (bias of F against the extrapolated midpoint; work = cells x steps / ref)")
    def work(c):
        p = os.path.join(RUNS, c)
        return (json.load(open(os.path.join(p, "state.json")))["info"]["mesh"]["meshes"]["M0"]["cells"] *
                summary(p)["solver_steps"])
    for lab, e in sorted(EXT.items(), key=lambda x: x[1]["re"]):
        c8, c12, c16 = e["cases"]
        levels = [("medium", MEDIUM.get(lab)), ("C8", c8), ("C10/12" if lab == "N10 1" else "C12", c12),
                  ("C12" if lab == "N10 1" else "C16", c16)]
        w8 = work(c8)
        row = "  %-12s Re %7.0f ext %.4g+-%.1f%% |" % (lab, e["re"], e["fm"], 100 * e["fh"] / e["fm"])
        for name, c in levels:
            if not done(c):
                row += " %s -- |" % name
                continue
            Fc = summary(os.path.join(RUNS, c))["F_damp_mean_N"]
            row += " %s %+.1f%% (%+.1f..%+.1f) w %.2f |" % (name, 100 * (Fc / e["fm"] - 1), 100 * (Fc / e["hi"] - 1),
                                                         100 * (Fc / e["lo"] - 1), work(c) / w8)
        print(row)


def fix3(x):
    """Three significant digits in fixed notation."""
    d = 2 - int(math.floor(math.log10(abs(x))))
    return "%.*f" % (d, x) if d >= 0 else "%d" % round(x, d)


def sec_uncertainty():
    """LaTeX rows of the uncertainty table: per condition the orifice Reynolds number, the damping force of the
    reference mesh, the estimate (midpoint of the two extrapolations) with its half-width, C_d,o likewise, and the
    bias of the medium mesh and of the levels C = 8, 12, 16 against the estimate; last row the work (cells of M0 x
    solver steps) relative to C = 8, range over the conditions."""
    if not EXT:
        sec_inplane()
    print("== uncertainty table rows")
    def work(c):
        p = os.path.join(RUNS, c)
        return (json.load(open(os.path.join(p, "state.json")))["info"]["mesh"]["meshes"]["M0"]["cells"] *
                summary(p)["solver_steps"])
    ws = {"medium": [], "C12": [], "C16": []}
    for lab, e in sorted(EXT.items(), key=lambda x: (" 20bar" in x[0], x[1]["re"])):
        n, v = lab.split()[:2]
        cond = "$N=%s$, $%s~\\mathrm{m/s}$%s" % (n[1:], "%g" % float(v), ", $20~\\mathrm{bar}$" if "20bar" in lab else "")
        c8, c12, c16 = e["cases"]
        def mark(q):
            return "\\textsuperscript{a}" if any(x.startswith(q) for x in e.get("ls1", [])) else ""
        def sig4(x):
            d = 3 - int(math.floor(math.log10(abs(x))))
            return "%.*f" % (max(d, 0), x)
        cells = [cond, "$%s$" % fix3(e["re"]) if e["re"] < 1e4 else "$%.1f\\times10^{%d}$" % (
            e["re"] / 10 ** int(math.log10(e["re"])), int(math.log10(e["re"]))),
            "$%s$%s" % (sig4(e["fm"]), mark("F")), "$%.1f$" % (100 * e["U"] / e["fm"]),
            "$%.3f$%s" % (e["cm"], mark("C")), "$%.1f$" % (100 * e["Uc"] / e["cm"])]
        levels = (("medium", MEDIUM.get(lab)), ("C8", c8), ("C12", c12), ("C16", c16))
        if lab == "N10 1":                     # levels C = 8, 10, 12: no C = 16
            levels = (("medium", MEDIUM.get(lab)), ("C8", c8), ("C12", c16), ("C16", None))
        for name, c in levels:
            if not done(c):
                cells.append("")
                continue
            cells.append("$%+.0f$" % (100 * (summary(os.path.join(RUNS, c))["F_damp_mean_N"] / e["fm"] - 1)))
            if name in ws:
                ws[name].append(work(c) / work(c8))
        print(" & ".join(cells) + " \\\\")
    print("work relative to C = 8: " + ", ".join("%s %.2f-%.2f" % (k, min(w), max(w)) for k, w in ws.items() if w))


SECTIONS = ("area", "kinematic", "branches", "window", "hyst", "gci", "meshstudy", "inplane", "levels", "further",
            "guidance", "uncertainty")
DEFAULT = ("area", "kinematic", "branches", "window", "hyst", "gci", "meshstudy", "inplane", "further", "guidance")


def main(argv=None):
    global RUNS, TABLES
    ap = argparse.ArgumentParser(description="Verification numbers on the reference mesh.")
    ap.add_argument("--runs", default="runs", help="directory holding the study directories (default: runs)")
    ap.add_argument("--tables", default="tables", help="directory of the derived CSV tables (default: tables)")
    ap.add_argument("sections", nargs="*", metavar="section",
                    help="any of %s (default: %s)" % (" ".join(SECTIONS), " ".join(DEFAULT)))
    args = ap.parse_args(argv)
    bad = [x for x in args.sections if x not in SECTIONS]
    if bad:
        ap.error("unknown section(s): %s" % " ".join(bad))
    RUNS, TABLES = args.runs, args.tables
    configure(runs=RUNS, tables=TABLES)
    for s in args.sections or DEFAULT:
        globals()["sec_" + s]()


if __name__ == "__main__":
    sys.exit(main())
