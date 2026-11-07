#!/usr/bin/env python3
"""Data files of the paper figures and derived case tables, computed from the stored case records.

usage: figdata.py [--runs DIR] [--tables DIR] [--out DIR] [output ...]

  --runs    directory holding the study directories (matrixC8ax, studiesC8ax, ...)        default: runs
  --tables  directory of the derived CSV tables                                            default: tables
            read:    orifice_benchmark.csv, orifice_benchmark_T2.csv (orifice_benchmark.py collect),
                     inplane_est.csv, inplane_cdo.csv (verification_numbers.py inplane)
            written: energy_balance.csv, orifice_convergence.csv, cases_table.csv
  --out     directory of the figure data files, one sub-folder per figure                  default: figures
  output    any of loops ratio cd benchmark mass spring history energy orifice cases       default: all

Every figure data file starts with a provenance line naming the source case directories relative to --runs.
The functions of this module are also used by paper_numbers.py, verification_numbers.py and bench_table.py.
"""

import argparse
import csv
import functools
import glob
import json
import math
import os
import sys

import numpy as np

RUNS = "runs"               # directory holding the study directories (set by configure)
TABLES = "tables"           # directory of the derived CSV tables (set by configure)
FIG = "figures"             # directory of the figure data, one sub-folder per figure (set by configure)
MATRIX = "matrixC8ax"       # reference matrix: orifice-refined cross-section (curvature size 8) with axial cells of
                            # 0.1 mm in the orifices; "matrixC8": same cross-section with 1 mm axial cells,
                            # "matrix": medium mesh
AUX = "studiesC8ax"         # verification and sensitivity studies on the reference mesh
REF1 = "matrixC8"           # matrix with 1 mm axial cells, for the effect of the axial refinement
FURTHER = "furtherC8ax"     # further cases on the reference mesh (17_further_ref.yaml)
FURTHER1 = "furtherC8"      # further cases with 1 mm axial cells (13_further_C8.yaml)
BATCH_LOG = "slurm-*.out"   # log of the batch scheduler, present in the cases solved on the compute node

GAMMA = 1040.0 / (1040.0 - 8314.462618 / 28.0134)
R_GAS = 8314.462618 / 28.0134


def configure(runs=None, tables=None, out=None):
    """Set the directory of the study directories, of the derived tables and of the figure data."""
    global RUNS, TABLES, FIG
    if runs is not None:
        RUNS = runs
    if tables is not None:
        TABLES = tables
    if out is not None:
        FIG = out


@functools.lru_cache(maxsize=None)
def read_ts(case):
    """Columns of results/timeseries.csv as float arrays (non-numeric entries NaN); cached, do not modify."""
    path = os.path.join(case, "results", "timeseries.csv")
    with open(path) as f:
        r = csv.reader(f)
        head = next(r)
        rows = list(r)
    cols = {}
    for j, k in enumerate(head):
        col = [row[j] if j < len(row) else "" for row in rows]
        try:
            cols[k] = np.array(col, dtype=float)
        except ValueError:
            vals = []
            for v in col:
                try:
                    vals.append(float(v))
                except ValueError:
                    vals.append(np.nan)
            cols[k] = np.array(vals)
    return cols


def summary(case):
    """Case summary written by the post stage, results/summary.json."""
    return json.load(open(os.path.join(case, "results", "summary.json")))


def resolved(case):
    """Resolved configuration of the case, config.resolved.json."""
    return json.load(open(os.path.join(case, "config.resolved.json")))


def finished(pattern):
    """Case directories matching pattern (relative to the runs directory) that have a case summary."""
    return sorted(c for c in glob.glob(os.path.join(RUNS, pattern))
                  if os.path.exists(os.path.join(c, "results", "summary.json")))


def write(fig, name, header, rows, sources):
    """Write one figure data file: provenance line, column names, rows."""
    d = os.path.join(FIG, fig)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name), "w") as f:
        f.write("# sources: %s\n" % " ".join(os.path.relpath(s, RUNS) for s in sources))
        f.write(" ".join(header) + "\n")
        for r in rows:
            f.write(" ".join("%.6g" % v for v in r) + "\n")


def thin(x, n=400):
    """Indices of at most n samples spread evenly over x."""
    if len(x) <= n:
        return np.arange(len(x))
    return np.unique(np.linspace(0, len(x) - 1, n).astype(int))


def fig_loops():
    """Force-displacement loops, N = 2, all speeds; force normalised by p0*A_rod."""
    for case in finished(MATRIX + "/N2_v*"):
        cfg = resolved(case)
        d = cfg["_derived"]
        T = d["stroke_time_s"]
        a_rod = d["areas_mm2"]["rod"] * 1e-6
        p0 = cfg["gas"]["initial_pressure_Pa"]
        ts = read_ts(case)
        for lab, sel in (("c", ts["time_s"] <= T), ("r", ts["time_s"] > T)):
            x, F = ts["piston_x_mm"][sel], ts["Fx"][sel] / (p0 * a_rod)
            i = thin(x)
            v = cfg["motion"]["mean_speed_mps"]
            write("forceLoops", "N2_v%g_%s.dat" % (v, lab), ["x", "F"],
                  np.column_stack([x[i], F[i]]), [case])
        # adiabatic gas spring, Eq. (gasSpring), normalised the same way
        g, a = cfg["geometry"], d["areas_mm2"]
        xs = np.linspace(d["x0_mm"] - cfg["motion"]["stroke_mm"], d["x0_mm"], 200)
        V = (a["cylinder"] * xs + a["orifices"] * g["piston_length_mm"]
             + a["annulus"] * (g["cylinder_length_mm"] - xs - g["piston_length_mm"]))
        write("forceLoops", "adiabatic.dat", ["x", "F"],
              np.column_stack([xs, (d["volume_initial_mm3"] / V) ** GAMMA]), [case])


def inplane_bias(re, high_pressure=False):
    """Estimated in-plane bias of the reference mesh in the damping force, F(C=8)/F(converged) - 1, at the orifice
    Reynolds number re: log-linear interpolation of inplane_bias.csv in the tables directory (written by
    verification_numbers.py inplane from the C = 8/12/16 triples) over the conditions at the same pre-charge pressure
    (1 or 20 bar), averaged over conditions at the same Reynolds number, constant beyond their ends; 0 if the table
    does not exist."""
    path = os.path.join(TABLES, "inplane_bias.csv")
    if not os.path.exists(path):
        return 0.0
    t = np.array([[float(x) for x in l.split(",")[:2]] for l in open(path)
                  if l[0].isdigit() and ("20bar" in l) == high_pressure])
    res = np.unique(t[:, 0])                     # conditions at the same Reynolds number: mean bias
    bias = [t[t[:, 0] == r, 1].mean() for r in res]
    return float(np.interp(math.log10(re), np.log10(res), bias))


def fig_ratio():
    """Damping-to-spring ratio against the ideal orifice Mach number (discrete orifice area), all finished cases,
    and its ratio to Eq. (ratio) with the kinematic factor kappa and the long-orifice coefficient of Lichtarowicz
    et al. at the orifice Reynolds number of the case (no fitted constant); est.dat, est20.dat: the in-plane
    estimates of the nine conditions of the uncertainty table over Eq. (ratio) on the full orifice area (a = 1),
    to which the in-plane refinement converges, with their uncertainty U as error bars."""
    rows, src = [], []
    for case in finished(MATRIX + "/N*_v*") + finished(AUX + "/p20_N*_v*") + finished(AUX + "/lam_*"):
        s = summary(case)
        if "damping_to_spring_ratio" not in s:
            continue
        st = case_state(case)
        cfg = st["cfg"]
        k, afam = kappa(cfg, st["a"])
        cl = cd_long_orifice(st["re"])
        model = 0.5 * GAMMA * afam * k * (st["ma"] / cl) ** 2
        dev = s["damping_to_spring_ratio"] / model
        lam = 1.0 if cfg["gas"]["turbulence"] == "laminar" else 0.0
        rows.append([cfg["geometry"]["orifice"]["count"], cfg["_derived"]["mean_speed_mps"],
                     cfg["gas"]["initial_pressure_Pa"] / 1e5, st["ma"], s["damping_to_spring_ratio"],
                     cd_orifice(case) / st["a"], afam, k, lam, st["re"], cl, model, dev, 0.0])
        src.append(case)
    rows.sort()
    hdr = ["N", "v", "p0bar", "Ma", "ratio", "Cdo", "AfAm", "kappa", "laminar", "Re", "CdL", "model", "dev", "err"]
    nan = [[np.nan] * len(hdr)]
    for fig in ("dampingRatio", "dampingDeviation"):
        write(fig, "cases.dat", hdr, rows, src)
        for n in (2, 4, 6, 8, 10):
            write(fig, "N%d.dat" % n, hdr, [r for r in rows if r[0] == n and r[2] == 1.0 and r[8] == 0.0] or nan, src)
        write(fig, "p20.dat", hdr, [r for r in rows if r[2] > 1.5] or nan, src)
        write(fig, "laminar.dat", hdr, [r for r in rows if r[8] == 1.0] or nan, src)
    # Eq. (ratio) with kappa, A_f/A_m of N = 2 and the long-orifice coefficient at 1 bar (Re_o proportional to Ma_o)
    ref = [r for r in rows if r[0] == 2 and r[2] == 1.0 and r[8] == 0.0][0]
    re_per_ma = ref[9] / ref[3]
    ma = np.logspace(-4, 0, 81)
    write("dampingRatio", "analytic.dat", ["Ma", "ratio"],
          np.column_stack([ma, [0.5 * GAMMA * ref[6] * ref[7] * (m / cd_long_orifice(m * re_per_ma)) ** 2
                                for m in ma]]), src)
    # converged estimates: finer in-plane levels (C = 16) of the same conditions
    conv = []
    for case in finished(AUX + "/orif_C16") + finished(FURTHER + "/orif*_C16_v*") + finished(FURTHER + "/p20_N*_C16"):
        s, st = summary(case), case_state(case)
        k, afam = kappa(st["cfg"], st["a"])
        model = 0.5 * GAMMA * afam * k * (st["ma"] / cd_long_orifice(st["re"])) ** 2
        conv.append([st["cfg"]["geometry"]["orifice"]["count"], st["ma"], st["re"],
                     s["damping_to_spring_ratio"] / model])
    write("dampingDeviation", "C16.dat", ["N", "Ma", "Re", "dev"], conv or [[np.nan] * 4], [])
    # in-plane estimates of the nine conditions (verification_numbers.py inplane): ratio to Eq. (ratio) on the full
    # orifice area, +-U
    est = os.path.join(TABLES, "inplane_est.csv")
    if os.path.exists(est):
        t = [l.strip().split(",") for l in open(est) if l[0].isdigit()]
        e1 = [[float(r[0]), float(r[1]), float(r[2]), float(r[3])] for r in t if r[4] == "0"]
        e20 = [[float(r[0]), float(r[1]), float(r[2]), float(r[3])] for r in t if r[4] == "1"]
        write("dampingDeviation", "est.dat", ["Ma", "Re", "dev", "err"], e1 or [[np.nan] * 4], [est])
        write("dampingDeviation", "est20.dat", ["Ma", "Re", "dev", "err"], e20 or [[np.nan] * 4], [est])
    with open(os.path.join(FIG, "dampingRatio", "cd_used.txt"), "w") as f:
        mid = [r for r in rows if 0.02 < r[3] < 0.3 and r[2] == 1.0 and r[8] == 0.0]
        f.write("model: Eq. (ratio) with kappa and C_d,o of Lichtarowicz et al. at Re_o; line for N = 2 "
                "(A_f/A_m %.4f, kappa %.4f); 16 cases 0.02<Ma<0.3: dev %.3f .. %.3f; fitted C_d,o %.4f\n"
                % (ref[6], ref[7], min(r[12] for r in mid), max(r[12] for r in mid),
                   float(np.mean([r[5] for r in mid]))))
    return rows


def fig_cd():
    """Orifice discharge coefficient on the discrete orifice area against the orifice Reynolds number, with the
    long-orifice correlation and the finer in-plane level C = 16."""
    rows, src = [], []
    for case in finished(MATRIX + "/N*_v*") + finished(AUX + "/p20_N*_v*") + finished(AUX + "/lam_*"):
        s = summary(case)
        if "Cd_effective" not in s:
            continue
        st = case_state(case)
        cfg = st["cfg"]
        rows.append([cfg["geometry"]["orifice"]["count"], cfg["_derived"]["mean_speed_mps"],
                     cfg["gas"]["initial_pressure_Pa"] / 1e5, st["re"], st["ma"], cd_orifice(case) / st["a"],
                     s["Cd_effective"] / st["a"], 1.0 if cfg["gas"]["turbulence"] == "laminar" else 0.0])
        src.append(case)
    rows.sort(key=lambda r: r[3])
    hdr = ["N", "v", "p0bar", "Re", "Ma", "Cdo", "Cd", "laminar"]
    nan = [[np.nan] * len(hdr)]
    write("dischargeCoefficient", "cases.dat", hdr, rows, src)
    for n in (2, 4, 6, 8, 10):
        write("dischargeCoefficient", "N%d.dat" % n, hdr,
              [r for r in rows if r[0] == n and r[2] == 1.0 and r[7] == 0.0] or nan, src)
    write("dischargeCoefficient", "p20.dat", hdr, [r for r in rows if r[2] > 1.5] or nan, src)
    write("dischargeCoefficient", "laminar.dat", hdr, [r for r in rows if r[7] == 1.0] or nan, src)
    c16, s16 = [], []
    for case in finished(AUX + "/orif_C16") + finished(FURTHER + "/orif*_C16_v*") + finished(FURTHER + "/p20_N*_C16"):
        st = case_state(case)
        c16.append([st["cfg"]["geometry"]["orifice"]["count"], st["re"], cd_orifice(case) / st["a"]])
        s16.append(case)
    c16.sort(key=lambda r: r[1])
    write("dischargeCoefficient", "C16.dat", ["N", "Re", "Cdo"], c16 or [[np.nan] * 3], s16)
    re = np.logspace(1, 6, 101)
    write("dischargeCoefficient", "lichtarowicz.dat", ["Re", "Cd"],
          np.column_stack([re, [cd_long_orifice(r) for r in re]]), [])
    return rows


def fig_mass():
    """Relative mass and volume error histories of selected cases."""
    for name in ("N2_v0.1", "N2_v1.0", "N2_v4.0", "N10_v0.1", "N10_v1.0", "N10_v4.0"):
        case = os.path.join(RUNS, MATRIX, name)
        if not os.path.exists(os.path.join(case, "results", "timeseries.csv")):
            continue
        cfg = resolved(case)
        T = cfg["_derived"]["stroke_time_s"]
        ts = read_ts(case)
        m = ts["m_total"] / ts["m_total"][0] - 1.0
        r = ts["V_cfd_m3"] / ts["V_analytic_m3"]
        r = r / r[0] - 1.0
        i = thin(ts["time_s"], 600)
        am = np.maximum(np.abs(m), 1e-12)  # log axis: floor at 1e-12
        write("massVolume", name + ".dat", ["tT", "dm", "absdm", "dV"],
              np.column_stack([ts["time_s"][i] / T, m[i], am[i], r[i]]), [case])


def fig_spring():
    """Gas-spring part of the force against the adiabatic model, low-speed cases."""
    for case in finished(MATRIX + "/N*_v0.025") + finished(MATRIX + "/N*_v0.1"):
        cfg = resolved(case)
        d = cfg["_derived"]
        T = d["stroke_time_s"]
        g = cfg["geometry"]
        a = d["areas_mm2"]
        p0 = cfg["gas"]["initial_pressure_Pa"]
        ts = read_ts(case)
        sel = ts["time_s"] <= T
        x = ts["piston_x_mm"][sel]
        V = (a["cylinder"] * x + a["orifices"] * g["piston_length_mm"]
             + a["annulus"] * (g["cylinder_length_mm"] - x - g["piston_length_mm"]))
        F_ad = p0 * a["rod"] * 1e-6 * (d["volume_initial_mm3"] / V) ** GAMMA
        i = thin(x)
        name = os.path.basename(case)
        write("gasSpring", name + ".dat", ["x", "F", "Fad", "err"],
              np.column_stack([x[i], ts["Fx"][sel][i], F_ad[i], ts["Fx"][sel][i] / F_ad[i] - 1.0]),
              [case])


def fig_history():
    """Chamber-averaged pressures and mean orifice Mach number over the cycle, 2 m/s, all orifice counts."""
    for case in finished(MATRIX + "/N*_v2.0"):
        cfg = resolved(case)
        T = cfg["_derived"]["stroke_time_s"]
        ts = read_ts(case)
        i = thin(ts["time_s"], 800)
        name = "N%d.dat" % cfg["geometry"]["orifice"]["count"]
        write("pressureHistory", name, ["tT", "p1", "p2"],
              np.column_stack([ts["time_s"][i] / T, ts["p1"][i] / 1e3, ts["p2"][i] / 1e3]), [case])
        write("machHistory", name, ["tT", "Ma"],
              np.column_stack([ts["time_s"][i] / T, ts["Ma_orifice_mean"][i]]), [case])


def table_energy():
    """Energy balance of one cycle (adiabatic walls): internal-energy rise at the initial volume against E_d."""
    rows = []
    for case in finished(MATRIX + "/N*_v*"):
        cfg, s = resolved(case), summary(case)
        ts = read_ts(case)
        p0 = cfg["gas"]["initial_pressure_Pa"]
        V0 = ts["V_analytic_m3"][0]
        dU = (ts["p_avg"][-1] - p0) * V0 / (GAMMA - 1.0)
        rows.append({"case": os.path.relpath(case, RUNS), "t_end_s": ts["time_s"][-1],
                     "p_avg_end_Pa": ts["p_avg"][-1], "dU_J": dU, "E_diss_J": s.get("E_diss_J")})
    with open(os.path.join(TABLES, "energy_balance.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def cd_orifice(case):
    """Orifice discharge coefficient of a damper case in the usual definition: actual mass flow, density of the
    upstream chamber 1 and the chamber-averaged pressure difference, over 0.4-0.6 T of the compression."""
    ts, cfg = read_ts(case), resolved(case)
    T = cfg["_derived"]["stroke_time_s"]
    sel = (ts["time_s"] > 0.4 * T) & (ts["time_s"] < 0.6 * T)
    p1, p2, r1, md = (float(np.mean(ts[k][sel])) for k in ("p1", "p2", "rho1", "mdot"))
    a_o = cfg["geometry"]["orifice"]["count"] * math.pi * (cfg["geometry"]["orifice"]["radius_mm"] * 1e-3) ** 2
    return abs(md) / (a_o * math.sqrt(2.0 * r1 * (p1 - p2))) if p1 > p2 else float("nan")


def lichtarowicz(re_l, ld=4.0):
    """Discharge coefficient of long orifices, Lichtarowicz et al. (1965); re_l based on the ideal jet velocity."""
    cdu = 0.827 - 0.0085 * ld
    return 1.0 / (1.0 / cdu + 20.0 / re_l * (1.0 + 2.25 * ld)
                  - 0.005 * ld / (1.0 + 7.5 * math.log10(0.00015 * re_l) ** 2))


def area_ratio(case):
    """Discrete orifice area over N pi r_o^2: the regular polygon of the mean number n of boundary segments per
    orifice (orifice-wall faces / (N * axial segments in the piston)), or the ratio recorded by the mesh stage."""
    cfg = resolved(case)
    N = cfg["geometry"]["orifice"]["count"]
    ns = max(2, int(round(cfg["geometry"]["piston_length_mm"] / cfg["mesh"]["piston_cell_size_mm"])))
    rep = json.load(open(os.path.join(case, "mesh", "mesh_report.json")))[0]
    n = rep["faces"]["orificeWall"] / (N * ns)
    if "orifice_area_ratio" in rep:
        return rep["orifice_area_ratio"], n
    return n / (2.0 * math.pi) * math.sin(2.0 * math.pi / n), n


def kappa(cfg, a=1.0, npts=200):
    """Kinematic factor of Eq. (ratio): stroke average (central 90 %) of (V0/V) q^2 (1 - A_m V1/(A_f V)) with
    q = (A_s - A_m V1/V)/(A_s - A_m); returns (kappa, A_f/A_m) for the discrete orifice area a A_o."""
    d, g = cfg["_derived"], cfg["geometry"]
    A = d["areas_mm2"]
    As, Am, Ao = A["cylinder"], A["rod"], A["orifices"] * a
    Af = As - Ao
    L, lp, S, x0 = g["cylinder_length_mm"], g["piston_length_mm"], cfg["motion"]["stroke_mm"], d["x0_mm"]

    def V(x):
        return As * x + Ao * lp + (As - Am) * (L - x - lp)
    V0 = V(x0)
    k = []
    for x in np.linspace(x0 - 0.95 * S, x0 - 0.05 * S, npts):
        v1 = As * x + 0.5 * Ao * lp
        q = (As - Am * v1 / V(x)) / (As - Am)
        k.append(V0 / V(x) * q * q * (1.0 - Am * v1 / (Af * V(x))))
    return float(np.mean(k)), Af / Am


def cd_long_orifice(re_o, ld=4.0):
    """Lichtarowicz et al. at the orifice Reynolds number of the mean orifice velocity (Re_L = Re_o / C_d)."""
    cl = 0.75
    for _ in range(40):
        cl = lichtarowicz(re_o / cl, ld)
    return cl


def case_state(case):
    """Ideal orifice velocity, Mach and Reynolds numbers on the discrete orifice area."""
    cfg = resolved(case)
    d, gas = cfg["_derived"], cfg["gas"]
    T0, p0 = gas["initial_temperature_K"], gas["initial_pressure_Pa"]
    c0 = math.sqrt(GAMMA * R_GAS * T0)
    mu0 = gas["sutherland_As"] * T0 ** 1.5 / (T0 + gas["sutherland_Ts_K"])
    a, _ = area_ratio(case)
    uo = d["mean_speed_mps"] * d["orifice_area_ratio"] / a
    re = p0 / (R_GAS * T0) * uo * 2e-3 * cfg["geometry"]["orifice"]["radius_mm"] / mu0
    return dict(cfg=cfg, a=a, uo=uo, ma=uo / c0, re=re, re_per_ma=re / (uo / c0))


def fig_benchmark():
    """Single-orifice benchmark (orifice_benchmark.csv from orifice_benchmark.py collect) against the
    long-orifice correlation and the damper cases in the same definition of the coefficient. Runs that stopped
    before their end time, or whose mean coefficient changed by more than 1% between the third and the last
    quarter of the run (unsettled), are left out."""
    path = os.path.join(TABLES, "orifice_benchmark.csv")
    rows = list(csv.DictReader(open(path)))
    # mesh-B runs repeated with Co_max = 0.5 (orifice_benchmark_T2.csv) replace the failed runs of the first set
    path2 = os.path.join(TABLES, "orifice_benchmark_T2.csv")
    if os.path.exists(path2):
        def good(r):
            return r["Cd"] not in ("", "nan") and r["complete"] == "1" and abs(float(r["drift"])) < 0.01
        rep = {r["case"]: r for r in csv.DictReader(open(path2)) if good(r)}
        rows = [rep.get(r["case"], r) if not good(r) else r for r in rows]
    hdr = ["Re", "Cd"]
    series = {"A": ("A", "kOmegaSST"), "B": ("B", "kOmegaSST"), "Blam": ("B", "laminar"),
              "Ax": ("Ax", "kOmegaSST"), "C": ("C", "kOmegaSST"), "X01": ("X0.1", "kOmegaSST")}
    for name, (mesh, turb) in series.items():
        sel = [[float(r["Re_o"]), float(r["Cd"])] for r in rows
               if r["mesh"] == mesh and r["turbulence"] == turb and r["Cd"] not in ("", "nan")
               and r["complete"] == "1" and abs(float(r["drift"])) < 0.01]
        write("orificeBenchmark", name + ".dat", hdr, sorted(sel) or [[np.nan] * 2], [path])
    # correlation, plotted against the Reynolds number of the mean orifice velocity, Re_o = Cd Re_L
    rel = np.logspace(1.5, 6.2, 200)
    write("orificeBenchmark", "lichtarowicz.dat", hdr,
          [[lichtarowicz(x) * x, lichtarowicz(x)] for x in rel], [path])
    dm, src = [], []
    for case in finished(MATRIX + "/N*_v*") + finished(AUX + "/p20_N*_v*"):
        cfg = resolved(case)
        d = cfg["_derived"]
        if cfg["gas"]["turbulence"] == "laminar":
            continue
        rho0 = cfg["gas"]["initial_pressure_Pa"] / (R_GAS * cfg["gas"]["initial_temperature_K"])
        T0 = cfg["gas"]["initial_temperature_K"]
        mu0 = cfg["gas"]["sutherland_As"] * T0 ** 1.5 / (T0 + cfg["gas"]["sutherland_Ts_K"])
        st = case_state(case)      # coefficients and Reynolds numbers on the discrete orifice area
        dm.append([st["re"], cd_orifice(case) / st["a"], cfg["gas"]["initial_pressure_Pa"] / 1e5])
        src.append(case)
    write("orificeBenchmark", "damper.dat", ["Re", "Cd", "p0bar"], sorted(dm), src)
    # matrix with the same cross-section and axial cells of 1 mm in the orifices
    dm, src = [], []
    for case in finished(REF1 + "/N*_v*"):
        cfg = resolved(case)
        d = cfg["_derived"]
        T0 = cfg["gas"]["initial_temperature_K"]
        rho0 = cfg["gas"]["initial_pressure_Pa"] / (R_GAS * T0)
        mu0 = cfg["gas"]["sutherland_As"] * T0 ** 1.5 / (T0 + cfg["gas"]["sutherland_Ts_K"])
        st = case_state(case)
        dm.append([st["re"], cd_orifice(case) / st["a"]])
        src.append(case)
    write("orificeBenchmark", "damper1mm.dat", hdr, sorted(dm) or [[np.nan] * 2], src)
    # in-plane estimates of the damper (Section 3.1, uncertainty table): inplane_cdo.csv from verification_numbers.py
    # inplane (all conditions with three levels); without it, the two 1 m/s conditions from the Richardson
    # extrapolations from C = 16, 12, 8 (N = 10: 12, 10, 8) with h ~ 1/C and h ~ n_q^-1/2, half their distance as
    # error bar
    est = os.path.join(TABLES, "inplane_cdo.csv")
    if os.path.exists(est):
        dm = [[float(x) for x in l.split(",")[:3]] for l in open(est) if l[0].isdigit()]
        write("orificeBenchmark", "damperExt.dat", ["Re", "Cd", "err"], sorted(dm), [est])
        return rows
    dm, src = [], []
    for levels in ((AUX + "/orif_C16", AUX + "/orif_C12", MATRIX + "/N2_v1.0"),
                   (AUX + "/orifN10_C12", AUX + "/orifN10_C10", MATRIX + "/N10_v1.0")):
        cases = [os.path.join(RUNS, c) for c in levels]
        if not all(os.path.exists(os.path.join(c, "results", "summary.json")) for c in cases):
            continue
        cdo = [cd_orifice(c) / area_ratio(c)[0] for c in cases]
        nq = [json.load(open(os.path.join(c, "mesh", "mesh_report.json")))[0]["faces"]["orificeInlet"] /
              resolved(c)["geometry"]["orifice"]["count"] for c in cases]
        ests = [gci(cdo, [1.0 / resolved(c)["mesh"]["section"]["curvature_elements"] for c in cases])[1],
                gci(cdo, [x ** -0.5 for x in nq])[1]]
        dm.append([case_state(cases[0])["re"], 0.5 * (max(ests) + min(ests)), 0.5 * (max(ests) - min(ests))])
        src += cases
    hdr_ext = ["Re", "Cd", "err"]
    write("orificeBenchmark", "damperExt.dat", hdr_ext, sorted(dm) or [[np.nan] * 3], src)
    return rows


def gci(f, h):
    """Celik et al. (2008) three-level estimate; f, h ordered fine -> coarse. Returns (p, f_ext, GCI_fine).
    The order solves the signed equation e32/e21 = r21^p (r32^p - 1) / (r21^p - 1) for monotone convergence;
    without a root (no power law fits the three levels) it raises ValueError."""
    r21, r32 = h[1] / h[0], h[2] / h[1]
    e21, e32 = f[1] - f[0], f[2] - f[1]
    ratio = e32 / e21

    def res(p):
        return r21 ** p * (r32 ** p - 1.0) / (r21 ** p - 1.0) - ratio
    lo, hi = 0.01, 20.0
    if res(lo) * res(hi) > 0:
        raise ValueError("no apparent order: e32/e21 = %.3f" % ratio)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if res(lo) * res(mid) <= 0:
            hi = mid
        else:
            lo = mid
    p = 0.5 * (lo + hi)
    f_ext = (r21 ** p * f[0] - f[1]) / (r21 ** p - 1.0)
    return p, f_ext, 1.25 * abs(e21 / f[0]) / (r21 ** p - 1.0)


def table_orifice():
    """Mesh studies at N = 2 and 10 (and 6), 1 m/s: in-plane (curvature size C) and axial (cell size h_f in the
    orifices and at the piston faces) resolution against F_d, C_d and E_d. Cases with 1 mm axial cells
    (piston_cell_size_mm = 1) are the studies of the cross-section alone. Writes orifice_convergence.csv into the
    tables directory and prints the three-level grid convergence indices."""
    names = [
        # reference mesh (axial 0.1 mm): in-plane C = 6, 8, 10, 12, 16 and axial 0.2, 0.1, 0.05 mm
        AUX + "/orif_C6", MATRIX + "/N2_v1.0", AUX + "/orif_C10", AUX + "/orif_C12", AUX + "/orif_C16",
        AUX + "/ax02_N2_v1.0", AUX + "/ax005_N2_v1.0",
        AUX + "/orifN10_C6", MATRIX + "/N10_v1.0", AUX + "/orifN10_C10", AUX + "/orifN10_C12",
        AUX + "/ax02_N10_v1.0", AUX + "/ax005_N10_v1.0",
        AUX + "/rep_N2_v1.0",
        # in-plane C = 12 against the reference at other speeds and at 20 bar
        MATRIX + "/N2_v0.1", AUX + "/orif_C12_v0.1", MATRIX + "/N2_v4.0", AUX + "/orif_C12_v4.0",
        MATRIX + "/N10_v4.0", AUX + "/orifN10_C12_v4.0",
        AUX + "/p20_N2_v1.0", AUX + "/p20_N2_C12", AUX + "/p20_N10_v1.0", AUX + "/p20_N10_C12",
        # 1 mm axial cells: global refinement and in-plane refinement of the orifices
        "studies/mesh_medium", "studies/orif_C6", "studies/orif_C8", FURTHER1 + "/orif_C10", "studies/orif_C12",
        "studies/orif_C12p", "studies/mesh_fine", "studies/mesh_xfine", "matrix/N6_v1.0", "studies/orifN6_C8",
        "matrix/N10_v1.0", FURTHER1 + "/orifN10_C6", "studies/orifN10_C8", FURTHER1 + "/orifN10_C10",
        "studiesC8/orifN10_C12"]
    rows = []
    for name in names:
        case = os.path.join(RUNS, name)
        if not os.path.exists(os.path.join(case, "results", "summary.json")):
            continue
        s = summary(case)
        cfg = resolved(case) if os.path.exists(os.path.join(case, "config.resolved.json")) else None
        m0 = json.load(open(os.path.join(case, "mesh", "mesh_report.json")))[0]
        rows.append({"case": name, "N": s["orifice_count"],
                     "curvature": cfg["mesh"]["section"].get("curvature_elements") if cfg else None,
                     "axial_mm": cfg["mesh"].get("face_cell_size_mm") or cfg["mesh"]["piston_cell_size_mm"]
                     if cfg else 1.0,
                     "quads_per_orifice": m0["faces"]["orificeInlet"] // s["orifice_count"],
                     "cells_M0": m0["hexahedra"], "orifice_cells": m0["cells"]["orifices"],
                     "F_damp": s["F_damp_mean_N"], "Cd": s["Cd_effective"], "E_diss": s["E_diss_J"],
                     "Cd_o": cd_orifice(case), "steps": s.get("solver_steps"),
                     "solver_h": s["solver_wall_s"] / 3600.0, "cores": cfg["run"]["n_procs"] if cfg else None})
    with open(os.path.join(TABLES, "orifice_convergence.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    seq = {r["case"]: r for r in rows}
    # in-plane: representative size 1/C at the orifice edge (the number of quadrilaterals per orifice does not
    # grow with C**2 and gives no power law); axial: the cell size h_f; three finest levels, fine -> coarse
    sets = (("N2 in-plane", (AUX + "/orif_C16", AUX + "/orif_C12", AUX + "/orif_C10"), "curvature"),
            ("N10 in-plane", (AUX + "/orifN10_C12", AUX + "/orifN10_C10", MATRIX + "/N10_v1.0"), "curvature"),
            ("N2 axial", (AUX + "/ax005_N2_v1.0", MATRIX + "/N2_v1.0", AUX + "/ax02_N2_v1.0"), "axial"),
            ("N10 axial", (AUX + "/ax005_N10_v1.0", MATRIX + "/N10_v1.0", AUX + "/ax02_N10_v1.0"), "axial"),
            ("N2 in-plane 1mm", ("studies/orif_C12", FURTHER1 + "/orif_C10", "studies/orif_C8"), "curvature"),
            ("N10 in-plane 1mm", ("studiesC8/orifN10_C12", FURTHER1 + "/orifN10_C10", "studies/orifN10_C8"),
             "curvature"))
    for label, names3, kind in sets:
        fine = [seq.get(n) for n in names3]
        if not all(fine):
            continue
        h = [1.0 / r["curvature"] if kind == "curvature" else r["axial_mm"] for r in fine]
        ref = MATRIX + "/N2_v1.0" if label.startswith("N2 ") else MATRIX + "/N10_v1.0"
        for key in ("F_damp", "Cd", "Cd_o", "E_diss"):
            try:
                p, ext, g = gci([r[key] for r in fine], h)
            except ValueError as e:
                print("GCI %-16s %-6s %s" % (label, key, e))
                continue
            dev = "  ref %+.1f%%" % (100 * (seq[ref][key] / ext - 1)) if ref in seq and "1mm" not in label else ""
            print("GCI %-16s %-6s p = %.2f  extrapolated %.5g  GCI_fine %.1f%%  coarsest-of-3 %+.1f%%%s"
                  % (label, key, p, ext, 100 * g, 100 * (fine[2][key] / ext - 1), dev))
    return rows


def table_cases():
    """One row per finished case of all studies under the runs directory, for the robustness and cost table;
    written to cases_table.csv in the tables directory."""
    rows = []
    for case in finished("*/*"):
        s, st = summary(case), json.load(open(os.path.join(case, "state.json")))
        stages = st.get("stages", {})
        cfg = resolved(case) if os.path.exists(os.path.join(case, "config.resolved.json")) else None
        rows.append({
            "case": os.path.relpath(case, RUNS), "N": s["orifice_count"], "v": s["mean_speed_mps"],
            # solve stage on the compute node when the log of the batch scheduler is present, else on the
            # workstation (the two machines of Section 3 of the paper)
            "machine": "compute node" if glob.glob(os.path.join(case, BATCH_LOG)) else "workstation",
            "cores": cfg["run"]["n_procs"] if cfg else None,
            "p0": s["initial_pressure_Pa"], "model": s["turbulence"],
            "steps": s.get("solver_steps"), "solver_h": s.get("solver_wall_s", 0) / 3600.0,
            "mesh_s": stages.get("mesh", {}).get("wall_s"),
            "total_h": sum(v.get("wall_s", 0) for v in stages.values()) / 3600.0,
            "mass_drift": s.get("mass_drift_max_rel"), "max_jump": max(map(abs, s.get("mass_jumps_rel") or [0])),
            "V_err": s.get("volume_error_max_rel"), "Co_max": s.get("Co_max"),
            "F_damp": s.get("F_damp_mean_N"), "E_diss": s.get("E_diss_J"),
            "n_poly": s.get("n_polytropic"), "Cd": s.get("Cd_effective"),
            "ratio": s.get("damping_to_spring_ratio"), "Ma_peak": s.get("Ma_orifice_mean_peak"),
        })
    with open(os.path.join(TABLES, "cases_table.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return rows


OUTPUTS = {"loops": fig_loops, "ratio": fig_ratio, "cd": fig_cd, "benchmark": fig_benchmark, "mass": fig_mass,
           "spring": fig_spring, "history": fig_history, "energy": table_energy, "orifice": table_orifice,
           "cases": table_cases}


def main(argv=None):
    ap = argparse.ArgumentParser(description="Data files of the paper figures and derived case tables.")
    ap.add_argument("--runs", default="runs", help="directory holding the study directories (default: runs)")
    ap.add_argument("--tables", default="tables", help="directory of the derived CSV tables (default: tables)")
    ap.add_argument("--out", default="figures", help="directory of the figure data (default: figures)")
    ap.add_argument("outputs", nargs="*", metavar="output",
                    help="any of %s (default: all)" % " ".join(OUTPUTS))
    args = ap.parse_args(argv)
    bad = [o for o in args.outputs if o not in OUTPUTS]
    if bad:
        ap.error("unknown output(s): %s" % " ".join(bad))
    configure(args.runs, args.tables, args.out)
    os.makedirs(TABLES, exist_ok=True)
    res = {}
    for name, fn in OUTPUTS.items():
        if not args.outputs or name in args.outputs:
            res[name] = fn()
    if "ratio" in res and "cases" in res:
        print("ratio points:", len(res["ratio"]), " cases in table:", len(res["cases"]))


if __name__ == "__main__":
    sys.exit(main())
