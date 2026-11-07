# SPDX-License-Identifier: GPL-3.0-or-later
"""Collect function-object output of all mesh segments into CSV files and a summary."""

import csv
import glob
import json
import math
import os
import re

import numpy as np

from . import config as C

trapezoid = getattr(np, "trapezoid", None) or np.trapz  # numpy < 2.0 has only trapz

FO = {  # function object -> (file, {column-substring: short name})
    "forces": ("force.dat", {"total_x": "Fx", "pressure_x": "Fx_p", "viscous_x": "Fx_v"}),
    "forcesFront": ("force.dat", {"total_x": "Fx_front"}),
    "forcesBack": ("force.dat", {"total_x": "Fx_back"}),
    # v2406 writes zone averages as volAverage(p) etc. (the zone name is not in the header)
    "chamber1": ("volFieldValue.dat", {"(p)": "p1", "(T)": "T1", "(rho)": "rho1"}),
    "chamber2": ("volFieldValue.dat", {"(p)": "p2", "(T)": "T2", "(rho)": "rho2"}),
    "orifices": ("volFieldValue.dat", {"(p)": "pO", "(T)": "TO", "(rho)": "rhoO"}),
    "all": ("volFieldValue.dat", {"(p)": "p_avg", "(T)": "T_avg", "(rho)": "rho_avg"}),
    "mass1": ("volFieldValue.dat", {"rho)": "m1"}),
    "mass2": ("volFieldValue.dat", {"rho)": "m2"}),
    "massOrifices": ("volFieldValue.dat", {"rho)": "mO"}),
    "massAll": ("volFieldValue.dat", {"rho)": "m_total"}),
    "orificeFlux": ("surfaceFieldValue.dat", {"phi)": "mdot"}),
    "orificeState": ("surfaceFieldValue.dat", {"(p)": "p_exit", "(T)": "T_exit", "(rho)": "rho_exit"}),
    "orificeUmax": ("volFieldValue.dat", {"mag(U))": "Umax_orifice"}),
}


def read_dat(path):
    header, rows = None, []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if "Time" in line:
                    header = line[1:].split()
                continue
            vals = line.replace("(", " ").replace(")", " ").split()
            if vals:
                rows.append([float(v) for v in vals])
    if header is None or not rows:
        return None, None
    return header, np.array(rows)


def force_columns(header):
    # v2406 force.dat: Time total_x total_y total_z pressure_x ... viscous_x ...
    return {h: i for i, h in enumerate(header)}


def load_segment(seg_dir, start_name):
    out = {}
    for fo, (fname, cols) in FO.items():
        files = glob.glob(os.path.join(seg_dir, "postProcessing", fo, "*", fname))
        if not files:
            continue
        header, data = read_dat(sorted(files)[0])
        if data is None:
            continue
        idx = {}
        for sub, short in cols.items():
            for i, h in enumerate(header):
                if sub in h:
                    idx[short] = i
                    break
        out.setdefault("time", data[:, 0])
        n = min(len(out["time"]), len(data))
        for short, i in idx.items():
            out[short] = data[:n, i]
    return out


def courant_from_log(log):
    t, dt, co = [], [], []
    cur_co, cur_dt = None, None
    with open(log) as f:
        for line in f:
            if line.startswith("Courant Number"):
                cur_co = float(line.split("max:")[1])
            elif line.startswith("deltaT ="):
                cur_dt = float(line.split("=")[1])
            elif line.startswith("Time ="):
                t.append(float(line.split("=")[1]))
                co.append(cur_co if cur_co is not None else np.nan)
                dt.append(cur_dt if cur_dt is not None else np.nan)
    return np.array(t), np.array(dt), np.array(co)


def collect(case):
    cfg, d = case.cfg, case.d
    g, m, gas = cfg["geometry"], cfg["motion"], cfg["gas"]
    T, S, x0 = d["stroke_time_s"], m["stroke_mm"], d["x0_mm"]
    prof = m["profile"]
    a = d["areas_mm2"]
    A_rod = a["rod"] * 1e-6
    A_orf = a["orifices"] * 1e-6
    p0, p_amb = gas["initial_pressure_Pa"], gas["ambient_pressure_Pa"]
    R, gam = d["R_J_per_kgK"], d["gamma"]
    rdir = os.path.join(case.dir, "results")
    os.makedirs(rdir, exist_ok=True)

    rows, seg_info = [], []
    for k, seg in enumerate(d["schedule"]["segments"]):
        sd = os.path.join(case.dir, "segments", seg["name"])
        s = load_segment(sd, seg["start_name"])
        if not s:
            continue
        lt, ldt, lco = courant_from_log(os.path.join(sd, "log.rhoPimpleFoam"))
        n = len(s["time"])
        for i in range(n):
            t = s["time"][i]
            r = {"time_s": t, "segment": k + 1, "mesh": seg["mesh"]}
            for key, arr in s.items():
                if key != "time" and i < len(arr):
                    r[key] = arr[i]
            j = np.searchsorted(lt, t - 1e-15)
            if j < len(lt) and abs(lt[j] - t) <= 1e-12 * max(1.0, t):
                r["dt_s"], r["Co_max"] = ldt[j], lco[j]
            rows.append(r)
        seg_info.append({"segment": k + 1, "name": seg["name"], "mesh": seg["mesh"],
                         "t_start": seg["start"], "t_end": seg["end"],
                         "m_first": float(s["m_total"][0]) if "m_total" in s else None,
                         "m_last": float(s["m_total"][-1]) if "m_total" in s else None})

    for r in rows:
        t = r["time_s"]
        xp = C.piston_x(prof, t, T, S, x0)
        r["piston_x_mm"] = xp
        r["piston_u_mps"] = C.piston_u(prof, t, T, S)
        Va = (a["cylinder"] * xp + a["orifices"] * g["piston_length_mm"]
              + a["annulus"] * (g["cylinder_length_mm"] - xp - g["piston_length_mm"])) * 1e-9
        r["V_analytic_m3"] = Va
        if "m_total" in r and "rho_avg" in r and r["rho_avg"] > 0:
            r["V_cfd_m3"] = r["m_total"] / r["rho_avg"]
        if "Fx_p" in r:
            r["Fx_gauge"] = r["Fx_p"] - p_amb * A_rod + r.get("Fx_v", 0.0)
            r["Fx_rel_p0"] = r["Fx_p"] - p0 * A_rod + r.get("Fx_v", 0.0)
        if "mdot" in r and "rho_exit" in r and "T_exit" in r and r["rho_exit"] > 0:
            u_mean = abs(r["mdot"]) / (r["rho_exit"] * A_orf)
            c = math.sqrt(gam * R * r["T_exit"])
            r["U_orifice_mean_mps"] = u_mean
            r["Ma_orifice_mean"] = u_mean / c
        if "Umax_orifice" in r and "TO" in r:
            r["Ma_orifice_max"] = r["Umax_orifice"] / math.sqrt(gam * R * r["TO"])
        if "p1" in r and "p2" in r:
            r["dp_12_Pa"] = r["p1"] - r["p2"]

    keys = ["time_s", "segment", "mesh", "piston_x_mm", "piston_u_mps", "Fx", "Fx_p", "Fx_v",
            "Fx_gauge", "Fx_rel_p0", "Fx_front", "Fx_back", "p_avg", "T_avg", "rho_avg",
            "p1", "T1", "rho1", "p2", "T2", "rho2", "pO", "TO", "rhoO", "dp_12_Pa",
            "m1", "m2", "mO", "m_total", "V_cfd_m3", "V_analytic_m3", "mdot",
            "p_exit", "T_exit", "rho_exit", "U_orifice_mean_mps", "Ma_orifice_mean",
            "Umax_orifice", "Ma_orifice_max", "dt_s", "Co_max"]
    with open(os.path.join(rdir, "timeseries.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for r in rows:
            w.writerow([("%.10g" % r[k]) if isinstance(r.get(k), (float, np.floating))
                        else r.get(k, "") for k in keys])

    summ = summarise(case, rows, seg_info)
    with open(os.path.join(rdir, "summary.json"), "w") as f:
        json.dump(summ, f, indent=1)
    with open(os.path.join(rdir, "segments.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(seg_info[0].keys()) + ["mass_jump_rel"])
        w.writeheader()
        for i, s in enumerate(seg_info):
            s = dict(s)
            if i and seg_info[i - 1]["m_last"] and s["m_first"]:
                s["mass_jump_rel"] = s["m_first"] / seg_info[i - 1]["m_last"] - 1.0
            w.writerow(s)
    try:
        plots(case, rows, rdir)
    except Exception as e:  # plotting is a convenience, never fail the case on it
        summ["plot_error"] = str(e)
    case.say("post: Fx_damp_mean = %.4g N, E_diss = %.4g J, n_poly = %.4f, max |dm/m0| = %.2e"
             % (summ.get("F_damp_mean_N", float("nan")), summ.get("E_diss_J", float("nan")),
                summ.get("n_polytropic", float("nan")), summ.get("mass_drift_max_rel", float("nan"))))
    return summ


def branch(rows, key, t0, t1, skip=()):
    """(x, value) of one stroke, sorted by x; time windows in `skip` are left out."""
    pts = [(r["piston_x_mm"], r[key]) for r in rows if t0 < r["time_s"] < t1 and key in r
           and not any(a <= r["time_s"] <= b for a, b in skip)]
    pts.sort()
    if len(pts) < 4:
        return None, None
    x, y = np.array(pts).T
    return x, y


def summarise(case, rows, seg_info):
    cfg, d = case.cfg, case.d
    T, S, x0 = d["stroke_time_s"], cfg["motion"]["stroke_mm"], d["x0_mm"]
    p0 = cfg["gas"]["initial_pressure_Pa"]
    gam = d["gamma"]
    A_rod = d["areas_mm2"]["rod"] * 1e-6
    out = {"case": cfg["case"]["name"], "orifice_count": cfg["geometry"]["orifice"]["count"],
           "mean_speed_mps": d["mean_speed_mps"], "stroke_time_s": T,
           "initial_pressure_Pa": p0, "profile": cfg["motion"]["profile"],
           "turbulence": cfg["gas"]["turbulence"], "n_rows": len(rows)}
    if not rows:
        return out
    m = np.array([r.get("m_total", np.nan) for r in rows])
    out["mass_initial_kg"] = float(m[0])
    out["mass_drift_max_rel"] = float(np.nanmax(np.abs(m / m[0] - 1.0)))
    jumps = []
    for i in range(1, len(seg_info)):
        if seg_info[i - 1]["m_last"] and seg_info[i]["m_first"]:
            jumps.append(seg_info[i]["m_first"] / seg_info[i - 1]["m_last"] - 1.0)
    out["mass_jumps_rel"] = jumps
    vc = np.array([r.get("V_cfd_m3", np.nan) for r in rows])
    va = np.array([r["V_analytic_m3"] for r in rows])
    out["volume_error_max_rel"] = float(np.nanmax(np.abs(vc / va - 1.0)))
    co = np.array([r.get("Co_max", np.nan) for r in rows], dtype=float)
    out["Co_max"] = float(np.nanmax(co)) if np.isfinite(co).any() else None
    for key in ("Ma_orifice_max", "Ma_orifice_mean", "Umax_orifice"):
        v = np.array([r.get(key, np.nan) for r in rows], dtype=float)
        if np.isfinite(v).any():
            out[key + "_peak"] = float(np.nanmax(v))

    # polytropic exponent of the compression stroke: p_avg * V^n = const
    sel = [(r["p_avg"], r["V_cfd_m3"]) for r in rows
           if r["time_s"] <= T and "p_avg" in r and "V_cfd_m3" in r]
    if len(sel) > 10:
        p, V = np.array(sel).T
        out["n_polytropic"] = float(-np.polyfit(np.log(V), np.log(p), 1)[0])
        out["p_end_compression_Pa"] = float(p[-1])
        out["p_end_adiabatic_Pa"] = float(p0 * (d["volume_initial_mm3"] / d["volume_min_mm3"]) ** gam)

    # spring / damping decomposition on the first cycle; the numerical transient that
    # follows each mesh switch (a few acoustic transit times) is left out
    c0 = math.sqrt(gam * d["R_J_per_kgK"] * cfg["gas"]["initial_temperature_K"])
    t_ac = cfg["geometry"]["cylinder_length_mm"] * 1e-3 / c0
    t_skip = 5.0 * t_ac
    switches = [s["start"] for s in d["schedule"]["segments"][1:]]
    skip = [(t, t + t_skip) for t in switches]
    out["acoustic_time_s"] = t_ac
    out["switch_skip_s"] = t_skip
    xc, Fc = branch(rows, "Fx", 0.0, T, skip)
    xr, Fr = branch(rows, "Fx", T, 2.0 * T, skip)
    if xc is not None and xr is not None:
        lo, hi = x0 - S, x0
        xs = np.linspace(lo + 0.05 * S, hi - 0.05 * S, 200)
        fc = np.interp(xs, xc, Fc)
        fr = np.interp(xs, xr, Fr)
        fd = 0.5 * (fc - fr)
        fs = 0.5 * (fc + fr)
        xf = np.linspace(lo + 1e-6, hi - 1e-6, 400)
        dE = np.interp(xf, xc, Fc) - np.interp(xf, xr, Fr)
        out["E_diss_J"] = float(trapezoid(dE, xf * 1e-3))
        out["F_damp_mean_N"] = float(np.mean(fd))
        out["F_damp_max_N"] = float(np.max(fd))
        out["F_spring_mid_N"] = float(np.interp(x0 - 0.5 * S, xs, fs))
        out["F_spring_range_N"] = float(fs.max() - fs.min())
        out["F_start_N"] = float(Fc[np.argmax(xc)])
        out["F_end_compression_N"] = float(Fc[np.argmin(xc)])
        out["damping_to_spring_ratio"] = out["F_damp_mean_N"] / (p0 * A_rod)
    # pressure drop across the piston at mid-compression
    mid = [r for r in rows if 0.4 * T < r["time_s"] < 0.6 * T and "dp_12_Pa" in r]
    if mid:
        dp = np.mean([r["dp_12_Pa"] for r in mid])
        rho = np.mean([r["rho2"] for r in mid])
        # ideal orifice velocity from the instantaneous piston speed (sinusoidal motion
        # is faster than the mean at mid-stroke)
        u_ideal = np.mean([abs(r["piston_u_mps"]) for r in mid]) * d["orifice_area_ratio"]
        out["dp_mid_compression_Pa"] = float(dp)
        if dp > 0:
            out["Cd_effective"] = float(u_ideal / math.sqrt(2.0 * dp / rho))
    st = case.state.get("segments", {})
    out["solver_wall_s"] = float(sum(v.get("wall_s", 0.0) for v in st.values()))
    out["solver_steps"] = int(sum(v.get("steps", 0) for v in st.values()))
    return out


def plots(case, rows, rdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    T = case.d["stroke_time_s"]
    t = np.array([r["time_s"] for r in rows])
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    for t0, t1, lab in ((0, T, "compression"), (T, 2 * T, "rebound")):
        sel = [(r["piston_x_mm"], r["Fx"]) for r in rows if t0 < r["time_s"] < t1 and "Fx" in r]
        if sel:
            x, f = np.array(sel).T
            ax[0].plot(x, f, label=lab)
    ax[0].set_xlabel("piston-front position x [mm]")
    ax[0].set_ylabel("F_x [N]")
    ax[0].legend()
    m = np.array([r.get("m_total", np.nan) for r in rows])
    ax[1].plot(t, m / m[0] - 1.0)
    ax[1].set_xlabel("t [s]")
    ax[1].set_ylabel("m/m0 - 1")
    for key in ("p1", "p2"):
        ax[2].plot(t, [r.get(key, np.nan) for r in rows], label=key)
    ax[2].set_xlabel("t [s]")
    ax[2].set_ylabel("p [Pa]")
    ax[2].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(rdir, "quicklook.png"), dpi=110)
    plt.close(fig)
