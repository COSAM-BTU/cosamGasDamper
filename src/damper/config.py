# SPDX-License-Identifier: GPL-3.0-or-later
"""Configuration loading, validation and derived quantities."""

import copy
import hashlib
import json
import math
import os
import re

import yaml

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULTS = os.path.join(REPO, "config", "defaults.yaml")


class ConfigError(ValueError):
    pass


def deep_merge(base, over):
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def set_dotted(cfg, key, value):
    d = cfg
    parts = key.split(".")
    for p in parts[:-1]:
        d = d.setdefault(p, {})
    d[parts[-1]] = value


def get_dotted(cfg, key):
    d = cfg
    for p in key.split("."):
        d = d[p]
    return d


_NUM = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")


def coerce_numbers(obj):
    """PyYAML (YAML 1.1) reads e.g. 1.0e5 as a string; turn such strings into floats."""
    if isinstance(obj, dict):
        return {k: coerce_numbers(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [coerce_numbers(v) for v in obj]
    if isinstance(obj, str) and _NUM.match(obj.strip()):
        return float(obj)
    return obj


def load(path_or_dict):
    with open(DEFAULTS) as f:
        cfg = yaml.safe_load(f)
    user = path_or_dict
    if isinstance(path_or_dict, str):
        with open(path_or_dict) as f:
            user = yaml.safe_load(f) or {}
    user = copy.deepcopy(user)
    user.pop("sweep", None)
    um = user.get("motion") or {}
    # the default mean speed yields to a user-given stroke time (and vice versa)
    if um.get("stroke_time_s") is not None and "mean_speed_mps" not in um:
        cfg["motion"]["mean_speed_mps"] = None
    return coerce_numbers(deep_merge(cfg, user))


def tname(t, precision):
    """OpenFOAM time name for timeFormat general and the given timePrecision."""
    return "%.*g" % (precision, t)


def tround(t, precision):
    return float(tname(t, precision))


def derive(cfg):
    """Validate the configuration and compute everything the stages need."""
    g, m, me, nu = cfg["geometry"], cfg["motion"], cfg["mesh"], cfg["numerics"]
    o = g["orifice"]
    rs, rm = g["cylinder_radius_mm"], g["rod_radius_mm"]
    L, lp = g["cylinder_length_mm"], g["piston_length_mm"]
    N, ro, rp = int(o["count"]), o["radius_mm"], o["pitch_radius_mm"]
    S = m["stroke_mm"]
    K = int(m["meshes_per_stroke"])

    errs = []
    if N < 1:
        errs.append("orifice.count must be >= 1")
    if not rm + 0.2 * ro < rp - ro:
        errs.append("orifices overlap the rod: rod_radius < pitch_radius - orifice_radius required")
    if not rp + ro < rs - 0.2 * ro:
        errs.append("orifices cut the cylinder wall: pitch_radius + orifice_radius < cylinder_radius")
    if N >= 2 and not 2.0 * rp * math.sin(math.pi / N) > 2.4 * ro:
        errs.append("neighbouring orifices overlap (2 r_p sin(pi/N) <= 2.4 r_o)")
    if (m.get("stroke_time_s") is None) == (m.get("mean_speed_mps") is None):
        errs.append("give exactly one of motion.stroke_time_s and motion.mean_speed_mps")
    if K < 1:
        errs.append("motion.meshes_per_stroke must be >= 1")
    if m["profile"] not in ("constant", "sinusoidal"):
        errs.append("motion.profile must be constant or sinusoidal")
    if nu.get("map_method", "interpolate") not in ("interpolate", "cellVolumeWeight"):
        errs.append("numerics.map_method must be interpolate or cellVolumeWeight")
    x0 = m.get("initial_piston_front_mm")
    if x0 is None:
        x0 = L / 2.0 + S / 2.0 - lp / 2.0
    hmin = 2.0 * me["chamber_cell_size_mm"]
    if x0 - S < hmin:
        errs.append("piston hits the bottom cap: x0 - stroke < 2 cells")
    if x0 + lp > L - hmin:
        errs.append("piston hits the top cap: x0 + piston_length > L - 2 cells")
    if errs:
        raise ConfigError("; ".join(errs))

    if m.get("stroke_time_s") is not None:
        T = float(m["stroke_time_s"])
        v = S * 1e-3 / T
    else:
        v = float(m["mean_speed_mps"])
        T = S * 1e-3 / v
    P = int(nu["time_precision"])
    T = tround(T, P)

    xk = [x0 - k * S / K for k in range(K + 1)]
    positions = []
    hf = me.get("face_cell_size_mm")
    for k in range(K):
        xf = xk[k]
        pos = {"name": "M%d" % k, "index": k, "piston_front_x": xf,
               "segments_chamber1": max(2, int(round(xf / me["chamber_cell_size_mm"]))),
               "segments_chamber2": max(2, int(round((L - xf - lp) / me["chamber_cell_size_mm"])))}
        if hf:
            q = float(me.get("face_cell_growth") or 1.15)
            pos["segments_chamber1"], pos["density_chamber1"] = graded_path(
                xf, float(hf), me["chamber_cell_size_mm"], q, True)
            pos["segments_chamber2"], pos["density_chamber2"] = graded_path(
                L - xf - lp, float(hf), me["chamber_cell_size_mm"], q, False)
        positions.append(pos)

    sched = schedule(T, S, K, int(m["cycles"]), m["profile"], P)
    a_rod = math.pi * rm ** 2
    a_cyl = math.pi * rs ** 2
    a_orf = N * math.pi * ro ** 2

    def volume(xf):  # mm^3, analytic gas volume for piston front at xf
        return a_cyl * xf + a_orf * lp + (a_cyl - a_rod) * (L - xf - lp)

    gas = cfg["gas"]
    R = 8314.462618 / gas["molar_mass_kg_per_kmol"]
    gamma = gas["cp_J_per_kgK"] / (gas["cp_J_per_kgK"] - R)
    max_dt = nu["max_dt_s"] if nu["max_dt_s"] is not None else T / 400.0
    write_dt = tround(T / nu["writes_per_stroke"], P)
    d = {
        "stroke_time_s": T, "mean_speed_mps": v, "x0_mm": x0, "grid_mm": xk,
        "positions": positions, "schedule": sched,
        "end_time_s": sched["end_time"], "max_dt_s": max_dt, "write_interval_s": write_dt,
        "areas_mm2": {"rod": a_rod, "cylinder": a_cyl, "orifices": a_orf,
                      "piston_front": a_cyl - a_orf, "piston_back": a_cyl - a_rod - a_orf,
                      "annulus": a_cyl - a_rod},
        "volume_initial_mm3": volume(x0), "volume_min_mm3": volume(x0 - S),
        "R_J_per_kgK": R, "gamma": gamma,
        "orifice_area_ratio": (a_cyl - a_rod) / a_orf,
    }
    return d


def graded_path(length, h_face, h_far, growth, face_at_end):
    """Axial distribution of one chamber with cells of h_face next to the piston face, growing geometrically by
    `growth` to h_far: number of segments and the density table (t, 1/h) of the SALOME 1D distribution."""
    def size(s):                                    # local cell size at distance s from the piston face
        return min(h_face + (growth - 1.0) * s, h_far)
    s_star = (h_far - h_face) / (growth - 1.0)
    if length <= s_star:
        n = math.log(size(length) / h_face) / (growth - 1.0)
    else:
        n = math.log(h_far / h_face) / (growth - 1.0) + (length - s_star) / h_far
    # table points clustered at the face, where the size changes, and uniform beyond
    near = min(length, 2.0 * s_star)
    dist = sorted(set([near * (i / 60.0) ** 2 for i in range(61)]
                      + [near + (length - near) * i / 20.0 for i in range(1, 21)]))
    pts = sorted(((1.0 - s / length) if face_at_end else s / length, 1.0 / size(s)) for s in dist)
    table = [v for pt in pts for v in pt]
    return max(2, int(round(n))), table


def displacement(profile, tau, T, S):
    """Stroke displacement s(tau) in mm for 0 <= tau <= T (monotonic)."""
    if profile == "constant":
        return S * tau / T
    return 0.5 * S * (1.0 - math.cos(math.pi * tau / T))


def speed(profile, tau, T, S):
    """ds/dtau in mm/s."""
    if profile == "constant":
        return S / T
    return 0.5 * S * math.pi / T * math.sin(math.pi * tau / T)


def grid_time(profile, k, K, T):
    """Time within a stroke at which s = k S / K."""
    if profile == "constant":
        return k * T / K
    return T / math.pi * math.acos(1.0 - 2.0 * k / K)


def piston_x(profile, t, T, S, x0):
    """Piston-front position (mm) at time t; compression first (towards x = 0)."""
    n, tau = divmod(t, 2.0 * T)
    if tau <= T:
        return x0 - displacement(profile, tau, T, S)
    return x0 - S + displacement(profile, tau - T, T, S)


def piston_u(profile, t, T, S):
    """Piston velocity (m/s, along +x) at time t."""
    n, tau = divmod(t, 2.0 * T)
    if tau < T:
        return -speed(profile, tau, T, S) * 1e-3
    return speed(profile, tau - T, T, S) * 1e-3


def schedule(T, S, K, cycles, profile, P):
    """Mesh segments: which sub-mesh is used over which time interval.

    Events are the instants at which the piston reaches a grid position
    x_k = x0 - k S/K.  Between two events the piston moves between adjacent
    grid points.  A sub-mesh M_j (natural position x_j) is kept as long as the
    piston stays within one grid spacing of x_j; otherwise the next segment
    starts on the sub-mesh whose natural position equals the current one.
    """
    events = []  # (time, grid index)
    for c in range(cycles):
        base = 2.0 * T * c
        for k in range(K + 1):
            events.append((base + grid_time(profile, k, K, T), k))
        for k in range(K - 1, -1, -1):
            events.append((base + 2.0 * T - grid_time(profile, k, K, T), k))
    ev = []
    for t, k in events:
        t = tround(t, P)
        if not ev or t > ev[-1][0]:
            ev.append((t, k))
    segments = []
    cur = None
    for (ta, ga), (tb, gb) in zip(ev[:-1], ev[1:]):
        if cur is not None and abs(ga - cur["mesh"]) <= 1 and abs(gb - cur["mesh"]) <= 1:
            cur["end"] = tb
            cur["grid_path"].append(gb)
            continue
        mesh = min(ga, K - 1)
        cur = {"mesh": mesh, "start": ta, "end": tb, "grid_path": [ga, gb],
               "init": "initial" if not segments else "map"}
        segments.append(cur)
    for i, s in enumerate(segments):
        s["name"] = "seg%02d_M%d" % (i + 1, s["mesh"])
        s["start_name"] = tname(s["start"], P)
        s["end_name"] = tname(s["end"], P)
    return {"segments": segments, "events": ev, "end_time": ev[-1][0]}


def velocity_table(cfg, d):
    """Piecewise-linear table of the piston velocity (m/s) over the whole run."""
    m = cfg["motion"]
    T, S = d["stroke_time_s"], m["stroke_mm"]
    P = int(cfg["numerics"]["time_precision"])
    eps = 1e-9 * T
    pts = []
    if m["profile"] == "constant":
        # left-continuous at every reversal: a time step that ends exactly at the
        # reversal still uses the old velocity, the next one the new velocity
        v = S * 1e-3 / T
        for c in range(int(m["cycles"])):
            b = 2.0 * T * c
            pts += [(b + (eps if c else 0.0), -v), (b + T, -v),
                    (b + T + eps, v), (b + 2 * T, v)]
    else:
        n = int(m["velocity_table_points_per_stroke"])
        for c in range(int(m["cycles"])):
            for i in range(2 * n + 1):
                t = 2.0 * T * c + i * T / n
                if pts and t <= pts[-1][0]:
                    continue
                pts.append((t, piston_u("sinusoidal", min(t, 2 * T * (c + 1) - eps), T, S)))
    return [(float("%.*g" % (P + 4, t)), u) for t, u in pts]


def stage_hash(cfg, keys):
    sub = {k: cfg.get(k) for k in keys}
    return hashlib.sha1(json.dumps(sub, sort_keys=True).encode()).hexdigest()[:12]
