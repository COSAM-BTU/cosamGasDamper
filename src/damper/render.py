# SPDX-License-Identifier: GPL-3.0-or-later
"""Headless rendering of mid-plane slices with ParaView's pvbatch."""

import json
import math
import os
import shutil
import subprocess

from . import config as C

SCRIPT = os.path.join(C.REPO, "src", "render", "render_slices.py")


def time_label_value(label, T):
    return {"compression_mid": 0.5 * T, "compression_end": T,
            "rebound_mid": 1.5 * T, "cycle_end": 2.0 * T}[label]


def run(case):
    cfg, d = case.cfg, case.d
    T = d["stroke_time_s"]
    P = int(cfg["numerics"]["time_precision"])
    segs = d["schedule"]["segments"]
    a = math.radians(cfg["geometry"]["orifice"]["start_angle_deg"])
    jobs = []
    for label in cfg["render"]["times"]:
        t = C.tround(time_label_value(label, T), P)
        seg = next((s for s in segs if s["start"] < t <= s["end"] + 1e-12), segs[0])
        jobs.append({"label": label, "time": t,
                     "case": os.path.join(case.dir, "segments", seg["name"], "case.foam")})
    out = os.path.join(case.dir, "frames")
    os.makedirs(out, exist_ok=True)
    spec = {"jobs": jobs, "fields": cfg["render"]["fields"], "out_dir": out,
            "resolution": cfg["render"]["resolution"],
            "normal": [0.0, -math.sin(a), math.cos(a)],
            "tangent": [0.0, math.cos(a), math.sin(a)],
            "radius_m": cfg["geometry"]["cylinder_radius_mm"] * 1e-3,
            "ranges": cfg["render"].get("ranges") or {},
            "colormap": cfg["render"].get("colormap") or "Viridis (matplotlib)",
            "scalar_bar": cfg["render"].get("scalar_bar", True),
            "length_m": cfg["geometry"]["cylinder_length_mm"] * 1e-3}
    spath = os.path.join(out, "render.json")
    with open(spath, "w") as f:
        json.dump(spec, f, indent=1)
    log = os.path.join(case.logdir, "render.log")
    env = dict(os.environ)
    env["PYTHONNOUSERSITE"] = "1"  # keep user-site VTK wheels away from pvbatch
    cmd = [cfg["run"]["pvbatch"], SCRIPT, spath]
    if not env.get("DISPLAY") and shutil.which("xvfb-run"):
        # no display: render through a virtual X server
        cmd = ["xvfb-run", "-a", "-s", "-screen 0 1920x1080x24"] + cmd
    with open(log, "w") as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
    pngs = sorted(p for p in os.listdir(out) if p.endswith(".png"))
    case.say("render: %d images (pvbatch exit %d)" % (len(pngs), r.returncode))
    if r.returncode != 0 or not pngs:
        from .pipeline import StageError
        raise StageError("render failed (pvbatch exit %d), see %s" % (r.returncode, log))
    return {"images": pngs, "exit": r.returncode}
