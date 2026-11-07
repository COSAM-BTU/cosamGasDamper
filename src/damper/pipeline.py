# SPDX-License-Identifier: GPL-3.0-or-later
"""Stages of one damper case: mesh -> solve (mesh segments + mapping) -> post -> render."""

import fcntl
import glob
import json
import os
import re
import shutil
import subprocess
import time

from . import config as C
from . import foamcase as F

REPO = C.REPO
STAGES = ("mesh", "solve", "post", "render")


class StageError(RuntimeError):
    pass


def num(text, pattern):
    """First float matched by pattern (OpenFOAM often ends numbers with a full stop)."""
    return float(re.search(pattern, text).group(1).rstrip("."))


class Case:
    def __init__(self, case_dir, cfg_path=None):
        self.dir = os.path.abspath(case_dir)
        os.makedirs(self.dir, exist_ok=True)
        src = cfg_path or os.path.join(self.dir, "config.yaml")
        if cfg_path and os.path.abspath(cfg_path) != os.path.join(self.dir, "config.yaml"):
            shutil.copy(cfg_path, os.path.join(self.dir, "config.yaml"))
        self.cfg = C.load(src)
        self.d = C.derive(self.cfg)
        self.cfg["_derived"] = self.d
        self.tab = C.velocity_table(self.cfg, self.d)
        self.logdir = os.path.join(self.dir, "logs")
        os.makedirs(self.logdir, exist_ok=True)
        self.state_path = os.path.join(self.dir, "state.json")
        self.state = json.load(open(self.state_path)) if os.path.exists(self.state_path) else {}
        with open(os.path.join(self.dir, "config.resolved.json"), "w") as f:
            json.dump(self.cfg, f, indent=1)

    # ------------------------------------------------------------------ utils
    def save_state(self):
        tmp = self.state_path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(self.state, f, indent=1)
        os.replace(tmp, self.state_path)

    def say(self, msg):
        line = "[%s %s] %s" % (time.strftime("%H:%M:%S"), self.cfg["case"]["name"], msg)
        print(line, flush=True)
        with open(os.path.join(self.dir, "run.log"), "a") as f:
            f.write(line + "\n")

    def foam(self, cmd, cwd, log):
        """Run an OpenFOAM command with the configured environment."""
        bashrc = self.cfg["run"]["openfoam_bashrc"]
        full = "source %s >/dev/null 2>&1; set -o pipefail; %s" % (bashrc, cmd)
        with open(log, "w") as f:
            r = subprocess.run(["bash", "-c", full], cwd=cwd, stdout=f, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            tail = open(log).read()[-3000:]
            raise StageError("command failed (%d): %s\n  in %s\n  log %s\n%s"
                             % (r.returncode, cmd, cwd, log, tail))

    # ------------------------------------------------------------------ mesh
    def stage_mesh(self):
        cfg, d = self.cfg, self.d
        g, me = cfg["geometry"], cfg["mesh"]
        mdir = os.path.join(self.dir, "mesh")
        os.makedirs(mdir, exist_ok=True)
        sec = os.path.join(mdir, "section.med")
        if me.get("section_med"):
            shutil.copy(me["section_med"], sec)
        job = {
            "cylinder_radius": g["cylinder_radius_mm"], "rod_radius": g["rod_radius_mm"],
            "orifice_radius": g["orifice"]["radius_mm"],
            "orifice_pitch_radius": g["orifice"]["pitch_radius_mm"],
            "orifice_count": g["orifice"]["count"],
            "orifice_start_angle_deg": g["orifice"]["start_angle_deg"],
            "cylinder_length": g["cylinder_length_mm"], "piston_length": g["piston_length_mm"],
            "segments_piston": max(2, int(round(g["piston_length_mm"] / me["piston_cell_size_mm"]))),
            "gmsh": {"algorithm_2d": me["section"]["algorithm_2d"],
                     "subdiv_algo": me["section"]["subdiv_algo"],
                     "smoothing_steps": me["section"]["smoothing_steps"],
                     "size_factor": me["section"]["size_factor"],
                     "min_size": me["section"]["min_size_mm"],
                     "max_size": me["section"]["max_size_mm"],
                     "curvature_elements": me["section"].get("curvature_elements")},
            "preserve_orifice_area": bool(me["section"].get("preserve_orifice_area")),
            "merge_tolerance": me["merge_tolerance_mm"],
            "section_med": sec, "reuse_section": bool(me.get("section_med")),
            "report_path": os.path.join(mdir, "mesh_report.json"),
            "positions": [dict(p, unv_path=os.path.join(mdir, p["name"] + ".unv"))
                          for p in d["positions"]],
        }
        jpath = os.path.join(mdir, "mesh.json")
        with open(jpath, "w") as f:
            json.dump(job, f, indent=1)
        script = os.path.join(REPO, "src", "salome", "damper_mesh.py")
        log = os.path.join(self.logdir, "salome.log")
        t0 = time.time()
        lock = open(os.path.join(os.path.expanduser("~"), ".damper_salome.lock"), "w")
        fcntl.flock(lock, fcntl.LOCK_EX)  # SALOME sessions are serialised
        try:
            with open(log, "w") as f:
                r = subprocess.run([cfg["run"]["salome"], "-t", script, "args:" + jpath],
                                   cwd=mdir, stdout=f, stderr=subprocess.STDOUT)
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)
            lock.close()
        txt = open(log, errors="replace").read()
        if "DAMPER_MESH_OK" not in txt:
            raise StageError("SALOME mesh generation failed, see " + log)
        t_salome = time.time() - t0
        reports = json.load(open(job["report_path"]))
        info = {"salome_s": t_salome, "meshes": {}}
        V_tol = 5e-3
        for p, rep in zip(d["positions"], reports):
            case = os.path.join(mdir, p["name"])
            if os.path.exists(case):
                shutil.rmtree(case)
            F.mesh_system(case)
            self.foam("ideasUnvToFoam %s" % os.path.join(mdir, p["name"] + ".unv"), case,
                      os.path.join(self.logdir, "ideasUnvToFoam_%s.log" % p["name"]))
            self.foam("transformPoints -scale 0.001", case,
                      os.path.join(self.logdir, "transformPoints_%s.log" % p["name"]))
            bfile = os.path.join(case, "constant", "polyMesh", "boundary")
            btxt = open(bfile).read()
            names = re.findall(r"^\s{4}(\w+)\s*\n\s{4}\{", btxt, re.M)
            if sorted(names) != sorted(F.WALLS):
                raise StageError("unexpected patches in %s: %s" % (p["name"], names))
            btxt = re.sub(r"type\s+patch;", "type            wall;", btxt)
            btxt = re.sub(r"\s*physicalType\s+\w+;", "", btxt)
            open(bfile, "w").write(btxt)
            clog = os.path.join(self.logdir, "checkMesh_%s.log" % p["name"])
            self.foam("checkMesh -constant", case, clog)
            ctxt = open(clog).read()
            if "Mesh OK." not in ctxt:
                raise StageError("checkMesh failed for %s, see %s" % (p["name"], clog))
            V = num(ctxt, r"Total volume = ([0-9.eE+-]+)")
            xf = p["piston_front_x"]
            a = d["areas_mm2"]
            Va = (a["cylinder"] * xf + a["orifices"] * g["piston_length_mm"]
                  + a["annulus"] * (g["cylinder_length_mm"] - xf - g["piston_length_mm"])) * 1e-9
            if abs(V / Va - 1.0) > V_tol:
                raise StageError("volume check failed for %s: %g vs analytic %g" % (p["name"], V, Va))
            info["meshes"][p["name"]] = {
                "cells": int(re.search(r"cells:\s+(\d+)", ctxt).group(1)),
                "volume_m3": V, "volume_analytic_m3": Va,
                "max_skewness": num(ctxt, r"Max skewness = ([0-9.eE+-]+)"),
                "max_non_orthogonality": num(ctxt, r"non-orthogonality Max: ([0-9.eE+-]+)"),
                "salome": rep,
            }
            self.say("mesh %s: %d cells, V/V_analytic-1 = %.2e"
                     % (p["name"], info["meshes"][p["name"]]["cells"], V / Va - 1.0))
        return info

    # ------------------------------------------------------------------ solve
    def segment_dir(self, seg):
        return os.path.join(self.dir, "segments", seg["name"])

    def setup_segment(self, seg, prev):
        case = self.segment_dir(seg)
        if os.path.exists(case):
            shutil.rmtree(case)
        os.makedirs(os.path.join(case, "constant"))
        shutil.copytree(os.path.join(self.dir, "mesh", "M%d" % seg["mesh"], "constant", "polyMesh"),
                        os.path.join(case, "constant", "polyMesh"))
        F.constant(case, self.cfg)
        F.system(case, self.cfg, seg, self.cfg["run"]["n_procs"])
        F.fields(case, seg["start_name"], self.cfg, self.tab)
        open(os.path.join(case, "case.foam"), "w").close()
        if seg["init"] == "map":
            src = self.segment_dir(prev)
            if not os.path.isdir(os.path.join(src, seg["start_name"])):
                raise StageError("source time %s missing in %s" % (seg["start_name"], src))
            method = self.cfg["numerics"].get("map_method", "interpolate")
            if method == "interpolate":
                cmd = "mapFields %s -sourceTime %s -consistent -noFunctionObjects"
            elif method == "cellVolumeWeight":
                cmd = "mapFieldsPar %s -sourceTime %s -consistent -mapMethod cellVolumeWeight -noFunctionObjects"
            else:
                raise StageError("unknown numerics.map_method: %s" % method)
            self.foam(cmd % (src, seg["start_name"]), case,
                      os.path.join(self.logdir, "mapFields_%s.log" % seg["name"]))
            # mapFields sets aside target fields it cannot map (the point field
            # pointMotionUx, which only carries the piston-velocity table): restore them
            tdir = os.path.join(case, seg["start_name"])
            for f in glob.glob(os.path.join(tdir, "*.unmapped")):
                os.replace(f, f[:-len(".unmapped")])
            for f in ("mag(U)", "mag(U).gz"):
                if os.path.exists(os.path.join(tdir, f)):
                    os.remove(os.path.join(tdir, f))
            if not glob.glob(os.path.join(tdir, "pointMotionUx*")):
                raise StageError("pointMotionUx missing after mapping in " + tdir)
        return case

    def run_segment(self, seg, case):
        n = int(self.cfg["run"]["n_procs"])
        name = seg["name"]
        if n > 1:
            self.foam("decomposePar -force", case, os.path.join(self.logdir, "decomposePar_%s.log" % name))
            cmd = "%s --bind-to none -np %d rhoPimpleFoam -parallel" % (self.cfg["run"]["mpirun"], n)
        else:
            cmd = "rhoPimpleFoam"
        slog = os.path.join(case, "log.rhoPimpleFoam")
        t0 = time.time()
        self.foam(cmd, case, slog)
        wall = time.time() - t0
        txt = open(slog).read()
        if not re.search(r"^End\s*$", txt, re.M):
            raise StageError("solver did not reach End in %s" % name)
        times = re.findall(r"^Time = ([0-9.eE+-]+)", txt, re.M)
        if not times or abs(float(times[-1]) - seg["end"]) > 1e-9 * max(1.0, seg["end"]):
            raise StageError("segment %s ended at %s, expected %s" % (name, times[-1:], seg["end"]))
        if n > 1:
            self.foam("reconstructPar", case, os.path.join(self.logdir, "reconstructPar_%s.log" % name))
            if not self.cfg["run"]["keep_processor_dirs"]:
                for pd in glob.glob(os.path.join(case, "processor*")):
                    shutil.rmtree(pd)
        if not os.path.isdir(os.path.join(case, seg["end_name"])):
            raise StageError("end time directory %s missing in %s" % (seg["end_name"], name))
        return {"wall_s": wall, "steps": len(times)}

    def stage_solve(self):
        segs = self.d["schedule"]["segments"]
        done = self.state.setdefault("segments", {})
        prev = None
        for seg in segs:
            if done.get(seg["name"], {}).get("status") == "done":
                prev = seg
                continue
            case = self.setup_segment(seg, prev)
            self.say("solve %s (M%d, t = %s -> %s, %s)" % (seg["name"], seg["mesh"], seg["start_name"],
                                                           seg["end_name"], seg["init"]))
            info = self.run_segment(seg, case)
            info["status"] = "done"
            done[seg["name"]] = info
            self.save_state()
            self.say("  %s done: %d steps, %.0f s" % (seg["name"], info["steps"], info["wall_s"]))
            prev = seg
        return {"segments": done}

    # ------------------------------------------------------------------ post/render
    def stage_post(self):
        from . import post
        return post.collect(self)

    def stage_render(self):
        if not self.cfg["render"]["enabled"]:
            return {"skipped": True}
        from . import render
        return render.run(self)

    # ------------------------------------------------------------------ driver
    def run(self, start=None, only=None, force=False):
        lock = open(os.path.join(self.dir, ".lock"), "w")
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            raise StageError("another driver is already running in " + self.dir)
        lock.write(str(os.getpid()))
        lock.flush()
        order = list(STAGES)
        stamps = self.state.setdefault("stages", {})
        if only:
            order = [only]
        elif start:
            order = order[order.index(start):]
            for s in order:
                stamps.pop(s, None)
                if s == "solve":
                    self.state.pop("segments", None)
        for st in order:
            if not force and stamps.get(st, {}).get("status") == "done" and not only:
                continue
            t0 = time.time()
            self.say("stage %s" % st)
            info = getattr(self, "stage_" + st)()
            stamps[st] = {"status": "done", "wall_s": time.time() - t0,
                          "finished": time.strftime("%Y-%m-%d %H:%M:%S")}
            self.state.setdefault("info", {})[st] = info
            self.save_state()
        self.say("case complete")
