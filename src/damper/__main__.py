# SPDX-License-Identifier: GPL-3.0-or-later
"""Command line interface.

  run    <case_dir> [--config cfg.yaml] [--from STAGE | --only STAGE]
  sweep  <cfg.yaml> --root DIR [--max-cores N]   expand `sweep:` and run all cases
  check  <cfg.yaml>                              validate a configuration and print the schedule
"""

import argparse
import itertools
import json
import os
import subprocess
import sys
import time

import yaml

from . import config as C
from .pipeline import STAGES, Case


def cmd_run(a):
    case = Case(a.case_dir, a.config)
    case.run(start=a.start, only=a.only, force=a.force)


def cmd_check(a):
    cfg = C.load(a.config)
    d = C.derive(cfg)
    print(json.dumps({k: d[k] for k in ("stroke_time_s", "mean_speed_mps", "x0_mm", "grid_mm",
                                        "end_time_s", "write_interval_s", "orifice_area_ratio")},
                     indent=1))
    for s in d["schedule"]["segments"]:
        print("%-10s M%d  %s -> %s  %s" % (s["name"], s["mesh"], s["start_name"], s["end_name"],
                                           s["init"]))


def resolve_paths(cfg):
    """A relative mesh.section_med is relative to the directory the sweep is started from; make it absolute,
    because every case runs in its own process with a different working directory."""
    sec = cfg.get("mesh", {}).get("section_med")
    if sec and not os.path.isabs(sec):
        cfg["mesh"]["section_med"] = os.path.abspath(sec)
    return cfg


def expand_sweep(path):
    raw = yaml.safe_load(open(path)) or {}
    variants = raw.pop("variants", None)
    if variants:  # named cases, each a set of dotted-key overrides
        cases = []
        for v in variants:
            cfg = json.loads(json.dumps(raw))
            for k, val in v.get("set", {}).items():
                C.set_dotted(cfg, k, val)
            cfg.setdefault("case", {})["name"] = v["name"]
            cfg["case"].pop("name_template", None)
            cases.append(resolve_paths(cfg))
        return cases
    sweep = raw.pop("sweep", None) or {}
    keys = list(sweep)
    cases = []
    for values in itertools.product(*[sweep[k] for k in keys]):
        cfg = json.loads(json.dumps(raw))
        tags = []
        for k, v in zip(keys, values):
            C.set_dotted(cfg, k, v)
            tags.append("%s%s" % (k.split(".")[-1][:1].upper() if k.split(".")[-1] else k,
                                  ("%g" % v).replace(".", "p") if isinstance(v, float) else v))
        name = cfg.get("case", {}).get("name", "case")
        tmpl = raw.get("case", {}).get("name_template")
        if tmpl:
            flat = {k.replace(".", "_"): v for k, v in zip(keys, values)}
            name = tmpl.format(**flat)
        else:
            name = name + "_" + "_".join(str(t) for t in tags)
        cfg.setdefault("case", {})["name"] = name
        sec = cfg.get("mesh", {}).get("section_med")
        if sec and "{name}" in sec:  # reuse the stored cross-section of the case of the same name
            cfg["mesh"]["section_med"] = sec.format(name=name)
        cases.append(resolve_paths(cfg))
    return cases


def cmd_sweep(a):
    cases = []
    for cfgfile in a.config:  # several study files share one queue and core budget
        cases += expand_sweep(cfgfile)
    names = [c["case"]["name"] for c in cases]
    if len(set(names)) != len(names):
        raise SystemExit("duplicate case names in the sweep: %s" % names)
    root = os.path.abspath(a.root)
    os.makedirs(root, exist_ok=True)
    jobs = []
    for cfg in cases:
        C.derive(C.deep_merge(C.load({}), cfg))  # validate before launching anything
        cdir = os.path.join(root, cfg["case"]["name"])
        os.makedirs(cdir, exist_ok=True)
        with open(os.path.join(cdir, "config.yaml"), "w") as f:
            yaml.safe_dump(cfg, f, sort_keys=False)
        n = int(C.load(cfg)["run"]["n_procs"])
        jobs.append((cdir, n))
    print("sweep: %d cases under %s, core budget %d" % (len(jobs), root, a.max_cores), flush=True)
    if a.dry_run:
        for cdir, n in jobs:
            print("  ", os.path.basename(cdir), n)
        return
    running, used, failed = [], 0, []
    queue = list(jobs)
    while queue or running:
        for p in list(running):
            if p[0].poll() is not None:
                running.remove(p)
                used -= p[2]
                status = "ok" if p[0].returncode == 0 else "FAILED(%d)" % p[0].returncode
                if p[0].returncode != 0:
                    failed.append(p[1])
                print("[%s] finished %s %s" % (time.strftime("%H:%M:%S"), os.path.basename(p[1]),
                                               status), flush=True)
        while queue and used + queue[0][1] <= a.max_cores:
            cdir, n = queue.pop(0)
            log = open(os.path.join(cdir, "driver.log"), "a")
            p = subprocess.Popen([sys.executable, "-m", "damper", "run", cdir],
                                 stdout=log, stderr=subprocess.STDOUT,
                                 cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            running.append((p, cdir, n))
            used += n
            print("[%s] started  %s (%d cores, %d in use)" % (time.strftime("%H:%M:%S"),
                                                               os.path.basename(cdir), n, used),
                  flush=True)
        time.sleep(5)
    print("sweep done: %d ok, %d failed %s" % (len(jobs) - len(failed), len(failed),
                                               [os.path.basename(f) for f in failed]), flush=True)


def cmd_collect(a):
    """One row per case from results/summary.json under the given roots."""
    import csv
    import glob
    rows = []
    for root in a.roots:
        for f in sorted(glob.glob(os.path.join(root, "*", "results", "summary.json"))):
            s = json.load(open(f))
            s["case_dir"] = os.path.basename(os.path.dirname(os.path.dirname(f)))
            s.pop("mass_jumps_rel", None)
            rows.append(s)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    w = csv.DictWriter(open(a.out, "w", newline="") if a.out else sys.stdout, fieldnames=keys)
    w.writeheader()
    for r in rows:
        w.writerow(r)


def main():
    ap = argparse.ArgumentParser(prog="damper")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("case_dir")
    r.add_argument("--config")
    g = r.add_mutually_exclusive_group()
    g.add_argument("--from", dest="start", choices=STAGES)
    g.add_argument("--only", choices=STAGES)
    r.add_argument("--force", action="store_true")
    r.set_defaults(func=cmd_run)
    s = sub.add_parser("sweep")
    s.add_argument("config", nargs="+")
    s.add_argument("--root", required=True)
    s.add_argument("--max-cores", type=int, default=8)
    s.add_argument("--dry-run", action="store_true")
    s.set_defaults(func=cmd_sweep)
    c = sub.add_parser("check")
    c.add_argument("config")
    c.set_defaults(func=cmd_check)
    k = sub.add_parser("collect")
    k.add_argument("roots", nargs="+")
    k.add_argument("--out")
    k.set_defaults(func=cmd_collect)
    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
