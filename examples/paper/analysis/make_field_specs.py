#!/usr/bin/env python3
"""Render specs (format of src/render/render_slices.py) for the field sequences of the paper:
one spec per orifice count N = 2, 6, 10, mean speed 2 m/s, t = T/4, T/2, 3T/4 and T of the compression stroke;
with --rebound one spec for N = 2 at t = 5T/4, 3T/2, 7T/4 and 2T of the rebound stroke.

usage: make_field_specs.py [--runs DIR] [--study NAME] [--out DIR] [--ranges FILE] [--rebound]

  --runs     directory holding the study directories                                     default: runs
  --study    study of the matrix cases N2_v2.0, N6_v2.0, N10_v2.0                         default: matrixC8ax
  --out      directory of the specs spec_N<N>.json (spec_R2.json with --rebound); the images go into its
             sub-folders N<N> (R2)                                                        default: field_specs
  --ranges   JSON file of fixed colour ranges per orifice count, e.g. {"2": {"U": [0, 150], "p": [1e5, 1.3e5]}};
             without it the specs carry no ranges (for field_ranges.py)

The case paths in the specs are formed from --runs as given; pass an absolute path if pvbatch runs elsewhere.
Then: pvbatch src/render/render_slices.py spec_N2.json
"""

import argparse
import json
import math
import os
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(description="Render specs for the field sequences of the paper.")
    ap.add_argument("--runs", default="runs", help="directory holding the study directories (default: runs)")
    ap.add_argument("--study", default="matrixC8ax", help="study of the matrix cases (default: matrixC8ax)")
    ap.add_argument("--out", default="field_specs", help="directory of the specs and images (default: field_specs)")
    ap.add_argument("--ranges", default=None, help="JSON file of fixed colour ranges per orifice count")
    ap.add_argument("--rebound", action="store_true", help="rebound stroke of N = 2 instead of the compression")
    args = ap.parse_args(argv)
    runs, out = args.runs, args.out
    ranges = json.load(open(args.ranges)) if args.ranges and args.ranges != "-" else {}
    rebound = args.rebound
    speed = "2.0"
    # the first four mesh segments end at T/4, T/2, 3T/4; T lies inside the fourth (M3)
    steps = [("t1", "seg01_M0", 0.25), ("t2", "seg02_M1", 0.5), ("t3", "seg03_M2", 0.75), ("t4", "seg04_M3", 1.0)]
    if rebound:
        # rebound stroke, N = 2 only: M3 spans 3T/4...3T/2, M2 3T/2...7T/4, M1 7T/4...2T
        steps = [("t5", "seg04_M3", 1.25), ("t6", "seg04_M3", 1.5), ("t7", "seg05_M2", 1.75), ("t8", "seg06_M1", 2.0)]
    os.makedirs(out, exist_ok=True)
    for n in ((2,) if rebound else (2, 6, 10)):
        case = os.path.join(runs, args.study, "N%d_v%s" % (n, speed))
        T = json.load(open(os.path.join(case, "config.resolved.json")))["_derived"]["stroke_time_s"]
        a = math.pi / 2.0          # first orifice on +z (as in the render stage of the framework)
        spec = {
            "jobs": [{"label": "N%d_%s" % (n, lab), "time": round(f * T, 12),
                      "case": os.path.join(case, "segments", seg, "case.foam")} for lab, seg, f in steps],
            "fields": ["U", "p"],
            "out_dir": os.path.join(out, ("R%d" if rebound else "N%d") % n),
            "resolution": [1800, 420],
            "normal": [0, -math.sin(a), math.cos(a)],
            "tangent": [0, math.cos(a), math.sin(a)],
            "radius_m": 0.01,
            "length_m": 0.12,
            "scalar_bar": False,
            "colormap": "Viridis (matplotlib)",
        }
        if str(n) in ranges:
            spec["ranges"] = ranges[str(n)]
        os.makedirs(spec["out_dir"], exist_ok=True)
        with open(os.path.join(out, ("spec_R%d.json" if rebound else "spec_N%d.json") % n), "w") as f:
            json.dump(spec, f, indent=1)
        print("spec N=%d, T=%g s" % (n, T))


if __name__ == "__main__":
    sys.exit(main())
