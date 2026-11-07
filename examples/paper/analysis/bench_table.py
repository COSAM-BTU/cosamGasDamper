#!/usr/bin/env python3
"""Rows of the benchmark results table (discharge coefficient of the single orifice per mesh and Reynolds number).

usage: bench_table.py [--tables DIR]

Reads orifice_benchmark.csv and, if present, orifice_benchmark_T2.csv (mesh-B runs repeated with Co_max = 0.5)
from the tables directory (default: tables), both written by orifice_benchmark.py collect; completed runs of the
second set replace the failed ones of the first set and are marked with superscript d; -- marks a failed,
unphysical or unsettled run (|drift| > 1 %) or one still running, an empty cell a run not attempted. The last row
is the correlation of Lichtarowicz et al., Eq. (lichtarowicz), solved with Re_L = Re_o / C_d,o at the nominal
Reynolds number of the column.
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figdata import cd_long_orifice  # noqa: E402

TABLES = "tables"           # directory of the derived CSV tables (set in main)
COLS = [("Re25", 25), ("Re100", 100), ("Re301", 301), ("Re1005", 1005), ("Re3014", 3014), ("Re5023", 5023),
        ("Re10046", 10046), ("Re30139", 30139), ("Re100463", 100463)]
ROWS = [("A", "A", "kOmegaSST"), ("X0.5", "X, $h_\\mathrm{f}=0.5$", "kOmegaSST"),
        ("X0.25", "X, $h_\\mathrm{f}=0.25$", "kOmegaSST"), ("X0.1", "X, $h_\\mathrm{f}=0.1$", "kOmegaSST"),
        ("X0.05", "X, $h_\\mathrm{f}=0.05$", "kOmegaSST"), ("X0.025", "X, $h_\\mathrm{f}=0.025$", "kOmegaSST"),
        ("Ax", "Ax", "kOmegaSST"), ("B", "B", "kOmegaSST"), ("B", "B, laminar", "laminar"), ("C", "C", "kOmegaSST")]


def load(name):
    """Rows of a benchmark CSV file of the tables directory by case name; empty if the file does not exist."""
    path = os.path.join(TABLES, name)
    if not os.path.exists(path):
        return {}
    return {r["case"]: r for r in csv.DictReader(open(path))}


def ok(r):
    """Completed run with a finite coefficient whose mean settled: |drift| (last minus third quarter) <= 1 %."""
    return (r is not None and r.get("complete") == "1" and r.get("Cd") not in ("", "nan")
            and abs(float(r.get("drift") or "nan")) <= 0.01)


def main(argv=None):
    global TABLES
    ap = argparse.ArgumentParser(description="Rows of the single-orifice benchmark table (LaTeX).")
    ap.add_argument("--tables", default="tables", help="directory of the benchmark CSV files (default: tables)")
    TABLES = ap.parse_args(argv).tables
    first, second = load("orifice_benchmark.csv"), load("orifice_benchmark_T2.csv")
    for mesh, label, turb in ROWS:
        cells = []
        for key, _ in COLS:
            name = "%s_%s_%s" % (mesh, "lam" if turb == "laminar" else "sst", key)
            r1, r2 = first.get(name), second.get(name)
            if ok(r1):
                cells.append("$%.3f$" % float(r1["Cd"]))
            elif ok(r2):
                cells.append("$%.3f$\\textsuperscript{d}" % float(r2["Cd"]))
            elif any(r is not None and r.get("complete") == "0" for r in (r1, r2)):
                cells.append("--")
            else:
                cells.append("")
        print(label + " & " + " & ".join(cells) + " \\\\")
    print("Eq.~\\eqref{eq:lichtarowicz} & " + " & ".join("$%.3f$" % cd_long_orifice(re) for _, re in COLS) + " \\\\")


if __name__ == "__main__":
    sys.exit(main())
