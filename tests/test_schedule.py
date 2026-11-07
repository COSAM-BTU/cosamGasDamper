# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests of the configuration, the mesh-segment schedule and the velocity table.

Run with:  PYTHONPATH=src python3 -m pytest tests   (or: PYTHONPATH=src python3 tests/test_schedule.py)
"""

import math

from damper import config as C


def make(**motion):
    cfg = C.load({"motion": motion})
    d = C.derive(cfg)
    cfg["_derived"] = d
    return cfg, d


def test_constant_schedule_k4():
    cfg, d = make(mean_speed_mps=1.0)
    segs = d["schedule"]["segments"]
    assert [s["mesh"] for s in segs] == [0, 1, 2, 3, 2, 1]
    T = d["stroke_time_s"]
    assert math.isclose(segs[3]["start"], 0.75 * T) and math.isclose(segs[3]["end"], 1.5 * T)
    assert segs[0]["init"] == "initial" and all(s["init"] == "map" for s in segs[1:])


def test_sub_mesh_deformation_bounded():
    """No sub-mesh is deformed by more than one grid spacing S/K."""
    for prof in ("constant", "sinusoidal"):
        cfg, d = make(mean_speed_mps=0.5, profile=prof, cycles=2)
        S, K, T, x0 = cfg["motion"]["stroke_mm"], cfg["motion"]["meshes_per_stroke"], d["stroke_time_s"], d["x0_mm"]
        for s in d["schedule"]["segments"]:
            xm = d["grid_mm"][s["mesh"]]
            for i in range(201):
                t = s["start"] + (s["end"] - s["start"]) * i / 200
                x = C.piston_x(prof, t, T, S, x0)
                assert abs(x - xm) <= S / K + 1e-9, (prof, s["name"], t, x, xm)


def test_sinusoidal_boundaries():
    cfg, d = make(mean_speed_mps=1.0, profile="sinusoidal")
    T = d["stroke_time_s"]
    starts = [s["start"] for s in d["schedule"]["segments"]]
    for got, want in zip(starts[1:4], (T / 3, T / 2, 2 * T / 3)):
        assert math.isclose(got, want, rel_tol=1e-8)


def test_velocity_table_left_continuous():
    cfg, d = make(mean_speed_mps=1.0, cycles=2)
    tab = dict(C.velocity_table(cfg, d))
    T = d["stroke_time_s"]
    assert tab[T] < 0 and tab[2 * T] > 0   # value AT the reversal is the old one


def test_displacement_integrates_to_stroke():
    for prof in ("constant", "sinusoidal"):
        cfg, d = make(mean_speed_mps=0.4, profile=prof)
        T, S = d["stroke_time_s"], cfg["motion"]["stroke_mm"]
        n = 20000
        disp = sum(-C.piston_u(prof, (i + 0.5) * T / n, T, S) * T / n for i in range(n)) * 1e3
        assert math.isclose(disp, S, rel_tol=1e-6)


def test_time_names_exact():
    cfg, d = make(stroke_time_s=0.72)
    for s in d["schedule"]["segments"]:
        assert float(s["end_name"]) == s["end"]


def test_invalid_geometry_rejected():
    try:
        C.derive(C.load({"geometry": {"orifice": {"count": 40}}}))
    except C.ConfigError:
        return
    raise AssertionError("overlapping orifices were accepted")


def test_invalid_map_method_rejected():
    try:
        C.derive(C.load({"numerics": {"map_method": "nearest"}}))
    except C.ConfigError:
        return
    raise AssertionError("an unknown mapping method was accepted")


def test_orifice_options_default_off():
    cfg = C.load({})
    assert cfg["mesh"]["section"]["preserve_orifice_area"] is False
    assert cfg["numerics"]["map_method"] == "interpolate"


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
