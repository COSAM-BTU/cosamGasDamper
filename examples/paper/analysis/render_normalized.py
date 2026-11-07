# pvbatch: velocity magnitude normalised by the ideal orifice velocity, |u|/u_o, at mid-compression
# (t = T/2, end of mesh segment 2) in the plane through the axis and the first orifice; one PNG per case,
# no scalar bar, same camera as src/render/render_slices.py.
# usage: pvbatch render_normalized.py out_dir vmax case_dir [case_dir ...]
#        vmax = 0 only prints the data range of |u|/u_o on the slice; otherwise writes <out_dir>/uStar_<case>.png
#        with the colour range 0 ... vmax
import argparse
import json
import math
import os

from paraview.simple import (AssignViewToLayout, Calculator, ColorBy, CreateLayout, CreateView,
                             GetColorTransferFunction, Hide, OpenFOAMReader, Render, SaveScreenshot,
                             Show, Slice)

ap = argparse.ArgumentParser(description="Velocity magnitude over the ideal orifice velocity at mid-compression.")
ap.add_argument("out_dir", help="directory of the images")
ap.add_argument("vmax", type=float, help="upper end of the colour range; 0: only print the data ranges")
ap.add_argument("cases", nargs="+", metavar="case_dir", help="case directory")
args = ap.parse_args()
out, vmax = args.out_dir, args.vmax
os.makedirs(out, exist_ok=True)
L, R = 0.12, 0.01
W, H = 1800, 420
a = math.pi / 2.0                       # first orifice on +z
n = [0.0, -math.sin(a), math.cos(a)]
tangent = [0.0, math.cos(a), math.sin(a)]
view = CreateView("RenderView")
AssignViewToLayout(view=view, layout=CreateLayout("n"), hint=0)
view.ViewSize = [W, H]
view.OrientationAxesVisibility = 0
view.Background = [1, 1, 1]
try:
    view.UseColorPaletteForBackground = 0
except AttributeError:
    pass

for case in args.cases:
    cfg = json.load(open(os.path.join(case, "config.resolved.json")))
    d = cfg["_derived"]
    # ideal orifice velocity on the discrete orifice area (polygon of n boundary segments per orifice)
    rep = json.load(open(os.path.join(case, "mesh", "mesh_report.json")))[0]
    nseg = rep["faces"]["orificeWall"] / (cfg["geometry"]["orifice"]["count"] * max(
        2, int(round(cfg["geometry"]["piston_length_mm"] / cfg["mesh"]["piston_cell_size_mm"]))))
    a_h = rep.get("orifice_area_ratio", nseg / (2 * math.pi) * math.sin(2 * math.pi / nseg))
    u_o = d["mean_speed_mps"] * d["orifice_area_ratio"] / a_h
    t = 0.5 * d["stroke_time_s"]
    foam = os.path.join(case, "segments", "seg02_M1", "case.foam")
    rd = OpenFOAMReader(FileName=foam)
    rd.MeshRegions = ["internalMesh"]
    rd.CellArrays = ["U"]
    rd.UpdatePipeline(t)
    sl = Slice(Input=rd)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = [0.5 * L, 0.0, 0.0]
    sl.SliceType.Normal = n
    calc = Calculator(Input=sl)
    calc.ResultArrayName = "uStar"
    calc.Function = "mag(U)/%.10g" % u_o
    calc.UpdatePipeline(t)
    rng = calc.PointData["uStar"].GetRange(0)
    label = os.path.basename(case.rstrip("/"))
    print("RANGE %s t=%g u_o=%.4g uStar %.4g %.4g" % (label, t, u_o, rng[0], rng[1]))
    if vmax <= 0:
        continue
    view.ViewTime = t
    disp = Show(calc, view)
    ColorBy(disp, ("POINTS", "uStar"))
    lut = GetColorTransferFunction("uStar")
    lut.ApplyPreset("Viridis (matplotlib)", True)
    lut.RescaleTransferFunction(0.0, vmax)
    disp.SetScalarBarVisibility(view, False)
    Render(view)          # the first render may reset the camera; set it afterwards
    view.CameraFocalPoint = [0.5 * L, 0.0, 0.0]
    view.CameraPosition = [0.5 * L, n[1] * 3 * L, n[2] * 3 * L]
    view.CameraViewUp = tangent
    view.CameraParallelProjection = 1
    view.CameraParallelScale = max(1.25 * R, 0.55 * L * H / W)
    Render(view)
    SaveScreenshot(os.path.join(out, "uStar_%s.png" % label), view, ImageResolution=[W, H])
    print("wrote", label)
    Hide(calc, view)
