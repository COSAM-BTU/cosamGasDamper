# SPDX-License-Identifier: GPL-3.0-or-later
"""pvbatch script: slice through the cylinder axis and the first orifice, one PNG
per (time, field).  Usage: pvbatch render_slices.py render.json"""

import json
import os
import sys

from paraview.simple import (ColorBy, CreateView, GetColorTransferFunction, GetScalarBar,
                             Hide, OpenFOAMReader, Render, SaveScreenshot, Show, Slice)
try:
    from paraview.simple import AssignViewToLayout, CreateLayout
except ImportError:  # older ParaView
    AssignViewToLayout = CreateLayout = None

spec = json.load(open(sys.argv[1]))
view = CreateView("RenderView")
if CreateLayout is not None:
    AssignViewToLayout(view=view, layout=CreateLayout("render"), hint=0)
view.ViewSize = spec["resolution"]
view.OrientationAxesVisibility = 0
view.Background = [1, 1, 1]
try:
    view.UseColorPaletteForBackground = 0
except AttributeError:
    pass
n = spec["normal"]
L = spec["length_m"]
labels = {"p": "p [Pa]", "U": "|U| [m/s]", "T": "T [K]"}

t = spec["tangent"]                     # in-plane direction to the first orifice
R = spec["radius_m"]
W, H = spec["resolution"]
ranges = spec.get("ranges", {})

for job in spec["jobs"]:
    rd = OpenFOAMReader(FileName=job["case"])
    rd.MeshRegions = ["internalMesh"]
    rd.CellArrays = spec["fields"]
    try:
        rd.CaseType = "Reconstructed Case"
    except AttributeError:
        pass
    rd.UpdatePipeline(job["time"])
    sl = Slice(Input=rd)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = [0.5 * L, 0.0, 0.0]
    sl.SliceType.Normal = n
    view.ViewTime = job["time"]
    disp = Show(sl, view)
    prev = None
    for fld in spec["fields"]:
        if prev is not None:
            GetScalarBar(prev, view).Visibility = 0
        if fld == "U":
            ColorBy(disp, ("POINTS", "U", "Magnitude"))
        else:
            ColorBy(disp, ("POINTS", fld))
        lut = GetColorTransferFunction(fld)
        lut.ApplyPreset(spec.get("colormap", "Viridis (matplotlib)"), True)
        if fld in ranges:
            lut.RescaleTransferFunction(*ranges[fld])
        else:
            disp.RescaleTransferFunctionToDataRange(True, False)
        disp.SetScalarBarVisibility(view, bool(spec.get("scalar_bar", True)))
        bar = GetScalarBar(lut, view)
        bar.Title = labels.get(fld, fld)
        bar.ComponentTitle = ""
        bar.TitleColor = [0, 0, 0]
        bar.LabelColor = [0, 0, 0]
        bar.TitleFontSize = 20
        bar.LabelFontSize = 18
        bar.Orientation = "Horizontal"
        bar.WindowLocation = "Any Location"
        bar.Position = [0.30, 0.04]
        bar.ScalarBarLength = 0.4
        bar.Visibility = 1 if spec.get("scalar_bar", True) else 0
        # camera: look along the slice normal, x to the right, first orifice upwards;
        # render once first, because the first render may reset the camera
        Render(view)
        view.CameraFocalPoint = [0.5 * L, 0.0, 0.0]
        view.CameraPosition = [0.5 * L, n[1] * 3 * L, n[2] * 3 * L]
        view.CameraViewUp = list(t)
        view.CameraParallelProjection = 1
        view.CameraParallelScale = max(1.25 * R, 0.55 * L * H / W)
        Render(view)
        out = os.path.join(spec["out_dir"], "%s_%s.png" % (fld, job["label"]))
        SaveScreenshot(out, view, ImageResolution=spec["resolution"])
        print("wrote", out)
        prev = lut
    Hide(sl, view)
