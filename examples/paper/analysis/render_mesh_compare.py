# pvbatch: cross-section mesh (bottom cap, x = 0) of several cases, the whole section and a close-up of the
# first orifice (on the z axis at the pitch radius), for a side-by-side comparison of mesh resolutions.
# usage: pvbatch render_mesh_compare.py out_dir label=case_dir [label=case_dir ...]
#        (reads the sub-mesh case_dir/mesh/M0; writes <out_dir>/<label>_full.png and <label>_zoom.png)
import argparse
import os

from paraview.simple import (AssignViewToLayout, CreateLayout, CreateView, Hide, OpenFOAMReader, Render,
                             SaveScreenshot, Show)

PITCH = 0.007                       # m, orifice pitch radius (first orifice at 90 deg, i.e. on +z)
ap = argparse.ArgumentParser(description="Cross-section mesh at the bottom cap, whole section and first orifice.")
ap.add_argument("out_dir", help="directory of the images")
ap.add_argument("cases", nargs="+", metavar="label=case_dir", help="image label and case directory")
args = ap.parse_args()
out = args.out_dir
os.makedirs(out, exist_ok=True)
view = CreateView("RenderView")
AssignViewToLayout(view=view, layout=CreateLayout("s"), hint=0)
view.ViewSize = [900, 900]
view.OrientationAxesVisibility = 0
view.Background = [1, 1, 1]
try:
    view.UseColorPaletteForBackground = 0
except AttributeError:
    pass
for arg in args.cases:
    label, case = arg.split("=", 1)
    foam = os.path.join(case, "mesh", "M0", "case.foam")
    open(foam, "a").close()
    rd = OpenFOAMReader(FileName=foam)
    rd.MeshRegions = ["patch/bottomCap"]
    rd.CellArrays = []
    rd.UpdatePipeline(0)
    d = Show(rd, view)
    d.SetRepresentationType("Surface With Edges")
    d.AmbientColor = d.DiffuseColor = [0.85, 0.88, 0.95]
    d.EdgeColor = [0.1, 0.1, 0.3]
    view.CameraParallelProjection = 1
    Render(view)
    for tag, focal, scale in (("full", [0.0, 0.0, 0.0], 0.0104), ("zoom", [0.0, 0.0, PITCH], 0.0016)):
        view.CameraFocalPoint = focal
        view.CameraPosition = [focal[0] - 0.1, focal[1], focal[2]]
        view.CameraViewUp = [0.0, 0.0, 1.0]
        view.CameraParallelScale = scale
        Render(view)
        SaveScreenshot(os.path.join(out, "%s_%s.png" % (label, tag)), view, ImageResolution=[900, 900])
    print("wrote", label)
    Hide(rd, view)
