# pvbatch: bottom-cap (x = 0) cross-section mesh of each case, surface with edges.
# usage: pvbatch render_sections.py out_dir case_dir1 [case_dir2 ...]
#        (reads the sub-mesh case_dir/mesh/M0; writes <out_dir>/<case name>.png)
import argparse
import os

from paraview.simple import (AssignViewToLayout, CreateLayout, CreateView, OpenFOAMReader,
                             Render, SaveScreenshot, Show, Hide)

ap = argparse.ArgumentParser(description="Bottom-cap cross-section mesh of each case.")
ap.add_argument("out_dir", help="directory of the images")
ap.add_argument("cases", nargs="+", metavar="case_dir", help="case directory")
args = ap.parse_args()
out = args.out_dir
os.makedirs(out, exist_ok=True)
view = CreateView("RenderView")
AssignViewToLayout(view=view, layout=CreateLayout("s"), hint=0)
view.ViewSize = [800, 800]
view.OrientationAxesVisibility = 0
view.Background = [1, 1, 1]
try:
    view.UseColorPaletteForBackground = 0
except AttributeError:
    pass
for case in args.cases:
    foam = os.path.join(case, "mesh", "M0", "case.foam")
    open(foam, "a").close()
    rd = OpenFOAMReader(FileName=foam)
    rd.MeshRegions = ["patch/bottomCap"]
    rd.CellArrays = []
    rd.UpdatePipeline(0)
    d = Show(rd, view)
    d.SetRepresentationType("Surface With Edges")
    d.AmbientColor = [0.85, 0.88, 0.95]
    d.DiffuseColor = [0.85, 0.88, 0.95]
    d.EdgeColor = [0.1, 0.1, 0.3]
    view.CameraParallelProjection = 1
    Render(view)          # the first render may reset the camera; set it afterwards
    view.CameraPosition = [-0.1, 0.0, 0.0]
    view.CameraFocalPoint = [0.0, 0.0, 0.0]
    view.CameraViewUp = [0.0, 0.0, 1.0]
    view.CameraParallelScale = 0.0104
    Render(view)
    name = os.path.basename(case.rstrip("/"))
    SaveScreenshot(os.path.join(out, name + ".png"), view, ImageResolution=[800, 800])
    print("wrote", name)
    Hide(rd, view)
