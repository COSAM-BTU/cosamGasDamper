# pvbatch: axial section (plane through the axis and the first orifice) of one sub-mesh,
# surface with edges; one full-length image and one close-up of the piston.
# usage: pvbatch render_side.py out_dir xmin xmax zmin zmax label=mesh_dir [label=mesh_dir ...]
#        (close-up window in metres; mesh_dir contains constant/polyMesh; writes <out_dir>/<label>_full.png and
#        <label>_zoom.png)
import argparse
import os

from paraview.simple import (AssignViewToLayout, CreateLayout, CreateView, Hide, OpenFOAMReader,
                             Render, SaveScreenshot, Show, Slice)

ap = argparse.ArgumentParser(description="Axial section of one sub-mesh per case, full length and close-up.")
ap.add_argument("out_dir", help="directory of the images")
for w in ("xmin", "xmax", "zmin", "zmax"):
    ap.add_argument(w, type=float, help="close-up window, m")
ap.add_argument("meshes", nargs="+", metavar="label=mesh_dir", help="image label and sub-mesh directory")
args = ap.parse_args()
out = args.out_dir
zx0, zx1, zz0, zz1 = args.xmin, args.xmax, args.zmin, args.zmax
os.makedirs(out, exist_ok=True)
L = 0.12
view = CreateView("RenderView")
AssignViewToLayout(view=view, layout=CreateLayout("s"), hint=0)
view.OrientationAxesVisibility = 0
view.Background = [1, 1, 1]
try:
    view.UseColorPaletteForBackground = 0
except AttributeError:
    pass


def shoot(name, focal, half_height, size):
    view.ViewSize = size
    view.CameraParallelProjection = 1
    Render(view)          # the first render may reset the camera; set it afterwards
    # look along +y: x to the right, z (first orifice) upwards, bottom cap on the left
    view.CameraFocalPoint = focal
    view.CameraPosition = [focal[0], -3 * L, focal[2]]
    view.CameraViewUp = [0.0, 0.0, 1.0]
    view.CameraParallelScale = half_height
    Render(view)
    SaveScreenshot(os.path.join(out, name), view, ImageResolution=size)
    print("wrote", name)


for arg in args.meshes:
    label, case = arg.split("=", 1)
    foam = os.path.join(case, "case.foam")
    open(foam, "a").close()
    rd = OpenFOAMReader(FileName=foam)
    rd.MeshRegions = ["internalMesh"]
    rd.CellArrays = []
    rd.UpdatePipeline(0)
    sl = Slice(Input=rd)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = [0.5 * L, 1e-7, 0.0]   # a hair off the plane of symmetry: no degenerate cuts
    sl.SliceType.Normal = [0.0, 1.0, 0.0]
    try:
        sl.Triangulatetheslice = 0                 # keep the cut polygons, i.e. the cell outlines
    except AttributeError:
        pass
    d = Show(sl, view)
    d.SetRepresentationType("Surface With Edges")
    d.AmbientColor = [0.85, 0.88, 0.95]
    d.DiffuseColor = [0.85, 0.88, 0.95]
    d.EdgeColor = [0.1, 0.1, 0.3]
    d.LineWidth = 1.0
    # full length: 3600 x 640 px, the 120 mm x 20 mm section fills the width
    shoot(label + "_full.png", [0.5 * L, 0.0, 0.0], 0.51 * L * 640 / 3600, [3600, 640])
    # close-up of the piston region
    w, h = zx1 - zx0, zz1 - zz0
    px = 2000
    shoot(label + "_zoom.png", [0.5 * (zx0 + zx1), 0.0, 0.5 * (zz0 + zz1)], 0.5 * h, [px, int(round(px * h / w))])
    Hide(sl, view)
