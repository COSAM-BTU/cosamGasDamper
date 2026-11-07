# pvbatch: perspective view of the boundary mesh of one sub-mesh, with the near half of the cylinder wall and
# of the caps cut away, so that the piston faces, the orifice walls and the rod are visible.
# usage: pvbatch render_3d.py out_dir label=mesh_dir [label=mesh_dir ...]   (mesh_dir contains constant/polyMesh)
# writes <out_dir>/<label>_3d.png with a transparent background (compose_view3d.py arranges them)
import argparse
import os

from paraview.simple import (AssignViewToLayout, Clip, CreateLayout, CreateView, Hide, OpenFOAMReader, Render,
                             SaveScreenshot, Show)

ap = argparse.ArgumentParser(description="Perspective view of the boundary mesh of one sub-mesh per case.")
ap.add_argument("out_dir", help="directory of the images")
ap.add_argument("meshes", nargs="+", metavar="label=mesh_dir", help="image label and sub-mesh directory")
args = ap.parse_args()
out = args.out_dir
os.makedirs(out, exist_ok=True)
L = 0.12
W, H = 1600, 1000
view = CreateView("RenderView")
AssignViewToLayout(view=view, layout=CreateLayout("p"), hint=0)
view.ViewSize = [W, H]
view.OrientationAxesVisibility = 0
view.Background = [1, 1, 1]
try:
    view.UseColorPaletteForBackground = 0
except AttributeError:
    pass

DARK = [0.1, 0.1, 0.3]
GROUPS = [  # patches, colour, edge colour, cut away the near half; on the finely meshed piston faces and orifice
    # walls the edges take a darker shade of the face colour, so that the colour stays visible
    (["cylinderWall", "bottomCap", "topCap"], [0.85, 0.88, 0.95], DARK, True),
    (["pistonFront", "pistonBack"], [0.75, 0.75, 0.75], [0.45, 0.45, 0.45], False),
    (["orificeWall"], [0.90, 0.55, 0.45], [0.60, 0.30, 0.22], False),
    (["rodWall"], [0.60, 0.60, 0.60], DARK, False),
]

for arg in args.meshes:
    label, case = arg.split("=", 1)
    foam = os.path.join(case, "case.foam")
    open(foam, "a").close()
    shown = []
    for patches, colour, edge, cut in GROUPS:
        rd = OpenFOAMReader(FileName=foam)
        rd.MeshRegions = ["patch/" + p for p in patches]
        rd.CellArrays = []
        rd.UpdatePipeline(0)
        src = rd
        if cut:
            src = Clip(Input=rd)
            src.ClipType = "Plane"
            src.ClipType.Origin = [0.5 * L, 0.0, 0.0]
            src.ClipType.Normal = [0.0, -1.0, 0.0]     # keep y > 0, the far half seen from the camera
        d = Show(src, view)
        d.SetRepresentationType("Surface With Edges")
        d.AmbientColor = colour
        d.DiffuseColor = colour
        d.EdgeColor = edge
        shown.append(src)
    Render(view)          # the first render may reset the camera; set it afterwards
    view.CameraParallelProjection = 0
    view.CameraFocalPoint = [0.066, 0.0, 0.0]
    view.CameraPosition = [0.004, -0.112, 0.056]
    view.CameraViewUp = [0.0, 0.0, 1.0]
    view.CameraViewAngle = 30
    Render(view)
    SaveScreenshot(os.path.join(out, label + "_3d.png"), view, ImageResolution=[W, H], TransparentBackground=1)
    print("wrote", label)
    for s in shown:
        Hide(s, view)
