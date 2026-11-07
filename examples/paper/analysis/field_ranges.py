# pvbatch: data ranges of |U| and p on the axial slice for the jobs of a render spec
# (spec format of src/render/render_slices.py, written by make_field_specs.py), used to choose fixed colour ranges.
# usage: pvbatch field_ranges.py spec.json
import argparse
import json

from paraview.simple import OpenFOAMReader, Slice

ap = argparse.ArgumentParser(description="Data ranges of |U| and p on the axial slice for the jobs of a render spec.")
ap.add_argument("spec", help="render spec (JSON) written by make_field_specs.py")
spec = json.load(open(ap.parse_args().spec))
L = spec["length_m"]
for job in spec["jobs"]:
    rd = OpenFOAMReader(FileName=job["case"])
    rd.MeshRegions = ["internalMesh"]
    rd.CellArrays = spec["fields"]
    rd.UpdatePipeline(job["time"])
    sl = Slice(Input=rd)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = [0.5 * L, 0.0, 0.0]
    sl.SliceType.Normal = spec["normal"]
    sl.UpdatePipeline(job["time"])
    pd = sl.PointData
    u = pd["U"].GetRange(-1)
    p = pd["p"].GetRange(0)
    print("RANGE %s t=%g |U| %.4g %.4g p %.6g %.6g" % (job["label"], job["time"], u[0], u[1], p[0], p[1]))
