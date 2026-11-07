#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Parametric hexahedral mesh of a gas-filled piston-cylinder damper section.

Runs inside SALOME (``salome -t damper_mesh.py args:<mesh.json>``).  The input
JSON is written by the driver; it holds the cross-section and meshing
parameters and the list of piston positions (one sub-mesh M_k per entry).
Lengths are in millimetres.

The quad cross-section is generated ONCE (or read from a previously saved MED
file) and swept for every piston position, so all sub-meshes of a case share
an identical cross-section and differ only in their axial node distribution.
GMSH's quasi-structured quad algorithm is not bit-reproducible between
sessions, therefore the section is stored with the case for provenance.

Topology (x is the cylinder axis, 0 <= x <= L):
  chamber1  : full disk of radius rs,               0         <= x <= xf
  orifices  : N disks of radius ro on pitch radius,  xf        <= x <= xf + lp
  chamber2  : annulus rm < r < rs (rod side),        xf + lp   <= x <= L

The cross-section is one GEOM partition (outer disk, rod disk, N orifice
disks) meshed once with the GMSH quasi-structured quad algorithm and swept
along x, so all three blocks share conforming nodes.  Boundary and zone groups
are assigned from element positions, not from topological names, which makes
the generator valid for any orifice count N >= 1.
"""

import json
import math
import sys

import salome

salome.salome_init()

import GEOM  # noqa: E402,F401
import SMESH  # noqa: E402
from salome.geom import geomBuilder  # noqa: E402
from salome.smesh import smeshBuilder  # noqa: E402


def read_input():
    path = None
    for arg in sys.argv[1:]:
        path = arg.split("args:", 1)[1] if arg.startswith("args:") else arg
    if path is None:
        raise SystemExit("usage: salome -t damper_mesh.py args:<mesh.json>")
    with open(path) as f:
        return json.load(f)


def orifice_centres(c):
    n = c["orifice_count"]
    a0 = math.radians(c.get("orifice_start_angle_deg", 0.0))
    return [(c["orifice_pitch_radius"] * math.cos(a0 + 2.0 * math.pi * k / n),
             c["orifice_pitch_radius"] * math.sin(a0 + 2.0 * math.pi * k / n))
            for k in range(n)]


def build_section(geompy, c):
    """Return (partition, annulus_faces, rod_faces, orifice_faces)."""
    O = geompy.MakeVertex(0, 0, 0)
    Vx = geompy.MakeVectorDXDYDZ(1, 0, 0)
    outer = geompy.MakeDiskPntVecR(O, Vx, c["cylinder_radius"])
    tools = [geompy.MakeDiskPntVecR(O, Vx, c["rod_radius"])]
    for (y, z) in orifice_centres(c):
        tools.append(geompy.MakeDiskPntVecR(geompy.MakeVertex(0, y, z), Vx,
                                            c["orifice_radius"]))
    part = geompy.MakePartition([outer], tools, [], [],
                                geompy.ShapeType["FACE"], 0, [], 0)
    geompy.addToStudy(part, "section")

    faces = geompy.ExtractShapes(part, geompy.ShapeType["FACE"], True)
    a_rod = math.pi * c["rod_radius"] ** 2
    a_orf = math.pi * c["orifice_radius"] ** 2
    annulus, rod, orif = [], [], []
    for f in faces:
        area = geompy.BasicProperties(f)[1]
        _, y, z = geompy.PointCoordinates(geompy.MakeCDG(f))
        r = math.hypot(y, z)
        if abs(area - a_rod) < 1e-3 * a_rod and r < 1e-6 * c["cylinder_radius"]:
            rod.append(f)
        elif abs(area - a_orf) < 1e-3 * a_orf and \
                abs(r - c["orifice_pitch_radius"]) < 1e-3 * c["orifice_radius"]:
            orif.append(f)
        else:
            annulus.append(f)
    if len(rod) != 1 or len(orif) != c["orifice_count"] or len(annulus) != 1:
        raise SystemExit("section classification failed: annulus=%d rod=%d orifices=%d"
                         % (len(annulus), len(rod), len(orif)))
    groups = {}
    for name, flist in (("annulus", annulus), ("rodDisk", rod), ("orifice", orif)):
        g = geompy.CreateGroup(part, geompy.ShapeType["FACE"])
        geompy.UnionList(g, flist)
        geompy.addToStudyInFather(part, g, name)
        groups[name] = g
    return part, groups


def path_mesh(smesh, geompy, x0, x1, nseg, density=None, name="path"):
    edge = geompy.MakeLineTwoPnt(geompy.MakeVertex(x0, 0, 0),
                                 geompy.MakeVertex(x1, 0, 0))
    geompy.addToStudy(edge, name)
    m = smesh.Mesh(edge, name)
    algo = m.Segment()
    if density is None:
        algo.NumberOfSegments(int(nseg))
    else:
        hyp = algo.NumberOfSegments(int(nseg), None, [])
        hyp.SetConversionMode(1)
        hyp.SetTableFunction(list(density))
    if not m.Compute():
        raise SystemExit("1D path mesh failed: " + name)
    return m


def groups_by_name(mesh, groups):
    out = {}
    for g in groups:
        out[g.GetName()] = g
    return out


def section_mesh(smesh, geompy, c):
    """Quad mesh of the cross-section at x = 0 with groups annulus/rodDisk/orifice."""
    sec_file = c.get("section_med")
    if sec_file and c.get("reuse_section") and __import__("os").path.exists(sec_file):
        meshes, status = smesh.CreateMeshesFromMED(sec_file)
        mesh = meshes[0]
        print("section read from", sec_file)
    else:
        part, ggroups = build_section(geompy, c)
        mesh = smesh.Mesh(part, "section")
        g2d = mesh.Triangle(algo=smeshBuilder.GMSH_2D)
        gp = g2d.Parameters()
        gm = c["gmsh"]
        gp.SetMaxSize(float(gm.get("max_size", 1e22)))
        gp.Set2DAlgo(int(gm.get("algorithm_2d", 6)))
        gp.SetSubdivAlgo(int(gm.get("subdiv_algo", 0)))
        gp.SetSmouthSteps(int(gm.get("smoothing_steps", 25)))
        gp.SetSizeFactor(float(gm.get("size_factor", 0.8)))
        gp.SetMinSize(float(gm.get("min_size", 0.15)))
        if gm.get("curvature_elements"):
            # elements per 2*pi of a curved edge: refines the small orifice circles far more than the bore
            gp.SetMeshCurvatureSize(float(gm["curvature_elements"]))
        gp.SetIs2d(1)
        for name, g in ggroups.items():
            mesh.GroupOnGeom(g, name, SMESH.FACE)
        if not mesh.Compute():
            raise SystemExit("2D section mesh failed")
        if sec_file:
            mesh.ExportMED(sec_file)
    nonquad = mesh.NbTriangles() + mesh.NbPolygons()
    if nonquad:
        raise SystemExit("2D section mesh is not all-quad (%d non-quad faces)" % nonquad)
    names = sorted(g.GetName() for g in mesh.GetGroups())
    if names != ["annulus", "orifice", "rodDisk"]:
        raise SystemExit("unexpected section groups: %s" % names)
    return mesh


def orifice_area_ratio(mesh, c):
    """Area of the discretised orifice discs over N*pi*r_o^2 (< 1: the circles are inscribed polygons)."""
    g = [g for g in mesh.GetGroups() if g.GetName() == "orifice"][0]
    area = sum(mesh.GetArea(f) for f in g.GetIDs())
    return area / (c["orifice_count"] * math.pi * c["orifice_radius"] ** 2)


def preserve_orifice_area(mesh, c):
    """Scale the nodes of every orifice disc radially about its centre so that the polygonal disc has the area
    pi*r_o^2; the neighbouring annulus cells absorb the displacement of the circle nodes, (s - 1) r_o."""
    s = 1.0 / math.sqrt(orifice_area_ratio(mesh, c))
    centres = orifice_centres(c)
    g = [g for g in mesh.GetGroups() if g.GetName() == "orifice"][0]
    nodes = set()
    for f in g.GetIDs():
        nodes.update(mesh.GetElemNodes(f))
    for n in nodes:
        x, y, z = mesh.GetNodeXYZ(n)
        cy, cz = min(centres, key=lambda q: math.hypot(y - q[0], z - q[1]))
        mesh.MoveNode(n, x, cy + s * (y - cy), cz + s * (z - cz))
    return s


def volume_mesh(smesh, geompy, section, c, pos):
    """Sweep the section for one piston position and export it as UNV."""
    L = c["cylinder_length"]
    lp = c["piston_length"]
    xf = pos["piston_front_x"]
    xb = xf + lp
    rs, rm, ro = c["cylinder_radius"], c["rod_radius"], c["orifice_radius"]
    ro *= c.get("orifice_scale", 1.0)  # discs scaled by preserve_orifice_area
    if not (0.0 < xf and xb < L):
        raise SystemExit("piston outside cylinder: xf=%g xb=%g L=%g" % (xf, xb, L))
    tag = pos["name"]

    mesh = smesh.CopyMesh(section, tag, 1, 0)
    base = {g.GetName(): g for g in mesh.GetGroups()}
    n2d = mesh.NbFaces()

    p1 = path_mesh(smesh, geompy, 0.0, xf, pos["segments_chamber1"],
                   pos.get("density_chamber1", c.get("density_chamber1")), tag + "_path1")
    p2 = path_mesh(smesh, geompy, xf, xb, c["segments_piston"], None, tag + "_path2")
    p3 = path_mesh(smesh, geompy, xb, L, pos["segments_chamber2"],
                   pos.get("density_chamber2", c.get("density_chamber2")), tag + "_path3")

    # 1) chamber1: sweep the whole section from x=0 to x=xf
    g1, err = mesh.ExtrusionAlongPathObjects([], [], [mesh], p1, None, 1, 0, [], 0,
                                             0, [0, 0, 0], 1, [], 0)
    if err != SMESH.SMESH_MeshEditor.EXTR_OK:
        raise SystemExit("extrusion chamber1 failed: %s" % err)
    g1 = groups_by_name(mesh, g1)
    # 2) orifices: sweep the orifice tops through the piston
    g2, err = mesh.ExtrusionAlongPathObjects([], [], [g1["orifice_top"]], p2, None, 1, 0,
                                             [], 0, 0, [0, 0, 0], 1, [], 0)
    if err != SMESH.SMESH_MeshEditor.EXTR_OK:
        raise SystemExit("extrusion orifices failed: %s" % err)
    g2 = groups_by_name(mesh, g2)
    # 3) chamber2: copy the bottom annulus to x=xb (piston back face), make it
    #    conform with the orifice exits, and sweep both up to x=L
    tg = mesh.TranslateObject(base["annulus"], smesh.MakeDirStruct(xb, 0, 0), True, True)
    tg = groups_by_name(mesh, tg)
    back = tg[[k for k in tg if k.startswith("annulus")][0]]
    otop = g2[[k for k in g2 if k.endswith("_top")][0]]
    tol = c.get("merge_tolerance", 1e-6)
    coinc = mesh.FindCoincidentNodes(tol)
    if coinc:
        mesh.MergeNodes(coinc)
    g3, err = mesh.ExtrusionAlongPathObjects([], [], [back, otop], p3, None, 1, 0, [], 0,
                                             0, [0, 0, 0], 1, [], 0)
    if err != SMESH.SMESH_MeshEditor.EXTR_OK:
        raise SystemExit("extrusion chamber2 failed: %s" % err)
    g3 = groups_by_name(mesh, g3)
    if mesh.FindCoincidentNodes(tol):
        raise SystemExit("coincident nodes remain after the chamber2 sweep")

    # cell zones from the sweep groups
    vol = {"chamber1": [], "orifices": [], "chamber2": []}
    for key, gs in (("chamber1", g1), ("orifices", g2), ("chamber2", g3)):
        for name, g in gs.items():
            if name.endswith("_extruded"):
                vol[key] += list(g.GetIDs())

    # drop every helper group (keep the elements) and add the side walls
    for g in mesh.GetGroups():
        mesh.RemoveGroup(g)
    mesh.MakeBoundaryElements(SMESH.BND_2DFROM3D, "", "", 0, [])

    # classify all faces by position
    centres = orifice_centres(c)
    tolx = 1e-6 * L + 1e-9
    cls = {k: [] for k in ("bottomCap", "topCap", "pistonFront", "pistonBack",
                           "orificeInlet", "orificeOutlet", "cylinderWall",
                           "rodWall", "orificeWall")}

    def in_orifice(y, z):
        return any(math.hypot(y - cy, z - cz) < ro for (cy, cz) in centres)

    free = set(mesh.GetIdsFromFilter(smesh.GetFilter(SMESH.FACE, SMESH.FT_FreeFaces)))
    for fid in mesh.GetElementsByType(SMESH.FACE):
        x, y, z = mesh.BaryCenter(fid)
        internal = fid not in free
        if internal and not (abs(x - xf) < tolx or abs(x - xb) < tolx):
            raise SystemExit("unexpected internal face at x=%g r=%g" % (x, math.hypot(y, z)))
        r = math.hypot(y, z)
        if abs(x) < tolx:
            cls["bottomCap"].append(fid)
        elif abs(x - L) < tolx:
            cls["topCap"].append(fid)
        elif abs(x - xf) < tolx or abs(x - xb) < tolx:
            front = abs(x - xf) < tolx
            if in_orifice(y, z) != internal:
                raise SystemExit("piston face/orifice mismatch at x=%g y=%g z=%g" % (x, y, z))
            if internal:
                cls["orificeInlet" if front else "orificeOutlet"].append(fid)
            else:
                cls["pistonFront" if front else "pistonBack"].append(fid)
        elif xf < x < xb:
            cls["orificeWall"].append(fid)
        elif x > xb and r < 0.5 * (rm + rs):
            cls["rodWall"].append(fid)
        else:
            cls["cylinderWall"].append(fid)
    for name, ids in cls.items():
        if not ids:
            raise SystemExit("empty face group: " + name)
        mesh.MakeGroupByIds(name, SMESH.FACE, ids)
    for name, ids in vol.items():
        mesh.MakeGroupByIds(name, SMESH.VOLUME, ids)

    report = {
        "name": tag, "piston_front_x": xf,
        "nodes": mesh.NbNodes(), "hexahedra": mesh.NbHexas(),
        "volumes": mesh.NbVolumes(), "section_quads": n2d,
        "faces": {k: len(v) for k, v in cls.items()},
        "cells": {k: len(v) for k, v in vol.items()},
        "merged_nodes": sum(len(s) - 1 for s in coinc),
        "orifice_area_ratio": c["orifice_area_ratio"],
        "orifice_area_ratio_polygon": c["orifice_area_ratio_polygon"],
    }
    if report["volumes"] != report["hexahedra"]:
        raise SystemExit("non-hexahedral cells present: %s" % report)
    mesh.ExportUNV(pos["unv_path"])
    return report


def main():
    c = read_input()
    geompy = geomBuilder.New()
    smesh = smeshBuilder.New()
    section = section_mesh(smesh, geompy, c)
    c["orifice_area_ratio_polygon"] = c["orifice_area_ratio"] = orifice_area_ratio(section, c)
    if c.get("preserve_orifice_area"):
        c["orifice_scale"] = preserve_orifice_area(section, c)
        c["orifice_area_ratio"] = orifice_area_ratio(section, c)
        print("orifice discs scaled by %.5f: area ratio %.5f -> %.5f"
              % (c["orifice_scale"], c["orifice_area_ratio_polygon"], c["orifice_area_ratio"]))
    reports = [volume_mesh(smesh, geompy, section, c, pos) for pos in c["positions"]]
    with open(c["report_path"], "w") as f:
        json.dump(reports, f, indent=2)
    print("DAMPER_MESH_OK", json.dumps([(r["name"], r["hexahedra"]) for r in reports]))


main()
