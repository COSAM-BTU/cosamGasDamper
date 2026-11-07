#!/usr/bin/env python3
"""Flow at a constant flow rate through a single sharp-edged orifice of length l = 4d between two plenums,
for comparison of the discharge coefficient with the long-orifice data of Lichtarowicz et al. (1965).

Same orifice as the damper (d = 2 mm, l = 8 mm), same solver (rhoPimpleFoam), time scheme, Courant limit,
gas model, turbulence model, wall functions and convection schemes as the damper cases. The run starts from
rest, ramps the inflow up over 5% of the run and lasts eight transit times of the domain at the mean orifice
velocity; the coefficient is the time average over the second half of the run (mode "transient", default).
Mode "steady" uses local time stepping instead. Two-dimensional axisymmetric wedge meshes written with blockMesh:

  A    in-plane cell 0.14 mm in the orifice, axial cell 1 mm (in-plane resolution of the damper reference mesh,
       with uniform axial cells of the chamber size)
  Ar   as A, but radially graded to 0.005 mm at the orifice wall
  Ax   as A, but axially graded to 0.01 mm at the orifice inlet edge
  B    radially and axially graded (0.005 mm at the wall, 0.01 mm at the inlet edge)
  C    B with all cell sizes halved
  Xh   as A, but uniform axial cells of h mm in the orifice, graded to 1 mm away from the plate
       (X0.1: axial cell size of the damper reference mesh in the orifice and at its faces)

usage: orifice_benchmark.py setup ROOT [transient|steady]   write all cases
       orifice_benchmark.py collect ROOT                    print the discharge coefficients (CSV)
"""
import math
import os
import sys

D_O, L_O, R_S = 2.0, 8.0, 10.0          # mm: orifice diameter and length, plenum radius
L_UP, L_DN = 20.0, 60.0                 # mm: plenum lengths upstream and downstream of the orifice
T0, M, CP, AS, TS = 300.0, 28.0134, 1040.0, 1.401e-6, 107.0
R_GAS = 8314.462618 / M
WEDGE = 2.5                             # deg, half angle

# (p0 [Pa], ideal orifice velocity [m/s]) for orifice Reynolds numbers 25 ... 1e5
POINTS = [(1e5, 0.2), (1e5, 0.8), (1e5, 2.4), (1e5, 8.0), (1e5, 24.0), (1e5, 40.0),
          (2e6, 4.0), (2e6, 12.0), (2e6, 40.0)]

# sections per direction: (length fraction of the block, first cell, last cell) in mm
MESHES = {
    "A": {"rc": [(1.0, 1 / 7, 1 / 7)], "ro": [(1.0, 0.15, 1.0)], "xu": [(1.0, 1.0, 1.0)],
          "xo": [(1.0, 1.0, 1.0)], "xd": [(1.0, 1.0, 2.0)]},
    "Ar": {"rc": [(1.0, 0.05, 0.005)], "ro": [(1.0, 0.005, 1.0)], "xu": [(1.0, 1.0, 1.0)],
           "xo": [(1.0, 1.0, 1.0)], "xd": [(1.0, 1.0, 2.0)]},
    "Ax": {"rc": [(1.0, 1 / 7, 1 / 7)], "ro": [(1.0, 0.15, 1.0)], "xu": [(1.0, 1.0, 0.01)],
           "xo": [(0.125, 0.01, 0.1), (0.75, 0.1, 0.1), (0.125, 0.1, 0.02)], "xd": [(1.0, 0.02, 2.0)]},
    "B": {"rc": [(1.0, 0.05, 0.005)], "ro": [(1.0, 0.005, 1.0)], "xu": [(1.0, 1.0, 0.01)],
          "xo": [(0.125, 0.01, 0.1), (0.75, 0.1, 0.1), (0.125, 0.1, 0.02)], "xd": [(1.0, 0.02, 2.0)]},
}
MESHES["C"] = {k: [(f, a / 2, b / 2) for f, a, b in v] for k, v in MESHES["B"].items()}
# axial resolution of the orifice: in-plane as A, uniform axial cell h in the orifice, graded from h at the plate
# faces to 1 mm upstream and 2 mm downstream
for _h in (0.5, 0.25, 0.1, 0.05, 0.025):
    MESHES["X%g" % _h] = {"rc": [(1.0, 1 / 7, 1 / 7)], "ro": [(1.0, 0.15, 1.0)], "xu": [(1.0, 1.0, _h)],
                          "xo": [(1.0, _h, _h)], "xd": [(1.0, _h, 2.0)]}
X_POINTS = [(1e5, 8.0), (1e5, 40.0), (2e6, 40.0)]
TURB = {"A": ["kOmegaSST"], "Ar": ["kOmegaSST"], "Ax": ["kOmegaSST"], "B": ["kOmegaSST", "laminar"],
        "C": ["kOmegaSST"]}
TURB.update({m: ["kOmegaSST"] for m in MESHES if m.startswith("X")})
LAMINAR_MAX_RE = 1500.0
RAMP = 1000                             # steady mode: pseudo-time iterations of the inflow ramp
FLOWS = 8.0                             # transient mode: run time in domain transit times L/u_o
MODE = "transient"


def mu(T):
    return AS * T ** 1.5 / (T + TS)


def section(L, a, b):
    """Number of cells and expansion ratio (last/first) of a geometric distribution over L."""
    if abs(a - b) < 1e-12 * max(a, b):
        return max(1, int(round(L / a))), 1.0
    q = (L - a) / (L - b)
    n = max(1, int(round(1 + math.log(b / a) / math.log(q))))
    return n, b / a


def grading(L, secs):
    """blockMesh multi-grading of one block edge and its total cell count."""
    parts = [(f, *section(f * L, a, b)) for f, a, b in secs]
    n = sum(p[1] for p in parts)
    if len(parts) == 1:
        return n, "%.6g" % parts[0][2]
    return n, "(" + " ".join("(%.6g %.6g %.6g)" % (f, c / n, e) for f, c, e in parts) + ")"


def header(cls, obj, loc):
    return ("FoamFile\n{\n    version 2.0;\n    format ascii;\n    class %s;\n    location \"%s\";\n"
            "    object %s;\n}\n\n" % (cls, loc, obj))


def write(case, rel, cls, body):
    path = os.path.join(case, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    loc, obj = os.path.split(rel)
    with open(path, "w") as f:
        f.write(header(cls, obj, loc) + body)


def block_mesh(m):
    xs = [-L_UP, 0.0, L_O, L_O + L_DN]
    rs = [0.0, D_O / 2, R_S]
    c, s = math.cos(math.radians(WEDGE)), math.sin(math.radians(WEDGE))
    verts, idx = [], {}

    def v(i, j, side):
        """vertex index of (x_i, r_j) on the front (-1) or back (+1) wedge plane; one vertex on the axis"""
        key = (i, j, 0 if j == 0 else side)
        if key not in idx:
            idx[key] = len(verts)
            verts.append((xs[i], rs[j] * c, side * rs[j] * s if j else 0.0))
        return idx[key]

    nrc, grc = grading(rs[1], m["rc"])
    nro, gro = grading(rs[2] - rs[1], m["ro"])
    nxu, gxu = grading(L_UP, m["xu"])
    nxo, gxo = grading(L_O, m["xo"])
    nxd, gxd = grading(L_DN, m["xd"])
    nx = [(nxu, gxu), (nxo, gxo), (nxd, gxd)]
    nr = [(nrc, grc), (nro, gro)]
    blocks = []
    cells = 0
    for i in range(3):
        for j in range(2):
            if i == 1 and j == 1:
                continue                       # the piston (plate) around the orifice
            h = [v(i, j, -1), v(i + 1, j, -1), v(i + 1, j + 1, -1), v(i, j + 1, -1),
                 v(i, j, 1), v(i + 1, j, 1), v(i + 1, j + 1, 1), v(i, j + 1, 1)]
            blocks.append("    hex (%s) (%d %d 1) simpleGrading (%s %s 1)"
                          % (" ".join(map(str, h)), nx[i][0], nr[j][0], nx[i][1], nr[j][1]))
            cells += nx[i][0] * nr[j][0]

    def face(i0, j0, i1, j1):
        return "(%d %d %d %d)" % (v(i0, j0, -1), v(i1, j1, -1), v(i1, j1, 1), v(i0, j0, 1))

    def wedge(i, j, side):
        q = [v(i, j, side), v(i + 1, j, side), v(i + 1, j + 1, side), v(i, j + 1, side)]
        return "(%s)" % " ".join(map(str, q if side > 0 else q[::-1]))

    patches = {
        "inlet": ("patch", [face(0, 0, 0, 1), face(0, 1, 0, 2)]),
        "outlet": ("patch", [face(3, 0, 3, 1), face(3, 1, 3, 2)]),
        "upstreamFace": ("wall", [face(1, 1, 1, 2)]),
        "downstreamFace": ("wall", [face(2, 1, 2, 2)]),
        "orificeWall": ("wall", [face(1, 1, 2, 1)]),
        "outerWall": ("wall", [face(0, 2, 1, 2), face(2, 2, 3, 2)]),
        "front": ("wedge", [wedge(i, j, -1) for i in range(3) for j in range(2) if not (i == 1 and j == 1)]),
        "back": ("wedge", [wedge(i, j, 1) for i in range(3) for j in range(2) if not (i == 1 and j == 1)]),
        "axis": ("empty", [face(i, 0, i + 1, 0) for i in range(3)]),
    }
    body = "convertToMeters 0.001;\n\nvertices\n(\n"
    body += "".join("    (%.9g %.9g %.9g)\n" % p for p in verts)
    body += ");\n\nblocks\n(\n" + "\n".join(blocks) + "\n);\n\nedges ();\n\nboundary\n(\n"
    for name, (typ, faces) in patches.items():
        body += "    %s\n    {\n        type %s;\n        faces\n        (\n%s\n        );\n    }\n" % (
            name, typ, "\n".join("            " + f for f in faces))
    body += ");\n"
    return body, cells


WALLS = ["upstreamFace", "downstreamFace", "orificeWall", "outerWall"]


def run_time(u_o):
    return FLOWS * (L_UP + L_O + L_DN) * 1e-3 / u_o


def fields(case, p0, u_o, turb):
    u_in = u_o * (D_O / 2 / R_S) ** 2
    ramp = RAMP if MODE == "steady" else 0.05 * run_time(u_o)
    wedge = "    front\n    {\n        type wedge;\n    }\n    back\n    {\n        type wedge;\n    }\n" \
            "    axis\n    {\n        type empty;\n    }\n"

    def fld(name, cls, dims, internal, inlet, outlet, wall):
        body = "dimensions %s;\ninternalField uniform %s;\nboundaryField\n{\n" % (dims, internal)
        body += "    inlet\n    {\n%s    }\n    outlet\n    {\n%s    }\n" % (inlet, outlet)
        body += "".join("    %s\n    {\n%s    }\n" % (w, wall) for w in WALLS) + wedge + "}\n"
        write(case, "0/" + name, cls, body)

    fld("U", "volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)",
        "        type uniformFixedValue;\n        uniformValue table ((0 (0 0 0)) (%.9g (%.9g 0 0)));\n"
        "        value uniform (0 0 0);\n" % (ramp, u_in),
        "        type inletOutlet;\n        inletValue uniform (0 0 0);\n        value uniform (0 0 0);\n",
        "        type noSlip;\n")
    fld("p", "volScalarField", "[1 -1 -2 0 0 0 0]", "%.9g" % p0,
        "        type zeroGradient;\n",
        "        type fixedValue;\n        value uniform %.9g;\n" % p0,
        "        type zeroGradient;\n")
    fld("T", "volScalarField", "[0 0 0 1 0 0 0]", "%.9g" % T0,
        "        type fixedValue;\n        value uniform %.9g;\n" % T0,
        "        type inletOutlet;\n        inletValue uniform %.9g;\n        value uniform %.9g;\n" % (T0, T0),
        "        type zeroGradient;\n")
    if turb == "laminar":
        return
    k0 = max(1.5 * (0.05 * u_o) ** 2, 1e-10)
    om0 = k0 ** 0.5 / (0.09 ** 0.25 * 0.2e-3)
    for name, val, wall in (("k", k0, "kqRWallFunction"), ("omega", om0, "omegaWallFunction")):
        fld(name, "volScalarField", "[0 2 -2 0 0 0 0]" if name == "k" else "[0 0 -1 0 0 0 0]", "%.9g" % val,
            "        type fixedValue;\n        value uniform %.9g;\n" % val,
            "        type inletOutlet;\n        inletValue uniform %.9g;\n        value uniform %.9g;\n" % (val, val),
            "        type %s;\n        value uniform %.9g;\n" % (wall, val))
    fld("nut", "volScalarField", "[0 2 -1 0 0 0 0]", "0",
        "        type calculated;\n        value uniform 0;\n", "        type calculated;\n        value uniform 0;\n",
        "        type nutUSpaldingWallFunction;\n        value uniform 0;\n")
    fld("alphat", "volScalarField", "[1 -1 -1 0 0 0 0]", "0",
        "        type calculated;\n        value uniform 0;\n", "        type calculated;\n        value uniform 0;\n",
        "        type compressible::alphatWallFunction;\n        Prt 0.85;\n        value uniform 0;\n")


def system(case, turb, iters, u_o):
    probes = ""
    for name, patch, op, flds in (("pUp", "upstreamFace", "areaAverage", "(p rho)"),
                                  ("pDown", "downstreamFace", "areaAverage", "(p rho)"),
                                  ("pIn", "inlet", "areaAverage", "(p rho)"),
                                  ("mIn", "inlet", "sum", "(phi)")):
        probes += ("    %s\n    {\n        type surfaceFieldValue;\n        libs (fieldFunctionObjects);\n"
                   "        writeControl timeStep;\n        writeInterval 20;\n        log false;\n"
                   "        writeFields false;\n        regionType patch;\n        name %s;\n"
                   "        operation %s;\n        fields %s;\n    }\n" % (name, patch, op, flds))
    fo = ("functions\n{\n%s    yPlus\n    {\n        type yPlus;\n"
          "        libs (fieldFunctionObjects);\n        writeControl writeTime;\n    }\n}\n" % probes)
    if MODE == "steady":
        write(case, "system/controlDict", "dictionary",
              "application rhoPimpleFoam;\nstartFrom latestTime;\nstartTime 0;\nstopAt endTime;\n"
              "endTime %d;\ndeltaT 1;\nwriteControl timeStep;\nwriteInterval %d;\npurgeWrite 1;\n"
              "writeFormat ascii;\nwritePrecision 10;\nwriteCompression on;\ntimeFormat general;\n"
              "runTimeModifiable no;\n\n%s" % (iters, iters // 4, fo))
    else:
        te = run_time(u_o)
        write(case, "system/controlDict", "dictionary",
              "application rhoPimpleFoam;\nstartFrom latestTime;\nstartTime 0;\nstopAt endTime;\n"
              "endTime %.9g;\ndeltaT 1e-9;\nwriteControl adjustable;\nwriteInterval %.9g;\npurgeWrite 1;\n"
              "writeFormat ascii;\nwritePrecision 10;\nwriteCompression on;\ntimeFormat general;\n"
              "timePrecision 10;\nrunTimeModifiable no;\nadjustTimeStep yes;\nmaxCo 1;\nmaxDeltaT %.9g;\n\n%s"
              % (te, te / 4, te / 2000, fo.replace("writeInterval 20;", "writeInterval 10;")))
    turbdiv = ("    div(phi,k)      Gauss limitedLinear 1;\n    div(phi,omega)  Gauss limitedLinear 1;\n"
               if turb != "laminar" else "")
    write(case, "system/fvSchemes", "dictionary",
          "ddtSchemes\n{\n    default %s;\n}\ngradSchemes\n{\n    default Gauss linear;\n}\n"
          "divSchemes\n{\n    default none;\n    div(phi,U)      Gauss linearUpwind grad(U);\n"
          "    div(phi,e)      Gauss linearUpwind grad(e);\n    div(phi,K)      Gauss linear;\n"
          "    div(phiv,p)     Gauss linear;\n%s"
          "    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;\n}\n"
          "laplacianSchemes\n{\n    default Gauss linear corrected;\n}\n"
          "interpolationSchemes\n{\n    default linear;\n}\nsnGradSchemes\n{\n    default corrected;\n}\n"
          "wallDist\n{\n    method meshWave;\n}\n" % ("localEuler" if MODE == "steady" else "Euler", turbdiv))
    eqs = "U|e|k|omega" if turb != "laminar" else "U|e"
    if MODE != "steady":                # as the damper cases
        write(case, "system/fvSolution", "dictionary",
              "solvers\n{\n    \"rho.*\"\n    {\n        solver diagonal;\n    }\n"
              "    p\n    {\n        solver GAMG;\n        smoother GaussSeidel;\n        tolerance 1e-8;\n"
              "        relTol 0.01;\n    }\n    pFinal\n    {\n        $p;\n        tolerance 1e-7;\n"
              "        relTol 0;\n        maxIter 200;\n    }\n"
              "    \"(%s)\"\n    {\n        solver smoothSolver;\n        smoother symGaussSeidel;\n"
              "        tolerance 1e-9;\n        relTol 0.1;\n    }\n"
              "    \"(%s)Final\"\n    {\n        $U;\n        relTol 0;\n    }\n}\n\n"
              "PIMPLE\n{\n    momentumPredictor yes;\n    transonic no;\n    nOuterCorrectors 2;\n"
              "    nCorrectors 2;\n    nNonOrthogonalCorrectors 0;\n    pMinFactor 0.2;\n    pMaxFactor 5;\n}\n\n"
              "relaxationFactors\n{\n    equations\n    {\n        \".*\" 1;\n    }\n}\n" % (eqs, eqs))
        return
    write(case, "system/fvSolution", "dictionary",
          "solvers\n{\n    \"rho.*\"\n    {\n        solver diagonal;\n    }\n"
          "    \"p.*\"\n    {\n        solver GAMG;\n        smoother GaussSeidel;\n        tolerance 1e-10;\n"
          "        relTol 0.01;\n    }\n"
          "    \"(%s).*\"\n    {\n        solver smoothSolver;\n        smoother symGaussSeidel;\n"
          "        tolerance 1e-12;\n        relTol 0.1;\n    }\n}\n\n"
          "PIMPLE\n{\n    momentumPredictor yes;\n    transonic no;\n    nOuterCorrectors 1;\n"
          "    nCorrectors 2;\n    nNonOrthogonalCorrectors 0;\n    maxCo 0.2;\n"
          "    rDeltaTSmoothingCoeff 0.02;\n    rDeltaTDampingCoeff 0.5;\n    maxDeltaT 1;\n"
          "    pMinFactor 0.5;\n    pMaxFactor 2;\n}\n\n"
          "relaxationFactors\n{\n    equations\n    {\n        \".*\" 1;\n    }\n}\n" % eqs)


def constant(case, turb):
    write(case, "constant/thermophysicalProperties", "dictionary",
          "thermoType\n{\n    type heRhoThermo;\n    mixture pureMixture;\n    transport sutherland;\n"
          "    thermo hConst;\n    equationOfState perfectGas;\n    specie specie;\n"
          "    energy sensibleInternalEnergy;\n}\n\nmixture\n{\n    specie\n    {\n        molWeight %.6g;\n"
          "    }\n    thermodynamics\n    {\n        Cp %.6g;\n        Hf 0;\n    }\n    transport\n    {\n"
          "        As %.6g;\n        Ts %.6g;\n    }\n}\n" % (M, CP, AS, TS))
    write(case, "constant/turbulenceProperties", "dictionary",
          "simulationType laminar;\n" if turb == "laminar" else
          "simulationType RAS;\n\nRAS\n{\n    RASModel %s;\n    turbulence on;\n    printCoeffs off;\n}\n" % turb)


def cases():
    for mesh in MESHES:
        for turb in TURB[mesh]:
            for p0, u_o in (X_POINTS if mesh.startswith("X") else POINTS):
                re = p0 / (R_GAS * T0) * u_o * D_O * 1e-3 / mu(T0)
                if turb == "laminar" and re > LAMINAR_MAX_RE:
                    continue
                name = "%s_%s_Re%d" % (mesh, "lam" if turb == "laminar" else "sst", round(re))
                yield name, mesh, turb, p0, u_o, re


def setup(root):
    for name, mesh, turb, p0, u_o, re in cases():
        case = os.path.join(root, name)
        body, cells = block_mesh(MESHES[mesh])
        write(case, "system/blockMeshDict", "dictionary", body)
        iters = 10000 if cells < 5000 else 30000 if cells < 40000 else 60000
        system(case, turb, iters, u_o)
        constant(case, turb)
        fields(case, p0, u_o, turb)
        open(os.path.join(case, "case.foam"), "w").close()
        print(name, cells, iters if MODE == "steady" else "%.4g s" % run_time(u_o))


def rows_of(path):
    return [list(map(float, l.split())) for l in open(path) if l.strip() and not l.startswith("#")]


def collect(root):
    """Discharge coefficient averaged over the second half of the run, weighted with the sampling intervals;
    spread: half the range of the instantaneous coefficient over that window, relative to the mean;
    drift: mean over the last quarter of the run minus mean over the third quarter, relative to the mean;
    complete: 1 if the run reached its end time."""
    import glob
    print("case,mesh,turbulence,p0_Pa,Re_o,Cd,spread,Re_L,Cd_Lichtarowicz,time,drift,complete")
    a_o = math.pi * (D_O / 2 * 1e-3) ** 2
    ld = L_O / D_O
    cdu = 0.827 - 0.0085 * ld
    for name, mesh, turb, p0, u_o, re in cases():
        pp = os.path.join(root, name, "postProcessing")
        try:
            up = rows_of(glob.glob(os.path.join(pp, "pUp", "*", "surfaceFieldValue.dat"))[0])
            dn = rows_of(glob.glob(os.path.join(pp, "pDown", "*", "surfaceFieldValue.dat"))[0])
            mi = rows_of(glob.glob(os.path.join(pp, "mIn", "*", "surfaceFieldValue.dat"))[0])
        except (IndexError, OSError):
            print("%s,%s,%s,%g,%.4g,,,,,,," % (name, mesh, turb, p0, re))
            continue
        n = min(len(up), len(dn), len(mi))
        if n < 2:
            print("%s,%s,%s,%g,%.4g,,,,,,," % (name, mesh, turb, p0, re))
            continue
        t_end = up[n - 1][0]
        win = [i for i in range(n - 1) if up[i][0] >= 0.5 * t_end]
        cds, res, wts = [], [], []
        for i in win:
            dp, rho = up[i][1] - dn[i][1], up[i][2]
            mdot = -mi[i][1] * 360.0 / (2 * WEDGE)
            u = mdot / (rho * a_o)
            cds.append(u / math.sqrt(2 * dp / rho) if dp > 0 else float("nan"))
            res.append(rho * u * D_O * 1e-3 / mu(T0))
            wts.append(up[i + 1][0] - up[i][0])
        if not wts or sum(wts) <= 0:
            print("%s,%s,%s,%g,%.4g,,,,,,," % (name, mesh, turb, p0, re))
            continue
        cd = sum(c * w for c, w in zip(cds, wts)) / sum(wts)
        re = sum(r * w for r, w in zip(res, wts)) / sum(wts)
        spread = (max(cds) - min(cds)) / 2 / cd if cd == cd else float("nan")
        k0 = next((k for k, i in enumerate(win) if up[i][0] >= 0.75 * t_end), len(win))
        third, last = (sum(cds[k] * wts[k] for k in ks) / sum(wts[k] for k in ks) if sum(wts[k] for k in ks) > 0
                       else float("nan") for ks in (range(k0), range(k0, len(win))))
        drift = (last - third) / cd
        complete = int(t_end >= 0.99 * run_time(u_o))
        re_l = re / cd if cd > 0 else float("nan")
        lich = float("nan")
        if re_l == re_l:
            lich = 1 / (1 / cdu + 20 / re_l * (1 + 2.25 * ld)
                        - 0.005 * ld / (1 + 7.5 * math.log10(0.00015 * re_l) ** 2))
        print("%s,%s,%s,%g,%.4g,%.4f,%.1e,%.4g,%.4f,%.4g,%.1e,%d" % (name, mesh, turb, p0, re, cd, spread, re_l, lich,
                                                                     t_end, drift, complete))


if __name__ == "__main__":
    if len(sys.argv) > 3:
        MODE = sys.argv[3]
    {"setup": setup, "collect": collect}[sys.argv[1]](sys.argv[2])
