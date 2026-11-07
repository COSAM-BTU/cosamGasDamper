# SPDX-License-Identifier: GPL-3.0-or-later
"""Writers for the OpenFOAM dictionaries of one mesh segment."""

import os

HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    location    "{loc}";
    object      {obj};
}}
"""

MOVING = ("pistonFront", "pistonBack")          # piston faces
CARRIED = ("orificeWall", "rodWall")            # walls that translate with the piston
FIXED = ("bottomCap", "topCap", "cylinderWall")  # stationary walls
WALLS = MOVING + CARRIED + FIXED


def write(path, cls, obj, body, loc=None):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if loc is None:
        loc = os.path.basename(os.path.dirname(path))
    with open(path, "w") as f:
        f.write(HEADER.format(cls=cls, loc=loc, obj=obj))
        f.write(body)


def table_scalar(tab, sign=1.0):
    return "table\n        (\n" + "".join(
        "            (%.12g %.12g)\n" % (t, sign * u) for t, u in tab) + "        )"


def table_vector(tab):
    return "table\n        (\n" + "".join(
        "            (%.12g (%.12g 0 0))\n" % (t, u) for t, u in tab) + "        )"


def fields(case, tdir, cfg, tab):
    """Initial/boundary fields in the start-time directory of a segment."""
    gas = cfg["gas"]
    d = os.path.join(case, tdir)
    p0, T0 = gas["initial_pressure_Pa"], gas["initial_temperature_K"]
    zg = "".join("    %s\n    {\n        type zeroGradient;\n    }\n" % w for w in WALLS)

    write(os.path.join(d, "p"), "volScalarField", "p",
          "dimensions [1 -1 -2 0 0 0 0];\ninternalField uniform %.10g;\nboundaryField\n{\n%s}\n"
          % (p0, zg), tdir)
    write(os.path.join(d, "T"), "volScalarField", "T",
          "dimensions [0 0 0 1 0 0 0];\ninternalField uniform %.10g;\nboundaryField\n{\n%s}\n"
          % (T0, zg), tdir)

    u = ""
    for w in MOVING:
        u += "    %s\n    {\n        type movingWallVelocity;\n        value uniform (0 0 0);\n    }\n" % w
    for w in CARRIED:
        u += ("    %s\n    {\n        type uniformFixedValue;\n        uniformValue %s;\n"
              "        value uniform (0 0 0);\n    }\n") % (w, table_vector(tab))
    for w in FIXED:
        u += "    %s\n    {\n        type noSlip;\n    }\n" % w
    write(os.path.join(d, "U"), "volVectorField", "U",
          "dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (0 0 0);\nboundaryField\n{\n%s}\n"
          % u, tdir)

    pm = ""
    for w in MOVING + ("orificeWall",):
        pm += ("    %s\n    {\n        type uniformFixedValue;\n        uniformValue %s;\n"
               "        value uniform 0;\n    }\n") % (w, table_scalar(tab))
    for w in ("bottomCap", "topCap"):
        pm += "    %s\n    {\n        type fixedValue;\n        value uniform 0;\n    }\n" % w
    for w in ("cylinderWall", "rodWall"):
        pm += "    %s\n    {\n        type slip;\n    }\n" % w
    write(os.path.join(d, "pointMotionUx"), "pointScalarField", "pointMotionUx",
          "dimensions [0 1 -1 0 0 0 0];\ninternalField uniform 0;\nboundaryField\n{\n%s}\n"
          % pm, tdir)

    if gas["turbulence"] != "laminar":
        I = gas["turbulence_intensity"]
        lt = gas["turbulence_length_scale_mm"] * 1e-3
        uref = max(abs(u_) for _, u_ in tab) * cfg["_derived"]["orifice_area_ratio"]
        k0 = max(1.5 * (I * uref) ** 2, 1e-8)
        om0 = k0 ** 0.5 / (0.09 ** 0.25 * lt)
        nu_wall = "".join("    %s\n    {\n        type nutUSpaldingWallFunction;\n"
                          "        value uniform 0;\n    }\n" % w for w in WALLS)
        k_wall = "".join("    %s\n    {\n        type kqRWallFunction;\n"
                         "        value uniform %.6g;\n    }\n" % (w, k0) for w in WALLS)
        om_wall = "".join("    %s\n    {\n        type omegaWallFunction;\n"
                          "        value uniform %.6g;\n    }\n" % (w, om0) for w in WALLS)
        at_wall = "".join("    %s\n    {\n        type compressible::alphatWallFunction;\n"
                          "        Prt 0.85;\n        value uniform 0;\n    }\n" % w for w in WALLS)
        write(os.path.join(d, "k"), "volScalarField", "k",
              "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform %.6g;\nboundaryField\n{\n%s}\n"
              % (k0, k_wall), tdir)
        write(os.path.join(d, "omega"), "volScalarField", "omega",
              "dimensions [0 0 -1 0 0 0 0];\ninternalField uniform %.6g;\nboundaryField\n{\n%s}\n"
              % (om0, om_wall), tdir)
        write(os.path.join(d, "nut"), "volScalarField", "nut",
              "dimensions [0 2 -1 0 0 0 0];\ninternalField uniform 0;\nboundaryField\n{\n%s}\n"
              % nu_wall, tdir)
        write(os.path.join(d, "alphat"), "volScalarField", "alphat",
              "dimensions [1 -1 -1 0 0 0 0];\ninternalField uniform 0;\nboundaryField\n{\n%s}\n"
              % at_wall, tdir)


def constant(case, cfg):
    gas, me = cfg["gas"], cfg["mesh"]
    d = os.path.join(case, "constant")
    write(os.path.join(d, "thermophysicalProperties"), "dictionary", "thermophysicalProperties", """
thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       sutherland;
    thermo          hConst;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleInternalEnergy;
}

mixture
{
    specie
    {
        molWeight   %.8g;
    }
    thermodynamics
    {
        Cp          %.8g;
        Hf          0;
    }
    transport
    {
        As          %.8g;
        Ts          %.8g;
    }
}
""" % (gas["molar_mass_kg_per_kmol"], gas["cp_J_per_kgK"], gas["sutherland_As"],
       gas["sutherland_Ts_K"]))
    if gas["turbulence"] == "laminar":
        turb = "simulationType laminar;\n"
    else:
        turb = ("simulationType RAS;\n\nRAS\n{\n    RASModel        %s;\n    turbulence      on;\n"
                "    printCoeffs     off;\n}\n" % gas["turbulence"])
    write(os.path.join(d, "turbulenceProperties"), "dictionary", "turbulenceProperties", turb)
    diff = me["motion_diffusivity"]
    write(os.path.join(d, "dynamicMeshDict"), "dictionary", "dynamicMeshDict", """
dynamicFvMesh       dynamicMotionSolverFvMesh;
motionSolverLibs    (fvMotionSolvers);
motionSolver        velocityComponentLaplacian;
component           x;
diffusivity         directional (%g %g %g);
""" % tuple(diff))


FUNCTIONS = """
functions
{
    forces
    {
        type            forces;
        libs            (forces);
        writeControl    timeStep;
        writeInterval   1;
        log             false;
        patches         (pistonFront pistonBack orificeWall rodWall);
        rho             rho;
        p               p;
        U               U;
        pRef            0;
        CofR            (0 0 0);
    }
    forcesFront
    {
        $forces;
        patches         (pistonFront);
    }
    forcesBack
    {
        $forces;
        patches         (pistonBack);
    }
    chamber1
    {
        type            volFieldValue;
        libs            (fieldFunctionObjects);
        writeControl    timeStep;
        writeInterval   1;
        log             false;
        writeFields     false;
        regionType      cellZone;
        name            chamber1;
        operation       volAverage;
        fields          (p T rho);
    }
    chamber2
    {
        $chamber1;
        name            chamber2;
    }
    orifices
    {
        $chamber1;
        name            orifices;
    }
    all
    {
        $chamber1;
        regionType      all;
    }
    mass1
    {
        $chamber1;
        operation       volIntegrate;
        fields          (rho);
    }
    mass2
    {
        $mass1;
        name            chamber2;
    }
    massOrifices
    {
        $mass1;
        name            orifices;
    }
    massAll
    {
        $mass1;
        regionType      all;
    }
    orificeFlux
    {
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        writeControl    timeStep;
        writeInterval   1;
        log             false;
        writeFields     false;
        surfaceFormat   none;
        regionType      faceZone;
        name            orificeOutlet;
        operation       sum;
        fields          (phi);
    }
    orificeState
    {
        $orificeFlux;
        operation       areaAverage;
        fields          (p T rho);
    }
    magU
    {
        type            mag;
        libs            (fieldFunctionObjects);
        field           U;
        executeControl  timeStep;
        writeControl    writeTime;
    }
    orificeUmax
    {
        $chamber1;
        name            orifices;
        operation       max;
        fields          (mag(U));
    }
}
"""


def system(case, cfg, seg, n_procs):
    nu, d = cfg["numerics"], cfg["_derived"]
    # write interval that divides the segment exactly: the adjustable time step
    # lands on write times (counted from startTime), so the segment end is hit exactly
    span = seg["end"] - seg["start"]
    n_w = max(1, int(round(span / d["write_interval_s"])))
    write_dt = span / n_w
    p0 = cfg["gas"]["initial_pressure_Pa"]
    sdir = os.path.join(case, "system")
    write(os.path.join(sdir, "controlDict"), "dictionary", "controlDict", """
application     rhoPimpleFoam;
startFrom       startTime;
startTime       %s;
stopAt          endTime;
endTime         %s;
deltaT          %.6g;
writeControl    adjustable;
writeInterval   %.10g;
purgeWrite      0;
writeFormat     ascii;
writePrecision  10;
writeCompression on;
timeFormat      general;
timePrecision   %d;
runTimeModifiable no;
adjustTimeStep  yes;
maxCo           %.6g;
maxDeltaT       %.6g;
%s""" % (seg["start_name"], seg["end_name"], nu["initial_dt_s"], write_dt,
         nu["time_precision"], nu["max_courant"], d["max_dt_s"], FUNCTIONS))

    turb = cfg["gas"]["turbulence"] != "laminar"
    write(os.path.join(sdir, "fvSchemes"), "dictionary", "fvSchemes", """
ddtSchemes
{
    default         Euler;
}
gradSchemes
{
    default         Gauss linear;
}
divSchemes
{
    default         none;
    div(phi,U)      Gauss linearUpwind grad(U);
    div(phi,e)      Gauss linearUpwind grad(e);
    div(phi,K)      Gauss linear;
    div(phiv,p)     Gauss linear;
    div(meshPhi,p)  Gauss linear;
    div(phi,k)      Gauss limitedLinear 1;
    div(phi,omega)  Gauss limitedLinear 1;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes
{
    default         Gauss linear corrected;
}
interpolationSchemes
{
    default         linear;
}
snGradSchemes
{
    default         corrected;
}
wallDist
{
    method          meshWave;
}
""")
    eqs = "U|e|k|omega" if turb else "U|e"
    write(os.path.join(sdir, "fvSolution"), "dictionary", "fvSolution", """
solvers
{
    "rho.*"
    {
        solver          diagonal;
    }
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          0.01;
    }
    pFinal
    {
        $p;
        tolerance       1e-7;
        relTol          0;
        maxIter         200;
    }
    "pcorr.*"
    {
        $p;
        tolerance       1e-6;
        relTol          0;
    }
    "(%s)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-9;
        relTol          0.1;
    }
    "(%s)Final"
    {
        $U;
        relTol          0;
    }
    cellMotionUx
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0;
    }
}

PIMPLE
{
    momentumPredictor        yes;
    transonic                no;
    nOuterCorrectors         %d;
    nCorrectors              %d;
    nNonOrthogonalCorrectors %d;
    correctPhi               no;
    pMin                     %.10g;
    pMax                     %.10g;
}

relaxationFactors
{
    equations
    {
        ".*"        1;
    }
}
""" % (eqs, eqs, nu["n_outer_correctors"], nu["n_correctors"], nu["n_non_orth_correctors"],
       nu["p_min_factor"] * p0, nu["p_max_factor"] * p0))
    write(os.path.join(sdir, "decomposeParDict"), "dictionary", "decomposeParDict", """
numberOfSubdomains %d;
method          scotch;
""" % n_procs)


def mesh_system(case):
    """Minimal system/ needed by the mesh utilities."""
    sdir = os.path.join(case, "system")
    write(os.path.join(sdir, "controlDict"), "dictionary", "controlDict", """
application     rhoPimpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         1;
deltaT          1;
writeControl    timeStep;
writeInterval   1;
writeFormat     ascii;
writePrecision  12;
writeCompression off;
timeFormat      general;
timePrecision   10;
""")
    write(os.path.join(sdir, "fvSchemes"), "dictionary", "fvSchemes",
          "ddtSchemes{default Euler;}\ngradSchemes{default Gauss linear;}\n"
          "divSchemes{default none;}\nlaplacianSchemes{default Gauss linear corrected;}\n"
          "interpolationSchemes{default linear;}\nsnGradSchemes{default corrected;}\n")
    write(os.path.join(sdir, "fvSolution"), "dictionary", "fvSolution", "\n")
