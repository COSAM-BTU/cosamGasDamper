# cosamGasDamper

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24042536.svg)](https://doi.org/10.5281/zenodo.24042536)

Automated, end-to-end CFD of large-stroke, gas-filled piston–cylinder dampers
and gas springs with one to ten orifices, with OpenFOAM and SALOME.

You describe a damper in one YAML file (bore, rod and orifice dimensions,
number of orifices, stroke, piston speed, pre-charge pressure, mesh density,
number of cores) and run one command. The framework builds the parametric
geometry and hexahedral meshes, splits the stroke into mesh segments, runs
`rhoPimpleFoam` with a moving mesh, transfers the solution between segments,
checks the gas mass and volume at every time step, and writes time-series CSV
files, a case summary with the spring and damping forces and the dissipated
energy, and rendered mid-plane images. A sweep file expands into many cases
that are executed under a core budget, so that a parameter study runs
unattended.

> This offering is not approved or endorsed by OpenCFD Limited, producer and
> distributor of the OpenFOAM software via www.openfoam.com, and owner of the
> OPENFOAM® and OpenCFD® trade marks.

## What it does

```
config.yaml ──► geometry & section mesh (SALOME, N-orifice parametric, Gmsh quad)
            ──► K sub-meshes swept along the axis (one per piston position)
            ──► mesh segments: rhoPimpleFoam + dynamic mesh (velocityComponentLaplacian)
                  │  piston velocity from one Function1 table (no restarts on reversal)
                  └─ mapFields between segments (each sub-mesh deforms by ≤ S/K);
                     optional volume-weighted, mass-conserving transfer (mapFieldsPar)
            ──► results/timeseries.csv, summary.json, segments.csv, quicklook.png
            ──► frames/*.png (pvbatch, headless)
```

* Any number of equally spaced orifices `N ≥ 1`; faces and zones are classified
  from their positions, so the generator has no per-N templates.
* All sub-meshes of a case share one cross-section, which is stored
  (`mesh/section.med`) and can be reused.
* Built-in checks: analytic gas volume of every mesh, `checkMesh`, solver end
  time, and, at every time step, total gas mass and the CFD gas volume against
  the analytic volume of the prescribed piston position.
* Restartable: a finished stage or segment is not repeated
  (`run --from STAGE`, `run --only STAGE`).

## Requirements

Tested with

| tool | version |
|---|---|
| OpenFOAM (ESI/OpenCFD) | v2406 |
| SALOME | 9.12 (with Gmsh 4.10) |
| ParaView (`pvbatch`) | 5.11 – 6.1 |
| Python | 3.10 with PyYAML, NumPy, Matplotlib |
| MPI | Open MPI 4 |

on Ubuntu 22.04. The locations of OpenFOAM, SALOME and `pvbatch` are keys of
the `run` section of the configuration (`config/defaults.yaml`).

## Installation

No build step is needed: clone the repository and make sure that the tools
above are installed.

```bash
git clone https://github.com/COSAM-BTU/cosamGasDamper.git
cd cosamGasDamper
PYTHONPATH=src python3 tests/test_schedule.py     # nine unit tests, no OpenFOAM needed
```

The tests cover the configuration validation, the mesh-segment schedule, the
velocity table and the defaults of the options; with pytest installed,
`PYTHONPATH=src python3 -m pytest tests` runs the same tests.

## Quick start

```bash
./run.sh check tutorials/01_quick/config.yaml          # validate, print the segment schedule
./run.sh run runs/quick --config tutorials/01_quick/config.yaml
```

The tutorial (ten orifices, 1 m/s, one compression–rebound cycle on the
medium mesh) runs from the configuration file to the images in about half an
hour on eight cores; `tutorials/01_quick/README.md` lists the expected output.

A parameter study:

```bash
./run.sh sweep examples/paper/01_matrix.yaml --root runs/matrix --max-cores 48
./run.sh collect runs/matrix --out matrix.csv           # one table of all case summaries
```

## Configuration

`config/defaults.yaml` lists every key with its default; a case file only
overrides what differs. Units are part of the key names.

```yaml
case:
  name: N4_v0.4
geometry:
  cylinder_radius_mm: 10.0
  rod_radius_mm: 4.0
  cylinder_length_mm: 120.0
  piston_length_mm: 8.0
  orifice: {count: 4, radius_mm: 1.0, pitch_radius_mm: 7.0}
motion:
  stroke_mm: 90.0
  mean_speed_mps: 0.4          # or stroke_time_s
  profile: constant            # constant | sinusoidal
  meshes_per_stroke: 4
gas:
  initial_pressure_Pa: 1.0e+5
  turbulence: kOmegaSST        # or laminar
mesh:
  piston_cell_size_mm: 1.0     # axial cell size inside the orifices
  face_cell_size_mm: null      # axial cell size next to the piston faces, growing to chamber_cell_size_mm;
                               # 0.1 resolves the separation at the orifice inlet; null = off
  section:
    curvature_elements: null   # local refinement of the orifices (8 for the reference mesh of the article); null = off
    preserve_orifice_area: false   # true: scale each orifice polygon to the area pi*r_o^2
numerics:
  map_method: interpolate      # or cellVolumeWeight: transfer between sub-meshes that conserves the gas mass
run:
  n_procs: 8
```

The discretised orifice circles are inscribed polygons whose area is
`(n/2π) sin(2π/n)` times `π r_o²` for `n` boundary segments (0.989 for the
reference mesh, 0.936 for the default section); `mesh/mesh_report.json`
records this ratio (`orifice_area_ratio`), and discharge coefficients should be
based on the discrete area unless `preserve_orifice_area` is set.

A sweep file adds either a Cartesian product

```yaml
sweep:
  geometry.orifice.count: [2, 4, 6, 8, 10]
  motion.mean_speed_mps: [0.025, 0.1, 0.4, 1.0, 2.0, 4.0]
```

or a list of named `variants` with dotted-key overrides (see
`examples/paper/02_mesh.yaml`). A relative `mesh.section_med` (a stored
cross-section to reuse) is relative to the directory the sweep is started from.

## Output of a case

```
runs/<case>/
  config.yaml  config.resolved.json  state.json  run.log
  mesh/        section.med, M0..M{K-1}/constant/polyMesh, mesh_report.json
  segments/    seg01_M0 ... (OpenFOAM cases, reconstructed)
  results/     timeseries.csv  summary.json  segments.csv  quicklook.png
  frames/      p_*.png  U_*.png  T_*.png
  logs/
```

`timeseries.csv` holds, per time step: piston position and velocity, the
axial force on the piston assembly (total, pressure, viscous, gauge), the
volume-averaged pressure, temperature and density of both chambers and the
orifices, gas masses, CFD and analytic gas volume, orifice mass flow, mean and
peak orifice Mach number, time step and Courant number. `summary.json`
holds the spring/damping decomposition, the energy dissipated per cycle, the
polytropic exponent of the compression, the largest mass and volume errors and
the solver cost.

## Reproducing the article

`examples/paper/` holds the configuration files of all studies of the
accompanying article, the single-orifice benchmark and, in
`examples/paper/analysis/`, the scripts that produce its figures and tables
from the case records; `examples/paper/README.md` explains the order of the
runs. The case summaries, mesh reports, stored cross-sections and derived
tables of the published runs are archived as a separate dataset on Zenodo.

## Repository layout

```
run.sh                      entry point
config/defaults.yaml        all keys and defaults
src/damper/                 driver (config, schedule, OpenFOAM case writer, stages, post, render)
src/salome/damper_mesh.py   parametric geometry and meshing (runs inside SALOME)
src/render/render_slices.py pvbatch script
tests/                      unit tests
tutorials/                  small case with its expected output
examples/paper/             studies, benchmark and analysis scripts of the accompanying article
```

## Citing

If you use this software, please cite the accompanying article and the
archived release (see `CITATION.cff`):

L. Aydınbakar, M. Gömlekli, An open-source multi-mesh OpenFOAM
framework for the automated simulation of large-stroke gas-filled
piston–cylinder dampers with multiple orifices, submitted to Advances in
Engineering Software (2026).

Software, all versions: https://doi.org/10.5281/zenodo.24042536.
Data of the article (case records, derived tables, figure data):
https://doi.org/10.5281/zenodo.24042564.

## Licence

GPL-3.0-or-later, see `LICENSE`.

## Authors

Levent Aydınbakar and Murat Gömlekli, Department of Mechanical Engineering,
Bursa Technical University, Türkiye.
