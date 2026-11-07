# Analysis scripts of the paper

These scripts turn the case records written by the framework into the figure data, tables and numbers of the
paper *An open-source multi-mesh OpenFOAM framework for the automated simulation of large-stroke gas-filled
piston–cylinder dampers with multiple orifices* and its supplementary material.

| Script | Purpose |
|---|---|
| `figdata.py` | data files of the line plots (one sub-folder per figure) and the derived case tables; also the shared helper functions |
| `paper_numbers.py` | numbers quoted in the text of Section 3 (printed) |
| `verification_numbers.py` | in-plane discretisation uncertainty, mesh-study tables, checks of the reduced model and of the force decomposition, further cases (printed; the section `inplane` also writes three CSV tables) |
| `bench_table.py` | LaTeX rows of the single-orifice benchmark table (printed) |
| `make_field_specs.py` | render specs of the velocity and pressure field sequences, for `src/render/render_slices.py` |
| `field_ranges.py` | data ranges of the fields of a render spec, to choose fixed colour ranges (pvbatch) |
| `render_normalized.py` | images of the velocity magnitude over the ideal orifice velocity (pvbatch) |
| `render_sections.py` | images of the cross-section mesh at the bottom cap (pvbatch) |
| `render_mesh_compare.py` | images of the cross-section mesh, whole section and close-up of the first orifice (pvbatch) |
| `render_side.py` | images of the axial section of a sub-mesh, full length and close-up (pvbatch) |
| `render_3d.py` | perspective images of the boundary mesh of a sub-mesh (pvbatch) |
| `compose_view3d.py` | arranges the images of `render_3d.py` into one figure (`view3d.tex`) |

## Requirements

* Python 3 with NumPy (`figdata.py`, `paper_numbers.py`, `verification_numbers.py`, `bench_table.py`,
  `make_field_specs.py`; checked with Python 3.14 and NumPy 2.4); Pillow in addition for `compose_view3d.py`.
* ParaView with `pvbatch` (5.13 for the paper) for `field_ranges.py` and the `render_*.py` scripts, together with
  `src/render/render_slices.py` of this repository for the field images. These scripts need the meshes
  (`mesh/M*/constant/polyMesh`) or the field files (`segments/seg*/`) of the cases, which the published dataset
  does not contain; they are produced by running the cases (for the mesh images the mesh stage,
  `run --only mesh`, is enough).

## Inputs and outputs

All paths are command-line options with defaults relative to the current working directory:

* `--runs` (default `runs`): directory holding the study directories, `<runs>/<study>/<case>/`. This is either the
  `runs` directory of your own sweeps or the `cases/` folder of the published dataset. The studies and the
  configuration files in `examples/paper` that produce them:

  | Study | Configuration | Content |
  |---|---|---|
  | `matrix` | `01_matrix.yaml` | parametric matrix, medium mesh |
  | `studies` | `02`–`10` | mesh, time-step, model and segment studies on the medium mesh; orifice refinement with axial cells of 1 mm |
  | `matrixC8` | `11_matrix_C8.yaml` | matrix with refined orifices (curvature size 8), axial cells of 1 mm |
  | `studiesC8` | `12_studies_C8.yaml` | studies on that mesh |
  | `furtherC8` | `13_further_C8.yaml` | further cases on that mesh |
  | `matrixC8ax` | `14_matrix_ref.yaml` | matrix on the reference mesh (curvature size 8, axial cells of 0.1 mm) |
  | `studiesC8ax` | `15_studies_ref.yaml` | verification and sensitivity studies on the reference mesh |
  | `sections` | `16_sections_ref.yaml` | meshes of N = 1 to 10 for the mesh figures (mesh stage only) |
  | `furtherC8ax` | `17_further_ref.yaml` | further cases on the reference mesh |

  The sweeps are started from the `runs` directory, for example
  `cd runs && ../run.sh sweep ../examples/paper/14_matrix_ref.yaml --root matrixC8ax`; see `examples/paper/README.md`
  for all commands and their order.
* `--tables` (default `tables`): derived CSV tables. Read: `orifice_benchmark.csv` and `orifice_benchmark_T2.csv`
  (the output of `examples/paper/orifice_benchmark.py collect` for the first set of benchmark runs and for the
  mesh-B runs repeated with `Co_max = 0.5`), `inplane_est.csv` and `inplane_cdo.csv`. Written:
  `inplane_bias.csv`, `inplane_cdo.csv`, `inplane_est.csv` (`verification_numbers.py inplane`) and
  `energy_balance.csv`, `orifice_convergence.csv`, `cases_table.csv` (`figdata.py`).
* `--out` (default `figures`, `figdata.py` only): one sub-folder per figure, named as the figure folders of the
  paper. Every data file starts with a line `# sources: ...` naming its source case directories relative to
  `--runs`.

The per-case files read by the scripts are `results/summary.json`, `config.resolved.json`,
`mesh/mesh_report.json`, `state.json` and `results/timeseries.csv`; in addition `figdata.py cases` checks for
the log of the batch scheduler (`slurm-*.out`) to label the machine of each case.

## Running

From the directory that holds `runs/` (paths below relative to the repository root):

```sh
A=examples/paper/analysis
python3 $A/verification_numbers.py inplane          # tables/inplane_{bias,cdo,est}.csv, needed by figdata.py
python3 $A/figdata.py                               # figures/<figure>/*.dat, tables/*.csv, GCI of the mesh studies
python3 $A/paper_numbers.py > paper_numbers.txt     # numbers of Section 3
python3 $A/verification_numbers.py > verification.txt
python3 $A/verification_numbers.py levels uncertainty
python3 $A/bench_table.py
```

* `figdata.py [--runs] [--tables] [--out] [output ...]` with outputs `loops ratio cd benchmark mass spring history
  energy orifice cases` (default: all). Without `inplane_est.csv` the in-plane estimates are left out of
  `dampingDeviation`; without `inplane_cdo.csv` `damperExt.dat` falls back to the two 1 m/s conditions.
* `paper_numbers.py [--runs] [--sections a,b,...] [matrix [aux]]` with sections `fields energy conservation
  segments gas_spring cd_ratio cd_details cost studies axial` (default: all) and the studies of the matrix and of
  the sensitivity studies (default `matrixC8ax studiesC8ax`; `matrixC8 studiesC8` gives the values with axial
  cells of 1 mm, `matrix studies` those of the medium mesh).
* `verification_numbers.py [--runs] [--tables] [section ...]` with sections `area kinematic branches window hyst
  gci meshstudy inplane levels further guidance uncertainty`; the default runs all except `levels` and
  `uncertainty`, which run `inplane` themselves. The script takes about one to two minutes.
* `bench_table.py [--tables]`.
* Field figures: `make_field_specs.py [--runs] [--study matrixC8ax] [--out field_specs] [--ranges ranges.json]
  [--rebound]` writes `spec_N2.json`, `spec_N6.json`, `spec_N10.json` (or `spec_R2.json`);
  `pvbatch field_ranges.py spec_N2.json` prints the data ranges on the slice, from which the fixed ranges of
  `ranges.json` were chosen (paper: `{"2": {"U": [0, 180], "p": [90000, 140000]}, "6": {"U": [0, 60],
  "p": [100000, 125000]}, "10": {"U": [0, 30], "p": [100000, 125000]}}`, the rebound with the ranges of N = 2);
  `pvbatch src/render/render_slices.py spec_N2.json` writes the images.
* `pvbatch render_normalized.py out_dir vmax case_dir ...` (paper: `vmax` = 1.8, cases
  `matrixC8ax/N{2,10}_v{0.4,1.0,2.0,4.0}`; `vmax` = 0 prints the data ranges only).
* `pvbatch render_sections.py out_dir case_dir ...` (paper: `sections/sec_N1` ... `sec_N10`).
* `pvbatch render_side.py out_dir xmin xmax zmin zmax label=mesh_dir ...` (paper: window
  `0.046 0.074 -0.0005 0.010`, `N1=runs/sections/sec_N1/mesh/M2`, `N10=runs/sections/sec_N10/mesh/M2`).
* `pvbatch render_mesh_compare.py out_dir label=case_dir ...` (paper: labels `N2med N2C8 N2C12 N10med N10C8
  N10C12` for the medium mesh, the reference cross-section C = 8 and the level C = 12 of N = 2 and 10).
* `pvbatch render_3d.py out_dir N1=runs/sections/sec_N1/mesh/M2 ... N10=runs/sections/sec_N10/mesh/M2`, then
  `compose_view3d.py out_dir`.

Equations are referred to in the code by their labels in the paper: `gasSpring` (Eq. 8), `ratio` (Eq. 9),
`kappa` (Eq. 10), `hysteresis` (Eq. 11), `lichtarowicz` (Eq. 12).

## Figures and tables of the paper

Data folder: sub-folder of `--out`. Time series: whether the output needs `results/timeseries.csv` of the cases,
which the published dataset does not contain (see below).

| Item | Label | Script, function / section | Data | Time series |
|---|---|---|---|---|
| Fig. 1 | `fig:geometry` | schematic (TikZ) | – | – |
| Fig. 2 | `fig:architecture` | schematic (TikZ) | – | – |
| Fig. 3 | `fig:sections` | `render_sections.py` | `sections/sec_N*` meshes | meshes |
| Fig. 4 | `fig:sideMesh` | `render_side.py` | `sections/sec_N1`, `sec_N10` meshes | meshes |
| Fig. 5 | `fig:meshCompare` | `render_mesh_compare.py` | meshes of the medium, C = 8 and C = 12 cross-sections | meshes |
| Fig. 6 | `fig:massVolume` | `figdata.py mass` (`fig_mass`) | `massVolume/` | yes |
| Fig. 7 | `fig:forceLoops` | `figdata.py loops` (`fig_loops`) | `forceLoops/` | yes |
| Fig. 8 | `fig:pressureHistory` | `figdata.py history` (`fig_history`) | `pressureHistory/` | yes |
| Fig. 9 | `fig:dischargeCoefficient` | `figdata.py cd` (`fig_cd`) | `dischargeCoefficient/` | yes |
| Fig. 10 | `fig:orificeBenchmark` | `figdata.py benchmark` (`fig_benchmark`); `verification_numbers.py inplane` | `orificeBenchmark/` | yes |
| Fig. 11 | `fig:fieldsSpeed` | `render_normalized.py` | `matrixC8ax` fields | fields |
| Fig. 12 | `fig:fieldsN2` | `make_field_specs.py`, `field_ranges.py`, `src/render/render_slices.py` | `matrixC8ax/N2_v2.0` fields | fields |
| Fig. 13 (a) | `fig:dampingRatio` | `figdata.py ratio` (`fig_ratio`) | `dampingRatio/` | yes |
| Fig. 13 (b) | `fig:dampingRatio` | `figdata.py ratio` (`fig_ratio`); `verification_numbers.py inplane` | `dampingDeviation/` | yes |
| Tables 1–4 | | text | – | – |
| Table 5 | `tab:studies` | configuration files `examples/paper/*.yaml` | – | – |
| Table 6 | `tab:uncertainty` | `verification_numbers.py uncertainty` (`sec_uncertainty`) | printed | yes |
| Table 7 | `tab:cost` | `paper_numbers.py --sections cost`; core-hours per case from `tables/cases_table.csv` (`solver_h` x `cores`) | printed | no |
| Fig. S1 | `fig:view3d` | `render_3d.py`, `compose_view3d.py` (the N labels were added to `view3d.tex` afterwards) | `sections/sec_N*` meshes | meshes |
| Figs. S2, S3 | `fig:fieldsN6`, `fig:fieldsN10` | as Fig. 12 | `matrixC8ax/N6_v2.0`, `N10_v2.0` fields | fields |
| Fig. S4 | `fig:fieldsR2` | as Fig. 12 with `make_field_specs.py --rebound` | `matrixC8ax/N2_v2.0` fields | fields |
| Fig. S5 | `fig:machHistory` | `figdata.py history` (`fig_history`) | `machHistory/` | yes |
| Table S1 | `tab:meshStudy` | `verification_numbers.py meshstudy`; estimate row from `inplane` | printed | yes |
| Table S2 | `tab:levels` | `verification_numbers.py levels` (orders); GCI per choice from `inplane` | printed | yes |
| Table S3 | `tab:segments` | `paper_numbers.py --sections segments` | printed | no |
| Table S4 | `tab:benchmark` | `examples/paper/orifice_benchmark.py` (meshes of the benchmark) | – | – |
| Table S5 | `tab:benchmarkCd` | `bench_table.py` (cells without a settled value, `--` in the output, are labelled s, u, n in the paper from the state of the run) | `tables/orifice_benchmark*.csv` | no |

`figdata.py spring` writes `gasSpring/` (force of the low-speed cases against the adiabatic gas spring), which is
not shown in the paper.

Numbers in the text: Section 3 from `paper_numbers.py`; the mesh-study grid convergence indices of Section 3.1
from `figdata.py orifice` (printed); the in-plane uncertainty (3.1), the area of the orifice polygons (3.1, S1),
the thermal hysteresis (3.3), the kinematic factor and the model comparison (3.5), the equal-area and 20 bar cases
(3.5) and the force decomposition (S3) from `verification_numbers.py`; the energy balance (3.3) from
`paper_numbers.py --sections energy` and `tables/energy_balance.csv`.

## What can be regenerated from the published dataset

The dataset holds every case record except `results/timeseries.csv` (and the meshes and fields). The time series
is needed wherever an output uses a history (force loops, pressures, mass) or the orifice discharge coefficient
`C_d,o`, which is averaged from the columns `p1`, `p2`, `rho1` and `mdot` over 0.4–0.6 of the compression stroke
(`figdata.cd_orifice`).

| Output | Per-case files read | From the dataset alone |
|---|---|---|
| `figdata.py loops`, `mass`, `spring`, `history` | `timeseries.csv`, `config.resolved.json` | no (`mass` skips cases without it) |
| `figdata.py ratio`, `cd` | `summary.json`, `config.resolved.json`, `mesh_report.json`, `timeseries.csv`; tables `inplane_est.csv` | no |
| `figdata.py benchmark` | `config.resolved.json`, `mesh_report.json`, `timeseries.csv`; tables `orifice_benchmark*.csv`, `inplane_cdo.csv` | no |
| `figdata.py energy` | `timeseries.csv`, `summary.json`, `config.resolved.json` | no |
| `figdata.py orifice` | `summary.json`, `config.resolved.json`, `mesh_report.json`, `timeseries.csv` | no |
| `figdata.py cases` | `summary.json`, `state.json`, `config.resolved.json`, presence of `slurm-*.out` | yes |
| `paper_numbers.py` `fields`, `energy`, `conservation`, `cd_details`, `studies`, `axial` | `timeseries.csv` with `summary.json`, `config.resolved.json`, `mesh_report.json` | no |
| `paper_numbers.py` `segments`, `gas_spring` | `summary.json` | yes |
| `paper_numbers.py` `cd_ratio` | `summary.json`, `config.resolved.json` | yes |
| `paper_numbers.py` `cost` | `summary.json`, `config.resolved.json`, `mesh_report.json`, `state.json` | yes |
| `verification_numbers.py` `area`, `hyst` | `summary.json`, `config.resolved.json`, `mesh_report.json` | yes |
| `verification_numbers.py` `kinematic`, `branches`, `window`, `gci`, `meshstudy`, `inplane`, `levels`, `further` | `timeseries.csv` with `summary.json`, `config.resolved.json`, `mesh_report.json` | no |
| `verification_numbers.py` `guidance`, `uncertainty` | as `inplane`, and `state.json` | no |
| `bench_table.py` | tables `orifice_benchmark*.csv` only | yes |
| `make_field_specs.py` | `config.resolved.json` (the rendering needs the fields) | specs only |
| `field_ranges.py`, `render_*.py`, `compose_view3d.py` | fields or meshes | no |

The figure data files and tables of the paper computed by these scripts are part of the dataset, so the outputs
marked "no" can be inspected there; regenerating them requires running the cases with the configuration files of
`examples/paper`.
