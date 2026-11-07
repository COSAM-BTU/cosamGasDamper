# Studies of the accompanying article

This directory holds everything needed to repeat the simulations of the article
"An open-source multi-mesh OpenFOAM framework for the automated simulation of
large-stroke gas-filled piston–cylinder dampers with multiple orifices":

* one configuration (sweep) file per study,
* `orifice_benchmark.py` and `run_par.sh`, the single-orifice benchmark against
  the long-orifice data of Lichtarowicz et al. (1965),
* `analysis/`, the scripts that turn the case records into the figures and
  tables of the article (see `analysis/README.md`).

The case summaries, segment tables, mesh reports, resolved configurations and
stored cross-sections of all cases, and the derived tables, are archived as a
separate Zenodo dataset (see the main `README.md`).

## Studies

Every file is a sweep file for `run.sh sweep`. The run-directory names in the
second column are the ones the analysis scripts expect.

| File | Run directory | Cases | Content |
|---|---|---|---|
| `01_matrix.yaml` | `matrix` | 30 | orifice count `N = 2...10` x mean speed `0.025...4 m/s`, medium mesh |
| `02_mesh.yaml`, `07_mesh_extra.yaml` | `studies` | 6 | global mesh refinement, `N = 2`, 1 m/s |
| `03_courant.yaml` ... `06_sinusoidal.yaml`, `08_segments.yaml` | `studies` | 9 | Courant number, pre-charge pressure, laminar flow, sinusoidal motion, sub-meshes per stroke (medium mesh) |
| `09_orifice_refinement.yaml`, `10_orifice_count.yaml` | `studies` | 6 | local refinement of the orifice cross-sections (curvature size `C`) |
| `11_matrix_C8.yaml` | `matrixC8` | 30 | matrix on the refined cross-section (`C = 8`) with axial cells of 1 mm |
| `12_studies_C8.yaml` | `studiesC8` | 10 | sensitivity studies on that mesh |
| `13_further_C8.yaml` | `furtherC8` | 9 | Courant number and further refinement levels on that mesh |
| `14_matrix_ref.yaml` | `matrixC8ax` | 30 | matrix on the reference mesh (`C = 8`, axial cells of 0.1 mm in and at the piston) |
| `15_studies_ref.yaml` | `studiesC8ax` | 27 | verification and sensitivity studies on the reference mesh |
| `16_sections_ref.yaml` | `sections` | 10 | reference mesh for `N = 1...10`, mesh stage only (mesh figures) |
| `17_further_ref.yaml` | `furtherC8ax` | 18 | further in-plane levels, orifice-area option, conservative transfer, equal total orifice area |

The core counts (`run.n_procs`) in the files are those used for the published
runs; adjust them to the machine at hand. On the reference mesh a case costs
about 100 to 1400 core-hours, most at the lowest speeds; the medium-mesh matrix
costs about 340 core-hours in total.

## Running the studies

Start every sweep from one working directory: the run directories are created
there, and the stored cross-sections that a file reuses (`mesh.section_med`) are
paths relative to it.

```bash
mkdir -p runs && cd runs
../run.sh sweep ../examples/paper/01_matrix.yaml --root matrix --max-cores 48
../run.sh sweep ../examples/paper/0[2-9]_*.yaml ../examples/paper/10_orifice_count.yaml \
    --root studies --max-cores 48
../run.sh sweep ../examples/paper/11_matrix_C8.yaml --root matrixC8 --max-cores 48
../run.sh sweep ../examples/paper/12_studies_C8.yaml --root studiesC8 --max-cores 48
../run.sh sweep ../examples/paper/13_further_C8.yaml --root furtherC8 --max-cores 48
../run.sh sweep ../examples/paper/14_matrix_ref.yaml --root matrixC8ax --max-cores 48
../run.sh sweep ../examples/paper/15_studies_ref.yaml --root studiesC8ax --max-cores 48
../run.sh sweep ../examples/paper/17_further_ref.yaml --root furtherC8ax --max-cores 48

# mesh stage only for the mesh figures
../run.sh sweep ../examples/paper/16_sections_ref.yaml --root sections --dry-run
for c in sections/sec_N*; do ../run.sh run "$c" --only mesh; done

# one table of all case summaries of a study
../run.sh collect matrixC8ax --out matrixC8ax.csv
```

The order matters where a file reuses the cross-section of another study, so
that only the studied parameter differs between a case and its reference:
`14`, `15` and `17` reuse cross-sections of `matrixC8`; `15` also reuses those
of `studies` (`orif_C6`, `orif_C12`), `studiesC8` and `furtherC8`; and `17` reuses
some of `studiesC8ax` and of its own earlier cases (those are listed first in
the file). Because the quadrilateral cross-section mesher is not
bit-reproducible, a cross-section generated anew differs from the published one
by a few per cent in the cell count, and the damping force by about 1 % on the
medium mesh. To repeat a study on exactly the published meshes, copy the
stored cross-sections of the dataset into the working directory first; they
are arranged in the same `<run directory>/<case>/mesh/section.med` layout:

```bash
rsync -a --include='*/' --include='section.med' --exclude='*' <dataset>/cases/ runs/
```

Copy only the cross-sections, not the whole case records: their `state.json`
marks the stages as finished, and the framework would skip the cases. The
files `13`, `14`, `15` and `17` then reuse the published cross-sections. The
other files generate new cross-sections, which overwrite the copied ones of
their run directory; to keep a published one, set `mesh.section_med` of the
case to it (for example `matrixC8/N2_v1.0/mesh/section.med`, copied to a
separate folder first).

## Orifice benchmark

```bash
python3 ../examples/paper/orifice_benchmark.py setup orificeBench     # writes all cases
for c in orificeBench/*/; do ../examples/paper/run_par.sh "$c" 8; done
python3 ../examples/paper/orifice_benchmark.py collect orificeBench > orifice_benchmark.csv
```

`run_par.sh` uses the OpenFOAM environment of the shell, or sources the
bashrc given in `FOAM_BASHRC`. The runs of mesh B that did not settle with the
default Courant limit were repeated with `maxCo 0.5` in `system/controlDict`
(all other settings unchanged); their coefficients are collected into
`orifice_benchmark_T2.csv` in the same way and replace the failed runs in the
benchmark table of the article.
