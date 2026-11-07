# Changelog

## 1.0.0 (2026-11-07)

First public release, accompanying the article "An open-source multi-mesh OpenFOAM framework for the automated
simulation of large-stroke gas-filled piston–cylinder dampers with multiple orifices".

- Configuration-driven driver (`run.sh`: `check`, `run`, `sweep`, `collect`) with restartable stages
  (mesh, solve, post, render) and a sweep runner under a core budget.
- Parametric SALOME/Gmsh generator for `N >= 1` equally spaced orifices: one quasi-structured quadrilateral
  cross-section swept into hexahedral sub-meshes; optional local refinement of the orifices and preservation of
  the orifice area.
- Segmented multi-mesh scheme for large strokes with field transfer by consistent interpolation or by the
  volume-weighted, mass-conserving `mapFieldsPar` transfer; gas mass and volume checked at every time step.
- Post-processing into time series, segment tables and a force and energy summary; headless pvbatch renderings.
- Tutorial case, unit tests, the configuration files of all studies of the article, the single-orifice benchmark
  and the analysis and plotting scripts of the article.
