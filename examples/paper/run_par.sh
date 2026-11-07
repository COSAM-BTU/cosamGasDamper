#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Run one orifice-benchmark case written by orifice_benchmark.py.
# usage: run_par.sh CASE_DIR NPROC
# OpenFOAM is taken from the current environment; if it is not loaded, the bashrc given in FOAM_BASHRC
# (default: the v2406 package location on Ubuntu) is sourced.
if [ -z "${WM_PROJECT_DIR:-}" ]; then
  source "${FOAM_BASHRC:-/usr/lib/openfoam/openfoam2406/etc/bashrc}" || exit 1
fi
cd "$1" || exit 1
blockMesh > log.blockMesh 2>&1 && checkMesh > log.checkMesh 2>&1 || exit 1
if [ "$2" -gt 1 ]; then
  printf 'FoamFile\n{\n    version 2.0;\n    format ascii;\n    class dictionary;\n    object decomposeParDict;\n}\n\nnumberOfSubdomains %d;\nmethod scotch;\n' "$2" > system/decomposeParDict
  decomposePar -force > log.decomposePar 2>&1 || exit 1
  mpirun --bind-to none -np "$2" rhoPimpleFoam -parallel > log.rhoPimpleFoam 2>&1
else
  rhoPimpleFoam > log.rhoPimpleFoam 2>&1
fi
echo "$(date +%H:%M:%S) $1 $?"
