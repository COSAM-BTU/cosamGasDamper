#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Entry point:  ./run.sh run <case_dir> --config my.yaml
#               ./run.sh sweep sweep.yaml --root runs --max-cores 48
#               ./run.sh check my.yaml
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="$here/src${PYTHONPATH:+:$PYTHONPATH}"
exec python3 -m damper "$@"
