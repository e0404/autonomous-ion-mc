#!/usr/bin/env bash
# Reproduce every measurement quoted in docs/research/warp-architecture.md.
# Usage (from the repository root, with the project environment synced):
#   bash validation/scripts/warp-architecture/run_all.sh <results-dir>
# Optional: STEPS="06 07" selects steps by their two-digit prefix; STEP_TIMEOUT (s, default 600)
# bounds each command (a timed-out step is archived with exit=124 instead of hanging the suite).
# Each command's stdout+stderr is archived verbatim in <results-dir>/NN-<name>.txt
# together with environment.txt. Raw outputs are committed as measurement provenance.
set -u
out="${1:?results directory required}"
mkdir -p "$out"
here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"
PY="${PYTHON:-python}"
{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "python=$($PY -c 'import sys;print(sys.version.split()[0])')"
  echo "warp=$($PY -c 'import warp;print(warp.__version__)' 2>/dev/null)"
  echo "numpy=$($PY -c 'import numpy;print(numpy.__version__)')"
  echo "cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')"
  echo "logical_cpus=$(nproc)"
  echo "kernel=$(uname -r)"
  echo "git_sha=$(git -C "$here" rev-parse HEAD 2>/dev/null)"
} > "$out/environment.txt"
STEPS="${STEPS:-}"; STEP_TIMEOUT="${STEP_TIMEOUT:-600}"
run() {
  local name="$1"; shift
  if [ -n "$STEPS" ] && ! [[ " $STEPS " == *" ${name%%-*} "* ]]; then return; fi
  echo "== $name: $*"
  { echo "# command: $*"; echo "# step_timeout_s: $STEP_TIMEOUT"; timeout "$STEP_TIMEOUT" "$@"; echo "# exit=$?"; } > "$out/$name.txt" 2>&1
  tail -n 3 "$out/$name.txt"
}
run 01-compile-cold-backward-on  $PY compile_probe.py 1
run 02-compile-cold-backward-off $PY compile_probe.py 0
run 03-toy-f32-1thread  $PY run_toy.py f32 20000 1
run 04-toy-f64-1thread  $PY run_toy.py f64 20000 1
run 05-toy-f32-4threads $PY run_toy.py f32 80000 4
run 06-toy-f32-8threads $PY run_toy.py f32 160000 8
run 07-toy-f32-16threads $PY run_toy.py f32 320000 16
run 08-toy-f32-1thread-repeat $PY run_toy.py f32 20000 1
run 09-precision $PY precision.py
run 10-rng-seed-dupes $PY ../rng/rng_seed_dupes.py
run 11-rng-overlap-1e5x1000 $PY ../rng/rng_overlap.py 100000 1000
run 12-rng-overlap-1e6x1000 $PY ../rng/rng_overlap.py 1000000 1000
run 13-rng-overlap-1e6x2000 $PY ../rng/rng_overlap.py 1000000 2000
run 14-philox-kat-and-python $PY ../rng/philox.py
echo "done: $out"
