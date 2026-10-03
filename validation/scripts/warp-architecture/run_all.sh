#!/usr/bin/env bash
# Reproduce every measurement quoted in docs/research/warp-architecture.md.
# Usage (from the repository root, with the project environment synced):
#   PYTHON="$PWD/.venv/bin/python" bash validation/scripts/warp-architecture/run_all.sh <results-dir>
# Optional: STEPS="06 07" selects steps by their two-digit prefix; STEP_TIMEOUT (seconds,
# default 600) bounds each command; a timed-out step is archived with exit=124.
# Every step's stdout+stderr is archived verbatim in <results-dir>/NN-<name>.txt together
# with environment.txt (versions, hardware, exact git SHA, dirty state and script hashes).
# The script exits non-zero if any step failed or timed out, and summarize.py writes
# SUMMARY.md from the raw files so that no number is transcribed by hand.
set -u
out="${1:?results directory required}"
mkdir -p "$out"
here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"
PY="${PYTHON:-python}"
repo="$(git -C "$here" rev-parse --show-toplevel)"
{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "python=$($PY -c 'import sys;print(sys.version.split()[0])')"
  echo "warp=$($PY -c 'import warp;print(warp.__version__)' 2>/dev/null)"
  echo "numpy=$($PY -c 'import numpy;print(numpy.__version__)')"
  echo "cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')"
  echo "logical_cpus=$(nproc)"
  echo "kernel=$(uname -r)"
  echo "git_sha=$(git -C "$here" rev-parse HEAD)"
  # The results directory being written is excluded; everything else under validation/scripts counts.
  dirty="$(git -C "$repo" status --porcelain -- validation/scripts ':!validation/scripts/warp-architecture/results')"
  echo "scripts_dirty=$([ -n "$dirty" ] && echo yes || echo no)"
  [ -n "$dirty" ] && echo "dirty_entries=$(echo "$dirty" | tr '\n' ';')"
  echo "step_timeout_s=${STEP_TIMEOUT:-600}"
  echo "steps_filter=${STEPS:-all}"
  echo "script_sha256:"
  sha256sum "$here"/*.py "$here"/*.sh "$here"/../rng/*.py | sed 's/^/  /'
} > "$out/environment.txt"
STEPS="${STEPS:-}"; STEP_TIMEOUT="${STEP_TIMEOUT:-600}"
failures=0
run() {
  local name="$1"; shift
  if [ -n "$STEPS" ] && ! [[ " $STEPS " == *" ${name%%-*} "* ]]; then return; fi
  echo "== $name: $*"
  {
    echo "# command: $*"
    echo "# step_timeout_s: $STEP_TIMEOUT"
    timeout "$STEP_TIMEOUT" "$@"
    echo "# exit=$?"
  } > "$out/$name.txt" 2>&1
  local code
  code="$(sed -n 's/^# exit=//p' "$out/$name.txt" | tail -n 1)"
  [ "$code" = "0" ] || failures=$((failures + 1))
  tail -n 2 "$out/$name.txt"
}
run 01-compile-cold-backward-on  "$PY" compile_probe.py 1
run 02-compile-cold-backward-off "$PY" compile_probe.py 0
run 03-toy-f32-1thread  "$PY" run_toy.py f32 20000 1
run 04-toy-f64-1thread  "$PY" run_toy.py f64 20000 1
run 05-toy-f32-4threads "$PY" run_toy.py f32 80000 4
run 06-toy-f32-8threads "$PY" run_toy.py f32 160000 8
run 07-toy-f32-16threads "$PY" run_toy.py f32 320000 16
run 08-toy-f32-1thread-repeat "$PY" run_toy.py f32 20000 1
run 09-precision "$PY" precision.py
run 10-rng-seed-dupes "$PY" ../rng/rng_seed_dupes.py
run 11-rng-overlap-1e5x1000 "$PY" ../rng/rng_overlap.py 100000 1000
run 12-rng-overlap-1e6x1000 "$PY" ../rng/rng_overlap.py 1000000 1000
run 13-rng-overlap-1e6x2000 "$PY" ../rng/rng_overlap.py 1000000 2000
run 14-philox-kat-and-python "$PY" ../rng/philox.py
"$PY" "$here/summarize.py" "$out" > "$out/SUMMARY.md"
echo "done: $out (failed_or_timed_out_steps=$failures)"
[ "$failures" -eq 0 ]
