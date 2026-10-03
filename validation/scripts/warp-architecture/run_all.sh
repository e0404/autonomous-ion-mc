#!/usr/bin/env bash
# Reproduce every measurement quoted in docs/research/warp-architecture.md.
# Usage (from the repository root, with the project environment synced):
#   PYTHON="$PWD/.venv/bin/python" bash validation/scripts/warp-architecture/run_all.sh <results-dir>
# Optional: STEPS="06 07" selects steps by their two-digit prefix; STEP_TIMEOUT (seconds,
# default 600) bounds each command; a timed-out step is archived with exit=124.
# Every step's stdout+stderr is archived verbatim in <results-dir>/NN[-rK]-<name>.txt together
# with environment.txt (versions, hardware, exact git SHA, dirty state and script hashes).
# Timing steps are repeated REPEATS times (default 3); summarize.py reports median and range.
# Any failure (environment probe, hashing, step, summary) makes the script exit non-zero.
set -euo pipefail
out="${1:?results directory required}"
mkdir -p "$out"
here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"
PY="${PYTHON:-python}"
REPEATS="${REPEATS:-3}"
STEPS="${STEPS:-}"
STEP_TIMEOUT="${STEP_TIMEOUT:-600}"
repo="$(git -C "$here" rev-parse --show-toplevel)"

probe() { local label="$1"; shift; local value; value="$("$@")" || { echo "environment probe failed: $label" >&2; exit 1; }; echo "$label=$value"; }
{
  probe date_utc date -u +%Y-%m-%dT%H:%M:%SZ
  probe python "$PY" -c 'import sys;print(sys.version.split()[0])'
  probe warp "$PY" -c 'import warp;print(warp.__version__)'
  probe numpy "$PY" -c 'import numpy;print(numpy.__version__)'
  probe cpu sh -c "grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //'"
  probe logical_cpus nproc
  probe kernel uname -r
  probe git_sha git -C "$here" rev-parse HEAD
  # The results directory being written is excluded; everything else under validation/scripts counts.
  dirty="$(git -C "$repo" status --porcelain -- validation/scripts ':!validation/scripts/warp-architecture/results')"
  echo "scripts_dirty=$([ -n "$dirty" ] && echo yes || echo no)"
  [ -z "$dirty" ] || echo "dirty_entries=$(echo "$dirty" | tr '\n' ';')"
  echo "step_timeout_s=$STEP_TIMEOUT"
  echo "repeats=$REPEATS"
  echo "steps_filter=${STEPS:-all}"
  echo "script_sha256:"
  hashes="$(sha256sum "$here"/*.py "$here"/*.sh "$here"/../rng/*.py)"
  echo "$hashes" | sed 's/^/  /'
} > "$out/environment.txt"

failures=0
run_once() {  # name, command...
  local name="$1"; shift
  echo "== $name: $*"
  local code=0
  {
    echo "# command: $*"
    echo "# step_timeout_s: $STEP_TIMEOUT"
    if timeout "$STEP_TIMEOUT" "$@"; then code=0; else code=$?; fi
    echo "# exit=$code"
  } > "$out/$name.txt" 2>&1 || true
  code="$(sed -n 's/^# exit=//p' "$out/$name.txt" | tail -n 1)"
  [ "$code" = "0" ] || failures=$((failures + 1))
  tail -n 2 "$out/$name.txt"
}
selected() { [ -z "$STEPS" ] || [[ " $STEPS " == *" $1 "* ]]; }
run() {        # single execution: NN-name
  local step="$1" name="$2"; shift 2
  selected "$step" || return 0
  run_once "$step-$name" "$@"
}
run_repeated() {  # REPEATS executions: NN-rK-name
  local step="$1" name="$2"; shift 2
  selected "$step" || return 0
  local k
  for k in $(seq 1 "$REPEATS"); do run_once "$step-r$k-$name" "$@"; done
}
run_repeated 01 compile-cold-backward-on  "$PY" compile_probe.py 1
run_repeated 02 compile-cold-backward-off "$PY" compile_probe.py 0
run_repeated 03 toy-f32-1thread  "$PY" run_toy.py f32 20000 1
run_repeated 04 toy-f64-1thread  "$PY" run_toy.py f64 20000 1
run_repeated 05 toy-f32-4threads "$PY" run_toy.py f32 80000 4
run_repeated 06 toy-f32-8threads "$PY" run_toy.py f32 160000 8
run_repeated 07 toy-f32-16threads "$PY" run_toy.py f32 320000 16
run 08 precision "$PY" precision.py
run 09 rng-seed-dupes "$PY" ../rng/rng_seed_dupes.py
run 10 rng-overlap-1e5x1000 "$PY" ../rng/rng_overlap.py 100000 1000
run 11 rng-overlap-1e6x1000 "$PY" ../rng/rng_overlap.py 1000000 1000
run 12 rng-overlap-1e6x2000 "$PY" ../rng/rng_overlap.py 1000000 2000
run_repeated 13 philox-kat-and-cost "$PY" ../rng/philox.py
"$PY" "$here/summarize.py" "$out" > "$out/SUMMARY.md"
echo "done: $out (failed_or_timed_out_steps=$failures)"
[ "$failures" -eq 0 ]
