#!/usr/bin/env bash
# Reproduce every measurement quoted in docs/research/warp-architecture.md.
# Usage (from the repository root, with the project environment synced):
#   PYTHON="$PWD/.venv/bin/python" bash validation/scripts/warp-architecture/run_all.sh <new-results-dir>
# The results directory must not exist yet: every archive is one complete execution, so
# evidence from different executions or SHAs can never be mixed. Optional environment:
# REPEATS (timing repeats, default 3) and STEP_TIMEOUT (seconds per command, default 600;
# a timed-out step is archived with exit=124).
# Every step's stdout+stderr is archived verbatim in <results-dir>/NN[-rK]-<name>.txt with a
# per-step provenance header (git SHA, UTC start time), next to environment.txt (versions,
# hardware, exact git SHA, dirty state and script hashes). summarize.py then writes SUMMARY.md.
# Any failure (environment probe, hashing, archive write, step, summary) exits non-zero.
set -euo pipefail
out="${1:?new results directory required}"
if [ -e "$out" ]; then
  echo "refusing to reuse an existing results directory: $out" >&2
  exit 1
fi
mkdir -p "$out"
here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"
PY="${PYTHON:-python}"
REPEATS="${REPEATS:-3}"
STEP_TIMEOUT="${STEP_TIMEOUT:-600}"
repo="$(git -C "$here" rev-parse --show-toplevel)"
sha="$(git -C "$here" rev-parse HEAD)"

probe() {  # label command... : the command must succeed and print a non-empty value
  local label="$1"; shift
  local value
  value="$("$@")" || { echo "environment probe failed: $label" >&2; exit 1; }
  [ -n "$value" ] || { echo "environment probe returned nothing: $label" >&2; exit 1; }
  echo "$label=$value"
}
cpu_model() { local line; line="$(grep -m1 'model name' /proc/cpuinfo)" || return 1; echo "${line#*: }"; }
{
  probe date_utc date -u +%Y-%m-%dT%H:%M:%SZ
  probe python "$PY" -c 'import sys;print(sys.version.split()[0])'
  probe warp "$PY" -c 'import warp;print(warp.__version__)'
  probe numpy "$PY" -c 'import numpy;print(numpy.__version__)'
  probe cpu cpu_model
  probe logical_cpus nproc
  probe kernel uname -r
  echo "git_sha=$sha"
  # The results directory being written is excluded; everything else under validation/scripts counts.
  dirty="$(git -C "$repo" status --porcelain -- validation/scripts ':!validation/scripts/warp-architecture/results')"
  echo "scripts_dirty=$([ -n "$dirty" ] && echo yes || echo no)"
  [ -z "$dirty" ] || echo "dirty_entries=$(echo "$dirty" | tr '\n' ';')"
  echo "step_timeout_s=$STEP_TIMEOUT"
  echo "repeats=$REPEATS"
  echo "script_sha256:"
  hashes="$(sha256sum "$here"/*.py "$here"/*.sh "$here"/../rng/*.py)"
  echo "$hashes" | sed 's/^/  /'
} > "$out/environment.txt"

failures=0
run_once() {  # name, command...
  local name="$1"; shift
  local file="$out/$name.txt"
  if [ -e "$file" ]; then echo "duplicate step name: $name" >&2; exit 1; fi
  echo "== $name: $*"
  local code=0 started
  started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  # Header lines are written with set -e active: an unwritable archive aborts the run.
  {
    echo "# command: $*"
    echo "# git_sha: $sha"
    echo "# started_utc: $started"
    echo "# step_timeout_s: $STEP_TIMEOUT"
  } > "$file"
  # The measured command itself may fail; its exit status is recorded, not masked.
  set +e
  timeout "$STEP_TIMEOUT" "$@" >> "$file" 2>&1
  code=$?
  set -e
  echo "# exit=$code" >> "$file"
  # Verify that this execution produced the trailer that is parsed below.
  [ "$(tail -n 1 "$file")" = "# exit=$code" ] || { echo "archive verification failed: $file" >&2; exit 1; }
  [ "$code" -eq 0 ] || failures=$((failures + 1))
  tail -n 2 "$file"
}
run() { run_once "$1-$2" "${@:3}"; }
run_repeated() {
  local step="$1" name="$2" k; shift 2
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
