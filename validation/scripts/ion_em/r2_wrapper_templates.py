# ruff: noqa: E501
"""V3-008 revision 2: the four reviewed chunk-wrapper sources of the sealed ionmc sample runs.

Every sample chunk was produced by a host-runner call ``argv == ["python", "-c", <source>]``
(Codex REVIEW-84382c07). The evaluator (``r2_io.sample_invocation``) accepts ONLY a source that is
byte-identical to one of the templates below rendered for the run's sealed chunks; any other signed
source fails closed. The templates are the verbatim sources of the committed signed host-run
records (``validation/results/transport/v8r2-samples/host-records/``) with the single
case/k/scale assignment line replaced by a placeholder line (:data:`PLACEHOLDER_K` /
:data:`PLACEHOLDER_KS`); :func:`render` substitutes only that line (the sources contain f-string
braces, so ``str.format`` is never applied to a whole source).

What every template does: it sets the single-process lv8 environment (``IONMC_RUN_SHA`` = the
checked-out HEAD, ``IONMC_RUN_SUITE`` = lv8), creates the out directory
``validation/generated/<subdir>/case<c>`` (published by the host runner to
``validation/generated/<run_id>/<subdir>/case<c>``), runs
``steps_v8.py td-tl-sample --case <c> --k <k> --scale <scale> [--seed-base 20621004] --out-dir
<out>`` as a subprocess once per k, writes that subprocess's stdout and stderr verbatim to
``case<c>-k<k>.stdout.txt`` / ``.stderr.txt``, prints a short summary of the step document and
exits non-zero if a step's return code is non-zero, its document is missing or not ``pass``
(the host record then says ``succeeded: false``). It writes nothing else.

- ``T1`` (1 run, main case 1 k 1): one k; assignment line after ``os.environ.update``.
- ``T2`` (30 runs, main cases 1-5): as T1 with the assignment line second.
- ``T3`` (45 runs, main cases 6-12): ``ks`` list, one step per k (stops early, exit 1, if a step
  fails; skips the remaining k after 820 s, so a sealed k without its step output fails closed in
  the evaluator).
- ``T4`` (13 runs, held-out cases 10-12): as T2 with ``--seed-base 20621004`` (amendment
  R2-STOP3) and the out directory ``r2-samples-heldout``.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass

PLACEHOLDER_K = "case, k, scale = {case}, {k}, '{scale}'"
PLACEHOLDER_KS = "case, ks, scale = {case}, {ks}, '{scale}'"
_SCALE = re.compile(r"^\d+(?:\.\d+)?$")


@dataclass(frozen=True)
class WrapperTemplate:
    id: str
    source: str  # verbatim, with the placeholder assignment line
    sha256: str  # sha256 of ``source`` (UTF-8)
    kind: str  # seal kind allowed: "main" | "heldout"
    subdir: str  # out directory below validation/generated
    seed_base: int  # the --seed-base it passes (20581004 = steps_v8.py default when absent)
    multi_k: bool  # True: ``ks`` list (PLACEHOLDER_KS); False: one ``k`` (PLACEHOLDER_K)
    example_run: str  # a committed host record with this source


T1_SOURCE = """\
import os, subprocess, sys, json, time, pathlib
sha = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
os.environ.update({'PYTHONPATH': '/workspace/src', 'IONMC_CACHE_DIR': '/workspace/.ionmc-cache/ionmc-data', 'IONMC_REQUIRE_DATA': '1', 'IONMC_SINGLE_PROCESS': '1', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1', 'PYTHONUNBUFFERED': '1', 'PYTHONFAULTHANDLER': '1', 'IONMC_RUN_SHA': sha, 'IONMC_RUN_SUITE': 'lv8'})
case, k, scale = {case}, {k}, '{scale}'
out = pathlib.Path(f'validation/generated/r2-samples/case{case}'); out.mkdir(parents=True, exist_ok=True)
cmd = [sys.executable, 'validation/scripts/transport/steps_v8.py', 'td-tl-sample', '--case', str(case), '--k', str(k), '--scale', scale, '--out-dir', str(out)]
t0 = time.time(); p = subprocess.run(cmd, text=True, capture_output=True)
(out / f'case{case}-k{k}.stdout.txt').write_text(p.stdout); (out / f'case{case}-k{k}.stderr.txt').write_text(p.stderr)
doc = json.loads(p.stdout.split('#JSON-BEGIN')[1].split('#JSON-END')[0]) if '#JSON-BEGIN' in p.stdout else None
keep = {k2: doc.get(k2) for k2 in ('case', 'k', 'seed', 'r_index', 'histories', 'npz', 'npz_sha256', 'valid', 'pass', 'wall_s', 'git_sha', 'suite', 'counters')} if doc else None
print('#JSON-BEGIN'); print(json.dumps(keep)); print('#JSON-END')
npz = out / f'ionmc-case{case}-k{k}.npz'
print(f'case {case} k {k} rc={p.returncode} elapsed={time.time()-t0:.1f}s npz_bytes={npz.stat().st_size if npz.exists() else None}')
if p.returncode != 0 or not doc or not doc.get('pass'):
    print(p.stdout[-3000:]); print(p.stderr[-3000:]); sys.exit(1)
"""
T1_SHA256 = "48a902ec74e7f62c0edc85747e7f4c507cbd65fc9fa74a7557fd83225f8e3b9b"

T2_SOURCE = """\
import os, subprocess, sys, json, time, pathlib
case, k, scale = {case}, {k}, '{scale}'
sha = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
os.environ.update({'PYTHONPATH': '/workspace/src', 'IONMC_CACHE_DIR': '/workspace/.ionmc-cache/ionmc-data', 'IONMC_REQUIRE_DATA': '1', 'IONMC_SINGLE_PROCESS': '1', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1', 'PYTHONUNBUFFERED': '1', 'PYTHONFAULTHANDLER': '1', 'IONMC_RUN_SHA': sha, 'IONMC_RUN_SUITE': 'lv8'})
out = pathlib.Path(f'validation/generated/r2-samples/case{case}'); out.mkdir(parents=True, exist_ok=True)
cmd = [sys.executable, 'validation/scripts/transport/steps_v8.py', 'td-tl-sample', '--case', str(case), '--k', str(k), '--scale', scale, '--out-dir', str(out)]
t0 = time.time(); p = subprocess.run(cmd, text=True, capture_output=True)
(out / f'case{case}-k{k}.stdout.txt').write_text(p.stdout); (out / f'case{case}-k{k}.stderr.txt').write_text(p.stderr)
doc = json.loads(p.stdout.split('#JSON-BEGIN')[1].split('#JSON-END')[0]) if '#JSON-BEGIN' in p.stdout else None
keep = {k2: doc.get(k2) for k2 in ('case', 'k', 'seed', 'r_index', 'histories', 'npz', 'npz_sha256', 'valid', 'pass', 'wall_s', 'git_sha', 'suite', 'counters')} if doc else None
print('#JSON-BEGIN'); print(json.dumps(keep)); print('#JSON-END')
npz = out / f'ionmc-case{case}-k{k}.npz'
print(f'case {case} k {k} rc={p.returncode} elapsed={time.time()-t0:.1f}s npz_bytes={npz.stat().st_size if npz.exists() else None}')
if p.returncode != 0 or not doc or not doc.get('pass'):
    print(p.stdout[-3000:]); print(p.stderr[-3000:]); sys.exit(1)
"""
T2_SHA256 = "2f3a3c376eefe4c737c71b37823e45cabd190825dfaaa79e5d1bd92bae3ca68a"

T3_SOURCE = """\
import os, subprocess, sys, json, time, pathlib
case, ks, scale = {case}, {ks}, '{scale}'
sha = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
os.environ.update({'PYTHONPATH': '/workspace/src', 'IONMC_CACHE_DIR': '/workspace/.ionmc-cache/ionmc-data', 'IONMC_REQUIRE_DATA': '1', 'IONMC_SINGLE_PROCESS': '1', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1', 'PYTHONUNBUFFERED': '1', 'PYTHONFAULTHANDLER': '1', 'IONMC_RUN_SHA': sha, 'IONMC_RUN_SUITE': 'lv8'})
out = pathlib.Path(f'validation/generated/r2-samples/case{case}'); out.mkdir(parents=True, exist_ok=True)
T0 = time.time(); rc_all = 0
for i, k in enumerate(ks):
    if i > 0 and time.time() - T0 > 820:
        print(f'SKIP k={k}: elapsed {time.time()-T0:.0f}s leaves no room'); break
    cmd = [sys.executable, 'validation/scripts/transport/steps_v8.py', 'td-tl-sample', '--case', str(case), '--k', str(k), '--scale', scale, '--out-dir', str(out)]
    t0 = time.time(); p = subprocess.run(cmd, text=True, capture_output=True)
    (out / f'case{case}-k{k}.stdout.txt').write_text(p.stdout); (out / f'case{case}-k{k}.stderr.txt').write_text(p.stderr)
    doc = json.loads(p.stdout.split('#JSON-BEGIN')[1].split('#JSON-END')[0]) if '#JSON-BEGIN' in p.stdout else None
    keep = {k2: doc.get(k2) for k2 in ('case', 'k', 'seed', 'r_index', 'histories', 'npz', 'npz_sha256', 'valid', 'pass', 'wall_s', 'git_sha', 'suite', 'counters')} if doc else None
    print('#JSON-BEGIN'); print(json.dumps(keep)); print('#JSON-END')
    npz = out / f'ionmc-case{case}-k{k}.npz'
    print(f'case {case} k {k} rc={p.returncode} elapsed={time.time()-t0:.1f}s npz_bytes={npz.stat().st_size if npz.exists() else None}', flush=True)
    if p.returncode != 0 or not doc or not doc.get('pass'):
        print(p.stdout[-3000:]); print(p.stderr[-3000:]); rc_all = 1; break
sys.exit(rc_all)
"""
T3_SHA256 = "9a0f9cd18736891da354ffd0cff38f60de5c0144bb8713d5e98ec8fd7a41fcc8"

T4_SOURCE = """\
import os, subprocess, sys, json, time, pathlib
case, k, scale = {case}, {k}, '{scale}'
sha = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
os.environ.update({'PYTHONPATH': '/workspace/src', 'IONMC_CACHE_DIR': '/workspace/.ionmc-cache/ionmc-data', 'IONMC_REQUIRE_DATA': '1', 'IONMC_SINGLE_PROCESS': '1', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1', 'PYTHONUNBUFFERED': '1', 'PYTHONFAULTHANDLER': '1', 'IONMC_RUN_SHA': sha, 'IONMC_RUN_SUITE': 'lv8'})
out = pathlib.Path(f'validation/generated/r2-samples-heldout/case{case}'); out.mkdir(parents=True, exist_ok=True)
cmd = [sys.executable, 'validation/scripts/transport/steps_v8.py', 'td-tl-sample', '--case', str(case), '--k', str(k), '--scale', scale, '--seed-base', '20621004', '--out-dir', str(out)]
t0 = time.time(); p = subprocess.run(cmd, text=True, capture_output=True)
(out / f'case{case}-k{k}.stdout.txt').write_text(p.stdout); (out / f'case{case}-k{k}.stderr.txt').write_text(p.stderr)
doc = json.loads(p.stdout.split('#JSON-BEGIN')[1].split('#JSON-END')[0]) if '#JSON-BEGIN' in p.stdout else None
keep = {k2: doc.get(k2) for k2 in ('case', 'k', 'seed', 'seed_base', 'r_index', 'histories', 'npz', 'npz_sha256', 'valid', 'pass', 'wall_s', 'git_sha', 'suite', 'counters', 'end_of_range_model')} if doc else None
print('#JSON-BEGIN'); print(json.dumps(keep)); print('#JSON-END')
npz = out / f'ionmc-case{case}-k{k}.npz'
print(f'case {case} k {k} rc={p.returncode} elapsed={time.time()-t0:.1f}s npz_bytes={npz.stat().st_size if npz.exists() else None}')
if p.returncode != 0 or not doc or not doc.get('pass'):
    print(p.stdout[-3000:]); print(p.stderr[-3000:]); sys.exit(1)
"""
T4_SHA256 = "6b86c32b12c6fcad2d72a671c68522a1cd99224178d8eada9ed41a8589acf9ea"


TEMPLATES: tuple[WrapperTemplate, ...] = (
    WrapperTemplate(
        "T1",
        T1_SOURCE,
        T1_SHA256,
        "main",
        "r2-samples",
        20581004,
        False,
        "RUN-20261108T051423Z-5c33a3e4",
    ),
    WrapperTemplate(
        "T2",
        T2_SOURCE,
        T2_SHA256,
        "main",
        "r2-samples",
        20581004,
        False,
        "RUN-20261108T053031Z-29586a53",
    ),
    WrapperTemplate(
        "T3",
        T3_SOURCE,
        T3_SHA256,
        "main",
        "r2-samples",
        20581004,
        True,
        "RUN-20261108T152342Z-d0231a16",
    ),
    WrapperTemplate(
        "T4",
        T4_SOURCE,
        T4_SHA256,
        "heldout",
        "r2-samples-heldout",
        20621004,
        False,
        "RUN-20261111T032108Z-0f2c091b",
    ),
)


def render(template: WrapperTemplate, case: int, ks: list[int], scale: str) -> str:
    """``template.source`` with its one placeholder line replaced by the assignment for ``case``,
    ``ks`` (T1/T2/T4: exactly one k, rendered as an int; T3: ``repr`` of the sorted list, e.g.
    ``[1, 2]``) and ``scale`` (a decimal string, rendered quoted). Only that line changes."""
    if type(case) is not int or not ks or any(type(k) is not int for k in ks):
        raise ValueError(f"case/ks must be ints: {case!r} {ks!r}")
    if not isinstance(scale, str) or not _SCALE.match(scale):
        raise ValueError(f"scale {scale!r} is not a decimal string")
    if template.multi_k:
        placeholder = PLACEHOLDER_KS
        line = f"case, ks, scale = {case}, {sorted(ks)!r}, '{scale}'"
    else:
        if len(ks) != 1:
            raise ValueError(f"template {template.id} runs exactly one k, not {ks!r}")
        placeholder = PLACEHOLDER_K
        line = f"case, k, scale = {case}, {ks[0]}, '{scale}'"
    lines = template.source.split("\n")
    if lines.count(placeholder) != 1:
        raise ValueError(f"template {template.id} has no unique placeholder line")
    lines[lines.index(placeholder)] = line
    return "\n".join(lines)


def check_digests() -> None:
    """Every template source hashes to its recorded sha256 (raises AssertionError otherwise)."""
    for t in TEMPLATES:
        got = hashlib.sha256(t.source.encode("utf-8")).hexdigest()
        if got != t.sha256:
            raise AssertionError(f"wrapper template {t.id}: sha256 {got} != {t.sha256}")
