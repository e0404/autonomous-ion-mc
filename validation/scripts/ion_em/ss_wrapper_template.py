# ruff: noqa: E501
"""V5-SS: the single reviewed host-wrapper source of the ionmc thin-slab sample runs.

Every V5-SS ionmc sample is one host-runner call ``argv == ["python", "-c", <source>]`` with
``<source> == render_wrapper(case, sha, out_dir, scale)`` (the form of the r2 templates T1-T4 of
``r2_wrapper_templates.py``, adapted to ``steps_v8.py ss-sample``). The seal
(``ss_evaluate.py seal``) accepts no other source. The template sets the single-process
environment (``IONMC_RUN_SHA`` = the checked-out HEAD, ``IONMC_RUN_SUITE`` = v5ss), refuses to run
unless HEAD equals the rendered ``sha`` (the implementation commit), creates
``validation/generated/<out_dir>/case<c>`` (published by the host runner under
``validation/generated/<run_id>/<out_dir>/case<c>``), runs ``steps_v8.py ss-sample --case <c>
--k 0 --scale <scale> --seed-base 20631004`` once, writes that subprocess's stdout/stderr verbatim
to ``case<c>-k0.stdout.txt`` / ``.stderr.txt``, prints a short summary of the step document and
exits non-zero if the step failed. :func:`render_wrapper` substitutes only the placeholder line.
"""

from __future__ import annotations

import hashlib
import re

PLACEHOLDER = "case, scale, expected_sha, subdir = {case}, '{scale}', '{sha}', '{subdir}'"
SS_WRAPPER_TEMPLATE = """\
import os, subprocess, sys, json, time, pathlib
case, scale, expected_sha, subdir = {case}, '{scale}', '{sha}', '{subdir}'
sha = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
if sha != expected_sha:
    print(f'HEAD {sha} is not the implementation commit {expected_sha}'); sys.exit(1)
os.environ.update({'PYTHONPATH': '/workspace/src', 'IONMC_CACHE_DIR': '/workspace/.ionmc-cache/ionmc-data', 'IONMC_REQUIRE_DATA': '1', 'IONMC_SINGLE_PROCESS': '1', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1', 'PYTHONUNBUFFERED': '1', 'PYTHONFAULTHANDLER': '1', 'IONMC_RUN_SHA': sha, 'IONMC_RUN_SUITE': 'v5ss'})
out = pathlib.Path(f'validation/generated/{subdir}/case{case}'); out.mkdir(parents=True, exist_ok=True)
cmd = [sys.executable, 'validation/scripts/transport/steps_v8.py', 'ss-sample', '--case', str(case), '--k', '0', '--scale', scale, '--seed-base', '20631004', '--out-dir', str(out)]
t0 = time.time(); p = subprocess.run(cmd, text=True, capture_output=True)
(out / f'case{case}-k0.stdout.txt').write_text(p.stdout); (out / f'case{case}-k0.stderr.txt').write_text(p.stderr)
doc = json.loads(p.stdout.split('#JSON-BEGIN')[1].split('#JSON-END')[0]) if '#JSON-BEGIN' in p.stdout else None
keep = {k2: doc.get(k2) for k2 in ('case', 'k', 'seed', 'seed_base', 'r_index', 'histories', 'npz', 'npz_sha256', 'valid', 'pass', 'wall_s', 'git_sha', 'suite', 'counters', 'n_escaped', 'theta_anchor_rad', 'theta_w_rad')} if doc else None
print('#JSON-BEGIN'); print(json.dumps(keep)); print('#JSON-END')
npz = out / f'ionmc-ss-case{case}-k0.npz'
print(f'case {case} k 0 rc={p.returncode} elapsed={time.time()-t0:.1f}s npz_bytes={npz.stat().st_size if npz.exists() else None}')
if p.returncode != 0 or not doc or not doc.get('pass'):
    print(p.stdout[-3000:]); print(p.stderr[-3000:]); sys.exit(1)
"""
SS_WRAPPER_SHA256 = "021088be16f79292e3a20a420160eff6e463dd599842255114d8b946d01f06eb"
_SCALE = re.compile(r"^\d+(?:\.\d+)?$")
_SHA = re.compile(r"^[0-9a-f]{40}$")
_SUBDIR = re.compile(r"^[a-z0-9][a-z0-9-]*$")


def check_digest() -> None:
    """The template hashes to its recorded sha256 (raises AssertionError otherwise)."""
    got = hashlib.sha256(SS_WRAPPER_TEMPLATE.encode("utf-8")).hexdigest()
    if got != SS_WRAPPER_SHA256:
        raise AssertionError(f"SS wrapper template sha256 {got} != {SS_WRAPPER_SHA256}")


def render_wrapper(case: int, sha: str, out_dir: str, scale: str) -> str:
    """The template with its placeholder line replaced for ``case`` (int), ``sha`` (40-hex
    implementation commit), ``out_dir`` (subdirectory name below ``validation/generated``) and
    ``scale`` (decimal string). Only that line changes; invalid arguments raise ValueError."""
    if type(case) is not int:
        raise ValueError(f"case must be an int: {case!r}")
    if not isinstance(sha, str) or not _SHA.match(sha):
        raise ValueError(f"sha {sha!r} is not a 40-hex commit SHA")
    if not isinstance(out_dir, str) or not _SUBDIR.match(out_dir):
        raise ValueError(f"out_dir {out_dir!r} is not a plain subdirectory name")
    if not isinstance(scale, str) or not _SCALE.match(scale):
        raise ValueError(f"scale {scale!r} is not a decimal string")
    check_digest()
    lines = SS_WRAPPER_TEMPLATE.split("\n")
    if lines.count(PLACEHOLDER) != 1:
        raise ValueError("the template has no unique placeholder line")
    lines[lines.index(PLACEHOLDER)] = (
        f"case, scale, expected_sha, subdir = {case}, '{scale}', '{sha}', '{out_dir}'"
    )
    return "\n".join(lines)
