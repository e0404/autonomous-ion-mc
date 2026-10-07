"""Regenerate the T1 trace fixture ``tests/ionmc/data/trace_baseline_20mev.npz`` (V3-003D, plan
amendment 24 of ``validation/plans/v3-003-acceptance.md``).

Usage (argv only; run on a clean tree at the commit to be recorded)::

    python validation/scripts/transport/make_trace_baseline.py \
        --out tests/ionmc/data/trace_baseline_20mev.npz --code-sha $(git rev-parse HEAD)

The configuration is the one of ``tests/ionmc/test_transport_diagnostic.py`` (``_config`` with the
python backend and the pinned legacy physics options: straggling model
``bohr_gauss_clamped_gamma_v1``, ``short_step_fraction`` 1e-3), imported from the test so that the
generator and the guard cannot drift apart. Fails closed if ``src``, ``pyproject.toml`` or ``uv.lock``
have uncommitted changes or ``--code-sha`` is not ``HEAD``. The 17 trace columns are written together with ``meta`` (a JSON string): generating
commit, ``TransportTables.sha256`` and table identity (including ``range_construction``), the
python, numpy, warp and ionmc versions, the environment flags, the sha256 of the previous fixture
(the file at ``--out`` before it is overwritten) and the comparison with it: row counts, the first
discrete difference, and the maximum absolute difference per column, for the rows before the first
discrete difference and for every history whose discrete columns are identical.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO = Path(__file__).resolve().parents[3]
DISCRETE = ("history", "step", "ix", "iy", "iz", "reason", "blocks", "attempts")


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(REPO), *args], capture_output=True, text=True, check=True
    ).stdout.strip()


def compare_with_previous(
    new: dict[str, np.ndarray], old: dict[str, np.ndarray], columns: tuple[str, ...]
) -> dict[str, Any]:
    """Row counts, first discrete difference and per-column max |difference| (see module doc)."""
    cont = [c for c in columns if c not in DISCRETE]
    n = min(len(new["history"]), len(old["history"]))
    bad = np.zeros(n, dtype=bool)
    for c in DISCRETE:
        bad |= np.asarray(new[c])[:n] != np.asarray(old[c])[:n]
    out: dict[str, Any] = {
        "rows_new": int(len(new["history"])),
        "rows_previous": int(len(old["history"])),
        "first_discrete_difference_row": int(np.argmax(bad)) if bad.any() else None,
    }

    def maxdiff(a: dict[str, np.ndarray], b: dict[str, np.ndarray], sa: Any, sb: Any) -> dict:
        return {
            c: float(np.max(np.abs(np.asarray(a[c])[sa] - np.asarray(b[c])[sb]))) for c in cont
        }

    k = out["first_discrete_difference_row"]
    out["max_abs_diff_rows_before_first_discrete_difference"] = maxdiff(
        new, old, slice(0, n if k is None else k), slice(0, n if k is None else k)
    )
    per_hist: dict[str, Any] = {}
    for h in sorted(set(np.unique(new["history"])) | set(np.unique(old["history"]))):
        mn, mo = np.asarray(new["history"]) == h, np.asarray(old["history"]) == h
        entry: dict[str, Any] = {"rows_new": int(mn.sum()), "rows_previous": int(mo.sum())}
        same = mn.sum() == mo.sum() and all(
            np.array_equal(np.asarray(new[c])[mn], np.asarray(old[c])[mo]) for c in DISCRETE
        )
        entry["discrete_identical"] = bool(same)
        if same:
            entry["max_abs_diff"] = maxdiff(new, old, mn, mo)
        per_hist[str(int(h))] = entry
    out["per_history"] = per_hist
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True, help="fixture to (re)write")
    ap.add_argument("--code-sha", required=True, help="must equal HEAD of a clean tree")
    ap.add_argument(
        "--note", default="", help="free text recorded in meta (for example the discrete flip)"
    )
    args = ap.parse_args(argv)
    head = _git("rev-parse", "HEAD")
    # the code that defines the trace (the package and the dependency pins) must be clean; the
    # generator, the test and the fixture itself may be uncommitted (they are committed with it)
    if _git("status", "--porcelain", "--", "src", "pyproject.toml", "uv.lock"):
        raise SystemExit("refusing to generate with uncommitted changes in src/ or the lock files")
    if args.code_sha != head:
        raise SystemExit(f"--code-sha {args.code_sha} is not HEAD {head}")
    sys.path.insert(0, str(REPO))
    import warp as wp

    import ionmc
    from ionmc.simulation import Simulation
    from ionmc.transport.tally import TRACE_COLUMNS
    from tests.ionmc.test_transport_diagnostic import _config

    out = Path(args.out)
    if not out.is_absolute():
        out = REPO / out
    prev: dict[str, np.ndarray] | None = None
    prev_sha = None
    if out.exists():
        prev_sha = hashlib.sha256(out.read_bytes()).hexdigest()
        with np.load(out) as z:
            prev = {c: z[c] for c in TRACE_COLUMNS}
    res = Simulation(_config("python", legacy=True)).run()
    trace = res.diagnostics["trace"]
    tables = res.effective_config.tables
    arrays = {c: np.asarray(trace[c], dtype=np.float64) for c in TRACE_COLUMNS}
    meta: dict[str, Any] = {
        "generating_commit": head,
        "generator": "validation/scripts/transport/make_trace_baseline.py",
        "backend": "python",
        "config": "tests/ionmc/test_transport_diagnostic.py::_config(python, legacy=True)",
        # informative and machine specific (float64 arrays; libm ulps): the portable identity is
        # table_identity, which the test asserts; this hash is only compared softly
        "transport_tables_sha256": tables.sha256,
        "table_identity": json.loads(json.dumps(_thaw(tables.identity), default=str)),
        "range_construction": tables.identity[0]["range_construction"],
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "warp": wp.config.version,
            "ionmc": getattr(ionmc, "__version__", "unknown"),
        },
        "environment": {
            k: os.environ.get(k)
            for k in ("IONMC_SINGLE_PROCESS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")
        },
        "previous_fixture_sha256": prev_sha,
        "comparison_with_previous": (
            compare_with_previous(arrays, prev, TRACE_COLUMNS) if prev is not None else None
        ),
        "note": args.note,
    }
    np.savez(out, **arrays, meta=np.array(json.dumps(meta, sort_keys=True)))
    print(json.dumps(meta["comparison_with_previous"], indent=1, sort_keys=True))
    print("wrote", out, "sha256", hashlib.sha256(out.read_bytes()).hexdigest())
    return 0


def _thaw(value: Any) -> Any:
    from ionmc.transport.tables import thaw

    return thaw(value)


if __name__ == "__main__":
    raise SystemExit(main())
