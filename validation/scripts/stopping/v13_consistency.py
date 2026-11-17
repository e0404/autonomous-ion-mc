"""V3-013 step S4: consistency checks of the model ``icru90-water-v1`` (report-only aggregates).

Usage::

    uv run python validation/scripts/stopping/v13_consistency.py \
        [--cache-dir CACHE] [--out-dir validation/results/stopping/v3-013]

Needs the cached objects ``geant4-icru90-stopping-11.4.2``, ``nist-pstar-water-2005`` and
``nist-astar-water-2005`` (read through the verified provenance cache; nothing is downloaded).
Writes two JSON files holding AGGREGATES ONLY (counts, windows, max and RMS statistics; no table
value, no per-node value, no stopping power or range, no ratio of a single node):

* ``cc-int.json``: rows CC-INT (interpolation error between nodes, log-log versus a monotone cubic
  in log-log space, both built from the array only; max and RMS per decade), the report-only parts
  of CC-NODE (chord error of the fine table at the ICRU nodes, max and RMS) and CC-I2 (start-range
  offset of the model against the NIST start of the I2 construction, one number per projectile;
  measured maximum relative deviation of the node arrays R(E_i) - R(E_0), the gated quantity),
  and CC-SPAN (node counts and coverage);
* ``v-c3.json``: V-C3, ``icru90-water-v1`` against PSTAR (proton) and ASTAR (alpha), stopping power
  and CSDA range over the common window (max and RMS relative deviation).

These are consistency results for implementation correctness, never evidence of physical accuracy
(plan ``validation/plans/v3-013-acceptance.md``, sections A.5, A.6 and B). The gated assertions of
CC-NODE (at the fine nodes) and CC-I2 are in ``tests/ionmc/test_stopping_icru90_model.py``, which
imports the functions of this module. The JSON records the SHA-256 of this script and of its
input set (the three cache objects, ``src/ionmc/data/icru90.py``,
``src/ionmc/physics/stopping_tabulated.py``, ``src/ionmc/physics/stopping_icru90.py``).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

sys.dont_write_bytecode = True

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
WORKTREE = HERE.parents[2]
sys.path.insert(0, str(WORKTREE / "src"))

from ionmc.data import cache  # noqa: E402
from ionmc.data.icru90 import Icru90Water, load_icru90_water  # noqa: E402
from ionmc.data.nist_star import StarTable, load_star_table  # noqa: E402
from ionmc.materials import water  # noqa: E402
from ionmc.physics.projectiles import ALPHA, PROTON, Projectile  # noqa: E402
from ionmc.physics.stopping_icru90 import Icru90WaterSource, icru90_start_range  # noqa: E402
from ionmc.physics.stopping_icru90 import load_icru90_model as _load_model  # noqa: E402

ICRU_ID = "geant4-icru90-stopping-11.4.2"
PSTAR_ID, ASTAR_ID = "nist-pstar-water-2005", "nist-astar-water-2005"
PROJECTILES: tuple[Projectile, ...] = (PROTON, ALPHA)
INPUT_FILES = (
    "src/ionmc/data/icru90.py",
    "src/ionmc/physics/stopping_tabulated.py",
    "src/ionmc/physics/stopping_icru90.py",
)
DECADES = ((1.0, 10.0), (10.0, 100.0), (100.0, 1000.0))  # MeV/u, clipped to the model window
SAMPLES_PER_DECADE = 400
OUT_DIR = WORKTREE / "validation/results/stopping/v3-013"
NOTE = (
    "Aggregates only: no ICRU 90, PSTAR or ASTAR table value, per-node value or per-energy ratio "
    "is written. Consistency results for implementation correctness, not evidence of physical "
    "accuracy."
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_inputs(cache_dir: Path) -> tuple[Icru90WaterSource, Icru90Water, StarTable, StarTable]:
    """The model, its parsed arrays and the PSTAR and ASTAR tables, all from the verified cache."""
    model = _load_model(str(cache_dir))
    data = load_icru90_water(cache.verify(ICRU_ID, cache_dir))
    pstar = load_star_table(cache.verify(PSTAR_ID, cache_dir))
    astar = load_star_table(cache.verify(ASTAR_ID, cache_dir))
    return model, data, pstar, astar


def arrays_per_u(data: Icru90Water, proj: Projectile) -> tuple[np.ndarray, np.ndarray]:
    """Parsed ICRU 90 array of ``proj`` as (energy per nucleon [MeV/u], S [MeV cm2/g])."""
    if proj.name == "proton":
        return data.proton_energy_mev / proj.a, data.proton_stopping
    return data.alpha_energy_mev / proj.a, data.alpha_stopping


def _loglog(e: np.ndarray, e_u: np.ndarray, s: np.ndarray) -> np.ndarray:
    return np.exp(np.interp(np.log(e), np.log(e_u), np.log(s)))


def _stats(dev: np.ndarray) -> dict[str, Any]:
    return {
        "n_points": int(dev.size),
        "max_abs_relative_deviation": float(np.max(np.abs(dev))),
        "rms_relative_deviation": float(np.sqrt(np.mean(dev**2))),
    }


# ---- CC-NODE ----


def cc_node_fine(model: Icru90WaterSource, data: Icru90Water, proj: Projectile) -> float:
    """Maximum relative deviation of the transport table's S at the fine nodes from the log-log
    interpolant of the parsed array evaluated there (gated at 1e-12)."""
    t = model.table(water(78.0), proj)
    e_u, s = arrays_per_u(data, proj)
    ref = _loglog(t.energy_per_u, e_u, s)
    return float(np.max(np.abs(t.s_el_mass / ref - 1.0)))


def cc_node_chord(model: Icru90WaterSource, data: Icru90Water, proj: Projectile) -> dict[str, Any]:
    """Chord error of the fine table at the ICRU nodes inside its window (off the fine grid):
    ``table.stopping_at(node) / array_value - 1``; max and RMS only."""
    t = model.table(water(78.0), proj)
    e_u, s = arrays_per_u(data, proj)
    lo, hi = float(t.energy_per_u[0]), float(t.energy_per_u[-1])
    m = (e_u >= lo * (1 + 1e-12)) & (e_u <= hi * (1 - 1e-12))
    return _stats(t.stopping_at(e_u[m]) / s[m] - 1.0)


# ---- CC-INT ----


def pchip_slopes(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Fritsch-Carlson monotone cubic Hermite slopes (end slopes by the shape-preserving
    three-point formula)."""
    h = np.diff(x)
    delta = np.diff(y) / h
    n = x.size
    d = np.zeros(n)
    for k in range(1, n - 1):
        if delta[k - 1] * delta[k] > 0.0:
            w1 = 2.0 * h[k] + h[k - 1]
            w2 = h[k] + 2.0 * h[k - 1]
            d[k] = (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k])

    def end(h0: float, h1: float, d0: float, d1: float) -> float:
        v = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
        if np.sign(v) != np.sign(d0):
            return 0.0
        if np.sign(d0) != np.sign(d1) and abs(v) > 3.0 * abs(d0):
            return 3.0 * d0
        return float(v)

    d[0] = end(h[0], h[1], delta[0], delta[1])
    d[-1] = end(h[-1], h[-2], delta[-1], delta[-2])
    return d


def pchip_eval(x: np.ndarray, y: np.ndarray, d: np.ndarray, xq: np.ndarray) -> np.ndarray:
    k = np.clip(np.searchsorted(x, xq, side="right") - 1, 0, x.size - 2)
    h = x[k + 1] - x[k]
    t = (xq - x[k]) / h
    h00 = (1 + 2 * t) * (1 - t) ** 2
    h10 = t * (1 - t) ** 2
    h01 = t**2 * (3 - 2 * t)
    h11 = t**2 * (t - 1)
    return h00 * y[k] + h10 * h * d[k] + h01 * y[k + 1] + h11 * h * d[k + 1]


def cc_int(
    data: Icru90Water, proj: Projectile, window: tuple[float, float]
) -> list[dict[str, Any]]:
    """Log-log interpolation against a monotone cubic in log-log space (both from the array only),
    sampled uniformly in ln E (``SAMPLES_PER_DECADE`` per decade), max and RMS per decade."""
    e_u, s = arrays_per_u(data, proj)
    x, y = np.log(e_u), np.log(s)
    d = pchip_slopes(x, y)
    out = []
    for lo_d, hi_d in DECADES:
        lo, hi = max(lo_d, window[0]), min(hi_d, window[1])
        if not lo < hi:
            continue
        n = max(2, int(math.ceil(SAMPLES_PER_DECADE * math.log10(hi / lo))) + 1)
        xq = np.linspace(math.log(lo), math.log(hi), n)
        dev = np.exp(np.interp(xq, x, y) - pchip_eval(x, y, d, xq)) - 1.0
        out.append({"decade_per_u_mev": [lo_d, hi_d], "window_per_u_mev": [lo, hi], **_stats(dev)})
    return out


# ---- CC-I2 ----


def load_compare_nist() -> Any:
    """The module ``compare_nist.py`` (its ``_icru_range_table`` is the I2 construction)."""
    path = HERE / "compare_nist.py"
    spec = importlib.util.spec_from_file_location("v13_compare_nist", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["v13_compare_nist"] = mod
    spec.loader.exec_module(mod)
    return mod


def i2_table(data: Icru90Water, star: StarTable, proj: Projectile) -> Any:
    """The I2 "integrated ICRU 90" range table (log-log fine grid started at the NIST range)."""
    mod = load_compare_nist()
    if proj.name == "proton":
        return mod._icru_range_table(proj, star, data.proton_energy_mev, data.proton_stopping)
    return mod._icru_range_table(proj, star, data.alpha_energy_mev, data.alpha_stopping)


def cc_i2_arrays(
    model: Icru90WaterSource, data: Icru90Water, star: StarTable, proj: Projectile
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``R(E_i) - R(E_0)`` of the model and of the I2 construction on the common fine grid, and
    the grid itself (the node energies are identical)."""
    t_model = model.table(water(78.0), proj)
    t_i2 = i2_table(data, star, proj)
    if not np.array_equal(t_model.energy_per_u, t_i2.energy_per_u):
        raise AssertionError("the fine grids of the model and of the I2 construction differ")
    r_m = t_model.csda_range_g_cm2 - t_model.csda_range_g_cm2[0]
    r_i = t_i2.csda_range_g_cm2 - t_i2.csda_range_g_cm2[0]
    return r_m, r_i, t_model.energy_per_u


def cc_i2_max_rel(
    model: Icru90WaterSource, data: Icru90Water, star: StarTable, proj: Projectile
) -> float:
    """Maximum relative deviation of the node arrays at the nodes after the first (gated 1e-12)."""
    r_m, r_i, _ = cc_i2_arrays(model, data, star, proj)
    return float(np.max(np.abs(r_m[1:] / r_i[1:] - 1.0)))


def cc_i2_start_offset(data: Icru90Water, star: StarTable, proj: Projectile) -> float:
    """Relative offset of the model's start range (own-array integral to 1 MeV/u) against the NIST
    start of the I2 construction: ``R_own / R_NIST - 1`` at 1 MeV/u."""
    e_u, s = arrays_per_u(data, proj)
    own = icru90_start_range(e_u, s, float(proj.a))
    nist = float(
        np.exp(np.interp(np.log(proj.a), np.log(star.energy_mev), np.log(star.csda_range)))
    )
    return own / nist - 1.0


# ---- V-C3 ----


def v_c3(model: Icru90WaterSource, star: StarTable, proj: Projectile) -> dict[str, Any]:
    """Model against the NIST program over the common window: stopping power and CSDA range at the
    NIST nodes inside the model window (max and RMS relative deviation)."""
    t = model.table(water(78.0), proj)
    e_u = star.energy_mev / proj.a
    lo, hi = float(t.energy_per_u[0]), float(t.energy_per_u[-1])
    m = (e_u >= lo * (1 + 1e-12)) & (e_u <= hi * (1 + 1e-12))
    e_c = np.clip(e_u[m], lo, hi)
    return {
        "reference": star.program,
        "reference_dataset_id": star.dataset_id,
        "window_per_u_mev": [lo, hi],  # the model window, not the reference grid energies
        "stopping": _stats(t.stopping_at(e_c) / star.s_electronic[m] - 1.0),
        "csda_range": _stats(t.range_at(e_c) / star.csda_range[m] - 1.0),
    }


# ---- aggregate documents ----


def _provenance(cache_dir: Path, worktree: Path) -> dict[str, Any]:
    objs = {}
    for ds_id in (ICRU_ID, PSTAR_ID, ASTAR_ID):
        m = cache.read_manifest(ds_id, cache_dir)
        if m is None:
            raise FileNotFoundError(f"dataset {ds_id!r} is not cached in {cache_dir}")
        objs[f"{ds_id} (cache object)"] = str(m["sha256"])
    files = {rel: sha256(worktree / rel) for rel in INPUT_FILES}
    return {
        "input_sha256": {**objs, **files},
        "script_sha256": sha256(Path(__file__).resolve()),
        "note": NOTE,
    }


def compute(cache_dir: Path, worktree: Path = WORKTREE) -> dict[str, dict[str, Any]]:
    """The two documents ``cc-int`` and ``v-c3`` (aggregates only)."""
    model, data, pstar, astar = load_inputs(cache_dir)
    star = {"proton": pstar, "alpha": astar}
    prov = _provenance(cache_dir, worktree)
    cc: dict[str, Any] = {"model_id": "icru90-water-v1", "rows": "CC-INT, CC-NODE, CC-I2, CC-SPAN"}
    cc["points_per_decade_fine_grid"] = 200
    cc["cc_int_samples_per_decade"] = SAMPLES_PER_DECADE
    per: dict[str, Any] = {}
    for proj in PROJECTILES:
        e_u, _ = arrays_per_u(data, proj)
        win = model.window_per_u(proj)
        per[proj.name] = {
            "cc_span": {
                "n_nodes": int(e_u.size),
                "window_per_u_mev": list(win),
                "window_inside_span": bool(e_u[0] <= win[0] and win[1] <= e_u[-1] * (1 + 1e-12)),
                "n_nodes_below_window": int(np.count_nonzero(e_u < win[0])),
                "n_nodes_in_window": int(np.count_nonzero((e_u >= win[0]) & (e_u <= win[1]))),
            },
            "cc_int_per_decade": cc_int(data, proj, win),
            "cc_node_chord_error_at_icru_nodes": cc_node_chord(model, data, proj),
            "cc_node_max_rel_dev_at_fine_nodes": cc_node_fine(model, data, proj),
            "cc_i2_max_rel_dev_node_arrays": cc_i2_max_rel(model, data, star[proj.name], proj),
            "cc_i2_start_range_relative_offset_vs_nist_start": cc_i2_start_offset(
                data, star[proj.name], proj
            ),
        }
    cc["species"] = per
    vc3: dict[str, Any] = {"model_id": "icru90-water-v1", "row": "V-C3"}
    vc3["species"] = {p.name: v_c3(model, star[p.name], p) for p in PROJECTILES}
    return {"cc-int": {**cc, **prov}, "v-c3": {**vc3, **prov}}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", type=Path, default=cache.resolve_cache_dir())
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    docs = compute(args.cache_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for name, doc in docs.items():
        (args.out_dir / f"{name}.json").write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
