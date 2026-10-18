"""D-E3 variant elastic tables (exploratory, report-only; Amendment 19 (b); plan section 3.1).

Loads the pinned elastic table ``1fd24cff...`` and writes the eight D-E3 variants into a scratch cache copy
(``IONMC_CACHE_DIR`` for the runs). Nothing under ``src/``, ``tests/ionmc``, ``validation/scripts/transport`` or
``benchmarks/transport`` is touched; the pinned cache is only read.

Variants (sigma unchanged except where stated; 256 equiprobable mu_CM bins per node as in the pinned table):

- ``all``          byte-identical arrays of the pinned table (base copy; the bit-for-bit check of ``--check-base``);
- ``pp-only``      sigma_el = 0 for every target other than H-1;
- ``pA-only``      sigma_pp = 0 (H-1);
- ``pA-sigma0.8``  sigma_el x 0.8 for every target other than H-1;
- ``pA-Rx1.2``     p+A shape = black disk at 1.2 R(E);
- ``pA-s0.5`` / ``pA-s0.9``  black disk times exp(-q^2 s^2), s in fm;
- ``pA-gauss``     exp(-q^2 R^2 / 4).

The parametric shapes are the objects of D-E0 (``de1_shapes.Shapes.mu_quantiles``, 16385-point shape grid as the
builder). Every variant carries an ``exploratory`` key, the recomputed npz sha256 and the id of
``ionmc.nuclear.build.table_id``; the qualification flags are copied from the pinned table unearned (disclosed in the
key). It is loaded only through ``PhysicsConfig(elastic_table_id=..., allow_exploratory_table=True)``.

Scratch cache: ``--scratch-cache DIR``, else the environment variable ``IONMC_DE3_SCRATCH``, else the portable default
``validation/generated/elastic_attribution/de3-cache/`` (git-ignored). The directory holds a marker file
``.de3-scratch-cache``; the script refuses to rebuild or clear a directory that is not empty and has no marker, and a
directory that is, contains or lies inside the pinned data cache. ``IONMC_CACHE_DIR`` is set only for the scope of the preflight
(printed) and then restored; an already-set value is never overwritten silently.

Run: ``python -m validation.scripts.elastic_attribution.make_variant_tables [--scratch-cache DIR]`` (writes the cache, prints ids,
normalisation statistics and the preflight).
"""

from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np

from . import common as c
from .de1_shapes import HBARC, family_weight

DEFAULT_SCRATCH = c.OUT / "de3-cache"  # inside the repository, git-ignored (validation/generated/)
MARKER = ".de3-scratch-cache"
SCRATCH_ENV = "IONMC_DE3_SCRATCH"


def resolve_scratch(arg: str | Path | None = None) -> Path:
    """The scratch cache directory: ``--scratch-cache`` argument, else the environment variable ``IONMC_DE3_SCRATCH``, else
    ``validation/generated/elastic_attribution/de3-cache/``."""
    return Path(arg or os.environ.get(SCRATCH_ENV) or DEFAULT_SCRATCH).expanduser()


SCRATCH = resolve_scratch()
REAL_CACHE = c.CACHE_DIR  # the configured pinned-data cache (IONMC_CACHE_DIR at import, else the repository cache); read only
NUCLEAR_TABLE_ID = "00e8031d5f67704f256e904ebab89b0c2a75845bcf1887c29e1518bffaafe5fb"
SHAPE_KIND = {"pA-Rx1.2": "Rx1.2", "pA-s0.5": "s0.5", "pA-s0.9": "s0.9", "pA-gauss": "gauss"}
VARIANTS = ("all", "pp-only", "pA-only", "pA-sigma0.8", *SHAPE_KIND)
PURPOSE = (
    "V3-005D D-E3 sensitivity probe (Amendment 19 (b)), report-only, never evidence; qualification flags copied from the "
    "pinned table unearned"
)
NORM_CHECK_STEP = 8  # every 8th valid node is re-integrated on a 4x finer shape grid


def _derived(root: Path) -> Path:
    return root / "derived"


def _build_arrays(
    variant: str, base: dict[str, np.ndarray], names: tuple[str, ...]
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Variant arrays and the normalisation statistics (max |bin probability - 1/256| on the fine grid)."""
    a = {k: v.copy() for k, v in base.items()}
    stats: dict[str, Any] = {}
    pa = [i for i, n in enumerate(names) if n != "H-1"]
    h = names.index("H-1")
    if variant == "pp-only":
        a["sigma_barn"][pa] = 0.0
    elif variant == "pA-only":
        a["sigma_barn"][h] = 0.0
    elif variant == "pA-sigma0.8":
        a["sigma_barn"][pa] *= 0.8
    elif variant in SHAPE_KIND:
        kind = SHAPE_KIND[variant]
        sh = _ShapeBuilder(base, names)
        n_edges = a["edges_mu"].shape[2]
        u = np.arange(n_edges) / (n_edges - 1)
        worst_fine = 0.0
        n_nodes = 0
        for it in pa:
            for k in range(a["edges_mu"].shape[1]):
                if not a["valid"][it, k]:
                    continue  # nodes below the declared domain keep the pinned (unused) edges
                edges = sh.edges(kind, it, k, u)
                a["edges_mu"][it, k] = edges
                n_nodes += 1
                if k % NORM_CHECK_STEP == 0:
                    worst_fine = max(worst_fine, sh.fine_dev(kind, it, k, edges))
        stats = {"nodes_rebuilt": n_nodes, "max_bin_probability_deviation_fine_grid": worst_fine,
                 "n_quantiles": n_edges - 1}  # fmt: skip
    return a, stats


class _ShapeBuilder:
    """Quantile edges of one family at one table node through ``de1_shapes`` (radius and p_CM of the node)."""

    def __init__(self, base: dict[str, np.ndarray], names: tuple[str, ...]) -> None:
        self.base, self.names = base, names

    def _wcdf(self, kind: str, it: int, k: int, n_s: int) -> tuple[np.ndarray, np.ndarray]:
        s = np.linspace(0.0, 1.0, n_s)
        r, pcm = float(self.base["radius_fm"][it, k]), float(self.base["p_cm_mev"][it, k])
        g = 4.0 * s * family_weight(kind, 2.0 * pcm / HBARC * s, r)
        f = np.concatenate(([0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(s))))
        return s, f / f[-1]

    def edges(self, kind: str, it: int, k: int, u: np.ndarray) -> np.ndarray:
        s, f = self._wcdf(kind, it, k, 16385)
        s_u = np.interp(1.0 - u, f, s)
        e = 1.0 - 2.0 * s_u * s_u
        e[0], e[-1] = -1.0, 1.0
        if not np.all(np.diff(e) > 0.0):
            raise AssertionError(f"edges not strictly ascending: {self.names[it]} node {k}")
        return e

    def fine_dev(self, kind: str, it: int, k: int, edges: np.ndarray) -> float:
        """Max |integral of the density over a bin - 1/n_q| on a 4x finer shape grid."""
        s, f = self._wcdf(kind, it, k, 4 * 16384 + 1)
        s_e = np.sqrt(np.clip((1.0 - edges) / 2.0, 0.0, 1.0))
        frac_forward = np.interp(s_e, s, f)  # fraction with s' < s(edge); mu descends with s
        p = -np.diff(frac_forward)  # events between consecutive ascending-mu edges
        return float(np.max(np.abs(p - 1.0 / (edges.size - 1))))


def _write_variant(
    variant: str, pinned_json: dict[str, Any], base: dict[str, np.ndarray], out_dir: Path
) -> dict[str, Any]:
    from ionmc.nuclear.build import table_id

    names = tuple(pinned_json["target_names"])
    arrays, stats = _build_arrays(variant, base, names)
    _assert_table(arrays, names)
    tmp = out_dir / "tmp-variant.npz"
    np.savez(tmp, **arrays)
    sha = c.sha256_file(tmp)
    side = json.loads(json.dumps(pinned_json))
    side["npz_sha256"] = sha
    side["exploratory"] = {"d_e3_variant": variant, "source_table_id": c.TABLE_ID, "purpose": PURPOSE,
                           "parameters": _parameters(variant), "normalisation": stats}  # fmt: skip
    side["table_id"] = table_id(side)
    tid = side["table_id"]
    tmp.rename(out_dir / f"elastic-proton-{tid}.npz")
    (out_dir / f"elastic-proton-{tid}.json").write_text(
        json.dumps(side, indent=1, sort_keys=True) + "\n"
    )
    return {"variant": variant, "table_id": tid, "npz_sha256": sha, "normalisation": stats}


def _parameters(variant: str) -> dict[str, Any]:
    return {
        "all": {},
        "pp-only": {"sigma_el_pA": 0.0},
        "pA-only": {"sigma_pp": 0.0},
        "pA-sigma0.8": {"sigma_scale_pA": 0.8},
        "pA-Rx1.2": {"R_factor": 1.2},
        "pA-s0.5": {"s_fm": 0.5},
        "pA-s0.9": {"s_fm": 0.9},
        "pA-gauss": {"form": "exp(-q^2 R^2/4)"},
    }[variant]


def _assert_table(
    a: dict[str, np.ndarray], names: tuple[str, ...], check_ends: bool = True
) -> None:
    """Normalisation (n_q equiprobable bins, endpoints, ascending), non-negative density, sigma = 0 below the domain."""
    edges = a["edges_mu"]
    for it, name in enumerate(names):
        lo, hi = (-0.96, 0.96) if name == "H-1" else (-1.0, 1.0)
        v = a["valid"][it].astype(bool)
        e = edges[it][v]
        if check_ends:
            assert np.allclose(e[:, 0], lo) and np.allclose(e[:, -1], hi), f"{name}: edge range"
        else:  # swapped nodes: the end edges are the sample extremes, inside the physical range
            assert e[:, 0].min() >= lo - 1e-12 and e[:, -1].max() <= hi + 1e-12, (
                f"{name}: edges outside the range"
            )
        assert np.all(np.diff(e, axis=1) > 0.0), (
            f"{name}: density not positive (edges not ascending)"
        )
        below = a["grid_e_mev"] < a["target_e_min_mev"][it]
        assert np.all(a["sigma_barn"][it][below] == 0.0) and not np.any(a["valid"][it][below]), (
            f"{name}: sigma != 0 below domain"
        )
        assert np.all(a["sigma_barn"][it] >= 0.0), f"{name}: negative sigma"


def check_scratch(scratch: Path) -> None:
    """Refuse a scratch directory that could be (or hold) the real cache, or that is non-empty without the marker."""
    sp, rp = scratch.resolve(), REAL_CACHE.resolve()
    if sp == rp or rp in sp.parents or sp in rp.parents:
        raise SystemExit(
            f"refusing to use {scratch}: it is, contains or lies inside the pinned data cache {REAL_CACHE}"
        )
    if scratch.exists() and any(scratch.iterdir()) and not (scratch / MARKER).is_file():
        raise SystemExit(
            f"refusing to rebuild {scratch}: it is not empty and has no {MARKER} marker (not a D-E3 scratch cache)"
        )


def build(scratch: Path = SCRATCH) -> dict[str, Any]:
    """(Re)build the scratch cache: manifests and objects (copied), the pinned nuclear and elastic tables, the variants."""
    src = REAL_CACHE
    check_scratch(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    (scratch / MARKER).write_text(
        "D-E3 scratch cache (variant elastic tables); safe to delete and rebuild\n"
    )
    for sub in ("manifests", "objects"):
        if not (scratch / sub).is_dir():
            shutil.copytree(src / sub, scratch / sub)
    der = _derived(scratch)
    kept: dict[str, Any] = {}
    if (scratch / "variants.json").is_file() and der.is_dir():  # installed builder tables survive a rebuild
        for v, e in _builder_entries(scratch).items():
            kept[v] = (e, {n: (der / n).read_bytes() for n in (f"elastic-proton-{e['table_id']}.json", f"elastic-proton-{e['table_id']}.npz")})
    if der.is_dir():
        shutil.rmtree(der)
    der.mkdir()
    for _e, files in kept.values():
        for n, b in files.items():
            (der / n).write_bytes(b)
    for stem in (f"nuclear-proton-{NUCLEAR_TABLE_ID}", f"elastic-proton-{c.TABLE_ID}"):
        for ext in (".json", ".npz"):
            shutil.copy2(src / "derived" / (stem + ext), der / (stem + ext))
    pinned = json.loads((src / "derived" / f"elastic-proton-{c.TABLE_ID}.json").read_text())
    base = {
        k: v
        for k, v in np.load(
            src / "derived" / f"elastic-proton-{c.TABLE_ID}.npz", allow_pickle=False
        ).items()
    }
    out = {v: _write_variant(v, pinned, base, der) for v in VARIANTS}
    doc = {
        "banner": c.BANNER,
        "pinned_table_id": c.TABLE_ID,
        "nuclear_table_id": NUCLEAR_TABLE_ID,
        "variants": out,
        **({"builders": {v: e for v, (e, _f) in kept.items()}} if kept else {}),
    }
    (scratch / "variants.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    return doc


BUILDER_VARIANTS = {
    "m2prime": {"model": "m2prime-pauli-blocked", "pauli_variant": "clementel-villi"},
    "m2": {"model": "m2-optical-glauber", "pauli_variant": None},
    "rb": {"model": "m2prime-pauli-blocked", "pauli_variant": "bdc-pn"},
}
"""Builder tables (``ionmc data build elastic-proton --model m2|m2prime``; V3-005D Amendments 20 (c), 21 (e)) installable as named variants:
``m2prime`` (R-C Clementel-Villi, the model under test), ``m2`` (report-only) and ``rb`` (Bertulani & De Conti reading R-B, report-only)."""


def _builder_entries(scratch: Path) -> dict[str, dict[str, Any]]:
    doc = json.loads((scratch / "variants.json").read_text())
    return dict(doc.get("builders", {}))


def install_builder_table(scratch: Path, variant: str, source: Path) -> dict[str, Any]:
    """Install the builder table pair (``source``: a directory holding exactly one ``elastic-proton-<id>.json/.npz`` pair, or the json) into
    the scratch cache as the named variant. Fail closed: the npz sha256 equals the sidecar's, the table id RECOMPUTED from the sidecar (with
    that sha256) equals the file name, and the exploratory block names the model and Pauli variant expected for ``variant``. The binding
    ``{table_id, npz_sha256, json_sha256}`` is recorded in ``variants.json`` (``builders``); installing the same pair again is a no-op, a
    different table for an installed variant is refused."""
    from ionmc.nuclear.build import table_id

    if variant not in BUILDER_VARIANTS:
        raise SystemExit(f"unknown builder variant {variant!r}; expected one of {sorted(BUILDER_VARIANTS)}")
    check_scratch(scratch)
    if not (scratch / "variants.json").is_file() or not (scratch / MARKER).is_file():
        raise SystemExit(f"{scratch}: not a built D-E3 scratch cache (run make_variant_tables first)")
    src = Path(source)
    jp = src if src.is_file() else _only(src, "elastic-proton-*.json")
    side = json.loads(jp.read_text())
    tid = side["table_id"]
    npz_p = jp.with_suffix(".npz")
    if jp.name != f"elastic-proton-{tid}.json" or not npz_p.is_file():
        raise SystemExit(f"{jp}: file name does not carry the recorded table id or the npz is missing")
    sha = c.sha256_file(npz_p)
    if sha != side.get("npz_sha256") or table_id({**side, "npz_sha256": sha}) != tid:
        raise SystemExit(f"{jp}: npz sha256 or the recomputed table id differs from the record")
    want, ex = BUILDER_VARIANTS[variant], side.get("exploratory") or {}
    if ex.get("model") != want["model"] or ex.get("pauli_variant") != want["pauli_variant"]:
        raise SystemExit(f"{jp}: exploratory block {ex.get('model')!r}/{ex.get('pauli_variant')!r} is not the {variant!r} construction")
    entries = _builder_entries(scratch)
    if variant in entries and entries[variant]["table_id"] != tid:
        raise SystemExit(f"variant {variant!r} is already installed with table {entries[variant]['table_id'][:12]}; use a fresh scratch cache")
    der = _derived(scratch)
    for p in (jp, npz_p):
        shutil.copy2(p, der / p.name)
    entry = {"table_id": tid, "npz_sha256": sha, "json_sha256": c.sha256_file(jp), "model": ex["model"], "pauli_variant": ex.get("pauli_variant")}
    doc = json.loads((scratch / "variants.json").read_text())
    doc.setdefault("builders", {})[variant] = entry
    (scratch / "variants.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    return entry


def _only(d: Path, pattern: str) -> Path:
    files = sorted(Path(d).glob(pattern))
    if len(files) != 1:
        raise SystemExit(f"{d}: expected exactly one {pattern}, found {len(files)}")
    return files[0]


def verify_variant(scratch: Path, variant: str) -> dict[str, str]:
    """``{table_id, npz_sha256}`` of an installed builder variant, re-verified against the files in the scratch cache (the binding that
    the run documents record)."""
    e = _builder_entries(scratch).get(variant)
    if e is None:
        raise SystemExit(f"builder variant {variant!r} is not installed in {scratch}")
    der = _derived(scratch)
    if c.sha256_file(der / f"elastic-proton-{e['table_id']}.npz") != e["npz_sha256"] or c.sha256_file(der / f"elastic-proton-{e['table_id']}.json") != e["json_sha256"]:
        raise SystemExit(f"builder variant {variant!r}: the installed files differ from the recorded binding")
    return {"table_id": e["table_id"], "npz_sha256": e["npz_sha256"]}


def variant_ids(scratch: Path = SCRATCH) -> dict[str, str]:
    doc = json.loads((scratch / "variants.json").read_text())
    return {**{v: r["table_id"] for v, r in doc["variants"].items()}, **{v: r["table_id"] for v, r in doc.get("builders", {}).items()}}


def preflight(scratch: Path = SCRATCH) -> dict[str, Any]:
    """With ``IONMC_CACHE_DIR`` = the scratch cache: the derived directory holds exactly the expected set, every table
    loads (variants with ``allow_exploratory``), and a Simulation of every variant loads exactly the expected ids."""
    with c.scoped_cache_env(scratch):
        return _preflight(scratch)


def _preflight(scratch: Path) -> dict[str, Any]:
    from ionmc.nuclear import elastic_tables as et
    from ionmc.nuclear import tables as nt

    ids = variant_ids(scratch)
    builders = _builder_entries(scratch)
    expected = {
        f"nuclear-proton-{NUCLEAR_TABLE_ID}",
        f"elastic-proton-{c.TABLE_ID}",
        *(f"elastic-proton-{t}" for t in ids.values()),
    }
    present = {p.stem for p in _derived(scratch).glob("*.npz")} | {
        p.stem for p in _derived(scratch).glob("*.json")
    }
    assert present == expected, (
        f"derived set differs: extra {sorted(present - expected)}, missing {sorted(expected - present)}"
    )
    for v, tid in ids.items():
        t = et.ElasticTable.load(scratch, tid, allow_exploratory=True)
        ex = t.info["exploratory"]
        if v in builders:
            verify_variant(scratch, v)
            assert t.exploratory and ex.get("model") == BUILDER_VARIANTS[v]["model"] and ex.get("pauli_variant") == BUILDER_VARIANTS[v]["pauli_variant"]
            continue
        assert t.exploratory and (ex.get("d_e3_variant") or ex.get("d_e6_variant")) == v
    loaded: list[str] = []
    orig_e, orig_n = et.ElasticTable.load.__func__, nt.NuclearTable.load.__func__  # type: ignore[attr-defined]

    def log_e(cls: Any, cache_dir: Any, table_id: str, allow_exploratory: bool = False) -> Any:
        loaded.append(f"elastic-proton-{table_id}")
        return orig_e(cls, cache_dir, table_id, allow_exploratory)

    def log_n(cls: Any, cache_dir: Any, table_id: str) -> Any:
        loaded.append(f"nuclear-proton-{table_id}")
        return orig_n(cls, cache_dir, table_id)

    et.ElasticTable.load = classmethod(log_e)  # type: ignore[method-assign,assignment]
    nt.NuclearTable.load = classmethod(log_n)  # type: ignore[method-assign,assignment]
    res: dict[str, Any] = {}
    try:
        from .de3_probes import build_config

        for v, tid in ids.items():
            loaded.clear()
            from ionmc.simulation import Simulation

            Simulation(build_config(150, "emel", 100, 1, tid))
            want = {f"nuclear-proton-{NUCLEAR_TABLE_ID}", f"elastic-proton-{tid}"}
            res[v] = {"loaded": sorted(set(loaded)), "ok": set(loaded) == want}
            assert set(loaded) == want, (
                f"{v}: loaded {sorted(set(loaded))}, expected {sorted(want)}"
            )
    finally:
        et.ElasticTable.load = classmethod(orig_e)  # type: ignore[method-assign,assignment]
        nt.NuclearTable.load = classmethod(orig_n)  # type: ignore[method-assign,assignment]
    return {"derived_files_expected_equals_present": True, "per_variant": res}


def main(argv: list[str] | None = None) -> None:
    import argparse

    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument(
        "--scratch-cache",
        type=Path,
        default=None,
        help=f"scratch cache directory (env {SCRATCH_ENV}; default {DEFAULT_SCRATCH})",
    )
    ap.add_argument("--install", nargs=2, action="append", metavar=("VARIANT", "SOURCE"), default=[],
                    help=f"install a builder table pair as a named variant ({', '.join(BUILDER_VARIANTS)}) into the existing scratch cache and run the preflight")  # fmt: skip
    a = ap.parse_args(argv)
    scratch = resolve_scratch(a.scratch_cache)
    print(f"scratch cache: {scratch}", flush=True)
    if a.install:
        for v, src in a.install:
            print(v, json.dumps(install_builder_table(scratch, v, Path(src))))
        print(json.dumps(preflight(scratch), indent=1))
        return
    doc = build(scratch)
    for v, r in doc["variants"].items():
        print(f"{v:12s} {r['table_id']}  npz {r['npz_sha256'][:12]}  {r['normalisation'] or ''}")
    print(json.dumps(preflight(scratch), indent=1))
    sys.stdout.flush()


if __name__ == "__main__":
    main()
