"""V3-013 seal material: the ionmc sample definitions and the model-only predictions.

``validation/scripts/stopping/steps_v13.py`` (plan ``validation/plans/v3-013-acceptance.md``
sections A.6, C, D steps 3-4). Pure checks (seeds, grid edges against the TOPAS bundle, the
precondition arithmetic with the plan literals, the configurations, the fresh-cache hook) run
everywhere; the recomputation of ``vc-predictions.json`` needs the ICRU 90 cache object and skips
with a record when it is absent (``IONMC_REQUIRE_DATA=1`` makes that a failure). No sample is run
here.
"""

from __future__ import annotations

import importlib.util
import json
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PLAN = REPO / "validation" / "plans" / "v3-013-acceptance.md"
TOPAS = REPO / "validation" / "reference_cases" / "topas"
PREDICTIONS = REPO / "validation" / "results" / "stopping" / "v3-013" / "vc-predictions.json"


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


S13 = _load("_steps_v13", REPO / "validation" / "scripts" / "stopping" / "steps_v13.py")
GEN = _load(
    "_gen_v13_for_steps", REPO / "validation" / "scripts" / "stopping" / "gen_v13_bundles.py"
)
FAM = _load("_v3_013_seeds_for_steps", HERE / "test_v3_013_seed_families.py")


def test_seed_formulas_and_family_membership():
    fam = FAM.ion13_seeds()
    assert S13.ION13_BASE == FAM.ION13_BASE == 20612004
    assert [S13.seed_of(c) for c in (1, 2, 3)] == [20612004, 20612104, 20612204]
    assert [S13.r_index_of(c) for c in (1, 2, 3)] == [0, 1, 2]
    for c in (1, 2, 3):
        for k in (0, 1, 99):
            assert S13.seed_of(c, k) == 20612004 + 100 * (c - 1) + k and S13.seed_of(c, k) in fam
            assert S13.d_seed_of(c, k) == 20612004 + 1500 + k and S13.d_seed_of(c, k) in fam
    assert S13.R_INDEX_VC2D == 15
    # r_index 0-2 only (the V-C1/V-C2 samples); 10 (A.2) and 11-14 (LV) are not touched
    used = {S13.seed_of(c, k) for c in (1, 2, 3) for k in range(100)}
    for r in (10, 11, 12, 13, 14, 15):
        assert not used & set(range(20612004 + 100 * r, 20612004 + 100 * r + 100))
    with pytest.raises(ValueError):
        S13.seed_of(1, 100)
    with pytest.raises(ValueError):
        S13.seed_of(4)
    # the plan text carries the formula and the r_index assignment
    c = PLAN.read_text().split("## C. Seeds", 1)[1].split("\n## D.", 1)[0]
    assert "seed = 20612004 + 100·r_index + k" in c
    assert "r_index 0–2 = He 75/125/200 V-C1/V-C2 samples" in c and "r_index 15 = V-C2-D" in c


def test_run_definitions_match_the_plan_and_the_bundle_generator():
    assert (S13.N_HISTORIES, S13.N_BATCHES, S13.N_HISTORIES_D) == (1_000_000, 40, 2000)
    a6 = PLAN.read_text().split("### A.6", 1)[1].split("## B.", 1)[0]
    assert "1e6 histories, 40 batches" in a6 and "2 000 histories per model" in a6
    for c, (e_u, length, hz, r_plan) in S13.CASES.items():
        g = next(x for x in GEN.CASES if x[0] == c)
        assert (g[1], g[3], g[4], g[2]) == (e_u, length, hz, r_plan)
        assert S13.A_NUCLEONS == GEN.A_NUCLEONS and S13.FIELD_MM == GEN.FIELD_MM


def test_idd_grid_edges_equal_the_topas_bundle():
    for c, (e_u, length, hz, _r) in S13.CASES.items():
        inp = (TOPAS / f"he-water-{e_u:g}mevu-v13-idd-em-part0" / "input.txt").read_text()
        nz = next(
            int(ln.split("=")[1]) for ln in inp.splitlines() if ln.startswith("i:Sc/IDD/ZBins")
        )
        hlz = next(
            float(ln.split("=")[1].split()[0])
            for ln in inp.splitlines()
            if ln.startswith("d:Ge/Water/HLZ")
        )
        assert (2 * hlz, nz) == (length, round(length / hz))
        origin, spacing, shape = S13.grid_spec(c)
        assert (
            shape == (1, 1, nz) and origin == (-60.0, -60.0, 0.0) and spacing == (120.0, 120.0, hz)
        )
        edges = S13.idd_edges_z_mm(c)
        assert (
            len(edges) == nz + 1
            and edges[0] == 0.0
            and math.isclose(edges[-1], length, rel_tol=1e-12)
        )
        from ionmc.scoring import ScoringGrid

        grid = ScoringGrid(origin, spacing, shape, name="idd")
        for ax in range(3):
            assert np.array_equal(
                grid.edges_mm(ax), S13.idd_edges_z_mm(c) if ax == 2 else grid.edges_mm(ax)
            )
        assert [nz for nz in (S13.grid_spec(i)[2][2] for i in (1, 2, 3))] == [2000, 4000, 3500]


def test_precondition_arithmetic_with_the_plan_literals():
    a6 = PLAN.read_text().split("### A.6", 1)[1].split("## B.", 1)[0]
    assert "σ_exp = √2 · SE_1e6" in a6 and "SE_1e6 = 0.01 mm · √(1e5/1e6)" in a6
    assert "|ΔR_pred| ≥ 2 · (0.02 mm + 1.645 σ_exp)" in a6
    assert "|ΔR_pred,D| ≥ 5e-4 · R_CSDA,default(E0)" in a6
    se = 0.01 * math.sqrt(1e5 / 1e6)
    sigma_exp = math.sqrt(2.0) * se
    assert math.isclose(S13.SIGMA_EXP_MM, sigma_exp, rel_tol=1e-15) and math.isclose(
        sigma_exp, 0.0044721359549995794, rel_tol=1e-12
    )
    thr = 2.0 * (0.02 + 1.645 * sigma_exp)
    assert S13.VC2_PRECONDITION_MM == thr and math.isclose(thr, 0.0547140, abs_tol=1e-6)
    assert (
        not S13.vc2_informative(0.0547)
        and S13.vc2_informative(0.0548)
        and S13.vc2_informative(-0.0548)
    )
    assert S13.vc2d_informative(0.05, 100.0) and not S13.vc2d_informative(0.0499, 100.0)
    assert S13.vc2_tolerance_mm(0.1) == 0.02 and S13.vc2_tolerance_mm(-0.5) == 0.05
    assert S13.TOL_FLOOR_MM == 0.02 and S13.TOL_FRACTION == 0.10 and S13.Z_TOST == 1.645


def test_gaussian_average():
    f = lambda e: 3.0 * e + 0.5 * e * e  # noqa: E731
    assert S13.gaussian_average(f, 10.0, 0.0) == f(10.0)
    assert math.isclose(S13.gaussian_average(f, 10.0, 2.0), f(10.0) + 0.5 * 4.0, rel_tol=1e-12)


def test_configurations_of_the_samples():
    for c, (e_u, length, hz, _r) in S13.CASES.items():
        cfg = S13.build_config(
            c,
            "bethe-ionmc-v1",
            n_histories=S13.N_HISTORIES,
            n_batches=S13.N_BATCHES,
            backend="warp-cpu",
            seed=S13.seed_of(c),
        )
        ph, run = cfg.physics, cfg.run
        assert (run.backend, run.precision, run.n_histories, run.n_batches, run.seed) == (
            "warp-cpu",
            "float64",
            1_000_000,
            40,
            S13.seed_of(c),
        )
        assert ph.stopping_model == "bethe-ionmc-v1" and ph.stopping is None and ph.nuclear is False
        assert ph.multiple_scattering and ph.straggling and ph.end_of_range_model == "default"
        assert cfg.source.projectile.name == "alpha" and cfg.source.kinetic_energy_mev == 4 * e_u
        assert cfg.geometry.material.mean_excitation_eV == 78.0
        (g,) = cfg.scoring
        assert g.name == "idd" and g.shape == (1, 1, round(length / hz))
        det = S13.build_config(
            c,
            "bethe-ionmc-v1",
            n_histories=2000,
            n_batches=20,
            backend="python",
            seed=S13.d_seed_of(c),
            deterministic=True,
            track_end=True,
        )
        assert not det.physics.multiple_scattering and not det.physics.straggling
        assert det.diagnostics.track_end_positions and det.run.seed == 20612004 + 1500
    with pytest.raises(ValueError):
        S13.build_config(
            1,
            "nist-pstar-astar-water-2005",
            n_histories=1000,
            n_batches=20,
            backend="python",
            seed=1,
        )


def test_pair_seeds_are_independent_of_the_model():
    import inspect

    assert "model" not in inspect.signature(S13.seed_of).parameters
    assert "model" not in inspect.signature(S13.d_seed_of).parameters


def test_fresh_warp_cache_hook(tmp_path, monkeypatch):
    monkeypatch.setenv("WARP_CACHE_PATH", "unchanged")
    busy = tmp_path / "warp-cache-t"
    busy.mkdir()
    (busy / "stale.so").write_bytes(b"x")
    with pytest.raises(SystemExit, match="not empty"):
        S13.prepare_warp_cache(tmp_path, "t")
    assert os.environ["WARP_CACHE_PATH"] == "unchanged"
    if "warp" in sys.modules:
        with pytest.raises(SystemExit, match="imported before"):
            S13.prepare_warp_cache(tmp_path, "fresh")
    monkeypatch.delitem(sys.modules, "warp", raising=False)
    d = S13.prepare_warp_cache(tmp_path, "fresh")
    assert d.is_dir() and not any(d.iterdir()) and os.environ["WARP_CACHE_PATH"] == str(d.resolve())
    assert S13.warp_modules(d) == {}
    (d / "wp_mod").mkdir()
    (d / "wp_mod" / "a.cpp").write_bytes(b"abc")
    (d / "b.so").write_bytes(b"")
    mods = S13.warp_modules(d)
    assert set(mods) == {"wp_mod/a.cpp", "b.so"}
    assert (
        mods["wp_mod/a.cpp"] == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    )


def test_environment_digest_hook_degrades_to_a_recorded_note():
    digests, note = S13.environment_digests()
    if (REPO / "validation" / "scripts" / "stopping" / "env_digest.py").is_file():
        assert (digests is None) != (note is None)
    else:
        assert digests is None and note == S13.HAVE_ENV_DIGEST_NOTE


# -- the committed prediction file --
def _doc() -> dict[str, Any]:
    return json.loads(PREDICTIONS.read_text())


def test_committed_predictions_are_consistent():
    d = _doc()
    assert d["schema"] == S13.PREDICTIONS_SCHEMA
    assert PREDICTIONS.read_bytes().endswith(b"\n") and not PREDICTIONS.read_bytes().endswith(
        b"\n\n"
    )
    assert d["constants"]["vc2_precondition_mm"] == S13.VC2_PRECONDITION_MM
    assert d["constants"]["sigma_exp_mm"] == S13.SIGMA_EXP_MM
    assert [e["case"] for e in d["energies"]] == [1, 2, 3]
    for e in d["energies"]:
        assert e["sigma_e_mev"] == 0.0
        assert math.isclose(e["delta_r_pred_mm"], e["delta_r_pred_d_mm"], rel_tol=0, abs_tol=1e-12)
        assert math.isclose(
            e["delta_r_pred_d_mm"],
            e["r_csda_icru90_mm"] - e["r_csda_default_mm"],
            rel_tol=0,
            abs_tol=1e-12,
        )
        assert (
            e["r_plan_mm"] == S13.CASES[e["case"]][3] and e["e0_total_mev"] == 4 * e["e_mev_per_u"]
        )
        assert e["vc2_precondition_holds"] == S13.vc2_informative(e["delta_r_pred_mm"])
        assert e["vc2d_informative"] == S13.vc2d_informative(
            e["delta_r_pred_d_mm"], e["r_csda_default_mm"]
        )
        assert e["vc2_gate"].startswith(
            "gating" if e["vc2_precondition_holds"] else "non-informative"
        )
        assert e["vc2d_gate"].startswith("gating" if e["vc2d_informative"] else "substitute")
        assert math.isclose(
            e["vc2_tolerance_mm"], S13.vc2_tolerance_mm(e["delta_r_pred_mm"]), rel_tol=1e-12
        )
    assert set(d["models"]) == set(S13.MODELS)
    for m, pin in d["models"].items():
        assert pin["model_id_effective"] == m and len(pin["spec_sha256"]) == 64
        assert pin["content_sha256"] is None or len(pin["content_sha256"]) == 64
    assert d["models"]["icru90-water-v1"]["content_sha256"] is not None

    # aggregates only: no array of values anywhere in the file
    def no_arrays(x: Any) -> None:
        if isinstance(x, dict):
            for v in x.values():
                no_arrays(v)
        elif isinstance(x, list):
            assert len(x) <= 8 and all(not isinstance(v, float) for v in x)
            for v in x:
                no_arrays(v)

    no_arrays(d)


V13_PREDICTIONS_DERIVATION_COMMIT = "53d779ae73d1520080b3350b0b979c32795604af"
V13_PREDICTIONS_TREE = "c5b33495898cbf45f78c72a3b7e416610d369e83"


def test_committed_predictions_pin_script_and_tree():
    """The committed predictions pin the script and the ``src/ionmc`` tree at their derivation.

    After the V3-013 integration (develop 95937553) the tree pin is a historical record of the
    derivation commit (V3-013 plan, post-evaluation closure of 2026-11-23): later tasks change
    ``src/ionmc`` without re-deriving the sealed predictions, so the pin is compared with the tree
    of the recorded derivation commit, not with HEAD.
    """
    d = _doc()
    assert d["script"]["path"] == S13.SELF_REL
    assert d["script"]["sha256"] == S13._sha256_file(REPO / S13.SELF_REL)
    assert d["src_ionmc_tree"] == V13_PREDICTIONS_TREE
    try:
        tree = subprocess.run(
            ["git", "-C", str(REPO), "rev-parse", f"{V13_PREDICTIONS_DERIVATION_COMMIT}:src/ionmc"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return  # the derivation commit is not reachable (shallow clone): the constant above pins it
    assert tree == V13_PREDICTIONS_TREE


def _close(a: Any, b: Any, path: str = "") -> None:
    if isinstance(a, dict):
        assert isinstance(b, dict) and set(a) == set(b), path
        for k in a:
            _close(a[k], b[k], f"{path}/{k}")
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            _close(x, y, f"{path}[{i}]")
    elif isinstance(a, float) and not isinstance(a, bool):
        assert isinstance(b, float) and math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12), (
            path,
            a,
            b,
        )
    else:
        assert a == b, (path, a, b)


def test_predictions_reproduce_from_the_committed_inputs():
    from ionmc.data import cache

    cdir = cache.resolve_cache_dir()
    try:
        cache.verify("geant4-icru90-stopping-11.4.2", cdir)
    except (FileNotFoundError, cache.IntegrityError) as exc:
        if os.environ.get("IONMC_REQUIRE_DATA") == "1":
            pytest.fail(f"IONMC_REQUIRE_DATA=1 but the ICRU 90 object is not cached: {exc}")
        pytest.skip(
            f"ICRU 90 object not cached in {cdir}: the committed predictions are not recomputed"
        )
    _close(_doc(), S13.compute_predictions())


def test_normalise_cache_dir_keeps_null_and_replaces_only_strings() -> None:
    """Codex REVIEW-f5d6f67e finding 1: an analytic model records a null cache directory and
    the sealed pin must keep it null; a data-backed model's directory string (top level and in
    every data_used entry) becomes the placeholder; nothing else changes."""
    analytic = {"model_id_effective": "bethe-ionmc-v1", "cache_dir": None, "data_used": []}
    out = S13.normalise_cache_dir(analytic)
    assert out["cache_dir"] is None and out["data_used"] == []
    assert out == analytic and out is not analytic
    backed = {
        "model_id_effective": "icru90-water-v1",
        "cache_dir": "/abs/cache/ionmc-data",
        "data_used": [
            {"key": "icru90", "cache_dir": "/abs/cache/ionmc-data"},
            {"key": "x", "cache_dir": None},
        ],
    }
    out = S13.normalise_cache_dir(backed)
    assert out["cache_dir"] == S13.CACHE_DIR_PLACEHOLDER
    assert out["data_used"] == [
        {"key": "icru90", "cache_dir": S13.CACHE_DIR_PLACEHOLDER},
        {"key": "x", "cache_dir": None},
    ]
    assert backed["cache_dir"] == "/abs/cache/ionmc-data"  # the input is not modified
    assert "cache_dir" not in S13.normalise_cache_dir({"model_id_effective": "m"})
