"""Rows CAP-1a, CAP-1b, CAP-1c, CAP-2..CAP-6 of plan v3-006 (step S4): the capability contract of
decision 0045 section 4 (``ionmc.contract``).

Declared == enforced over the enumerated dimensions (sweep against ``validate``), the targeted
positive/negative pair of every other ``validate`` rule, the derived declaration, the frozen
requested -> effective reconciliation table, the discriminating dispatch and the honoured check.
Nothing here runs more than a dozen histories (single process, python backend); sweeps only call
``validate``. Cache-backed cases use the ``_require_datasets`` convention of the other V3-006 tests.
"""

from __future__ import annotations

import ast
import copy
import dataclasses
import inspect
import json
import os
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_v3_006_configio as cio  # noqa: E402
from test_v3_006_results_format import _mutate  # noqa: E402

from ionmc import config, contract  # noqa: E402
from ionmc.config import (  # noqa: E402
    DiagnosticsOptions,
    EffectiveConfig,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
    validate,
)
from ionmc.configio import load_toml, to_dict  # noqa: E402
from ionmc.errors import BackendUnavailableError, UnsupportedCombinationError  # noqa: E402
from ionmc.geometry import BoxPhantom, VoxelGeometry  # noqa: E402
from ionmc.materials import PMMA, WATER  # noqa: E402
from ionmc.physics.projectiles import ALPHA, CARBON12, PROTON  # noqa: E402
from ionmc.physics.stopping import BetheStoppingSource  # noqa: E402
from ionmc.results import ResultFormatError, load, save  # noqa: E402
from ionmc.scoring import ScoringGrid, TallyRequest  # noqa: E402
from ionmc.simulation import (  # noqa: E402
    NuclearTransportCounters,
    Simulation,
    capabilities,
)
from ionmc.sources import PencilBeamSource  # noqa: E402

BETHE = "bethe-ionmc-v1"
LS_EFF = "bethe-ionmc-v1+ls-eff"
HERE = Path(__file__).resolve().parent
SLAB_TOML = HERE / "data" / "v3_006" / "voxel_slab_lookup.toml"


def _require_datasets(*dataset_ids: str) -> None:
    """Skip a cache-backed case without the dataset (CI tier); fail when the dataset is absent
    under IONMC_REQUIRE_DATA=1 or IONMC_LV_TIER=1. A corrupt object, malformed manifest or pin
    mismatch (``cache.IntegrityError``) is never a skip."""
    from ionmc.data import cache

    cdir = cache.resolve_cache_dir(None)
    required = "1" in (os.environ.get("IONMC_REQUIRE_DATA"), os.environ.get("IONMC_LV_TIER"))
    for dataset_id in dataset_ids:
        try:
            cache.verify(dataset_id, cdir)
        except FileNotFoundError as exc:
            if required:
                pytest.fail(f"{dataset_id} is required here but not cached: {exc}")
            pytest.skip(f"{dataset_id} is not in the cache (IONMC_CACHE_DIR)")


def _nuclear_table_id_or_none() -> str | None:
    """A built nuclear table of the cache, or None (the nuclear python cases are then reported
    as not evaluated, never as agreeing)."""
    from ionmc.data import cache
    from ionmc.nuclear.tables import NuclearTable, NuclearTableError

    try:
        cdir = cache.resolve_cache_dir(None)
        for p in sorted((cdir / "derived").glob("nuclear-proton-*.json")):
            tid = p.stem.removeprefix("nuclear-proton-")
            try:
                if "transport_path_bound_terms" in NuclearTable.load(None, tid).info:
                    return tid
            except NuclearTableError:
                continue
    except (OSError, FileNotFoundError):
        pass
    return None


class _SyntheticNuclearTable:
    """A cache-free stand-in for ``NuclearTable`` that carries exactly what ``validate`` and the
    effective summary read (energy bounds, path-bound terms, identity fields, per-material rows).
    It is installed with ``NuclearTable.load`` patched, so every nuclear cell of the sweep runs the
    real ``validate`` nuclear branch (species/material/energy rules, deuteron tables, summary)."""

    table_id = "synthetic-nuclear-ci"

    def __init__(self) -> None:
        term = {"n_max": 1.0, "t_lab_max_mev": 250.0}
        self.info = {
            "schema": "synthetic", "builder_version": "synthetic", "npz_sha256": "0" * 64,
            "sources": {}, "elements": {}, "q_plus_table_mev": 0.0, "recoil_t_max_mev": 1.0,
            "transport_energy_bound_mev": 250.0, "history_energy_bound_mev": 250.0,
            "transport_path_bound_terms": {"p": dict(term), "d": dict(term)},
            "grid": {"n_points": 2}, "multiplicity": {"model": "synthetic"},
        }  # fmt: skip

    def material_rows(self, material: Any) -> Any:
        import types

        return types.SimpleNamespace(f_e=1.0, material=material)


SYNTHETIC_NUCLEAR_ID = _SyntheticNuclearTable.table_id


def _patch_synthetic_nuclear(monkeypatch: pytest.MonkeyPatch) -> str:
    from ionmc.nuclear.tables import NuclearTable

    table = _SyntheticNuclearTable()
    monkeypatch.setattr(NuclearTable, "load", classmethod(lambda cls, cache_dir, tid: table))
    return SYNTHETIC_NUCLEAR_ID


def _nuclear_boundary(monkeypatch: pytest.MonkeyPatch) -> str:
    """The nuclear table the cells are evaluated with. LV tier (IONMC_LV_TIER=1): the real cached
    table, whose absence FAILS. Lightweight CI: the synthetic boundary above. Never 'not
    evaluated'."""
    if os.environ.get("IONMC_LV_TIER") == "1":
        tid = _nuclear_table_id_or_none()
        if tid is None:
            pytest.fail("a built nuclear table is required in the local-validation tier")
        return _memoize_real_nuclear(monkeypatch, tid)
    return _patch_synthetic_nuclear(monkeypatch)


def _memoize_real_nuclear(monkeypatch: pytest.MonkeyPatch, tid: str) -> str:
    """The real cached table, loaded once for the sweep instead of once per nuclear cell (the
    cells read the same object ``validate`` would load; any other table id loads as usual)."""
    from ionmc.nuclear.tables import NuclearTable

    real_load = NuclearTable.load
    table = real_load(None, tid)

    def load(cls: Any, cache_dir: Any, table_id: str) -> Any:
        return table if table_id == tid else real_load(cache_dir, table_id)

    monkeypatch.setattr(NuclearTable, "load", classmethod(load))
    return tid


def _cfg(
    species: Any = PROTON,
    energy: float = 20.0,
    sigma: float = 0.0,
    geometry: Any = None,
    scoring: tuple[ScoringGrid, ...] | None = None,
    physics: dict[str, Any] | None = None,
    run: dict[str, Any] | None = None,
    diagnostics: DiagnosticsOptions | None = None,
    tallies: tuple[TallyRequest, ...] = (),
    direction: tuple[float, float, float] = (0.0, 0.0, 1.0),
) -> SimulationConfig:
    ph = {"nuclear": False, "stopping_model": BETHE, "multiple_scattering": species.z < 2}
    rn = {"backend": "python", "precision": "float64", "seed": 5, "n_histories": 40,
          "n_batches": 4}  # fmt: skip
    return SimulationConfig(
        PencilBeamSource(species, (0.0, 0.0, -5.0), direction, energy, sigma),
        geometry or BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 20.0), WATER),
        scoring or (ScoringGrid((-10.0, -10.0, 0.0), (5.0, 5.0, 5.0), (4, 4, 4)),),
        PhysicsOptions(**(ph | (physics or {}))),
        RunOptions(**(rn | (run or {}))),
        diagnostics or DiagnosticsOptions(),
        tallies=tallies,
    )


LET = (TallyRequest("let", "dose", "let_d"),)


# ---- CAP-1a: the sweep -----------------------------------------------------------------------


def _bad(outcomes: list[contract.Outcome]) -> list[contract.Outcome]:
    return [o for o in outcomes if o.verdict in ("third", "disagree")]


def test_cap1a_sweep_declared_equals_enforced(monkeypatch: pytest.MonkeyPatch) -> None:
    tid = _nuclear_boundary(monkeypatch)
    matrix = contract.Matrix(nuclear_table_id=tid)
    assert matrix.size() == 4 * 4 * 3 * 2 * 2 * 2 * 2 * 3 == 2304  # not vacuous
    out = contract.sweep(matrix)
    assert len(out) == matrix.size()
    assert _bad(out) == []  # zero third outcomes and zero disagreements
    counts = contract.summarise(out)
    print("CAP-1a sweep", len(out), counts, "nuclear table:", tid)
    supported = [o for o in out if o.actual == contract.OK]
    refused = [o for o in out if o.actual == contract.UNSUPPORTED]
    unavailable = [o for o in out if o.actual == contract.UNAVAILABLE]
    assert len(supported) > 100 and len(refused) > 500
    # every energy-domain class is exercised on both sides for a proton
    for dc, expect in (("below", False), ("inside", True), ("upper-edge", True), ("above", False)):
        cases = [o for o in out if o.species == "proton" and o.domain_class == dc
                 and o.backend == "python" and o.precision == "float64" and not o.nuclear
                 and o.tallies == "none" and not o.trace and o.model == BETHE]  # fmt: skip
        assert len(cases) == 1 and (cases[0].actual == contract.OK) is expect, (dc, cases)
    # an ion at the edges (domain of decision 0043 section 5)
    edge = [o for o in out if o.species == "alpha" and o.domain_class == "upper-edge"
            and o.backend == "python" and o.precision == "float64" and not o.nuclear
            and o.tallies == "none" and not o.trace and o.model == BETHE]  # fmt: skip
    assert edge[0].actual == contract.OK
    if config.cuda_available():
        assert not unavailable
    else:  # CUDA cases are refused by BackendUnavailableError, never skipped
        cuda_ok = [
            o for o in out if o.backend == "warp-cuda" and o.declared == contract.UNAVAILABLE
        ]
        assert cuda_ok and all(o.actual == contract.UNAVAILABLE for o in cuda_ok)
        assert len(unavailable) == len(cuda_ok)
    assert counts.get("not-evaluated", 0) == 0  # every nuclear cell is evaluated, in every tier
    nuclear_ok = [o for o in out if o.nuclear and o.declared == contract.OK]
    assert nuclear_ok and all(o.actual == contract.OK for o in nuclear_ok)


def test_cap1a_nuclear_cells_on_the_synthetic_boundary(monkeypatch: pytest.MonkeyPatch) -> None:
    """Cache-free in every tier: the nuclear cells (proton, alpha, carbon) run the real validate
    nuclear branch against the synthetic table, so declared-vs-enforced of nuclear is checked in
    the lightweight CI as well."""
    tid = _patch_synthetic_nuclear(monkeypatch)
    matrix = contract.Matrix(nuclear=(True,), backends=("python", "warp-cpu"), nuclear_table_id=tid)
    out = contract.sweep(matrix)
    assert len(out) == matrix.size() == 4 * 4 * 2 * 2 * 1 * 2 * 2 * 3
    assert _bad(out) == []
    assert not [o for o in out if o.verdict == "not-evaluated"]
    ok = [o for o in out if o.declared == contract.OK]
    assert ok and all(o.actual == contract.OK and o.backend == "python" for o in ok)
    refused = {(o.species == "proton", o.backend) for o in out if o.actual == contract.UNSUPPORTED}
    assert (True, "warp-cpu") in refused and (False, "python") in refused  # ion + backend rules


def test_cap1a_declared_run_combinations_equal_validate_including_workers(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tid = _nuclear_boundary(monkeypatch)
    declared = contract.declared_capabilities()["run_combinations"]
    assert len(declared) == 3 * 2 * 2 * 2 * 2 * 2
    seen: set[str] = set()
    for d in declared:
        w = 2 if d["workers_above_one"] else 1
        cfg = _cfg(
            energy=100.0,
            physics={"nuclear": d["nuclear"], "nuclear_table_id": tid if d["nuclear"] else None},
            run={"backend": d["backend"], "precision": d["precision"], "cpu_workers": w},
            diagnostics=DiagnosticsOptions(trace_histories=1 if d["trace"] else 0),
            tallies=LET if d["tallies"] else (),
        )
        try:
            validate(cfg)
            actual = contract.OK
        except (UnsupportedCombinationError, BackendUnavailableError) as exc:
            actual = type(exc).__name__
        assert actual == d["status"], d
        seen.add(actual)
    assert seen >= {contract.OK, contract.UNSUPPORTED}
    assert (contract.UNAVAILABLE in seen) is (not config.cuda_available())


# ---- CAP-1b: cache-backed sweep (LV; the CUDA cases belong to HR) -----------------------------


def test_cap1b_cache_backed_tabulated_models_sweep() -> None:
    """The frozen matrix: species x domain class x backend x precision x nuclear x tally mode
    {none, LET, full set with a synthetic lookup} x trace x tabulated model x material {water,
    PMMA}. LV tier: a missing cache object FAILS (``_require_datasets``), never skips there."""
    _require_datasets("nist-pstar-water-2005", "geant4-icru90-stopping-11.4.2")
    tid = _nuclear_table_id_or_none()
    if tid is None:
        pytest.fail("a built nuclear table is required for the cache-backed matrix")
    matrix = contract.Matrix(
        species=(PROTON, ALPHA, CARBON12),
        models=("nist-pstar-astar-water-2005", "icru90-water-v1"),
        nuclear=(False, True),
        tallies=("none", "let", "full"),
        trace=(False, True),
        materials=(WATER, PMMA),
        nuclear_table_id=tid,
    )
    out = contract.sweep(matrix)
    assert matrix.size() == 3 * 4 * 3 * 2 * 2 * 3 * 2 * 2 * 2 == 3456
    assert len(out) == matrix.size()
    assert _bad(out) == []
    assert not [o for o in out if o.verdict == "not-evaluated"]  # a missing object is a failure
    counts = contract.summarise(out)
    print("CAP-1b sweep", len(out), counts)
    assert counts.get("agree", 0) == len(out)
    assert {o.tallies for o in out} == {"none", "let", "full"}
    assert any(o.actual == contract.OK and o.tallies == "full" for o in out)
    assert any(o.actual == contract.OK and o.material == "water" for o in out)
    assert any(o.material == PMMA.name and o.actual == contract.UNSUPPORTED for o in out)
    assert {o.material for o in out} == {WATER.name, PMMA.name}


DEFAULT_SELECTION = BETHE
PMMA_BOX = BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 20.0), PMMA)
WATER_BOX = BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 20.0), WATER)
WATER_PMMA = VoxelGeometry(
    (-10.0, -10.0, 0.0), (20.0, 20.0, 10.0), (1, 1, 2), (WATER, PMMA),
    material_index=np.array([[[0, 1]]], dtype=np.int32),
)  # fmt: skip


def test_cap1b_the_water_row_names_the_selected_model_with_and_without_water() -> None:
    """The LET water row follows the selected stopping model, whether or not the geometry holds
    water (config.py: the LET medium is part of the tally definition), for each tabulated model
    and for the default."""
    _require_datasets("nist-pstar-water-2005", "geant4-icru90-stopping-11.4.2")
    rows: dict[tuple[str, str], dict[str, Any]] = {}
    phantom = dataclasses.replace(WATER, name="phantom")  # water composition, not named water
    for model in (DEFAULT_SELECTION, "nist-pstar-astar-water-2005", "icru90-water-v1"):
        tabulated = STOPPING_MODELS_KIND(model) == "tabulated"
        cases = [("water", WATER_BOX), ("no-water", PMMA_BOX), ("mixed", WATER_PMMA)]
        for label, geometry in cases:
            cfg = _cfg(geometry=geometry, physics={"stopping_model": model}, tallies=LET)
            if tabulated and label != "water":  # a tabulated model has water data only
                with pytest.raises(UnsupportedCombinationError, match="pmma"):
                    validate(cfg)
                continue
            eff = validate(cfg)
            assert eff.stopping_identity()["model_id_effective"] == model
            row = eff.summary()["tables"]["water_row"]
            assert row["material"] == "water" and row["projectile"] == "proton"
            rows[(model, label)] = row
        # a geometry whose material is not NAMED water: the row is still the selected model's
        eff = validate(_cfg(geometry=BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 20.0), phantom),
                            physics={"stopping_model": model}, tallies=LET))  # fmt: skip
        assert eff.summary()["tables"]["water_row"]["material"] == "water"
        rows[(model, "renamed-water")] = eff.summary()["tables"]["water_row"]
    # the row names the dataset of the selected model, independent of the geometry
    for model, dataset in (("nist-pstar-astar-water-2005", "nist-pstar-water-2005"),
                           ("icru90-water-v1", "geant4-icru90-stopping-11.4.2")):  # fmt: skip
        assert rows[(model, "water")]["dataset_id"] == dataset
        assert rows[(model, "renamed-water")]["dataset_id"] == dataset
    for label in ("water", "no-water", "mixed", "renamed-water"):
        assert (
            rows[(DEFAULT_SELECTION, label)]["source"]
            == rows[(DEFAULT_SELECTION, "water")]["source"]
        )
    assert (
        rows[(DEFAULT_SELECTION, "water")]["source"] != rows[("icru90-water-v1", "water")]["source"]
    )


def STOPPING_MODELS_KIND(model: str) -> str:  # noqa: N802 - reads like the registry lookup
    from ionmc.physics.stopping_models import STOPPING_MODELS

    return str(STOPPING_MODELS[model].kind)


# ---- CAP-1c: every remaining validate rule has a positive/negative pair -----------------------

W = "workers"


def _grids(n: int, names: bool = True, **kw: Any) -> tuple[ScoringGrid, ...]:
    sp = kw.get("spacing", (5.0, 5.0, 5.0))
    return tuple(
        ScoringGrid((-10.0, -10.0, 0.0), sp, (4, 4, 4), f"g{i}" if names else "same")
        for i in range(n)
    )


def _voxel(shape_ok: bool = True) -> VoxelGeometry:
    import numpy as np

    return VoxelGeometry(
        (0.0, 0.0, 0.0), (5.0, 5.0, 5.0), (2, 2, 2), (WATER,),
        np.zeros((2, 2, 2) if shape_ok else (2, 2, 3), dtype=np.int16),
    )  # fmt: skip


# id -> (needle of the validate source site or None, make(bad) -> config, error type)
Make = Callable[[bool], SimulationConfig]
U, B = UnsupportedCombinationError, BackendUnavailableError


def _alpha(bad: bool, **kw: Any) -> SimulationConfig:
    return _cfg(ALPHA, 400.0, **kw)


RULES: list[tuple[str, str | None, Make, type[Exception]]] = [
    # --- the sites of validate, in source order (raise fail / raise BackendUnavailableError)
    (
        "config-type",
        "config must be a SimulationConfig",
        lambda bad: "not a config" if bad else _cfg(),
        U,
    ),  # type: ignore[arg-type,return-value]
    (
        "projectile-registry",
        "is not supported: the transport engine implements",
        lambda bad: _cfg(cio.CUSTOM_PROJECTILE.source.projectile) if bad else _cfg(),
        U,
    ),
    (
        "ion-nuclear",
        "nuclear=True is not available for the ion source",
        lambda bad: _alpha(bad, physics={"nuclear": bad, "nuclear_table_id": "x" if bad else None}),
        U,
    ),
    (
        "ion-e-cut-mev",
        "e_cut_mev applies to protons only",
        lambda bad: _alpha(bad, physics={"e_cut_mev": 3.0 if bad else 2.0}),
        U,
    ),
    (
        "ion-energy-domain",
        "is outside the supported domain",
        lambda bad: _cfg(ALPHA, 100.0 if bad else 400.0, physics={"multiple_scattering": False}),
        U,
    ),
    (
        "straggling-model-choice",
        'choice("straggling_model"',
        lambda bad: _cfg(physics={"straggling_model": "bogus" if bad else "bohr_gamma_v1"}),
        U,
    ),
    (
        "mcs-model-choice",
        'choice("mcs_model"',
        lambda bad: _cfg(physics={"mcs_model": "bogus" if bad else "differential_moliere"}),
        U,
    ),
    (
        "eor-ion-only",
        "'continuous_residual_v1' is available for ion sources only",
        lambda bad: _cfg(
            physics={"end_of_range_model": "continuous_residual_v1" if bad else "point_v1"}
        ),
        U,
    ),  # fmt: skip
    (
        "delta-electrons-choice",
        'choice("delta_electrons"',
        lambda bad: _cfg(physics={"delta_electrons": "bogus" if bad else "local"}),
        U,
    ),
    (
        "python-float32",
        "runs in float64 only",
        lambda bad: _cfg(run={"precision": "float32" if bad else "float64"}),
        U,
    ),
    (
        "cuda-workers",
        "is meaningless on backend 'warp-cuda'",
        lambda bad: _cfg(run={"backend": "warp-cuda" if bad else "warp-cpu", "cpu_workers": 2}),
        U,
    ),
    (
        "seed-range",
        "seed must be in [0, 2**64)",
        lambda bad: _cfg(run={"seed": 2**64 if bad else 5}),
        U,
    ),
    (
        "n-histories-range",
        "n_histories must be in [1, 2**32]",
        lambda bad: _cfg(run={"n_histories": 2**32 + 4 if bad else 40}),
        U,
    ),
    (
        "n-batches-min",
        "n_batches must be >= 2",
        lambda bad: _cfg(run={"n_batches": 1 if bad else 4}),
        U,
    ),
    (
        "n-batches-multiple",
        "must be a multiple of n_batches",
        lambda bad: _cfg(run={"n_histories": 42 if bad else 40}),
        U,
    ),
    (
        "workers-max",
        "cpu_workers must be <=",
        lambda bad: _cfg(
            run={"cpu_workers": config.MAX_CPU_WORKERS + 1 if bad else 2, "n_histories": 400}
        ),
        U,
    ),  # fmt: skip
    (
        "workers-vs-histories",
        "must not exceed n_histories",
        lambda bad: _cfg(run={"cpu_workers": 5 if bad else 2, "n_histories": 4, "n_batches": 2}),
        U,
    ),
    (
        "grids-count",
        "scoring grids are required",
        lambda bad: _cfg(scoring=_grids(config.MAX_SCORING_GRIDS + 1 if bad else 2)),
        U,
    ),
    (
        "grid-names-unique",
        "scoring grid names must be unique",
        lambda bad: _cfg(scoring=_grids(2, names=not bad)),
        U,
    ),
    (
        "memory-budget",
        "above the memory budget",
        lambda bad: _cfg(run={"memory_budget_bytes": 1 if bad else 2**31}),
        U,
    ),
    (
        "scoring-pieces",
        "voxel pieces per leg",
        lambda bad: _cfg(
            scoring=(ScoringGrid((-10.0, -10.0, 0.0), (0.01 if bad else 5.0,) * 3, (1, 1, 1)),),
            physics={"max_step_mm": 100.0},
        ),
        U,
    ),  # fmt: skip
    (
        "table-floor-vs-cutoff",
        "the stopping tables start at",
        lambda bad: _cfg(physics={"e_cut_mev": 1.0 if bad else 2.0}),
        U,
    ),
    (
        "source-energy-window",
        "outside [2 e_cut, table maximum]",
        lambda bad: _cfg(energy=3.0 if bad else 20.0),
        U,
    ),
    (
        "source-sigma-fraction",
        "exceeds",
        lambda bad: _cfg(energy=20.0, sigma=1.5 if bad else 0.5),
        U,
    ),
    (
        "fixed-point-capacity",
        "exceed the capacity",
        lambda bad: _cfg(
            energy=200.0, run={"n_histories": 2**32 if bad else 40, "n_batches": 2 if bad else 4}
        ),
        U,
    ),  # fmt: skip
    (
        "max-steps-counter",
        "block-counter bound",
        lambda bad: _cfg(run={"max_steps": 2**26 if bad else 1000}),
        U,
    ),
    (
        "trace-float64",
        "the per-step trace is only available in float64",
        lambda bad: _cfg(
            run={"backend": "warp-cpu", "precision": "float32" if bad else "float64"},
            diagnostics=DiagnosticsOptions(trace_histories=2),
        ),
        U,
    ),
    (
        "trace-vs-histories",
        "exceeds n_histories",
        lambda bad: _cfg(diagnostics=DiagnosticsOptions(trace_histories=50 if bad else 2)),
        U,
    ),
    (
        "trace-buffer",
        "trace buffers",
        lambda bad: _cfg(
            run={"max_steps": 2**24 if bad else 1000},
            diagnostics=DiagnosticsOptions(trace_histories=40 if bad else 2),
        ),
        U,
    ),
    (
        "channel-backends",
        "does not implement scoring channels yet",
        lambda bad: _cfg(run={"backend": "python" if bad else "warp-cpu"}, tallies=LET),
        U,
    ),  # patched
    (
        "cuda-device",
        "requires a CUDA device",
        lambda bad: _cfg(run={"backend": "warp-cuda" if bad else "warp-cpu"}),
        B,
    ),
    # --- construction-time rules of the config dataclasses (geometry, grid, source, cutoffs, ...)
    (
        "geometry-kind",
        None,
        lambda bad: dataclasses.replace(_cfg(), geometry="box") if bad else _cfg(),
        U,
    ),
    (
        "voxel-material-index-shape",
        None,
        lambda bad: _cfg(
            geometry=_voxel(not bad), scoring=(ScoringGrid((0, 0, 0), (5.0,) * 3, (2, 2, 2)),)
        ),
        U,
    ),  # fmt: skip
    (
        "scoring-grid-shape",
        None,
        lambda bad: _cfg(
            scoring=(ScoringGrid((-10.0, -10.0, 0.0), (5.0,) * 3, (0, 4, 4) if bad else (4, 4, 4)),)
        ),
        U,
    ),  # fmt: skip
    ("source-sigma-negative", None, lambda bad: _cfg(sigma=-1.0 if bad else 0.0), U),
    ("cutoff-positive", None, lambda bad: _cfg(physics={"e_cut_mev": -1.0 if bad else 2.0}), U),
    (
        "energy-loss-fraction-limit",
        None,
        lambda bad: _cfg(physics={"max_energy_loss_fraction": 0.5 if bad else 0.02}),
        U,
    ),
    (
        "end-of-range-mode",
        None,
        lambda bad: _cfg(physics={"end_of_range_model": "bogus" if bad else "default"}),
        U,
    ),
    (
        "mcs-mapping-eval-mode",
        None,
        lambda bad: _alpha(
            bad,
            physics={"mcs_mapping_eval": "bogus" if bad else "root", "multiple_scattering": False},
        ),
        U,
    ),  # fmt: skip
    (
        "chunk-power-of-two",
        None,
        lambda bad: _cfg(run={"chunk_histories": 3000 if bad else 2048}),
        U,
    ),
    ("chunk-minimum", None, lambda bad: _cfg(run={"chunk_histories": 512 if bad else 1024}), U),
    (
        "trace-negative",
        None,
        lambda bad: _cfg(diagnostics=DiagnosticsOptions(trace_histories=-1 if bad else 0)),
        U,
    ),
    ("workers-minimum", None, lambda bad: _cfg(run={"cpu_workers": 0 if bad else 1}), U),
]


def _validate_site_sources() -> list[str]:
    """Source text of every ``raise fail(...)``, ``raise BackendUnavailableError(...)`` and
    ``choice(...)`` statement of ``config.validate``, in source order."""
    src = inspect.getsource(config.validate)
    tree = ast.parse(src)
    sites: list[tuple[int, str]] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Raise) and isinstance(node.exc, ast.Call):
            f = node.exc.func
            if isinstance(f, ast.Name) and f.id in ("fail", "BackendUnavailableError"):
                sites.append((node.lineno, ast.get_source_segment(src, node) or ""))
        elif isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            f = node.value.func
            if isinstance(f, ast.Name) and f.id == "choice":
                sites.append((node.lineno, ast.get_source_segment(src, node) or ""))
    return [s for _, s in sorted(sites)]


def test_cap1c_checklist_covers_every_validate_raise_site_in_source_order() -> None:
    sites = _validate_site_sources()
    listed = [(rid, needle) for rid, needle, _, _ in RULES if needle is not None]
    n_raise = sum(1 for s in sites if s.startswith("raise"))
    print("CAP-1c validate sites", len(sites), "(raise sites:", n_raise, ") checklist", len(RULES))
    assert n_raise == 28 and len(sites) == 31  # a new rule must be added to the checklist
    assert len(listed) == len(sites)
    for (rid, needle), site in zip(listed, sites, strict=True):
        assert needle in site, (rid, needle, site[:80])
    ids = [r[0] for r in RULES]
    assert len(set(ids)) == len(ids)
    assert len(RULES) >= len(sites)


@pytest.mark.parametrize("rule", RULES, ids=[r[0] for r in RULES])
def test_cap1c_positive_and_negative_pair(
    rule: tuple[str, str | None, Make, type[Exception]], monkeypatch: pytest.MonkeyPatch
) -> None:
    rid, _, make, err = rule
    if rid == "channel-backends":
        monkeypatch.setattr(config, "CHANNEL_BACKENDS", ("warp-cpu",))
    if rid == "cuda-device":
        if config.cuda_available():  # the negative case needs an absent device
            validate(make(True))
            return
    # positive: passes validate (the unavailable device is its own positive on warp-cpu)
    validate(make(False))
    with pytest.raises(err):
        validate(make(True))


# ---- CAP-2: the declaration is derived --------------------------------------------------------


def _paths(a: Any, b: Any, prefix: str = "") -> set[str]:
    if isinstance(a, dict) and isinstance(b, dict):
        out: set[str] = set()
        for k in a.keys() | b.keys():
            if k not in a or k not in b:
                out.add(f"{prefix}/{k}")
            else:
                out |= _paths(a[k], b[k], f"{prefix}/{k}")
        return out
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        return {p for i, (x, y) in enumerate(zip(a, b, strict=True))
                for p in _paths(x, y, f"{prefix}[{i}]")}  # fmt: skip
    return set() if a == b else {prefix}


def test_cap2_patching_an_enforcing_constant_changes_the_declaration_identically(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    base = capabilities()
    monkeypatch.setattr(config, "MAX_CPU_WORKERS", 2)
    got = _paths(base, capabilities())
    assert got == {"/max_cpu_workers", "/declared/limits/max_cpu_workers"}
    validate(_cfg(run={"cpu_workers": 2, "n_histories": 400}))  # enforcement moved with it
    with pytest.raises(UnsupportedCombinationError, match="cpu_workers must be <= 2"):
        validate(_cfg(run={"cpu_workers": 3, "n_histories": 400}))
    monkeypatch.undo()

    monkeypatch.setattr(config, "MAX_SCORING_GRIDS", 1)
    got = _paths(base, capabilities())
    assert got == {"/max_scoring_grids", "/declared/limits/max_scoring_grids"}
    with pytest.raises(UnsupportedCombinationError, match="between 1 and 1 scoring grids"):
        validate(_cfg(scoring=_grids(2)))
    monkeypatch.undo()

    monkeypatch.setattr(config, "NUCLEAR_MAX_ENERGY_MEV", 123.0)
    got = _paths(base, capabilities())
    assert "/declared/source_energy_rules/nuclear_max_total_mev" in got
    assert any("window_total_mev_nuclear" in p for p in got)
    assert all(p.startswith(("/declared/", "/species_details")) for p in got), got
    monkeypatch.undo()

    monkeypatch.setitem(config.ION_DOMAIN_MEV_PER_U, "alpha", (60.0, 200.0))
    got = _paths(base, capabilities())
    assert any("/species_details/alpha/energy_domain" in p for p in got)
    assert any("/declared/source_energy_rules/ion_domain_mev_per_u/alpha" in p for p in got)
    validate(_cfg(ALPHA, 4 * 200.0, physics={"multiple_scattering": False}))
    monkeypatch.undo()

    # a channel backend removed from the enforcing tuple changes declaration and enforcement
    monkeypatch.setattr(config, "CHANNEL_BACKENDS", ("python",))
    got = _paths(base, capabilities())
    assert "/declared/limits/channel_backends" in got
    assert all(p.startswith("/declared/") for p in got)
    rows = capabilities()["declared"]["run_combinations"]
    cell = [r for r in rows if r["backend"] == "warp-cpu" and r["precision"] == "float64"
            and r["tallies"] and not r["nuclear"] and not r["trace"]
            and not r["workers_above_one"]]  # fmt: skip
    assert cell[0]["status"] == contract.UNSUPPORTED
    with pytest.raises(UnsupportedCombinationError, match="scoring channels"):
        validate(_cfg(run={"backend": "warp-cpu"}, tallies=LET))


def test_cap2_free_text_only_under_informative() -> None:
    d = capabilities()["declared"]

    def walk(x: Any, path: str) -> None:
        if isinstance(x, dict):
            for k, v in x.items():
                if k != "informative":
                    walk(v, f"{path}/{k}")
        elif isinstance(x, list):
            for i, v in enumerate(x):
                walk(v, f"{path}[{i}]")
        elif isinstance(x, str):
            assert " " not in x and len(x) < 60 or path.endswith("/rule"), (path, x)

    walk(d, "declared")
    assert "informative" in d and isinstance(capabilities()["informative"], dict)
    for key in capabilities()["informative"]["free_text_keys"]:
        assert key.split(".")[0] in capabilities()


# ---- CAP-3: reconciliation --------------------------------------------------------------------


def _summary_of(cfg: SimulationConfig) -> tuple[dict[str, Any], dict[str, Any]]:
    return to_dict(cfg), validate(cfg).summary()


def _seen(rq: Any, eff: Any, into: set[tuple[str, str, str]], ptrs: set[str]) -> None:
    for r in contract.reconcile(rq, eff):
        into.add((r.row, r.cls, r.rule))
        for p in r.pointers:
            assert contract.pointer_resolves(eff, p), (r.row, p)
            ptrs.add(p)


def _nuclear_surgery(rq: dict[str, Any], eff: dict[str, Any]) -> tuple[Any, Any]:
    """A nuclear configuration as dicts (the nuclear setup needs the cache)."""
    rq, eff = copy.deepcopy(rq), copy.deepcopy(eff)
    rq["physics"].update(nuclear=True, nuclear_table_id="tid", e_cut_deuteron_mev=4.0)
    eff["physics"].update(nuclear=True, e_cut_deuteron_mev=4.0)
    eff["nuclear"] = {"table_id": "tid"}
    eff["stopping"]["species_transported"] = ["proton", "deuteron"]
    return rq, eff


def test_cap3a_every_row_class_and_rule_is_exercised_and_recorded() -> None:
    seen: set[tuple[str, str, str]] = set()
    ptrs: set[str] = set()
    slab, _ = load_toml(SLAB_TOML)
    cfgs = [
        _cfg(
            tallies=LET,
            direction=(0.0, 0.0, 2.0),
            physics={"end_of_range_model": "point_v1"},
            run={"max_steps": 1000},
        ),
        _cfg(physics={"straggling": False, "multiple_scattering": False, "stopping_model": LS_EFF}),
        _cfg(run={"backend": "warp-cpu", "precision": "float32"}),
        _cfg(ALPHA, 400.0, physics={"stopping_model": LS_EFF, "multiple_scattering": True}),
        _cfg(ALPHA, 400.0, physics={"multiple_scattering": False}),
        cio.USER_STOPPING,
        slab,
    ]
    for cfg in cfgs:
        rq, eff = _summary_of(cfg)
        _seen(rq, eff, seen, ptrs)
        assert contract.reconcile(cfg, validate(cfg)) == contract.reconcile(rq, eff)
    rq, eff = _summary_of(cfgs[0])
    _seen(*_nuclear_surgery(rq, eff), seen, ptrs)
    table = {(r.id, c, rule) for r in contract.RECONCILIATION_TABLE for c, rule in r.outcomes}
    print("CAP-3a rows", len(contract.RECONCILIATION_TABLE), "outcomes", len(table))
    assert len(contract.RECONCILIATION_TABLE) == 24
    assert table - seen == set(), sorted(table - seen)
    assert seen - table == set()
    assert {c for _, c, _ in seen} == set(contract.CLASSES)
    every = {p for r in contract.RECONCILIATION_TABLE for p in r.pointers}
    assert every - ptrs == set(), sorted(every - ptrs)  # every effective pointer resolves


def test_cap3b_cfg1_generated_set_has_zero_unlisted_differences() -> None:
    slab, _ = load_toml(SLAB_TOML)
    box, _ = load_toml(HERE / "data" / "v3_006" / "proton_box.toml")
    generated = [*cio.API_ROUND_TRIP.values(), cio.USER_STOPPING, cio.ALPHA_BOX,
                 cio._inmemory_lookup_config(), slab, box]  # fmt: skip
    checked = 0
    for cfg in generated:
        try:
            eff = validate(cfg)
        except UnsupportedCombinationError:
            continue  # does not pass validate (custom projectile, missing cache object)
        res = contract.reconcile(cfg, eff)
        assert len(res) == len(contract.RECONCILIATION_TABLE) - 1  # box xor voxel row
        for r in res:
            assert all(contract.pointer_resolves(eff.summary(), p) for p in r.pointers)
        checked += 1
    print("CAP-3b CFG-1 configurations validated and reconciled:", checked, "of", len(generated))
    assert checked >= 6


@pytest.mark.parametrize(
    ("pointer", "value"),
    [
        ("/source/kinetic_energy_mev", 21.0),
        ("/source/position_mm", [0.0, 0.0, -4.0]),
        ("/unit_direction", [0.0, 1.0, 0.0]),
        ("/geometry/shape", [2, 1, 1]),
        ("/scoring/0/spacing_mm", [1.0, 1.0, 1.0]),
        ("/physics/max_step_mm", 0.5),
        ("/physics/e_cut_effective_mev", 3.0),
        ("/physics/end_of_range_model_effective", "continuous_residual_v1"),
        ("/physics/straggling", False),
        ("/physics/ls_effective", "eff"),
        ("/stopping/model_id_effective", "icru90-water-v1"),
        ("/backend", "warp-cpu"),
        ("/production", True),
        ("/max_steps_origin", "user"),
        ("/chunk_histories", 2048),
        ("/seed", 6),
        ("/tallies/quantities/0/name", "other"),
    ],
)
def test_cap3c_a_patched_effective_field_raises(pointer: str, value: Any) -> None:
    rq, eff = _summary_of(_cfg(tallies=LET))
    contract.reconcile(rq, eff)
    node: Any = eff
    parts = pointer.split("/")[1:]
    for p in parts[:-1]:
        node = node[int(p)] if isinstance(node, list) else node[p]
    node[int(parts[-1]) if isinstance(node, list) else parts[-1]] = value
    with pytest.raises(UnsupportedCombinationError, match="requested .*, effective"):
        contract.reconcile(rq, eff)


def test_cap3c_unrequested_feature_objects_and_a_changed_request_raise() -> None:
    rq, eff = _summary_of(_cfg(tallies=LET))
    bad = copy.deepcopy(rq)
    bad["physics"]["multiple_scattering"] = False  # request changed, effective not
    with pytest.raises(UnsupportedCombinationError, match="requested False, effective True"):
        contract.reconcile(bad, eff)
    nuc_rq, nuc_eff = _nuclear_surgery(rq, eff)
    contract.reconcile(nuc_rq, nuc_eff)
    dropped = copy.deepcopy(nuc_eff)
    del dropped["nuclear"]
    with pytest.raises(UnsupportedCombinationError, match="effective <absent>"):
        contract.reconcile(nuc_rq, dropped)
    with pytest.raises(UnsupportedCombinationError, match="nuclear"):
        contract.reconcile(rq, nuc_eff)  # an effective nuclear setup nobody requested


def test_cap3d_decision_0044_selection_is_unchanged() -> None:
    with pytest.raises(UnsupportedCombinationError, match="exactly one of"):
        PhysicsOptions(nuclear=False)
    with pytest.raises(UnsupportedCombinationError, match="exactly one of"):
        PhysicsOptions(nuclear=False, stopping=BetheStoppingSource(), stopping_model=BETHE)
    assert PhysicsOptions(nuclear=False, stopping_model=BETHE).selection.requested_id == BETHE


# ---- CAP-4: discriminating dispatch -----------------------------------------------------------


def _drop_in_summary(
    monkeypatch: pytest.MonkeyPatch, edit: Callable[[dict[str, Any]], None]
) -> None:
    orig = EffectiveConfig.summary

    def summary(self: EffectiveConfig) -> dict[str, Any]:
        out = orig(self)
        edit(out)
        return out

    monkeypatch.setattr(EffectiveConfig, "summary", summary)


def test_cap4_a_consumer_that_drops_the_tallies_makes_simulation_raise(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    Simulation(_cfg(tallies=LET))  # control: constructs
    monkeypatch.setattr(config, "_compile_tallies", lambda *a, **k: None)
    with pytest.raises(UnsupportedCombinationError, match="no channel plan"):
        Simulation(_cfg(tallies=LET))


@pytest.mark.parametrize(
    ("what", "cfg", "edit", "match"),
    [
        ("lookups", "slab", lambda d: d["tallies"]["lookups"].clear(), "lookups"),
        (
            "multiple scattering",
            "box",
            lambda d: d.pop("scattering_length_g_cm2"),
            "multiple_scattering",
        ),
        (
            "straggling",
            "box",
            lambda d: d["physics"].update(straggling_model=None),
            "straggling_model",
        ),
        (
            "stopping model",
            "box",
            lambda d: d["stopping"].pop("model_id_effective"),
            "stopping.model_id_effective",
        ),
        ("tallies", "tallies", lambda d: d.update(tallies=None), "no channel plan"),
    ],
)
def test_cap4_dropping_a_feature_object_from_the_effective_config_raises(
    monkeypatch: pytest.MonkeyPatch, what: str, cfg: str, edit: Callable[[Any], None], match: str
) -> None:
    config_ = {"slab": load_toml(SLAB_TOML)[0], "box": _cfg(), "tallies": _cfg(tallies=LET)}[cfg]
    Simulation(config_)
    _drop_in_summary(monkeypatch, edit)
    with pytest.raises(UnsupportedCombinationError, match=match):
        Simulation(config_)


def test_cap4_a_dropped_nuclear_setup_raises_in_reconcile_and_in_simulation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rq, eff = _summary_of(_cfg())
    nrq, neff = _nuclear_surgery(rq, eff)
    del neff["nuclear"]
    with pytest.raises(UnsupportedCombinationError, match="effective <absent>"):
        contract.reconcile(nrq, neff)
    tid = _patch_synthetic_nuclear(monkeypatch)  # cache-free: runs in the lightweight CI too
    cfg = _cfg(energy=100.0, physics={"nuclear": True, "nuclear_table_id": tid})
    Simulation(cfg)
    _drop_in_summary(monkeypatch, lambda d: d.pop("nuclear"))
    with pytest.raises(UnsupportedCombinationError, match="nuclear"):
        Simulation(cfg)


# ---- CAP-5 / CAP-6: honoured check and persistence --------------------------------------------


@pytest.fixture(scope="module")
def tally_result() -> Any:
    cfg = _cfg(
        tallies=LET,
        physics={"stopping_model": LS_EFF},
        run={"n_histories": 12, "n_batches": 3},
    )
    return Simulation(cfg).run()


@pytest.fixture(scope="module")
def trace_result() -> Any:
    cfg = _cfg(
        run={"n_histories": 12, "n_batches": 3}, diagnostics=DiagnosticsOptions(trace_histories=2)
    )
    return Simulation(cfg).run()


def test_cap5_result_lacking_a_requested_tally_is_refused(
    tmp_path: Path, tally_result: Any
) -> None:
    r = tally_result
    save(r, tmp_path / "ok")  # control
    grid = r.grids[0]
    thinner = dataclasses.replace(grid, quantities={})
    with pytest.raises(ResultFormatError, match="requested tally 'let'.*left no array"):
        save(dataclasses.replace(r, grids=(thinner, *r.grids[1:])), tmp_path / "thin")
    assert not (tmp_path / "thin").exists()


def test_cap5_unrequested_tally_array_and_nuclear_tallies_are_refused(
    tmp_path: Path, tally_result: Any, trace_result: Any
) -> None:
    r = tally_result
    # the same result presented for a request without tallies: reconcile refuses it
    cfg = dataclasses.replace(r.requested_config, tallies=())
    with pytest.raises(UnsupportedCombinationError, match="LET water row"):
        save(dataclasses.replace(r, requested_config=cfg), tmp_path / "extra")
    # nuclear counters on a non-nuclear request
    base = {f.name: getattr(r.counters, f.name) for f in dataclasses.fields(r.counters)}
    nuc = NuclearTransportCounters(**base)
    with pytest.raises(ResultFormatError, match="nuclear request and nuclear counters"):
        save(dataclasses.replace(r, counters=nuc), tmp_path / "nuc")
    # a nuclear request on a result without nuclear tallies
    req = dataclasses.replace(
        r.requested_config,
        physics=dataclasses.replace(r.requested_config.physics, nuclear=True, nuclear_table_id="x"),
    )
    with pytest.raises(ResultFormatError, match="nuclear"):
        save(dataclasses.replace(r, requested_config=req), tmp_path / "nuc2")
    # trace: requested but missing, and present but not requested
    t = trace_result
    save(t, tmp_path / "trace_ok")
    diag = {k: v for k, v in t.diagnostics.items() if k != "trace"}
    with pytest.raises(ResultFormatError, match="trace diagnostics left no array"):
        save(dataclasses.replace(t, diagnostics=diag), tmp_path / "notrace")
    cfg0 = dataclasses.replace(t.requested_config, diagnostics=DiagnosticsOptions())
    with pytest.raises(ResultFormatError, match="no trace was requested"):
        save(dataclasses.replace(t, requested_config=cfg0), tmp_path / "unasked")


def test_cap5_the_reader_refuses_the_same_dishonoured_bundles(
    tmp_path: Path, tally_result: Any, trace_result: Any
) -> None:
    b = save(tally_result, tmp_path / "b")
    load(b)
    _mutate(b, lambda d: d["requested"]["tallies"].clear())
    with pytest.raises(ResultFormatError):
        load(b)
    t = save(trace_result, tmp_path / "t")
    _mutate(t, lambda d: d["requested"]["diagnostics"].update(trace_histories=0))
    with pytest.raises(ResultFormatError, match="trace"):
        load(t)


def test_cap6_requested_and_effective_are_persisted_separately(
    tmp_path: Path, tally_result: Any
) -> None:
    b = save(tally_result, tmp_path / "b")
    doc = json.loads((b / "result.json").read_text(encoding="utf-8"))
    rq, eff = doc["requested"], doc["effective"]
    assert rq["physics"]["end_of_range_model"] == "default"
    assert eff["physics"]["end_of_range_model"] == "default"
    assert eff["physics"]["end_of_range_model_effective"] == "point_v1"
    assert rq["physics"]["stopping_model"] == LS_EFF
    assert eff["stopping"]["model_id_requested"] == LS_EFF
    assert eff["stopping"]["model_id_effective"] == BETHE  # the z = 1 LS no-op
    assert eff["physics"]["ls_requested"] == "eff" and eff["physics"]["ls_effective"] == "off"
    assert rq["run"]["max_steps"] is None
    assert isinstance(eff["max_steps"], int) and eff["max_steps_origin"] == "computed"
    res = {r["row"]: r for r in doc["effective_extra"]["resolutions"]}
    assert len(doc["effective_extra"]["resolutions"]) == len(contract.RECONCILIATION_TABLE) - 1
    assert (res["end-of-range"]["class"], res["end-of-range"]["rule"]) == (
        "derived",
        "default-by-z",
    )
    assert (res["max-steps"]["class"], res["max-steps"]["rule"]) == ("derived", "computed")
    assert (res["ls"]["class"], res["ls"]["rule"]) == ("derived", "ls-z1-off")
    assert res["stopping-model"]["rule"] == "ls-variant-z1-default"
    assert res["run-only"]["class"] == "requested-only"
    assert doc["effective_extra"]["resolutions_status"].startswith("computed")
    # the persisted resolutions are what the single implementation computes, and load binds them
    assert doc["effective_extra"]["resolutions"] == [
        r.to_dict() for r in contract.reconcile(rq, eff)
    ]
    load(b)
    _mutate(b, lambda d: d["effective_extra"]["resolutions"].pop())
    with pytest.raises(ResultFormatError, match="resolutions"):
        load(b)


def test_cap6_simulation_keeps_the_resolutions_and_result_gains_no_field() -> None:
    sim = Simulation(_cfg(tallies=LET))
    assert [r.row for r in sim.resolutions] == [
        r.id for r in contract.RECONCILIATION_TABLE if r.id != "geometry-voxel"
    ]
    from ionmc.simulation import Result

    assert not [f.name for f in dataclasses.fields(Result) if "resolution" in f.name]


# ---- CAP-3c: every request field class is compared field by field (Codex REVIEW-57661746) -----


def _mutation_base() -> tuple[SimulationConfig, dict[str, Any], dict[str, Any], dict[str, Any]]:
    lk = (contract._synthetic_lookup(),)
    grids = (
        ScoringGrid((-10.0, -10.0, 0.0), (5.0, 5.0, 5.0), (4, 4, 4), "dose"),
        ScoringGrid((-10.0, -10.0, 0.0), (10.0, 10.0, 10.0), (2, 2, 2), "other"),
    )
    tallies = (
        TallyRequest("let", "dose", "let_d"),
        TallyRequest("ls", "dose", "lookup_sum", lookup=lk[0].name),
        TallyRequest("spec", "other", "fluence_spectrum",
                     energy_edges_mev_per_u=(1.0, 101.0, 201.0, 301.0, 401.0)),
    )  # fmt: skip
    cfg = dataclasses.replace(_cfg(scoring=grids, tallies=tallies), lookups=lk)
    eff = validate(cfg)
    return cfg, to_dict(cfg), eff.summary(), contract.evidence_from_effective(eff)


def _mutations() -> list[
    tuple[str, Callable[[dict[str, Any], dict[str, Any], dict[str, Any]], None]]
]:
    from ionmc.configio import _material_dict

    return [
        ("tally quantity", lambda rq, eff, ev: rq["tallies"][0].update(quantity="let_t")),
        ("tally species", lambda rq, eff, ev: rq["tallies"][0].update(species=["alpha"])),
        ("tally generation", lambda rq, eff, ev: rq["tallies"][0].update(generation="primary")),
        ("tally grid", lambda rq, eff, ev: rq["tallies"][0].update(grid="other")),
        ("tally lookup binding",
         lambda rq, eff, ev: rq["tallies"][1].update(lookup="elsewhere")),
        ("tally spectrum edges",
         lambda rq, eff, ev: rq["tallies"][2].update(
             energy_edges_mev_per_u=[2.0, 102.0, 202.0, 302.0, 402.0])),
        ("tally removed and replaced",
         lambda rq, eff, ev: rq["tallies"].__setitem__(
             0, {**rq["tallies"][0], "name": "renamed"})),
        ("lookup units", lambda rq, eff, ev: rq["lookups"][0].update(units="2")),
        ("lookup content", lambda rq, eff, ev: rq["lookups"][0].update(content_sha256="f" * 64)),
        ("effective lookup content",
         lambda rq, eff, ev: eff["tallies"]["lookups"][0].update(content_sha256="e" * 64)),
        ("effective channel generation",
         lambda rq, eff, ev: eff["tallies"]["channels"][0].update(generation=[1, 1])),
        ("effective channel species",
         lambda rq, eff, ev: eff["tallies"]["channels"][0].update(species=[])),
        ("effective quantity units",
         lambda rq, eff, ev: eff["tallies"]["quantities"][0].update(units="Gy")),
        ("diagnostics end positions",
         lambda rq, eff, ev: rq["diagnostics"].update(track_end_positions=True)),
        ("diagnostics escape records",
         lambda rq, eff, ev: rq["diagnostics"].update(escape_records=True)),
        ("diagnostics trace", lambda rq, eff, ev: rq["diagnostics"].update(trace_histories=2)),
        ("realised trace count (evidence)",
         lambda rq, eff, ev: ev["diagnostics"].update(trace_count=1)),
        ("tally let_medium", lambda rq, eff, ev: rq["tallies"][0].update(let_medium="bone")),
        ("tally dose_reference",
         lambda rq, eff, ev: rq["tallies"][0].update(dose_reference="water")),
        ("effective LET water row medium",
         lambda rq, eff, ev: eff["tables"]["water_row"].update(material="pmma")),
        ("realised diagnostics (evidence)",
         lambda rq, eff, ev: ev["diagnostics"].update(escape_records=True)),
        ("material definition",
         lambda rq, eff, ev: rq["materials"].__setitem__(
             0, _material_dict(dataclasses.replace(WATER, density_g_cm3=1.1)))),
        ("effective material identity",
         lambda rq, eff, ev: eff["tables"]["materials"][0].update(material_sha256="d" * 64)),
    ]  # fmt: skip


def test_cap3c_reconcile_compares_every_request_field_class() -> None:
    _cfg_, rq0, eff0, ev0 = _mutation_base()
    contract.reconcile(rq0, eff0, ev0)  # the control passes
    assert {r.row for r in contract.reconcile(rq0, eff0, ev0)} >= {"tallies", "lookups"}
    for label, mutate in _mutations():
        rq, eff, ev = copy.deepcopy(rq0), copy.deepcopy(eff0), copy.deepcopy(ev0)
        mutate(rq, eff, ev)
        with pytest.raises(UnsupportedCombinationError):
            contract.reconcile(rq, eff, ev)
        print("CAP-3c refused:", label)


def test_cap3c_trace_count_is_reconciled_exactly() -> None:
    """A request between positive trace counts (1 -> 2) is detected: the evidence carries the
    realised count min(trace_histories, n_histories), compared exactly, never as a Boolean."""
    cfg = _cfg(diagnostics=DiagnosticsOptions(trace_histories=1), run={"n_histories": 40})
    eff = validate(cfg)
    rq, summary, ev = to_dict(cfg), eff.summary(), contract.evidence_from_effective(eff)
    assert ev["diagnostics"]["trace_count"] == 1
    contract.reconcile(rq, summary, ev)
    rq2 = copy.deepcopy(rq)
    rq2["diagnostics"]["trace_histories"] = 2
    with pytest.raises(UnsupportedCombinationError, match="trace_count: requested 2, effective 1"):
        contract.reconcile(rq2, summary, ev)
    ev2 = copy.deepcopy(ev)
    ev2["diagnostics"]["trace_count"] = 2
    with pytest.raises(UnsupportedCombinationError, match="trace_count: requested 1, effective 2"):
        contract.reconcile(rq, summary, ev2)
    # more traces than histories: validate refuses the request, and so does the reconciliation
    # (uncapped: a request can never grow past the realised count)
    cfg3 = _cfg(diagnostics=DiagnosticsOptions(trace_histories=100), run={"n_histories": 40})
    with pytest.raises(UnsupportedCombinationError, match="exceeds n_histories"):
        validate(cfg3)
    with pytest.raises(UnsupportedCombinationError, match="at most n_histories 40"):
        contract.requested_trace_count(41, 40)
    assert contract.requested_trace_count(40, 40) == 40
    assert contract.requested_trace_count(0, 40) is None
    rq4 = copy.deepcopy(rq)
    rq4["diagnostics"]["trace_histories"] = 41
    with pytest.raises(UnsupportedCombinationError, match="at most n_histories 40"):
        contract.reconcile(rq4, summary, ev)
    # a negative count (DiagnosticsOptions forbids it) is refused, never read as "no trace"
    with pytest.raises(UnsupportedCombinationError, match="requested -1, effective 'a count >= 0'"):
        contract.requested_trace_count(-1, 40)
    cfg0 = _cfg(run={"n_histories": 40})
    eff0 = validate(cfg0)
    rq5, ev5 = to_dict(cfg0), contract.evidence_from_effective(eff0)
    contract.reconcile(rq5, eff0.summary(), ev5)
    rq5["diagnostics"]["trace_histories"] = -1
    with pytest.raises(UnsupportedCombinationError, match="trace_histories: requested -1"):
        contract.reconcile(rq5, eff0.summary(), ev5)
    # from persisted origins: the length of trace_end_history
    ev3 = contract.evidence_from_origins(
        [
            ({"kind": "diagnostic", "array": "trace_end_history"}, None, (7,)),
            ({"kind": "diagnostic", "array": "trace_end_code"}, None, (7,)),
        ]
    )
    assert ev3["diagnostics"] == {"track_end_positions": False, "escape_records": False,
                                  "trace_histories": True, "trace_count": 7}  # fmt: skip


def test_cap5_writer_and_reader_refuse_a_changed_trace_count_and_tally_semantics(
    tmp_path: Path, tally_result: Any, trace_result: Any
) -> None:
    """Regression mutations for REVIEW-d77d34e0: the exact trace count and the let_medium /
    dose_reference semantics are checked by the writer (save) and the reader (load)."""
    t = trace_result  # trace_histories=2 of 12 histories
    more = dataclasses.replace(
        t.requested_config, diagnostics=DiagnosticsOptions(trace_histories=3)
    )
    with pytest.raises((ResultFormatError, UnsupportedCombinationError), match="trace"):
        save(dataclasses.replace(t, requested_config=more), tmp_path / "trace3")
    assert not (tmp_path / "trace3").exists()
    b = save(t, tmp_path / "trace_ok")
    _mutate(b, lambda d: d["requested"]["diagnostics"].update(trace_histories=3))
    with pytest.raises(ResultFormatError, match="trace"):
        load(b)
    _mutate(b, lambda d: d["requested"]["diagnostics"].update(trace_histories=1))
    with pytest.raises(ResultFormatError, match="trace"):
        load(b)
    # every history traced (N of N): the request N + 1 is refused by the writer and the reader
    full = Simulation(_cfg(run={"n_histories": 4, "n_batches": 2},
                           diagnostics=DiagnosticsOptions(trace_histories=4))).run()  # fmt: skip
    assert full.diagnostics["trace_end_history"].shape == (4,)
    plus = dataclasses.replace(
        full.requested_config, diagnostics=DiagnosticsOptions(trace_histories=5)
    )
    with pytest.raises((ResultFormatError, UnsupportedCombinationError), match="trace"):
        save(dataclasses.replace(full, requested_config=plus), tmp_path / "trace5")
    assert not (tmp_path / "trace5").exists()
    bf = save(full, tmp_path / "trace_full")
    load(bf)
    _mutate(bf, lambda d: d["requested"]["diagnostics"].update(trace_histories=5))
    with pytest.raises(ResultFormatError, match="trace"):
        load(bf)
    # 0 -> -1 on a bundle without a trace: refused by the reader
    r = tally_result
    b0 = save(r, tmp_path / "notrace")
    load(b0)
    _mutate(b0, lambda d: d["requested"]["diagnostics"].update(trace_histories=-1))
    with pytest.raises(ResultFormatError, match="trace"):
        load(b0)
    for field, value in (("let_medium", "bone"), ("dose_reference", "water")):
        req = dataclasses.replace(
            r.requested_config,
            tallies=tuple(
                dataclasses.replace(q, **{field: value}) for q in r.requested_config.tallies
            ),
        )
        with pytest.raises((ResultFormatError, UnsupportedCombinationError), match=field):
            save(dataclasses.replace(r, requested_config=req), tmp_path / f"bad_{field}")
        assert not (tmp_path / f"bad_{field}").exists()
        b2 = save(r, tmp_path / f"ok_{field}")
        _mutate(b2, lambda d, f=field, v=value: d["requested"]["tallies"][0].update({f: v}))
        with pytest.raises(ResultFormatError, match=field):
            load(b2)


def test_cap3c_voxel_array_digests_are_compared() -> None:
    cfg, _ = load_toml(SLAB_TOML)
    eff = validate(cfg)
    rq, summary, ev = to_dict(cfg), eff.summary(), contract.evidence_from_effective(eff)
    contract.reconcile(rq, summary, ev)
    assert ev["geometry_digests"]["material_index"] and ev["geometry_digests"]["density"]
    for role in ("material_index", "density"):
        bad = copy.deepcopy(ev)
        bad["geometry_digests"][role] = "c" * 64  # the arrays the engine holds differ
        with pytest.raises(UnsupportedCombinationError, match=role):
            contract.reconcile(rq, summary, bad)
        rq2 = copy.deepcopy(rq)
        rq2["geometry"][role]["content_sha256"] = "b" * 64  # the request names other arrays
        with pytest.raises(UnsupportedCombinationError, match=role):
            contract.reconcile(rq2, summary, ev)


def test_cap5_the_writer_compares_the_full_quantity_identity(
    tmp_path: Path, tally_result: Any
) -> None:
    r = tally_result
    grid = r.grids[0]
    q = grid.quantities["let"]
    for field, value in (("quantity", "let_t"), ("units", "Gy"), ("definition", "other text")):
        altered = dataclasses.replace(q, **{field: value})
        bad = dataclasses.replace(grid, quantities={**grid.quantities, "let": altered})
        with pytest.raises(ResultFormatError, match="not the requested quantity|left no array"):
            save(dataclasses.replace(r, grids=(bad, *r.grids[1:])), tmp_path / f"bad_{field}")
        assert not (tmp_path / f"bad_{field}").exists()


# ---- the LS memo of the Bethe tables (pure memoisation, bit-identical values) -----------------


def test_ls_memo_returns_bit_identical_private_copies() -> None:
    """``stopping.delta_l_ls_memo`` returns the uncached value bit for bit, a fresh array on every
    call (a caller's mutation never reaches the memo), stays bounded, and is bypassed for
    unequal shapes (the broadcast path of ``delta_l_ls``)."""
    from ionmc.physics import lindhard_sorensen, stopping

    memo = stopping._LsMemo(entries=2)
    q = np.full(5, 6.0)
    beta = np.linspace(0.2, 0.6, 5)
    ref = lindhard_sorensen.delta_l_ls(q, beta)
    first = memo(q, beta)
    second = memo(q, beta)
    assert np.array_equal(first, ref) and np.array_equal(second, ref)
    assert first is not second and (memo.hits, memo.misses) == (1, 1)
    first[0] = 0.0
    assert np.array_equal(memo(q, beta), ref)  # the memoised array is private
    memo(np.full(5, 2.0), beta)
    memo(np.full(5, 8.0), beta)
    assert len(memo._values) == 2  # bounded
    one = memo(np.array([6.0]), beta)  # unequal shapes: broadcast by the uncached function
    assert np.array_equal(one, lindhard_sorensen.delta_l_ls(6.0, beta))
    assert memo.misses == 3 and len(memo._values) == 2  # the bypass is neither hit nor miss
