"""Consumer consistency of the stopping model selection (plan v3-013 A.2; decision 0044).

A test-only model ``test-bethe-x1.05`` (tabulated kind: S = 1.05 times the default's S on the
default grid, start range = the default's / 1.05) is selected by ``stopping_model`` through the
test hook :func:`~ionmc.physics.stopping_models.register_test_model`. If every consumer takes its
stopping power from the one selected source, every length-scaled quantity of the test run is the
default's divided by 1.05 (rows CC-R, CC-LET, CC-RES, CC-T); the identity rows CC-ID, CC-NEG and
CC-HASH check that every table identity names the test model and that a consumer bypassing the
selection is detected. Seeds: plan v3-013 section C, ``ion13`` r_index 10: 20613004 + k.
No reference data and no table values: the test model derives from the analytic default.
"""

from __future__ import annotations

import json
from collections.abc import Iterator
from contextlib import contextmanager
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from ionmc import config as config_module
from ionmc.config import (
    DiagnosticsOptions,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
    validate,
)
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import VoxelGeometry
from ionmc.materials import PMMA, WATER
from ionmc.physics.projectiles import ALPHA, CARBON12, DEUTERON, OXYGEN16, PROTON, Projectile
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.physics.stopping_evidence import STOPPING_EVIDENCE
from ionmc.physics.stopping_models import (
    DEFAULT_STOPPING_MODEL_ID,
    STOPPING_MODELS,
    register_test_model,
    resolve,
    spec_of,
    spec_sha256,
)
from ionmc.physics.stopping_selection import (
    effective_model_id,
    select,
    table_model_id,
)
from ionmc.physics.stopping_tabulated import TabulatedStoppingSource
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.simulation import Simulation
from ionmc.transport.reference import END_CUTOFF, _Reference
from ionmc.transport.tables import TransportTables, thaw

TEST_ID = "test-bethe-x1.05"
FACTOR = 1.05
SEED0 = 20613004  # ion13 family, r_index 10 (the A.2 test seed): 20612004 + 100 * 10
SPECIES = (PROTON, DEUTERON, ALPHA, CARBON12, OXYGEN16)
ALL_LS = ("bethe-ionmc-v1", "bethe-ionmc-v1+ls-eff", "bethe-ionmc-v1+ls-bare_mixed")


def make_test_source() -> TabulatedStoppingSource:
    """S = 1.05 x the default's S at the default grid nodes (water; the five default species),
    grid bounds and density of the default, start range = the default's constant-S start range
    divided by 1.05."""
    default = resolve(DEFAULT_STOPPING_MODEL_ID)
    water = {p.name: default.table(WATER, p) for p in SPECIES}
    arrays = {n: (t.energy_per_u, FACTOR * t.s_el_mass) for n, t in water.items()}
    r0 = {n: float(t.csda_range_g_cm2[0]) / FACTOR for n, t in water.items()}
    return TabulatedStoppingSource(
        name=TEST_ID,
        metadata={"model_id": TEST_ID, "version": "1", "decision": "test"},
        arrays=arrays,
        e_min_per_u=1.0,
        e_max_per_u=500.0,
        start_range_id="default constant-S start range / 1.05",
        start_range=r0,
    )


@contextmanager
def _registered_test_model() -> Iterator[None]:
    with register_test_model(make_test_source().instance_spec(), make_test_source):
        yield


@pytest.fixture
def registered() -> Iterator[None]:
    with _registered_test_model():
        yield


# ---- configuration builders ----


def _cfg(
    proj: Projectile,
    model: str,
    scale: float = 1.0,
    *,
    backend: str = "python",
    e_per_u: float | None = None,
    n_histories: int = 2000,
    seed: int = SEED0,
    mcs: bool = False,
    straggling: bool = False,
    tallies: tuple[TallyRequest, ...] = (),
    diagnostics: bool = False,
) -> SimulationConfig:
    """Small water box; EVERY length of the request is multiplied by ``scale`` (voxel spacing, box
    and world length, source position, ``z_clip`` = ``z_exit_mm``, ``max_step_mm``,
    ``range_rho_f_mm``, the scoring grid); energies and dimensionless fractions are not."""
    e_u = (100.0 if e_per_u is None else e_per_u) if proj.z > 1 else 100.0
    shape = (3, 3, 20)
    geo = VoxelGeometry(
        (-7.5 * scale, -7.5 * scale, 0.0),
        (5.0 * scale, 5.0 * scale, 5.0 * scale),
        shape,
        (WATER,),
        np.zeros(shape, dtype=np.int32),
        z_exit_mm=95.0 * scale,
    )
    from ionmc.sources import PencilBeamSource  # noqa: PLC0415

    return SimulationConfig(
        source=PencilBeamSource(
            proj, (0.0, 0.0, 2.0 * scale), (0.0, 0.0, 1.0), e_u * proj.a, 0.0, 0.0
        ),
        geometry=geo,
        scoring=(
            ScoringGrid(
                (-7.5 * scale, -7.5 * scale, 0.0),
                (2.5 * scale, 2.5 * scale, 5.0 * scale),
                (6, 6, 20),
            ),
        ),
        physics=PhysicsOptions(
            nuclear=False,
            stopping_model=model,
            multiple_scattering=mcs,
            straggling=straggling,
            max_step_mm=1.0 * scale,
            range_rho_f_mm=0.1 * scale,
        ),
        run=RunOptions(
            backend=backend,  # type: ignore[arg-type]
            precision="float64",
            seed=seed,
            n_histories=n_histories,
            n_batches=20 if n_histories >= 20 else 2,
        ),
        tallies=tallies,
        diagnostics=DiagnosticsOptions(track_end_positions=diagnostics),
    )


# ---- the generic source and the hook ----


def test_hook_registers_and_restores_and_refuses_collisions() -> None:
    before = (dict(STOPPING_MODELS), dict(STOPPING_EVIDENCE))
    spec = make_test_source().instance_spec()
    with register_test_model(spec, make_test_source):
        assert TEST_ID in STOPPING_MODELS
        digest = spec_sha256(STOPPING_MODELS[TEST_ID])
        assert [r.status for r in STOPPING_EVIDENCE[digest]] == ["unvalidated"]
        assert resolve(TEST_ID).name == TEST_ID
        with pytest.raises(ValueError, match="already registered"):
            with register_test_model(spec, make_test_source):
                pass
        assert TEST_ID in STOPPING_MODELS  # the refused nested registration did not remove it
    assert (dict(STOPPING_MODELS), dict(STOPPING_EVIDENCE)) == before
    with pytest.raises(ValueError, match="unknown stopping model id"):
        resolve(TEST_ID)
    default = STOPPING_MODELS[DEFAULT_STOPPING_MODEL_ID]
    with pytest.raises(ValueError, match="start with 'test-'"):  # a registered id cannot be reused
        with register_test_model(default, make_test_source):
            pass
    from dataclasses import replace  # noqa: PLC0415

    with pytest.raises(ValueError, match="start with 'test-'"):
        with register_test_model(replace(spec, model_id="x-bethe"), make_test_source):
            pass
    try:  # the registry is restored also when the body raises
        with register_test_model(spec, make_test_source):
            raise RuntimeError("body")
    except RuntimeError:
        pass
    assert TEST_ID not in STOPPING_MODELS


def test_tabulated_source_identity_and_fail_closed_domain(registered: None) -> None:
    src = resolve(TEST_ID)
    mid, spec = spec_of(src)
    assert mid == TEST_ID and spec is STOPPING_MODELS[TEST_ID]
    assert spec.kind == "tabulated" and spec.data == ()
    # a source with different content is never reported under the registered id
    other = make_test_source()
    arrays = {n: (e, 1.01 * s) for n, (e, s) in other._arrays.items()}
    changed = TabulatedStoppingSource(
        other.name, other.metadata, arrays, 1.0, 500.0, other.start_range_id, other.start_range
    )
    assert spec_of(changed)[0] == "unregistered"
    assert spec_of(changed)[1].decision == "user"
    narrower = TabulatedStoppingSource(
        other.name,
        other.metadata,
        other._arrays,
        2.0,
        500.0,
        other.start_range_id,
        other.start_range,
    )
    assert spec_of(narrower)[0] == "unregistered"
    # no extrapolation, no data for the projectile, structural material rule
    with pytest.raises(ValueError, match="outside the array span"):
        TabulatedStoppingSource(
            other.name, other.metadata, other._arrays, 0.5, 500.0, "x", other.start_range
        )
    with pytest.raises(ValueError, match="no data for projectile"):
        TabulatedStoppingSource(
            other.name,
            other.metadata,
            {"proton": other._arrays["proton"]},
            1.0,
            500.0,
            "x",
            other.start_range,
        ).table(WATER, ALPHA)
    with pytest.raises(ValueError, match="liquid water only"):
        src.table(PMMA, PROTON)
    with pytest.raises(ValueError, match="model_id"):
        TabulatedStoppingSource("n", {"version": "1"}, other._arrays, 1.0, 500.0, "x", {})


# ---- CC-R ----


@pytest.mark.parametrize("proj", [PROTON, ALPHA], ids=lambda p: p.name)
def test_cc_r_csda_range_is_the_defaults_over_the_factor(
    proj: Projectile, registered: None
) -> None:
    e = np.geomspace(2.0, 480.0, 10)  # ten energies per nucleon
    r_def = resolve(DEFAULT_STOPPING_MODEL_ID).table(WATER, proj).range_at(e)
    r_new = resolve(TEST_ID).table(WATER, proj).range_at(e)
    assert np.max(np.abs(r_new * FACTOR / r_def - 1.0)) <= 1e-12
    # the start range is part of the table, not only the increments
    t_def = resolve(DEFAULT_STOPPING_MODEL_ID).table(WATER, proj)
    t_new = resolve(TEST_ID).table(WATER, proj)
    np.testing.assert_allclose(t_new.csda_range_g_cm2 * FACTOR, t_def.csda_range_g_cm2, rtol=1e-12)
    assert t_new.metadata["range_construction"] == "exact-loglog-quadrature-v1"


# ---- CC-LET, CC-ID ----


def test_cc_let_and_cc_id_water_row_and_identities(registered: None) -> None:
    tallies = (TallyRequest("lt", "dose", "let_t"),)
    for proj in (PROTON, ALPHA):
        d = validate(_cfg(proj, DEFAULT_STOPPING_MODEL_ID, tallies=tallies))
        t = validate(_cfg(proj, TEST_ID, tallies=tallies))
        assert t.tables.has_water and d.tables.has_water
        ratio = np.exp(t.tables.water_ln_s_mass - d.tables.water_ln_s_mass)
        assert np.max(np.abs(ratio / FACTOR - 1.0)) <= 1e-12  # CC-LET at every node
        ratio_m = np.exp(t.tables.ln_s_mass - d.tables.ln_s_mass)
        assert np.max(np.abs(ratio_m / FACTOR - 1.0)) <= 1e-12
        blk = t.summary()["stopping"]
        assert blk["model_id_requested"] == blk["model_id_effective"] == TEST_ID
        assert blk["kind"] == "tabulated" and blk["version"] == "1"
        assert blk["spec_sha256"] == spec_sha256(STOPPING_MODELS[TEST_ID])
        assert [r["status"] for r in blk["evidence_at_run"]] == ["unvalidated"]
        entries = [*thaw(t.tables.identity), thaw(t.tables.water_identity)]
        assert len(entries) == 2
        for entry in entries:
            assert entry["metadata"]["model_id"] == TEST_ID
            assert entry["metadata"]["version"] == "1"
            assert table_model_id(entry, proj) == TEST_ID
        for entry in [*thaw(d.tables.identity), thaw(d.tables.water_identity)]:
            assert "model_id" not in entry["metadata"]  # pre-existing metadata unchanged
        assert t.tables.sha256 != d.tables.sha256


# ---- CC-RES ----


@pytest.mark.parametrize("proj", [PROTON, ALPHA], ids=lambda p: p.name)
def test_cc_res_residual_range_segment(proj: Projectile, registered: None) -> None:
    d = _Reference(validate(_cfg(proj, DEFAULT_STOPPING_MODEL_ID, n_histories=20)))
    t = _Reference(validate(_cfg(proj, TEST_ID, n_histories=20)))
    e_cut = d.e_cut
    assert e_cut == t.e_cut
    floor = float(np.exp(d.tab.ln_e0[0]))
    energies = [e_cut, 0.9 * e_cut, 0.7 * e_cut, 0.5 * e_cut, 0.25 * e_cut]
    assert any(e < floor for e in energies) or proj.z > 1  # the constant-floor part is exercised
    for e in energies:
        r_d, r_t = d._residual_range(0, e), t._residual_range(0, e)
        assert abs(r_t * FACTOR / r_d - 1.0) <= 1e-12, e
        # the inverse of the segment: the energy at the scaled remaining range is the same
        assert abs(t._residual_energy(0, r_d / FACTOR) / d._residual_energy(0, r_d) - 1.0) <= 1e-12


# ---- CC-T ----


def _mean_cutoff_depth(cfg: SimulationConfig) -> tuple[float, float]:
    res = Simulation(cfg).run()
    diag = res.diagnostics
    code, z = diag["end_code"], diag["end_position_mm"][:, 2]
    assert np.all(code == END_CUTOFF)
    assert res.effective_config.summary()["stopping"]["model_id_effective"] in {
        DEFAULT_STOPPING_MODEL_ID,
        TEST_ID,
    }
    return float(np.mean(z)), float(np.max(np.abs(z - z[0])))


@pytest.mark.parametrize("backend", ["python", "warp-cpu"])
@pytest.mark.parametrize(("proj", "k"), [(PROTON, 0), (ALPHA, 1)], ids=lambda x: str(x))
def test_cc_t_mean_cutoff_depth_scales_with_the_range(
    backend: str, proj: Projectile, k: int, registered: None
) -> None:
    seed = SEED0 + k
    z_def, spread_def = _mean_cutoff_depth(
        _cfg(proj, DEFAULT_STOPPING_MODEL_ID, backend=backend, seed=seed, diagnostics=True)
    )
    z_new, spread_new = _mean_cutoff_depth(
        _cfg(proj, TEST_ID, 1.0 / FACTOR, backend=backend, seed=seed, diagnostics=True)
    )
    residual = abs(z_new * FACTOR / z_def - 1.0)
    print(
        f"CC-T {backend} {proj.name}: z_default={z_def!r} z_test={z_new!r} residual={residual:.3e}"
    )
    assert z_def > 50.0  # the cutoff point is deep inside the box, past many steps
    assert residual <= 1e-9


# ---- CC-NEG ----


def _default_tables(proj: Projectile) -> TransportTables:
    src = BetheStoppingSource()
    return TransportTables.from_stopping_tables(
        [src.table(WATER, proj)], water=src.table(WATER, proj), projectile=proj
    )


def test_cc_neg_each_consumer_bypassing_the_selection_is_detected(
    monkeypatch: pytest.MonkeyPatch, registered: None
) -> None:
    tallies = (TallyRequest("lt", "dose", "let_t"),)
    ok = _cfg(PROTON, TEST_ID, tallies=tallies)
    validate(ok)  # control: the unpatched run validates
    default_src = BetheStoppingSource()
    real_tables, real_water = config_module._stopping_tables, config_module._water_table
    real_from = TransportTables.from_stopping_tables

    # 1. material tables from the default source (no water row, no ion rows involved)
    with monkeypatch.context() as m:
        m.setattr(
            config_module,
            "_stopping_tables",
            lambda source, materials, projectile=PROTON: real_tables(
                default_src, materials, projectile
            ),
        )
        with pytest.raises(UnsupportedCombinationError, match=r"material tables"):
            validate(ok)
    # 2. the LET water row from the default source
    with monkeypatch.context() as m:
        m.setattr(
            config_module,
            "_water_table",
            lambda config, projectile=None: default_src.table(
                WATER, config.source.projectile if projectile is None else projectile
            ),
        )
        with pytest.raises(UnsupportedCombinationError, match=r"LET water row"):
            validate(ok)
    # 3. the exact-Omega ion rows (alpha, multiple scattering on) from the default source: the
    # material tables are repaired to the test source so that only that group is wrong
    ion_cfg = _cfg(ALPHA, TEST_ID, mcs=True)
    test_tables = real_tables(resolve(TEST_ID), (WATER,), ALPHA)

    def right_transport_tables(tables: Any, *args: Any, **kwargs: Any) -> TransportTables:
        return real_from(test_tables, *args, **kwargs)

    with monkeypatch.context() as m:
        m.setattr(
            config_module,
            "_stopping_tables",
            lambda source, materials, projectile=PROTON: real_tables(
                default_src, materials, projectile
            ),
        )
        m.setattr(TransportTables, "from_stopping_tables", right_transport_tables)
        with pytest.raises(UnsupportedCombinationError, match=r"exact-Omega ion tables"):
            validate(ion_cfg)
    assert real_water is config_module._water_table  # patches were undone
    # 4. deuteron tables of a nuclear run (the nuclear table cache is not needed: the assertion is
    # applied to the deuteron tables a nuclear setup would hand over)
    sel = select(None, TEST_ID)
    eff = validate(ok)
    deut_test = real_from(
        [resolve(TEST_ID).table(WATER, DEUTERON)],
        water=resolve(TEST_ID).table(WATER, DEUTERON),
        projectile=DEUTERON,
    )
    good = SimpleNamespace(deuteron_tables=deut_test)
    config_module._assert_single_model(sel, PROTON, eff.tables, good)  # type: ignore[arg-type]
    bad = SimpleNamespace(deuteron_tables=_default_tables(DEUTERON))
    with pytest.raises(UnsupportedCombinationError, match=r"deuteron"):
        config_module._assert_single_model(sel, PROTON, eff.tables, bad)  # type: ignore[arg-type]


# ---- CC-HASH part 2 ----


def _key(entry: dict[str, Any], proj: Projectile) -> str:
    return json.dumps(
        [entry["metadata"], entry["e_min_mev"], entry["e_max_mev"], proj.name], sort_keys=True
    )


def test_cc_hash_metadata_and_injective_identity_map(registered: None) -> None:
    default_meta = resolve(DEFAULT_STOPPING_MODEL_ID).table(WATER, PROTON).metadata
    assert "model_id" not in default_meta and "version" not in default_meta
    t_def = _default_tables(PROTON)
    new_src = resolve(TEST_ID)
    t_new = TransportTables.from_stopping_tables(
        [new_src.table(WATER, PROTON)], water=new_src.table(WATER, PROTON), projectile=PROTON
    )
    assert t_new.sha256 != t_def.sha256
    assert thaw(t_new.identity)[0]["metadata"]["model_id"] == TEST_ID
    assert thaw(t_new.identity)[0]["metadata"]["version"] == "1"

    key_to_id: dict[str, str] = {}
    pair_to_key: dict[tuple[str, str], str] = {}
    for mid in (*ALL_LS, TEST_ID):
        src = resolve(mid)
        for proj in SPECIES:
            if proj.name not in STOPPING_MODELS[mid].domain.species:
                continue
            tab = src.table(WATER, proj)
            entry = {
                "metadata": json.loads(json.dumps(tab.metadata, default=str)),
                "e_min_mev": float(tab.energy_per_u[0] * proj.a),
                "e_max_mev": float(tab.energy_per_u[-1] * proj.a),
            }
            eff_id = effective_model_id(mid, proj)
            assert table_model_id(entry, proj) == eff_id, (mid, proj.name)
            key = _key(entry, proj)
            assert key_to_id.setdefault(key, eff_id) == eff_id  # a function of (metadata, species)
            pair = (eff_id, proj.name)
            assert pair_to_key.setdefault(pair, key) == key
    # injective on the effective (model, species) pairs: distinct pairs, distinct keys
    assert len(set(pair_to_key.values())) == len(pair_to_key)
    # z = 1 variants of the Lindhard-Sorensen options share the default's key (they resolve to it)
    assert effective_model_id("bethe-ionmc-v1+ls-eff", PROTON) == DEFAULT_STOPPING_MODEL_ID
    assert (DEFAULT_STOPPING_MODEL_ID, "proton") in pair_to_key
    assert (TEST_ID, "alpha") in pair_to_key
