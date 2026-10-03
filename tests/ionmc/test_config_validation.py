"""C1 (fail-closed configuration rules) and checks of geometry, source and scoring types."""

from collections.abc import Callable
from dataclasses import replace
from typing import Any

import numpy as np
import pytest

from ionmc.config import (
    DiagnosticsOptions,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
    validate,
)
from ionmc.errors import (
    TransportLimitError,
    UnsupportedCombinationError,
)
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.materials import AIR, PMMA, WATER
from ionmc.physics.projectiles import DEUTERON, PROTON
from ionmc.physics.stopping import BetheStoppingSource, StoppingSource, StoppingTable
from ionmc.scoring import (
    MEV_PER_G_TO_GY,
    ScoringGrid,
    overlap_matrix,
    reduce_batches,
    voxel_mass_g,
)
from ionmc.simulation import Simulation
from ionmc.sources import PencilBeamSource

MakeConfig = Callable[..., SimulationConfig]


def _rejects(cfg_or_callable: Any, exc: type[Exception] = UnsupportedCombinationError) -> None:
    with pytest.raises(exc):
        if callable(cfg_or_callable):
            cfg_or_callable()
        else:
            validate(cfg_or_callable)


class _NoTableSource:
    """A stopping source that has no table for any material."""

    name = "no-tables"

    def table(self, material: Any, projectile: Any) -> StoppingTable:
        raise ValueError("this source has no table for the material")


def test_c1_valid_configuration_records_requested_and_effective(
    make_config: MakeConfig,
) -> None:
    cfg = make_config(direction=(0.0, 3.0, 4.0))
    eff = validate(cfg)
    assert eff.requested is cfg
    assert eff.unit_direction == pytest.approx((0.0, 0.6, 0.8))
    assert eff.backend == "python" and eff.precision == "float64" and not eff.production
    assert eff.max_steps > 100 and eff.max_steps_origin == "computed"
    assert eff.rng["u01"] == "((w >> 8) + 0.5) * 2**-24"
    s = eff.summary()
    assert s["tables"]["sha256"] == eff.tables.sha256
    assert s["scattering_length_g_cm2"][0] == pytest.approx(46.88, rel=3e-3)
    assert s["physics"]["scattering_E_s_mev"] == 15.0
    pinned = validate(make_config(max_steps=5000))
    assert pinned.max_steps == 5000 and pinned.max_steps_origin == "user"


def test_c1_nuclear_has_no_default_and_true_is_rejected(
    make_config: MakeConfig, bethe: BetheStoppingSource
) -> None:
    with pytest.raises(TypeError):
        PhysicsOptions(stopping=bethe)  # type: ignore[call-arg]
    with pytest.raises(UnsupportedCombinationError):
        PhysicsOptions(nuclear=1, stopping=bethe)  # type: ignore[arg-type]
    cfg = make_config()
    _rejects(replace(cfg, physics=replace(cfg.physics, nuclear=True)))
    # the same rule through Simulation, before any transport
    _rejects(lambda: Simulation(replace(cfg, physics=replace(cfg.physics, nuclear=True))))


def test_c1_projectile_must_be_a_proton(make_config: MakeConfig) -> None:
    cfg = make_config()
    _rejects(replace(cfg, source=replace(cfg.source, projectile=DEUTERON)))


@pytest.mark.parametrize("field", ["straggling_model", "mcs_model", "delta_electrons"])
def test_c1_unknown_models_are_rejected(make_config: MakeConfig, field: str) -> None:
    cfg = make_config()
    _rejects(replace(cfg, physics=replace(cfg.physics, **{field: "highland"})))


def test_c1_backend_and_precision_rules(make_config: MakeConfig) -> None:
    cfg = make_config()
    run = cfg.run
    _rejects(replace(cfg, run=replace(run, precision="float32")))  # python is float64 only
    with pytest.raises(UnsupportedCombinationError):
        RunOptions(backend="gpu", precision="float64", seed=1, n_histories=2)  # type: ignore[arg-type]
    with pytest.raises(UnsupportedCombinationError):
        RunOptions(backend="python", precision="float16", seed=1, n_histories=2)  # type: ignore[arg-type]
    # warp-cpu is available in both precisions and validates
    for precision in ("float32", "float64"):
        eff = validate(replace(cfg, run=replace(run, backend="warp-cpu", precision=precision)))
        assert eff.backend == "warp-cpu" and eff.production == (precision == "float32")
    # a warp-cuda request never falls back when no device exists; cpu_workers on CUDA is rejected
    _rejects(
        replace(cfg, run=replace(run, backend="warp-cuda", precision="float32", cpu_workers=2))
    )
    # multi-process execution is supported on python and warp-cpu
    assert validate(replace(cfg, run=replace(run, n_histories=4, cpu_workers=2))).requested
    # float32 with a trace, oversized traces, chunk and worker bounds
    _rejects(
        replace(
            cfg,
            run=replace(run, backend="warp-cpu", precision="float32"),
            diagnostics=DiagnosticsOptions(trace_histories=1),
        )
    )
    _rejects(replace(cfg, diagnostics=DiagnosticsOptions(trace_histories=3)))  # > n_histories
    _rejects(
        replace(
            cfg,
            run=replace(run, n_histories=2, n_batches=2, max_steps=10**8),
            diagnostics=DiagnosticsOptions(trace_histories=2),
        )
    )  # 1 GiB trace limit
    _rejects(replace(cfg, run=replace(run, cpu_workers=3)))  # more workers than histories
    _rejects(replace(cfg, run=replace(run, n_histories=2000, n_batches=2, cpu_workers=1000)))
    _rejects(
        replace(
            cfg,
            run=replace(
                run, memory_budget_bytes=8 * 30 * 30 * 30 * 2 * 3 - 1, n_histories=6, cpu_workers=3
            ),
        )
    )  # one copy per worker
    with pytest.raises(UnsupportedCombinationError):
        RunOptions(backend="python", precision="float64", seed=1, n_histories=2, cpu_workers=0)


@pytest.mark.parametrize(
    "run_kwargs",
    [
        {"n_histories": 0},
        {"n_histories": 2**32 + 20},
        {"n_histories": 7, "n_batches": 2},
        {"n_histories": 4, "n_batches": 1},
        {"seed": -1},
        {"seed": 2**64},
    ],
)
def test_c1_run_size_and_seed_bounds(make_config: MakeConfig, run_kwargs: dict[str, int]) -> None:
    cfg = make_config()
    _rejects(replace(cfg, run=replace(cfg.run, **run_kwargs)))
    for bad in (1.5, True, "3"):
        with pytest.raises(UnsupportedCombinationError):
            RunOptions(backend="python", precision="float64", seed=bad, n_histories=2)  # type: ignore[arg-type]


def test_c1_source_energy_rules(make_config: MakeConfig) -> None:
    cfg = make_config()
    for energy in (3.9, 600.0):  # below 2 e_cut, above the table maximum (500 MeV)
        _rejects(replace(cfg, source=replace(cfg.source, kinetic_energy_mev=energy)))
    _rejects(replace(cfg, source=replace(cfg.source, energy_sigma_mev=1.01)))  # > 5 % of 20 MeV
    validate(replace(cfg, source=replace(cfg.source, energy_sigma_mev=1.0)))
    with pytest.raises(UnsupportedCombinationError):
        PencilBeamSource(PROTON, (0, 0, 0), (0, 0, 0), 20.0)  # zero direction
    for bad in (float("nan"), float("inf"), -1.0, 0.0):
        with pytest.raises(UnsupportedCombinationError):
            PencilBeamSource(PROTON, (0, 0, 0), (0, 0, 1), bad)
    with pytest.raises(UnsupportedCombinationError):
        PencilBeamSource(PROTON, (0, 0), (0, 0, 1), 20.0)  # type: ignore[arg-type]
    with pytest.raises(UnsupportedCombinationError):
        PencilBeamSource(PROTON, (0, 0, 0), (0, 0, 1), 20.0, energy_sigma_mev=-1.0)
    with pytest.raises(UnsupportedCombinationError):
        PencilBeamSource("proton", (0, 0, 0), (0, 0, 1), 20.0)  # type: ignore[arg-type]


def test_c1_table_rules(make_config: MakeConfig, bethe: BetheStoppingSource) -> None:
    cfg = make_config()
    # table floor 1.5 MeV is above half of e_cut = 2 MeV
    floor = BetheStoppingSource(e_min_per_u=1.5)
    _rejects(replace(cfg, physics=replace(cfg.physics, stopping=floor)))
    # a material without a table
    _rejects(replace(cfg, physics=replace(cfg.physics, stopping=_NoTableSource())))
    with pytest.raises(UnsupportedCombinationError):
        PhysicsOptions(nuclear=False, stopping=object())  # type: ignore[arg-type]
    src: StoppingSource = bethe
    assert src.name == "bethe"


def _geo(**kw: Any) -> VoxelGeometry:
    args: dict[str, Any] = dict(
        origin_mm=(0.0, 0.0, 0.0),
        spacing_mm=(1.0, 1.0, 1.0),
        shape=(2, 2, 2),
        materials=(WATER,),
        material_index=np.zeros((2, 2, 2), dtype=np.int32),
    )
    args.update(kw)
    return VoxelGeometry(**args)


def test_c1_geometry_rules() -> None:
    _geo()
    _rejects(lambda: _geo(material_index=np.full((2, 2, 2), 1, dtype=np.int32)))  # out of range
    _rejects(lambda: _geo(material_index=np.full((2, 2, 2), -1, dtype=np.int32)))
    _rejects(lambda: _geo(material_index=np.zeros((2, 2, 3), dtype=np.int32)))  # shape
    _rejects(lambda: _geo(material_index=np.zeros((2, 2, 2), dtype=np.float64)))
    for bad in (0.0, -1.0, float("nan")):
        _rejects(lambda bad=bad: _geo(spacing_mm=(1.0, bad, 1.0)))
    _rejects(lambda: _geo(shape=(2, 0, 2)))
    _rejects(lambda: _geo(origin_mm=(0.0, float("inf"), 0.0)))
    _rejects(lambda: _geo(materials=()))
    _rejects(lambda: _geo(materials=("water",)))
    dens = np.ones((2, 2, 2))
    _geo(density_g_cm3=dens)
    for bad_density in (0.0, -1.0, float("nan")):
        d = dens.copy()
        d[1, 1, 1] = bad_density
        _rejects(lambda d=d: _geo(density_g_cm3=d))  # vacuum voxels are unsupported
    _rejects(lambda: _geo(density_g_cm3=np.ones((2, 2, 3))))
    # a BoxPhantom is the one-voxel geometry
    box = BoxPhantom((1.0, 2.0, 3.0), (4.0, 5.0, 6.0), WATER)
    g = box.to_geometry()
    assert g.shape == (1, 1, 1) and g.upper_mm == (5.0, 7.0, 9.0) and g.n_voxels == 1
    assert g.densities_g_cm3()[0, 0, 0] == 1.0
    _rejects(lambda: BoxPhantom((0, 0, 0), (1.0, 0.0, 1.0), WATER))
    _rejects(lambda: BoxPhantom((0, 0, 0), (1.0, 1.0, 1.0), "water"))  # type: ignore[arg-type]
    assert not g.material_index.flags.writeable


def test_c1_scoring_rules(make_config: MakeConfig) -> None:
    cfg = make_config()
    grid = cfg.scoring[0]
    five = tuple(replace(grid, name=f"g{i}") for i in range(5))
    _rejects(replace(cfg, scoring=five))  # more than four grids
    _rejects(replace(cfg, scoring=()))
    _rejects(replace(cfg, scoring=(grid, grid)))  # duplicate names
    # memory budget: n_batches * voxels * 8 B for the float64 python backend
    tight = replace(cfg.run, memory_budget_bytes=2 * grid.n_voxels * 8 - 1)
    _rejects(replace(cfg, run=tight))
    validate(replace(cfg, run=replace(cfg.run, memory_budget_bytes=2 * grid.n_voxels * 8)))
    # max_step above the smallest scoring spacing: rejected, not clamped
    _rejects(replace(cfg, physics=replace(cfg.physics, max_step_mm=2.5)))
    fine = ScoringGrid((0.0, 0.0, 0.0), (2.0, 2.0, 0.5), (4, 4, 4))
    _rejects(replace(cfg, scoring=(grid, fine)))
    for bad in (0.0, -1.0, float("nan")):
        with pytest.raises(UnsupportedCombinationError):
            ScoringGrid((0, 0, 0), (1.0, bad, 1.0), (2, 2, 2))
    with pytest.raises(UnsupportedCombinationError):
        ScoringGrid((0, 0, 0), (1, 1, 1), (2, 2))  # type: ignore[arg-type]
    with pytest.raises(UnsupportedCombinationError):
        ScoringGrid((0, 0, 0), (1, 1, 1), (2, 2, 2.5))  # type: ignore[arg-type]
    with pytest.raises(UnsupportedCombinationError):
        ScoringGrid((0, 0, 0), (1, 1, 1), (2, 2, 2), name="")


def test_c1_step_counter_bound(make_config: MakeConfig) -> None:
    cfg = make_config()
    ok = 2**32 // 65  # largest step bound with 65 * max_steps < 2**32
    validate(replace(cfg, run=replace(cfg.run, max_steps=ok)))
    _rejects(replace(cfg, run=replace(cfg.run, max_steps=ok + 1)))  # (1 + 64) * max_steps >= 2^32


def test_c1_numeric_option_bounds(bethe: BetheStoppingSource) -> None:
    ok = dict(nuclear=False, stopping=bethe)
    PhysicsOptions(**ok)
    for kwargs in (
        {"e_cut_mev": 0.0},
        {"e_cut_mev": float("nan")},
        {"max_step_mm": -1.0},
        {"max_energy_loss_fraction": 0.0},
        {"max_energy_loss_fraction": 0.3},
        {"range_alpha": 1.5},
        {"range_rho_f_mm": 0.0},
        {"short_step_fraction": 0.5},
        {"straggling": "yes"},
        {"straggling_model": 3},
    ):
        with pytest.raises(UnsupportedCombinationError):
            PhysicsOptions(**ok, **kwargs)  # type: ignore[arg-type]
    for dkw in ({"trace_histories": -1}, {"track_end_positions": 1}, {"escape_records": "x"}):
        with pytest.raises(UnsupportedCombinationError):
            DiagnosticsOptions(**dkw)  # type: ignore[arg-type]
    with pytest.raises(UnsupportedCombinationError):
        RunOptions("python", "float64", 1, 2, max_steps=0)
    with pytest.raises(UnsupportedCombinationError):
        RunOptions("python", "float64", 1, 2, worker_timeout_s=-1.0)
    with pytest.raises(UnsupportedCombinationError):
        RunOptions("python", "float64", 1, 2, allow_invalid_result=1)  # type: ignore[arg-type]
    with pytest.raises(UnsupportedCombinationError):
        SimulationConfig(  # type: ignore[arg-type]
            source=None,
            geometry=BoxPhantom((0, 0, 0), (1, 1, 1), WATER),
            scoring=(),
            physics=PhysicsOptions(**ok),
            run=RunOptions("python", "float64", 1, 2),
        )


def test_transport_limit_error_carries_result(make_config: MakeConfig) -> None:
    err = TransportLimitError("x")
    assert err.result is None and isinstance(err, RuntimeError)


# ---------------------------------------------------------------------------------------------
# scoring: voxel masses and the batch estimator
# ---------------------------------------------------------------------------------------------
def test_voxel_mass_is_the_exact_overlap_integral() -> None:
    idx = np.array([[[0, 1]]], dtype=np.int32)  # two voxels along z: water, PMMA
    geo = VoxelGeometry((0, 0, 0), (2.0, 2.0, 3.0), (1, 1, 2), (WATER, PMMA), idx)
    # a grid shifted by half a voxel in z, finer in x/y and larger than the geometry in x
    grid = ScoringGrid((-1.0, 0.0, 1.5), (1.0, 1.0, 3.0), (4, 2, 2))
    mass = voxel_mass_g(grid, geo)
    assert mass.shape == (4, 2, 2)
    # voxel (1, 0, 0): x 0..1 (inside), y 0..1, z 1.5..4.5 -> 1.5 mm water + 1.5 mm PMMA
    expected = (1.0 * 1.0 * (1.5 * 1.0 + 1.5 * 1.19)) * 1e-3
    assert mass[1, 0, 0] == pytest.approx(expected, rel=1e-14)
    assert mass[0, :, :].sum() == 0.0  # x -1..0 lies outside the geometry
    # the total mass inside the grid equals the mass of the covered part of the geometry
    covered = 2.0 * 2.0 * (1.5 * 1.0 + 3.0 * 1.19) * 1e-3  # z 1.5..6 of the geometry
    assert mass.sum() == pytest.approx(covered, rel=1e-14)
    full = ScoringGrid((0.0, 0.0, 0.0), (2.0, 2.0, 3.0), (1, 1, 2))
    assert voxel_mass_g(full, geo).sum() == pytest.approx(
        2.0 * 2.0 * (3.0 * 1.0 + 3.0 * 1.19) * 1e-3, rel=1e-14
    )
    # density override
    dens = np.array([[[1.0, 0.5]]])
    geo2 = VoxelGeometry((0, 0, 0), (2.0, 2.0, 3.0), (1, 1, 2), (WATER, AIR), idx, dens)
    assert voxel_mass_g(full, geo2)[0, 0, 1] == pytest.approx(2.0 * 2.0 * 3.0 * 0.5e-3)
    ov = overlap_matrix(np.array([0.0, 1.0, 2.0]), np.array([0.5, 1.5]))
    assert np.allclose(ov, [[0.5], [0.5]])
    assert MEV_PER_G_TO_GY == 1.602176634e-10


def test_batch_reducer() -> None:
    sums = np.array([[2.0, 0.0, 4.0], [4.0, 0.0, 0.0], [6.0, 0.0, 2.0], [0.0, 0.0, 2.0]])
    st = reduce_batches(sums, 2)
    x = sums / 2.0
    assert np.allclose(st.mean, x.mean(axis=0))
    assert np.allclose(st.variance_of_mean, x.var(axis=0, ddof=1) / 4.0)
    assert list(st.n_nonzero) == [3, 0, 3]
    with pytest.raises(ValueError):
        reduce_batches(sums[:1], 2)
    with pytest.raises(ValueError):
        reduce_batches(sums, 0)
    with pytest.raises(ValueError):
        reduce_batches(np.full((2, 2), np.nan), 1)


# ---------------------------------------------------------------------------------------------
# review fixes: projectile identity, table identity, effective table provenance
# ---------------------------------------------------------------------------------------------
def test_forged_proton_like_projectiles_are_rejected(make_config: MakeConfig) -> None:
    from ionmc.physics.projectiles import Projectile

    cfg = make_config()
    # wrong mass or name: constructible (consistent fields) but not the canonical proton
    for forged in (
        Projectile("proton", "p", 1, 1, 940.0),
        Projectile("antiproton", "p", 1, 1, PROTON.mass_mev),
        Projectile("proton", "H", 1, 1, PROTON.mass_mev),
    ):
        _rejects(replace(cfg, source=replace(cfg.source, projectile=forged)))
    # invalid fields are rejected when the projectile is built
    for args in (
        ("proton", "p", 1, 1, float("nan")),
        ("proton", "p", 1, 1, float("inf")),
        ("proton", "p", 1, 1, -938.0),
        ("proton", "p", 1, 1, 0.0),
        ("proton", "p", 1, 1, 5000.0),  # mass inconsistent with the mass number
        ("", "p", 1, 1, PROTON.mass_mev),
        ("proton", " ", 1, 1, PROTON.mass_mev),
        ("proton", "p", 0, 1, PROTON.mass_mev),
        ("proton", "p", 2, 1, PROTON.mass_mev),
        ("proton", "p", 1.0, 1, PROTON.mass_mev),
        ("proton", "p", True, 1, PROTON.mass_mev),
    ):
        with pytest.raises(ValueError):
            Projectile(*args)  # type: ignore[arg-type]
    _rejects(lambda: Simulation(replace(cfg, source=replace(cfg.source, projectile=DEUTERON))))


class _WrongTableSource:
    """Returns a valid table, but not for the requested material or projectile."""

    def __init__(self, material: Any = AIR, projectile: Any = PROTON) -> None:
        self.name = "wrong-table"
        self._material, self._projectile = material, projectile

    def table(self, material: Any, projectile: Any) -> StoppingTable:
        return BetheStoppingSource().table(self._material, self._projectile)


def test_stopping_source_results_are_checked_against_the_request(make_config: MakeConfig) -> None:
    from ionmc.physics.projectiles import ALPHA

    cfg = make_config()
    for source in (_WrongTableSource(AIR), _WrongTableSource(WATER, ALPHA)):
        _rejects(replace(cfg, physics=replace(cfg.physics, stopping=source)))
    # same name but a different I value or density is also a different material
    for forged in (WATER.with_I(75.0, "water"), replace_density(WATER, 1.1)):
        _rejects(replace(cfg, physics=replace(cfg.physics, stopping=_WrongTableSource(forged))))
    validate(replace(cfg, physics=replace(cfg.physics, stopping=_WrongTableSource(WATER))))
    _rejects(
        lambda: validate(replace(cfg, physics=replace(cfg.physics, stopping=_NotATableSource())))
    )


class _NotATableSource:
    name = "not-a-table"

    def table(self, material: Any, projectile: Any) -> Any:
        return {"energy": [1.0]}


def replace_density(material: Any, density: float) -> Any:
    from ionmc.materials import Material

    return Material(
        material.name,
        density,
        dict(material.mass_fractions),
        material.I_eV,
        material.sternheimer,
        material.source,
    )


def _synthetic_pstar() -> Any:
    """Smooth STAR-like proton table from the analytic model at I = 75 eV (plumbing only)."""
    import numpy as np

    from ionmc.data.nist_star import parse_star_text
    from ionmc.physics.stopping import bethe_mass_stopping

    e = np.geomspace(1.0, 600.0, 120)
    s = bethe_mass_stopping(e, PROTON, WATER.with_I(75.0))
    r = np.cumsum(np.concatenate(([0.0], np.diff(e) / s[1:]))) + 1e-3
    lines = ["PSTAR: x", "WATER, LIQUID", "", "h", "h", "h", ""]
    for ei, si, ri in zip(e, s, r, strict=True):
        lines.append(f"{ei:.6E} {si:.6E} 1.0E-3 {si:.6E} {ri:.6E} {ri * 0.9:.6E} 1.0")
    return parse_star_text("\n".join(lines))


def test_effective_table_provenance_reaches_the_summary(make_config: MakeConfig) -> None:
    """A NIST STAR table is at I = 75 eV and carries a dataset identity: the effective I, the
    dataset id and the hashes are in the table identity, the table hash and the summary, and
    they differ from the requested material's 78 eV and from the analytic tables."""
    from ionmc.physics.stopping import NistStarStoppingSource

    star = replace(
        _synthetic_pstar(), dataset_id="nist-pstar-water-2005", sha256="ab" * 32, version="v-test"
    )
    cfg = make_config()
    star_cfg = replace(cfg, physics=replace(cfg.physics, stopping=NistStarStoppingSource(star)))
    eff, ref = validate(star_cfg), validate(cfg)
    ident = eff.summary()["tables"]["materials"][0]
    assert ident["I_eV_requested"] == 78.0 and ident["I_eV_effective"] == 75.0
    assert ident["dataset_id"] == "nist-pstar-water-2005" and ident["source_sha256"] == "ab" * 32
    assert ident["content_sha256"] == star.content_sha256 and ident["source"] == "nist-star"
    assert ident["metadata"]["effective_material"]["I_eV"] == 75.0
    assert ident["metadata"]["requested_material"]["I_eV"] == 78.0
    ref_ident = ref.summary()["tables"]["materials"][0]
    assert ref_ident["I_eV_effective"] == 78.0 and ref_ident["dataset_id"] is None
    assert eff.tables.sha256 != ref.tables.sha256
    other = replace(star, sha256="cd" * 32)
    other_cfg = replace(cfg, physics=replace(cfg.physics, stopping=NistStarStoppingSource(other)))
    assert validate(other_cfg).tables.sha256 != eff.tables.sha256  # the dataset hash is hashed
    import json

    json.dumps(eff.summary())  # JSON-serialisable


def test_material_identity_is_the_full_definition(make_config: MakeConfig) -> None:
    """Same name, density, Z/A and I but different density-effect parameters or composition
    must not be accepted as the requested material."""
    from ionmc.materials import Material, SternheimerParameters

    cfg = make_config()
    sp = WATER.sternheimer
    assert sp is not None
    other_sp = SternheimerParameters(sp.x0, sp.x1, sp.cbar + 0.01, sp.a, sp.m)
    same_aggregates = Material(
        "water", 1.0, dict(WATER.mass_fractions), WATER.I_eV, other_sp, WATER.source
    )
    assert (same_aggregates.density_g_cm3, same_aggregates.z_over_a, same_aggregates.I_eV) == (
        WATER.density_g_cm3,
        WATER.z_over_a,
        WATER.I_eV,
    )
    recomposed = Material(
        "water",
        1.0,
        {"H": 0.1119, "O": 0.8881},  # a different composition
        WATER.I_eV,
        WATER.sternheimer,
        WATER.source,
    )
    for forged in (same_aggregates, recomposed):
        _rejects(replace(cfg, physics=replace(cfg.physics, stopping=_WrongTableSource(forged))))
    validate(replace(cfg, physics=replace(cfg.physics, stopping=_WrongTableSource(WATER))))
    from ionmc.transport.tables import material_fingerprint

    assert material_fingerprint(WATER) != material_fingerprint(same_aggregates)
    assert material_fingerprint(WATER) != material_fingerprint(recomposed)
    assert validate(cfg).tables.identity[0]["material_sha256"] == material_fingerprint(WATER)


def test_transport_tables_are_deeply_immutable(make_config: MakeConfig) -> None:
    import copy

    eff = validate(make_config())
    t = eff.tables
    for name in (
        "ln_s_mass",
        "ln_r_mass",
        "ln_e_of_r",
        "e_min_mev",
        "z_over_a",
        "inv_rho_xs_cm2_g",
        "nominal_density_g_cm3",
        "ln_e0",
    ):
        arr = getattr(t, name)
        assert not arr.flags.writeable
        with pytest.raises(ValueError):
            arr[...] = 0.0
        with pytest.raises(ValueError):
            arr.setflags(write=True)  # cannot be made writable again
        before = arr.tobytes()
        with pytest.raises(ValueError):
            arr[...] = 1.0
        with pytest.raises(ValueError):
            arr.setflags(write=True)
        assert arr.dtype == np.float64 and arr.tobytes() == before
        root = arr
        while isinstance(root, np.ndarray):
            assert not root.flags.writeable
            with pytest.raises(ValueError):
                root.setflags(write=True)
            root = root.base
        assert isinstance(root, bytes)  # immutable backing buffer
    assert t.sha256 == eff.tables.sha256 == validate(make_config()).tables.sha256
    with pytest.raises(TypeError):
        t.identity[0]["I_eV_effective"] = 1.0  # type: ignore[index]
    with pytest.raises(TypeError):
        t.identity[0]["metadata"]["source"] = "x"  # type: ignore[index]
    with pytest.raises(AttributeError):
        t.sha256 = "x"  # type: ignore[misc]
    # any constructor path re-freezes: the caller's array stays writable, ours does not
    base = t.ln_s_mass.copy()
    view = base[:, :]
    from dataclasses import fields

    kwargs = {f.name: getattr(t, f.name) for f in fields(t)}
    kwargs["ln_s_mass"] = view
    t2 = type(t)(**kwargs)
    assert not t2.ln_s_mass.flags.writeable and base.flags.writeable
    # the summary is a deep copy: mutating it changes neither identity nor hash
    s = eff.summary()
    before = copy.deepcopy(s)
    s["tables"]["materials"][0]["I_eV_effective"] = -1.0
    s["tables"]["materials"][0]["metadata"]["source"] = "tampered"
    s["rng"]["generator"] = "x"
    assert eff.summary() == before
    assert t.identity[0]["I_eV_effective"] == 78.0 and eff.rng["generator"] == "philox4x32-10"


def test_voxel_geometry_arrays_are_deeply_immutable() -> None:
    mi = np.zeros((2, 1, 1), dtype=np.int64)
    dens = np.array([[[1.0]], [[1.1]]])
    g = VoxelGeometry((0, 0, 0), (1, 1, 1), (2, 1, 1), (WATER,), mi, dens)
    assert mi.flags.writeable and dens.flags.writeable  # caller's arrays untouched
    for arr, dtype in ((g.material_index, np.int32), (g.density_g_cm3, np.float64)):
        assert arr is not None and arr.dtype == dtype and arr.shape == (2, 1, 1)
        before = arr.tobytes()
        with pytest.raises(ValueError):
            arr[...] = 0
        root: object = arr
        while isinstance(root, np.ndarray):
            assert not root.flags.writeable
            with pytest.raises(ValueError):
                root.setflags(write=True)
            root = root.base
        assert isinstance(root, bytes)
        assert arr.tobytes() == before
    assert g.densities_g_cm3().tolist() == [[[1.0]], [[1.1]]]


def test_c1_cuda_backend_unavailable_raises_without_fallback(
    make_config: MakeConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    from ionmc.errors import BackendUnavailableError

    monkeypatch.setattr("ionmc.config.cuda_available", lambda: False)
    cfg = make_config(backend="warp-cuda", precision="float32")
    with pytest.raises(BackendUnavailableError, match="no fallback|nothing falls back"):
        validate(cfg)
    monkeypatch.setattr("ionmc.config.cuda_available", lambda: True)
    assert validate(cfg).backend == "warp-cuda"
    assert "requested" in validate(cfg).__dict__ and validate(cfg).requested is cfg
