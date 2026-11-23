"""Rows CFG-1..CFG-7 and the manifest part of PROV-8 (plan v3-006, step S1).

``ionmc.configio``: the canonical dict form, its strict inverse, the TOML loader, the
``InputManifest`` and the dataset checks. Nothing here runs transport: configurations are built,
serialised, rebuilt and (CFG-1, CFG-2) validated. Fixtures are synthetic (``data/v3_006``: tiny
geometry arrays and the mathematical test lookup of ``tests/data``); no table values, no seed of the
``out6`` family and no ``cpu_workers`` above 1 appear (REG-4, REG-5).
"""

from __future__ import annotations

import copy
import dataclasses
import json
import shutil
from dataclasses import fields, replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ionmc import configio
from ionmc.config import (
    DiagnosticsOptions,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
    validate,
)
from ionmc.configio import (
    ConfigError,
    DatasetError,
    InputManifest,
    from_dict,
    load_toml,
    load_toml_dict,
    load_toml_full,
    manifest_for_api_config,
    to_dict,
    verify_datasets,
)
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.lookup import LookupTable, resample_uniform
from ionmc.materials import MATERIALS, PMMA, WATER, Material
from ionmc.physics.projectiles import ALPHA, PROTON, Projectile
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.sources import PencilBeamSource

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "v3_006"
LOOKUP_FILE = HERE.parent / "data" / "synthetic_lookup.json"
BETHE = "bethe-ionmc-v1"


# ---- builders --------------------------------------------------------------------------------


def _json(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, allow_nan=False)


def _summary(config: SimulationConfig) -> str:
    return _json(validate(config).summary())


def _proton_box(**physics: Any) -> SimulationConfig:
    """An API configuration without arrays (the reconstructible API subset)."""
    return SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 20.0),
        BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 20.0), WATER),
        (ScoringGrid((-10.0, -10.0, 0.0), (2.0, 2.0, 1.0), (10, 10, 20)),),
        PhysicsOptions(nuclear=False, stopping_model=BETHE, **physics),
        RunOptions("python", "float64", 11, 8, n_batches=4),
        tallies=(TallyRequest("let", "dose", "let_d"),),
    )


FULL = SimulationConfig(  # every field of every non-array dataclass at a non-default value
    PencilBeamSource(PROTON, (1.0, 2.0, -5.0), (0.0, 0.6, 0.8), 20.0, 0.1, 1.0),
    BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 20.0), WATER),
    (ScoringGrid((-10.0, -10.0, 0.0), (2.0, 2.0, 1.0), (10, 10, 20), "grid_a"),),
    PhysicsOptions(
        nuclear=True,
        stopping_model="bethe-ionmc-v1+ls-eff",
        straggling=False,
        straggling_model="bohr_gauss_clamped_gamma_v1",
        multiple_scattering=False,
        mcs_model="other_mcs",
        delta_electrons="other_delta",
        e_cut_mev=3.0,
        max_step_mm=0.5,
        max_energy_loss_fraction=0.01,
        range_alpha=0.1,
        range_rho_f_mm=0.2,
        short_step_fraction=0.005,
        truncated_hinge_diagnostic=True,
        e_cut_deuteron_mev=5.0,
        nuclear_table_id="table-id",
        e_cut_ion_mev_per_u=3.0,
        mcs_mapping_eval="root",
        end_of_range_model="point_v1",
    ),
    RunOptions(
        "warp-cpu",
        "float32",
        11,
        16,
        n_batches=4,
        cpu_workers=1,  # REG-5: never above 1 in a test
        max_steps=1000,
        allow_invalid_result=True,
        worker_timeout_s=5.0,
        memory_budget_bytes=2**30,
        chunk_histories=2048,
    ),
    DiagnosticsOptions(True, True, 2),
    tallies=(
        TallyRequest(
            "spectrum",
            "grid_a",
            "fluence_spectrum",
            species=("proton",),
            generation="primary",
            energy_edges_mev_per_u=(1.0, 2.0, 4.0),
            let_medium="medium_a",
            dose_reference="water",
        ),
        TallyRequest("weighted", "grid_a", "lookup_sum", lookup="table_a"),
    ),
)
CUSTOM_PROJECTILE = replace(
    _proton_box(),
    source=PencilBeamSource(
        Projectile("custom-ion", "cx", 1, 1, 938.0), (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 20.0
    ),
)
ALPHA_BOX = replace(
    _proton_box(),
    source=PencilBeamSource(ALPHA, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 400.0),
    geometry=BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), PMMA),
    tallies=(),
)
CUSTOM_MATERIAL_BOX = replace(
    _proton_box(),
    geometry=BoxPhantom(
        (-10.0, -10.0, 0.0),
        (20.0, 20.0, 20.0),
        Material("synthetic_gel", 1.05, {"H": 0.1, "C": 0.2, "O": 0.7}, 70.0, None, "synthetic"),
    ),
)
USER_STOPPING = replace(
    _proton_box(),
    physics=PhysicsOptions(nuclear=False, stopping=BetheStoppingSource()),
)


def _inmemory_lookup_config() -> SimulationConfig:
    axis, values, record = resample_uniform(
        [0.0, 1.0, 3.0, 4.0], {"proton": [1.0, 1.1, 1.3, 1.4]}, n_points=5, spacing="linear"
    )
    table = LookupTable(
        "resampled",
        "synthetic factor",
        "1",
        "let_water_kev_um",
        "linear",
        axis,
        values,
        "synthetic citation",
        "CC0-1.0",
        "built in the test",
        True,
        resampling=record,
    )
    return replace(_proton_box(), lookups=(table,))


def _copy_inputs(tmp_path: Path) -> Path:
    """The synthetic input files in ``tmp_path/inputs`` (the TOML goes to ``tmp_path/cfg``)."""
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for name in ("slab_material_index.npy", "slab_density.npy"):
        shutil.copy(DATA / name, inputs / name)
    shutil.copy(LOOKUP_FILE, inputs / "lookup.json")
    (tmp_path / "cfg").mkdir()
    return tmp_path / "cfg"


def _sha(path: Path) -> str:
    return configio.hashlib.sha256(path.read_bytes()).hexdigest()


def _slab_toml(
    tmp_path: Path,
    *,
    z_exit: bool = True,
    density: bool = True,
    mi_sha: str | None = None,
    lookup_sha: str | None = None,
    with_lookup: bool = True,
) -> Path:
    """A TOML with a voxel slab (``..`` paths into ``tmp_path/inputs``), a lookup tally and every
    optional geometry key; returns the path of the TOML file."""
    cfg_dir = tmp_path / "cfg" if (tmp_path / "cfg").exists() else _copy_inputs(tmp_path)
    inputs = tmp_path / "inputs"
    lines = [
        'schema = "ionmc.config"',
        'schema_version = "1.0"',
        "[source]",
        'projectile = "proton"',
        "position_mm = [0.0, 0.0, -5.0]",
        "direction = [0.0, 0.0, 1.0]",
        "kinetic_energy_mev = 30.0",
        "[geometry]",
        'type = "voxel"',
        "origin_mm = [-6.0, -6.0, 0.0]",
        "spacing_mm = [4.0, 4.0, 4.0]",
        "shape = [3, 3, 8]",
        'materials = ["water", "pmma"]',
        f'material_index = {{ file = "../inputs/slab_material_index.npy", '
        f'sha256 = "{mi_sha or _sha(inputs / "slab_material_index.npy")}" }}',
    ]
    if density:
        lines.append(
            'density = { file = "../inputs/slab_density.npy", '
            f'sha256 = "{_sha(inputs / "slab_density.npy")}" }}'
        )
    if z_exit:
        lines.append("z_exit_mm = 30.0")
    lines += [
        "[[scoring]]",
        'name = "dose"',
        "origin_mm = [-6.0, -6.0, 0.0]",
        "spacing_mm = [4.0, 4.0, 4.0]",
        "shape = [3, 3, 8]",
        "[physics]",
        "nuclear = false",
        'stopping_model = "bethe-ionmc-v1"',
        "[run]",
        'backend = "warp-cpu"',
        'precision = "float64"',
        "seed = 13",
        "n_histories = 8",
        "n_batches = 4",
    ]
    if with_lookup:
        lines += [
            "[[tallies]]",
            'name = "let"',
            'grid = "dose"',
            'quantity = "let_d"',
            "[[tallies]]",
            'name = "weighted"',
            'grid = "dose"',
            'quantity = "lookup_sum"',
            'lookup = "synthetic_let_linear"',
            "[[lookups]]",
            'file = "../inputs/lookup.json"',
            f'sha256 = "{lookup_sha or _sha(inputs / "lookup.json")}"',
        ]
    path = cfg_dir / "slab.toml"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _api_slab_from_files(tmp_path: Path) -> SimulationConfig:
    """The API construction equal to :func:`_slab_toml` (same files, loaded by hand)."""
    inputs = tmp_path / "inputs"
    table = LookupTable.from_file(inputs / "lookup.json")
    return SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 30.0),
        VoxelGeometry(
            (-6.0, -6.0, 0.0),
            (4.0, 4.0, 4.0),
            (3, 3, 8),
            (WATER, PMMA),
            np.load(inputs / "slab_material_index.npy"),
            np.load(inputs / "slab_density.npy"),
            30.0,
        ),
        (ScoringGrid((-6.0, -6.0, 0.0), (4.0, 4.0, 4.0), (3, 3, 8)),),
        PhysicsOptions(nuclear=False, stopping_model=BETHE),
        RunOptions("warp-cpu", "float64", 13, 8, n_batches=4),
        tallies=(
            TallyRequest("let", "dose", "let_d"),
            TallyRequest("weighted", "dose", "lookup_sum", lookup=table.name),
        ),
        lookups=(table,),
    )


def _set(d: dict[str, Any], path: tuple[Any, ...], value: Any) -> dict[str, Any]:
    out = copy.deepcopy(d)
    node: Any = out
    for key in path[:-1]:
        node = node[key]
    node[path[-1]] = value
    return out


def _drop(d: dict[str, Any], path: tuple[Any, ...]) -> dict[str, Any]:
    out = copy.deepcopy(d)
    node: Any = out
    for key in path[:-1]:
        node = node[key]
    del node[path[-1]]
    return out


# ---- CFG-1 -----------------------------------------------------------------------------------

API_ROUND_TRIP = {
    "proton_box": _proton_box(),
    "alpha_pmma_box": ALPHA_BOX,
    "custom_material": CUSTOM_MATERIAL_BOX,
    "custom_projectile": CUSTOM_PROJECTILE,
    "ls_variant_with_options": _proton_box(end_of_range_model="point_v1", max_step_mm=0.25),
}


@pytest.mark.parametrize("name", sorted(API_ROUND_TRIP))
def test_cfg1_api_config_round_trip_without_manifest(name: str) -> None:
    c = API_ROUND_TRIP[name]
    d = to_dict(c)
    _json(d)  # canonical and JSON-serialisable (allow_nan=False)
    assert d["schema"] == "ionmc.config" and d["schema_version"] == "1.0"
    rebuilt = from_dict(json.loads(json.dumps(d)))
    assert to_dict(rebuilt) == d
    if name != "custom_projectile":  # a custom projectile is not a supported species
        assert _summary(rebuilt) == _summary(c)


def test_cfg1_every_field_of_every_config_dataclass_is_generated_and_round_trips(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # the nuclear dataset check needs a cache table; its mapping is tested separately (CFG-3)
    monkeypatch.setattr(configio, "_dataset_ids", lambda *a, **k: ())
    loaded = load_toml_full(_slab_toml(tmp_path))
    generated: list[SimulationConfig] = [
        FULL,
        loaded.config,
        USER_STOPPING,
        _inmemory_lookup_config(),
        CUSTOM_PROJECTILE,
        CUSTOM_MATERIAL_BOX,
        ALPHA_BOX,
    ]
    expected_classes = {
        "SimulationConfig",
        "PencilBeamSource",
        "Projectile",
        "VoxelGeometry",
        "BoxPhantom",
        "Material",
        "SternheimerParameters",
        "ScoringGrid",
        "PhysicsOptions",
        "RunOptions",
        "DiagnosticsOptions",
        "TallyRequest",
        "LookupTable",
    }
    seen: dict[str, dict[str, bool]] = {}

    def walk(obj: Any) -> None:
        if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
            table = seen.setdefault(type(obj).__name__, {})
            for f in fields(obj):
                if not f.init:
                    continue
                value = getattr(obj, f.name)
                touched = f.default is dataclasses.MISSING or _differs(value, f.default)
                table[f.name] = table.get(f.name, False) or touched
                if f.name != "stopping":  # a user object is not part of the schema
                    walk(value)
        elif isinstance(obj, tuple | list):
            for x in obj:
                walk(x)

    for c in generated:
        walk(c)
    assert (
        set(seen) == expected_classes
    )  # a new config dataclass must be added here and to configio
    not_generated = {
        (cls, name) for cls, table in seen.items() for name, touched in table.items() if not touched
    }
    # RunOptions.cpu_workers > 1 is excluded by REG-5 (single process): its value 1 is serialised
    assert not_generated - {("RunOptions", "cpu_workers")} == set()
    # the dataclass field lists the generator saw are the ones configio's schemas cover
    for cls in (PhysicsOptions, RunOptions, DiagnosticsOptions, SimulationConfig):
        assert {f.name for f in fields(cls) if f.init} == set(seen[cls.__name__])

    d = {
        "full": FULL,
        "inmemory_lookup": generated[3],
        "user_stopping": USER_STOPPING,
    }
    full_dict = to_dict(FULL)
    _json(full_dict)
    assert to_dict(from_dict(full_dict)) == full_dict  # every non-array field, exactly
    for name, c in d.items():
        _json(to_dict(c))  # every config serialises, including the not-reconstructible ones
        assert name


def _differs(value: Any, default: Any) -> bool:
    if isinstance(value, np.ndarray) or isinstance(default, np.ndarray):
        return True
    return bool(value != default)


def test_cfg1_toml_loaded_config_round_trips_with_the_manifest(tmp_path: Path) -> None:
    loaded = load_toml_full(_slab_toml(tmp_path))
    c, m = loaded.config, loaded.manifest
    d = to_dict(c, manifest=m)
    _json(d)
    assert d["geometry"]["material_index"]["ref"] == "material_index"
    assert d["geometry"]["density"]["ref"] == "density"
    assert d["lookups"][0]["ref"] == "lookup:synthetic_let_linear"
    assert set(d["geometry"]["material_index"]) == {"ref", "sha256", "content_sha256"}
    rebuilt = from_dict(json.loads(json.dumps(d)), manifest=m)
    assert to_dict(rebuilt, manifest=m) == d  # the exact round trip of the reconstructible subset
    # arrays compared by content digest
    assert configio.array_content_sha256(rebuilt.geometry.material_index) == (
        configio.array_content_sha256(c.geometry.material_index)
    )
    assert configio.array_content_sha256(rebuilt.geometry.density_g_cm3) == (
        configio.array_content_sha256(c.geometry.density_g_cm3)
    )
    assert rebuilt.lookups[0].content_sha256 == c.lookups[0].content_sha256
    assert rebuilt.geometry.z_exit_mm == 30.0
    assert _summary(rebuilt) == _summary(c)
    m.verify(rebuilt)  # the rebuilt configuration is the one the manifest describes


def test_cfg1_api_config_equal_to_the_manifest_files_serialises_to_the_same_references(
    tmp_path: Path,
) -> None:
    loaded = load_toml_full(_slab_toml(tmp_path))
    api = _api_slab_from_files(tmp_path)
    assert to_dict(api, manifest=loaded.manifest) == to_dict(
        loaded.config, manifest=loaded.manifest
    )
    assert _summary(api) == _summary(loaded.config)
    loaded.manifest.verify(api)  # ... and the manifest describes it too


def test_cfg1_summary_equal_after_round_trip_for_a_tabulated_free_set() -> None:
    for c in (_proton_box(), ALPHA_BOX, CUSTOM_MATERIAL_BOX):
        assert _summary(from_dict(to_dict(c))) == _summary(c)


def test_cfg1_not_reconstructible_nuclear_dataset_check_is_reached(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Without the stub the loader's nuclear table check runs (no cache table: DatasetError)."""
    empty = Path(__file__).parent / "does-not-exist-v3-006"
    monkeypatch.setenv("IONMC_CACHE_DIR", str(empty))
    with pytest.raises(DatasetError):
        from_dict(to_dict(FULL))


# ---- CFG-2 -----------------------------------------------------------------------------------


def _api_equivalents() -> dict[str, SimulationConfig]:
    gel = Material(
        "synthetic_gel",
        1.05,
        {"H": 0.10, "C": 0.20, "O": 0.70},
        70.0,
        None,
        "synthetic composition for the configuration tests",
    )
    grid = ScoringGrid((-10.0, -10.0, 0.0), (2.0, 2.0, 1.0), (10, 10, 20))
    proton_box = SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 20.0, 0.1, 1.0),
        BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 20.0), WATER),
        (grid,),
        PhysicsOptions(nuclear=False, stopping_model=BETHE, max_step_mm=0.5),
        RunOptions("python", "float64", 11, 8, n_batches=4),
    )
    warp_gel = SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 25.0),
        BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 20.0), gel),
        (ScoringGrid((-10.0, -10.0, 0.0), (4.0, 4.0, 2.0), (5, 5, 10)),),
        PhysicsOptions(nuclear=False, stopping_model=BETHE, straggling=False),
        RunOptions("warp-cpu", "float32", 12, 8, n_batches=4, chunk_histories=1024),
        tallies=(TallyRequest("let", "dose", "let_d"),),
    )
    inputs = DATA
    table = LookupTable.from_file(LOOKUP_FILE)
    slab = SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 30.0),
        VoxelGeometry(
            (-6.0, -6.0, 0.0),
            (4.0, 4.0, 4.0),
            (3, 3, 8),
            (WATER, PMMA),
            np.load(inputs / "slab_material_index.npy"),
            np.load(inputs / "slab_density.npy"),
        ),
        (ScoringGrid((-6.0, -6.0, 0.0), (4.0, 4.0, 4.0), (3, 3, 8)),),
        PhysicsOptions(nuclear=False, stopping_model=BETHE),
        RunOptions("warp-cpu", "float64", 13, 8, n_batches=4),
        tallies=(
            TallyRequest("let", "dose", "let_d"),
            TallyRequest("weighted", "dose", "lookup_sum", lookup=table.name),
        ),
        lookups=(table,),
    )
    carbon = SimulationConfig(
        PencilBeamSource(
            __import__("ionmc.physics.projectiles", fromlist=["CARBON12"]).CARBON12,
            (0.0, 0.0, -5.0),
            (0.0, 0.0, 2.0),
            1200.0,
        ),
        BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), WATER),
        (ScoringGrid((-10.0, -10.0, 0.0), (5.0, 5.0, 2.0), (4, 4, 20)),),
        PhysicsOptions(nuclear=False, stopping_model=BETHE, multiple_scattering=False),
        RunOptions("warp-cpu", "float64", 14, 8, n_batches=4),
    )
    return {
        "proton_box.toml": proton_box,
        "proton_box_warp_material.toml": warp_gel,
        "voxel_slab_lookup.toml": slab,
        "carbon_box.toml": carbon,
    }


@pytest.mark.parametrize("name", sorted(_api_equivalents()))
def test_cfg2_documented_toml_examples_equal_the_api_construction(name: str) -> None:
    api = _api_equivalents()[name]
    config, manifest = load_toml(DATA / name)
    assert _summary(config) == _summary(api)
    assert manifest.source == "toml"
    manifest.verify(config)
    # the same canonical dict (arrays of the TOML config referenced, those of the API config not)
    assert to_dict(config, manifest=manifest) == to_dict(api, manifest=manifest)


def test_cfg2_example_files_are_documented_and_at_least_three() -> None:
    examples = sorted(DATA.glob("*.toml"))
    assert len(examples) >= 3
    for p in examples:
        assert p.read_text("utf-8").lstrip().startswith("# Documented example")


def test_cfg2_output_table_is_returned_with_absolute_path() -> None:
    loaded = load_toml_full(DATA / "proton_box.toml")
    assert loaded.output.persist_batches is True
    assert loaded.output.path == str((DATA / "out" / "proton_box").resolve())
    assert load_toml_full(DATA / "voxel_slab_lookup.toml").output.persist_batches is False
    assert load_toml_full(DATA / "carbon_box.toml").output == configio.OutputOptions()


def test_cfg2_load_toml_dict_is_the_strictly_validated_document() -> None:
    d = load_toml_dict(DATA / "proton_box.toml")
    assert d["physics"]["stopping_model"] == BETHE
    assert _summary(from_dict(d, base_dir=DATA)) == _summary(load_toml(DATA / "proton_box.toml")[0])


# ---- CFG-3 -----------------------------------------------------------------------------------


def _base_dict() -> dict[str, Any]:
    return to_dict(_proton_box())


def _refuses(d: dict[str, Any], key_path: str, **kwargs: Any) -> ConfigError:
    with pytest.raises(ConfigError) as info:
        from_dict(d, **kwargs)
    assert info.value.key_path == key_path
    assert key_path in str(info.value)
    assert isinstance(info.value, UnsupportedCombinationError)
    return info.value


@pytest.mark.parametrize(
    ("path", "key_path"),
    [
        ((), "bogus"),
        (("source",), "source.bogus"),
        (("geometry",), "geometry.bogus"),
        (("scoring", 0), "scoring[0].bogus"),
        (("physics",), "physics.bogus"),
        (("run",), "run.bogus"),
        (("diagnostics",), "diagnostics.bogus"),
        (("tallies", 0), "tallies[0].bogus"),
        (("materials", 0), "materials[0].bogus"),
    ],
)
def test_cfg3_unknown_key_at_every_level(path: tuple[Any, ...], key_path: str) -> None:
    d = _set(_base_dict(), (*path, "bogus"), 1)
    _refuses(d, key_path)


def test_cfg3_unknown_key_in_nested_inline_tables() -> None:
    d = to_dict(CUSTOM_MATERIAL_BOX)
    assert "bogus" in str(_refuses(_set(d, ("materials", 0, "bogus"), 1), "materials[0].bogus"))
    gel = _set(d, ("materials", 0, "sternheimer"), {"x0": 1.0, "bogus": 2.0})
    _refuses(gel, "materials[0].sternheimer.bogus")


def test_cfg3_misspelt_table_in_toml_text(tmp_path: Path) -> None:
    text = (DATA / "proton_box.toml").read_text("utf-8").replace("[source]", "[sourse]")
    p = tmp_path / "bad.toml"
    p.write_text(text, encoding="utf-8")
    with pytest.raises(ConfigError) as info:
        load_toml(p)
    assert info.value.key_path == "sourse"


@pytest.mark.parametrize(
    ("path", "value", "key_path"),
    [
        (("run", "n_histories"), "8", "run.n_histories"),
        (("run", "n_histories"), 8.0, "run.n_histories"),
        (("run", "seed"), True, "run.seed"),
        (("source", "kinetic_energy_mev"), "20", "source.kinetic_energy_mev"),
        (("source", "position_mm"), [0.0, 1.0], "source.position_mm"),
        (("source", "position_mm"), [0.0, 1.0, "a"], "source.position_mm[2]"),
        (("physics", "nuclear"), "false", "physics.nuclear"),
        (("physics", "stopping_model"), 3, "physics.stopping_model"),
        (("physics", "max_step_mm"), float("nan"), "physics.max_step_mm"),
        (("geometry",), [1, 2], "geometry"),
        (("scoring",), {"a": 1}, "scoring"),
        (("scoring", 0, "shape"), [1, 2, 3.5], "scoring[0].shape[2]"),
        (("tallies", 0, "species"), "proton", "tallies[0].species"),
    ],
)
def test_cfg3_wrong_type(path: tuple[Any, ...], value: Any, key_path: str) -> None:
    _refuses(_set(_base_dict(), path, value), key_path)


def test_cfg3_semantic_errors_name_the_table() -> None:
    err = _refuses(_set(_base_dict(), ("run", "backend"), "gpu"), "run")
    assert "backend" in str(err)
    _refuses(_set(_base_dict(), ("source", "kinetic_energy_mev"), -1.0), "source")


def test_cfg3_missing_nuclear_and_missing_stopping_model() -> None:
    _refuses(_drop(_base_dict(), ("physics", "nuclear")), "physics.nuclear")
    err = _refuses(_drop(_base_dict(), ("physics", "stopping_model")), "physics.stopping_model")
    assert "required" in str(err)
    for table in ("source", "geometry", "scoring", "physics", "run"):
        _refuses(_drop(_base_dict(), (table,)), table)


def test_cfg3_missing_nuclear_in_toml_text(tmp_path: Path) -> None:
    text = (DATA / "proton_box.toml").read_text("utf-8").replace("nuclear = false\n", "")
    p = tmp_path / "x.toml"
    p.write_text(text, encoding="utf-8")
    with pytest.raises(ConfigError) as info:
        load_toml(p)
    assert info.value.key_path == "physics.nuclear"
    text2 = (
        (DATA / "proton_box.toml")
        .read_text("utf-8")
        .replace('stopping_model = "bethe-ionmc-v1"\n', "")
    )
    p.write_text(text2, encoding="utf-8")
    with pytest.raises(ConfigError) as info2:
        load_toml(p)
    assert info2.value.key_path == "physics.stopping_model"


def test_cfg3_schema_header_is_required(tmp_path: Path) -> None:
    base = (DATA / "proton_box.toml").read_text("utf-8")
    p = tmp_path / "x.toml"
    for old, new, key in (
        ('schema = "ionmc.config"\n', "", "schema"),
        ('schema = "ionmc.config"', 'schema = "other"', "schema"),
        ('schema_version = "1.0"\n', "", "schema_version"),
        ('schema_version = "1.0"', 'schema_version = "1.1"', "schema_version"),
        ('schema_version = "1.0"', 'schema_version = "2.0"', "schema_version"),
    ):
        p.write_text(base.replace(old, new), encoding="utf-8")
        with pytest.raises(ConfigError) as info:
            load_toml(p)
        assert info.value.key_path == key


def test_cfg3_toml_syntax_error_and_unreadable_file(tmp_path: Path) -> None:
    p = tmp_path / "x.toml"
    p.write_text("schema = [unterminated\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="not valid TOML"):
        load_toml(p)
    p.write_bytes(b"\xff\xfe")
    with pytest.raises(ConfigError, match="not valid TOML"):
        load_toml(p)
    with pytest.raises(ConfigError, match="cannot read"):
        load_toml(tmp_path / "missing.toml")


def test_cfg3_missing_sha256(tmp_path: Path) -> None:
    p = _slab_toml(tmp_path)
    text = p.read_text("utf-8")
    inputs = tmp_path / "inputs"
    p.write_text(
        text.replace(f', sha256 = "{_sha(inputs / "slab_material_index.npy")}"', ""),
        encoding="utf-8",
    )
    with pytest.raises(ConfigError) as info:
        load_toml(p)
    assert info.value.key_path == "geometry.material_index.sha256"
    p.write_text(text.replace(f'sha256 = "{_sha(inputs / "lookup.json")}"', ""), encoding="utf-8")
    with pytest.raises(ConfigError) as info2:
        load_toml(p)
    assert info2.value.key_path == "lookups[0].sha256"


def test_cfg3_sha256_mismatch(tmp_path: Path) -> None:
    p = _slab_toml(tmp_path, mi_sha="0" * 64)
    with pytest.raises(ConfigError) as info:
        load_toml(p)
    assert info.value.key_path == "geometry.material_index.sha256"
    assert "expected" in str(info.value)
    p = _slab_toml(tmp_path, lookup_sha="1" * 64)
    with pytest.raises(ConfigError) as info2:
        load_toml(p)
    assert info2.value.key_path == "lookups[0].sha256"
    p = _slab_toml(tmp_path, mi_sha="xyz")
    with pytest.raises(ConfigError, match="hexadecimal"):
        load_toml(p)


def test_cfg3_dotdot_paths_resolve_and_are_recorded_absolute_with_sha256(tmp_path: Path) -> None:
    p = _slab_toml(tmp_path)
    config, m = load_toml(p)
    assert m.source == "toml"
    assert m.config_file is not None
    assert m.config_file.path == str(p.resolve()) and m.config_file.sha256 == _sha(p)
    assert m.base_dir == str(p.resolve().parent)
    roles = {f.role: f for f in m.files}
    assert set(roles) == {"material_index", "density", "lookup:synthetic_let_linear"}
    for role, name in (
        ("material_index", "slab_material_index.npy"),
        ("density", "slab_density.npy"),
        ("lookup:synthetic_let_linear", "lookup.json"),
    ):
        target = (tmp_path / "inputs" / name).resolve()
        assert roles[role].path == str(target) and Path(roles[role].path).is_absolute()
        assert ".." not in Path(roles[role].path).parts
        assert roles[role].sha256 == _sha(target)
    j = m.to_dict()
    assert j["source"] == "toml" and j["base_dir"] == m.base_dir
    assert j["config_file"] == {"path": m.config_file.path, "sha256": m.config_file.sha256}
    assert [f["role"] for f in j["files"]] == [f.role for f in m.files]
    assert config.geometry.z_exit_mm == 30.0


def test_cfg3_array_item_forms(tmp_path: Path) -> None:
    d = to_dict(_proton_box())
    voxel = _set(
        d, ("geometry",), {**to_dict(_api_slab_from_files_for_forms(tmp_path))["geometry"]}
    )
    err = _refuses(
        _set(voxel, ("geometry", "material_index"), {"file": "x.npy", "sha256": "0" * 64, "z": 1}),
        "geometry.material_index.z",
    )
    assert "unknown key" in str(err)
    _refuses(_set(voxel, ("geometry", "material_index"), {}), "geometry.material_index")
    _refuses(
        _set(voxel, ("geometry", "material_index"), {"file": "x.npy", "sha256": "0" * 64}),
        "geometry.material_index.file",
    )  # a relative path needs a base_dir


def _api_slab_from_files_for_forms(tmp_path: Path) -> SimulationConfig:
    _copy_inputs(tmp_path)
    return _api_slab_from_files(tmp_path)


# ---- CFG-4 -----------------------------------------------------------------------------------


def test_cfg4_stopping_by_model_id_only() -> None:
    c = from_dict(to_dict(_proton_box()))
    assert c.physics.stopping is None and c.physics.stopping_model == BETHE
    assert c.physics.selection.selected_by == "stopping_model"


def test_cfg4_a_stopping_table_is_rejected(tmp_path: Path) -> None:
    d = _set(_base_dict(), ("physics", "stopping"), {"name": "mine"})
    err = _refuses(d, "physics.stopping")
    assert "stopping_model" in str(err)
    text = (DATA / "proton_box.toml").read_text("utf-8") + "\n[physics.stopping]\nname = 'mine'\n"
    p = tmp_path / "x.toml"
    p.write_text(text, encoding="utf-8")
    with pytest.raises(ConfigError) as info:
        load_toml(p)
    assert info.value.key_path == "physics.stopping"


def test_cfg4_unknown_id_raises_before_physics_options_is_built(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[int] = []
    original = PhysicsOptions.__init__

    def counting(self: PhysicsOptions, *a: Any, **k: Any) -> None:
        calls.append(1)
        original(self, *a, **k)

    monkeypatch.setattr(PhysicsOptions, "__init__", counting)
    from_dict(_base_dict())
    assert calls, "the monkeypatch must observe a valid construction"
    bad_id = _set(_base_dict(), ("physics", "stopping_model"), "no-such-model")
    stop = _set(_base_dict(), ("physics", "stopping"), {"x": 1})
    calls.clear()

    def boom(self: PhysicsOptions, *a: Any, **k: Any) -> None:
        raise AssertionError("PhysicsOptions must not be constructed")

    monkeypatch.setattr(PhysicsOptions, "__init__", boom)
    err = _refuses(bad_id, "physics.stopping_model")
    assert "no-such-model" in str(err) and BETHE in str(err)
    _refuses(stop, "physics.stopping")  # also before construction
    assert not calls


# ---- CFG-5 -----------------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(MATERIALS))
def test_cfg5_builtin_material_names_resolve_to_the_identical_material(name: str) -> None:
    d = _set(_base_dict(), ("geometry", "material"), name)
    d = _set(d, ("materials",), [{"name": name}])
    c = from_dict(d)
    assert isinstance(c.geometry, BoxPhantom)
    assert c.geometry.material is MATERIALS[name]
    # the materials table may be omitted: built-in names resolve without it
    d2 = _drop(d, ("materials",))
    assert from_dict(d2).geometry.material is MATERIALS[name]


def test_cfg5_unknown_builtin_name_and_unknown_reference() -> None:
    d = _set(_base_dict(), ("materials",), [{"name": "unobtainium"}])
    _refuses(d, "materials[0].name")
    _refuses(_set(_base_dict(), ("geometry", "material"), "unobtainium"), "geometry.material")


def test_cfg5_inline_definition_equals_the_api_material() -> None:
    inline = {
        "name": "synthetic_gel",
        "density_g_cm3": 1.05,
        "mass_fractions": {"H": 0.1, "C": 0.2, "O": 0.7},
        "I_eV": 70.0,
        "sternheimer": {"x0": 0.24, "x1": 2.8, "cbar": 3.5, "a": 0.09, "m": 3.4},
        "source": "synthetic",
    }
    d = _set(
        _set(_base_dict(), ("geometry", "material"), "synthetic_gel"), ("materials",), [inline]
    )
    c = from_dict(d)
    assert isinstance(c.geometry, BoxPhantom)
    expected = Material(
        "synthetic_gel",
        1.05,
        {"H": 0.1, "C": 0.2, "O": 0.7},
        70.0,
        configio.SternheimerParameters(0.24, 2.8, 3.5, 0.09, 3.4),
        "synthetic",
    )
    assert c.geometry.material == expected
    assert to_dict(c)["materials"] == [inline]  # a non-built-in material is written in full


@pytest.mark.parametrize(
    "defect",
    [
        {"mass_fractions": {"H": 0.5, "O": 0.4}},  # does not sum to 1
        {"mass_fractions": {"Xx": 1.0}},  # unknown element
        {"mass_fractions": {"H": 1.0}, "density_g_cm3": -1.0},  # negative density
        {"mass_fractions": {"H": 1.0}, "I_eV": 0.0},  # non-positive I
        {"mass_fractions": {"H": -0.2, "O": 1.2}},  # negative fraction
    ],
)
def test_cfg5_inline_defects_give_the_same_error_as_the_api(defect: dict[str, Any]) -> None:
    fields_ = {"name": "bad", "density_g_cm3": 1.0, "mass_fractions": {"H": 1.0}, **defect}
    with pytest.raises(ValueError) as api:
        Material(
            "bad",
            fields_["density_g_cm3"],
            fields_["mass_fractions"],
            fields_.get("I_eV"),
            None,
            "",
        )
    d = _set(_set(_base_dict(), ("geometry", "material"), "bad"), ("materials",), [fields_])
    err = _refuses(d, "materials[0]")
    assert str(api.value) in str(err)


def test_cfg5_duplicate_material_names_and_name_clash_on_output() -> None:
    d = _set(_base_dict(), ("materials",), [{"name": "water"}, {"name": "water"}])
    _refuses(d, "materials[1].name")
    clash = replace(
        _proton_box(),
        geometry=VoxelGeometry(
            (0.0, 0.0, 0.0),
            (1.0, 1.0, 1.0),
            (2, 1, 1),
            (WATER, WATER.with_I(75.0, "water")),
            np.array([[[0]], [[1]]]),
        ),
    )
    with pytest.raises(ConfigError, match="two different materials are named"):
        to_dict(clash)


# ---- CFG-6 -----------------------------------------------------------------------------------


def test_cfg6_no_environment_variable_overrides_a_physics_or_run_field(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    before = to_dict(load_toml(DATA / "proton_box.toml")[0])
    for name, value in {
        "IONMC_BACKEND": "warp-cuda",
        "IONMC_PRECISION": "float32",
        "IONMC_SEED": "999",
        "IONMC_N_HISTORIES": "1000000",
        "IONMC_N_BATCHES": "2",
        "IONMC_CPU_WORKERS": "8",
        "IONMC_NUCLEAR": "1",
        "IONMC_STOPPING_MODEL": "icru90-water-v1",
        "IONMC_MAX_STEP_MM": "9",
        "IONMC_E_CUT_MEV": "9",
        "IONMC_OFFLINE": "0",
        "IONMC_REQUIRE_DATA": "0",
    }.items():
        monkeypatch.setenv(name, value)
    config, manifest = load_toml(DATA / "proton_box.toml")
    assert to_dict(config) == before
    assert config.run.backend == "python" and config.run.seed == 11
    assert config.run.cpu_workers == 1 and config.physics.stopping_model == BETHE
    assert from_dict(before).run.seed == 11
    assert manifest.to_dict()["files"] == []


# ---- CFG-7 -----------------------------------------------------------------------------------


def test_cfg7_user_stopping_object_is_not_reconstructible() -> None:
    d = to_dict(USER_STOPPING)
    rec = d["physics"]["stopping"]
    assert rec["status"] == "not-reconstructible" and rec["reason"]
    assert rec["reconstructible"] is False
    ident = rec["stopping_object"]
    assert ident["model_id"] == BETHE and len(ident["spec_sha256"]) == 64
    assert ident["class"] == "BetheStoppingSource"
    assert "object at 0x" not in _json(d)  # no repr
    assert d["physics"]["stopping_model"] is None
    err = _refuses(d, "physics.stopping")
    assert "not-reconstructible" in str(err)


def test_cfg7_api_created_lookup_table_is_not_reconstructible() -> None:
    c = _inmemory_lookup_config()
    d = to_dict(c)
    rec = d["lookups"][0]
    assert set(rec) == {"status", "reason", "content_sha256"}
    assert rec["status"] == "not-reconstructible"
    assert rec["content_sha256"] == c.lookups[0].content_sha256
    _refuses(d, "lookups[0]")
    # a manifest that does not list it changes nothing; the API manifest never references it
    api = manifest_for_api_config(c)
    assert to_dict(c, manifest=api)["lookups"][0] == rec
    _refuses(to_dict(c, manifest=api), "lookups[0]", manifest=api)


def test_cfg7_arrays_without_a_manifest_or_matching_entry(tmp_path: Path) -> None:
    loaded = load_toml_full(_slab_toml(tmp_path))
    c, m = loaded.config, loaded.manifest
    d = to_dict(c)  # no manifest: no reference
    for ptr in (("geometry", "material_index"), ("geometry", "density")):
        rec = d[ptr[0]][ptr[1]]
        assert set(rec) == {"status", "reason", "content_sha256"}
        assert rec["status"] == "not-reconstructible"
        _refuses(
            _set(d, ptr[:1] + ("material_index",), rec) if ptr[1] == "material_index" else d,
            "geometry.material_index",
        )
    assert d["lookups"][0]["status"] == "not-reconstructible"
    # a manifest with no matching entry (the loader's manifest of another configuration)
    other = load_toml(DATA / "voxel_slab_lookup.toml")[1]
    d2 = to_dict(c, manifest=other)
    assert d2["geometry"]["material_index"]["status"] == "not-reconstructible"
    assert (
        d["geometry"]["material_index"]["content_sha256"]
        == (to_dict(c, manifest=m)["geometry"]["material_index"]["content_sha256"])
    )


def test_cfg7_reference_failures_name_the_pointer(tmp_path: Path) -> None:
    loaded = load_toml_full(_slab_toml(tmp_path))
    c, m = loaded.config, loaded.manifest
    d = to_dict(c, manifest=m)
    # no manifest at all
    err = _refuses(d, "geometry.material_index")
    assert "no manifest" in str(err)
    # a manifest lacking the role
    no_density = replace(m, files=tuple(f for f in m.files if f.role != "density"))
    err = _refuses(d, "geometry.density", manifest=no_density)
    assert "role 'density'" in str(err)
    api = manifest_for_api_config(c)  # entries without files
    err = _refuses(d, "geometry.material_index", manifest=api)
    assert "unavailable" in str(err)
    # a wrong digest in the reference
    _refuses(
        _set(d, ("geometry", "material_index", "content_sha256"), "2" * 64),
        "geometry.material_index.content_sha256",
        manifest=m,
    )
    _refuses(
        _set(d, ("geometry", "density", "sha256"), "3" * 64), "geometry.density.sha256", manifest=m
    )
    _refuses(
        _set(d, ("geometry", "material_index", "ref"), "density"),
        "geometry.material_index.ref",
        manifest=m,
    )
    _refuses(_set(d, ("lookups", 0, "ref"), "lookup:other"), "lookups[0]", manifest=m)
    # a file whose bytes changed after the manifest was made
    np.save(tmp_path / "inputs" / "slab_material_index.npy", np.ones((3, 3, 8), dtype=np.int32))
    err = _refuses(d, "geometry.material_index.sha256", manifest=m)
    assert "expected" in str(err)


def test_cfg7_changed_content_digest_with_matching_file_hash(tmp_path: Path) -> None:
    """The file hash is right but the manifest entry records another content digest."""
    loaded = load_toml_full(_slab_toml(tmp_path))
    m = loaded.manifest
    bad = replace(
        m,
        files=tuple(
            replace(f, content_sha256="4" * 64) if f.role == "material_index" else f
            for f in m.files
        ),
    )
    d = to_dict(loaded.config, manifest=m)
    _refuses(d, "geometry.material_index.content_sha256", manifest=bad)


@pytest.mark.parametrize(
    "change",
    [
        {"citation": "another citation"},
        {"license": "MIT"},
        {"source": "somewhere else"},
        {"resampling": {"method": "other"}},
    ],
)
def test_cfg7_lookup_with_other_provenance_emits_no_reference(
    tmp_path: Path, change: dict[str, Any]
) -> None:
    loaded = load_toml_full(_slab_toml(tmp_path))
    c, m = loaded.config, loaded.manifest
    table = c.lookups[0]
    twin = replace(table, **change)
    assert twin.content_sha256 == table.content_sha256 or "resampling" in change
    assert twin.content_sha256 == table.content_sha256  # numerically (and content-hash) equal
    c2 = replace(c, lookups=(twin,))
    rec = to_dict(c2, manifest=m)["lookups"][0]
    assert rec["status"] == "not-reconstructible" and rec["content_sha256"] == table.content_sha256
    assert to_dict(c, manifest=m)["lookups"][0]["ref"] == "lookup:synthetic_let_linear"
    # PROV-8: different provenance gives a different config_sha256 and the manifest refuses it
    assert configio.config_sha256(c2, m) != m.config_sha256
    with pytest.raises(ConfigError):
        m.verify(c2)
    _refuses(to_dict(c2, manifest=m), "lookups[0]", manifest=m)


def test_cfg7_no_default_is_substituted_for_a_missing_array(tmp_path: Path) -> None:
    d = to_dict(_api_slab_from_files_for_forms(tmp_path))
    _refuses(_drop(d, ("geometry", "material_index")), "geometry.material_index")
    _refuses(_set(d, ("geometry", "material_index"), None), "geometry.material_index")


# ---- PROV-8 (manifest part) ------------------------------------------------------------------


def test_prov8_toml_manifest_records_files_and_binds_the_configuration(tmp_path: Path) -> None:
    p = _slab_toml(tmp_path)
    c, m = load_toml(p)
    assert isinstance(m, InputManifest) and m.source == "toml"
    j = m.to_dict()
    _json(j)
    assert set(j) == {"source", "base_dir", "config_file", "files", "config_sha256"}
    for f in j["files"]:
        assert set(f) == {"role", "path", "sha256", "content_sha256", "item_sha256"}
        assert all(isinstance(f[k], str) for k in f)
    assert m.config_sha256 == configio.config_sha256(c, m)
    assert m.config_sha256 == configio.sha256_canonical(to_dict(c, manifest=m))
    m.verify(c)
    # determinism: the loader, verify and a second computation agree
    assert load_toml(p)[1].config_sha256 == m.config_sha256
    assert configio.config_sha256(c, m) == configio.config_sha256(c, m)


def test_prov8_api_manifest_has_unavailable_file_pointers(tmp_path: Path) -> None:
    c = _api_slab_from_files_for_forms(tmp_path)
    m = manifest_for_api_config(c)
    assert m.source == "api" and m.base_dir is None and m.config_file is None
    j = m.to_dict()
    unavailable = {"status": "unavailable", "reason": "constructed in process"}
    assert j["base_dir"] == unavailable and j["config_file"] == unavailable
    assert [f["role"] for f in j["files"]] == [
        "material_index",
        "density",
        "lookup:synthetic_let_linear",
    ]
    assert all(f["path"] == unavailable and f["sha256"] == unavailable for f in j["files"])
    m.verify(c)
    # an API config without arrays has no files at all
    plain = manifest_for_api_config(_proton_box())
    assert plain.files == ()
    plain.verify(_proton_box())
    assert to_dict(from_dict(to_dict(_proton_box()), manifest=plain), manifest=plain) == to_dict(
        _proton_box()
    )


def test_prov8_editing_a_referenced_file_after_load_is_refused(tmp_path: Path) -> None:
    c, m = load_toml(_slab_toml(tmp_path))
    m.verify(c)
    target = tmp_path / "inputs" / "slab_density.npy"
    np.save(target, np.full((3, 3, 8), 2.0))
    with pytest.raises(ConfigError, match="changed since it was loaded"):
        m.verify(c)
    target.unlink()
    with pytest.raises(ConfigError, match="cannot read"):
        m.verify(c)


def test_prov8_editing_the_toml_after_load_is_refused(tmp_path: Path) -> None:
    p = _slab_toml(tmp_path)
    c, m = load_toml(p)
    m.verify(c)
    p.write_text(p.read_text("utf-8") + "# edited after load\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="changed since it was loaded") as info:
        m.verify(c)
    assert info.value.key_path == "input_manifest.config_file"


def test_prov8_content_digest_that_disagrees_with_the_config_arrays_is_refused(
    tmp_path: Path,
) -> None:
    c, m = load_toml(_slab_toml(tmp_path))
    bad = replace(
        m,
        files=tuple(
            replace(f, content_sha256="5" * 64) if f.role == "density" else f for f in m.files
        ),
    )
    with pytest.raises(ConfigError, match="content digest"):
        bad.verify(c)
    # other arrays than those the manifest recorded
    other_geometry = replace(
        c.geometry,  # type: ignore[type-var]
        density_g_cm3=np.full((3, 3, 8), 1.5),
    )
    with pytest.raises(ConfigError, match="content digest"):
        m.verify(replace(c, geometry=other_geometry))
    with pytest.raises(ConfigError, match="no array of this role"):
        m.verify(replace(c, geometry=BoxPhantom((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), WATER)))


def test_prov8_manifest_of_another_configuration_is_refused_through_config_sha256(
    tmp_path: Path,
) -> None:
    c, m = load_toml(_slab_toml(tmp_path))
    # (i) config B has the same arrays but another scalar
    b = replace(c, run=replace(c.run, seed=99))
    assert to_dict(b, manifest=m)["geometry"] == to_dict(c, manifest=m)["geometry"]
    with pytest.raises(ConfigError, match="does not describe this configuration") as info:
        m.verify(b)
    assert info.value.key_path == "input_manifest.config_sha256"
    # (ii) a configuration mutated in one scalar after loading
    mutated = replace(c, source=replace(c.source, kinetic_energy_mev=31.0))
    with pytest.raises(ConfigError, match="does not describe this configuration"):
        m.verify(mutated)
    # a tally, a physics switch and a scoring-grid change are refused as well
    for changed in (
        replace(c, physics=replace(c.physics, max_step_mm=0.75)),
        replace(c, tallies=c.tallies[:1], lookups=()),
        replace(c, scoring=(replace(c.scoring[0], name="other"),)),
    ):
        with pytest.raises(ConfigError):
            m.verify(changed)
    m.verify(c)


def test_prov8_manifest_is_immutable_and_validates_itself() -> None:
    m = manifest_for_api_config(_proton_box())
    with pytest.raises(dataclasses.FrozenInstanceError):
        m.source = "toml"  # type: ignore[misc]
    with pytest.raises(ConfigError):
        InputManifest("ftp", None, None, (), "")
    with pytest.raises(ConfigError):
        InputManifest("toml", None, None, (), "")
    f = configio.InputFile("density", None, None, "a", "b")
    with pytest.raises(ConfigError, match="duplicate roles"):
        InputManifest("api", None, None, (f, f), "")


# ---- dataset checks (decision 0045 section 7) ------------------------------------------------


def test_datasets_missing_cache_object_is_a_dataset_error_not_a_config_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("IONMC_CACHE_DIR", str(tmp_path / "empty-cache"))
    d = _set(_base_dict(), ("physics", "stopping_model"), "nist-pstar-astar-water-2005")
    with pytest.raises(DatasetError) as info:
        from_dict(d)
    assert not isinstance(info.value, UnsupportedCombinationError)
    assert isinstance(info.value.__cause__, FileNotFoundError)
    # an unknown id is still a ConfigError, whatever the cache holds
    with pytest.raises(ConfigError):
        from_dict(_set(_base_dict(), ("physics", "stopping_model"), "nope"))
    # a nuclear run needs its table
    nuclear = _set(
        _set(_base_dict(), ("physics", "nuclear"), True), ("physics", "nuclear_table_id"), "0" * 64
    )
    with pytest.raises(DatasetError, match="nuclear table"):
        from_dict(nuclear)


def test_datasets_verify_datasets_of_an_analytic_config_needs_nothing(tmp_path: Path) -> None:
    assert verify_datasets(_proton_box(), tmp_path) == ()
    assert verify_datasets(USER_STOPPING, tmp_path) == ()


def test_datasets_verify_datasets_names_the_pinned_objects() -> None:
    """With the cache present the NIST model's pins verify; without it ``DatasetError``."""
    from ionmc.data import cache

    cdir = cache.resolve_cache_dir(None)
    ids = ("nist-pstar-water-2005", "nist-astar-water-2005")
    present = all((cdir / "manifests" / f"{i}.json").is_file() for i in ids)
    if not present:
        pytest.skip("the NIST datasets are not in the cache (IONMC_CACHE_DIR)")
    c = from_dict(_set(_base_dict(), ("physics", "stopping_model"), "nist-pstar-astar-water-2005"))
    assert set(verify_datasets(c)) == set(ids)
    with pytest.raises(DatasetError) as info:
        verify_datasets(c, cdir / "nowhere")
    assert info.value.__cause__ is not None


def test_datasets_loader_reads_no_physics_or_run_field_from_the_environment() -> None:
    src = (Path(configio.__file__)).read_text("utf-8")
    assert "os.environ" not in src and "getenv" not in src
