"""Rows OUT-0, OUT-1, OUT-2 (CI part), OUT-6, OUT-8 (bitwise part), PROV-8 (bundle part), REG-4 and
REG-5 of plan v3-006: real tiny transport runs written with ``ionmc.results.save`` and read back.

Every run here is at most a few dozen histories with a fixed small seed (never the ``out6`` seeds),
single process (``cpu_workers = 1``), python or warp-cpu float64. Rows whose generators belong to
later steps (full report, plots, ``ionmc verify`` CLI, three-backend LV/HR runs) are not here.
"""

from __future__ import annotations

import dataclasses
import json
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ionmc import results
from ionmc.config import PhysicsOptions, RunOptions, SimulationConfig
from ionmc.configio import (
    UNAVAILABLE,
    ConfigError,
    InputManifest,
    array_content_sha256,
    config_sha256,
    load_toml,
    manifest_for_api_config,
    sha256_canonical,
    to_dict,
)
from ionmc.errors import TransportLimitError
from ionmc.geometry import BoxPhantom
from ionmc.materials import WATER
from ionmc.physics.projectiles import ALPHA, PROTON
from ionmc.results import LoadedResult, load, save
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.simulation import Result, Simulation
from ionmc.sources import PencilBeamSource

DATA = Path(__file__).resolve().parent / "data" / "v3_006"
SRC = Path(__file__).resolve().parents[2] / "src"
SLAB_TOML = DATA / "voxel_slab_lookup.toml"
BOX_TOML = DATA / "proton_box.toml"

_RAN: list[SimulationConfig] = []  # every configuration run by this module (REG-5)


@dataclass
class Case:
    cfg: SimulationConfig
    manifest: InputManifest | None
    result: Result
    bundle: Path
    loaded: LoadedResult


def _run(cfg: SimulationConfig) -> Result:
    assert cfg.run.cpu_workers == 1
    _RAN.append(cfg)
    return Simulation(cfg).run()


def _case(root: Path, cfg: SimulationConfig, manifest: InputManifest | None) -> Case:
    r = _run(cfg)
    b = save(r, root / "bundle", manifest=manifest)
    return Case(cfg, manifest, r, b, load(b))


def _box_config(**run_kw: Any) -> SimulationConfig:
    run = dict(n_batches=3, seed=5) | run_kw
    return SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 20.0),
        BoxPhantom((-3.0, -3.0, 0.0), (6.0, 6.0, 6.0), WATER),
        # non-representable spacings: the edges are the interesting bitwise object of OUT-6
        (ScoringGrid((-1.2, -1.2, 0.1), (0.3, 0.3, 0.7), (8, 8, 6)),),
        PhysicsOptions(nuclear=False, stopping_model="bethe-ionmc-v1"),
        RunOptions(
            run.pop("backend", "python"), "float64", run.pop("seed"), run.pop("n", 12), **run
        ),
    )


@pytest.fixture(scope="module")
def case_a(tmp_path_factory: pytest.TempPathFactory) -> Case:
    """OUT-1 (a): proton box, python."""
    return _case(tmp_path_factory.mktemp("a"), _box_config(), None)


@pytest.fixture(scope="module")
def case_b(tmp_path_factory: pytest.TempPathFactory) -> Case:
    """OUT-1 (b): voxel slab, warp-cpu float64, LET + spectrum + synthetic lookup tallies; API."""
    cfg, _ = load_toml(SLAB_TOML)
    spectrum = TallyRequest(
        "spectrum", "dose", "fluence_spectrum", energy_edges_mev_per_u=(1.0, 11.0, 21.0, 31.0)
    )
    cfg = dataclasses.replace(cfg, tallies=(*cfg.tallies, spectrum))
    return _case(tmp_path_factory.mktemp("b"), cfg, None)


@pytest.fixture(scope="module")
def case_c(tmp_path_factory: pytest.TempPathFactory) -> Case:
    """OUT-1 (c): a real run with an exhausted step bound, allowed to return its invalid result."""
    cfg = _box_config(max_steps=2, allow_invalid_result=True, n=6)
    return _case(tmp_path_factory.mktemp("c"), cfg, None)


@pytest.fixture(scope="module")
def case_d(tmp_path_factory: pytest.TempPathFactory) -> Case:
    """OUT-1 (d): the TOML slab with lookup files and the loader's manifest."""
    cfg, man = load_toml(SLAB_TOML)
    return _case(tmp_path_factory.mktemp("d"), cfg, man)


@pytest.fixture(scope="module")
def case_e(tmp_path_factory: pytest.TempPathFactory) -> Case:
    """OUT-1 (e): the configuration of a TOML file built through the API, no manifest."""
    cfg, _ = load_toml(BOX_TOML)
    return _case(tmp_path_factory.mktemp("e"), cfg, None)


def _doc_json(b: Path) -> dict[str, Any]:
    return json.loads((b / "result.json").read_text(encoding="utf-8"))  # type: ignore[no-any-return]


def _assert_bitwise(c: Case) -> None:
    """Catalogue arrays equal the writer's inputs bit for bit (dtype, shape, bytes: NaN positions
    and payloads included) and the per-grid base arrays equal the ``GridResult`` fields."""
    persisted = c.loaded.arrays
    exp = {a.name: a.data for a in results._collect_arrays(c.result, True)}
    batches = bool(_doc_json(c.bundle)["statistics"]["batch_arrays"]["status"] == "persisted")
    assert batches
    assert set(persisted) == set(exp) == {e["name"] for e in c.loaded.doc["arrays"]}
    for n, a in persisted.items():
        e = np.ascontiguousarray(exp[n])
        assert (a.dtype, a.shape) == (e.dtype, e.shape), n
        assert a.tobytes() == e.tobytes(), n
    for g in c.result.grids:
        for base, field in (
            ("energy_mev", "energy_mev"),
            ("dose_gy", "dose_gy"),
            ("relative_uncertainty", "relative_uncertainty"),
            ("defined_mask", "defined_mask"),
            ("mass_g", "mass_g"),
        ):
            got, want = c.loaded.grid_array(g.name, base), np.asarray(getattr(g, field))
            assert (got.dtype, got.shape) == (want.dtype, want.shape)
            assert got.tobytes() == want.tobytes(), (g.name, base)
            assert np.array_equal(np.isnan(got), np.isnan(want)) if got.dtype.kind == "f" else True


def _assert_json(c: Case, source: str) -> None:
    doc, r = c.loaded.doc, c.result
    man = c.manifest or manifest_for_api_config(c.cfg)
    assert doc["requested"] == to_dict(c.cfg, manifest=man)
    assert doc["effective"] == json.loads(results.canonical_bytes(r.effective_config.summary()))
    assert doc["input_manifest"] == man.to_dict() and doc["input_manifest"]["source"] == source
    assert doc["valid"] is r.valid and doc["accounting"]["valid"] is r.valid
    run = doc["run"]
    assert (run["backend"], run["device"], run["precision"], run["seed"]) == (
        r.backend, r.device, r.precision, r.seed,
    )  # fmt: skip
    assert (run["n_histories"], run["n_batches"], run["rng"]) == (r.n_histories, r.n_batches, r.rng)
    assert [g["valid"] for g in doc["statistics"]["grids"]] == [g.valid for g in r.grids]
    assert doc["accounting"]["counters"] == json.loads(json.dumps(r.counters.as_dict()))
    assert c.loaded.doc == _doc_json(c.bundle)  # the reader returns the file's document unchanged


# ---- OUT-0 ------------------------------------------------------------------------------------

#: V2-OUTPUT clause -> JSON pointer ("*" = every list element) or ``array:<grid-array base name>``.
CLAUSES: tuple[tuple[str, str], ...] = (
    ("code version / SHA / dirty state", "/provenance/code"),
    ("requested configuration", "/requested"),
    ("effective configuration", "/effective"),
    ("requested physics", "/requested/physics"),
    ("effective physics", "/effective/physics"),
    ("species and source: requested", "/requested/source/projectile"),
    ("species and source: effective", "/effective/source"),
    ("geometry definition", "/requested/geometry"),
    ("material definitions", "/effective_extra/materials"),
    ("transport grid / geometry extents", "/effective/geometry"),
    ("scoring grids with coordinates", "/effective/scoring"),
    ("grid origin", "/effective/scoring/*/origin_mm"),
    ("grid spacing", "/effective/scoring/*/spacing_mm"),
    ("grid convention", "/effective_extra/grid_convention"),
    ("units table", "/units"),
    ("per-array units", "/arrays/*/units"),
    ("external data identifiers and hashes", "/provenance/data"),
    ("table identity and sha256", "/provenance/data/table_sha256"),
    ("model and reference versions", "/provenance/models"),
    ("stopping model version", "/provenance/models/stopping_version"),
    ("actual histories", "/run/n_histories"),
    ("batches", "/run/n_batches"),
    ("histories per batch", "/run/histories_per_batch"),
    ("seed", "/run/seed"),
    ("RNG identity", "/run/rng"),
    ("RNG generator", "/run/rng/generator"),
    ("normalisation", "/statistics/normalisation"),
    ("uncertainty definition", "/statistics/uncertainty"),
    ("validity mask definition", "/statistics/defined_mask"),
    ("validity mask array", "array:defined_mask"),
    ("undefined-value policy (relative uncertainty NaN + mask)", "array:relative_uncertainty"),
    ("escape: escaped energy", "/accounting/energy_balance/escaped_mev"),
    ("truncation: truncated energy", "/accounting/energy_balance/truncated_mev"),
    ("energy closure", "/accounting/energy_balance/closure_residual_mev"),
    ("transport counters", "/accounting/counters"),
    ("Result.valid", "/valid"),
    ("Result.valid (accounting)", "/accounting/valid"),
    ("invalidity reasons", "/accounting/reasons"),
    ("GridResult.valid (every grid)", "/statistics/grids/*/valid"),
    ("low_batch_count (every grid)", "/statistics/grids/*/low_batch_count"),
    ("backend", "/run/backend"),
    ("device", "/run/device"),
    ("precision", "/run/precision"),
    ("relevant environment", "/provenance/environment"),
    ("timings (labelled non-deterministic)", "/run/timings"),
    ("input manifest", "/input_manifest"),
    ("machine-readable data: catalogue", "/arrays"),
    ("machine-readable data: files", "/files/arrays.npz"),
    ("human-readable summary", "/files/report.md"),
    ("dose array", "array:dose_gy"),
)


def _resolve(node: Any, parts: list[str]) -> list[Any]:
    if not parts:
        return [node]
    head, rest = parts[0], parts[1:]
    if head == "*":
        assert isinstance(node, list) and node, "wildcard over an empty or non-list node"
        return [v for item in node for v in _resolve(item, rest)]
    assert isinstance(node, dict) and head in node, f"missing {head!r}"
    return _resolve(node[head], rest)


def _explicit(v: Any) -> bool:
    """Non-empty, or carrying an explicit status (decision 0040: undefined is explicit)."""
    if isinstance(v, dict):
        return bool(v) and ("status" in v or any(_explicit(x) for x in v.values()))
    if isinstance(v, list | str):
        return len(v) > 0
    return v is not None  # numbers and booleans (0 and False are real values)


def test_out0_clause_map(case_d: Case, case_b: Case, case_c: Case) -> None:
    assert len({c for c, _ in CLAUSES}) == len(CLAUSES)
    for case in (case_d, case_b, case_c):
        doc, loaded = case.loaded.doc, case.loaded
        for clause, ptr in CLAUSES:
            if ptr.startswith("array:"):
                for g in case.result.grids:
                    a = loaded.grid_array(g.name, ptr[6:])
                    assert a.size > 0, clause
                continue
            assert ptr.startswith("/"), clause
            vals = _resolve(doc, ptr[1:].split("/"))
            if ptr == "/accounting/reasons":  # explicit: empty exactly for a valid result
                assert (vals[0] == []) is case.result.valid, clause
                continue
            assert vals and all(_explicit(v) for v in vals), f"{clause}: {ptr} is empty"
        # the clauses whose values are the run's own
        r = case.result
        assert doc["valid"] is r.valid and doc["accounting"]["valid"] is r.valid
        assert [g["valid"] for g in doc["statistics"]["grids"]] == [g.valid for g in r.grids]
        assert (doc["run"]["backend"], doc["run"]["device"]) == (r.backend, r.device)
        assert doc["run"]["rng"] == r.rng and doc["run"]["rng"]["generator"] == "philox4x32-10"
        assert (case.bundle / "report.md").stat().st_size > 0
        assert doc["provenance"]["code"].get("status") or doc["provenance"]["code"]["git_commit"]
    assert case_d.result.backend == "warp-cpu" and case_d.result.device == "cpu"
    assert case_c.result.valid is False and case_d.result.valid is True


# ---- OUT-1 ------------------------------------------------------------------------------------


def test_out1a_proton_box_python(case_a: Case) -> None:
    assert case_a.result.backend == "python" and case_a.result.valid
    _assert_bitwise(case_a)
    _assert_json(case_a, "api")
    nan = np.isnan(case_a.loaded.grid_array("dose", "relative_uncertainty"))
    assert nan.any() and not case_a.loaded.grid_array("dose", "defined_mask")[nan].any()


def test_out1b_voxel_slab_let_spectrum_lookup(case_b: Case) -> None:
    assert (case_b.result.backend, case_b.result.precision) == ("warp-cpu", "float64")
    kinds = {e["origin"].get("quantity") for e in case_b.loaded.doc["arrays"]}
    assert {"let_d", "lookup_sum", "fluence_spectrum"} <= kinds
    _assert_bitwise(case_b)
    _assert_json(case_b, "api")
    files = case_b.loaded.doc["input_manifest"]["files"]
    assert {f["role"] for f in files} == {
        "material_index",
        "density",
        "lookup:synthetic_let_linear",
    }
    for f in files:
        assert f["path"] == f["sha256"] == f["bytes"] == UNAVAILABLE


def test_out1c_invalid_result_with_allow_invalid_result(case_c: Case, tmp_path: Path) -> None:
    r = case_c.result
    assert r.valid is False and r.counters.step_truncation > 0 and not r.grids[0].valid
    with pytest.raises(TransportLimitError):
        _run(_box_config(max_steps=2, n=6))
    _assert_bitwise(case_c)
    _assert_json(case_c, "api")
    doc = case_c.loaded.doc
    assert doc["valid"] is False and doc["accounting"]["reasons"] == ["step_truncation"]
    assert doc["statistics"]["grids"][0]["valid"] is False
    assert "INVALID RESULT" in (case_c.bundle / "report.md").read_text(encoding="utf-8")

    # decision 0045 section 1: an invalid result is written only when the run allowed it
    cfg = case_c.cfg
    strict = dataclasses.replace(cfg, run=dataclasses.replace(cfg.run, allow_invalid_result=False))
    with pytest.raises(results.ResultFormatError, match="allow_invalid_result"):
        save(dataclasses.replace(r, requested_config=strict), tmp_path / "strict")
    assert not (tmp_path / "strict").exists()


def test_out1d_toml_run_with_manifest(case_d: Case) -> None:
    _assert_bitwise(case_d)
    _assert_json(case_d, "toml")
    m = case_d.loaded.doc["input_manifest"]
    assert m["config_file"]["path"] == str(SLAB_TOML.resolve())
    assert m["base_dir"] == str(DATA.resolve())
    assert {f["role"] for f in m["files"]} == {
        "material_index", "density", "lookup:synthetic_let_linear",
    }  # fmt: skip
    for f in m["files"]:
        assert re.fullmatch(r"[0-9a-f]{64}", f["sha256"]) and Path(f["path"]).is_file()
        assert f["bytes"] == Path(f["path"]).stat().st_size


def test_out1e_same_config_through_the_api_without_manifest(case_e: Case) -> None:
    _assert_bitwise(case_e)
    _assert_json(case_e, "api")
    m = case_e.loaded.doc["input_manifest"]
    assert m["base_dir"] == m["config_file"] == UNAVAILABLE
    assert m["files"] == []  # a box has no array files


# ---- OUT-2 (CI part) --------------------------------------------------------------------------

_CHILD = r"""
import hashlib, json, sys, zipfile
sys.path.insert(0, sys.argv[1])
import numpy as np
from ionmc import results
from ionmc.configio import array_content_sha256
loaded = results.load(sys.argv[2])
doc = loaded.doc
for name, rec in doc["files"].items():  # every file digest, recomputed
    assert hashlib.sha256(open(sys.argv[2] + "/" + name, "rb").read()).hexdigest() == rec["sha256"]
with zipfile.ZipFile(sys.argv[2] + "/arrays.npz") as zf:  # every content digest, recomputed
    for e in doc["arrays"]:
        with zf.open(e["name"] + ".npy") as fh:
            a = np.lib.format.read_array(fh, allow_pickle=False)
        assert array_content_sha256(a) == e["content_sha256"], e["name"]
print("DUMP:" + json.dumps(results.canonical_bytes(loaded.dump()).decode("utf-8")))
"""


def test_out2_clean_process_round_trip(case_d: Case, case_b: Case) -> None:
    env = {k: v for k, v in os.environ.items() if not k.startswith("PYTHON")}
    env.update(IONMC_SINGLE_PROCESS="1", OMP_NUM_THREADS="1")
    for case in (case_d, case_b):
        p = subprocess.run(
            [sys.executable, "-I", "-c", _CHILD, str(SRC), str(case.bundle)],
            capture_output=True, text=True, timeout=300, env=env, check=False, cwd=case.bundle,
        )  # fmt: skip
        assert p.returncode == 0, p.stderr[-2000:]
        lines = [ln for ln in p.stdout.splitlines() if ln.startswith("DUMP:")]
        assert len(lines) == 1
        writer = {  # the writer's own document and array digests, from the in-memory result
            "doc": _doc_json(case.bundle),
            "content_sha256": {
                a.name: array_content_sha256(a.data)
                for a in sorted(results._collect_arrays(case.result, True), key=lambda a: a.name)
            },
        }
        assert json.loads(lines[0][5:]) == results.canonical_bytes(writer).decode("utf-8")
    assert case_d.bundle != case_b.bundle


# ---- OUT-6 ------------------------------------------------------------------------------------


def test_out6_coordinates(case_a: Case, case_d: Case) -> None:
    for case in (case_a, case_d):
        doc = case.loaded.doc
        assert doc["effective_extra"]["grid_convention"] == results.GRID_CONVENTION
        assert "corner of voxel (0,0,0)" in results.GRID_CONVENTION and "C order" in (
            results.GRID_CONVENTION
        )
        for rec, grid in zip(doc["effective"]["scoring"], case.cfg.scoring, strict=True):
            for ax in range(3):
                # edges recomputed from the persisted origin, spacing and shape alone
                edges = rec["origin_mm"][ax] + rec["spacing_mm"][ax] * np.arange(
                    rec["shape"][ax] + 1, dtype=np.float64
                )
                want = grid.edges_mm(ax)
                assert edges.dtype == want.dtype and edges.tobytes() == want.tobytes()
    # the input arrays of the voxel geometry come back bit for bit
    geom = case_d.cfg.geometry
    for name, want in (
        ("material_index", geom.material_index),
        ("density_g_cm3", geom.density_g_cm3),
    ):
        got = case_d.loaded.arrays["geometry__" + name]
        want = np.asarray(want)
        assert (got.dtype, got.shape) == (want.dtype, want.shape)
        assert got.tobytes() == np.ascontiguousarray(want).tobytes()
    assert case_d.loaded.doc["effective_extra"]["geometry_arrays"] == {
        "material_index": "geometry__material_index", "density_g_cm3": "geometry__density_g_cm3",
    }  # fmt: skip


# ---- OUT-8 (bitwise part) ---------------------------------------------------------------------


def _table_arrays(obj: Any, seen: set[int] | None = None, depth: int = 0) -> Iterator[Any]:
    """Every non-constant numpy array of at least 8 elements reachable from ``obj``."""
    seen = set() if seen is None else seen
    if id(obj) in seen or depth > 8:
        return
    seen.add(id(obj))
    if isinstance(obj, np.ndarray):
        if obj.size >= 8 and obj.dtype.kind in "fiu" and np.ptp(obj) != 0:
            yield obj
    elif isinstance(obj, dict):
        for v in obj.values():
            yield from _table_arrays(v, seen, depth + 1)
    elif isinstance(obj, list | tuple | frozenset | set):
        for v in obj:
            yield from _table_arrays(v, seen, depth + 1)
    elif hasattr(obj, "__dict__") or hasattr(obj, "__slots__"):
        if getattr(type(obj), "__module__", "").startswith("ionmc"):
            fields = vars(obj) if hasattr(obj, "__dict__") else {}
            for v in fields.values():
                yield from _table_arrays(v, seen, depth + 1)


def _assert_no_table_values(case: Case, minimum: int) -> None:
    tables = list(_table_arrays(case.result.effective_config.tables)) + list(
        _table_arrays(case.result.effective_config.nuclear)
    )
    assert len(tables) >= minimum, "the table collector found too few arrays (vacuous check)"
    table_digests = {array_content_sha256(a) for a in tables}
    with np.load(case.bundle / "arrays.npz", allow_pickle=False) as z:
        members = {n: z[n] for n in z.files}
    assert members
    for n, a in members.items():
        assert array_content_sha256(a) not in table_digests, n  # dtype, shape and bytes equal
    for a in members.values():  # and no float member holds a table row verbatim
        for t in tables:
            if a.dtype == t.dtype and a.shape == t.shape:
                assert a.tobytes() != t.tobytes()


def test_out8_no_table_value_equals_a_persisted_member(case_d: Case, case_a: Case) -> None:
    _assert_no_table_values(case_d, 6)
    _assert_no_table_values(case_a, 6)
    names = {e["name"] for e in case_d.loaded.doc["arrays"]}
    for banned in ("range_mm", "stopping", "let_water", "ln_s_mass", "nuclear_table"):
        assert not any(banned in n for n in names)


def test_out8_cache_backed_icru90_alpha_and_nuclear_tables(tmp_path: Path) -> None:
    cfg = SimulationConfig(
        PencilBeamSource(ALPHA, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 400.0),
        BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), WATER),
        (ScoringGrid((-10.0, -10.0, 0.0), (5.0, 5.0, 4.0), (4, 4, 10)),),
        PhysicsOptions(nuclear=False, stopping_model="icru90-water-v1", multiple_scattering=False),
        RunOptions("warp-cpu", "float64", 17, 16, n_batches=4),
    )
    _assert_no_table_values(_case(tmp_path, cfg, None), 4)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from test_v3_006_results_format import _nuclear_table_id

    ncfg = SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 150.0, 0.0, 0.0),
        BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 200.0), WATER),
        (ScoringGrid((-10.0, -10.0, 0.0), (5.0, 5.0, 10.0), (4, 4, 20)),),
        PhysicsOptions(
            nuclear=True, nuclear_table_id=_nuclear_table_id(), stopping_model="bethe-ionmc-v1"
        ),
        RunOptions("python", "float64", 16, 40, n_batches=4),
    )
    nuc = tmp_path / "nuc"
    nuc.mkdir()
    _assert_no_table_values(_case(nuc, ncfg, None), 20)


# ---- PROV-8 (bundle part) ---------------------------------------------------------------------


def _local_copy(tmp_path: Path) -> Path:
    """The slab example with its three referenced files in one directory (editable)."""
    d = tmp_path / "cfg"
    d.mkdir()
    for name in ("slab_material_index.npy", "slab_density.npy"):
        shutil.copy(DATA / name, d / name)
    shutil.copy(
        SLAB_TOML.parents[3] / "data" / "synthetic_lookup.json", d / "synthetic_lookup.json"
    )
    text = SLAB_TOML.read_text(encoding="utf-8")
    text = text.replace("../../../data/synthetic_lookup.json", "synthetic_lookup.json")
    text = text.replace("../v3_006/slab_density.npy", "slab_density.npy")
    (d / "slab.toml").write_text(text, encoding="utf-8")
    return d


def test_prov8_toml_and_api_forms(case_d: Case, case_b: Case) -> None:
    man, doc = case_d.manifest, case_d.loaded.doc
    assert man is not None
    m = doc["input_manifest"]
    assert (m["source"], m["config_sha256"]) == ("toml", man.config_sha256)
    assert m["config_sha256"] == sha256_canonical(doc["requested"])
    assert m["config_sha256"] == config_sha256(case_d.cfg, man)  # loader == verify == save
    assert (
        to_dict(case_d.cfg, manifest=man) == doc["requested"] == to_dict(case_d.cfg, manifest=man)
    )
    for f in m["files"]:
        assert re.fullmatch(r"[0-9a-f]{64}", f["sha256"])
        assert sha256_of(Path(f["path"])) == f["sha256"]
    assert sha256_of(Path(m["config_file"]["path"])) == m["config_file"]["sha256"]
    assert doc["requested"]["geometry"]["material_index"]["ref"] == "material_index"
    # api form: unavailable pointers, and no reference to a file in `requested`
    a = case_b.loaded.doc
    assert a["input_manifest"]["source"] == "api"
    assert a["requested"]["geometry"]["material_index"]["status"] == "not-reconstructible"
    assert all(f["path"] == UNAVAILABLE for f in a["input_manifest"]["files"])
    # pointer set P: absolute paths appear only at the toml pointers, never in the api bundle
    assert str(DATA.resolve()) not in json.dumps({k: v for k, v in a.items() if k != "run"})


def sha256_of(p: Path) -> str:
    import hashlib

    return hashlib.sha256(p.read_bytes()).hexdigest()


def test_prov8_swapped_and_stale_manifests_are_refused_by_save(
    case_d: Case, tmp_path: Path
) -> None:
    cfg, man, r = case_d.cfg, case_d.manifest, case_d.result
    assert man is not None
    other = dataclasses.replace(cfg, run=dataclasses.replace(cfg.run, seed=cfg.run.seed + 1))
    # (i) the manifest of config A for config B with the very same arrays
    assert other.geometry is cfg.geometry
    swapped = dataclasses.replace(r, requested_config=other)
    with pytest.raises(ConfigError):
        save(swapped, tmp_path / "i", manifest=man)
    with pytest.raises(ConfigError):
        man.verify(other)
    # (ii) a config mutated in one scalar after loading
    src = dataclasses.replace(cfg.source, kinetic_energy_mev=cfg.source.kinetic_energy_mev + 1.0)
    mutated = dataclasses.replace(cfg, source=src)
    with pytest.raises(ConfigError):
        save(dataclasses.replace(r, requested_config=mutated), tmp_path / "ii", manifest=man)
    with pytest.raises(ConfigError):
        man.verify(mutated)
    # a manifest entry whose content_sha256 disagrees with the config arrays
    bad = dataclasses.replace(man.files[0], content_sha256="1" * 64)
    forged = dataclasses.replace(man, files=(bad, *man.files[1:]))
    with pytest.raises(ConfigError):
        save(r, tmp_path / "iii", manifest=forged)
    assert not any(p.name in ("i", "ii", "iii") for p in tmp_path.iterdir())
    assert not [p for p in tmp_path.iterdir() if p.name.startswith(".ionmc-")]  # no partial bundle
    man.verify(cfg)  # the genuine pair is accepted
    save(r, tmp_path / "ok", manifest=man)


def test_prov8_edited_input_or_config_file_between_load_and_save(
    case_d: Case, tmp_path: Path
) -> None:
    d = _local_copy(tmp_path)
    cfg, man = load_toml(d / "slab.toml")
    r = dataclasses.replace(case_d.result, requested_config=cfg)
    save(r, tmp_path / "before", manifest=man)  # unchanged files: accepted
    dens = np.load(d / "slab_density.npy")
    np.save(d / "slab_density.npy", dens * 1.0000001)  # edit a referenced file
    with pytest.raises(ConfigError):
        save(r, tmp_path / "after", manifest=man)
    np.save(d / "slab_density.npy", dens)
    (d / "slab.toml").write_text(
        (d / "slab.toml").read_text(encoding="utf-8") + "\n# edited\n", encoding="utf-8"
    )  # edit the configuration file
    with pytest.raises(ConfigError):
        save(r, tmp_path / "after2", manifest=man)
    assert not (tmp_path / "after").exists() and not (tmp_path / "after2").exists()


# ---- REG-4 / REG-5 ----------------------------------------------------------------------------


def _added_files() -> list[Path]:
    """Files added by V3-006: from git when the base is reachable, else the known set."""
    root = Path(__file__).resolve().parents[2]
    for base in ("origin/v3/develop", "v3/develop"):
        try:
            out = subprocess.run(
                ["git", "-C", str(root), "diff", "--name-only", "--diff-filter=A",
                 f"{base}...HEAD"],
                capture_output=True, text=True, timeout=60, check=True,
            ).stdout  # fmt: skip
        except (OSError, subprocess.SubprocessError):
            continue
        files = [root / ln for ln in out.splitlines() if ln and (root / ln).is_file()]
        if files:
            return files
    known = [
        *sorted((root / "src" / "ionmc").glob("configio.py")),
        *sorted((root / "src" / "ionmc").glob("results.py")),
        *sorted((root / "src" / "ionmc").glob("_result_s*.py")),
        *sorted((root / "tests" / "ionmc").glob("test_v3_006_*.py")),
    ]
    return [f for f in known if f.is_file()]


def test_reg4_no_table_values_in_files_added_by_this_task(case_d: Case, case_a: Case) -> None:
    added = _added_files() + sorted(DATA.glob("*"))
    assert added
    tables = [*_table_arrays(case_d.result.effective_config.tables)]
    tables += [*_table_arrays(case_a.result.effective_config.tables)]
    assert tables
    digests = {array_content_sha256(t) for t in tables}
    for p in added:
        assert p.suffix != ".npz", p  # no bundle committed
        if p.suffix == ".npy":
            a = np.load(p, allow_pickle=False)
            assert a.size < 2000 and array_content_sha256(a) not in digests, p
        if p.suffix in (".json", ".toml", ".csv", ".dat", ".txt"):
            assert p.stat().st_size < 600_000, p
            text = p.read_text(encoding="utf-8", errors="replace")
            for t in tables[:40]:  # a table row printed in text: first 3 values at full precision
                row = ", ".join(repr(float(x)) for x in t.ravel()[:3])
                assert row not in text, p
    assert not [p for p in added if p.suffix in (".npz", ".h5", ".hdf5", ".pkl")]


def test_reg5_single_process_only() -> None:
    assert os.environ.get("IONMC_SINGLE_PROCESS") == "1", (
        "run V3-006 tests with IONMC_SINGLE_PROCESS=1"
    )
    assert _RAN and all(c.run.cpu_workers == 1 for c in _RAN)
    for p in sorted(DATA.glob("*.toml")):
        assert not re.search(r"cpu_workers\s*=\s*([2-9]|\d\d)", p.read_text(encoding="utf-8")), p
    assert multiprocessing.active_children() == []
    pool = (
        "Pro" + "cessPoolExecutor",
        "mp." + "Pool(",
        "multiprocessing." + "Pool(",
        "set_start_" + "method",
    )
    for p in sorted(Path(__file__).parent.glob("test_v3_006_*.py")):
        text = p.read_text(encoding="utf-8")
        assert not any(w in text for w in pool), p
