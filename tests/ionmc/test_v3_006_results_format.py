"""Rows OUT-4, OUT-5, OUT-7, OUT-8 (allowlist part), OUT-9, OUT-10, OUT-11 (plan v3-006, step S2a).

``ionmc.results`` writer and reader on synthetic results: a real ``EffectiveConfig`` from a
validated tiny configuration and hand-made arrays; no transport runs (those rows are step
S2b). No physics table value, no ``out6`` seed and no worker process appears (REG-4, REG-5).
"""

from __future__ import annotations

import io
import json
import os
import resource
import shutil
import struct
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ionmc import results
from ionmc.config import PhysicsOptions, RunOptions, SimulationConfig, validate
from ionmc.configio import array_content_sha256
from ionmc.geometry import BoxPhantom
from ionmc.materials import PMMA, WATER
from ionmc.physics.projectiles import PROTON
from ionmc.results import ResultFormatError, load, save
from ionmc.scoring import QuantityResult, ScoringGrid, TallyRequest
from ionmc.simulation import EnergyBalance, GridResult, Result, TransportCounters
from ionmc.sources import PencilBeamSource

SHAPE = (3, 2, 4)
NB = 4


def _config(names: tuple[str, ...] = ("dose",), tally: bool = False) -> SimulationConfig:
    return SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 20.0),
        # the last z layer of the grid lies outside the box: its voxel mass is 0 (dose is NaN)
        BoxPhantom((-3.0, -2.0, 0.0), (6.0, 4.0, 6.0), WATER),
        tuple(ScoringGrid((-3.0, -2.0, 0.0), (2.0, 2.0, 2.0), SHAPE, n) for n in names),
        PhysicsOptions(nuclear=False, stopping_model="bethe-ionmc-v1"),
        RunOptions("python", "float64", 11, 8, n_batches=NB),
        tallies=(TallyRequest("let", names[0], "let_d"),) if tally else (),
    )


def _batch_sums(seed: int = 0) -> np.ndarray[Any, Any]:
    """Synthetic per-batch sums ``[NB, nvox]``; voxel 0 has no deposit in any batch."""
    rng = np.random.default_rng(seed)  # synthetic arrays only
    sums = rng.uniform(0.5, 2.0, (NB, int(np.prod(SHAPE)))) * 2.0
    sums[:, 0] = 0.0
    sums[0, 1] = 0.0  # one voxel with a zero batch
    return sums


def _grid(cfg: SimulationConfig, eff: Any, gi: int = 0) -> GridResult:
    """A grid result built by the engine's own reduction (consistent derived arrays)."""
    from ionmc.simulation import _grid_result

    return _grid_result(cfg.scoring[gi], _batch_sums(gi), eff, True)


def _quantity(shape: tuple[int, ...] = SHAPE) -> QuantityResult:
    mean = np.linspace(1.0, 2.0, int(np.prod(shape))).reshape(shape)
    mask = np.ones(shape, dtype=bool)
    mean[0, 0, 0], mask[0, 0, 0] = np.nan, False
    return QuantityResult(
        "let", "let_d", "ratio", "MeV/mm", "sum(l S^2)/sum(l S), S in water",
        mean, 0.1 * mean, mask,
        np.where(mask, 4, 0).astype(np.int32), np.where(mask, 1e-9, np.nan), (3, 3),
    )  # fmt: skip


def _pieces() -> dict[str, QuantityResult]:
    """The automatic channels of a tally grid: the excluded-energy channel (the compiled plan of a
    LET tally has a local-class energy channel) and the two piece counts."""
    spec = (
        ("edep_excluded_from_let", "edep", "MeV", results._EXCLUDED_DEFINITION),
        ("scoring_pieces", "count", "1", results._PIECE_DEFINITIONS["scoring_pieces"]),
        ("scoring_pieces_local", "count", "1", results._PIECE_DEFINITIONS["scoring_pieces_local"]),
    )
    out = {}
    for n, q, u, d in spec:
        one = np.ones(SHAPE)
        out[n] = QuantityResult(
            n, q, "linear", u, d, one, 0.0 * one, np.ones(SHAPE, dtype=bool),
            np.full(SHAPE, 4, dtype=np.int32), 0.0 * one, (3,),
        )  # fmt: skip
    return out


def make_result(
    names: tuple[str, ...] = ("dose",),
    tally: bool = False,
    diagnostics: dict[str, Any] | None = None,
    cfg: SimulationConfig | None = None,
) -> Result:
    cfg = cfg or _config(names, tally)
    eff = validate(cfg)
    grids = []
    for i, _n in enumerate(names):
        g = _grid(cfg, eff, i)
        if tally and i == 0:
            g = GridResult(*[getattr(g, f) for f in g.__dataclass_fields__ if f != "quantities"],
                           quantities={"let": _quantity(), **_pieces()})  # fmt: skip
        grids.append(g)
    in_grid = tuple(float(np.sum(g.energy_mev)) * 8 for g in grids)  # bound to the arrays
    total = max(in_grid) + 1.0  # in_grid + outside = step_deposit + cutoff on every grid
    bal = EnergyBalance(
        total + 0.25, total - 0.5, 0.5, 0.25, 0.0, 0.0, in_grid,
        tuple(total - x for x in in_grid),
    )  # fmt: skip
    return Result(
        True, cfg, eff, "python", "cpu", "float64", 11, dict(eff.rng), 8, NB,
        tuple(grids), bal, TransportCounters(), {"total": 0.1},
        {"ionmc_version": "0", "python_version": "3", "numpy_version": np.__version__,
         "warp_version": None, "warp_error": None, "warp_devices": [],
         "git_sha": None, "git_dirty": None},
        diagnostics or {},
    )  # fmt: skip


# ---- helpers on bundles ----------------------------------------------------------------------


def _doc(b: Path) -> dict[str, Any]:
    return json.loads((b / "result.json").read_text(encoding="utf-8"))


def _put_doc(b: Path, doc: dict[str, Any]) -> None:
    (b / "result.json").write_bytes(results.canonical_bytes(doc))


def _refile(b: Path) -> None:
    """Recompute the recorded file hashes (an adversary who also fixes ``files``)."""
    doc = _doc(b)
    for name in doc["files"]:
        sha, size = results._sha256_path(b / name)
        doc["files"][name] = {"sha256": sha, "size": size}
    _put_doc(b, doc)


def _members(b: Path) -> dict[str, bytes]:
    with zipfile.ZipFile(b / "arrays.npz") as zf:
        return {i.filename: zf.read(i) for i in zf.infolist()}


def _npy(arr: np.ndarray[Any, Any]) -> bytes:
    buf = io.BytesIO()
    np.lib.format.write_array(buf, arr, allow_pickle=False)
    return buf.getvalue()


def _write_zip(b: Path, members: list[tuple[str, bytes]]) -> None:
    with zipfile.ZipFile(b / "arrays.npz", "w", zipfile.ZIP_DEFLATED) as zf:
        for n, data in members:
            zf.writestr(zipfile.ZipInfo(n, (1980, 1, 1, 0, 0, 0)), data)
    _refile(b)


@pytest.fixture
def bundle(tmp_path: Path) -> Path:
    return save(make_result(tally=True), tmp_path / "b")


def _mutate(b: Path, fn: Any) -> None:
    d = _doc(b)
    fn(d)
    _put_doc(b, d)


def _peak_kb() -> int:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


# ---- positive control and OUT-4 ---------------------------------------------------------------


def test_round_trip_bitwise(tmp_path: Path) -> None:
    r = make_result(tally=True)
    got = load(save(r, tmp_path / "b"))
    g = r.grids[0]
    for base in (
        "energy_mev",
        "dose_gy",
        "relative_uncertainty",
        "defined_mask",
        "batch_energy_mev",
    ):
        a, b = getattr(g, base), got.grid_array("dose", base)
        assert a.dtype == b.dtype and a.shape == b.shape
        assert np.array_equal(a, b, equal_nan=True)
    assert got.doc["schema"] == "ionmc.result" and got.doc["valid"] is True
    assert (tmp_path / "b" / "result.json").read_bytes().endswith(b"}\n")


def _unknown_top(d: dict[str, Any]) -> None:
    d["surprise"] = 1


def _unknown_nested(d: dict[str, Any]) -> None:
    d["run"]["surprise"] = 1


def _unknown_entry(d: dict[str, Any]) -> None:
    d["arrays"][0]["surprise"] = 1


def _missing_key(d: dict[str, Any]) -> None:
    del d["run"]["seed"]


def _v11(d: dict[str, Any]) -> None:
    d["schema_version"] = "1.1"


def _v20(d: dict[str, Any]) -> None:
    d["schema_version"] = "2.0"


def _missing_hash(d: dict[str, Any]) -> None:
    del d["files"]["report.md"]


@pytest.mark.parametrize(
    ("fn", "match"),
    [
        (_unknown_top, "unknown field"),
        (_unknown_nested, "unknown field"),
        (_unknown_entry, "unknown field"),
        (_missing_key, "missing required"),
        (_v11, "unknown schema_version"),
        (_v20, "unknown schema_version"),
        (_missing_hash, "unexpected entry"),
    ],
)
def test_out4_document_negative_controls(bundle: Path, fn: Any, match: str) -> None:
    _mutate(bundle, fn)
    with pytest.raises(ResultFormatError, match=match):
        load(bundle)


def test_out4_flipped_byte(bundle: Path) -> None:
    p = bundle / "arrays.npz"
    raw = bytearray(p.read_bytes())
    raw[len(raw) // 2] ^= 0xFF
    p.write_bytes(bytes(raw))
    with pytest.raises(ResultFormatError, match="sha256"):
        load(bundle)


def test_out4_altered_member_with_recomputed_file_hash(bundle: Path) -> None:
    m = _members(bundle)
    name = "grid__dose__energy_mev.npy"
    arr = np.load(io.BytesIO(m[name]), allow_pickle=False).copy()
    arr[1, 1, 1] += 1.0
    m[name] = _npy(arr)
    _write_zip(bundle, sorted(m.items()))
    with pytest.raises(ResultFormatError, match="content_sha256"):
        load(bundle)


def test_out4_nan_literal(bundle: Path) -> None:
    p = bundle / "result.json"
    p.write_bytes(p.read_bytes().replace(b'"seed": 11', b'"seed": NaN', 1))
    with pytest.raises(ResultFormatError, match="NaN"):
        load(bundle)


def test_out4_array_shape_disagrees_with_grid(tmp_path: Path) -> None:
    r = make_result()
    g = r.grids[0]
    bad = GridResult(*[getattr(g, f) for f in g.__dataclass_fields__ if f != "energy_mev"][:0] or [
        g.name, np.zeros((3, 2, 5)), g.energy_std_mev, g.dose_gy, g.dose_std_gy,
        g.relative_uncertainty, g.defined_mask, g.mass_g, g.n_nonzero_batches, g.sparse_mask,
        g.low_batch_count, g.batch_energy_mev, g.valid])  # fmt: skip
    r2 = Result(**{**r.__dict__, "grids": (bad,)})
    with pytest.raises(ResultFormatError, match="disagrees with grid|does not match"):
        load(save(r2, tmp_path / "b"))


def test_out4_extra_npz_member(bundle: Path) -> None:
    m = _members(bundle)
    m["extra.npy"] = _npy(np.zeros(3))
    _write_zip(bundle, sorted(m.items()))
    with pytest.raises(ResultFormatError, match="entries|not in the catalogue"):
        load(bundle)


# ---- OUT-5 ------------------------------------------------------------------------------------


def test_out5_units(bundle: Path) -> None:
    doc = _doc(bundle)
    for e in doc["arrays"]:
        assert e["units"] in results.UNIT_VOCABULARY, e["name"]
    for g in doc["statistics"]["grids"]:
        for k in ("energy_total_mev", "dose_max_gy", "relative_uncertainty_max"):
            assert g[k]["units"] in results.UNIT_VOCABULARY
    for k in doc["accounting"]["energy_balance"]:
        assert k.endswith("_mev") or k in ("in_grid_mev", "outside_mev"), k
    _mutate(bundle, lambda d: d["arrays"][0].update(units="parsec"))
    with pytest.raises(ResultFormatError, match="vocabulary"):
        load(bundle)


# ---- OUT-7 ------------------------------------------------------------------------------------


def test_out7_undefined_is_nan_and_scalars_carry_status(tmp_path: Path) -> None:
    got = load(save(make_result(), tmp_path / "b"))
    mask = got.grid_array("dose", "defined_mask")
    assert not mask[0, 0, 0] and np.isnan(got.grid_array("dose", "relative_uncertainty")[0, 0, 0])
    assert np.isnan(got.grid_array("dose", "dose_gy")[-1, -1, -1])
    r = make_result()
    g = r.grids[0]
    z = g.energy_mev * 0
    d0 = np.where(g.mass_g > 0, 0.0, np.nan)
    none = GridResult(
        g.name, z, z, d0, d0, np.full(SHAPE, np.nan), mask & False, g.mass_g,
        np.zeros(SHAPE, dtype=np.int32), g.sparse_mask & False, NB < 10, g.batch_energy_mev * 0,
        True,
    )  # fmt: skip
    eb = EnergyBalance(8.0, 7.0, 0.5, 0.25, 0.0, 0.0, (0.0,), (7.5,))  # no energy in the grid
    doc = load(
        save(Result(**{**r.__dict__, "grids": (none,), "energy_balance": eb}), tmp_path / "c")
    ).doc
    rec = doc["statistics"]["grids"][0]["dose_max_gy"]
    assert rec["value"] is None and rec["status"] == "undefined" and rec["reason"]
    _mutate(
        tmp_path / "c",
        lambda d: d["statistics"]["grids"][0].update(dose_max_gy={"value": 0.0, "units": "Gy"}),
    )
    with pytest.raises(ResultFormatError, match="undefined record"):
        load(tmp_path / "c")


def test_out7_patched_writer_writing_zero_is_rejected(tmp_path: Path, monkeypatch: Any) -> None:
    real = results._collect_arrays

    def zeroing(result: Result, persist: bool) -> list[Any]:
        out = []
        for a in real(result, persist):
            if a.origin.get("array") == "relative_uncertainty":
                a = results._Arr(a.name, np.nan_to_num(a.data, nan=0.0), a.axes, a.grid,
                                 a.units, a.definition, a.origin)  # fmt: skip
            out.append(a)
        return out

    monkeypatch.setattr(results, "_collect_arrays", zeroing)
    with pytest.raises(ResultFormatError, match="NaN, never 0"):
        load(save(make_result(), tmp_path / "b"))


# ---- OUT-8 (allowlist and injection) ----------------------------------------------------------


def test_out8_names_inside_allowlist_and_injection_refused(
    tmp_path: Path, bundle: Path, monkeypatch: Any
) -> None:
    with zipfile.ZipFile(bundle / "arrays.npz") as zf:
        names = [i.filename[:-4] for i in zf.infolist()]
    entries = _doc(bundle)["arrays"]
    assert sorted(names) == sorted(e["name"] for e in entries)
    results.check_allowlist(entries)
    for key in ("tables.range_mm", "range_mm", "stopping_power"):
        r = make_result(diagnostics={key: np.zeros(3)})
        with pytest.raises(ResultFormatError, match="allowlist"):
            save(r, tmp_path / "x")
    assert not (tmp_path / "x").exists()
    real = results._collect_arrays

    def inject(result: Result, persist: bool) -> list[Any]:
        extra = results._Arr("tables__range_mm", np.zeros(3), ("e",), None, "mm", "x",
                             {"kind": "diagnostic", "array": "range_mm"})  # fmt: skip
        return [*real(result, persist), extra]

    monkeypatch.setattr(results, "_collect_arrays", inject)
    with pytest.raises(ResultFormatError, match="allowlist"):
        save(make_result(), tmp_path / "y")
    assert sorted(p.name for p in tmp_path.iterdir()) == ["b"]  # no partial bundle


# ---- OUT-9 ------------------------------------------------------------------------------------


def test_out9_deterministic_bytes_and_level(tmp_path: Path) -> None:
    r = make_result(tally=True)
    b = save(r, tmp_path / "b")
    first = {n: (b / n).read_bytes() for n in ("result.json", "arrays.npz", "report.md")}
    save(r, b, overwrite=True)
    assert first == {n: (b / n).read_bytes() for n in first}
    with pytest.raises(FileExistsError):
        save(r, b)
    c = save(r, tmp_path / "c", zlib_level=0)
    assert (c / "arrays.npz").read_bytes() != first["arrays.npz"]
    da, dc = _doc(b), _doc(c)
    assert da["files"]["arrays.npz"]["sha256"] != dc["files"]["arrays.npz"]["sha256"]
    assert [e["content_sha256"] for e in da["arrays"]] == [
        e["content_sha256"] for e in dc["arrays"]
    ]
    assert load(c).dump()["content_sha256"] == load(b).dump()["content_sha256"]


# ---- OUT-10 -----------------------------------------------------------------------------------


def _forbid(monkeypatch: Any, target: Any, attr: str) -> None:
    def boom(*a: Any, **k: Any) -> Any:
        raise AssertionError(f"{attr} must not be called")

    monkeypatch.setattr(target, attr, boom)


def _npy_header_bytes(arr: np.ndarray[Any, Any], pad_to: int) -> bytes:
    hdr = f"{{'descr': '{arr.dtype.str}', 'fortran_order': False, 'shape': {arr.shape}, }}"
    hdr = hdr.ljust(pad_to - 11) + "\n"
    return b"\x93NUMPY\x01\x00" + struct.pack("<H", len(hdr)) + hdr.encode() + arr.tobytes()


def _swap_member(b: Path, name: str, data: bytes) -> None:
    m = _members(b)
    m[name] = data
    _write_zip(b, sorted(m.items()))


def test_out10_member_cases(bundle: Path) -> None:
    name = "grid__dose__energy_mev.npy"
    arr = np.load(io.BytesIO(_members(bundle)[name]))
    pristine = (bundle / "arrays.npz").read_bytes()
    cases: dict[str, bytes] = {
        "long header": _npy_header_bytes(arr, 4000),
        "disagreeing header": _npy(arr.astype(np.float32)),
        "wrong shape header": _npy(arr.reshape(-1)),
        "huge declared member": _npy(arr) + bytes(3_000_000),
        "truncated member": _npy(arr)[:-8],
    }
    for label, data in cases.items():
        (bundle / "arrays.npz").write_bytes(pristine)
        _swap_member(bundle, name, data)
        with pytest.raises(ResultFormatError):
            load(bundle)
        assert label


def test_out10_duplicate_member_and_extra_members(bundle: Path) -> None:
    m = list(_members(bundle).items())
    with pytest.warns(UserWarning):
        _write_zip(bundle, [*m, m[0]])
    with pytest.raises(ResultFormatError):
        load(bundle)
    _write_zip(bundle, [*m, *[(f"x{i}.npy", b"") for i in range(len(m) + 1)]])
    with pytest.raises(ResultFormatError, match="declares|catalogue"):
        load(bundle)


def test_out10_declared_shape_beyond_budget(bundle: Path) -> None:
    before = _peak_kb()
    _mutate(bundle, lambda d: d["arrays"][0].update(shape=[2**40, 1, 1]))
    with pytest.raises(ResultFormatError, match="max_array_bytes"):
        load(bundle)
    assert _peak_kb() - before < 64 * 1024
    _mutate(bundle, lambda d: d["arrays"][0].update(shape=[1, 1, 1, 1, 1]))
    with pytest.raises(ResultFormatError, match="invalid shape"):
        load(bundle)


def test_out10_oversized_json_and_budget_keywords(bundle: Path) -> None:
    with pytest.raises(ResultFormatError, match="max_json_bytes"):
        load(bundle, max_json_bytes=100)
    with pytest.raises(ResultFormatError, match="max_array_bytes"):
        load(bundle, max_array_bytes=8)
    with pytest.raises(ResultFormatError, match="max_total_bytes"):
        load(bundle, max_total_bytes=100)


def _zip_eocd(n: int, cd_size: int, cd_off: int) -> bytes:
    return struct.pack("<4s4H2LH", b"PK\x05\x06", 0, 0, n, n, cd_size, cd_off, 0)


def _zip64(n: int, cd_size: int, cd_off: int) -> bytes:
    rec = struct.pack("<4sQ2H2L4Q", b"PK\x06\x06", 44, 45, 45, 0, 0, n, n, cd_size, cd_off)
    loc = struct.pack("<4sLQL", b"PK\x06\x07", 0, 0, 1)
    return rec + loc + _zip_eocd(0xFFFF, 0xFFFFFFFF, 0xFFFFFFFF)


def test_out10_directory_bombs_refused_before_zipfile(bundle: Path, monkeypatch: Any) -> None:
    n_cat = len(_doc(bundle)["arrays"])
    bombs = {
        "entry count": _zip_eocd(60000, 46, 0),
        "directory size": _zip_eocd(n_cat, 2**30, 0),
        "offset": _zip_eocd(n_cat, 46, 2**31),
        "zip64 entry count": _zip64(2**40, 100, 0),
        "zip64 directory size": _zip64(1, 2**40, 0),
        "no eocd": b"PK" * 40,
    }
    _forbid(monkeypatch, zipfile, "ZipFile")
    for label, data in bombs.items():
        (bundle / "arrays.npz").write_bytes(data)
        _refile(bundle)
        before = _peak_kb()
        with pytest.raises(ResultFormatError):
            load(bundle)
        assert _peak_kb() - before < 64 * 1024, label


def _json_bomb(kind: str) -> bytes:
    if kind == "objects":
        return b"[" + b'{"a":1},' * 600_000 + b"{}]"
    if kind == "empties":
        return b"[" + b"{}," * 1_900_000 + b"{}]"
    if kind == "depth":
        return b"[" * 65 + b"]" * 65
    if kind == "numbers":
        return b"[" + b"1," * 2_100_000 + b"1]"
    return b"[" + b"true,false,null," * 700_000 + b"null]"


@pytest.mark.parametrize("kind", ["objects", "depth", "numbers", "literals", "empties"])
def test_out10_json_bombs_refused_before_json_loads(
    bundle: Path, monkeypatch: Any, kind: str
) -> None:
    (bundle / "result.json").write_bytes(_json_bomb(kind))
    _forbid(monkeypatch, json, "loads")
    before = _peak_kb()
    with pytest.raises(ResultFormatError, match="token|depth"):
        load(bundle)
    assert _peak_kb() - before < 64 * 1024


def test_out10_positive_control_at_default_budgets(bundle: Path) -> None:
    assert load(bundle).doc["kind"] == "transport"


# ---- OUT-11 -----------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "cfg_names",
    [("../up", "a/b", "A", "a"), ("a\\b", "UPPER", "x" * 80, "nul\0name")],
)
def test_out11_writer_names(tmp_path: Path, cfg_names: tuple[str, ...]) -> None:
    b = save(make_result(cfg_names), tmp_path / "out")
    doc = _doc(b)
    members = [e["name"] for e in doc["arrays"]]
    assert len(set(members)) == len(members)
    for m in members:
        assert results._MEMBER_RE.fullmatch(m + ".npy"), m
    assert {e["origin"]["grid"] for e in doc["arrays"] if e["grid"]} == set(cfg_names)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["out"]
    assert sorted(p.name for p in b.iterdir()) == ["arrays.npz", "report.md", "result.json"]
    assert len(load(b).arrays) == len(members)


def test_out11_writer_refuses_unsafe_entries() -> None:
    bad = {"name": "../x", "origin": {"kind": "grid", "array": "energy_mev"}, "units": "MeV"}
    with pytest.raises(ResultFormatError):
        results.check_allowlist([bad])


@pytest.mark.parametrize(
    "name", ["../evil.txt", "/abs/evil.txt", "./a.txt", "a/b.txt", "a\0.txt", "plots/../x.svg"]
)
def test_out11_reader_files_entries(bundle: Path, tmp_path: Path, name: str) -> None:
    outside = tmp_path / "evil.txt"
    outside.write_text("x")
    _mutate(bundle, lambda d: d["files"].update({name: {"sha256": "0" * 64, "size": 1}}))
    with pytest.raises(ResultFormatError, match="unsafe file name"):
        load(bundle)


@pytest.mark.parametrize("member", ["../x.npy", "/x.npy", "./x.npy", "a/b.npy", "a\0.npy"])
def test_out11_reader_member_names(bundle: Path, member: str) -> None:
    m = list(_members(bundle).items())
    _write_zip(bundle, [*m, (member, b"")])
    with pytest.raises(ResultFormatError):
        load(bundle)


def test_out11_reader_refuses_planted_symlinks(bundle: Path, tmp_path: Path) -> None:
    report = bundle / "report.md"
    outside = tmp_path / "outside.md"
    outside.write_bytes(report.read_bytes())  # same bytes, hence same sha256
    report.unlink()
    os.symlink(outside, report)
    with pytest.raises(ResultFormatError, match="symbolic link|regular"):
        load(bundle)
    report.unlink()
    os.symlink(bundle / "arrays.npz", report)  # a link to a file inside the bundle
    with pytest.raises(ResultFormatError, match="symbolic link|regular"):
        load(bundle)
    report.unlink()
    outside.rename(report)
    assert load(bundle).doc["schema"] == "ionmc.result"


# ---- strict blocks: requested, input_manifest, effective --------------------------------------


def _nested_unknown(path: tuple[str, ...]) -> Any:
    def fn(d: dict[str, Any]) -> None:
        node: Any = d
        for k in path:
            node = node[k]
        node["surprise"] = 1

    return fn


@pytest.mark.parametrize(
    ("path", "match"),
    [
        (("requested", "run"), "requested: .*run.surprise"),
        (("requested",), "requested: .*surprise"),
        (("input_manifest",), "unknown field"),
        (("effective",), "effective: unknown field"),
    ],
)
def test_strict_blocks_refuse_unknown_keys(bundle: Path, path: tuple[str, ...], match: str) -> None:
    _mutate(bundle, _nested_unknown(path))
    with pytest.raises(ResultFormatError, match=match):
        load(bundle)


def test_manifest_file_entry_unknown_key(bundle: Path) -> None:
    _mutate(bundle, lambda d: d["input_manifest"]["files"].append({"surprise": 1}))
    with pytest.raises(ResultFormatError, match="input_manifest.files"):
        load(bundle)


def test_effective_top_keys_equal_the_summaries_fixture() -> None:
    fx = json.loads(
        (Path(__file__).parent / "data/v3_006/summaries-95937553.json").read_text(encoding="utf-8")
    )
    union = {k for c in fx["cases"].values() for k in c["summary"]}
    assert union == set(results.EFFECTIVE_TOP_KEYS)


# ---- review REVIEW-741411bd regressions -------------------------------------------------------


def _invalid(r: Result, **kw: Any) -> Result:
    """``r`` as an invalid result: the engine marks every grid invalid with the run."""
    import dataclasses

    grids = tuple(dataclasses.replace(g, valid=False) for g in r.grids)
    return Result(**{**r.__dict__, "valid": False, "grids": grids, **kw})


def _with_channels(lookup: int, path: int) -> Result:
    from ionmc.transport.tally import ChannelRaw

    r = make_result()
    raw = ChannelRaw(np.zeros((NB, 1), dtype=np.int64), [], lookup, path)
    return _invalid(r, channel_raw=raw)


@pytest.mark.parametrize("which", ["lookup_out_of_domain", "path_bound_exceeded"])
def test_item1_invalid_by_channel_counter_round_trips(tmp_path: Path, which: str) -> None:
    r = _with_channels(*((2, 0) if which == "lookup_out_of_domain" else (0, 3)))
    doc = load(save(r, tmp_path / "b")).doc
    assert doc["valid"] is False and doc["accounting"]["reasons"] == [which]
    assert doc["accounting"]["channel_counters"][which] > 0
    _mutate(tmp_path / "b", lambda d: d["accounting"]["channel_counters"].update({which: 0}))
    with pytest.raises(ResultFormatError, match="valid disagrees"):
        load(tmp_path / "b")


def test_item2_swapped_grid_labels_and_renamed_array(tmp_path: Path) -> None:
    b = save(make_result(("g1", "g2")), tmp_path / "b")
    pristine = _doc(b)

    def swap(d: dict[str, Any]) -> None:
        other = {"g1": "g2", "g2": "g1"}
        for e in d["arrays"]:
            if e["grid"]:
                e["grid"] = e["origin"]["grid"] = other[e["grid"]]
        for g in d["statistics"]["grids"]:
            g["name"] = other[g["name"]]

    _mutate(b, swap)
    with pytest.raises(ResultFormatError, match="semantic origin"):
        load(b)
    _put_doc(b, pristine)

    def rename(d: dict[str, Any]) -> None:
        for e in d["arrays"]:
            if e["origin"].get("array") == "energy_mev" and e["grid"] == "g1":
                e["origin"]["array"] = "energy_std_mev"

    _mutate(b, rename)
    with pytest.raises(ResultFormatError, match="allowlist|duplicate semantic origin"):
        load(b)


def _drop_member(b: Path, name: str) -> None:
    m = {k: v for k, v in _members(b).items() if k != name + ".npy"}
    _write_zip(b, sorted(m.items()))
    _mutate(b, lambda d: d.update(arrays=[e for e in d["arrays"] if e["name"] != name]))


@pytest.mark.parametrize(
    "name",
    [
        "grid__dose__energy_std_mev",
        "grid__dose__n_nonzero_batches",
        "grid__dose__sparse_mask",
        "grid__dose__batch_energy_mev",
        "grid__dose__q__let__rounding_bound",
        "geometry__material_index",
    ],
)
def test_item3_omitted_arrays_are_refused(bundle: Path, name: str) -> None:
    _drop_member(bundle, name)
    with pytest.raises(ResultFormatError, match="must be present|missing|geometry"):
        load(bundle)


def test_item3_batches_declared_not_persisted_must_be_absent(tmp_path: Path) -> None:
    b = save(make_result(), tmp_path / "b", persist_batches=False)
    assert load(b).doc["statistics"]["batch_arrays"]["status"] == "not-persisted"
    _mutate(b, lambda d: d["statistics"].update(batch_arrays={"status": "persisted"}))
    with pytest.raises(ResultFormatError, match="batch_energy_mev"):
        load(b)


def test_item4_file_tree(bundle: Path) -> None:
    def register(name: str, data: bytes) -> None:
        (bundle / name).parent.mkdir(exist_ok=True)
        (bundle / name).write_bytes(data)
        _refile_add(bundle, name)

    def _refile_add(b: Path, name: str) -> None:
        d = _doc(b)
        sha, size = results._sha256_path(b / name)
        d["files"][name] = {"sha256": sha, "size": size}
        _put_doc(b, d)

    register("plots/p.svg", b"<svg/>")
    assert "plots/p.svg" in load(bundle).doc["files"]  # a listed plot is fine
    (bundle / "plots" / "stray.svg").write_bytes(b"x")
    with pytest.raises(ResultFormatError, match="unexpected entry"):
        load(bundle)
    (bundle / "plots" / "stray.svg").unlink()
    os.symlink(bundle / "arrays.npz", bundle / "plots" / "link.svg")
    with pytest.raises(ResultFormatError, match="symbolic link"):
        load(bundle)
    (bundle / "plots" / "link.svg").unlink()
    register("extra.txt", b"x")
    with pytest.raises(ResultFormatError, match="unsafe file name"):
        load(bundle)


def test_item4_report_is_required(bundle: Path) -> None:
    (bundle / "report.md").unlink()
    _mutate(bundle, lambda d: d["files"].pop("report.md"))
    with pytest.raises(ResultFormatError, match="report.md"):
        load(bundle)


@pytest.mark.parametrize(
    "fn",
    [
        lambda d: d["arrays"][0].update(axes=5),
        lambda d: d["arrays"][0].update(axes=[1, 2, 3]),
        lambda d: d["accounting"].update(counters="x"),
        lambda d: d["accounting"]["counters"].update(stall="many"),
        lambda d: d["accounting"].update(reasons=7),
        lambda d: d["run"].update(n_batches="four"),
        lambda d: d["statistics"].update(batch_arrays=[]),
        lambda d: d["effective"].update(scoring=5),
        lambda d: d["effective"].update(scoring=[{"name": 1, "shape": 2}]),
    ],
)
def test_item5_malformed_types_are_format_errors(bundle: Path, fn: Any) -> None:
    _mutate(bundle, fn)
    with pytest.raises(ResultFormatError):
        load(bundle)


def test_item6_output_and_stopping_record_schemas(bundle: Path) -> None:
    from ionmc.configio import ConfigError, to_dict, validate_dict_structure

    good = to_dict(_config())
    good["output"] = {"path": None, "persist_batches": True}
    validate_dict_structure(good)
    rec = {
        "status": "not-reconstructible",
        "reason": "r",
        "reconstructible": False,
        "stopping_object": {"class": "C", "model_id": "m", "spec_sha256": "0" * 64},
    }
    user = {**good["physics"], "stopping_model": None}
    validate_dict_structure({**good, "physics": {**user, "stopping": rec}})
    with pytest.raises(ConfigError, match="stopping_model"):
        validate_dict_structure({**good, "physics": {**good["physics"], "stopping": rec}})
    bad_rec = {**rec, "stopping_object": {**rec["stopping_object"], "surprise": "x"}}
    for physics in (
        {"stopping": {"surprise": 1}},
        {"stopping": bad_rec},
        {"stopping": {**rec, "x": 1}},
    ):
        with pytest.raises(ConfigError, match="physics.stopping"):
            validate_dict_structure({**good, "physics": {**user, **physics}})
    with pytest.raises(ConfigError, match=r"output\.surprise"):
        validate_dict_structure({**good, "output": {"surprise": 1}})
    _mutate(bundle, lambda d: d["requested"].update(output={"surprise": 1}))
    with pytest.raises(ResultFormatError, match="output"):
        load(bundle)
    _mutate(bundle, lambda d: d["requested"].update(output={}))
    _mutate(
        bundle,
        lambda d: d["requested"]["physics"].update(stopping={"surprise": 1}, stopping_model=None),
    )
    with pytest.raises(ResultFormatError, match="physics.stopping"):
        load(bundle)


def test_duplicate_semantic_origin_is_refused(bundle: Path) -> None:
    m = _members(bundle)
    src = "grid__dose__energy_mev"
    m["grid__zz__energy_mev.npy"] = m[src + ".npy"]  # same bytes under a second allowlisted name
    _write_zip(bundle, sorted(m.items()))

    def dup(d: dict[str, Any]) -> None:
        e = next(x for x in d["arrays"] if x["name"] == src)
        d["arrays"] = sorted(
            [*d["arrays"], {**e, "name": "grid__zz__energy_mev"}], key=lambda x: x["name"]
        )

    _mutate(bundle, dup)
    with pytest.raises(ResultFormatError, match="duplicate semantic origin"):
        load(bundle)


def test_unknown_budget_keyword_is_a_type_error(bundle: Path) -> None:
    with pytest.raises(TypeError):
        load(bundle, max_everything=1)  # type: ignore[call-arg]


# ---- review REVIEW-9469eb55 regressions -------------------------------------------------------

import dataclasses  # noqa: E402

from ionmc.configio import sha256_canonical  # noqa: E402
from ionmc.geometry import VoxelGeometry  # noqa: E402
from ionmc.physics.projectiles import CARBON12  # noqa: E402
from ionmc.simulation import (  # noqa: E402
    IonTransportCounters,
    NuclearEnergyBalance,
    Simulation,
)

_RUN_CACHE = os.environ.get("IONMC_CACHE_DIR")


def _ion_config(backend: str = "warp-cpu") -> SimulationConfig:
    return SimulationConfig(
        PencilBeamSource(CARBON12, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 1200.0),
        BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), WATER),
        (ScoringGrid((-10.0, -10.0, 0.0), (5.0, 5.0, 4.0), (4, 4, 10)),),
        PhysicsOptions(nuclear=False, stopping_model="bethe-ionmc-v1", multiple_scattering=False),
        RunOptions(backend, "float64", 15, 8, n_batches=4),
    )


def test_item1_ion_run_round_trip(tmp_path: Path) -> None:
    r = Simulation(_ion_config()).run()
    assert isinstance(r.counters, IonTransportCounters)
    got = load(save(r, tmp_path / "b"))
    doc = got.doc
    assert set(doc["accounting"]["counters"]) == set(r.counters.as_dict())
    assert doc["valid"] is True and got.effective_binding == "bundle-local"


def test_item1_bracket_miss_is_diagnostic_only(tmp_path: Path) -> None:
    cfg = dataclasses.replace(
        _config(),
        source=PencilBeamSource(CARBON12, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 1200.0),
        physics=PhysicsOptions(
            nuclear=False, stopping_model="bethe-ionmc-v1", multiple_scattering=False
        ),
    )
    r = make_result(cfg=cfg)
    r = Result(**{**r.__dict__, "counters": IonTransportCounters(mcs_table_bracket_miss=5)})
    doc = load(save(r, tmp_path / "b")).doc
    assert doc["valid"] is True and doc["accounting"]["reasons"] == []
    assert doc["accounting"]["counters"]["mcs_table_bracket_miss"] == 5
    bad = _invalid(r, counters=IonTransportCounters(mcs_mapping_failure=2))
    d2 = load(save(bad, tmp_path / "c")).doc
    assert d2["accounting"]["reasons"] == ["mcs_mapping_failure"]
    # the base counter schema is refused for an ion run
    _mutate(tmp_path / "b", lambda d: d["accounting"]["counters"].pop("mcs_mapping_failure"))
    with pytest.raises(ResultFormatError, match="keys must be exactly"):
        load(tmp_path / "b")


def _nuclear_table_id() -> str:
    from ionmc.data import cache
    from ionmc.nuclear.tables import NuclearTable, NuclearTableError

    cdir = cache.resolve_cache_dir(None)
    for p in sorted((cdir / "derived").glob("nuclear-proton-*.json")):
        tid = p.stem.removeprefix("nuclear-proton-")
        try:
            info = NuclearTable.load(None, tid).info
        except NuclearTableError:
            continue
        if "transport_path_bound_terms" in info:
            return tid
    pytest.skip("no built nuclear table in the cache (IONMC_CACHE_DIR)")


def test_item1_nuclear_run_round_trip(tmp_path: Path) -> None:
    tid = _nuclear_table_id()
    cfg = SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 150.0, 0.0, 0.0),
        BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 200.0), WATER),
        (ScoringGrid((-10.0, -10.0, 0.0), (5.0, 5.0, 10.0), (4, 4, 20)),),
        PhysicsOptions(nuclear=True, nuclear_table_id=tid, stopping_model="bethe-ionmc-v1"),
        RunOptions("python", "float64", 16, 40, n_batches=4),
    )
    r = Simulation(cfg).run()
    assert isinstance(r.energy_balance, NuclearEnergyBalance)
    got = load(save(r, tmp_path / "b"))
    doc = got.doc
    assert "majorant_violation" in doc["accounting"]["counters"]
    assert "nuclear_mev" in doc["accounting"]["energy_balance"]
    assert got.effective_binding == "bundle-local"
    # A-1 item 4: the exact nuclear key set and the signed nuclear terms in the closure
    b = tmp_path / "b"
    _mutate(
        b, lambda d: d["accounting"]["energy_balance"]["nuclear_mev"].pop("nuclear_alpha_local")
    )
    with pytest.raises(ResultFormatError, match="nuclear_mev"):
        load(b)
    _put_doc(b, doc)
    _mutate(b, lambda d: d["accounting"]["energy_balance"]["nuclear_mev"].update(
        nuclear_binding=doc["accounting"]["energy_balance"]["nuclear_mev"]["nuclear_binding"] + 1.0
    ))  # fmt: skip
    with pytest.raises(ResultFormatError, match="closure_residual"):
        load(b)


def _put_members(b: Path, extra: dict[str, bytes]) -> None:
    m = {**_members(b), **extra}
    _write_zip(b, sorted(m.items()))


def test_item2_quantity_relabelled(bundle: Path) -> None:
    def relabel(d: dict[str, Any]) -> None:
        for e in d["arrays"]:
            if e["origin"].get("name") == "let":
                e["origin"]["quantity"] = "let_t"

    _mutate(bundle, relabel)
    with pytest.raises(ResultFormatError, match="semantic origin"):
        load(bundle)


def test_item2_unit_and_axes_changes_refused(bundle: Path) -> None:
    _mutate(
        bundle,
        lambda d: d["arrays"][0].update(units="MeV" if d["arrays"][0]["units"] != "MeV" else "Gy"),
    )
    with pytest.raises(ResultFormatError):
        load(bundle)


def test_item2_unrequested_quantity_added(bundle: Path) -> None:
    m = _members(bundle)
    extra = {}
    for base in results.QUANTITY_ARRAYS:
        extra[f"grid__dose__q__other__{base}.npy"] = m[f"grid__dose__q__let__{base}.npy"]
    _put_members(bundle, extra)

    def add(d: dict[str, Any]) -> None:
        new = []
        for e in d["arrays"]:
            if e["origin"].get("name") == "let":
                c = json.loads(json.dumps(e))
                c["name"] = e["name"].replace("__q__let__", "__q__other__")
                c["origin"].update(name="other", index=3)
                new.append(c)
        d["arrays"] = sorted(d["arrays"] + new, key=lambda x: x["name"])

    _mutate(bundle, add)
    with pytest.raises(ResultFormatError, match="inventory differs"):
        load(bundle)


def test_item2_geometry_shape_in_effective_changed(bundle: Path) -> None:
    _mutate(bundle, lambda d: d["effective"]["geometry"].update(shape=[3, 2, 5]))
    with pytest.raises(ResultFormatError):
        load(bundle)


@pytest.mark.parametrize(
    "fn",
    [
        lambda d: d["run"].update(seed="x"),
        lambda d: d["run"].update(production="yes"),
        lambda d: d["statistics"]["grids"][0].update(valid="yes"),
        lambda d: d["statistics"]["grids"][0].update(shape=[1, "a"]),
        lambda d: d["accounting"].update(energy_balance={"x": 1}),
        lambda d: d["accounting"]["energy_balance"].update(initial_mev="a"),
        lambda d: d["statistics"].update(batch_arrays={"status": "persisted", "reason": "x"}),
        lambda d: d["input_manifest"].update(config_sha256="zz"),
        lambda d: d["effective_extra"].update(materials=5),
        lambda d: d["effective_extra"].update(profile_axis=7),
        lambda d: d["units"].update(extra="x"),
        lambda d: d["provenance"].update(code=[]),
        lambda d: d["statistics"].update(normalisation=5),
    ],
)
def test_item3_malformed_fields_are_refused(bundle: Path, fn: Any) -> None:
    _mutate(bundle, fn)
    with pytest.raises(ResultFormatError):
        load(bundle)


def test_item4_user_stopping_object_round_trip(tmp_path: Path) -> None:
    from ionmc.configio import to_dict, validate_dict_structure
    from ionmc.physics.stopping import BetheStoppingSource

    base = _config()
    cfg = dataclasses.replace(
        base, physics=PhysicsOptions(nuclear=False, stopping=BetheStoppingSource())
    )
    d = to_dict(cfg)
    assert d["physics"]["stopping_model"] is None and "stopping" in d["physics"]
    validate_dict_structure(d)
    doc = load(save(make_result(cfg=cfg), tmp_path / "b")).doc
    assert doc["requested"]["physics"]["stopping"]["reconstructible"] is False


def _voxel_config() -> SimulationConfig:
    base = _config()
    geom = VoxelGeometry(
        (-3.0, -2.0, 0.0), (2.0, 2.0, 2.0), SHAPE, (WATER,), np.zeros(SHAPE, dtype=np.int32)
    )
    return dataclasses.replace(base, geometry=geom)


def test_item5_stale_config_digest(bundle: Path) -> None:
    _mutate(bundle, lambda d: d["requested"]["run"].update(allow_invalid_result=True))
    with pytest.raises(ResultFormatError, match="config_sha256"):
        load(bundle)


def test_item5_reference_manifest_binding(tmp_path: Path) -> None:
    b = save(make_result(cfg=_voxel_config()), tmp_path / "b")
    pristine = _doc(b)
    mi = pristine["requested"]["geometry"]["material_index"]
    assert mi["status"] == "not-reconstructible"
    # api manifest with a different content digest than the requested record
    _mutate(b, lambda d: d["input_manifest"]["files"][0].update(content_sha256="1" * 64))
    with pytest.raises(ResultFormatError, match="digests differ"):
        load(b)
    _put_doc(b, pristine)

    def as_toml(d: dict[str, Any], ref_sha: str) -> None:
        content = d["requested"]["geometry"]["material_index"]["content_sha256"]
        d["requested"]["geometry"]["material_index"] = {
            "ref": "material_index", "sha256": ref_sha, "content_sha256": content,
        }  # fmt: skip
        m = d["input_manifest"]
        m.update(
            source="toml", base_dir="/x", config_file={"path": "/x/c.toml", "sha256": "3" * 64}
        )
        f = m["files"][0]
        f.update(path="/x/m.npy", sha256="2" * 64, bytes=123)
        m["config_sha256"] = sha256_canonical(d["requested"])

    _mutate(b, lambda d: as_toml(d, "2" * 64))
    assert load(b).doc["input_manifest"]["source"] == "toml"  # consistent: accepted
    _mutate(b, lambda d: as_toml(d, "4" * 64))
    with pytest.raises(ResultFormatError, match="disagrees with the manifest entry"):
        load(b)
    _mutate(b, lambda d: d["input_manifest"]["files"][0].update(role="density"))
    with pytest.raises(ResultFormatError):
        load(b)


# ---- review REVIEW-89227ea7 regressions -------------------------------------------------------


def test_item1_validly_typed_semantic_mismatches_are_refused(bundle: Path) -> None:
    pristine = _doc(bundle)
    cases = [
        lambda d: d["effective"].update(physics="corrupt"),
        lambda d: d["effective"]["physics"].pop("nuclear"),
        lambda d: d["effective"]["physics"].update(surprise=1),
        lambda d: d["provenance"].update(code={}),
        lambda d: d["provenance"]["code"].update(git_dirty="no"),
        lambda d: d["run"].update(seed=12),
        lambda d: d["effective"].update(seed=12),
        lambda d: d["run"].update(n_histories=16),
        lambda d: d["statistics"]["grids"][0].update(shape=[3, 2, 5]),
        lambda d: d["statistics"]["grids"][0].update(
            energy_total_mev={"value": 1.5, "units": "MeV"}
        ),
        lambda d: d["effective_extra"].update(profile_axis=0),
        lambda d: d["provenance"]["models"].update(stopping_model_id_effective="other"),
    ]
    for fn in cases:
        _put_doc(bundle, pristine)
        _mutate(bundle, fn)
        with pytest.raises(ResultFormatError):
            load(bundle)


def test_item1_array_record_types() -> None:
    from ionmc.configio import ConfigError, to_dict, validate_dict_structure

    good = to_dict(_config())
    geo = {**good["geometry"]}
    vox = {"type": "voxel", "origin_mm": [0, 0, 0], "spacing_mm": [1, 1, 1], "shape": [1, 1, 1],
           "materials": ["water"], "z_exit_mm": None, "density": None}  # fmt: skip
    rec = {"status": "not-reconstructible", "reason": "r", "content_sha256": "a" * 64}
    validate_dict_structure({**good, "geometry": {**vox, "material_index": rec}})
    for bad in (
        {**rec, "content_sha256": "xyz"},
        {**rec, "reason": 5},
        {**rec, "status": "other"},
    ):
        with pytest.raises(ConfigError):
            validate_dict_structure({**good, "geometry": {**vox, "material_index": bad}})
    assert geo


def test_item2_automatic_channel_removed(bundle: Path) -> None:
    gone = {
        e["name"] + ".npy"
        for e in _doc(bundle)["arrays"]
        if e["origin"].get("name") == "edep_excluded_from_let"
    }
    assert len(gone) == 5
    m = {k: v for k, v in _members(bundle).items() if k not in gone}
    _write_zip(bundle, sorted(m.items()))
    _mutate(
        bundle,
        lambda d: d.update(arrays=[e for e in d["arrays"] if e["name"] + ".npy" not in gone]),
    )
    with pytest.raises(ResultFormatError, match="inventory differs"):
        load(bundle)


def test_item2_grid_origin_changed_in_effective(bundle: Path) -> None:
    _mutate(bundle, lambda d: d["effective"]["scoring"][0].update(origin_mm=[-3.0, -2.0, 1.0]))
    with pytest.raises(ResultFormatError, match="effective.scoring disagrees"):
        load(bundle)


@pytest.fixture(scope="module")
def lookup_bundle(tmp_path_factory: pytest.TempPathFactory) -> Path:
    from ionmc.lookup import LookupTable

    data = Path(__file__).resolve().parent / "data" / "v3_006"
    table = LookupTable.from_file(
        Path(__file__).resolve().parents[1] / "data" / "synthetic_lookup.json"
    )
    cfg = SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 30.0),
        VoxelGeometry(
            (-6.0, -6.0, 0.0), (4.0, 4.0, 4.0), (3, 3, 8), (WATER, PMMA),
            np.load(data / "slab_material_index.npy"), np.load(data / "slab_density.npy"),
        ),
        (ScoringGrid((-6.0, -6.0, 0.0), (4.0, 4.0, 4.0), (3, 3, 8)),),
        PhysicsOptions(nuclear=False, stopping_model="bethe-ionmc-v1"),
        RunOptions("warp-cpu", "float64", 13, 8, n_batches=4),
        tallies=(
            TallyRequest("let", "dose", "let_d"),
            TallyRequest("weighted", "dose", "lookup_sum", lookup=table.name),
        ),
        lookups=(table,),
    )  # fmt: skip
    r = Simulation(cfg).run()
    return save(r, tmp_path_factory.mktemp("lk") / "b")


def test_item2_lookup_tally_round_trip_and_unit_change(lookup_bundle: Path, tmp_path: Path) -> None:
    got = load(lookup_bundle)
    names = {
        e["origin"].get("name") for e in got.doc["arrays"] if e["origin"]["kind"] == "quantity"
    }
    assert {"weighted", "let", "scoring_pieces"} <= names
    assert [e["role"] for e in got.doc["input_manifest"]["files"]]
    b = tmp_path / "copy"
    shutil.copytree(lookup_bundle, b)

    def unit(d: dict[str, Any]) -> None:
        for e in d["arrays"]:
            if e["origin"].get("name") == "weighted" and e["origin"]["array"] == "mean":
                e["units"] = "Gy"

    _mutate(b, unit)
    with pytest.raises(ResultFormatError, match="semantic origin"):
        load(b)


def test_item3_role_and_item_digest_mutations(lookup_bundle: Path, tmp_path: Path) -> None:
    b = tmp_path / "copy"
    shutil.copytree(lookup_bundle, b)
    pristine = _doc(b)

    def swap_roles(d: dict[str, Any]) -> None:
        f = d["input_manifest"]["files"]
        f[0]["role"], f[1]["role"] = f[1]["role"], f[0]["role"]

    def one_role(d: dict[str, Any]) -> None:
        d["input_manifest"]["files"][0]["role"] = "lookup:other"

    def item(d: dict[str, Any]) -> None:
        d["input_manifest"]["files"][0]["item_sha256"] = "5" * 64

    for fn, match in (
        (swap_roles, "digests differ|roles"),
        (one_role, "roles"),
        (item, "item_sha256"),
    ):
        _put_doc(b, pristine)
        _mutate(b, fn)
        with pytest.raises(ResultFormatError, match=match):
            load(b)


def test_fixture_summaries_conform_to_the_declared_schema() -> None:
    """The schema is declared by hand from config.py; every observed summary must conform (the
    schema is not derived from the fixture)."""
    from ionmc._result_shapes import check_effective

    fx = json.loads(
        (Path(__file__).parent / "data/v3_006/summaries-95937553.json").read_text(encoding="utf-8")
    )
    for name, c in fx["cases"].items():
        check_effective(c["summary"])
        assert name


def test_declared_schema_accepts_a_numeric_z_exit_and_rejects_a_string() -> None:
    from ionmc._result_shapes import SchemaError, check

    spec = {"z": ("or", "num", "null")}
    geo = ("rec", spec, frozenset())
    check({"z": 30.0}, geo)
    check({"z": None}, geo)
    check({"z": 30}, geo)
    with pytest.raises(SchemaError):
        check({"z": "30"}, geo)
    with pytest.raises(SchemaError):
        check({"z": True}, geo)


# ---- review REVIEW-b79a1dc7 regressions (fix round 4) -----------------------------------------


def _voxel_result(z_exit: float | None) -> Result:
    from ionmc.geometry import VoxelGeometry as VG

    geo = VG(
        (-3.0, -2.0, 0.0), (2.0, 2.0, 10.0), SHAPE, (WATER,), np.zeros(SHAPE, dtype=np.int32),
        z_exit_mm=z_exit,
    )  # fmt: skip
    cfg = SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 20.0),
        geo,
        (ScoringGrid((-3.0, -2.0, 0.0), (2.0, 2.0, 2.0), SHAPE, "dose"),),
        PhysicsOptions(nuclear=False, stopping_model="bethe-ionmc-v1"),
        RunOptions("python", "float64", 11, 8, n_batches=NB),
    )
    return make_result(cfg=cfg)


def test_numeric_z_exit_round_trip_and_binding(tmp_path: Path) -> None:
    b = save(_voxel_result(30.0), tmp_path / "b")
    got = load(b)
    assert got.doc["effective"]["geometry"]["z_exit_mm"] == 30.0
    assert got.doc["requested"]["geometry"]["z_exit_mm"] == 30.0
    for bad in (31.0, None, "30"):
        _mutate(b, lambda d, v=bad: d["effective"]["geometry"].update(z_exit_mm=v))
        with pytest.raises(ResultFormatError):
            load(b)
        _mutate(b, lambda d: d["effective"]["geometry"].update(z_exit_mm=30.0))
    load(b)


@pytest.mark.parametrize(
    ("key", "value"),
    [("backend", "warp-cpu"), ("precision", "float32"), ("seed", 12), ("n_histories", 16),
     ("n_batches", 2), ("cpu_workers", 2), ("chunk_histories", 7)],
)  # fmt: skip
def test_run_and_effective_changed_together_but_not_requested(
    bundle: Path, key: str, value: Any
) -> None:
    def fn(d: dict[str, Any]) -> None:
        d["run"][key] = value
        d["effective"][key] = value
        if key in ("n_histories", "n_batches"):
            d["run"]["histories_per_batch"] = d["run"]["n_histories"] // d["run"]["n_batches"]

    _mutate(bundle, fn)
    with pytest.raises(ResultFormatError, match="requested|multiple|Not|histories"):
        load(bundle)


@pytest.mark.parametrize(
    "fn",
    [
        lambda e: e["physics"].update(max_step_mm=5.0),  # differs from requested
        lambda e: e["physics"].update(max_step_mm=-1.0),
        lambda e: e["physics"].update(max_step_mm="1"),
        lambda e: e["physics"].update(max_rejection_attempts=3.5),
        lambda e: e["physics"].update(nuclear=True),
        lambda e: e["geometry"].update(origin_mm=[0.0, 0.0]),
        lambda e: e["geometry"].update(shape=[1.0, 1.0, 1.0]),
        lambda e: e["geometry"].update(spacing_mm=[1.0, 1.0, 1.0]),
        lambda e: e["geometry"].update(materials=["unobtainium"]),
        lambda e: e["geometry"].update(z_exit_mm=4.0),
        lambda e: e["source"].update(kinetic_energy_mev=21.0),
        lambda e: e["source"].update(z=2),
        lambda e: e.update(unit_direction=[0.0, 0.6, 0.8]),
        lambda e: e["scoring"][0].update(shape=[3, 2, 5]),
        lambda e: e["scoring"][0].update(spacing_mm=[2.0, 2.0, 3.0]),
        lambda e: e["rng"].update(key="x"),
        lambda e: e["stopping"].update(model_id_requested="other"),
        lambda e: e.update(max_steps_origin="user"),
    ],
)
def test_nested_effective_fields_are_bound_and_typed(bundle: Path, fn: Any) -> None:
    _mutate(bundle, lambda d: fn(d["effective"]))
    with pytest.raises(ResultFormatError):
        load(bundle)


@pytest.mark.parametrize(
    "fn",
    [
        lambda d: d["statistics"]["grids"][0].update(valid=False),
        lambda d: d["statistics"]["grids"][0].update(low_batch_count=False),
        lambda d: d["statistics"].update(normalisation="per batch"),
        lambda d: d["statistics"]["grids"][0]["energy_total_mev"].update(units="Gy"),
        lambda d: d["accounting"]["energy_balance"]["in_grid_mev"].__setitem__(0, 2.0),
        lambda d: d["accounting"]["energy_balance"].update(closure_residual_mev=0.5),
        lambda d: d["accounting"]["energy_balance"]["outside_mev"].append(1.0),
        lambda d: d["provenance"]["models"].update(mcs_model="other"),
        lambda d: d["provenance"]["data"].update(table_stopping_source="other"),
        lambda d: d["provenance"]["code"].update(git_dirty="no"),
        lambda d: d["effective_extra"]["materials"][0].update(name="other"),
    ],
)
def test_statistics_accounting_and_provenance_are_bound(bundle: Path, fn: Any) -> None:
    _mutate(bundle, fn)
    with pytest.raises(ResultFormatError):
        load(bundle)


def test_catalogue_units_and_effective_tallies_forged_together(bundle: Path) -> None:
    def fn(d: dict[str, Any]) -> None:
        d["effective"]["tallies"]["quantities"][0]["units"] = "Gy"
        for e in d["arrays"]:
            o = e["origin"]
            if o["kind"] == "quantity" and o.get("name") == "let" and o["array"] == "mean":
                e["units"] = "Gy"

    _mutate(bundle, fn)
    with pytest.raises(ResultFormatError, match="plan"):
        load(bundle)


def test_local_class_energy_channel_removed_with_catalogue_and_archive(bundle: Path) -> None:
    pre = "grid__dose__q__edep_excluded_from_let__"

    def fn(d: dict[str, Any]) -> None:
        t = d["effective"]["tallies"]
        t["channels"] = [
            c for c in t["channels"] if not (c["kind"] == "E" and c["class_mask"] == 2)
        ]
        t["n_channels"] = len(t["channels"])
        t["quantities"] = [q for q in t["quantities"] if q["name"] != "edep_excluded_from_let"]
        d["arrays"] = [e for e in d["arrays"] if not e["name"].startswith(pre)]

    _mutate(bundle, fn)
    m = {k: v for k, v in _members(bundle).items() if not k.startswith(pre)}
    _write_zip(bundle, sorted(m.items()))  # also recomputes the recorded file hashes
    with pytest.raises(ResultFormatError, match="plan"):
        load(bundle)


def test_derived_channel_plan_equals_the_engine_plan_for_every_quantity() -> None:
    """The reader-side plan builder agrees with ``compile_channels`` (EffectiveConfig.summary())."""
    from ionmc.configio import to_dict
    from ionmc.lookup import LookupTable

    table = LookupTable.from_file(
        Path(__file__).resolve().parents[1] / "data" / "synthetic_lookup.json"
    )
    names = ("edep", "dose", "fluence", "let_t", "let_d", "let_d_eps")
    tallies = [TallyRequest(f"t_{q}", "dose", q) for q in names]
    tallies += [
        TallyRequest("t_sum", "dose", "lookup_sum", lookup=table.name),
        TallyRequest("t_avg", "dose", "lookup_dose_avg", lookup=table.name, species=("proton",)),
        TallyRequest(
            "t_spec", "dose", "fluence_spectrum", energy_edges_mev_per_u=(1.0, 2.0, 4.0, 8.0, 16.0)
        ),
    ]
    cfg = SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 30.0),
        BoxPhantom((-3.0, -2.0, 0.0), (6.0, 4.0, 8.0), WATER),
        (ScoringGrid((-3.0, -2.0, 0.0), (2.0, 2.0, 2.0), SHAPE, "dose"),),
        PhysicsOptions(nuclear=False, stopping_model="bethe-ionmc-v1"),
        RunOptions("warp-cpu", "float64", 5, 8, n_batches=NB),
        tallies=tuple(tallies),
        lookups=(table,),
    )
    eff = validate(cfg)
    doc = {"requested": to_dict(cfg), "provenance": {"data": {"lookups": [table.provenance()]}}}
    ch, qs, total, n_res = results._derive_channel_plan(doc)
    plan = eff.summary()["tallies"]
    assert [{k: v for k, v in c.items() if k not in ("k", "bound_per_history")} for c in
            plan["channels"]] == ch  # fmt: skip
    assert plan["quantities"] == qs
    assert (plan["total_size"], plan["n_residual_columns"]) == (total, n_res)


# ---- amendment A-1: re-derivation of effective, accounting identities, array invariants -------


def _copy(b: Path, tmp_path: Path, name: str = "c") -> Path:
    dst = tmp_path / name
    shutil.copytree(b, dst)
    return dst


def test_a1_default_box_and_voxel_round_trips_are_re_derived(tmp_path: Path) -> None:
    for i, r in enumerate((make_result(tally=True), make_result(cfg=_voxel_config()),
                           _voxel_result(30.0))):  # fmt: skip
        got = load(save(r, tmp_path / f"b{i}"))
        assert got.effective_binding == "bundle-local" and len(got.binding_limits) >= 3


def _bump_max_steps(d: dict[str, Any]) -> None:
    d["effective"]["max_steps"] += 1


def _bump_scattering(d: dict[str, Any]) -> None:
    d["effective"]["scattering_length_g_cm2"][0] *= 1.01


def _bump_table_n_e(d: dict[str, Any]) -> None:
    d["effective"]["tables"]["n_e"] += 1


@pytest.mark.parametrize("fn", [_bump_scattering, _bump_table_n_e])
def test_a2_table_dependent_effective_fields_are_a_declared_tier1_limit(
    bundle: Path, fn: Any, tmp_path: Path, monkeypatch: Any
) -> None:
    """Tier 1 never opens the cache: a changed table-derived field is accepted (hash-protected
    only, named in binding_limits) and is found by verify_effective (tier 2)."""
    _mutate(bundle, fn)
    monkeypatch.setenv("IONMC_CACHE_DIR", str(tmp_path / "empty"))  # tier 1 must not need it
    got = load(bundle)
    assert got.effective_binding == "bundle-local" and "effective.tables.*" in got.binding_limits[0]


def test_fix8_max_steps_is_bound_through_the_piece_count_channel_bound(bundle: Path) -> None:
    _mutate(bundle, _bump_max_steps)  # the N bound is max_steps 2 scoring_pieces
    with pytest.raises(ResultFormatError, match="recomputed plan"):
        load(bundle)


def _phys(key: str, value: Any) -> Any:
    return lambda d: d["effective"]["physics"].__setitem__(key, value)


@pytest.mark.parametrize(
    "fn",
    [
        _phys("e_cut_effective_mev", 3.0),
        _phys("end_of_range_model_effective", "continuous_residual_v1"),
        _phys("ls_effective", "eff"),
        _phys("scattering_E_s_mev", 20.0),
        lambda d: d["effective"].__setitem__("production", True),
        lambda d: d["effective"]["stopping"].__setitem__("species_transported", ["alpha"]),
        lambda d: d["effective_extra"]["materials"][0].__setitem__("density_g_cm3", 2.0),
        lambda d: d["effective_extra"]["materials"][0]["mass_fractions"].__setitem__("H", 0.5),
        lambda d: d["effective_extra"]["materials"][0].__setitem__("source", "other"),
        lambda d: d["effective_extra"]["materials"][0].__setitem__("I_eV_effective", 99.0),
        lambda d: d["effective_extra"].__setitem__("materials", []),
    ],
)
def test_a2_tier1_reconciliation_and_material_negatives(bundle: Path, fn: Any) -> None:
    _mutate(bundle, fn)
    with pytest.raises(ResultFormatError):
        load(bundle)


def test_a2_lookup_plan_is_bound_to_the_persisted_f_max(
    lookup_bundle: Path, tmp_path: Path
) -> None:
    got = load(lookup_bundle)
    assert any("f_max" in x for x in got.binding_limits)
    d = _doc(lookup_bundle)
    ch = d["effective"]["tallies"]["channels"]
    fe = next(i for i, c in enumerate(ch) if c["kind"] == "FE")
    other = next(i for i, c in enumerate(ch) if c["kind"] not in ("FE", "N"))
    assert d["effective_extra"]["lookup_plan"][0]["channel"] == fe

    def edit(tag: str, fn: Any) -> Path:
        b = _copy(lookup_bundle, tmp_path, tag)
        _mutate(b, fn)
        return b

    tl = lambda d: d["effective"]["tallies"]  # noqa: E731
    cases = {
        "k": lambda d: tl(d)["channels"][fe].update(k=tl(d)["channels"][fe]["k"] + 1),
        "bound": lambda d: tl(d)["channels"][fe].update(
            bound_per_history=tl(d)["channels"][fe]["bound_per_history"] * 1.01),
        "fmax_plan": lambda d: d["effective_extra"]["lookup_plan"][0].update(
            f_max=d["effective_extra"]["lookup_plan"][0]["f_max"] * 4.0),
        "bounds_fmax": lambda d: tl(d)["bounds"].update(f_max=tl(d)["bounds"]["f_max"] * 2.0),
        "plan_empty": lambda d: d["effective_extra"].__setitem__("lookup_plan", []),
        "plan_name": lambda d: d["effective_extra"]["lookup_plan"][0].update(lookup="other"),
        }  # fmt: skip
    for tag, fn in cases.items():
        with pytest.raises(ResultFormatError):
            load(edit(tag, fn))
    # table-derived numerical bounds are a declared tier-1 limit (found by verify_effective)
    other_k = lambda d: tl(d)["channels"][other].update(k=tl(d)["channels"][other]["k"] + 1)  # noqa: E731
    bl = lambda d: tl(d)["bounds"].update(B_L_mm=tl(d)["bounds"]["B_L_mm"] * 1.01)  # noqa: E731
    # the quantum exponent of a channel that feeds a quantity is bound through rounding_bound
    # (fix round 7); the table-derived bound of the plan stays a declared limit
    with pytest.raises(ResultFormatError, match="quantum exponent"):
        load(edit("other_k", other_k))
    with pytest.raises(ResultFormatError, match="recomputed plan|bounds"):
        load(edit("bl", bl))  # B_L is bound to the channel bounds and the derived bounds (fix 8)
    assert load(lookup_bundle).binding_limits
    d["effective_extra"]["lookup_plan"]  # the unchanged bundle loads (done above)


def test_a1_stopping_identity_block_is_validated(bundle: Path, monkeypatch: Any) -> None:
    def refuse(block: Any, transported: Any = None) -> None:
        raise ValueError("planted")

    monkeypatch.setattr(results, "validate_identity_block", refuse)
    with pytest.raises(ResultFormatError, match="identity block"):
        load(bundle)


def _eb(key: str, fn: Any) -> Any:
    return lambda d: d["accounting"]["energy_balance"].__setitem__(
        key, fn(d["accounting"]["energy_balance"][key])
    )


@pytest.mark.parametrize(
    ("fn", "match"),
    [
        (_eb("closure_residual_mev", lambda v: v + 1.0), "closure_residual"),
        (_eb("initial_mev", lambda v: v + 1.0), "closure_residual"),
        (_eb("escaped_mev", lambda v: v + 1.0), "closure_residual"),
        (_eb("outside_mev", lambda v: [v[0] + 1.0]), "deposits"),
        (_eb("quantization_mev", lambda v: [0.5]), "deposits"),
        (_eb("step_deposit_mev", lambda v: v + 1.0), "closure_residual|deposits"),
        (_eb("in_grid_mev", lambda v: [v[0] * 1.5]), "deposits|in_grid"),
        (_eb("in_grid_mev", lambda v: [-v[0]]), "negative|deposits"),
    ],
)
def test_a1_accounting_identities_negative_controls(bundle: Path, fn: Any, match: str) -> None:
    _mutate(bundle, fn)
    with pytest.raises(ResultFormatError, match=match):
        load(bundle)


def _alter(b: Path, member: str, fn: Any) -> None:
    """Change one array in place and fix the catalogue digest and the file hashes (an adversary)."""
    m = _members(b)
    arr = np.load(io.BytesIO(m[member + ".npy"]), allow_pickle=False).copy()
    fn(arr)
    m[member + ".npy"] = _npy(arr)
    _write_zip(b, sorted(m.items()))
    d = _doc(b)
    for e in d["arrays"]:
        if e["name"] == member:
            e["content_sha256"] = array_content_sha256(arr)
    _put_doc(b, d)
    _refile(b)


def _set(idx: tuple[int, ...], value: Any) -> Any:
    def f(a: np.ndarray[Any, Any]) -> None:
        a[idx] = value

    return f


def _scale(idx: Any, factor: float) -> Any:
    def f(a: np.ndarray[Any, Any]) -> None:
        a[idx] *= factor

    return f


_G, _Q = "grid__dose__", "grid__dose__q__let__"


@pytest.mark.parametrize(
    ("member", "fn"),
    [
        (_G + "dose_gy", _scale((1, 1, 1), 1.01)),
        (_G + "dose_std_gy", _scale((1, 1, 1), 1.01)),
        (_G + "relative_uncertainty", _scale((1, 1, 1), 1.1)),
        (_G + "energy_std_mev", _set((1, 1, 1), -1.0)),
        (_G + "energy_std_mev", _set((1, 1, 1), np.inf)),
        (_G + "mass_g", _scale((1, 1, 1), 2.0)),
        (_G + "defined_mask", _set((1, 1, 1), False)),
        (_G + "n_nonzero_batches", _set((1, 1, 1), 0)),
        (_G + "n_nonzero_batches", _set((1, 1, 1), 99)),
        (_G + "sparse_mask", _set((1, 1, 1), True)),
        (_G + "batch_energy_mev", _scale(Ellipsis, 1.01)),
        (_G + "batch_energy_mev", _scale((0,), 1.5)),
        (_Q + "std", _set((1, 1, 1), -1.0)),
        (_Q + "n_nonzero", _set((1, 1, 1), 1)),
        (_Q + "mean", _set((0, 0, 0), 0.0)),
        (_Q + "rounding_bound", _set((0, 0, 0), 0.0)),
        (_Q + "rounding_bound", _set((1, 1, 1), -1.0)),
        ("grid__dose__q__scoring_pieces__defined_mask", _set((1, 1, 1), False)),
    ],
)
def test_a1_array_invariants_negative_controls(bundle: Path, member: str, fn: Any) -> None:
    _alter(bundle, member, fn)
    with pytest.raises(ResultFormatError):
        load(bundle)


def test_a1_array_invariants_positive_control_unchanged_arrays(bundle: Path) -> None:
    _alter(bundle, _G + "dose_gy", lambda a: None)  # re-written identically: still accepted
    assert load(bundle).effective_binding == "bundle-local"


# ---- amendment A-2: ratio masks, nuclear destinations, per-quantity batch arrays --------------


def test_a2_ratio_defined_mask_is_n_nonzero_threshold_in_both_directions(bundle: Path) -> None:
    # a bin that is undefined although its denominator is nonzero in every batch (all NaN, so
    # only the threshold rule can refuse it)
    for base in ("mean", "std", "rounding_bound"):
        _alter(bundle, _Q + base, _set((1, 1, 1), np.nan))
    _alter(bundle, _Q + "defined_mask", _set((1, 1, 1), False))
    with pytest.raises(ResultFormatError, match="n_nonzero >="):
        load(bundle)


def _nuclear_doc(**over: float) -> dict[str, Any]:
    nuc = {
        k: 0.0
        for k in (
            "nuclear_local",
            "nuclear_alpha_local",
            "nuclear_escaped_neutron",
            "nuclear_escaped_gamma",
            "nuclear_binding",
            "nuclear_imbalance",
        )
    }
    nuc.update(nuclear_local=2.0, nuclear_alpha_local=1.0, **over)
    eb = {
        "initial_mev": 5.0,
        "step_deposit_mev": 3.0,
        "cutoff_mev": 0.0,
        "escaped_mev": 0.0,
        "truncated_mev": 0.0,
        "unaccounted_mev": 0.0,
        "closure_residual_mev": 0.0,
        "in_grid_mev": [5.0],
        "outside_mev": [0.0],
        "quantization_mev": [],
        "nuclear_mev": nuc,
    }
    return {"accounting": {"energy_balance": eb}, "effective": {"scoring": [{}]}}  # fmt: skip


def test_a2_unsigned_nuclear_destinations_and_alpha_local_part_of_local() -> None:
    results._check_energy_structure(_nuclear_doc())
    for bad, match in (({"nuclear_alpha_local": 3.0}, "alpha_local"),):
        d = _nuclear_doc()
        d["accounting"]["energy_balance"]["nuclear_mev"].update(bad)
        with pytest.raises(ResultFormatError, match=match):
            results._check_energy_structure(d)
    d = _nuclear_doc()
    nm = d["accounting"]["energy_balance"]["nuclear_mev"]
    nm.update(nuclear_escaped_gamma=-1.0, nuclear_local=3.0)  # keeps the closure, breaks the sign
    with pytest.raises(ResultFormatError, match="non-negative"):
        results._check_energy_structure(d)


def test_fix7_nuclear_binding_and_imbalance_are_signed() -> None:
    d = _nuclear_doc(nuclear_binding=-1.0, nuclear_imbalance=1.0)  # exothermic: negative binding
    results._check_energy_structure(d)
    d = _nuclear_doc(nuclear_binding=-0.5, nuclear_imbalance=0.25)  # breaks only the closure
    with pytest.raises(ResultFormatError, match="closure"):
        results._check_energy_structure(d)


def test_a2_quantity_batch_arrays_are_persisted_and_bind_means(
    lookup_bundle: Path, tmp_path: Path
) -> None:
    doc = _doc(lookup_bundle)
    assert doc["statistics"]["quantity_batch_arrays"] == {"status": "persisted"}
    names = [e["name"] for e in doc["arrays"] if e["origin"].get("array", "").startswith("batch_")]
    assert any(n.endswith("__batch_den") for n in names)  # the let_d ratio
    assert any(n.endswith("__batch_num") for n in names)
    assert not any("table" in n for n in names)
    ratio_q = next(e for e in doc["arrays"] if e["origin"].get("quantity") == "let_d")
    stem = ratio_q["name"].rsplit("__", 1)[0]
    lin_q = next(e for e in doc["arrays"] if e["origin"].get("quantity") == "lookup_sum")
    lstem = lin_q["name"].rsplit("__", 1)[0]
    for tag, member, fn in (
        ("num", stem + "__batch_num", _scale((0,), 1.5)),
        ("den", stem + "__batch_den", _scale((1,), 0.5)),
        ("lin", lstem + "__batch_num", _scale(Ellipsis, 1.01)),
    ):
        b = _copy(lookup_bundle, tmp_path, tag)
        _alter(b, member, fn)
        with pytest.raises(ResultFormatError, match="batch arrays|disagrees"):
            load(b)


def test_a2_persist_batches_false_records_the_status(tmp_path: Path) -> None:
    r = Simulation(_ion_config()).run()
    got = load(save(r, tmp_path / "b", persist_batches=False))
    assert got.doc["statistics"]["quantity_batch_arrays"]["status"] == "not-persisted"
    assert any("not recomputed" in x for x in got.binding_limits)


# ---- fix round 7: frozen schema-1.0 rules, RNG identity, stopping pins, rounding bound ----------


def test_fix7_frozen_snapshot_equals_live_definitions() -> None:
    import math as _m

    from ionmc import _result_schema_1_0 as fz
    from ionmc.config import MAX_REJECTION_ATTEMPTS, U01_MAPPING
    from ionmc.materials import ELEMENTS, MATERIALS
    from ionmc.physics.em import E_S_MEV
    from ionmc.physics.projectiles import PROJECTILES
    from ionmc.scoring import MEV_PER_G_TO_GY, min_defined_batches
    from ionmc.simulation import ratio_rounding_bound
    from ionmc.transport import channels as ch
    from ionmc.transport.tally import NUCLEAR_TALLY_NAMES

    assert fz.PROJECTILES == {n: (p.z, p.a, p.mass_mev) for n, p in PROJECTILES.items()}
    assert fz.ELEMENTS == {s: (e.Z, e.A_g_mol, e.I_eV) for s, e in ELEMENTS.items()}
    assert set(fz.MATERIALS) == set(MATERIALS)
    for n, m in MATERIALS.items():
        assert fz.material_definition({"name": n}) == results._material_def(m)
    inline = {"name": "x", "density_g_cm3": 1.1, "mass_fractions": {"H": 0.1, "O": 0.9}}
    from ionmc.configio import material_from_record

    want = results._material_def(material_from_record(inline))
    got = fz.material_definition(inline)
    assert got["I_eV_effective"] == pytest.approx(want["I_eV_effective"], rel=1e-14)
    assert {k: v for k, v in got.items() if k != "I_eV_effective"} == {
        k: v for k, v in want.items() if k != "I_eV_effective"
    }
    assert (fz.SCATTERING_E_S_MEV, fz.MAX_REJECTION_ATTEMPTS) == (E_S_MEV, MAX_REJECTION_ATTEMPTS)
    assert fz.MEV_PER_G_TO_GY == MEV_PER_G_TO_GY and fz.U01_MAPPING == U01_MAPPING
    assert (fz.K_MAX, fz.CAPACITY_EXPONENT, fz.E_QUANTUM_EXPONENT) == (
        ch.K_MAX, ch.CAPACITY_EXPONENT, ch.E_QUANTUM_EXPONENT)  # fmt: skip
    assert (fz.FE_F_MIN_EXPONENT, fz.FE_F_MAX_EXPONENT, fz.GEN_MAX) == (
        ch.FE_F_MIN_EXPONENT, ch.FE_F_MAX_EXPONENT, ch.GEN_MAX)  # fmt: skip
    assert (fz.CLASS_STEP, fz.CLASS_LOCAL) == (ch.CLASS_STEP, ch.CLASS_LOCAL)
    assert fz.NUCLEAR_TALLY_NAMES == NUCLEAR_TALLY_NAMES
    from ionmc.physics.stopping import RANGE_CONSTRUCTION

    assert fz.NUCLEAR_MAX_PARTICLES == ch.MAX_PARTICLES
    assert fz.RANGE_CONSTRUCTION == RANGE_CONSTRUCTION
    for f in (1e-9, 0.37, 1.0, 5.5, 1e6):
        assert fz.fe_exponent(f) == ch.fe_exponent(f)
    for hpb, bound in ((1, 1.0), (1000, 3.7), (10**6, 1e-3), (50, 2.0**30)):
        assert fz.adaptive_exponent(hpb, bound) == ch.adaptive_exponent(hpb, bound)
    assert all(fz.min_defined_batches(n) == min_defined_batches(n) for n in range(1, 40))
    x, y = np.array([1.0, 2.0, 0.0]), np.array([4.0, 0.1, 0.0])
    a, b, r = np.array([0.1, 0.1, 0.1]), np.array([0.2, 0.2, 0.2]), np.array([0.25, 20.0, 0.5])
    assert np.array_equal(
        fz.ratio_rounding_bound(x, a, y, b, r), ratio_rounding_bound(x, a, y, b, r)
    )
    cfg = _config()
    assert fz.rng_identity(cfg.run.precision) == validate(cfg).rng
    assert _m.isfinite(fz.MEV_PER_G_TO_GY)


def test_fix7_tier1_does_not_use_the_live_registries(bundle: Path, monkeypatch: Any) -> None:
    import ast
    import inspect

    import ionmc.config
    import ionmc.materials
    import ionmc.physics.em
    import ionmc.physics.projectiles
    import ionmc.scoring
    import ionmc.transport.channels as chm

    forbidden = {
        "E_S_MEV", "MAX_REJECTION_ATTEMPTS", "PROJECTILES", "DEUTERON", "MEV_PER_G_TO_GY",
        "min_defined_batches", "fe_exponent", "NUCLEAR_TALLY_NAMES", "material_from_record",
        "_material_def", "CLASS_STEP", "CLASS_LOCAL", "GEN_MAX", "MATERIALS", "ELEMENTS",
    }  # fmt: skip
    tier1 = (
        results._check_run_binding, results._check_effective_binding, results._check_materials,
        results._check_reconciliation, results._check_energy_structure,
        results._derive_channel_plan, results._check_channel_plan, results._check_lookup_plan,
        results._check_quantity_batches, results._check_rounding_bound,
        results._quantity_channels, results._check_array_invariants,
        results._mass_from_geometry, results._check_invariants, results._check_stopping_pins,
    )  # fmt: skip
    for fn in tier1:
        names = {
            n.id if isinstance(n, ast.Name) else n.attr
            for n in ast.walk(ast.parse(inspect.getsource(fn).lstrip()))
            if isinstance(n, ast.Name | ast.Attribute)
        }
        # a frozen.<name> attribute is the allowed route
        src = inspect.getsource(fn)
        live = {n for n in forbidden & names if f"frozen.{n}" not in src or n in ("MATERIALS",)}
        assert not live, (fn.__name__, live)

    # drift of every live constant and helper does not change the verdict of load()
    def boom(*a: Any, **k: Any) -> Any:
        raise AssertionError("live helper used by tier 1")

    monkeypatch.setattr(results, "fe_exponent", boom, raising=False)
    monkeypatch.setattr(results, "CLASS_STEP", 99, raising=False)
    monkeypatch.setattr(results, "CLASS_LOCAL", 98, raising=False)
    monkeypatch.setattr(ionmc.physics.em, "E_S_MEV", 99.0)
    monkeypatch.setattr(ionmc.config, "MAX_REJECTION_ATTEMPTS", 1)
    monkeypatch.setattr(ionmc.scoring, "MEV_PER_G_TO_GY", 1.0)
    monkeypatch.setattr(ionmc.scoring, "min_defined_batches", boom)
    monkeypatch.setattr(chm, "fe_exponent", boom)
    monkeypatch.setitem(
        ionmc.physics.projectiles.PROJECTILES, "proton",
        dataclasses.replace(PROTON, mass_mev=PROTON.mass_mev * 1.001, name="proton"),
    )  # fmt: skip
    monkeypatch.setitem(
        ionmc.materials.MATERIALS, "water", dataclasses.replace(WATER, density_g_cm3=2.0)
    )
    assert load(bundle).doc["schema"] == "ionmc.result"


def test_fix7_rng_identity_is_the_frozen_record(bundle: Path) -> None:
    def both(d: dict[str, Any]) -> None:  # an adversary that keeps run.rng == effective.rng
        d["run"]["rng"]["generator"] = "philox4x32-7"
        d["effective"]["rng"]["generator"] = "philox4x32-7"

    _mutate(bundle, both)
    with pytest.raises(ResultFormatError, match="RNG identity"):
        load(bundle)


def _pin_doc(**over: Any) -> dict[str, Any]:
    from ionmc.physics.stopping_models import (
        NIST_CONTENT_SHA256,
        NIST_MODEL_ID,
        STOPPING_MODELS,
        _canonical,
        spec_sha256,
    )

    spec = STOPPING_MODELS[NIST_MODEL_ID]
    pin = next(_canonical(p) for p in spec.data if p.dataset_id == "nist-pstar-water-2005")
    st: dict[str, Any] = {
        "model_id_requested": NIST_MODEL_ID,
        "model_id_effective": NIST_MODEL_ID,
        "spec_sha256": spec_sha256(spec),
        "version": spec.version,
        "kind": spec.kind,
        "interpolation": spec.interpolation,
        "range_start": spec.range_start,
        "decision": spec.decision,
        "domain": _canonical(spec.domain),
        "options": {**_canonical(spec.options), "class": "NistStarCombinedSource",
                    "ls_requested": "off", "ls_effective": "off"},
        "selected_by": "stopping_model",
        "range_construction": "exact-loglog-quadrature-v1",
        "content_sha256": NIST_CONTENT_SHA256["proton"],
        "domain_resolved": {"proton": {"e_min_per_u": 1.0, "e_max_per_u": 500.0}},
        "species_transported": ["proton"],
        "data_declared": [_canonical(p) for p in spec.data],
        "data_used": [{"pin": pin, "role_in_this_run": "construction",
                       "table": ["proton", "water"], "cache_dir": None}],
    }  # fmt: skip
    st.update(over)
    return {
        "effective": {"stopping": st, "geometry": {"materials": ["water"]},
                      "source": {"projectile": "proton", "z": 1, "a": 1,
                                 "kinetic_energy_mev": 100.0, "energy_sigma_mev": 0.0},
                      "physics": {"e_cut_effective_mev": 1.0}},
        "effective_extra": {"materials": [{
            "name": "water", "density_g_cm3": 1.0, "mass_fractions": {"H": 0.111894, "O": 0.888106},
            "I_eV_requested": 78.0, "I_eV_effective": 78.0}]},
        "requested": {"physics": {"stopping_model": NIST_MODEL_ID}, "tallies": []},
    }  # fmt: skip


def test_fix7_stopping_pins_are_registry_backed() -> None:
    results._check_stopping_pins(_pin_doc())
    with pytest.raises(ResultFormatError, match="data_used"):  # pin removed
        results._check_stopping_pins(_pin_doc(data_used=[]))
    d = _pin_doc()  # rewritten pin with mirrored provenance (declared and used agree)
    st = d["effective"]["stopping"]
    for p in (st["data_used"][0]["pin"], *st["data_declared"]):
        if p["dataset_id"] == "nist-pstar-water-2005":
            p["archive_sha256"] = "0" * 64
    with pytest.raises(ResultFormatError, match="data_declared|pin"):
        results._check_stopping_pins(d)
    d = _pin_doc()  # a different source projectile needs different pins
    d["effective"]["source"].update(projectile="alpha", z=2)
    with pytest.raises(ResultFormatError, match="data_used|domain_resolved|content_sha256"):
        results._check_stopping_pins(d)


def test_fix8_stopping_identity_is_complete_and_data_used_is_an_exact_multiset() -> None:
    results._check_stopping_pins(_pin_doc())
    for key, bad in (("version", "9"), ("kind", "analytic"), ("interpolation", "x"),
                     ("range_start", "x"), ("decision", "x"), ("domain", {}),
                     ("options", {"class": "C"}), ("selected_by", "stopping-object"),
                     ("range_construction", "other"), ("content_sha256", "a" * 64),
                     ("domain_resolved", None)):  # fmt: skip
        with pytest.raises(ResultFormatError, match=key):
            results._check_stopping_pins(_pin_doc(**{key: bad}))
    d = _pin_doc()  # a duplicated use is not the multiset of the run
    st = d["effective"]["stopping"]
    st["data_used"] = [st["data_used"][0], dict(st["data_used"][0])]
    with pytest.raises(ResultFormatError, match="data_used"):
        results._check_stopping_pins(d)
    d = _pin_doc()  # a projectile the model does not support, a disallowed material
    d["effective"]["source"].update(projectile="carbon12", z=6, a=12)
    d["effective"]["stopping"]["species_transported"] = ["carbon12"]
    with pytest.raises(ResultFormatError, match="does not support"):
        results._check_stopping_pins(d)
    d = _pin_doc()
    d["effective_extra"]["materials"][0]["density_g_cm3"] = 1.2
    with pytest.raises(ResultFormatError, match="material rule"):
        results._check_stopping_pins(d)
    d = _pin_doc()  # a forged nuclear run: transported deuterons are not a NIST projectile
    d["effective"]["stopping"]["species_transported"] = ["proton", "deuteron"]
    with pytest.raises(ResultFormatError, match="does not support"):
        results._check_stopping_pins(d)
    d = _pin_doc()  # below-window cutoff and above-window source
    d["effective"]["physics"] = {"e_cut_effective_mev": 0.5}
    d["effective"]["source"].update(kinetic_energy_mev=100.0, energy_sigma_mev=0.0)
    with pytest.raises(ResultFormatError, match="energy interval"):
        results._check_stopping_pins(d)
    d["effective"]["physics"]["e_cut_effective_mev"] = 1.0
    results._check_stopping_pins(d)
    d["effective"]["source"]["kinetic_energy_mev"] = 600.0
    with pytest.raises(ResultFormatError, match="energy interval"):
        results._check_stopping_pins(d)
    d = _pin_doc()  # a material outside the geometry
    d["effective"]["stopping"]["data_used"][0]["table"] = ["proton", "bone"]
    with pytest.raises(ResultFormatError, match="data_used"):
        results._check_stopping_pins(d)
    d = _pin_doc()  # tallies add the LET water row: a missing use is detected
    d["effective"]["geometry"]["materials"] = ["w2"]
    d["effective_extra"]["materials"].append(dict(d["effective_extra"]["materials"][0], name="w2"))
    d["effective"]["stopping"]["data_used"][0]["table"] = ["proton", "w2"]
    d["requested"]["tallies"] = [{"name": "t"}]
    with pytest.raises(ResultFormatError, match="data_used"):
        results._check_stopping_pins(d)


def test_fix7_user_stopping_object_binds_what_it_can_and_records_the_limit(
    tmp_path: Path,
) -> None:
    from ionmc.physics.stopping import BetheStoppingSource

    cfg = dataclasses.replace(
        _config(), physics=PhysicsOptions(nuclear=False, stopping=BetheStoppingSource())
    )
    b = save(make_result(cfg=cfg), tmp_path / "b")
    got = load(b)
    assert any("user object" in x for x in got.binding_limits)  # tier 1: the declared limit
    _mutate(b, lambda d: d["effective"]["stopping"].update(species_transported=["alpha"]))
    with pytest.raises(ResultFormatError):  # every bindable field is still bound
        load(b)


def test_fix7_rounding_bound_is_recomputed_from_the_plan_and_piece_counts(
    lookup_bundle: Path, tmp_path: Path
) -> None:
    assert load(lookup_bundle).doc["statistics"]["quantity_batch_arrays"]["status"] == "persisted"
    for q in ("let", "weighted", "scoring_pieces"):
        b = _copy(lookup_bundle, tmp_path, q)
        _alter(b, f"grid__dose__q__{q}__rounding_bound", _scale(Ellipsis, 1.5))
        with pytest.raises(ResultFormatError, match="rounding_bound"):
            load(b)
    b = _copy(lookup_bundle, tmp_path, "k")  # the persisted quantum exponent is part of the plan
    _mutate(b, lambda d: d["effective"]["tallies"]["channels"][0].update(
        k=d["effective"]["tallies"]["channels"][0]["k"] + 1))  # fmt: skip
    with pytest.raises(ResultFormatError):
        load(b)
    b = _copy(lookup_bundle, tmp_path, "cnt")  # piece counts feed the bound
    _alter(b, "grid__dose__q__scoring_pieces__batch_num", _scale(Ellipsis, 2.0))
    with pytest.raises(ResultFormatError):
        load(b)


@pytest.mark.parametrize("which", ["batch_num", "batch_den"])
def test_fix7_batch_sums_must_be_non_negative(
    lookup_bundle: Path, tmp_path: Path, which: str
) -> None:
    b = _copy(lookup_bundle, tmp_path, which)
    _alter(b, f"grid__dose__q__let__{which}", _set((0, 0, 0, 0), -1.0))
    with pytest.raises(ResultFormatError, match="non-negative"):
        load(b)


def test_fix7_stopping_spec_sha256_is_registry_backed() -> None:
    results._check_stopping_pins(_pin_doc())
    with pytest.raises(ResultFormatError, match="spec_sha256"):
        results._check_stopping_pins(_pin_doc(spec_sha256="0" * 64))


def test_fix8_every_channel_exponent_and_bound_is_recomputed(
    lookup_bundle: Path, tmp_path: Path
) -> None:
    d = _doc(lookup_bundle)
    chans = d["effective"]["tallies"]["channels"]
    assert {c["kind"] for c in chans} >= {"E", "N"}
    e_idx = [i for i, c in enumerate(chans) if c["kind"] == "E"]
    assert all(chans[i]["k"] == 30 for i in e_idx)
    # an energy-channel exponent 30 -> 29 with the dependent rounding-bound arrays and hashes
    # coherently rewritten (rounding half-quantum doubles): only the plan recomputation can object
    b = _copy(lookup_bundle, tmp_path, "e29")
    for q in d["effective"]["tallies"]["quantities"]:
        if q["kind"] == "linear" and chans[q["numerator"]]["kind"] == "E":
            _alter(
                b,
                f"grid__{d['effective']['scoring'][chans[q['numerator']]['grid']]['name']}"
                f"__q__{q['name']}__rounding_bound",
                _scale(Ellipsis, 2.0),
            )
    _mutate(b, lambda dd: [dd["effective"]["tallies"]["channels"][i].update(k=29) for i in e_idx])
    with pytest.raises(ResultFormatError, match="quantum exponent"):
        load(b)
    # every other kind: the bound per history and the exponent follow the persisted bounds
    for tag, fn in {
        "n_bound": lambda dd: [
            c.update(bound_per_history=c["bound_per_history"] * 2)
            for c in dd["effective"]["tallies"]["channels"]
            if c["kind"] == "N"
        ],
        "e_bound": lambda dd: [
            c.update(bound_per_history=c["bound_per_history"] * 2)
            for c in dd["effective"]["tallies"]["channels"]
            if c["kind"] == "E"
        ],
        "n_k": lambda dd: [
            c.update(k=1) for c in dd["effective"]["tallies"]["channels"] if c["kind"] == "N"
        ],
        "b_es": lambda dd: dd["effective"]["tallies"]["bounds"].update(
            B_ES=dd["effective"]["tallies"]["bounds"]["B_ES"] * 2
        ),
    }.items():
        bb = _copy(lookup_bundle, tmp_path, tag)
        _mutate(bb, fn)
        with pytest.raises(ResultFormatError):
            load(bb)


def test_fix8_voxel_material_digests_ignore_live_registry_drift(
    tmp_path: Path, monkeypatch: Any
) -> None:
    from ionmc import materials

    b = save(_voxel_result(30.0), tmp_path / "b")
    assert load(b).doc["requested"]["materials"] == [{"name": "water"}]
    drifted = dataclasses.replace(materials.MATERIALS["water"], density_g_cm3=1.5)
    monkeypatch.setitem(materials.MATERIALS, "water", drifted)
    from ionmc import configio

    monkeypatch.setitem(configio.MATERIALS, "water", drifted)
    assert load(b).doc["input_manifest"]["files"]  # frozen definitions, not the live registry
    _mutate(b, lambda d: d["input_manifest"]["files"][0].update(item_sha256="0" * 64))
    with pytest.raises(ResultFormatError, match="item_sha256"):
        load(b)


def test_fix8_user_stopping_object_is_bound_to_the_requested_record(tmp_path: Path) -> None:
    from ionmc.physics.stopping import BetheStoppingSource

    cfg = dataclasses.replace(
        _config(), physics=PhysicsOptions(nuclear=False, stopping=BetheStoppingSource())
    )
    b = save(make_result(cfg=cfg), tmp_path / "b")
    load(b)
    for tag, fn in {
        "class": lambda d: d["effective"]["stopping"]["options"].update(**{"class": "Other"}),
        "id": lambda d: d["effective"]["stopping"].update(model_id_requested="other"),
        "sha": lambda d: d["effective"]["stopping"].update(spec_sha256="0" * 64),
    }.items():
        bb = _copy(b, tmp_path, tag)
        _mutate(bb, fn)
        with pytest.raises(ResultFormatError):
            load(bb)


def test_fix9_stopping_options_class_selection_and_content_are_bound(tmp_path: Path) -> None:
    b = save(make_result(tally=True), tmp_path / "b")
    load(b)
    mirror = "a" * 64

    def content(d: dict[str, Any]) -> None:  # mirrored into provenance: still not an analytic's
        d["effective"]["stopping"]["content_sha256"] = mirror
        d["provenance"]["data"]["stopping_content_sha256"] = mirror

    for tag, fn in {
        "selected": lambda d: d["effective"]["stopping"].update(selected_by="stopping-object"),
        "class": lambda d: d["effective"]["stopping"]["options"].update(
            **{"class": "NistStarCombinedSource"}
        ),
        "range": lambda d: d["effective"]["stopping"].update(range_construction="other"),
        "content": content,
    }.items():
        bb = _copy(b, tmp_path, tag)
        _mutate(bb, fn)
        with pytest.raises(ResultFormatError):
            load(bb)


def test_fix9_history_energy_and_particle_bounds_are_derived_from_the_request(
    bundle: Path, tmp_path: Path
) -> None:
    d = _doc(bundle)
    tl = d["effective"]["tallies"]
    assert not any("nuclear" in x and "max_particles" in x for x in load(bundle).binding_limits)

    def e_side(dd: dict[str, Any]) -> None:  # every E bound, E_hi_mev and B_ES scaled together
        t = dd["effective"]["tallies"]
        for c in t["channels"]:
            if c["kind"] == "E":
                c["bound_per_history"] *= 2.0
        t["bounds"]["E_hi_mev"] *= 2.0
        t["bounds"]["B_ES"] *= 2.0

    bb = _copy(bundle, tmp_path, "e")
    _mutate(bb, e_side)
    with pytest.raises(ResultFormatError, match="E_hi_mev"):
        load(bb)
    # sigma != 0: E_hi must be E0 + 6 sigma (the model's table cap is far above): every E bound,
    # E_hi_mev and B_ES moved coherently to another in-range value is refused
    d = _doc(bundle)
    d["effective"]["source"]["energy_sigma_mev"] = 1.0
    tl = d["effective"]["tallies"]

    def set_e(e: float) -> None:
        for c in tl["channels"]:
            if c["kind"] == "E":
                c["bound_per_history"] = e
        tl["bounds"]["B_ES"] = tl["bounds"]["S_w_max_mev_per_mm"] * e
        tl["bounds"]["E_hi_mev"] = e

    e0 = d["effective"]["source"]["kinetic_energy_mev"]
    set_e(e0 + 3.0)
    with pytest.raises(ResultFormatError, match="E_hi_mev"):
        results._check_channel_numerics(d, tl["channels"])
    # above the table maximum (500 MeV/u, within the window tolerance): E_hi is exactly the cap
    d["effective"]["source"].update(kinetic_energy_mev=499.9999999)
    d["effective"]["source"]["energy_sigma_mev"] = (500.0 - 499.9999999 + 2.5e-10) / 6.0
    set_e(500.0)
    results._check_channel_numerics(d, tl["channels"])
    set_e(
        d["effective"]["source"]["kinetic_energy_mev"]
        + 6.0 * d["effective"]["source"]["energy_sigma_mev"]
    )  # the uncapped value, coherent everywhere
    with pytest.raises(ResultFormatError, match="E_hi_mev"):
        results._check_channel_numerics(d, tl["channels"])
    d["effective"]["source"].update(kinetic_energy_mev=e0, energy_sigma_mev=1.0)
    set_e(e0 + 6.0)
    # nuclear run: E_bound = max(E_hi, persisted history bound), max_particles = 32 (frozen)
    set_e(e0 + 6.0)
    d["effective"]["nuclear"] = {"history_energy_bound_mev": 2.0 * e0}
    nuc = e0 * 2.0
    for c in tl["channels"]:
        if c["kind"] == "E":
            c["bound_per_history"] = nuc
    tl["bounds"].update(
        E_bound_mev=nuc, max_particles=32.0, B_ES=tl["bounds"]["S_w_max_mev_per_mm"] * nuc
    )
    for c in tl["channels"]:
        if c["kind"] == "N":
            c["bound_per_history"] *= 32.0
    results._check_channel_numerics(d, tl["channels"])  # coherent and derivable: accepted
    for c in tl["channels"]:  # E side moved coherently away from the persisted history bound
        if c["kind"] == "E":
            c["bound_per_history"] = 3.0 * e0
    tl["bounds"].update(E_bound_mev=3.0 * e0, B_ES=tl["bounds"]["S_w_max_mev_per_mm"] * 3.0 * e0)
    with pytest.raises(ResultFormatError, match="particle|energy bound"):
        results._check_channel_numerics(d, tl["channels"])
    for c in tl["channels"]:  # N side: max_particles and every N bound moved together
        if c["kind"] == "E":
            c["bound_per_history"] = nuc
        if c["kind"] == "N":
            c["bound_per_history"] *= 2.0
    tl["bounds"].update(
        E_bound_mev=nuc, max_particles=64.0, B_ES=tl["bounds"]["S_w_max_mev_per_mm"] * nuc
    )
    with pytest.raises(ResultFormatError, match="particle|recomputed plan"):
        results._check_channel_numerics(d, tl["channels"])


def test_fix10_tabulated_bundle_content_digest_is_registry_bound(tmp_path: Path) -> None:
    from ionmc.physics.stopping_models import NIST_MODEL_ID

    cfg = dataclasses.replace(
        _config(), physics=PhysicsOptions(nuclear=False, stopping_model=NIST_MODEL_ID)
    )
    b = save(make_result(cfg=cfg), tmp_path / "b")
    load(b)
    mirror = "a" * 64

    def content(d: dict[str, Any]) -> None:  # mirrored into provenance: still not the registry's
        d["effective"]["stopping"]["content_sha256"] = mirror
        d["provenance"]["data"]["stopping_content_sha256"] = mirror

    _mutate(b, content)
    with pytest.raises(ResultFormatError, match="NIST digest"):
        load(b)
