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
import struct
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ionmc import results
from ionmc.config import PhysicsOptions, RunOptions, SimulationConfig, validate
from ionmc.geometry import BoxPhantom
from ionmc.materials import WATER
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
        BoxPhantom((-3.0, -2.0, 0.0), (6.0, 4.0, 8.0), WATER),
        tuple(ScoringGrid((-3.0, -2.0, 0.0), (2.0, 2.0, 2.0), SHAPE, n) for n in names),
        PhysicsOptions(nuclear=False, stopping_model="bethe-ionmc-v1"),
        RunOptions("python", "float64", 11, 8, n_batches=NB),
        tallies=(TallyRequest("let", names[0], "let_d"),) if tally else (),
    )


def _grid(name: str, shape: tuple[int, ...] = SHAPE, seed: int = 0) -> GridResult:
    rng = np.random.default_rng(seed)  # synthetic arrays only
    en = rng.uniform(0.5, 2.0, shape)
    en[0, 0, 0] = 0.0
    std = 0.1 * en
    mass = np.full(shape, 0.01)
    mass[-1, -1, -1] = 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(en > 0, std / en, np.nan)
        dose = np.where(mass > 0, en * 3.0, np.nan)
        dstd = np.where(mass > 0, std * 3.0, np.nan)
    mask = (en > 0) & (mass > 0)
    return GridResult(
        name, en, std, dose, dstd, rel, mask, mass,
        np.full(shape, 3, dtype=np.int32), mask & False, False,
        rng.uniform(0.5, 2.0, (NB, *shape)), True,
    )  # fmt: skip


def _quantity(shape: tuple[int, ...] = SHAPE) -> QuantityResult:
    mean = np.linspace(1.0, 2.0, int(np.prod(shape))).reshape(shape)
    mask = np.ones(shape, dtype=bool)
    mean[0, 0, 0], mask[0, 0, 0] = np.nan, False
    return QuantityResult(
        "let", "let_d", "ratio", "keV/um", "synthetic ratio", mean, 0.1 * mean, mask,
        np.full(shape, 4, dtype=np.int32), np.where(mask, 1e-9, np.nan), (3, 3),
    )  # fmt: skip


def _pieces() -> dict[str, QuantityResult]:
    out = {}
    for n in ("scoring_pieces", "scoring_pieces_local"):
        one = np.ones(SHAPE)
        out[n] = QuantityResult(
            n, "count", "linear", "1", "synthetic piece count", one, 0.0 * one,
            np.ones(SHAPE, dtype=bool), np.full(SHAPE, 4, dtype=np.int32), 0.0 * one, (3,),
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
    for i, n in enumerate(names):
        g = _grid(n, seed=i)
        if tally and i == 0:
            g = GridResult(*[getattr(g, f) for f in g.__dataclass_fields__ if f != "quantities"],
                           quantities={"let": _quantity(), **_pieces()})  # fmt: skip
        grids.append(g)
    bal = EnergyBalance(
        8.0, 7.0, 0.5, 0.25, 0.0, 0.0, tuple(1.0 for _ in names), (6.0,) * len(names)
    )
    return Result(
        True, cfg, eff, "python", "cpu", "float64", 11, {"name": "philox4x32-10"}, 8, NB,
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
        g.name, z, z, d0, d0, np.full(SHAPE, np.nan), mask & False, g.mass_g, g.n_nonzero_batches,
        g.sparse_mask, False, g.batch_energy_mev, True,
    )  # fmt: skip
    doc = load(save(Result(**{**r.__dict__, "grids": (none,)}), tmp_path / "c")).doc
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
    if kind == "depth":
        return b"[" * 65 + b"]" * 65
    if kind == "numbers":
        return b"[" + b"1," * 2_100_000 + b"1]"
    return b"[" + b"true,false,null," * 700_000 + b"null]"


@pytest.mark.parametrize("kind", ["objects", "depth", "numbers", "literals"])
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


def _with_channels(lookup: int, path: int) -> Result:
    from ionmc.transport.tally import ChannelRaw

    r = make_result()
    raw = ChannelRaw(np.zeros((NB, 1), dtype=np.int64), [], lookup, path)
    return Result(**{**r.__dict__, "valid": False, "channel_raw": raw})


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
    doc = load(save(r, tmp_path / "b")).doc
    assert set(doc["accounting"]["counters"]) == set(r.counters.as_dict())
    assert doc["valid"] is True


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
    bad = Result(
        **{**r.__dict__, "valid": False, "counters": IonTransportCounters(mcs_mapping_failure=2)}
    )
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
    doc = load(save(r, tmp_path / "b")).doc
    assert "majorant_violation" in doc["accounting"]["counters"]
    assert "nuclear_mev" in doc["accounting"]["energy_balance"]


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
    _mutate(bundle, lambda d: d["requested"]["run"].update(seed=12))
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
        f.update(path="/x/m.npy", sha256="2" * 64)
        m["config_sha256"] = sha256_canonical(d["requested"])

    _mutate(b, lambda d: as_toml(d, "2" * 64))
    assert load(b).doc["input_manifest"]["source"] == "toml"  # consistent: accepted
    _mutate(b, lambda d: as_toml(d, "4" * 64))
    with pytest.raises(ResultFormatError, match="disagrees with the manifest entry"):
        load(b)
    _mutate(b, lambda d: d["input_manifest"]["files"][0].update(role="density"))
    with pytest.raises(ResultFormatError):
        load(b)
