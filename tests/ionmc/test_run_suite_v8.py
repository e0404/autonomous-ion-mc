"""V3-008 suites ``lv8`` and ``hr8`` of ``run_suite.py`` and their steps in ``steps_v8.py``: step
lists and names, seed derivation (plan section C), the hashed source set, the summarizer tags and
a tiny-scale smoke of the cheap CPU steps (CI seeds: the steps' own derivation from the lv8 base at
``--scale`` far below 1, labelled ``reduced``; nothing here is qualification evidence)."""

from __future__ import annotations

import importlib.util
import io
import json
import sys
from contextlib import redirect_stdout
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "validation" / "scripts" / "transport"
RESERVED = range(20572004, 20632004)


def _load(name: str) -> ModuleType:
    if str(SCRIPTS) not in sys.path:
        sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


RS = _load("run_suite")
V8 = _load("steps_v8")
SUM = _load("summarize")

LV8 = [
    "01-i6-lv-deterministic-end-depth",
    "02-i8-range-straggling-he150",
    "03-i8-range-straggling-c290",
    "04-i8-range-straggling-o350",
    *[f"{i + 5:02d}-i13-r80-vs-csda-{n}" for i, n in enumerate(
        ("he50", "he150", "he250", "c80", "c290", "c430", "o100", "o350", "o430"))],
    "14-t1-ion-trace-parity-c290-o350",
    "15-r1-ion-repeatability-c290",
    "16-r2-a16-trace-fixture",
]  # fmt: skip
SHARDED_LV8 = ("t12-ion-lv", "i9", "i12", "i14", "i15")
SHARD_COUNTS = {"t12-ion-lv": 4, "i9": 6, "i12": 36, "i14": 48, "i15": 20, "t12-ion-hr": 10}


def _sharded_names(rows: tuple[str, ...], start: int) -> list[str]:
    out = []
    for row in rows:
        for unit, n_sh, _n, _b, _m in V8.shard_plan(row):
            out += [f"{row}-{unit}-s{j}" for j in range(n_sh)]
        out.append(f"{row}-combine")
    return [f"{i:02d}-{n}" for i, n in enumerate(out, start)]


def test_suite_step_names_and_seed_bases() -> None:
    assert {"lv8", "hr8"} <= set(RS.SUITES)
    assert RS.full_step_names("lv8", 2) == LV8 + _sharded_names(SHARDED_LV8, 17)
    assert len(RS.full_step_names("lv8", 2)) == 16 + sum(SHARD_COUNTS[r] + 1 for r in SHARDED_LV8)
    assert RS.full_step_names("hr8", 2) == ["01-t13-ion-chunks-cuda-c290",
                                            *_sharded_names(("t12-ion-hr",), 2)]  # fmt: skip
    assert RS.DEFAULT_SEED_BASES["lv8"] == 20581004 == V8.LV8_BASE
    assert RS.DEFAULT_SEED_BASES["hr8"] == 20591004 == V8.HR8_BASE
    assert RS.DEFERRED_STEPS["lv8"] == () and RS.DEFERRED_STEPS["hr8"] == ()
    assert RS.step_timeout_s("lv8", "16-r2-a16-trace-fixture", 1500) == 3300
    for name in RS.full_step_names("lv8", 2) + RS.full_step_names("hr8", 2):
        assert SUM.expected_tag(name) not in (None, False), name  # every step has its document tag
    for row in (*SHARDED_LV8, "t12-ion-hr"):
        names = [n for n in RS.full_step_names("lv8", 2) + RS.full_step_names("hr8", 2)
                 if n.split("-", 1)[1].startswith(row + "-")]  # fmt: skip
        assert {SUM.expected_tag(n) for n in names} == {f"{row}-shard", f"{row}-combine"}
        assert sum(n.endswith("-combine") for n in names) == 1
    assert SUM.STEP_FILE.match("135-i15-combine.txt") and SUM.STEP_FILE.match("01-x.txt")
    assert not SUM.STEP_FILE.match("environment.txt")
    assert SUM.expected_tag("01-t13-ion-chunks-cuda-c290") == "t13-ion"
    assert SUM.expected_tag("14-t1-ion-trace-parity-c290-o350") == "t1-ion"
    assert SUM.expected_tag("15-r1-ion-repeatability-c290") == "r1-ion"
    assert SUM.seed_blockers(20581004, "lv8") == [] and SUM.seed_blockers(20591004, "hr8") == []
    assert SUM.seed_blockers(20601004, "lv8")
    # the third re-run base (plan section C) qualifies lv8 only
    assert SUM.seed_blockers(20572004, "lv8") == [] and SUM.seed_blockers(20572004, "hr8")
    assert SUM.seed_blockers(20621004, "lv8")  # not declared for run_suite re-runs


def test_lv8_accepts_the_third_rerun_base() -> None:
    """``--seed-base 20572004`` reaches every lv8 step unchanged (no base guard in run_suite or
    steps_v8; every seed is ``base + 100 r_index + k``)."""
    for _name, cmd, _env in RS.suite_steps("lv8", 1, 0.01, seed_base=20572004):
        assert cmd[cmd.index("--seed-base") + 1] == "20572004"
    assert V8.row_seed("i9", 0, 20572004) == 20572004 + 100 * V8.ROW_R_INDEX["i9"]


def test_argv_carries_seed_base_timeout_and_cuda_env() -> None:
    steps = RS.suite_steps("lv8", 1, 0.01, seed_base=20601004)
    for _name, cmd, _env in steps:
        assert cmd[1].endswith("steps_v8.py")
        assert cmd[cmd.index("--seed-base") + 1] == "20601004"
        assert "--timeout" in cmd
    by = {n: c for n, c, _ in steps}
    assert by["03-i8-range-straggling-c290"][2:5] == ["i8", "--point", "1"]
    assert "--scale" not in by["16-r2-a16-trace-fixture"]
    (name, cmd, env), *hr = RS.suite_steps("hr8", 1, 1.0)  # fmt: skip
    assert env == {"IONMC_REQUIRE_CUDA": "1"} and cmd[2] == "t13-ion"
    assert cmd[cmd.index("--seed-base") + 1] == "20591004"
    assert all(e == {"IONMC_REQUIRE_CUDA": "1"} and c[2] in ("shard", "combine") for _n, c, e in hr)
    assert all(c[c.index("--row") + 1] == "t12-ion-hr" for _n, c, _e in hr)


def test_shard_and_combine_wiring(tmp_path: Path) -> None:
    """Every shard writes into ``OUT/samples`` with its unit and index; the combine step of the
    row reads ``--dirs OUT [imported ...]`` (imports need the partials manifest, as lv5)."""
    out = tmp_path / "OUT"
    by = {n: c for n, c, _ in RS.suite_steps("lv8", 1, 0.5, out=out)}
    sh = by["22-i9-he150-s0"]
    assert sh[2:9] == ["shard", "--row", "i9", "--unit", "he150", "--shard", "0"]
    assert sh[sh.index("--out-dir") + 1] == str(out / "samples")
    assert sh[sh.index("--scale") + 1] == "0.5"
    comb = by["28-i9-combine"]
    assert comb[2:5] == ["combine", "--row", "i9"]
    assert comb[comb.index("--dirs") + 1] == str(out) and "--partials-manifest" not in comb
    with pytest.raises(SystemExit):
        RS.suite_steps("lv8", 1, 1.0, out=out, import_dirs=["OLD"])
    by = {n: c for n, c, _ in RS.suite_steps("lv8", 1, 1.0, out=out, import_dirs=["OLD"],
                                             partials_manifest="m.json")}  # fmt: skip
    comb = by["28-i9-combine"]
    assert comb[comb.index("--dirs") + 1 : comb.index("--dirs") + 5] == [
        str(out),
        "OLD",
        "--partials-manifest",
        "m.json",
    ]
    for row, count in SHARD_COUNTS.items():  # k = position in shard_plan; all < 100, distinct
        plan = V8.shard_plan(row)
        ks = [V8.shard_k(row, u, j) for u, n_sh, _n, _b, _m in plan for j in range(n_sh)]
        assert ks == list(range(count))
        for u, n_sh, n, _b, _m in plan:
            assert n_sh * n == {"python": 4000}.get(u, 1_000_000 if "cuda" in u else 200_000)
            frozen_s = V8.planned_shard_seconds(row, u)
            assert frozen_s < 1100, (row, u, frozen_s)  # <= 1500 s step timeout with margin
    assert V8.row_seed("i9", 5) == 20581004 + 4400 + 5
    assert V8.row_seed("t12-ion-hr", 9, V8.HR8_BASE) == 20591004 + 5400 + 9
    with pytest.raises(SystemExit):
        V8.shard_k("i9", "c290", 2)
    with pytest.raises(SystemExit):
        V8.shard_k("i12", "theta-s2", 0)


def test_hashed_source_set() -> None:
    v8 = RS.source_file_list("lv8")
    assert v8 == RS.source_file_list("hr8")
    assert v8[: len(RS.SOURCE_FILES_V5)] == RS.SOURCE_FILES_V5
    for f in ("validation/plans/v3-008-acceptance.md", "decisions/0043-ion-em-transport.md",
              "validation/scripts/transport/steps_v8.py"):  # fmt: skip
        assert f in v8 and (ROOT / f).is_file()


def test_row_r_index_and_seeds() -> None:
    idx = {**V8.ROW_R_INDEX, **V8.RESERVED_R_INDEX}
    assert len(set(idx.values())) == len(idx)
    assert all(42 <= r <= 60 for r in idx.values())  # 30-41 are the td-tl-sample rows
    assert V8.RESERVED_R_INDEX == {}  # every reserved row is implemented
    assert {
        r: V8.ROW_R_INDEX[r] for r in ("i9", "i12", "i14", "i15", "t12-ion-lv", "t12-ion-hr")
    } == {  # fmt: skip
        "i9": 44,
        "i12": 45,
        "i14": 46,
        "i15": 47,
        "t12-ion-lv": 49,
        "t12-ion-hr": 54,
    }
    assert not set(idx.values()) & {V8.r_index_of(c) for c in V8.CASES}
    assert V8.row_seed("i8", 1) == 20581004 + 4200 + 1
    assert V8.row_seed("t13-ion", 1, V8.HR8_BASE) == 20591004 + 5500 + 1
    seeds = {V8.row_seed(r, k, b) for r in V8.ROW_R_INDEX for k in range(100)
             for b in (V8.LV8_BASE, V8.HR8_BASE)}  # fmt: skip
    assert all(s in RESERVED for s in seeds)
    td = {V8.seed_of(c, k) for c in V8.CASES for k in range(100)}
    assert not seeds & td
    with pytest.raises(ValueError):
        V8.row_seed("i8", 100)


def _run(argv: list[str]) -> tuple[int, dict]:  # type: ignore[type-arg]
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = V8.main(argv)
    text = buf.getvalue()
    return rc, json.loads(text.split("#JSON-BEGIN", 1)[1].split("#JSON-END", 1)[0])


def test_i8_quadrature_is_finite_and_scales_like_the_proton() -> None:
    """sigma_R/R(C 290) with corrections on is close to sigma_R/R(proton at equal velocity) times
    sqrt(m_p/M) (the I4b identity holds within 1e-3 without corrections; a few % with them)."""
    proj = V8.CARBON12
    sq, rq = V8.sigma_range_quadrature(proj, 290.0 * proj.a, 2.0 * proj.a)
    assert 0.001 < sq / rq < 0.006
    r = V8.csda_mm(proj, 290.0)
    assert abs(rq / r - 1.0) < 0.01


def test_smoke_i6_and_i8_tiny_scale() -> None:
    rc, doc = _run(["i6-lv"])
    assert rc == 0 and doc["pass"] and doc["step"] == "i6-lv" and len(doc["cases"]) == 18
    rc, doc = _run(["i8", "--point", "0", "--scale", "0.002"])
    assert doc["reduced"] and doc["histories"] == 200 and doc["seed"] == 20581004 + 4200
    assert doc["counters"] and doc["valid"] and doc["end_codes_all_cutoff"]
    assert "i8b_report_only" in doc and rc in (0, 1)  # 200 histories: the 5 % bound may fail


@pytest.mark.parametrize(("row", "unit"), [
    ("i9", "c290"), ("i12", "theta-s1"), ("i12", "sigma-s5"), ("i14", "s1-f002"),
    ("i15", "theta-v2-sh1"), ("i15", "sigma-v1-sh1"), ("t12-ion-lv", "python"),
    ("t12-ion-lv", "cpu64"),
])  # fmt: skip
def test_smoke_shard_tiny_scale(tmp_path: Path, row: str, unit: str) -> None:
    """One shard of each sharded lv8 row at <= 200 histories: sealed partial, seed, PARTIAL."""
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = V8.main(["shard", "--row", row, "--unit", unit, "--shard", "1", "--scale", "0.0004",
                      "--out-dir", str(tmp_path)])  # fmt: skip
    text = buf.getvalue()
    doc = json.loads(text.split("#JSON-BEGIN", 1)[1].split("#JSON-END", 1)[0])
    assert rc == 0 and doc["pass"] and doc["step"] == f"{row}-shard" and doc["reduced"]
    assert 20 <= doc["histories"] <= 200
    assert doc["seed"] == V8.row_seed(row, V8.shard_k(row, unit, 1))
    name = f"{row}-{unit}-s1.json"
    part = json.loads((tmp_path / name).read_text())
    assert f"PARTIAL {name} {part['content_sha256']}" in text
    assert part["valid"] is True and part["n"] == doc["histories"] and part["seed"] == doc["seed"]


def _write_synthetic(tmp_path: Path, row: str, scale: float, payload) -> None:  # type: ignore[no-untyped-def]
    import argparse

    import numpy as np

    rng = np.random.default_rng(7)
    base = V8.HR8_BASE if row == "t12-ion-hr" else V8.LV8_BASE
    a = argparse.Namespace(out_dir=str(tmp_path / "samples"), scale=scale, seed_base=base)
    with redirect_stdout(io.StringIO()):
        for unit, n_sh, n_nom, b, _m in V8.shard_plan(row):
            for j in range(n_sh):
                k = V8.shard_k(row, unit, j)
                n = V8.scaled(n_nom, scale, b)
                V8.write_partial8(a, f"{row}-{unit}-s{j}", {
                    "row": row, "unit": unit, "shard": j, "k": k, "n": n, "batches": b,
                    "seed": V8.row_seed(row, k, a.seed_base), "valid": True, "counters": {},
                    **payload(rng, unit, n, b)})  # fmt: skip


def _combine(tmp_path: Path, row: str, scale: float) -> tuple[int, dict]:  # type: ignore[type-arg]
    return _run(["combine", "--row", row, "--dirs", str(tmp_path), "--scale", str(scale)])


def _t12_payload(rng, unit, n, b):  # type: ignore[no-untyped-def]
    import numpy as np

    z = np.arange(200.0)
    idd = np.exp(-(((z - 150.0) / 20.0) ** 2)) + 0.2
    lat = np.exp(-(((np.arange(200.0) - 100.0) / 10.0) ** 2))
    noise = lambda m: m * (1.0 + 0.01 * rng.standard_normal((b, len(m))))  # noqa: E731
    return {"arrays": {"idd": noise(idd).tolist(), "lat25": noise(lat).tolist(),
                       "lat50": noise(lat).tolist(), "lat90": noise(lat).tolist()},
            "scalars": {"total_deposit_mev": (100 + rng.standard_normal(b)).tolist(),
                        "r80_mm": (160 + 0.1 * rng.standard_normal(b)).tolist(),
                        "sigma_lat_05R_mm": (2 + 0.01 * rng.standard_normal(b)).tolist()},
            "n_batches": b, "histories_per_batch": n // b, "precision": "float32"
            if "cuda32" in unit else "float64", "backend": "synthetic"}  # fmt: skip


def test_combine_t12_ion_hr_synthetic_and_fail_closed(tmp_path: Path) -> None:
    """The combine step pools every shard, compares the ``T12_PAIRS`` (f32-vs-f64 probe included)
    and fails closed on a missing, altered or invalid shard (exit != 0 through SystemExit)."""
    _write_synthetic(tmp_path, "t12-ion-hr", 0.001, _t12_payload)
    rc, doc = _combine(tmp_path, "t12-ion-hr", 0.001)
    assert doc["step"] == "t12-ion-hr-combine" and doc["reduced"]
    assert set(doc["t12_ion_hr"]) == {f"{x}_vs_{y}" for x, y in V8.T12_PAIRS["t12-ion-hr"]}
    assert "c290-cuda32_vs_c290-cuda64" in doc["t12_ion_hr"]
    assert doc["attestation"]["partials"] and rc == (0 if doc["pass"] else 1)
    f = tmp_path / "samples" / "t12-ion-hr-o350-cpu64-s1.json"
    good = f.read_text()
    f.write_text(good.replace('"valid": true', '"valid": false'))
    with pytest.raises(SystemExit, match="content_sha256"):
        _combine(tmp_path, "t12-ion-hr", 0.001)
    f.unlink()
    with pytest.raises(SystemExit, match="found 0 copies"):
        _combine(tmp_path, "t12-ion-hr", 0.001)
    f.write_text(good)
    with pytest.raises(SystemExit, match="scale"):  # bound to the run's scale
        _combine(tmp_path, "t12-ion-hr", 0.002)


def test_combine_i14_synthetic_r80_and_invalid_shard(tmp_path: Path) -> None:
    """I14 combine on identical synthetic IDDs: rebatched to 40 batches, Delta R80 = 0, every
    comparison against the reference present; a shard recorded invalid stops the combine."""
    import numpy as np

    z = np.arange(220.0)
    idd = np.where(z < 163, 1.0 + z / 300.0, np.exp(-(((z - 163) / 2.0) ** 2)))

    def payload(rng, unit, n, b):  # type: ignore[no-untyped-def]
        return {"idd_batches": (idd * (1 + 0.01 * rng.standard_normal((b, len(z))))).tolist(),
                "csda_range_mm": 163.0}  # fmt: skip

    _write_synthetic(tmp_path, "i14", 0.001, payload)
    rc, doc = _combine(tmp_path, "i14", 0.001)
    body = doc["i14"]
    assert body["reference"] == "s1-f002" and len(body["comparisons"]) == 5
    assert body["max_pairwise_delta_r80_mm"] < 0.05 and rc == (0 if doc["pass"] else 1)
    assert V8.rebatch(np.arange(8.0).reshape(4, 2), 2).tolist() == [[1.0, 2.0], [5.0, 6.0]]
    import argparse

    a = argparse.Namespace(out_dir=str(tmp_path / "samples"), scale=0.001, seed_base=V8.LV8_BASE)
    with redirect_stdout(io.StringIO()):  # an invalid shard, correctly sealed
        V8.write_partial8(a, "i14-s1-f002-s1", {
            "row": "i14", "unit": "s1-f002", "shard": 1, "k": 47, "n": 100, "batches": 20,
            "seed": V8.row_seed("i14", 47), "valid": False, "counters": {}, **payload(
                np.random.default_rng(1), "", 100, 20)})  # fmt: skip
    with pytest.raises(SystemExit, match="validity"):
        _combine(tmp_path, "i14", 0.001)


def test_ion_path_carries_the_material_scattering_length() -> None:
    """``steps_v8.ion_path`` (A2 of I9 and the lateral grid, theta^2 quadrature of I12 / I15) is
    ``ionmc.transport.ion_mcs_checks.IonPath`` with the transport's 1/X_S of water, not the former
    placeholder 0.0123 cm2/g."""
    from ionmc.materials import WATER
    from ionmc.physics.projectiles import CARBON12
    from ionmc.physics.scattering import inverse_scattering_length_cm2_per_g
    from ionmc.transport.ion_mcs_checks import IonPath

    path = V8.ion_path(CARBON12, 290.0)
    assert isinstance(path, IonPath)
    assert path.inv_rho_xs_cm2_g == inverse_scattering_length_cm2_per_g(WATER)
    assert path.inv_rho_xs_cm2_g != 0.0123
    assert path.rho_g_cm3 == WATER.density_g_cm3
