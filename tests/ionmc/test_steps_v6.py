"""Steps of row V6 / V6-sens (V3-005C C8): producers, the comparator dry run on reduced shards, the
fail-closed rules (V6_SHARDS, consumed bases, mismatched shards) and the registration."""

# ruff: noqa: E501

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

from tests.ionmc.test_v7r_infra import _load

REPO = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def s6() -> ModuleType:
    sys.path.insert(0, str(REPO / "validation" / "scripts" / "transport"))
    sys.modules.pop("steps_v6", None)
    return _load("steps_v6")


def _json_of(text: str) -> dict[str, Any]:
    return json.loads(text.split("#JSON-BEGIN", 1)[1].split("#JSON-END", 1)[0])


def test_v6_shards_fail_closed_until_14a(
    s6: ModuleType, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    assert s6.V6_SHARDS is None
    for step in ("v6-compare", "v6-cuda-compare", "v6-unsized"):
        with pytest.raises(SystemExit, match="V6_SHARDS"):
            s6.main([step, "--dirs", str(tmp_path)])
    with pytest.raises(SystemExit, match="consumed"):
        s6.main(["v6-shard", "--seed-base", "20471004", "--out-dir", str(tmp_path)])


def test_seeds_and_names(s6: ModuleType) -> None:
    s6.base.SEED_BASE = 20491004
    assert s6.seed_of(s6.V6_CUDA_R_INDEX, 2) == 20491004 + 3000 + 2
    s6.base.SEED_BASE = 20481004
    assert (
        s6.seed_of(s6.V6_R_INDEX, 5) == 20481004 + 12005
        and s6.seed_of(s6.V6_SENS_R_INDEX, 1) == 20481004 + 20001
    )
    assert s6.sens_names(2)[0] == "v6-sens-nominal-s0.json" and len(s6.sens_names(2)) == 18
    assert s6.cuda_names(3) == ("v6-cuda-s0.json", "v6-cuda-s1.json", "v6-cuda-s2.json")


def test_e0_step_record_and_cross_checks(
    s6: ModuleType, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    rc = s6.main(["v6-e0", "--out-dir", str(tmp_path), "--seed-base", "20505000"])
    doc = _json_of(capsys.readouterr().out)
    rec = json.loads((tmp_path / "e0_record.json").read_text())
    assert rec["E0_MeV"] == doc["E0_MeV"] and 176.0 <= rec["E0_MeV"] <= 178.0
    assert rec["inside_recording_band"] is True and rec["inside_expectation_band"] is False
    assert rec["finding_outside_176_178"] is False and rec["dose_table_read"] is False
    assert (
        rec["tbl_parms_match_frozen"] is True
        and rec["data"]["halo_tex_ok"]
        and rec["data"]["tarball_ok"]
    )
    assert rec["cross_check_exact_inverse"]["ok"] is True
    tr = rec["transport_table_finding"]
    assert (
        tr["rel_diff"] > 0.0
        and tr["within_1e-6"] == (tr["rel_diff"] <= 1e-6)
        and tr["within_v3_003d_u2"] is True
    )
    assert (
        tr["status"].startswith("finding (report-only)")
        and "orchestrator decision 2026-10-11" in tr["reclassified_after_value_observed"]
    )
    assert (
        tr["transport_table_range_g_cm2"] > 21.321 and 5e-7 < tr["energy_equivalent_rel"] < 1.2e-6
    )
    assert (
        rec["table_identity"]["transport_tables_sha256"]
        and "effective_transport_table" in rec["table_identity"]
    )
    assert rc == 0 and rec["pass"] is True and set(rec["beam_table"]) == set(s6.gv6().SENS_NAMES)


def test_producers_compare_and_mismatch_refusal(
    s6: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    g = s6.gv6()
    monkeypatch.setattr(s6, "V6_SHARDS", 1)
    out = tmp_path / "samples"
    common = ["--shard", "0", "--scale", "0.002", "--out-dir", str(out), "--seed-base", "20505000"]
    assert s6.main(["v6-shard", *common]) == 0
    for v in g.SENS_NAMES:
        assert s6.main(["v6-sens-shard", "--variant", v, *common]) == 0
    capsys.readouterr()
    p = json.loads((out / "v6-s0.json").read_text())
    assert (
        p["seed"] == 20505000 + 12000
        and p["r_index"] == 12
        and p["backend"] == "warp-cpu"
        and p["precision"] == "float64"
    )
    assert (
        p["grid"]["kind"] == "annular"
        and p["grid"]["shape"] == [210, 1, 520]
        and p["weights_digest"]
    )
    assert (
        p["beam"]["E0_MeV"] == s6.beam_table()["nominal"]["E0_MeV"] and p["beam"]["sigma_E_MeV"] > 0
    )
    assert p["sphere"] == {"diameter_mm": 4.0, "subsamples": 16} and len(p["dose_batches"]) == 20
    assert p["valid"] is True and p["outside_over_initial"] < 1e-3
    sens = {v: json.loads((out / f"v6-sens-{v}-s0.json").read_text()) for v in g.SENS_NAMES}
    assert {d["seed"] for d in sens.values()} == {20505000 + 20000}  # common random numbers
    assert sens["rho-lo"]["beam"]["rho_xtheta"] != sens["nominal"]["beam"]["rho_xtheta"]
    assert sens["sx-hi"]["weights_digest"] == p["weights_digest"]
    rc = s6.main(
        [
            "v6-compare",
            "--dirs",
            str(tmp_path),
            "--scale",
            "0.002",
            "--seed-base",
            "20505000",
            "--out-dir",
            str(out),
        ]
    )
    doc = _json_of(capsys.readouterr().out)
    assert (
        doc["row"] == "V6"
        and doc["lineage"] == g.LINEAGE
        and doc["beam"]["roles"]["tbl:parms"] == g.BEAM_ROLE
    )
    assert doc["data"]["n_points"] == 295 and doc["data"]["n_missing"] == 165
    assert doc["data"]["sha256_halo_tex"] == g.HALO_TEX_SHA256 and doc["status"] in (
        "gating",
        "report-only",
    )
    assert doc["r0"]["n_points"] > 0 and "points" in doc["halo"] and "variants" in doc["v6_sens"]
    assert set(doc["v6_sens"]["variants"]) == set(g.VARIANTS) and set(
        doc["v6_sens"]["r0_variants_report"]
    ) == set(g.R0_VARIANTS)
    assert (
        doc["beam"]["A1_cm_rad"] == pytest.approx(-9.770e-3, rel=2e-4)
        and doc["monitor_uncertainty_note"]
    )
    assert (
        doc["n_histories"] == 2000 and rc in (0, 1) and doc["weights_digest"] == p["weights_digest"]
    )
    # mismatched shard: a resealed partial with a different beam is refused
    path = out / "v6-s0.json"
    bad = json.loads(path.read_text())
    bad["beam"]["sigma_x_mm"] *= 1.01
    del bad["content_sha256"]
    bad["content_sha256"] = s6.v5.content_digest(bad)
    path.write_text(json.dumps(bad, sort_keys=True))
    with pytest.raises(SystemExit, match="beam differ"):
        s6.main(
            [
                "v6-compare",
                "--dirs",
                str(tmp_path),
                "--scale",
                "0.002",
                "--seed-base",
                "20505000",
                "--out-dir",
                str(out),
            ]
        )
    # an altered digest-less edit is refused by the seal
    bad["beam"]["sigma_x_mm"] = 1.0
    path.write_text(json.dumps(bad, sort_keys=True))
    with pytest.raises(SystemExit, match="content_sha256"):
        s6.main(
            [
                "v6-compare",
                "--dirs",
                str(tmp_path),
                "--scale",
                "0.002",
                "--seed-base",
                "20505000",
                "--out-dir",
                str(out),
            ]
        )


def test_producers_do_not_read_the_dose_table_and_reuse_the_beam_builder(s6: ModuleType) -> None:
    cfg = s6.v6_config("rho-hi", 2000, 1, "warp-cpu", "float64")
    assert (
        cfg.source.angular_sigma_mrad == 31.55
        and cfg.source.x_theta_correlation == pytest.approx(-0.6247, abs=1e-4)
    )
    assert cfg.physics.nuclear and cfg.physics.elastic and cfg.scoring[0].kind == "annular"
    cuda = s6.v6_config("nominal", 2000, 1, "warp-cuda", "float32")
    assert cuda.run.chunk_histories == 2**14 and cuda.run.precision == "float32"
    assert math.isclose(cfg.source.kinetic_energy_mev, 177.84, abs_tol=0.01)


def test_constants_match_run_suite(s6: ModuleType) -> None:
    rs = _load("run_suite")
    assert s6.V6_SHARDS == rs.V6_SHARDS and tuple(s6.gv6().SENS_NAMES) == tuple(rs.V6_SENS_VARIANTS)
    assert (
        s6.HR5C_SEED_BASE == rs.HR5C_QUALIFICATION_SEED_BASE
        and s6.LV5C_SEED_BASE == rs.V5C_QUALIFICATION_SEED_BASE
    )


# ---- review fix: the sens-free CUDA path under the sealed lv5c status (end to end, reduced scale) ----
@pytest.fixture(scope="module")
def e2e(s6: ModuleType, tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """lv5c shards (V6 + nine V6-sens variants) and two emulated hr5c sample sets on warp-cpu float64 (the
    test-only backend override): ``same`` has the SAME seed as the lv5c V6 sample (identical batches), ``diff``
    another one. 2000 histories per shard."""
    mp = pytest.MonkeyPatch()
    root = tmp_path_factory.mktemp("v6e2e")
    mp.setattr(s6, "V6_SHARDS", 1)
    out: dict[str, Path] = {"root": root}
    try:
        mp.setenv("IONMC_RUN_SUITE", "lv5c")
        cpu = root / "cpu"
        common = [
            "--shard",
            "0",
            "--scale",
            "0.002",
            "--out-dir",
            str(cpu / "samples"),
            "--seed-base",
            "20505000",
        ]
        assert s6.main(["v6-shard", *common]) == 0
        for v in s6.gv6().SENS_NAMES:
            assert s6.main(["v6-sens-shard", "--variant", v, *common]) == 0
        out["cpu"] = cpu
        mp.setenv("IONMC_RUN_SUITE", "hr5c")
        mp.setenv("IONMC_V6_CUDA_TEST_BACKEND", "warp-cpu")
        for name, base_ in (
            ("same", 20514000),
            ("diff", 20515000),
        ):  # 20514000 + 3000 = 20505000 + 12000
            d = root / name
            assert s6.main(["v6-cuda-shard", "--shard", "0", "--scale", "0.002", "--out-dir", str(d / "samples"),
                            "--seed-base", str(base_)]) == 0  # fmt: skip
            out[name] = d
        yield_dirs = out
    finally:
        mp.undo()
    return yield_dirs


def _run_compare(
    s6: ModuleType,
    mp: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    dirs: Path,
    cuda: bool,
    base_: int,
    suite: str,
    sens_holds: bool | None,
    meas_from: Path,
) -> tuple[int, dict[str, Any]]:
    g = s6.gv6()
    orig = s6._data
    cpu_mean = np.array(
        json.loads((meas_from / "samples" / "v6-s0.json").read_text())["dose_batches"]
    ).mean(axis=0)

    def fake(gg: Any) -> Any:
        z, r, _m, data = orig(gg)
        m = np.where(cpu_mean > 0, cpu_mean, np.nan)
        m[0, 1] = np.nan
        return z, r, m, data

    mp.setattr(s6, "_data", fake)
    mp.setattr(s6, "V6_SHARDS", 1)
    if sens_holds is not None:
        mp.setattr(
            g, "_sens_status", lambda *a, **k: {"condition_holds": sens_holds, "variants": {}}
        )
    mp.setenv("IONMC_RUN_SUITE", suite)
    step = "v6-cuda-compare" if cuda else "v6-compare"
    rc = s6.main(
        [
            step,
            "--dirs",
            str(dirs),
            "--scale",
            "0.002",
            "--seed-base",
            str(base_),
            "--out-dir",
            str(dirs / "samples"),
        ]
    )
    return rc, _json_of(capsys.readouterr().out)


@pytest.mark.parametrize("sens_holds", [True, False])
def test_cuda_compare_inherits_lv5c_status_and_enforces_agreement(
    s6: ModuleType,
    e2e: dict[str, Path],
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
    sens_holds: bool,
) -> None:
    import shutil

    status = "gating" if sens_holds else "report-only"
    for name, agree in (("same", True), ("diff", False)):
        with monkeypatch.context() as mp:
            rc_cpu, cpu = _run_compare(
                s6, mp, capsys, e2e["cpu"], False, 20505000, "lv5c", sens_holds, e2e["cpu"]
            )
        assert cpu["status"] == status and cpu["r0"]["pass"] and cpu["halo"]["pass"]
        d = e2e[name]
        shutil.copy(
            e2e["cpu"] / "samples" / "v6-lv5c-verdict.json", d / "samples" / "v6-lv5c-verdict.json"
        )
        base_ = 20514000 if name == "same" else 20515000
        with monkeypatch.context() as mp:
            mp.setenv("IONMC_V6_CUDA_TEST_BACKEND", "warp-cpu")
            rc, doc = _run_compare(s6, mp, capsys, d, True, base_, "hr5c", sens_holds, e2e["cpu"])
        comb = doc["backend_combination"]
        # the CUDA sample is judged under the lv5c status; its gate is a boolean, never pending
        assert doc["status"] == status and comb["evaluated_status_both"] == [status, status]
        assert (
            doc["v6_sens"]["inherited_from"] == "lv5c"
            and doc["v6_sens"]["condition_holds"] is sens_holds
        )
        assert doc["halo"]["gating"] is sens_holds and doc["r0"]["pass"] in (True, False)
        if sens_holds:
            assert isinstance(doc["gate_pass"], bool) is True and comb["pass"] is not None
        if agree:
            assert (
                comb["disagreement_components"] == []
                and comb["pass"] is True
                and doc["pass"] is True
                and rc == 0
            )
        else:
            assert (
                comb["disagreement_components"]
                and comb["pass"] is False
                and doc["pass"] is False
                and rc == 1
            )
            assert any("backend disagreement" in r for r in doc["reasons"])
            if not sens_holds:  # report-only status: the disagreement is still a recorded failure
                assert doc["status"] == "report-only" and any(
                    "backend disagreement" in r for r in doc["reasons"]
                )


def test_cuda_compare_without_sealed_lv5c_verdict_fails_closed(
    s6: ModuleType,
    e2e: dict[str, Path],
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    d = e2e["root"] / "noverdict"
    import shutil

    shutil.copytree(e2e["same"], d)
    (d / "samples" / "v6-lv5c-verdict.json").unlink(missing_ok=True)
    with monkeypatch.context() as mp:
        mp.setenv("IONMC_V6_CUDA_TEST_BACKEND", "warp-cpu")
        rc, doc = _run_compare(s6, mp, capsys, d, True, 20514000, "hr5c", None, e2e["cpu"])
    assert (
        rc == 1
        and doc["pass"] is False
        and doc["status"] == "pending"
        and "backend_combination" not in doc
    )
    assert any("fail closed" in r for r in doc["reasons"])


def test_cuda_test_backend_override_refused_outside_pytest(
    s6: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("IONMC_V6_CUDA_TEST_BACKEND", "warp-cpu")
    monkeypatch.delenv("PYTEST_CURRENT_TEST")
    with pytest.raises(SystemExit, match="test-only"):
        s6.cuda_target()


def test_combine_backends_pending_cpu_and_sensfree_cuda_regression() -> None:
    """The defect of the review: a sens-free CUDA verdict (status pending, gate None) against a gating CPU one
    can never pass; evaluate(fixed_status=...) gives the CUDA sample the CPU status and a boolean gate."""
    import importlib.util as iu

    spec = iu.spec_from_file_location(
        "gv6_t", REPO / "validation" / "scripts" / "reference" / "gottschalk_v6.py"
    )
    assert spec and spec.loader
    g = iu.module_from_spec(spec)
    spec.loader.exec_module(g)
    z, r = np.array(g.DEPTHS_CM), np.array(g.RADII_CM, dtype=float)
    meas = 10.0 ** (-1.0 - 0.01 * np.arange(46)[:, None] - 0.1 * np.arange(10)[None, :])
    b = np.broadcast_to(meas, (20, 46, 10)) * (
        1.0 + 1e-9 * np.random.default_rng(1).standard_normal((20, 46, 10))
    )
    pending = g.evaluate(z, r, meas, b)
    assert pending["status"] == "pending" and pending["gate_pass"] is None
    gating = g.evaluate(z, r, meas, b, fixed_status="gating")
    rep = g.evaluate(z, r, meas, b, fixed_status="report-only")
    assert gating["gate_pass"] is True and gating["status"] == "gating" and rep["gate_pass"] is None
    assert g.combine_backends(gating, gating)["pass"] is True
    assert (
        g.combine_backends(gating, pending)["pass"] is False
    )  # status disagreement, never a silent pass
    assert "status" in g.combine_backends(gating, pending)["disagreement_components"]
    with pytest.raises(ValueError):
        g.evaluate(z, r, meas, b, fixed_status="pending")
    bad = g.evaluate(
        z,
        r,
        meas,
        b * np.where(np.arange(10) == 0, 1.05, 1.0)[None, None, :],
        fixed_status="report-only",
    )
    assert bad["r0"]["pass"] is False
    c = g.combine_backends(rep, bad)
    assert c["pass"] is False and "r0" in c["disagreement_components"]
    assert g.combine_backends(rep, rep)["pass"] is True
