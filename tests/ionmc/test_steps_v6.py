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
