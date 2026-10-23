"""Tests of the ENDF-6 identity / coverage verifier (validation/scripts/reference_data).

All inputs are synthetic ENDF-6 text built here; no data file is needed.
"""

from __future__ import annotations

import hashlib
import json
import math
import zipfile
from pathlib import Path
from typing import Any

import pytest

from validation.scripts.reference_data import endf_verify as ev

MAT = 625
ZA = 6012.0
AWR = 11.89365
LEGENDRE_MARK = 0.7777  # a Legendre coefficient that must never reach the record


# ---------------------------------------------------------------- synthetic ENDF-6 builders
def f(x: float) -> str:
    """ENDF-6 real in an 11-character field."""
    if x == 0.0:
        return "0.000000+0".rjust(11)
    mant, exp = f"{x:.6e}".split("e")
    return f"{mant}{int(exp):+d}".rjust(11)


def i(n: int) -> str:
    return str(n).rjust(11)


def cont(c1: float, c2: float, l1: int, l2: int, n1: int, n2: int) -> str:
    return f(c1) + f(c2) + i(l1) + i(l2) + i(n1) + i(n2)


def values(vals: list[float]) -> list[str]:
    return ["".join(f(v) for v in vals[k : k + 6]).ljust(66) for k in range(0, len(vals), 6)]


def interp(pairs: list[tuple[int, int]]) -> list[str]:
    flat = [x for p in pairs for x in p]
    return ["".join(i(x) for x in flat[k : k + 6]).ljust(66) for k in range(0, len(flat), 6)]


def tab1(
    c1: float, c2: float, pairs: list[tuple[int, int]], xs: list[float], ys: list[float]
) -> list[str]:
    inter = [v for xy in zip(xs, ys, strict=True) for v in xy]
    return [cont(c1, c2, 0, 0, len(pairs), len(xs)), *interp(pairs), *values(inter)]


def tab2(ne: int) -> list[str]:
    return [cont(0, 0, 0, 0, 1, ne), *interp([(ne, 2)])]


def section(mf: int, mt: int, rows: list[str], mat: int = MAT) -> list[str]:
    out = [r.ljust(66)[:66] + f"{mat:4d}{mf:2d}{mt:3d}{n:5d}" for n, r in enumerate(rows, 1)]
    out.append(" " * 66 + f"{mat:4d}{mf:2d}{0:3d}{0:5d}")  # SEND
    return out


def mu_a() -> tuple[list[float], list[float]]:
    return [-1.0, -0.5, 0.0, 0.5, 1.0], [1.0, 0.6, 0.8, 1.2, 1.6]


def mu_b() -> tuple[list[float], list[float]]:
    return [-1.0, 0.0, 1.0], [1.0, 2.0, 4.0]  # f = 2 * 2**mu, INT 4 (log-y, lin-mu)


def endf_text(law_a: int = 2, law_b: int = 4) -> str:
    comments = [
        "The extension to 150 MeV was not retained in this file.",
        "Second comment line.",
    ]
    directory = [cont(0, 0, mf, mt, 5, 0) for mf, mt in ((2, 151), (3, 2), (4, 2), (34, 2))]
    mf1 = [
        cont(ZA, AWR, 0, 0, 0, 0),
        cont(0, 0, 0, 0, 0, 6),
        cont(1.0, 2.0e7, 0, 0, 10, 8),
        cont(0, 0, 0, 0, len(comments), len(directory)),
        *[c.ljust(66) for c in comments],
        *directory,
    ]
    mf2 = [cont(ZA, AWR, 0, 0, 1, 0), cont(ZA, 1.0, 0, 0, 1, 0), cont(1e-5, 1e5, 0, 0, 0, 0)]
    xs = [1e-5, 50e6, 100e6, 150e6]
    mf3_1 = [cont(ZA, AWR, 0, 0, 0, 0), *tab1(0, 0, [(4, 2)], xs, [4.0, 3.0, 2.0, 1.0])]
    mf3_2 = [cont(ZA, AWR, 0, 0, 0, 0), *tab1(0, 0, [(4, 2)], xs, [3.0, 2.0, 1.0, 0.4])]
    mua, fa = mu_a()
    mub, fb = mu_b()
    mf4 = [
        cont(ZA, AWR, 0, 3, 0, 0),
        cont(0, AWR, 0, 2, 0, 0),
        *tab2(2),
        cont(0, 10e6, 0, 0, 3, 0),
        *values([0.1, LEGENDRE_MARK, 0.3]),
        cont(0, 20e6, 0, 0, 2, 0),
        *values([0.2, 0.4]),
        *tab2(2),
        *tab1(0, 30e6, [(len(mua), law_a)], mua, fa),
        *tab1(0, 40e6, [(len(mub), law_b)], mub, fb),
    ]
    mf34 = [
        cont(ZA, AWR, 0, 1, 0, 1),
        cont(0, 0, 0, 2, 1, 1),
        cont(0, 0, 1, 1, 0, 1),
        cont(0, 0, 1, 5, 3, 2),
        *values([1e-5, 7.0e6, 0.01]),
    ]
    lines = [" Synthetic tape label".ljust(66) + f"{0:4d}{0:2d}{0:3d}{0:5d}"]
    for mf, mt, rows in (
        (1, 451, mf1),
        (2, 151, mf2),
        (3, 1, mf3_1),
        (3, 2, mf3_2),
        (4, 2, mf4),
        (34, 2, mf34),
    ):
        lines += section(mf, mt, rows)
    lines.append(" " * 66 + f"{MAT:4d}{0:2d}{0:3d}{0:5d}")  # FEND-like
    lines.append(" " * 66 + f"{-1:4d}{0:2d}{0:3d}{0:5d}")  # TEND
    return "\n".join(lines) + "\n"


CHANGELOG = (
    "line one\n"
    "    * n-006_C_012.endf: natC evaluation restored from 20-150 MeV.\n"
    "unrelated\n"
    "* bug fixes in n-008_O_016.endf and others\n"
)


def build_zip(path: Path) -> None:
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("SUB/n-006_C_012.endf", endf_text())
        zf.writestr("SUB/CHANGELOG.txt", CHANGELOG)


def meta_for(zip_path: Path) -> dict[str, Any]:
    raw = zip_path.read_bytes()
    return {
        "record": "synthetic",
        "operator_authorisation": ["synthetic"],
        "acquisition": {
            "url": "https://example.invalid/x.zip",
            "bytes": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "retrieved_at_utc": "2026-01-01T00:00:00Z",
            "acquisition_id": "abc",
            "role": "exploratory",
            "acquisition_code_sha": "0" * 40,
            "citation": "cite",
            "license_basis": "public",
        },
        "findings": ["finding one"],
        "evaluation_comment_quotes": {"n-006_C_012.endf": "quote"},
        "git_policy": "decision 0041",
        "status": "VERIFIED",
    }


@pytest.fixture(autouse=True)
def _stub_external(monkeypatch: pytest.MonkeyPatch) -> None:
    """Fixed provenance (no git call) and a known Barashenkov stand-in."""
    monkeypatch.setattr(
        ev,
        "provenance",
        lambda argv: {
            "verifier_sha256": "v" * 64,
            "git_head": "h" * 40,
            "git_clean": True,
            "argv": ["python", "-m", ev.MODULE_NAME, *argv],
            "verifier": "x",
        },
    )
    monkeypatch.setattr(ev, "barashenkov_elastic_mb", lambda e, z: {100.0: 500.0, 150.0: 800.0}[e])


def parsed() -> dict[str, Any]:
    return ev.verify_member(endf_text().encode("ascii"), "n-006_C_012.endf", "SUB/n-006_C_012.endf")


# ---------------------------------------------------------------- float and record layer
@pytest.mark.parametrize(
    ("text", "value"),
    [
        ("1.0+5", 1.0e5),
        ("-1.5-3", -1.5e-3),
        ("  2.5000000+0", 2.5),
        ("", 0.0),
        ("   ", 0.0),
        ("1.234567+10", 1.234567e10),
        ("-1.000000-5", -1.0e-5),
        ("1.5e3", 1500.0),
        ("7", 7.0),
    ],
)
def test_endf_float(text: str, value: float) -> None:
    assert ev.endf_float(text) == pytest.approx(value, rel=1e-12)


def test_endf_float_rejects_garbage() -> None:
    for bad in ("abc", "nan", "1.0++5"):
        with pytest.raises(ValueError):
            ev.endf_float(bad)


def test_endf_int_blank_is_zero() -> None:
    assert ev.endf_int("          ") == 0
    assert ev.endf_int("        451") == 451


def test_sectioning_and_header() -> None:
    lines = endf_text().splitlines()
    secs = ev.split_sections(lines)
    assert {k[1:] for k in secs} == {(1, 451), (2, 151), (3, 1), (3, 2), (4, 2), (34, 2)}
    rec = parsed()
    assert rec["MAT"] == MAT
    h = rec["mf1_mt451"]
    assert (h["ZA"], h["NFOR"], h["NSUB"], h["NVER"], h["NWD"], h["NXC"]) == (ZA, 6, 10, 8, 2, 4)
    assert h["EMAX_eV"] == 2.0e7
    assert h["AWI"] == 1.0
    assert h["comment_lines_1_to_8"][1] == "Second comment line."
    assert [q["comment_line"] for q in h["comment_lines_mentioning_retained"]] == [1]
    assert h["directory_mf_mt"] == [[2, 151], [3, 2], [4, 2], [34, 2]]
    assert h["directory_has_mf4_mt2"] and not h["directory_has_mf6_mt2"]
    assert rec["mf2_mt151"] == {
        "EL_eV": 1e-5,
        "EH_eV": 1e5,
        "LRU": 0,
        "LRF": 0,
        "NRO": 0,
        "NAPS": 0,
    }
    assert rec["mf_present"] == [1, 2, 3, 4, 34]
    assert rec["mt_in_mf3"] == [1, 2]
    assert rec["covariance_mf"] == [34]
    assert rec["mf6_mt2_present"] is False
    assert rec["sha256"] == hashlib.sha256(endf_text().encode()).hexdigest()


def test_unparsable_identifiers_fail_closed() -> None:
    lines = endf_text().splitlines()
    lines[5] = lines[5][:66] + "XXXX" + lines[5][70:]
    with pytest.raises(ValueError, match="unparsable"):
        ev.split_sections(lines)


def test_mf3_tab1() -> None:
    rec = parsed()
    m = rec["mf3_mt2"]
    assert (m["NR"], m["NP"], m["E_min_eV"], m["E_max_eV"]) == (1, 4, 1e-5, 150e6)
    assert m["interp"] == [[4, 2]]
    assert m["n_points_ge_20MeV"] == 3
    assert rec["mf3_mt1"]["NP"] == 4


def test_mf4_separates_legendre_and_tabulated_energies() -> None:
    f4 = parsed()["mf4_mt2"]
    assert (f4["LTT"], f4["LI"], f4["LCT"]) == (3, 0, 2)
    leg, tab = f4["legendre"], f4["tabulated"]
    assert (leg["NE"], leg["E_min_eV"], leg["E_max_eV"], leg["NL_max"]) == (2, 10e6, 20e6, 3)
    assert leg["incident_energies_ge_20MeV_eV"] == [20e6]
    assert (tab["NE"], tab["E_min_eV"], tab["E_max_eV"], tab["NP_max"]) == (2, 30e6, 40e6, 5)
    assert tab["incident_energies_eV"] == [30e6, 40e6]
    assert tab["interp_laws_present"] == [2, 4]


def test_legendre_coefficients_are_not_recorded() -> None:
    assert str(LEGENDRE_MARK) not in ev.dumps({"r": parsed()})


def test_mf34_header_parse() -> None:
    f34 = parsed()["mf34_mt2"]
    assert f34["NMT1"] == 1
    (block,) = f34["blocks"]
    assert (block["MT1"], block["NL"], block["NL1"]) == (2, 1, 1)
    (sub,) = block["sub"]
    assert (sub["L"], sub["L1"], sub["LB"]) == (1, 1, [5])
    assert (sub["E_min_eV"], sub["E_max_eV"]) == (1e-5, 7.0e6)


def test_trailing_records_fail_closed() -> None:
    text = endf_text().replace(
        cont(0, 20e6, 0, 0, 2, 0), cont(0, 20e6, 0, 0, 1, 0)
    )  # declared NL shorter than the data: the section no longer parses cleanly
    with pytest.raises(ValueError):
        ev.verify_member(text.encode(), "x", "x")


# ---------------------------------------------------------------- interpolation laws, moments
def test_segment_laws_mixed_regions() -> None:
    t = ev.Tab1((0, 0, 0, 0), [(2, 2), (4, 4)], [0, 1, 2, 3], [1, 1, 1, 1])
    assert t.segment_laws() == [2, 4, 4]


def test_interpolate_all_laws() -> None:
    xs, ys = [1.0, 4.0], [2.0, 8.0]
    for law, x, expect in (
        (1, 2.0, 2.0),
        (2, 2.5, 5.0),
        (3, 2.0, 2.0 + 6.0 * 0.5),
        (4, 2.5, 2.0 * 4.0**0.5),
        (5, 2.0, 2.0 * 4.0**0.5),
    ):
        t = ev.Tab1((0, 0, 0, 0), [(2, law)], xs, ys)
        assert ev.interpolate(t, x) == pytest.approx(expect)
    with pytest.raises(ValueError, match="outside"):
        ev.interpolate(ev.Tab1((0, 0, 0, 0), [(2, 2)], xs, ys), 5.0)


def test_interpolate_unsupported_law() -> None:
    with pytest.raises(ValueError, match="unsupported"):
        ev.interpolate(ev.Tab1((0, 0, 0, 0), [(2, 6)], [1.0, 4.0], [2.0, 8.0]), 2.0)


def stats_by_energy() -> dict[float, dict[str, Any]]:
    return {s["E_eV"]: s for s in parsed()["mf4_mt2"]["tabulated_statistics"]}


def test_stats_int2_hand_computed() -> None:
    s = stats_by_energy()[30e6]
    norm = 0.4 + 0.35 + 0.5 + 0.7  # trapezoids of the four lin-lin segments
    first = -0.308333333 - 0.083333333 + 0.133333333 + 0.533333333  # int mu f, per segment
    assert s["norm_integral"] == pytest.approx(norm, rel=1e-6)
    assert s["mean_one_minus_mu"] == pytest.approx(1.0 - first / norm, rel=1e-6)
    assert s["interp_laws"] == [2]
    assert s["minima_deg_forward_of_130"] == [120.0]  # f(-0.5) = 0.6 < both neighbours


@pytest.mark.parametrize(("theta", "key"), [(10.0, "p_theta_gt_10deg"), (20.0, "p_theta_gt_20deg")])
def test_exceedance_interpolated_to_threshold_int2(theta: float, key: str) -> None:
    """The threshold lies between the last two tabulated points (mu = 0.5 and 1)."""
    mu0 = math.cos(math.radians(theta))
    f0 = 1.2 + 0.4 * (mu0 - 0.5) / 0.5
    partial = 0.5 * (1.2 + f0) * (mu0 - 0.5)
    expect = (0.4 + 0.35 + 0.5 + partial) / 1.95
    assert stats_by_energy()[30e6][key] == pytest.approx(expect, rel=1e-6)
    truncated = (0.4 + 0.35 + 0.5) / 1.95  # the biased value (stopping at the last node below)
    assert abs(expect - truncated) > 0.05


def test_stats_int4_exponential_exact() -> None:
    s = stats_by_energy()[40e6]
    k = math.log(2.0)
    norm = 2.0 * (2.0 - 0.5) / k
    first = 2.0 * (2.5 / k - 1.5 / k**2)
    assert s["interp_laws"] == [4]
    assert s["norm_integral"] == pytest.approx(norm, rel=1e-7)
    assert s["mean_one_minus_mu"] == pytest.approx(1.0 - first / norm, rel=1e-6)
    for theta, key in ((10.0, "p_theta_gt_10deg"), (20.0, "p_theta_gt_20deg")):
        mu0 = math.cos(math.radians(theta))
        assert s[key] == pytest.approx(2.0 * (2.0**mu0 - 0.5) / k / norm, rel=1e-7)
    assert s["minima_deg_forward_of_130"] == []


def test_local_minima_detector() -> None:
    mu = [math.cos(math.radians(a)) for a in (180, 170, 160, 150, 100, 60, 40, 30, 20, 10, 0)]
    mu = sorted(mu)  # ascending mu
    ang = [round(math.degrees(math.acos(m))) for m in mu]
    fmu = [1.0] * len(mu)
    # minimum at 160 deg (backward, excluded), strict minimum at 40 deg, plateau at 100/150 deg
    dip = {160: 0.5, 40: 0.5}
    fmu = [dip.get(a, 1.0 + 0.01 * k) for k, a in enumerate(ang)]
    tab = ev.Tab1((0, 30e6, 0, 0), [(len(mu), 2)], mu, fmu)
    assert ev.energy_statistics(tab)["minima_deg_forward_of_130"] == [40.0]
    flat = ev.Tab1((0, 30e6, 0, 0), [(3, 2)], [-1.0, 0.0, 1.0], [1.0, 1.0, 1.0])
    assert ev.energy_statistics(flat)["minima_deg_forward_of_130"] == []


def test_unsupported_tabulated_law_fails_closed() -> None:
    for law in (1, 3, 5):
        text = endf_text(law_a=law)
        with pytest.raises(ValueError, match="unsupported"):
            ev.verify_member(text.encode(), "x", "x")


def test_log_segment_with_zero_ordinate_fails_closed() -> None:
    t = ev.Tab1((0, 30e6, 0, 0), [(3, 4)], [-1.0, 0.0, 1.0], [0.0, 1.0, 2.0])
    with pytest.raises(ValueError, match="log undefined"):
        ev.energy_statistics(t)


# ---------------------------------------------------------------- sigma_el ratio, changelog
def test_sigma_ratio_uses_barashenkov_and_barn_to_mb() -> None:
    r = parsed()["sigma_el_over_barashenkov"]
    # MF3/MT2: 1.0 b at 100 MeV and 0.4 b at 150 MeV; stand-in formula: 500 and 800 mb
    assert r == {"100MeV": 2.0, "150MeV": 0.5}
    assert parsed()["Z"] == 6


def test_sigma_ratio_interpolates_with_own_law() -> None:
    tab = ev.Tab1((0, 0, 0, 0), [(3, 2)], [90e6, 100e6, 150e6], [2.0, 1.0, 0.4])
    assert ev.sigma_el_ratios(tab, 6, (100.0,)) == {"100MeV": 2.0}


def test_changelog_record() -> None:
    c = ev.changelog_record(CHANGELOG.encode(), ["n-006_C_012.endf", "n-008_O_016.endf"])
    assert c["sha256"] == hashlib.sha256(CHANGELOG.encode()).hexdigest()
    assert c["bytes"] == len(CHANGELOG) and c["lines"] == 4
    assert c["lines_mentioning_members"] == [
        {
            "line": 2,
            "member": "n-006_C_012.endf",
            "text": "    * n-006_C_012.endf: natC evaluation restored from 20-150 MeV.",
        },
        {
            "line": 4,
            "member": "n-008_O_016.endf",
            "text": "* bug fixes in n-008_O_016.endf and others",
        },
    ]


# ---------------------------------------------------------------- record, renderer, CLI
def build(tmp_path: Path) -> tuple[dict[str, Any], Path]:
    z = tmp_path / "src.zip"
    build_zip(z)
    rec = ev.build_record(
        z, ["SUB/n-006_C_012.endf"], "SUB/CHANGELOG.txt", meta_for(z), ["--zip", str(z)]
    )
    return rec, z


def all_keys(obj: Any) -> list[str]:
    if isinstance(obj, dict):
        return [str(k) for k in obj] + [k for v in obj.values() for k in all_keys(v)]
    if isinstance(obj, list):
        return [k for v in obj for k in all_keys(v)]
    return []


def test_render_is_pure_function_of_json(tmp_path: Path) -> None:
    rec, _ = build(tmp_path)
    again = json.loads(ev.dumps(rec))
    md1, md2 = ev.render_markdown(rec), ev.render_markdown(again)
    assert md1 == md2 == ev.render_markdown(again)
    assert "Legendre-part incident energies >= 20 MeV (MeV): 20" in md1
    assert "tabulated-part incident energies (MeV): 30, 40" in md1
    assert "n-006_C_012.endf` | 625" in md1
    assert "was not retained in this file" in md1
    assert "| 30 | 5 | [2] |" in md1 and "| 40 | 3 | [4] |" in md1


def test_record_has_no_absolute_cross_section_key(tmp_path: Path) -> None:
    rec, _ = build(tmp_path)
    for key in all_keys(rec):
        low = key.lower()
        assert not (low.endswith("_mb") or "barn" in low or low in ("sigma_el", "xs")), key
    assert set(rec["members"]["n-006_C_012.endf"]["sigma_el_over_barashenkov"]) == {
        "100MeV",
        "150MeV",
    }


def test_acquisition_mismatch_fails_closed(tmp_path: Path) -> None:
    z = tmp_path / "src.zip"
    build_zip(z)
    meta = meta_for(z)
    meta["acquisition"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="does not match"):
        ev.build_record(z, ["SUB/n-006_C_012.endf"], "SUB/CHANGELOG.txt", meta, [])


def test_cli_writes_deterministic_files(tmp_path: Path) -> None:
    z = tmp_path / "src.zip"
    build_zip(z)
    meta = tmp_path / "meta.json"
    meta.write_text(json.dumps(meta_for(z)), encoding="utf-8")
    args = [
        "--zip",
        str(z),
        "--member",
        "SUB/n-006_C_012.endf",
        "--changelog",
        "SUB/CHANGELOG.txt",
        "--acquisition-meta",
        str(meta),
    ]
    out1, out2 = tmp_path / "o1", tmp_path / "o2"
    assert ev.main([*args, "--out", str(out1)]) == 0
    assert ev.main([*args, "--out", str(out2)]) == 0
    j1 = (out1 / "acquisition.json").read_text(encoding="utf-8")
    assert j1 == (out2 / "acquisition.json").read_text(encoding="utf-8").replace(
        str(out2), str(out1)
    )
    rec = json.loads(j1)
    assert rec["acquisition"]["zip_entries"] == 2
    assert rec["provenance"]["zip_sha256"] == rec["acquisition"]["sha256"]
    assert (out1 / "ACQUISITION.md").read_text(encoding="utf-8") == ev.render_markdown(rec)
    assert j1.endswith("\n") and list(rec) == sorted(rec)


def test_script_hash_is_of_the_file_itself() -> None:
    assert ev.file_sha256(ev.SCRIPT_PATH) == hashlib.sha256(ev.SCRIPT_PATH.read_bytes()).hexdigest()
