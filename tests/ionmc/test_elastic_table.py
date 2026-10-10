"""V3-005C elastic table: shape, R rule, kinematics (P6), builder fixtures, loader fail-closed cases
(C1-ext), packed device, O-16 copy assertion, X-ENDF table presence and the P7 harness.
Synthetic cases
run without data; table cases build the real table once (single process, ~15 s) from the cache and
fail under ``IONMC_REQUIRE_DATA=1`` when the data are missing."""

from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ionmc.data import cache, endf6
from ionmc.data.registry import DATASETS
from ionmc.errors import UnsupportedCombinationError
from ionmc.materials import ELEMENTS, WATER, Material
from ionmc.nuclear import bgg, law5
from ionmc.nuclear import elastic_build as EB
from ionmc.nuclear.build import TARGETS, BuildError, table_id, write_npz_deterministic
from ionmc.nuclear.elastic_kin import two_body
from ionmc.nuclear.elastic_tables import (
    ElasticTable,
    ElasticTableMissingError,
    ElasticTablePinError,
    ElasticTableStaleError,
    ElasticTableUnqualifiedError,
    chi2_equiprobable,
    gamma_q,
)

MP = 938.27208816


# ---------------------------------------------------------------------------------------------
# shape, R rule, J1
# ---------------------------------------------------------------------------------------------
def test_bessel_j1_against_tabulated_values_and_series() -> None:
    x = np.array([1.0, 10.0, 3.831705970207512])
    j = EB.bessel_j1(x)
    assert abs(j[0] - 0.44005058574493355) < 1e-13
    assert abs(j[1] - 0.04347274616886144) < 1e-13
    assert abs(j[2]) < 1e-13  # first zero of J1
    assert abs(EB.bessel_j1(np.array([25.0]))[0] + 0.1253502) < 1e-6

    def series(v: float) -> float:
        return sum(
            (-1) ** m / (math.factorial(m) * math.factorial(m + 1)) * (v / 2) ** (2 * m + 1)
            for m in range(40)
        )

    assert abs(EB.bessel_j1(np.array([2.5]))[0] - series(2.5)) < 1e-14
    assert EB.disk_amplitude_sq(np.array([0.0]))[0] == 1.0


def test_radius_rule_inverts_sigma_nonel() -> None:
    for sig_b, t in ((0.3, 50.0), (0.46, 150.0), (0.30, 250.0)):
        r, lam = EB.disk_radius_fm(sig_b, t, 14899.17)
        assert math.pi * (r + lam) ** 2 == pytest.approx(100.0 * sig_b, rel=1e-14)
        pcm, _ = EB.kinematics_cm(t, 14899.17)
        assert lam == pytest.approx(EB.HBARC / pcm, rel=1e-14)


def test_disk_quantiles_round_trip_and_normalisation() -> None:
    r, pcm = 2.68, 560.0
    edges, _ = EB.disk_quantiles(r, pcm, 16385, 256)
    assert edges[0] == -1.0 and edges[-1] == 1.0 and np.all(np.diff(edges) >= 0.0)
    # independent CDF on a 4x finer grid: F(mu_i) = 1 - (1 - i/n) ... equiprobable bins
    s = np.linspace(0.0, 1.0, 65537)
    g = 4.0 * s * EB.disk_amplitude_sq(2.0 * pcm / EB.HBARC * s * r)
    cum = np.concatenate(([0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(s))))
    cum /= cum[-1]
    s_edges = np.sqrt((1.0 - edges) / 2.0)
    f = np.interp(s_edges, s, cum)  # cumulative from the forward direction
    assert np.max(np.abs(f - (1.0 - np.arange(257) / 256.0))) < 1e-5
    # dsigma/dOmega integrates to sigma_el (P6 (6))
    th = np.linspace(0.0, 180.0, 200001)
    d = EB.disk_density_mb_sr(0.159, r, pcm, th)
    integral = 2.0 * math.pi * np.trapezoid(d * np.sin(np.radians(th)), np.radians(th)) * 1e-3
    assert integral == pytest.approx(0.159, rel=1e-5)
    # first diffraction minimum at q R = 3.8317
    q_min = 3.831705970207512 / r
    th_min = 2.0 * math.degrees(math.asin(q_min * EB.HBARC / (2.0 * pcm)))
    dd = EB.disk_density_mb_sr(0.159, r, pcm, np.array([th_min - 0.5, th_min, th_min + 0.5]))
    assert dd[1] < 1e-6 * dd.max() + 1e-9 and dd[1] < dd[0] and dd[1] < dd[2]


def test_symmetric_quantiles_of_a_known_density() -> None:
    mu = np.linspace(-1.0, 1.0, 20001)
    edges = EB.symmetric_quantiles(mu, 1.0 + mu, 8)  # density (1 + mu)/2: CDF (1 + mu)^2 / 4
    assert np.allclose(edges, 2.0 * np.sqrt(np.arange(9) / 8.0) - 1.0, atol=2e-4)


# ---------------------------------------------------------------------------------------------
# kinematics
# ---------------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("m2", "t"), [(MP, 20.0), (MP, 150.0), (14899.17, 100.0), (37262.0, 250.0)]
)
def test_two_body_closure_and_equal_mass_recoil(m2: float, t: float) -> None:
    mu = np.linspace(-1.0, 1.0, 20001)
    k = two_body(t, MP, m2, mu)
    de, dp = k.closure(t, MP, m2)
    assert de <= 1e-9 and dp <= 1e-9
    # on-shell: E^2 - p^2 = m^2 for both outgoing particles
    assert np.max(np.abs(k.e1**2 - k.px1**2 - k.pz1**2 - MP**2)) <= 1e-6
    assert np.max(np.abs(k.e2**2 - k.px2**2 - k.pz2**2 - m2**2)) <= 1e-6
    if m2 == MP:
        assert np.max(np.abs(k.t2 - t * (1.0 - mu) / 2.0)) <= 1e-9
        assert np.max(np.abs(k.t1 + k.t2 - t)) <= 1e-9


# ---------------------------------------------------------------------------------------------
# grid / refinement / chi-square utility
# ---------------------------------------------------------------------------------------------
def test_grid_elastic_contains_anchor_14_mev_barashenkov_nodes_and_ends_at_250() -> None:
    g = EB.build_grid_elastic(np.array([1.0, 10.0, 150.0]), 50)
    assert g[0] == 1.0 and g[-1] == 250.0 and np.all(np.diff(g) > 0.0)
    for e in (14.0, 15.0, 17.0, 150.0, 160.0, 180.0, 200.0, 250.0, 100.0):
        assert np.any(np.abs(g - e) <= 1e-9), e


def test_refine_grid_meets_the_midpoint_tolerance() -> None:
    def f(e: float) -> list[float]:
        return [max(e - 3.0, 0.0) ** 2]

    g0 = np.linspace(1.0, 14.0, 14)
    g, added = EB.refine_grid(g0, f, 0.0, 14.0)
    assert added > 0
    for a, b in zip(g[:-1], g[1:], strict=True):
        m = 0.5 * (a + b)
        assert (
            abs(0.5 * (f(a)[0] + f(b)[0]) - f(m)[0])
            <= EB.REFINE_TOL * max(f(m)[0], 1e-3 * f(14.0)[0]) * 1.0001
        )


def test_chi2_utility_p_values() -> None:
    assert gamma_q(1.0, 1.0) == pytest.approx(math.exp(-1.0), rel=1e-12)
    assert 0.45 < gamma_q(31.5, 31.5) < 0.5
    rng = np.random.default_rng(1)
    edges = np.linspace(-1.0, 1.0, 257)
    mu = rng.uniform(-1.0, 1.0, 100000)
    assert chi2_equiprobable(mu, edges)[1] > 1e-3
    assert chi2_equiprobable(np.clip(mu * 0.9, -1, 1), edges)[1] < 1e-6  # wrong density detected


# ---------------------------------------------------------------------------------------------
# synthetic sealed table (CI-fast): loader, C1-ext fail-closed cases, composition, device
# ---------------------------------------------------------------------------------------------
def _synthetic(root: Path, *, qual: tuple[int, ...] = (1, 1, 1)) -> str:
    grid = np.array([1.0, 2.0, 5.0, 10.0, 14.0, 50.0, 100.0, 150.0, 250.0])
    n, n_t = grid.size, len(EB.TARGET_NAMES)
    sigma = np.tile(0.1 * (1 + np.arange(n_t))[:, None], (1, n)) * (1 + 0.1 * np.sin(grid))
    edges = np.tile(np.linspace(-1.0, 1.0, 5), (n_t, n, 1))
    e_min = np.array([12.5, 6.0, 3.0, 6.5, 2.5, 4.0, 2.5, 5.7])
    valid = (grid[None, :] >= e_min[:, None]).astype(np.int8)
    sigma = sigma * valid  # sigma = 0 and valid = 0 below the domain
    arrays: dict[str, Any] = {
        "grid_e_mev": grid, "sigma_barn": sigma, "edges_mu": edges,
        "qualification": np.array(qual, dtype=np.int8),
        "target_mass_mev": np.arange(1.0, n_t + 1.0),
        "target_e_min_mev": e_min, "valid": valid,
    }  # fmt: skip
    (root / "derived").mkdir(parents=True)
    tmp = root / "derived" / "x.npz"
    sha = write_npz_deterministic(tmp, arrays)
    info = {
        "schema": EB.SCHEMA, "builder_version": EB.BUILDER_VERSION, "options": {},
        "sources": {
            s: {"sha256": DATASETS[s].sha256, "version": DATASETS[s].version}
            for s in EB.SOURCE_IDS
        },
        "npz_sha256": sha, "target_names": list(EB.TARGET_NAMES),
        "elements": EB.element_rows_elastic(),
    }  # fmt: skip
    tid = table_id(info)
    info["table_id"] = tid
    tmp.replace(root / "derived" / f"elastic-proton-{tid}.npz")
    (root / "derived" / f"elastic-proton-{tid}.json").write_text(json.dumps(info), encoding="utf-8")
    return tid


def _reseal(src: Path, tid: str, dst: Path, edit_a: Any = None, edit_i: Any = None) -> str:
    (dst / "derived").mkdir(parents=True)
    npz = (src / "derived" / f"elastic-proton-{tid}.npz").read_bytes()
    info = json.loads((src / "derived" / f"elastic-proton-{tid}.json").read_text())
    if edit_a is not None:
        import io

        with np.load(io.BytesIO(npz), allow_pickle=False) as z:
            arr = {k: z[k].copy() for k in z.files}
        edit_a(arr)
        tmp = dst / "derived" / "t.npz"
        info["npz_sha256"] = write_npz_deterministic(tmp, arr)
        npz = tmp.read_bytes()
        tmp.unlink()
    if edit_i is not None:
        edit_i(info)
    info.pop("table_id")
    new = table_id(info)
    info["table_id"] = new
    (dst / "derived" / f"elastic-proton-{new}.npz").write_bytes(npz)
    (dst / "derived" / f"elastic-proton-{new}.json").write_text(json.dumps(info), encoding="utf-8")
    return new


def test_synthetic_table_loads_and_composes_water(tmp_path: Path) -> None:
    tid = _synthetic(tmp_path)
    tab = ElasticTable.load(tmp_path, tid)
    rows = tab.material_rows(WATER)
    assert set(rows.target_index) == {0, 3}  # H-1 and O-16
    assert rows.cum_fraction[-1].tolist() == [1.0] * tab.arrays["grid_e_mev"].size
    # Sigma_mass = N_A sum w scale sigma / A [cm2/g] by hand at the 14 MeV node
    k = 4
    h = WATER.mass_fractions["H"] * tab.arrays["sigma_barn"][0, k] / ELEMENTS["H"].A_g_mol
    o = WATER.mass_fractions["O"] * tab.arrays["sigma_barn"][3, k] / ELEMENTS["O"].A_g_mol
    assert rows.sigma_mass_cm2_g[k] == pytest.approx(6.02214076e23 * 1e-24 * (h + o), rel=1e-6)
    assert np.all(rows.sigma_hat_window >= rows.sigma_mass_cm2_g)
    with pytest.raises(UnsupportedCombinationError, match="Cu"):
        tab.material_rows(_copper())
    with pytest.raises(UnsupportedCombinationError):
        tab.sigma_barn(3, 251.0)
    with pytest.raises(UnsupportedCombinationError):
        tab.sigma_barn(3, 0.5)


def _copper() -> Material:
    from ionmc.materials import fractions_from_atom_counts

    return Material("copper", 8.96, fractions_from_atom_counts({"Cu": 1}), 322.0, None, "test")


def test_loader_fail_closed_cases_c1_ext(tmp_path: Path) -> None:
    root = tmp_path / "a"
    tid = _synthetic(root)
    with pytest.raises(ElasticTableMissingError):
        ElasticTable.load(tmp_path / "empty", tid)
    # stale: npz bytes changed
    npz = root / "derived" / f"elastic-proton-{tid}.npz"
    raw = bytearray(npz.read_bytes())
    raw[-3] ^= 0xFF
    bad = tmp_path / "stale"
    shutil.copytree(root, bad)
    (bad / "derived" / npz.name).write_bytes(bytes(raw))
    with pytest.raises(ElasticTableStaleError):
        ElasticTable.load(bad, tid)
    # stale: id not reproduced (sidecar field edited without resealing)
    ed = tmp_path / "edited"
    shutil.copytree(root, ed)
    jp = ed / "derived" / f"elastic-proton-{tid}.json"
    info = json.loads(jp.read_text())
    info["elements"]["O"]["sigma_scale"] = 2.0
    jp.write_text(json.dumps(info))
    with pytest.raises(ElasticTableStaleError):
        ElasticTable.load(ed, tid)

    # mis-pinned source (resealed): pin error
    def mispin(i: dict[str, Any]) -> None:
        i["sources"]["geant4-g4barashenkovdata-hh-11.4.2"]["sha256"] = "0" * 64

    new = _reseal(root, tid, tmp_path / "pin", None, mispin)
    with pytest.raises(ElasticTablePinError):
        ElasticTable.load(tmp_path / "pin", new)

    # wrong schema / builder (resealed): stale
    def schema(i: dict[str, Any]) -> None:
        i["builder_version"] = "old"

    new = _reseal(root, tid, tmp_path / "schema", None, schema)
    with pytest.raises(ElasticTableStaleError):
        ElasticTable.load(tmp_path / "schema", new)
    # unqualified: every qualification flag in turn
    for i, name in enumerate(EB.QUALIFICATION_FIELDS):

        def clear(a: dict[str, Any], i: int = i) -> None:
            a["qualification"][i] = 0

        new = _reseal(root, tid, tmp_path / f"q{i}", clear, None)
        with pytest.raises(ElasticTableUnqualifiedError, match=name):
            ElasticTable.load(tmp_path / f"q{i}", new)


def test_loader_rejects_nonzero_sigma_below_the_domain(tmp_path: Path) -> None:
    root = tmp_path / "d"
    tid = _synthetic(root)

    def bad(a: dict[str, Any]) -> None:
        a["sigma_barn"][3, 0] = 0.5  # O-16 below e_min = 6.5 MeV

    new = _reseal(root, tid, tmp_path / "bad", bad, None)
    with pytest.raises(ElasticTableUnqualifiedError, match="below"):
        ElasticTable.load(tmp_path / "bad", new)


def test_unsupported_error_hierarchy_is_fail_closed() -> None:
    assert issubclass(ElasticTableMissingError, UnsupportedCombinationError)
    assert issubclass(ElasticTableUnqualifiedError, UnsupportedCombinationError)


def test_endf_mt2_is_never_a_construction_input() -> None:
    EB.reject_endf_mt2_construction("H-1")
    for name in ("O-16", "C-12", "N-14", "Ca-40"):
        with pytest.raises(UnsupportedCombinationError):
            EB.reject_endf_mt2_construction(name)


def test_device_pack_hash_and_process_cache(tmp_path: Path) -> None:
    wp = pytest.importorskip("warp")
    from ionmc.transport import elastic_device as D

    tab = ElasticTable.load(tmp_path, _synthetic(tmp_path))
    host = D.pack_elastic(tab, (WATER,))
    assert host.n_targets == 8 and host.n_edges == 5 and host.arrays["edges"].size == 8 * 9 * 5
    dev = D.ElasticDevice(host, wp.float64, "cpu")
    assert dev.sha256 == D.host_sha256(host, "float64") == dev.device_sha256()
    dev32 = D.ElasticDevice(host, wp.float32, "cpu")
    assert dev32.sha256 == D.host_sha256(host, "float32") != dev.sha256
    assert np.array_equal(dev.readback()["mat_target"][:2], [0, 3])
    rb = dev.readback()
    assert np.array_equal(rb["e_min_shape"][1:], tab.arrays["target_e_min_mev"][1:])
    assert rb["e_min_shape"][0] == 0.0 and rb["e_min_pp"].tolist() == [12.5]
    assert tab.elastic_domain()["O-16"] == (6.5, 250.0)
    D.clear_elastic_device_cache()
    a = D.cached_elastic_device(tab, (WATER,), real=wp.float64, device="cpu")
    assert D.cached_elastic_device(tab, (WATER,), real=wp.float64, device="cpu") is a
    assert D.cached_elastic_device(tab, (WATER,), real=wp.float32, device="cpu") is not a
    D.clear_elastic_device_cache()


# ---------------------------------------------------------------------------------------------
# the real table (data-backed)
# ---------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def real() -> tuple[EB.ElasticBuildResult, ElasticTable]:
    cdir = cache.resolve_cache_dir(None)
    try:
        res = EB.build_elastic_proton(cdir)
    except FileNotFoundError:
        if os.environ.get("IONMC_REQUIRE_DATA") == "1":
            pytest.fail("data for the elastic table are not cached")
        pytest.skip("data not cached")
    return res, ElasticTable.load(cdir, res.table_id)


def test_real_table_nodes_midpoints_and_records(real: Any) -> None:
    res, tab = real
    info, a = res.info, tab.arrays
    grid = a["grid_e_mev"]
    assert grid[0] == 1.0 and grid[-1] == 250.0
    assert info["grid"]["max_midpoint_sigma_el_interpolation_error"] <= 1e-3  # P7 sigma check
    assert info["targets"][0]["max_midpoint_interpolation_error"] <= 1e-3
    for it, spec in enumerate(TARGETS, start=1):
        k_min = info["targets"][it]["e_min_shape_index"]
        for k in range(k_min, grid.size):
            assert a["sigma_barn"][it, k] == pytest.approx(
                bgg.bgg_elastic_mb(float(grid[k]), spec.z) * 1e-3,
                rel=(1e-3 if grid[k] < 14.0 else 1e-12),
            )
        assert np.all(a["sigma_barn"][it, :k_min] == 0.0) and np.all(
            a["radius_fm"][it, k_min:] > 0.0
        )
        assert float(a["inversion_residual"][it].max()) < 1e-12
        assert info["sigma_el_10_mev_mb"][spec.name] > 800.0  # recorded for the 0041 comparison
    assert info["targets"][3]["sigma_el_150_mev_mb"] == 159.0
    # R rule uses the transport's own sigma_nonel
    o = 3
    r_o = np.interp(150.0, grid, a["radius_fm"][o])
    lam = np.interp(150.0, grid, a["lambda_bar_fm"][o])
    assert math.pi * (r_o + lam) ** 2 == pytest.approx(
        100.0 * np.interp(150.0, grid, a["sigma_nonel_barn"][o]), rel=1e-12
    )
    assert 2.5 < r_o < 2.9
    assert (
        len(info["sigma_nonel_sha256"]) == 7 and info["o16_finding"]["o16_mt2_is_c12_copy"] is True
    )
    assert info["transcription_verified"] and all(info["transcription_verified"].values())
    # p-p: E_min_pp with first-negative diagnostics; sigma = 0 below
    h = info["targets"][0]
    assert h["negative_density_below_e_min_pp"] and h["first_negative_nodes"]
    assert 12.0 < h["e_min_pp_mev"] < 13.5 and h["e_min_pp"] == h["e_min_pp_mev"]
    assert h["no_negative_density_in_domain"] is True
    k_pp = h["e_min_pp_index"]
    assert np.all(a["sigma_barn"][0, :k_pp] == 0.0) and np.all(a["sigma_barn"][0, k_pp:] > 0.0)
    k150 = int(np.searchsorted(grid, 150.0))
    assert a["sigma_barn"][0, k150] * 1e3 == pytest.approx(26.262, abs=0.01)
    s_ratio = bgg.pp_bgg_mb(200.0) / bgg.pp_bgg_mb(150.0)
    k200 = int(np.searchsorted(grid, 200.0 - 1e-9))
    assert a["sigma_barn"][0, k200] == pytest.approx(
        a["sigma_barn"][0, k150] * bgg.pp_bgg_mb(float(grid[k200])) / bgg.pp_bgg_mb(150.0),
        rel=1e-12,
    )
    assert 0.9 < s_ratio < 1.0
    assert np.array_equal(a["edges_mu"][0, k200], a["edges_mu"][0, k150])  # shape fixed above 150
    # X-ENDF table present (report-only)
    assert a["xendf_la150_ni_mb_sr"].shape == (3, 3, 9) == a["xendf_model_mb_sr"].shape
    assert np.all(np.isfinite(a["xendf_la150_ni_mb_sr"])) and np.all(a["xendf_model_mb_sr"] > 0.0)
    assert info["xendf"]["report_only"] is True


def test_real_table_h1_p6_failure_recorded_and_positive_in_domain(real: Any) -> None:
    """The frozen row P6 fails for H-1 below E_min,pp (recorded, not converted into a pass);
    the density is non-negative at and above E_min,pp (revised domain)."""
    res, tab = real
    grid = tab.arrays["grid_e_mev"]
    h = res.info["targets"][0]
    k_pp = h["e_min_pp_index"]
    # recorded failure: every negative point lies below E_min,pp and has a negative value
    pts = h["first_negative_nodes"]
    assert pts and all(p["min_density_mb_sr"] < 0.0 for p in pts)
    assert max(p["e_mev"] for p in pts) < h["e_min_pp"]
    assert res.info["negative_density_below_e_min_pp"] is True
    assert res.info["o16_mt2_copy_recorded"] is True
    assert res.info["p6"]["no_negative_density_in_domain"] is True
    m = _material("p-001_H_001")
    sec, awi = law5.parse_law5(m.sections[(6, 2)]), law5.projectile_awi(m)
    mu = np.linspace(0.0, law5.MU_CUT_PP, 2001)
    # independent re-evaluation: negative at the recorded points, positive at/above E_min,pp
    for p in pts[-3:]:
        assert law5.ni_density_ltp1(sec, p["e_mev"] * 1e6, mu, awi, 1, 1).min() < 0.0
    for k in range(k_pp, int(np.searchsorted(grid, 150.0))):
        for e in (grid[k], 0.5 * (grid[k] + grid[k + 1])):
            assert law5.ni_density_ltp1(sec, float(e) * 1e6, mu, awi, 1, 1).min() >= 0.0
    # domain block, model revisions, NI above the cut, omitted-event estimates
    assert res.info["elastic_domain"]["H-1"] == [h["e_min_pp"], 250.0]
    dom = tab.elastic_domain()
    assert dom["H-1"] == (h["e_min_pp"], 250.0) and len(dom) == 8
    assert 12.0 < h["e_min_pp"] <= EB.E_MIN_PP_LIMIT_MEV
    for it, spec in enumerate(TARGETS, start=1):
        lo = res.info["targets"][it]["e_min_shape_mev"]
        assert res.info["elastic_domain"][spec.name] == [lo, 250.0] and lo <= 10.0
        assert dom[spec.name] == (lo, 250.0)
        assert res.info["targets"][it]["sigma_bgg_below_domain_max_mb"] > 0.0
    assert "Amendment 15" in res.info["model_revisions"]
    ni = res.info["sigma_ni_above_cut_mb"]
    assert ni["15_mev"] > 0.0 and ni["20_mev"] > 0.0
    po = res.info["pp_omitted_ni_correction"]
    assert (
        res.info["pp_omitted_ni_correction_events_per_history_150mev"] == po["events_per_history"]
    )
    assert res.info["pp_omitted_ni_correction_energy_mev"] == po["energy_mev"]
    assert 0.0 < po["events_per_history"] < 0.1 and 0.0 < po["energy_mev"] < 1.0
    assert po["nodes_mev"][0] == 1.0 and po["nodes_mev"][-1] == po["e_cut_mev"]
    assert len(po["m_barn"]) == len(po["nodes_mev"]) and min(po["m_barn"]) >= 0.0
    assert po["max_mean_recoil_fraction_of_t"] <= 0.5 + 1e-12  # (1 - mu)/2 <= 1/2 on mu >= 0
    assert 0.0 < po["residual_range_below_1mev_g_cm2"] < 0.01
    assert 0.18 < po["residual_range_at_cut_g_cm2"] < 0.21 and "method" in po
    assert "bound" not in " ".join(k for k in res.info if k.startswith("pp_omitted"))
    pa = res.info["pa_omitted_events_per_history_150mev"]
    assert 0.0 < pa["events_per_history"] < 0.01 and 0.05 < pa["residual_range_g_cm2"] < 0.07
    assert 0.0 < pa["residual_range_below_1mev_g_cm2"] < 0.01


def _material(member: str) -> endf6.EndfMaterial:
    z = cache.verify("endf-b8.0-protons", cache.resolve_cache_dir(None))
    return endf6.parse_endf(endf6.read_member(z, f"ENDF-B-VIII.0_protons/{member}.endf"))


def test_o16_copy_assertion_passes_on_la150_and_fails_closed_otherwise(real: Any) -> None:
    c12, o16, n14 = _material("p-006_C_012"), _material("p-008_O_016"), _material("p-007_N_014")
    f = EB.o16_copy_finding(c12, o16)
    assert (
        f["o16_mt2_is_c12_copy"]
        and f["from_mev"] == 24.0
        and f["max_relative_difference_mf3"] < 1e-4
    )
    with pytest.raises(BuildError):
        EB.o16_copy_finding(c12, n14)  # a genuine second evaluation: the assertion fails closed


class _Mat:
    """Fixture material: MF3/MT5 sigma and MF6/MT5 yields built from a base material."""

    def __init__(
        self,
        base: Any,
        sigma: Any = lambda e, y: y,
        yields: Any = lambda e, y: y,
    ) -> None:
        self.base, self.sigma, self.yields = base, sigma, yields

    def cross_section(self, mt: int) -> Any:
        t = self.base.cross_section(mt)
        return endf6.Tab1(t.nbt, t.interp, t.x.copy(), self.sigma(t.x, t.y))

    def products(self, mt: int) -> Any:
        import types

        out = []
        for p in self.base.products(mt).products:
            y = p.yield_
            out.append(types.SimpleNamespace(
                zap=p.zap,
                yield_=endf6.Tab1(y.nbt, y.interp, y.x.copy(), self.yields(y.x, y.y)),
            ))  # fmt: skip
        return types.SimpleNamespace(products=out)


def test_o16_mt5_not_a_copy_real_data(real: Any) -> None:
    c12, o16 = _material("p-006_C_012"), _material("p-008_O_016")
    f = EB.o16_mt5_copy_check(c12, o16)
    assert f["o16_mt5_not_c12_copy"] and f["fraction_rel_diff_gt_1e-3"] > 0.5
    assert f["sigma_not_copy"] and f["mf6_yields_not_copy"]
    assert f["ratio_o16_over_c12_cv"] > 1e-3
    assert f["max_relative_difference"] >= f["min_relative_difference"] >= 0.0
    prods = {k: v for k, v in f["mf6_yields_per_product"].items() if v is not None}
    assert prods and all(v["differs_from_c12"] for v in prods.values())
    assert {"n", "p"} <= set(prods)


@pytest.mark.parametrize(
    "case",
    ["exact", "copy_1e-6", "scaled", "copied_yields_different_sigma"],
)
def test_o16_mt5_copy_variants_fail_closed(real: Any, case: str) -> None:
    c12 = _material("p-006_C_012")
    wiggle = lambda e, y: y * (1.0 + 0.2 * np.sin(np.arange(y.size)))  # noqa: E731
    if case == "exact":
        fake = _Mat(c12)
    elif case == "copy_1e-6":
        fake = _Mat(c12, sigma=lambda e, y: y * (1.0 + 1e-6), yields=lambda e, y: y * (1.0 + 1e-6))
    elif case == "scaled":
        k = (16.0 / 12.0) ** (2.0 / 3.0)
        fake = _Mat(c12, sigma=lambda e, y: y * k, yields=lambda e, y: y * k)
    else:  # sigma genuinely different (non-constant ratio), yields copied
        fake = _Mat(c12, sigma=wiggle)
    with pytest.raises(BuildError, match="copy"):
        EB.o16_mt5_copy_check(c12, fake)  # type: ignore[arg-type]
    # distinct yields but copied (scaled) sigma also fail: the sigma half is gated independently
    if case == "scaled":
        with pytest.raises(BuildError, match="copy"):
            EB.o16_mt5_copy_check(c12, _Mat(c12, sigma=lambda e, y: y * k, yields=wiggle))  # type: ignore[arg-type]


def test_o16_mt5_parser_errors_propagate_as_build_error(real: Any) -> None:
    c12 = _material("p-006_C_012")

    class NoProducts(_Mat):
        def products(self, mt: int) -> Any:
            raise KeyError(mt)

    with pytest.raises(BuildError, match="cannot be evaluated"):
        EB.o16_mt5_copy_check(c12, NoProducts(c12, sigma=lambda e, y: y * 2.0))  # type: ignore[arg-type]


def _syn(xs: list[float], sig: list[float], ys: dict[int, list[float]]) -> Any:
    import types

    def tab(x: list[float], y: list[float]) -> Any:
        return endf6.Tab1(np.array([len(x)]), np.array([2]), np.array(x), np.array(y))

    prods = [types.SimpleNamespace(zap=z, yield_=tab(xs, y)) for z, y in ys.items()]
    return types.SimpleNamespace(
        cross_section=lambda mt: tab(xs, sig),
        products=lambda mt: types.SimpleNamespace(products=prods),
    )


def test_o16_mt5_union_grid_fallback_with_fewer_than_ten_shared_nodes() -> None:
    xc = [1e6 * x for x in (5.0, 10.0, 20.0, 40.0, 80.0, 150.0)]
    xo = [1e6 * x for x in (5.0, 12.0, 25.0, 50.0, 100.0, 150.0)]  # 2 shared nodes
    sc = [0.1, 0.2, 0.3, 0.25, 0.2, 0.15]
    so = [0.1, 0.25, 0.3, 0.35, 0.22, 0.18]
    yc = {1: [1.0, 1.1, 1.2, 1.3, 1.4, 1.5], 1001: [0.5, 0.6, 0.7, 0.6, 0.5, 0.4]}
    yo = {
        1: [1.0, 1.4, 1.1, 1.0, 1.6, 1.2],
        1001: [0.5, 0.3, 0.9, 0.8, 0.3, 0.6],
        2004: [0.0, 0.1, 0.2, 0.2, 0.1, 0.1],
    }
    f = EB.o16_mt5_copy_check(_syn(xc, sc, yc), _syn(xo, so, yo))
    assert f["compared_grid"] == "union_linear_interpolation" and f["n_compared"] == 10
    assert f["o16_mt5_not_c12_copy"]
    # the same fallback still fails closed for a copy
    with pytest.raises(BuildError, match="copy"):
        EB.o16_mt5_copy_check(_syn(xc, sc, yc), _syn(xc, sc, yc))


def test_pp_domain_sign_fixtures_p6d() -> None:
    grid = np.array([1.0, 5.0, 10.0, 12.0, 13.0, 14.0, 20.0, 150.0, 250.0])
    i150 = 7
    # a density negative only below the domain: recorded, the build proceeds
    node = np.array([-3.0, -2.0, -1.0, -0.1, 0.5, 0.6, 0.7, 0.8, np.inf])
    mid = np.array([-2.5, -1.5, -0.5, 0.2, 0.5, 0.6, 0.7, np.inf, np.inf])
    k_pp, pts = EB.pp_scan(grid, node, mid, i150)
    assert k_pp == 4 and len(pts) == 7 and all(p["min_density_mb_sr"] < 0 for p in pts)
    ok = EB.verify_pp_domain(grid, k_pp, i150, lambda e: 0.1)
    assert ok["n_points_checked"] == 4 + 3
    # the separate pass fails the build on a negative in-domain midpoint
    with pytest.raises(BuildError, match="P6-D"):
        EB.verify_pp_domain(grid, k_pp, i150, lambda e: -1e-9 if abs(e - 13.5) < 1e-9 else 0.1)
    # E_min,pp above the 15 MeV bound fails closed
    node2 = np.array([-3.0] * 6 + [0.7, 0.8, np.inf])
    with pytest.raises(BuildError, match="exceeds"):
        EB.pp_scan(grid, node2, np.full(9, np.inf), i150)


def test_domain_limit_fails_closed_above_10_mev() -> None:
    EB.check_domain_limit("O-16", 6.58)
    EB.check_domain_limit("X", 10.0)
    with pytest.raises(BuildError, match="exceeds"):
        EB.check_domain_limit("X", 10.5)


def test_pdg_is_not_a_construction_source() -> None:
    assert "pdg-rpp2022-pp-elastic" not in EB.SOURCE_IDS
    src = (
        Path(EB.__file__).read_text() + (Path(EB.__file__).parent / "elastic_tables.py").read_text()
    )
    assert "pdg-rpp2022" not in src
    pdg = DATASETS["pdg-rpp2022-pp-elastic"]
    assert pdg.role == "evaluation" and "NOT a construction input" in pdg.lineage


@pytest.mark.parametrize(
    ("name", "target", "energies"),
    [("H-1", 0, (20.0, 50.0, 100.0, 150.0, 200.0)),
     ("C-12", 1, (20.0, 50.0, 100.0, 150.0, 200.0, 250.0)),
     ("O-16", 3, (20.0, 50.0, 100.0, 150.0, 200.0, 250.0))],
)  # fmt: skip
def test_p7_harness_sampler_vs_table(
    real: Any, name: str, target: int, energies: tuple[float, ...]
) -> None:
    _, tab = real
    rng = np.random.default_rng(20450715 + target)
    for e in energies:
        mu = tab.sample_mu(target, e, rng.random(100_000))
        _, p, dof = chi2_equiprobable(mu, tab.edges_at(target, e), 64)
        assert dof == 63 and p > 0.001, (name, e, p)


def test_pp_recoil_energy_mean_is_a_quarter_of_t(real: Any) -> None:
    """T/4: the mean energy of the slower proton under the transported density at 150 MeV (a
    symmetric mu-density: slower proton has |mu| mapped to the half sphere)."""
    _, tab = real
    rng = np.random.default_rng(7)
    mu = tab.sample_mu(0, 150.0, rng.random(200_000))
    k = two_body(150.0, MP, MP, np.abs(mu))
    t_slow = np.minimum(k.t1, k.t2)
    assert t_slow.mean() == pytest.approx(150.0 / 4.0, rel=0.05)
    assert np.all(np.abs(k.t1 + k.t2 - 150.0) < 1e-9)


def test_real_table_missing_element_and_energy_fail_closed(real: Any) -> None:
    _, tab = real
    with pytest.raises(UnsupportedCombinationError):
        tab.material_rows(_copper())
    with pytest.raises(UnsupportedCombinationError):
        tab.edges_at(3, 250.5)
    rows = tab.material_rows(WATER)
    assert rows.sigma_at(100.0) > 0.0 and rows.sigma_at(150.0) < rows.sigma_at(100.0)


def test_real_build_is_deterministic(real: Any) -> None:
    res, _ = real
    again = EB.build_elastic_proton(cache.resolve_cache_dir(None))
    assert again.table_id == res.table_id
    assert hashlib.sha256(again.npz_path.read_bytes()).hexdigest() == res.info["npz_sha256"]


def test_pin_file_matches_the_built_table(real: Any) -> None:
    """The committed pin (ids and hashes only) equals the deterministic real build."""
    path = Path(__file__).resolve().parents[2] / "src" / "ionmc" / "data" / "elastic_table_pin.json"
    pin = json.loads(path.read_text(encoding="utf-8"))
    res, _ = real
    assert pin["builder_version"] == EB.BUILDER_VERSION and pin["table_schema"] == EB.SCHEMA
    assert pin["table_id"] == res.table_id and pin["npz_sha256"] == res.info["npz_sha256"]
    assert pin["grid_points"] == res.info["grid"]["n_points"]


def test_omitted_events_never_evaluate_sigma_below_the_table_floor() -> None:
    seen: list[float] = []

    def sigma(e: float) -> float:
        seen.append(e)
        return 100.0

    n = EB.omitted_events_per_history(6.0, 1.0e22, sigma)
    assert min(seen) >= 1.0 and n > 0.0
    # integral only: n * sigma * (R(6 MeV) - R(1 MeV)) with the project range table
    r = EB.residual_range_g_cm2
    assert n == pytest.approx(1.0e22 * 100.0e-27 * (r(6.0) - r(1.0)), rel=2e-3)
