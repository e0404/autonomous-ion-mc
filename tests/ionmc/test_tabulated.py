"""Tests of the tabulated stopping-power layer (reference paths).

Acceptance criteria are those fixed in decision 0008 before the comparisons
were run; the NIST PSTAR CSDA ranges quoted here were retrieved from
https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html (liquid water,
material 276) on 2026-09-10 and are a de-minimis quotation for validation.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc import materials, particles
from ionmc.backend import reference
from ionmc.data import MCSQUARE_G4_WATER, MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import prepare_stopping_table
from ionmc.reference_data.pstar_water import PSTAR_WATER_STOPPING_POWER
from ionmc.stopping_power import AnalyticStoppingPower
from ionmc.tabulated_stopping_power import TabulatedStoppingPower

#: NIST PSTAR CSDA ranges of protons in liquid water, g/cm^2 (4 significant figures).
PSTAR_CSDA_RANGE = {100.0: 7.718, 150.0: 15.77, 200.0: 25.96, 250.0: 37.94}
#: PSTAR CSDA range at the table floor (0.5 MeV) and at the analytic floor (1 MeV).
PSTAR_CSDA_AT_0_5_MEV = 8.869e-4
PSTAR_CSDA_AT_1_MEV = 2.458e-3

TOL_TABLE_VS_ANALYTIC_ABOVE_10_MEV = 0.010  # decision 0008
TOL_TABLE_RANGE_VS_PSTAR = 0.001  # decision 0008


@pytest.fixture(scope="module")
def pstar_table(pstar_cache_root):
    return TabulatedStoppingPower.from_dataset(
        MCSQUARE_PSTAR_WATER, materials.WATER, particles.PROTON, pstar_cache_root
    )


@pytest.fixture(scope="module")
def pstar_table_python(pstar_cache_root):
    return TabulatedStoppingPower.from_dataset(
        MCSQUARE_PSTAR_WATER,
        materials.WATER,
        particles.PROTON,
        pstar_cache_root,
        path="python",
    )


def test_table_reproduces_committed_pstar_subset_exactly(pstar_table) -> None:
    energies = np.array(sorted(PSTAR_WATER_STOPPING_POWER))
    expected = np.array([PSTAR_WATER_STOPPING_POWER[e] for e in energies])
    np.testing.assert_allclose(
        pstar_table.mass_stopping_power(energies), expected, rtol=1e-12
    )
    assert pstar_table.energy_range_mev == (0.5, 400.0)
    assert pstar_table.table.size == 800


def test_numpy_path_is_bitwise_identical_to_python_path(
    pstar_table, pstar_table_python
):
    grid = np.linspace(0.5, 400.0, 1601)
    np.testing.assert_array_equal(
        pstar_table.mass_stopping_power(grid),
        pstar_table_python.mass_stopping_power(grid),
    )
    np.testing.assert_array_equal(
        pstar_table.csda_range(grid[::100]), pstar_table_python.csda_range(grid[::100])
    )


def test_out_of_range_energies_are_rejected_unless_extrapolating(pstar_cache_root):
    strict = TabulatedStoppingPower.from_dataset(
        MCSQUARE_PSTAR_WATER, materials.WATER, particles.PROTON, pstar_cache_root
    )
    with pytest.raises(ValueError):
        strict.mass_stopping_power(0.25)
    with pytest.raises(ValueError):
        strict.mass_stopping_power(401.0)
    with pytest.raises(ValueError):
        strict.mass_stopping_power(-1.0)
    loose = TabulatedStoppingPower.from_dataset(
        MCSQUARE_PSTAR_WATER,
        materials.WATER,
        particles.PROTON,
        pstar_cache_root,
        allow_extrapolation=True,
    )
    # power-law continuation of the boundary segments
    s = loose.mass_stopping_power([0.25, 401.0])
    assert 0.0 < s[1] < s[0]


def test_tabulated_matches_analytic_within_criterion_above_10_mev(pstar_table) -> None:
    grid = np.linspace(10.0, 400.0, 781)
    analytic = AnalyticStoppingPower(materials.WATER, particles.PROTON, path="numpy")
    ratio = pstar_table.mass_stopping_power(grid) / analytic.mass_stopping_power(grid)
    assert np.max(np.abs(ratio - 1.0)) <= TOL_TABLE_VS_ANALYTIC_ABOVE_10_MEV


def test_tabulated_csda_range_matches_pstar(pstar_table) -> None:
    energies = np.array(sorted(PSTAR_CSDA_RANGE))
    expected = np.array([PSTAR_CSDA_RANGE[e] for e in energies])
    got = pstar_table.csda_range(energies) + PSTAR_CSDA_AT_0_5_MEV
    assert np.all(np.abs(got / expected - 1.0) <= TOL_TABLE_RANGE_VS_PSTAR), got
    # monotonic and consistent with the cumulative column at grid points
    grid = pstar_table.table.energy_mev
    np.testing.assert_allclose(
        pstar_table.csda_range(grid),
        pstar_table.table.csda_range_g_per_cm2,
        rtol=1e-12,
        atol=1e-15,
    )
    assert np.all(np.diff(pstar_table.csda_range(np.linspace(0.5, 400.0, 400))) > 0.0)


def test_analytic_csda_range_matches_pstar_with_floor_residual() -> None:
    analytic = AnalyticStoppingPower(materials.WATER, particles.PROTON, path="numpy")
    energies = np.array(sorted(PSTAR_CSDA_RANGE))
    expected = np.array([PSTAR_CSDA_RANGE[e] for e in energies])
    got = analytic.csda_range(energies) + PSTAR_CSDA_AT_1_MEV
    assert np.all(np.abs(got / expected - 1.0) <= 0.005)


def test_interpolation_is_exact_for_cubic_data() -> None:
    # PCHIP slopes are not the exact derivatives, so exactness holds for the
    # Hermite evaluation given exact slopes; check that path via the shared
    # source directly, then the full chain on a quadratic (exact at nodes).
    py = reference.python_reference("ionmc.physics.tabulated")
    e = np.array([1.0, 2.0, 3.0, 4.0])
    poly = np.polynomial.Polynomial([10.0, -1.0, 0.25, -0.01])
    s_nodes, d_nodes = poly(e), poly.deriv()(e)
    for x in (1.0, 1.3, 2.5, 3.999, 4.0):
        got = py.tabulated_mass_stopping_power(x, e, s_nodes, d_nodes, 4, 2)
        assert got == pytest.approx(poly(x), rel=1e-13)
    table = prepare_stopping_table(e, s_nodes)
    model = TabulatedStoppingPower(table, materials.WATER, particles.PROTON)
    np.testing.assert_allclose(model.mass_stopping_power(e), s_nodes, rtol=1e-14)


def test_partial_segment_quadrature_matches_cumulative_column(pstar_table) -> None:
    # R at a grid point computed through the partial-segment path (t = 1 of the
    # previous segment) must equal the cumulative column to quadrature accuracy.
    # the cumulative column and the on-the-fly query use the same 4-point
    # Gauss-Legendre rule, so at a grid node the query equals the column
    # exactly and the range is continuous across nodes.
    t = pstar_table.table
    np.testing.assert_allclose(
        pstar_table.csda_range(t.energy_mev),
        t.csda_range_g_per_cm2,
        rtol=1e-12,
        atol=1e-15,
    )
    # across a node the range is continuous to the 4-point Gauss-Legendre
    # single-segment error (~3e-7 g/cm^2 = a few nm of water), which is set by
    # the rule, not a coding discontinuity.
    e = 100.0
    left, right = pstar_table.csda_range([e - 1e-6, e + 1e-6])
    assert abs(right - left) < 1e-6


def test_pchip_reproduces_mcsquare_intermediate_rows_from_pstar_nodes(
    pstar_table,
) -> None:
    """The MCsquare authors generated the 0.5 MeV rows from the native PSTAR
    grid with a PCHIP interpolant (decision 0008); rebuilding them from the
    51 native nodes must reproduce the file below 300 MeV to ~1e-8."""
    e = pstar_table.table.energy_mev
    s = pstar_table.table.stopping_mev_cm2_per_g
    native = [x for x in np.arange(0.5, 10.01, 0.5)] + [
        12.5,
        15.0,
        17.5,
        20.0,
        25.0,
        27.5,
        30.0,
    ]
    native += (
        list(np.arange(35.0, 100.01, 5.0))
        + list(np.arange(125.0, 300.01, 25.0))
        + [350.0, 400.0]
    )
    native_e = np.array(sorted(set(np.round(native, 6))))
    mask = np.isin(np.round(e, 6), native_e)
    assert mask.sum() == 51
    rebuilt = TabulatedStoppingPower(
        prepare_stopping_table(e[mask], s[mask]), materials.WATER, particles.PROTON
    )
    below = e < 300.0
    rel = rebuilt.mass_stopping_power(e[below]) / s[below] - 1.0
    assert np.max(np.abs(rel)) < 1e-7


def test_locate_segment_covers_all_indices() -> None:
    py = reference.python_reference("ionmc.physics.tabulated")
    log_e = np.log(np.array([1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0]))
    n, steps = 7, 3
    for value, expected in (
        (0.5, 0),
        (1.0, 0),
        (1.5, 0),
        (2.0, 1),
        (7.9, 2),
        (8.0, 3),
        (63.0, 5),
        (64.0, 5),
        (100.0, 5),
    ):
        assert py.locate_segment(np.log(value), log_e, n, steps) == expected, value


def test_g4_table_differs_from_pstar_like_the_icru90_offset(
    pstar_cache_root, g4_cache_root, pstar_table
):
    g4 = TabulatedStoppingPower.from_dataset(
        MCSQUARE_G4_WATER, materials.WATER, particles.PROTON, g4_cache_root
    )
    energies = np.array([20.0, 50.0, 100.0, 200.0])
    table_ratio = g4.mass_stopping_power(energies) / pstar_table.mass_stopping_power(
        energies
    )
    a49 = AnalyticStoppingPower(materials.WATER, particles.PROTON, path="numpy")
    a90 = AnalyticStoppingPower(materials.WATER_ICRU90, particles.PROTON, path="numpy")
    model_ratio = a90.mass_stopping_power(energies) / a49.mass_stopping_power(energies)
    # both are about -0.5 %; they agree with each other to a few 1e-4
    assert np.all(table_ratio < 1.0)
    np.testing.assert_allclose(table_ratio, model_ratio, atol=1.5e-3)


def test_provenance_records_dataset_and_path(pstar_table) -> None:
    p = pstar_table.provenance()
    assert p["model"] == "tabulated-pchip" and p["path"] == "numpy"
    assert p["table"]["dataset"] == MCSQUARE_PSTAR_WATER.name
    assert p["table"]["sha256"] == MCSQUARE_PSTAR_WATER.sha256
    assert p["range_floor_energy_mev"] == 0.5


@pytest.mark.warp
def test_warp_paths_match_reference(warp_module, pstar_cache_root):
    """Warp CPU (float32) vs the float64 reference for the tabulated layer.

    Tolerances are those of decision 0005 (reference vs float32 Warp), the same
    criteria the analytic layer uses. CUDA is exercised by the host validation.
    """
    ref = TabulatedStoppingPower.from_dataset(
        MCSQUARE_PSTAR_WATER,
        materials.WATER,
        particles.PROTON,
        pstar_cache_root,
        path="python",
    )
    warp_model = TabulatedStoppingPower.from_dataset(
        MCSQUARE_PSTAR_WATER,
        materials.WATER,
        particles.PROTON,
        pstar_cache_root,
        path="warp",
        device="cpu",
    )
    grid = np.linspace(0.5, 400.0, 800)
    np.testing.assert_allclose(
        warp_model.mass_stopping_power(grid),
        ref.mass_stopping_power(grid),
        rtol=1e-5,
        atol=0.0,
    )
    e = np.array([50.0, 100.0, 150.0, 200.0, 250.0])
    np.testing.assert_allclose(
        warp_model.csda_range(e),
        ref.csda_range(e),
        rtol=2e-5,
        atol=0.0,
    )


@pytest.mark.cuda
def test_warp_cuda_matches_cpu(warp_module, cuda_available, pstar_cache_root):
    cpu = TabulatedStoppingPower.from_dataset(
        MCSQUARE_PSTAR_WATER,
        materials.WATER,
        particles.PROTON,
        pstar_cache_root,
        path="warp",
        device="cpu",
    )
    cuda = TabulatedStoppingPower.from_dataset(
        MCSQUARE_PSTAR_WATER,
        materials.WATER,
        particles.PROTON,
        pstar_cache_root,
        path="warp",
        device="cuda:0",
    )
    grid = np.linspace(0.5, 400.0, 800)
    np.testing.assert_allclose(
        cuda.mass_stopping_power(grid),
        cpu.mass_stopping_power(grid),
        rtol=4e-6,
        atol=1e-6,
    )
