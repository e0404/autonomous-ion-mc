"""Integral depth-dose scoring grid (decision ``0009``).

A uniform grid of ``n_bins`` depth bins over ``[0, depth_mm]`` accumulating the
energy deposited per bin (MeV). For an on-axis pencil beam with no lateral
transport this integral depth dose (energy per depth interval) is the
physically meaningful 1-D quantity; volumetric dose needs a lateral extent
introduced with multiple scattering (a later task). The grid is decoupled from
the transport geometry, as ``REQUIREMENTS.md`` requires.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DepthDoseGrid:
    """Uniform depth bins over ``[origin_mm, origin_mm + depth_mm]``.

    ``origin_mm`` (default 0) shifts the grid along the beam axis independently
    of the transport geometry, so the scoring grid's resolution *and* alignment
    are decoupled from the transport grid (decision ``0019``); the default
    reproduces the historical grid anchored at ``z = 0``.
    """

    depth_mm: float
    n_bins: int
    origin_mm: float = 0.0

    def __post_init__(self) -> None:
        if self.depth_mm <= 0.0 or self.n_bins <= 0:
            raise ValueError("depth and bin count must be positive")

    @property
    def bin_width_mm(self) -> float:
        return self.depth_mm / self.n_bins

    @property
    def edges_mm(self) -> np.ndarray:
        return np.linspace(
            self.origin_mm, self.origin_mm + self.depth_mm, self.n_bins + 1
        )

    @property
    def centers_mm(self) -> np.ndarray:
        e = self.edges_mm
        return 0.5 * (e[:-1] + e[1:])

    def empty(self) -> np.ndarray:
        """A zeroed energy-deposition array (MeV per bin)."""
        return np.zeros(self.n_bins, dtype=np.float64)


@dataclass(frozen=True)
class DepthLateralGrid:
    """2-D scoring grid: depth ``z`` in ``[depth_origin_mm, depth_origin_mm +
    depth_mm]`` x lateral ``x`` in ``[lateral_center_mm +/- half_width_mm]``
    (decisions 0011, 0019).

    Energy is accumulated per (depth, lateral-x) bin, summed over the second
    transverse axis ``y`` (a marginal projection). The depth marginal (sum over
    ``x``) is the integral depth dose; the second moment in ``x`` per depth
    slice gives the lateral spread ``sigma_x(z)``. ``depth_origin_mm`` and
    ``lateral_center_mm`` (both default 0) shift the grid independently of the
    transport geometry, decoupling scoring-grid alignment from the transport
    grid; the defaults reproduce the historical grid anchored at ``z = 0`` and
    centred at ``x = 0``.
    """

    depth_mm: float
    n_depth: int
    half_width_mm: float
    n_lateral: int
    depth_origin_mm: float = 0.0
    lateral_center_mm: float = 0.0

    def __post_init__(self) -> None:
        if self.depth_mm <= 0.0 or self.n_depth <= 0:
            raise ValueError("depth and depth-bin count must be positive")
        if self.half_width_mm <= 0.0 or self.n_lateral <= 0:
            raise ValueError("half width and lateral-bin count must be positive")

    @property
    def depth_bin_mm(self) -> float:
        return self.depth_mm / self.n_depth

    @property
    def lateral_bin_mm(self) -> float:
        return 2.0 * self.half_width_mm / self.n_lateral

    @property
    def lateral_lo_mm(self) -> float:
        """Lower lateral edge ``lateral_center_mm - half_width_mm``."""
        return self.lateral_center_mm - self.half_width_mm

    @property
    def depth_centers_mm(self) -> np.ndarray:
        edges = np.linspace(
            self.depth_origin_mm, self.depth_origin_mm + self.depth_mm, self.n_depth + 1
        )
        return 0.5 * (edges[:-1] + edges[1:])

    @property
    def lateral_centers_mm(self) -> np.ndarray:
        edges = np.linspace(
            self.lateral_center_mm - self.half_width_mm,
            self.lateral_center_mm + self.half_width_mm,
            self.n_lateral + 1,
        )
        return 0.5 * (edges[:-1] + edges[1:])

    def empty(self) -> np.ndarray:
        """A zeroed ``(n_depth, n_lateral)`` energy array [MeV]."""
        return np.zeros((self.n_depth, self.n_lateral), dtype=np.float64)

    def depth_dose(self, edep_zx: np.ndarray) -> np.ndarray:
        """Integral depth dose [MeV per depth bin] = sum over the lateral axis."""
        return edep_zx.sum(axis=1)

    def sigma_x_mm(self, edep_zx: np.ndarray) -> np.ndarray:
        """Lateral RMS ``sigma_x`` [mm] per depth slice (nan where no dose).

        ``sigma_x^2 = sum_x x^2 E / sum_x E - (sum_x x E / sum_x E)^2`` using the
        lateral bin centres; symmetric beams have mean ~0.
        """
        x = self.lateral_centers_mm
        w = edep_zx
        total = w.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = (w * x).sum(axis=1) / total
            mean_sq = (w * x * x).sum(axis=1) / total
            var = mean_sq - mean * mean
        sigma = np.sqrt(np.clip(var, 0.0, None))
        sigma[total <= 0.0] = np.nan
        return sigma


@dataclass(frozen=True)
class DoseGrid3D:
    """A lab-frame 3-D scoring grid accumulating deposited energy per voxel [MeV]
    (decision ``0021``).

    An ``Nx*Ny*Nz`` uniform grid with its minimum corner at ``origin_mm`` and
    voxel edge lengths ``spacing_mm``; it is independent of the transport grid
    (its own resolution and alignment). The 3-D voxel-grid transport path deposits
    each step's energy at the step's lab midpoint into the containing voxel. Flat
    index (shared with the deposition code): ``flat = (i*ny + j)*nz + k``.
    """

    shape: tuple[int, int, int]
    origin_mm: tuple[float, float, float] = (0.0, 0.0, 0.0)
    spacing_mm: tuple[float, float, float] = (2.0, 2.0, 2.0)

    def __post_init__(self) -> None:
        if any(int(n) < 1 for n in self.shape):
            raise ValueError("each grid dimension must have at least one voxel")
        if any(s <= 0.0 for s in self.spacing_mm):
            raise ValueError("voxel spacing must be positive on every axis")

    @property
    def nx(self) -> int:
        return int(self.shape[0])

    @property
    def ny(self) -> int:
        return int(self.shape[1])

    @property
    def nz(self) -> int:
        return int(self.shape[2])

    @property
    def n_voxels(self) -> int:
        return self.nx * self.ny * self.nz

    @property
    def voxel_volume_cm3(self) -> float:
        sx, sy, sz = self.spacing_mm
        return (sx * sy * sz) * 1.0e-3  # mm^3 -> cm^3

    def empty(self) -> np.ndarray:
        """A zeroed ``(nx, ny, nz)`` energy array [MeV]."""
        return np.zeros((self.nx, self.ny, self.nz), dtype=np.float64)

    def dose_gy(self, edep_mev: np.ndarray, density_g_per_cm3: float) -> np.ndarray:
        """Convert deposited energy [MeV] to absorbed dose [Gy] for a uniform
        voxel mass ``density * voxel_volume`` (1 MeV = 1.602176634e-13 J,
        1 g = 1e-3 kg)."""
        mass_kg = density_g_per_cm3 * self.voxel_volume_cm3 * 1.0e-3
        return edep_mev * (1.602176634e-13 / mass_kg)

    def axis_marginals(self, edep: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(z_centers_mm, depth_profile)`` -- energy summed over x and y
        per z-plane, with z-voxel centres in the lab frame."""
        oz = self.origin_mm[2]
        sz = self.spacing_mm[2]
        centers = oz + (np.arange(self.nz) + 0.5) * sz
        return centers, edep.sum(axis=(0, 1))

    def let_d_kev_um(
        self,
        let_num_mev_per_mm: np.ndarray,
        edep_mev: np.ndarray,
        min_dose_frac: float = 0.0,
    ) -> np.ndarray:
        """Dose-averaged LET ``LET_d = num / den`` [keV/um] (decision ``0023``).

        ``let_num_mev_per_mm`` is the per-voxel ``Sum eps_i * L_i`` numerator
        (MeV/mm, numerically keV/um since ``1 MeV/mm == 1 keV/um``) and
        ``edep_mev`` the co-registered dose energy ``Sum eps_i`` (the
        denominator). LET_d is returned only where the voxel dose exceeds
        ``min_dose_frac`` of the peak dose (a low-dose reporting mask, because the
        distal falloff where LET_d is highest is exactly where the dose -- and so
        the ratio's statistics -- is weakest); masked and zero-dose voxels are 0.
        """
        num = np.asarray(let_num_mev_per_mm, dtype=np.float64)
        den = np.asarray(edep_mev, dtype=np.float64)
        peak = float(den.max()) if den.size else 0.0
        threshold = max(0.0, min_dose_frac) * peak
        out = np.zeros_like(den)
        mask = den > threshold
        if threshold <= 0.0:
            mask = den > 0.0
        out[mask] = num[mask] / den[mask]
        return out


@dataclass(frozen=True)
class FluenceSpectrum:
    """A 1-D proton fluence-vs-energy histogram over a single scoring region
    (decision ``0026``). Linear energy bins over ``[e_lo_mev, e_hi_mev]``; the
    scored array is the raw track-length histogram ``counts[k] = Sum w_i*l_i`` [mm]
    (weighted step length ``l = s``) binned by the step-mean energy ``E_mid``. The
    default 160 bins over 0-160 MeV suit a therapeutic proton beam in water."""

    n_bins: int = 160
    e_lo_mev: float = 0.0
    e_hi_mev: float = 160.0

    def __post_init__(self) -> None:
        if self.n_bins <= 0:
            raise ValueError("need at least one energy bin")
        if self.e_hi_mev <= self.e_lo_mev:
            raise ValueError("e_hi_mev must exceed e_lo_mev")

    @property
    def bin_width_mev(self) -> float:
        return (self.e_hi_mev - self.e_lo_mev) / self.n_bins

    @property
    def edges_mev(self) -> np.ndarray:
        return np.linspace(self.e_lo_mev, self.e_hi_mev, self.n_bins + 1)

    @property
    def centers_mev(self) -> np.ndarray:
        e = self.edges_mev
        return 0.5 * (e[:-1] + e[1:])

    def empty(self) -> np.ndarray:
        """A zeroed raw track-length histogram [mm per bin]."""
        return np.zeros(self.n_bins, dtype=np.float64)

    def differential_fluence(
        self, counts_mm: np.ndarray, volume_mm3: float, n_histories: int = 1
    ) -> np.ndarray:
        """Differential fluence ``Phi(E)`` [protons cm^-2 MeV^-1] from the raw
        track-length histogram: ``counts * 100 / (V_mm3 * dE)`` (100 converts mm
        track length over mm^3 volume to cm^-2), divided by ``n_histories`` for the
        per-primary fluence."""
        scale = 100.0 / (volume_mm3 * self.bin_width_mev * max(1, n_histories))
        return np.asarray(counts_mm, dtype=np.float64) * scale

    def postprocess_lookup(
        self, counts_mm: np.ndarray, lookup_table: np.ndarray
    ) -> float:
        """Offline post-processing of the scored spectrum through a per-bin lookup
        table: ``Sum_k counts[k] * lookup_table[k]`` (decision 0026). With the same
        table and bin rule the on-the-fly accumulator reproduces this to round-off.
        """
        return float(
            np.dot(
                np.asarray(counts_mm, dtype=np.float64),
                np.asarray(lookup_table, dtype=np.float64),
            )
        )
