"""Physical constants used by ionmc.

All values are CODATA 2018 recommended values (Tiesinga et al., Rev. Mod.
Phys. 93, 025010 (2021); https://physics.nist.gov/cuu/Constants/), expressed
in the internal unit system documented in :mod:`ionmc.units`
(energy in MeV, length in mm, mass density in g/cm^3).

The constants are plain Python floats so that they can be used unchanged in
reference-Python code and inside Warp functions and kernels (Warp treats a
module-level float referenced from kernel code as a compile-time constant).
"""

from __future__ import annotations

# --- particle rest energies, MeV ------------------------------------------

#: Electron rest energy m_e c^2 [MeV].
ELECTRON_MASS_MEV: float = 0.51099895000

#: Proton rest energy m_p c^2 [MeV].
PROTON_MASS_MEV: float = 938.27208816

#: Atomic mass constant m_u c^2 [MeV] (1 u).
ATOMIC_MASS_UNIT_MEV: float = 931.49410242

# --- dimensionless and atomic constants ------------------------------------

#: Fine-structure constant alpha.
FINE_STRUCTURE_CONSTANT: float = 7.2973525693e-3

#: Avogadro constant [1/mol].
AVOGADRO: float = 6.02214076e23

#: Classical electron radius r_e [cm].
CLASSICAL_ELECTRON_RADIUS_CM: float = 2.8179403262e-13

#: Bethe-formula prefactor K = 4 pi N_A r_e^2 m_e c^2 [MeV cm^2 / mol].
#: Computed from the constants above; equals 0.307075... MeV cm^2/mol
#: (PDG value 0.307 075 MeV mol^-1 cm^2).
BETHE_K_MEV_CM2_PER_MOL: float = (
    4.0
    * 3.141592653589793
    * AVOGADRO
    * CLASSICAL_ELECTRON_RADIUS_CM**2
    * ELECTRON_MASS_MEV
)

# --- unit conversion helpers as constants -----------------------------------

#: Number of eV in one MeV.
EV_PER_MEV: float = 1.0e6

#: Number of mm in one cm.
MM_PER_CM: float = 10.0
