import math

import numpy as np

from ionmc.materials import water
from ionmc.physics import projectiles as P
from ionmc.physics import stopping as S

W = water(78.0)
ions = [P.CARBON12] + [getattr(P, n) for n in dir(P) if n.startswith("OXYGEN")][:1]
E = np.array([2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0])
for ion in ions:
    beta2, bg2, tmax = S.kinematics(E, ion)
    beta = np.sqrt(beta2)
    q = S.pierce_blann_charge(beta, ion)
    z = ion.z
    L0 = (
        0.5 * np.log(2 * S.ELECTRON_MASS_MEV * bg2 * tmax / (78e-6) ** 2)
        - beta2
        - 0.5 * S.density_effect(bg2, W)
    )
    sh = -S.shell_correction_over_z(bg2, W)
    bk = S.barkas_term(beta2, q, W)
    bl = S.bloch_term(beta2, q)
    mo = 0.5 * math.pi * S.ALPHA_FS * beta * q
    Lt = L0 + sh + bk + bl + mo
    full = S.bethe_mass_stopping(E, ion, W)
    nq = S.bethe_mass_stopping(E, ion, W, effective_charge=False)
    lsb = S.bethe_mass_stopping(E, ion, W, ls="bare_mixed")
    print(ion.name, "z", z)
    print(" E   q/z  1-q2/z2%  L0    shell%  barkas%  bloch%  mott%  S(noq)/S-1%  LSbare/S-1%")
    for i, e in enumerate(E):
        print(
            f"{e:5.0f} {q[i] / z:.4f} {100 * (1 - (q[i] / z) ** 2):6.2f} {L0[i]:.3f} {100 * sh[i] / Lt[i]:6.2f} {100 * bk[i] / Lt[i]:6.2f} {100 * bl[i] / Lt[i]:6.2f} {100 * mo[i] / Lt[i]:5.2f} {100 * (nq[i] / full[i] - 1):6.2f} {100 * (lsb[i] / full[i] - 1):6.2f}"
        )
    # alternative charges
    vv0 = (
        np.sqrt(E * 1000 * P.PROTON.mass_mev / ion.mass_mev * ion.a / 25.0)
        if hasattr(ion, "mass_mev")
        else None
    )
