"""ionmc: research-grade ion-therapy Monte Carlo (Python-first, Warp accelerated).

Implemented so far (Stage 0, task DEV-002):

* :mod:`ionmc.units`, :mod:`ionmc.constants` - unit conventions and constants;
* :mod:`ionmc.materials`, :mod:`ionmc.particles` - materials and projectiles;
* :mod:`ionmc.physics.stopping` - shared-source analytical stopping power;
* :mod:`ionmc.stopping_power` - Python API evaluating it on the reference
  Python, numpy and Warp execution paths;
* :mod:`ionmc.rng` - pure-Python mirror of Warp's random-number generator;
* :mod:`ionmc.backend` - the shared-source execution mechanism.
"""

from __future__ import annotations

from ionmc.materials import WATER, WATER_ICRU49, WATER_ICRU90, Material
from ionmc.particles import PROTON, Particle
from ionmc.stopping_power import AnalyticStoppingPower, Corrections, available_paths

__version__ = "0.1.0.dev0"

__all__ = [
    "PROTON",
    "WATER",
    "WATER_ICRU49",
    "WATER_ICRU90",
    "AnalyticStoppingPower",
    "Corrections",
    "Material",
    "Particle",
    "__version__",
    "available_paths",
]
