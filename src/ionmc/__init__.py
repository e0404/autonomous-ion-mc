"""IonMC: ion-therapy Monte Carlo transport.

The package is developed autonomously under the experiment protocol recorded in
``EXPERIMENT.md`` and ``experiment/v2/``. Public functionality is added task by
task; anything not exported here is not yet a supported public interface.
"""

from ionmc.provenance import code_identity

__version__ = "0.1.0.dev0"

__all__ = ["__version__", "code_identity"]
