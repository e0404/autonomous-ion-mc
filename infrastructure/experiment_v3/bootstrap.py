"""Historical bootstrap identity; v3 is prepared from the recorded v2 infra SHA.

The v1-overlay bootstrap is intentionally not executable in the v3 runtime.
Reproduction instructions live in experiment/v3/SETUP.md.
"""

IDENTITY = [
    "-c",
    "user.name=Autonomous IonMC Agent",
    "-c",
    "user.email=autonomous-ionmc-agent@users.noreply.github.com",
]
