"""Reader for geant-val experiment-curve JSON exports (EXFOR-derived inelastic cross sections).

The export is a JSON list of records. Per record this reader returns the target, beam particle,
observable, axis names and the x/y arrays with statistical and systematic y errors. Records
whose beam particle is not ``proton`` (the export includes ``C12`` projectiles) are flagged: the
x axis is labelled ``E (MeV)`` for them as well, but whether it is the total or the per-nucleon
energy is not stated in the file, so such records must not be used for proton comparisons
without that being resolved (``energy_caveat``). Registered parser name: ``geant_val_json``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

_CAVEAT = (
    "non-proton projectile; the export labels the energy axis 'E (MeV)' without stating whether "
    "it is total or per-nucleon kinetic energy"
)


class GeantValError(ValueError):
    """Malformed geant-val JSON."""


@dataclass(frozen=True)
class GeantValCurve:
    """One experiment curve: arrays are float64, units as named by ``x_axis``/``y_axis``."""

    record_id: int
    target: str
    beam: str
    observable: str
    x_axis: str
    y_axis: str
    x: NDArray[np.float64]
    y: NDArray[np.float64]
    y_stat_plus: NDArray[np.float64]
    y_stat_minus: NDArray[np.float64]
    y_sys_plus: NDArray[np.float64]
    y_sys_minus: NDArray[np.float64]
    energy_caveat: str | None


def parse_geant_val(text: str) -> list[GeantValCurve]:
    """Parse the JSON text into curves in file order; missing keys or length mismatches raise."""
    data = json.loads(text)
    if not isinstance(data, list):
        raise GeantValError("expected a JSON list of records")
    out = []
    for rec in data:
        try:
            meta, chart = rec["metadata"], rec["chart"]
            arrays = {
                k: np.asarray(chart[k], dtype=np.float64)
                for k in (
                    "xValues",
                    "yValues",
                    "yStatErrorsPlus",
                    "yStatErrorsMinus",
                    "ySysErrorsPlus",
                    "ySysErrorsMinus",
                )
            }
            beam = str(meta["beamParticle"])
            curve = GeantValCurve(
                record_id=int(rec["id"]),
                target=str(meta["targetName"]),
                beam=beam,
                observable=str(meta["observableName"]),
                x_axis=str(chart["xAxisName"]),
                y_axis=str(chart["yAxisName"]),
                x=arrays["xValues"],
                y=arrays["yValues"],
                y_stat_plus=arrays["yStatErrorsPlus"],
                y_stat_minus=arrays["yStatErrorsMinus"],
                y_sys_plus=arrays["ySysErrorsPlus"],
                y_sys_minus=arrays["ySysErrorsMinus"],
                energy_caveat=None if beam.lower() == "proton" else _CAVEAT,
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise GeantValError(f"malformed geant-val record: {exc!r}") from None
        n = curve.x.size
        if any(a.size != n for a in arrays.values()):
            raise GeantValError(f"record {curve.record_id}: array lengths differ")
        out.append(curve)
    return out
