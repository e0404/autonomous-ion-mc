"""V3-005C (F4b): the multiple-Coulomb-scattering power is purely electromagnetic, so hadronic
elastic scattering added by the elastic channel does not double count it."""

from __future__ import annotations

import ast
import importlib
import inspect
import math
from pathlib import Path

import ionmc.physics.em as em_mod
import ionmc.physics.scattering as sc_mod
from ionmc._wpfunc import python_twin
from ionmc.materials import WATER
from ionmc.physics.em import E_S_MEV, FDM_COEFFICIENTS, make_em
from ionmc.physics.scattering import inverse_scattering_length_cm2_per_g

FORBIDDEN = ("nuclear", "elastic", "bgg", "law5", "tripathi", "endf")


def _imports(mod: object) -> set[str]:
    tree = ast.parse(Path(inspect.getsourcefile(mod)).read_text(encoding="utf-8"))  # type: ignore[arg-type]
    out: set[str] = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            out.update(a.name for a in n.names)
        elif isinstance(n, ast.ImportFrom) and n.level == 0 and n.module:
            out.add(n.module)
            out.update(f"{n.module}.{a.name}" for a in n.names)
    return {i for i in out if i.startswith("ionmc")}


def _closure(*roots: str) -> set[str]:
    """Transitive in-package import closure (package ``__init__`` re-exports are not followed)."""
    seen: set[str] = set()
    todo = list(roots)
    while todo:
        name = todo.pop()
        if name in seen:
            continue
        try:
            mod = importlib.import_module(name)
        except ImportError:  # a name imported from a module, not a module
            continue
        seen.add(name)
        if Path(inspect.getsourcefile(mod) or "").name != "__init__.py":
            todo.extend(_imports(mod))
    return seen


# ionmc.physics.em takes a unit constant from ionmc.physics.stopping, which reads the EM stopping
# data (NIST); these ionmc.data modules are EM-only. No ENDF / AME / nuclear module may appear.
ALLOWED_DATA = {"ionmc.data.cache", "ionmc.data.nist_star", "ionmc.data.registry"}


def test_mcs_scattering_power_has_no_hadronic_term() -> None:
    # imports of the MCS modules contain no nuclear / elastic table dependency
    for mod in (em_mod, sc_mod):
        bad = [i for i in _imports(mod) if any(f in i.lower() for f in FORBIDDEN)]
        assert not bad, (mod.__name__, bad)
    closure = _closure("ionmc.physics.scattering", "ionmc.physics.em")
    assert not [m for m in closure if m.startswith("ionmc.nuclear")], closure
    data = {m for m in closure if m.startswith("ionmc.data")}
    assert data <= ALLOWED_DATA, data - ALLOWED_DATA
    # the inputs are pv, p1v1, charge, 1/X_S (Z, A of the material) and the density only
    em = python_twin(make_em)
    names = list(inspect.signature(em.scattering_power_dm).parameters)
    assert names == ["pv_mev", "p1v1_mev", "z", "inv_rho_xs_cm2_g", "rho_g_cm3"]
    # value at a reference energy equals the documented EM formula (Gottschalk 2010)
    e, m_p, rho = 100.0, 938.27208816, WATER.density_g_cm3
    pv = e * (e + 2.0 * m_p) / (e + m_p)
    e1 = 150.0
    p1v1 = e1 * (e1 + 2.0 * m_p) / (e1 + m_p)
    got = em.scattering_power_dm(pv, p1v1, 1.0, inverse_scattering_length_cm2_per_g(WATER), rho)
    c0, c1, c2, c3 = FDM_COEFFICIENTS
    l1, l2 = math.log10(1.0 - (pv / p1v1) ** 2), math.log10(pv)
    fdm = c0 + c1 * l1 + c2 * l2 + c3 * l1 * l2
    alpha, n_a, r_e = 1.0 / 137.035999084, 6.02214076e23, 2.8179403262e-13
    inv_xs = 0.0
    for el, w in WATER._fractions:
        a, z = el.A_g_mol, float(el.Z)
        inv_xs += (
            w
            * alpha
            * n_a
            * r_e**2
            * z
            * z
            / a
            * (2.0 * math.log(33219.0 * (a * z) ** (-1 / 3)) - 1.0)
        )
    want = fdm * (E_S_MEV / pv) ** 2 * rho * inv_xs / 10.0
    assert abs(got / want - 1.0) < 1e-6
