"""ENDF-6 structural identity / coverage verifier and record writer (V3-005D, exploratory role).

Parses, read-only and with plain code (no ``exec``, no external ENDF library), the parts of an
ENDF-6 neutron evaluation that identify it and bound its energy coverage, and writes a
deterministic record (``acquisition.json`` plus ``ACQUISITION.md`` rendered from that JSON).

Git policy (decision 0041, ``validation/results/elastic_attribution/README.md``): the record holds
hashes, headers, grid metadata (incident energies of the tabulation), counts, flags, RATIOS and
aggregate statistics (``<1 - mu_CM>``, exceedance probabilities, angles of local minima). It never
holds an absolute cross-section value, a Legendre coefficient or a per-point distribution; the
sigma_el comparison is stored as a ratio to the committed BGG/Barashenkov formula only.

Usage (repository root)::

    python -m validation.scripts.reference_data.endf_verify --zip ZIP \\
        --member ENDF-B-VIII.0_neutrons/n-006_C_012.endf --member ... \\
        --changelog ENDF-B-VIII.0_neutrons/CHANGELOG.txt --acquisition-meta META.json --out DIR

Interpolation laws (ENDF-6 INT): 1 histogram, 2 lin-lin, 3 lin-y/log-x, 4 log-y/lin-x,
5 log-log. The per-energy tabulated f(mu) supports INT 2 and 4 (exact segment integrals,
including partial segments up to an angular threshold); any other law raises ``ValueError``.
"""

from __future__ import annotations

import argparse
import bisect
import hashlib
import json
import math
import re
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Any

from ionmc.nuclear.bgg import barashenkov_elastic_mb

SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
MODULE_NAME = "validation.scripts.reference_data.endf_verify"
FORWARD_MINIMA_ANGLE_DEG = 130.0
EXCEEDANCE_ANGLES_DEG = (10.0, 20.0)
SIGMA_RATIO_ENERGIES_MEV = (100.0, 150.0)
BARN_TO_MB = 1.0e3
LOG_EXP_SWITCH = 1.0e-6  # |k dx| below this: the exponential segment is treated as linear


# --------------------------------------------------------------------------------------------
# ENDF-6 record layer
# --------------------------------------------------------------------------------------------
def endf_float(field: str) -> float:
    """ENDF-6 real: ``1.0+5`` = 1e5, ``-1.5-3`` = -1.5e-3; blank = 0; plain ``1.5e3`` accepted."""
    s = field.strip()
    if not s:
        return 0.0
    try:
        v = float(s)
    except ValueError:
        v = float(re.sub(r"(?<=[0-9.])([+-])(?=[0-9])", r"e\1", s))
    if not math.isfinite(v):
        raise ValueError(f"non-finite ENDF real {field!r}")
    return v


def endf_int(field: str) -> int:
    s = field.strip()
    return int(s) if s else 0


def line_ids(line: str) -> tuple[int, int, int]:
    """(MAT, MF, MT) from columns 67-80."""
    return endf_int(line[66:70]), endf_int(line[70:72]), endf_int(line[72:75])


def split_sections(lines: list[str]) -> dict[tuple[int, int, int], list[str]]:
    """Map (MAT, MF, MT) to the 66-column data lines of the section (SEND/FEND/MEND/TEND dropped).

    The first line is the tape label. Any other line whose identifiers cannot be parsed raises.
    """
    sections: dict[tuple[int, int, int], list[str]] = {}
    for n, ln in enumerate(lines[1:], start=2):
        if not ln.strip():
            continue
        try:
            mat, mf, mt = line_ids(ln)
        except ValueError as exc:
            raise ValueError(f"line {n}: unparsable MAT/MF/MT: {exc}") from exc
        if mat <= 0 or mf == 0 or mt == 0:
            continue
        sections.setdefault((mat, mf, mt), []).append(ln[:66].ljust(66))
    return sections


def cont(line: str) -> tuple[float, float, int, int, int, int]:
    return (
        endf_float(line[0:11]),
        endf_float(line[11:22]),
        endf_int(line[22:33]),
        endf_int(line[33:44]),
        endf_int(line[44:55]),
        endf_int(line[55:66]),
    )


def read_values(lines: list[str], start: int, count: int) -> tuple[list[float], int]:
    """``count`` reals from 6-per-line records starting at line ``start``; returns next index."""
    nlines = (count + 5) // 6
    if start + nlines > len(lines):
        raise ValueError("record runs past the end of the section")
    vals: list[float] = []
    for ln in lines[start : start + nlines]:
        for k in range(6):
            field = ln[11 * k : 11 * k + 11]
            if len(vals) < count:
                vals.append(endf_float(field))
            elif field.strip():
                raise ValueError("non-blank field beyond the declared record length")
    return vals, start + nlines


def read_interp(lines: list[str], start: int, nr: int) -> tuple[list[tuple[int, int]], int]:
    """Interpolation table (NBT, INT) pairs, 3 per line."""
    nlines = (nr + 2) // 3
    if start + nlines > len(lines):
        raise ValueError("interpolation table runs past the end of the section")
    pairs: list[tuple[int, int]] = []
    for j in range(nr):
        ln = lines[start + j // 3]
        o = 22 * (j % 3)
        pairs.append((endf_int(ln[o : o + 11]), endf_int(ln[o + 11 : o + 22])))
    return pairs, start + nlines


class Tab1:
    """ENDF TAB1 record: head constants, interpolation table, x and y arrays."""

    def __init__(
        self, head: tuple[float, float, int, int], interp: list[tuple[int, int]],
        x: list[float], y: list[float],
    ) -> None:
        self.c1, self.c2, self.l1, self.l2 = head
        self.interp = interp
        self.x = x
        self.y = y
        if len(x) != len(y) or len(x) < 2:
            raise ValueError("TAB1 needs at least two points")
        if any(b < a for a, b in zip(x, x[1:], strict=False)):
            raise ValueError("TAB1 abscissae are not non-decreasing")
        if not interp or interp[-1][0] != len(x):
            raise ValueError("TAB1 interpolation table does not end at the last point")

    def segment_laws(self) -> list[int]:
        """INT law of each of the NP-1 segments (the region containing the segment's end point)."""
        laws = []
        for j in range(len(self.x) - 1):
            end_point = j + 2
            for nbt, law in self.interp:
                if nbt >= end_point:
                    laws.append(law)
                    break
        return laws


def read_tab1(lines: list[str], start: int) -> tuple[Tab1, int]:
    c1, c2, l1, l2, nr, np_ = cont(lines[start])
    interp, i = read_interp(lines, start + 1, nr)
    vals, i = read_values(lines, i, 2 * np_)
    return Tab1((c1, c2, l1, l2), interp, vals[0::2], vals[1::2]), i


def read_tab2(lines: list[str], start: int) -> tuple[int, int, list[tuple[int, int]], int]:
    _c1, _c2, _l1, _l2, nr, ne = cont(lines[start])
    interp, i = read_interp(lines, start + 1, nr)
    return nr, ne, interp, i


# --------------------------------------------------------------------------------------------
# Interpolation and exact moments
# --------------------------------------------------------------------------------------------
def interpolate(tab: Tab1, x: float) -> float:
    """TAB1 value at ``x`` with the table's own (per-region) interpolation laws 1-5."""
    xs, ys = tab.x, tab.y
    if x < xs[0] or x > xs[-1]:
        raise ValueError("abscissa outside the tabulated range")
    if x == xs[-1]:
        return ys[-1]
    j = bisect.bisect_right(xs, x) - 1
    x1, x2, y1, y2 = xs[j], xs[j + 1], ys[j], ys[j + 1]
    law = tab.segment_laws()[j]
    if law == 1:
        return y1
    if law == 2:
        return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
    if law == 3:
        return y1 + (y2 - y1) * math.log(x / x1) / math.log(x2 / x1)
    if law in (4, 5):
        if y1 <= 0.0 or y2 <= 0.0:
            raise ValueError("log interpolation of a non-positive ordinate")
        if law == 4:
            return y1 * math.exp(math.log(y2 / y1) * (x - x1) / (x2 - x1))
        return y1 * math.exp(math.log(y2 / y1) * math.log(x / x1) / math.log(x2 / x1))
    raise ValueError(f"unsupported interpolation law INT={law}")


def segment_moments(
    x1: float, x2: float, y1: float, y2: float, law: int, a: float, b: float
) -> tuple[float, float]:
    """(int y dx, int x y dx) over [a, b] within the segment [x1, x2]; INT 2 and 4, exact."""
    if law not in (2, 4):
        raise ValueError(f"unsupported interpolation law INT={law} for f(mu)")
    if b <= a or x2 <= x1:
        return 0.0, 0.0
    dx = x2 - x1
    ta, tb = a - x1, b - x1
    k = 0.0
    if law == 4:
        if y1 < 0.0 or y2 < 0.0 or (y1 == 0.0) != (y2 == 0.0):
            raise ValueError("INT 4 segment with a zero/negative ordinate: log undefined")
        if y1 == 0.0:
            return 0.0, 0.0
        k = math.log(y2 / y1) / dx
    if law == 2 or abs(k * dx) < LOG_EXP_SWITCH:
        s = (y2 - y1) / dx
        m0 = y1 * (tb - ta) + s * (tb**2 - ta**2) / 2.0
        mt = y1 * (tb**2 - ta**2) / 2.0 + s * (tb**3 - ta**3) / 3.0
    else:
        ea, eb = math.exp(k * ta), math.exp(k * tb)
        m0 = y1 * (eb - ea) / k
        mt = y1 * (eb * (tb / k - 1.0 / k**2) - ea * (ta / k - 1.0 / k**2))
    return m0, mt + x1 * m0


def tab_moments(tab: Tab1, a: float, b: float) -> tuple[float, float]:
    """(int f dmu, int mu f dmu) over [a, b] with each segment's own law, partial ends exact."""
    m0 = m1 = 0.0
    for j, law in enumerate(tab.segment_laws()):
        x1, x2 = tab.x[j], tab.x[j + 1]
        lo, hi = max(a, x1), min(b, x2)
        if hi <= lo:
            continue
        d0, d1 = segment_moments(x1, x2, tab.y[j], tab.y[j + 1], law, lo, hi)
        m0 += d0
        m1 += d1
    return m0, m1


def _r(x: float) -> float:
    """Fixed significant-digit rounding of every derived statistic (deterministic output)."""
    return float(f"{x:.8g}")


def energy_statistics(tab: Tab1) -> dict[str, Any]:
    """Aggregate statistics of one tabulated f(mu) (CM): norm, <1-mu>, P(theta>10/20 deg), minima."""
    lo, hi = tab.x[0], tab.x[-1]
    norm, first = tab_moments(tab, lo, hi)
    if norm <= 0.0:
        raise ValueError("non-positive normalisation integral")
    stats: dict[str, Any] = {
        "E_eV": tab.c2,
        "NP": len(tab.x),
        "interp_laws": sorted(set(tab.segment_laws())),
        "mu_min": lo,
        "mu_max": hi,
        "norm_integral": _r(norm),
        "mean_one_minus_mu": _r((norm - first) / norm),
    }
    for theta in EXCEEDANCE_ANGLES_DEG:
        mu0 = math.cos(math.radians(theta))
        below, _ = tab_moments(tab, lo, min(mu0, hi))  # mu < cos(theta) <=> angle > theta
        stats[f"p_theta_gt_{theta:g}deg"] = _r(below / norm)
    minima = []
    for j in range(1, len(tab.x) - 1):
        f0, f1, f2 = tab.y[j - 1], tab.y[j], tab.y[j + 1]
        if f1 < f0 and f1 < f2:
            ang = math.degrees(math.acos(max(-1.0, min(1.0, tab.x[j]))))
            if ang < FORWARD_MINIMA_ANGLE_DEG:
                minima.append(round(ang, 3))
    stats[f"minima_deg_forward_of_{FORWARD_MINIMA_ANGLE_DEG:g}"] = sorted(minima)
    return stats


# --------------------------------------------------------------------------------------------
# Section parsers
# --------------------------------------------------------------------------------------------
def parse_mf1_mt451(sec: list[str]) -> dict[str, Any]:
    za, awr, lrp, lfi, nlib, nmod = cont(sec[0])
    elis, sta, lis, liso, _, nfor = cont(sec[1])
    awi, emax, lrel, _, nsub, nver = cont(sec[2])
    temp, _, ldrv, _, nwd, nxc = cont(sec[3])
    comments = [x.rstrip() for x in sec[4 : 4 + nwd]]
    directory = [cont(x) for x in sec[4 + nwd : 4 + nwd + nxc]]
    if len(directory) != nxc or 4 + nwd + nxc > len(sec):
        raise ValueError("MF1/MT451 directory is shorter than NXC")
    dir_pairs = sorted({(d[2], d[3]) for d in directory})
    return {
        "ZA": za, "AWR": awr, "LRP": lrp, "LFI": lfi, "NLIB": nlib, "NMOD": nmod,
        "ELIS": elis, "STA": sta, "LIS": lis, "LISO": liso, "NFOR": nfor,
        "AWI": awi, "EMAX_eV": emax, "LREL": lrel, "NSUB": nsub, "NVER": nver,
        "TEMP": temp, "LDRV": ldrv, "NWD": nwd, "NXC": nxc,
        "comment_lines_1_to_8": [x[:66].rstrip() for x in comments[:8]],
        "comment_lines_mentioning_retained": [
            {"comment_line": j + 1, "text": x[:66].rstrip()}
            for j, x in enumerate(comments)
            if "retain" in x.lower()
        ],
        "directory_mf_mt": [list(p) for p in dir_pairs],
        "directory_has_mf4_mt2": (4, 2) in dir_pairs,
        "directory_has_mf6_mt2": (6, 2) in dir_pairs,
    }


def parse_mf2_mt151(sec: list[str]) -> dict[str, Any]:
    # HEAD, then the isotope CONT (ZAI, ABN, 0, LFW, NER, 0), then the range CONT.
    el, eh, lru, lrf, nro, naps = cont(sec[2])
    return {"EL_eV": el, "EH_eV": eh, "LRU": lru, "LRF": lrf, "NRO": nro, "NAPS": naps}


def parse_mf3(sec: list[str]) -> tuple[Tab1, dict[str, Any]]:
    tab, end = read_tab1(sec, 1)
    if end != len(sec):
        raise ValueError("MF3 section has trailing records")
    meta = {
        "NR": len(tab.interp), "NP": len(tab.x),
        "E_min_eV": tab.x[0], "E_max_eV": tab.x[-1],
        "interp": [list(p) for p in tab.interp],
        "n_points_ge_20MeV": sum(1 for e in tab.x if e >= 20.0e6),
    }
    return tab, meta


def parse_mf4_mt2(sec: list[str]) -> tuple[dict[str, Any], list[Tab1]]:
    _za, _awr, _l1, ltt, _n1, _n2 = cont(sec[0])
    _c1, _c2, li, lct, _n1, _n2 = cont(sec[1])
    out: dict[str, Any] = {"LTT": ltt, "LI": li, "LCT": lct, "lines": len(sec)}
    if ltt not in (0, 1, 2, 3):
        raise ValueError(f"unsupported LTT={ltt}")
    leg_e: list[float] = []
    nl_max = 0
    tabs: list[Tab1] = []
    i = 2
    if ltt in (1, 3):
        _nr, ne, _interp, i = read_tab2(sec, i)
        for _ in range(ne):
            _t, e, _lt, _l2, nl, _n2 = cont(sec[i])
            _vals, i = read_values(sec, i + 1, nl)  # coefficients are read, never recorded
            leg_e.append(e)
            nl_max = max(nl_max, nl)
    if ltt in (2, 3):
        _nr, ne, _interp, i = read_tab2(sec, i)
        for _ in range(ne):
            tab, i = read_tab1(sec, i)
            tabs.append(tab)
    if i != len(sec):
        raise ValueError("MF4/MT2 section has trailing or missing records")
    tab_e = [t.c2 for t in tabs]
    out["legendre"] = {
        "NE": len(leg_e),
        "E_min_eV": leg_e[0] if leg_e else None,
        "E_max_eV": leg_e[-1] if leg_e else None,
        "NL_max": nl_max,
        "incident_energies_ge_20MeV_eV": [e for e in leg_e if e >= 20.0e6],
    }
    out["tabulated"] = {
        "NE": len(tabs),
        "E_min_eV": tab_e[0] if tab_e else None,
        "E_max_eV": tab_e[-1] if tab_e else None,
        "NP_max": max((len(t.x) for t in tabs), default=None),
        "interp_laws_present": sorted({law for t in tabs for law in t.segment_laws()}),
        "incident_energies_eV": tab_e,
    }
    return out, tabs


def parse_mf34_mt2(sec: list[str], mat: int) -> dict[str, Any]:
    """MF34/MT2: NMT1 blocks, each NL sub-sections (L, L1, LCT, NI) with NI LIST records."""
    _za, _awr, _l1, ltt, _n1, nmt1 = cont(sec[0])
    out: dict[str, Any] = {"LTT": ltt, "NMT1": nmt1, "blocks": []}
    i = 1
    for _ in range(nmt1):
        _c1, _c2, mat1, mt1, nl, nl1 = cont(sec[i])
        i += 1
        block: dict[str, Any] = {"MAT1": mat1, "MT1": mt1, "NL": nl, "NL1": nl1, "sub": []}
        same = mt1 == 2 and mat1 in (0, mat)  # same-reaction block: upper triangle of (L, L1)
        n_sections = nl * (nl + 1) // 2 if same else nl * nl1
        for _ in range(n_sections):
            _c1, _c2, lval, l1val, lct, ni = cont(sec[i])
            i += 1
            e_all: list[float] = []
            lbs: list[int] = []
            for _ in range(ni):
                _a, _b, ls, lb, nt, ne = cont(sec[i])
                vals, i = read_values(sec, i + 1, nt)
                lbs.append(lb)
                e_all.extend(vals[: 2 * ne : 2] if lb in (0, 1, 2, 3, 4) else vals[:ne])
            block["sub"].append(
                {
                    "L": lval, "L1": l1val, "LCT": lct, "NI": ni, "LB": sorted(set(lbs)),
                    "E_min_eV": min(e_all) if e_all else None,
                    "E_max_eV": max(e_all) if e_all else None,
                }
            )
        out["blocks"].append(block)
    if i != len(sec):
        raise ValueError("MF34/MT2 section has trailing or missing records")
    return out


def sigma_el_ratios(
    tab: Tab1, z: int, energies_mev: tuple[float, ...] = SIGMA_RATIO_ENERGIES_MEV
) -> dict[str, float]:
    """MF3/MT2 (own interpolation law) over ``barashenkov_elastic_mb(E, Z)``; ratios only."""
    ratios = {}
    for e in energies_mev:
        sigma_mb = interpolate(tab, e * 1.0e6) * BARN_TO_MB  # MF3 cross sections are in barns
        ratios[f"{e:g}MeV"] = _r(sigma_mb / barashenkov_elastic_mb(e, z))
    return ratios


def verify_member(data: bytes, name: str, zip_member: str) -> dict[str, Any]:
    lines = data.decode("ascii").splitlines()
    sections = split_sections(lines)
    mats = sorted({k[0] for k in sections})
    if len(mats) != 1:
        raise ValueError(f"{name}: expected exactly one MAT, found {mats}")
    mat = mats[0]
    rec: dict[str, Any] = {
        "zip_member": zip_member, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(),
        "lines": len(lines), "tape_label": lines[0][:66].rstrip(), "MAT": mat,
    }
    rec["mf1_mt451"] = parse_mf1_mt451(sections[(mat, 1, 451)])
    if (mat, 2, 151) in sections:
        rec["mf2_mt151"] = parse_mf2_mt151(sections[(mat, 2, 151)])
    mt2: Tab1 | None = None
    for mt in (1, 2):
        if (mat, 3, mt) in sections:
            tab, meta = parse_mf3(sections[(mat, 3, mt)])
            rec[f"mf3_mt{mt}"] = meta
            if mt == 2:
                mt2 = tab
    if mt2 is not None:
        z = int(rec["mf1_mt451"]["ZA"]) // 1000
        rec["Z"] = z
        rec["sigma_el_over_barashenkov"] = sigma_el_ratios(mt2, z)
    if (mat, 4, 2) in sections:
        mf4, tabs = parse_mf4_mt2(sections[(mat, 4, 2)])
        mf4["tabulated_statistics"] = [energy_statistics(t) for t in tabs]
        rec["mf4_mt2"] = mf4
    rec["mf6_mt2_present"] = (mat, 6, 2) in sections
    rec["mf34_mt2"] = (
        parse_mf34_mt2(sections[(mat, 34, 2)], mat) if (mat, 34, 2) in sections else None
    )
    rec["mf_present"] = sorted({k[1] for k in sections})
    rec["mt_in_mf3"] = sorted({k[2] for k in sections if k[1] == 3})
    rec["covariance_mf"] = sorted({k[1] for k in sections if 30 <= k[1] <= 40})
    rec["mt_in_mf33"] = sorted({k[2] for k in sections if k[1] == 33})
    rec["mt_in_mf34"] = sorted({k[2] for k in sections if k[1] == 34})
    return rec


def changelog_record(data: bytes, member_basenames: list[str]) -> dict[str, Any]:
    lines = data.decode("utf-8", errors="replace").splitlines()
    quotes = [
        {"line": n, "member": b, "text": ln}
        for n, ln in enumerate(lines, start=1)
        for b in member_basenames
        if b in ln
    ]
    return {
        "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(), "lines": len(lines),
        "lines_mentioning_members": quotes,
    }


# --------------------------------------------------------------------------------------------
# Provenance and CLI
# --------------------------------------------------------------------------------------------
def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, check=True
    ).stdout


def provenance(argv: list[str]) -> dict[str, Any]:
    return {
        "verifier": "validation/scripts/reference_data/endf_verify.py",
        "verifier_sha256": file_sha256(SCRIPT_PATH),
        "git_head": _git("rev-parse", "HEAD").strip(),
        "git_clean": _git("status", "--porcelain").strip() == "",
        "argv": ["python", "-m", MODULE_NAME, *argv],
    }


def build_record(
    zip_path: Path, members: list[str], changelog: str, meta: dict[str, Any], argv: list[str]
) -> dict[str, Any]:
    prov = provenance(argv)
    zip_sha = file_sha256(zip_path)
    acq = meta["acquisition"]
    if zip_sha != acq["sha256"] or zip_path.stat().st_size != acq["bytes"]:
        raise ValueError("source archive does not match the authored acquisition sha256/bytes")
    out: dict[str, Any] = {k: v for k, v in meta.items() if k != "findings"}
    with zipfile.ZipFile(zip_path) as zf:
        acq = dict(acq)
        acq["zip_entries"] = len(zf.infolist())
        out["acquisition"] = acq
        recs = {}
        for m in members:
            recs[Path(m).name] = verify_member(zf.read(m), Path(m).name, m)
        out["changelog"] = {"zip_member": changelog} | changelog_record(
            zf.read(changelog), sorted(recs)
        )
    out["provenance"] = prov | {"zip_sha256": zip_sha}
    out["members"] = recs
    out["findings"] = meta["findings"]
    return out


def dumps(rec: dict[str, Any]) -> str:
    return json.dumps(rec, indent=1, sort_keys=True, ensure_ascii=True) + "\n"


# --------------------------------------------------------------------------------------------
# Markdown (a pure function of the JSON record)
# --------------------------------------------------------------------------------------------
def _mev(e: float) -> str:
    return f"{e / 1.0e6:g}"


def _ene(e: float) -> str:
    """Energy [eV] with a readable unit."""
    if e < 1.0e3:
        return f"{e:g} eV"
    return f"{e / 1.0e3:g} keV" if e < 1.0e6 else f"{e / 1.0e6:g} MeV"


def render_markdown(rec: dict[str, Any]) -> str:
    a, p = rec["acquisition"], rec["provenance"]
    members = rec["members"]
    L: list[str] = []
    L.append("# ENDF/B-VIII.0 neutron sublibrary: acquisition, identity and coverage record "
             "(C-12 MAT 625, O-16 MAT 825)\n")
    L.append("Exploratory acquisition for the V3-005D M3 candidate. **No construction, calibration "
             "or evaluation use.** This page is rendered from `acquisition.json` by "
             "`validation/scripts/reference_data/endf_verify.py`; no absolute cross-section value, "
             "Legendre coefficient or per-point distribution is recorded (decision 0041).\n")
    L.append("## Acquisition\n")
    L.append(f"- URL `{a['url']}`; {a['bytes']:,} bytes; sha256 `{a['sha256']}`; retrieved "
             f"{a['retrieved_at_utc']}; acquisition id `{a['acquisition_id']}`; role "
             f"`{a['role']}`; acquisition code at `{a['acquisition_code_sha'][:8]}`; "
             f"{a['zip_entries']} zip entries.")
    for t in rec.get("operator_authorisation", []):
        L.append(f"- Authorisation: {t}")
    L.append(f"- Citation: {a['citation']}; license basis: {a['license_basis']}.\n")
    L.append("## Provenance of this record\n")
    L.append(f"- Verifier `{p['verifier']}` sha256 `{p['verifier_sha256']}`; git HEAD "
             f"`{p['git_head']}`; working tree clean: {p['git_clean']}.")
    L.append(f"- Command: `{' '.join(p['argv'])}`")
    L.append(f"- Source archive sha256 `{p['zip_sha256']}`.\n")
    L.append("## Identity of the members\n")
    L.append("| member | MAT | ZA | AWR | bytes | sha256 | evaluation (MF1/MT451) |")
    L.append("|---|---|---|---|---|---|---|")
    for name, m in members.items():
        h = m["mf1_mt451"]
        ev = " / ".join(x.strip() for x in h["comment_lines_1_to_8"][:2])
        L.append(f"| `{name}` | {m['MAT']} | {h['ZA']:g} | {h['AWR']:g} | {m['bytes']:,} | "
                 f"`{m['sha256']}` | {ev} |")
    L.append("")
    L.append("| member | NLIB | NVER | NFOR | NSUB | NMOD | LREL | EMAX (MeV) | NWD | NXC |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for name, m in members.items():
        h = m["mf1_mt451"]
        L.append(f"| `{name}` | {h['NLIB']} | {h['NVER']} | {h['NFOR']} | {h['NSUB']} | "
                 f"{h['NMOD']} | {h['LREL']} | {_mev(h['EMAX_eV'])} | {h['NWD']} | {h['NXC']} |")
    L.append("")
    L.append("## Coverage (MF3/MT2 and MF4/MT2)\n")
    L.append("| member | MF3/MT2 range (MeV), points | LTT / LI / LCT | Legendre part | "
             "tabulated part |")
    L.append("|---|---|---|---|---|")
    for name, m in members.items():
        f3, f4 = m["mf3_mt2"], m["mf4_mt2"]
        lg, tb = f4["legendre"], f4["tabulated"]
        L.append(f"| `{name}` | {_ene(f3['E_min_eV'])} to {_ene(f3['E_max_eV'])}, {f3['NP']} | "
                 f"{f4['LTT']} / {f4['LI']} / {f4['LCT']} | {_ene(lg['E_min_eV'])}"
                 f" to {_ene(lg['E_max_eV'])}, NE {lg['NE']}, "
                 f"NL <= {lg['NL_max']} | {_ene(tb['E_min_eV'])} to {_ene(tb['E_max_eV'])}, "
                 f"NE {tb['NE']}, NP <= {tb['NP_max']}, INT {tb['interp_laws_present']} |")
    L.append("")
    for name, m in members.items():
        f4 = m["mf4_mt2"]
        leg = ", ".join(_mev(e) for e in f4["legendre"]["incident_energies_ge_20MeV_eV"])
        tab = ", ".join(_mev(e) for e in f4["tabulated"]["incident_energies_eV"])
        L.append(f"- `{name}` Legendre-part incident energies >= 20 MeV (MeV): {leg}")
        L.append(f"- `{name}` tabulated-part incident energies (MeV): {tab}")
    L.append("")
    L.append("## Per-energy statistics of the tabulated f(mu) (CM; aggregates only)\n")
    for name, m in members.items():
        L.append(f"### `{name}`\n")
        L.append("| E (MeV) | NP | INT | norm | <1 - mu> | P(theta > 10 deg) | "
                 "P(theta > 20 deg) | local minima forward of 130 deg (deg) |")
        L.append("|---|---|---|---|---|---|---|---|")
        for s in m["mf4_mt2"]["tabulated_statistics"]:
            mins = s["minima_deg_forward_of_130"]
            L.append(f"| {_mev(s['E_eV'])} | {s['NP']} | {s['interp_laws']} | "
                     f"{s['norm_integral']:.6f} | {s['mean_one_minus_mu']:.5f} | "
                     f"{s['p_theta_gt_10deg']:.4f} | {s['p_theta_gt_20deg']:.4f} | "
                     f"{', '.join(f'{x:g}' for x in mins) if mins else 'none'} |")
        L.append("")
    L.append("## sigma_el ratio (MF3/MT2 over the committed BGG/Barashenkov elastic formula)\n")
    ens = sorted({k for m in members.values() for k in m["sigma_el_over_barashenkov"]},
                 key=lambda k: float(k[:-3]))
    L.append("| member | Z | " + " | ".join(ens) + " |")
    L.append("|---|---|" + "---|" * len(ens))
    for name, m in members.items():
        r = m["sigma_el_over_barashenkov"]
        L.append(f"| `{name}` | {m['Z']} | " + " | ".join(f"{r[k]:.4f}" for k in ens) + " |")
    L.append("")
    L.append("## Covariance and other sections\n")
    for name, m in members.items():
        f34 = m["mf34_mt2"]
        if f34 is None:
            c = "no MF34/MT2"
        else:
            parts = []
            for b in f34["blocks"]:
                subs = ", ".join(
                    f"L={s['L']}/L1={s['L1']} LB {s['LB']} E {_ene(s['E_min_eV'])}-"
                    f"{_ene(s['E_max_eV'])}" for s in b["sub"]
                )
                parts.append(f"MT1={b['MT1']} NL={b['NL']} NL1={b['NL1']}: {subs}")
            c = f"MF34/MT2 NMT1={f34['NMT1']}; " + "; ".join(parts)
        L.append(f"- `{name}`: MF {m['mf_present']}; covariance MF {m['covariance_mf']}; "
                 f"MF33 MT {m['mt_in_mf33']}; MF34 MT {m['mt_in_mf34']}; MF6/MT2 present: "
                 f"{m['mf6_mt2_present']}; {c}.")
    L.append("")
    L.append("## MF1 comment lines mentioning 'retain' (verbatim)\n")
    for name, m in members.items():
        for q in m["mf1_mt451"]["comment_lines_mentioning_retained"]:
            L.append(f"- `{name}` comment line {q['comment_line']}: `{q['text']}`")
    L.append("")
    cl = rec["changelog"]
    L.append("## CHANGELOG quotes\n")
    L.append(f"`{cl['zip_member']}`: {cl['bytes']:,} bytes, sha256 `{cl['sha256']}`, "
             f"{cl['lines']} lines. Verbatim lines mentioning the members:\n")
    for q in cl["lines_mentioning_members"]:
        L.append(f"- line {q['line']} ({q['member']}): `{q['text']}`")
    L.append("")
    L.append("## Findings\n")
    for f in rec["findings"]:
        L.append(f"- {f}")
    L.append("")
    for name, q in rec.get("evaluation_comment_quotes", {}).items():
        L.append(f"- Evaluation comment, `{name}`: {q}")
    L.append("")
    L.append(f"**Status:** {rec['status']}\n")
    return "\n".join(L)


def main(argv: list[str] | None = None) -> int:
    args_in = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--zip", required=True, type=Path)
    ap.add_argument("--member", action="append", required=True)
    ap.add_argument("--changelog", required=True)
    ap.add_argument("--acquisition-meta", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ns = ap.parse_args(args_in)
    meta = json.loads(ns.acquisition_meta.read_text(encoding="utf-8"))
    rec = build_record(ns.zip, ns.member, ns.changelog, meta, list(args_in))
    ns.out.mkdir(parents=True, exist_ok=True)
    (ns.out / "acquisition.json").write_text(dumps(rec), encoding="utf-8")
    (ns.out / "ACQUISITION.md").write_text(render_markdown(rec), encoding="utf-8")
    print(f"wrote {ns.out}/acquisition.json and ACQUISITION.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
