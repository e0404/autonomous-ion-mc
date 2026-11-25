"""Minimal deterministic SVG line-plot writer (decision 0045 section 5; decision 0037: SVG by an
internal writer, no plotting dependency).

Fixed canvas, axes with ticks and labels, polyline series, filled bands (for a +-1 sigma band),
and gaps: a non-finite value (NaN or infinity) breaks the path, it is never drawn as a coordinate.
Every number is written with ``format(x, ".6g")``; there are no timestamps, random ids, embedded
fonts or external references; text is XML-escaped; the output ends with one newline and parses
with ``xml.etree.ElementTree``. The output is a pure function of the arguments.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from xml.sax.saxutils import escape

WIDTH = 640
HEIGHT = 400
_LEFT, _RIGHT, _TOP, _BOTTOM = 72, 20, 44, 56
_PALETTE = ("#1f4e79", "#b45f06", "#38761d", "#7f3f98")
_BAND = "#9fc5e8"


def num(x: float) -> str:
    """The one number format of the writer (``format(x, ".6g")``); refuses a non-finite value."""
    if not math.isfinite(x):
        raise ValueError(f"non-finite number {x!r} in an SVG coordinate")
    s = format(float(x), ".6g")
    return "0" if s in ("-0", "-0.0") else s


def nice_ticks(lo: float, hi: float, target: int = 6) -> list[float]:
    """Round tick values covering ``[lo, hi]`` (deterministic; 1, 2, 5 times a power of ten)."""
    if not hi > lo:
        return [lo]
    raw = (hi - lo) / max(target - 1, 1)
    mag = 10.0 ** math.floor(math.log10(raw))
    step = mag * next(m for m in (1.0, 2.0, 5.0, 10.0) if raw <= m * mag * (1 + 1e-12))
    first = math.ceil(lo / step - 1e-9)
    last = math.floor(hi / step + 1e-9)
    return [round(i * step, 12) for i in range(first, last + 1)]


def _segments(xs: Iterable[float], ys: Iterable[float]) -> list[list[tuple[float, float]]]:
    """Maximal runs of points with both coordinates finite (a gap ends a run)."""
    out: list[list[tuple[float, float]]] = []
    run: list[tuple[float, float]] = []
    for x, y in zip(xs, ys, strict=True):
        if math.isfinite(x) and math.isfinite(y):
            run.append((float(x), float(y)))
        elif run:
            out.append(run)
            run = []
    if run:
        out.append(run)
    return out


class Figure:
    """One plot: add series, then :meth:`render`."""

    def __init__(self, title: str, xlabel: str, ylabel: str, note: str = "") -> None:
        self.title, self.xlabel, self.ylabel, self.note = title, xlabel, ylabel, note
        self._lines: list[tuple[str, list[list[tuple[float, float]]]]] = []
        self._bands: list[list[list[tuple[float, float]]]] = []
        self._bands_lo: list[list[list[tuple[float, float]]]] = []
        self._xs: list[float] = []  # every finite abscissa, also where the ordinate is undefined

    def line(self, label: str, xs: Iterable[float], ys: Iterable[float]) -> None:
        xs = list(xs)
        self._xs += [x for x in xs if math.isfinite(x)]
        self._lines.append((label, _segments(xs, ys)))

    def band(self, xs: Iterable[float], lo: Iterable[float], hi: Iterable[float]) -> None:
        """A filled band between ``lo`` and ``hi``; a point with an undefined bound is a gap."""
        xs = list(xs)
        self._xs += [x for x in xs if math.isfinite(x)]
        runs: list[list[tuple[float, float, float]]] = []
        run: list[tuple[float, float, float]] = []
        for x, a, b in zip(xs, lo, hi, strict=True):
            if math.isfinite(x) and math.isfinite(a) and math.isfinite(b):
                run.append((float(x), float(a), float(b)))
            elif run:
                runs.append(run)
                run = []
        if run:
            runs.append(run)
        self._bands.append([[(x, b) for x, _a, b in r] for r in runs])
        self._bands_lo.append([[(x, a) for x, a, _b in r] for r in runs])

    def _ranges(self) -> tuple[float, float, float, float]:
        pts = [p for _l, segs in self._lines for s in segs for p in s]
        pts += [p for b in (*self._bands, *self._bands_lo) for s in b for p in s]
        if not pts and not self._xs:
            return 0.0, 1.0, 0.0, 1.0
        x0, x1 = min([*self._xs, *(p[0] for p in pts)]), max([*self._xs, *(p[0] for p in pts)])
        y0, y1 = (min(p[1] for p in pts), max(p[1] for p in pts)) if pts else (0.0, 1.0)
        if not x1 > x0:
            x0, x1 = x0 - 0.5, x1 + 0.5
        if not y1 > y0:
            y0, y1 = y0 - 0.5, y1 + 0.5
        return x0, x1, y0, y1

    def render(self) -> str:
        x0, x1, y0, y1 = self._ranges()
        xt, yt = nice_ticks(x0, x1), nice_ticks(y0, y1)
        if xt and len(xt) > 1:
            x0, x1 = min(x0, xt[0]), max(x1, xt[-1])
        if yt and len(yt) > 1:
            y0, y1 = min(y0, yt[0]), max(y1, yt[-1])
        pw, ph = WIDTH - _LEFT - _RIGHT, HEIGHT - _TOP - _BOTTOM

        def px(x: float) -> str:
            return num(_LEFT + (x - x0) / (x1 - x0) * pw)

        def py(y: float) -> str:
            return num(_TOP + ph - (y - y0) / (y1 - y0) * ph)

        o: list[str] = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
            f'viewBox="0 0 {WIDTH} {HEIGHT}" font-family="sans-serif" font-size="11">',
            f"<title>{escape(self.title)}</title>",
            f'<rect x="0" y="0" width="{WIDTH}" height="{HEIGHT}" fill="#ffffff"/>',
            f'<text x="{num(_LEFT)}" y="22" font-size="14">{escape(self.title)}</text>',
        ]
        if self.note:
            o.append(f'<text x="{num(_LEFT)}" y="37" fill="#444444">{escape(self.note)}</text>')
        for t in xt:
            o.append(
                f'<line x1="{px(t)}" y1="{num(_TOP + ph)}" x2="{px(t)}" y2="{num(_TOP + ph + 4)}" '
                'stroke="#000000"/>'
            )
            o.append(
                f'<text x="{px(t)}" y="{num(_TOP + ph + 17)}" text-anchor="middle">'
                f"{escape(num(t))}</text>"
            )
        for t in yt:
            o.append(
                f'<line x1="{num(_LEFT - 4)}" y1="{py(t)}" x2="{num(_LEFT)}" y2="{py(t)}" '
                'stroke="#000000"/>'
            )
            o.append(
                f'<line x1="{num(_LEFT)}" y1="{py(t)}" x2="{num(_LEFT + pw)}" y2="{py(t)}" '
                'stroke="#e0e0e0"/>'
            )
            o.append(
                f'<text x="{num(_LEFT - 7)}" y="{py(t)}" text-anchor="end" '
                f'dominant-baseline="middle">{escape(num(t))}</text>'
            )
        o.append(
            f'<rect x="{num(_LEFT)}" y="{num(_TOP)}" width="{num(pw)}" height="{num(ph)}" '
            'fill="none" stroke="#000000"/>'
        )
        o.append(
            f'<text x="{num(_LEFT + pw / 2)}" y="{num(HEIGHT - 12)}" text-anchor="middle">'
            f"{escape(self.xlabel)}</text>"
        )
        cy = _TOP + ph / 2
        o.append(
            f'<text x="14" y="{num(cy)}" text-anchor="middle" '
            f'transform="rotate(-90 14 {num(cy)})">{escape(self.ylabel)}</text>'
        )
        for hi_runs, lo_runs in zip(self._bands, self._bands_lo, strict=True):
            for hi, lo in zip(hi_runs, lo_runs, strict=True):
                poly = [*hi, *reversed(lo)]  # one polygon per run of defined voxels
                pts = " ".join(f"{px(x)},{py(y)}" for x, y in poly)
                o.append(
                    f'<polygon points="{pts}" fill="{_BAND}" fill-opacity="0.6" stroke="none"/>'
                )
        for i, (label, segs) in enumerate(self._lines):
            colour = _PALETTE[i % len(_PALETTE)]
            for s in segs:
                if len(s) == 1:
                    o.append(
                        f'<circle cx="{px(s[0][0])}" cy="{py(s[0][1])}" r="2" fill="{colour}"/>'
                    )
                else:
                    pts = " ".join(f"{px(x)},{py(y)}" for x, y in s)
                    o.append(
                        f'<polyline points="{pts}" fill="none" stroke="{colour}" '
                        'stroke-width="1.5"/>'
                    )
            ly = _TOP + 14 + 14 * i
            o.append(
                f'<line x1="{num(_LEFT + pw - 150)}" y1="{num(ly - 4)}" '
                f'x2="{num(_LEFT + pw - 132)}" y2="{num(ly - 4)}" stroke="{colour}" '
                'stroke-width="1.5"/>'
            )
            o.append(f'<text x="{num(_LEFT + pw - 128)}" y="{num(ly)}">{escape(label)}</text>')
        if not any(segs for _l, segs in self._lines):
            o.append(
                f'<text x="{num(_LEFT + pw / 2)}" y="{num(_TOP + ph / 2)}" text-anchor="middle">'
                "no defined values</text>"
            )
        o.append("</svg>")
        return "\n".join(o) + "\n"
