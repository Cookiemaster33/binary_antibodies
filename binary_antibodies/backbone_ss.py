"""
backbone_ss.py
--------------
Dependency-free secondary-structure assignment from a CA trace.

DSSP is unavailable inside the RFdiffusion3/foundry container and on the GPU
hosts, and RFd3 backbones have no hydrogens, so hydrogen-bond based assignment
is not an option during the design loop. This module implements the P-SEA
CA-trace criteria (Labesse et al., CABIOS 1997), which classify each residue
from CA(i)->CA(i+n) distances and the virtual CA torsion. It is accurate enough
to rank designs by helix content, which is all Stage A needs.
"""

from __future__ import annotations

import numpy as np

# P-SEA distance windows (Angstrom) for d(i, i+2), d(i, i+3), d(i, i+4).
_HELIX_D2 = (5.0, 6.4)
_HELIX_D3 = (4.6, 6.3)
_HELIX_D4 = (5.6, 7.0)
_STRAND_D2 = (6.4, 7.4)
_STRAND_D3 = (9.4, 10.4)
_STRAND_D4 = (11.6, 13.2)

# Virtual CA torsion windows (degrees).
_HELIX_TAU = (30.0, 70.0)
_STRAND_TAU = (-180.0, -125.0)

_MIN_HELIX_RUN = 4
_MIN_STRAND_RUN = 3


def _in_range(value: float, window: tuple[float, float]) -> bool:
    return window[0] <= value <= window[1]


def _ca_torsion(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    b0, b1, b2 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
    m1 = np.cross(n1, b1 / np.linalg.norm(b1))
    x, y = float(np.dot(n1, n2)), float(np.dot(m1, n2))
    return float(np.degrees(np.arctan2(y, x)))


def _enforce_min_run(labels: list[str], code: str, min_run: int) -> None:
    """Drop runs of `code` shorter than min_run (P-SEA post-filter)."""
    n = len(labels)
    i = 0
    while i < n:
        if labels[i] != code:
            i += 1
            continue
        j = i
        while j < n and labels[j] == code:
            j += 1
        if j - i < min_run:
            for k in range(i, j):
                labels[k] = "L"
        i = j


def assign_ss(ca_coords: np.ndarray) -> str:
    """Return a per-residue SS string over {H, E, L} for a CA trace.

    H = alpha helix, E = beta strand / extended, L = loop / coil.
    """
    ca = np.asarray(ca_coords, dtype=float)
    n = len(ca)
    if n < 5:
        return "L" * n

    def dist(i: int, j: int) -> float:
        return float(np.linalg.norm(ca[j] - ca[i]))

    labels = ["L"] * n
    for i in range(n):
        # Helix test: residue i is helical if any 5-residue window covering it
        # satisfies the helix distance + torsion criteria.
        helix = False
        for start in range(max(0, i - 4), min(i + 1, n - 4)):
            if start + 4 >= n:
                break
            d2_ok = _in_range(dist(start, start + 2), _HELIX_D2)
            d3_ok = _in_range(dist(start, start + 3), _HELIX_D3)
            d4_ok = _in_range(dist(start, start + 4), _HELIX_D4)
            tau_ok = _in_range(
                _ca_torsion(ca[start], ca[start + 1], ca[start + 2], ca[start + 3]),
                _HELIX_TAU,
            )
            if (d2_ok and d3_ok and d4_ok) or (tau_ok and d3_ok):
                helix = True
                break
        if helix:
            labels[i] = "H"
            continue

        strand = False
        for start in range(max(0, i - 4), min(i + 1, n - 4)):
            if start + 4 >= n:
                break
            d2_ok = _in_range(dist(start, start + 2), _STRAND_D2)
            d3_ok = _in_range(dist(start, start + 3), _STRAND_D3)
            d4_ok = _in_range(dist(start, start + 4), _STRAND_D4)
            tau_ok = _in_range(
                _ca_torsion(ca[start], ca[start + 1], ca[start + 2], ca[start + 3]),
                _STRAND_TAU,
            )
            if (d2_ok and d3_ok and d4_ok) or (tau_ok and d3_ok):
                strand = True
                break
        if strand:
            labels[i] = "E"

    _enforce_min_run(labels, "H", _MIN_HELIX_RUN)
    _enforce_min_run(labels, "E", _MIN_STRAND_RUN)
    return "".join(labels)


def ss_fractions(ss: str) -> dict[str, float]:
    n = max(1, len(ss))
    return {
        "helix_frac": ss.count("H") / n,
        "strand_frac": ss.count("E") / n,
        "loop_frac": ss.count("L") / n,
    }


def ss_segments(ss: str, code: str, min_len: int = 1) -> list[tuple[int, int]]:
    """Inclusive 0-based (start, end) spans of `code` runs of at least min_len."""
    spans: list[tuple[int, int]] = []
    i = 0
    while i < len(ss):
        if ss[i] != code:
            i += 1
            continue
        j = i
        while j < len(ss) and ss[j] == code:
            j += 1
        if j - i >= min_len:
            spans.append((i, j - 1))
        i = j
    return spans


def summarize(ca_coords: np.ndarray, min_helix_len: int = 6) -> dict[str, object]:
    """Helix/strand/loop content plus helix-topology descriptors for one chain."""
    ss = assign_ss(ca_coords)
    helices = ss_segments(ss, "H", min_helix_len)
    strands = ss_segments(ss, "E", 3)
    helix_lengths = [end - start + 1 for start, end in helices]
    return {
        "ss": ss,
        "length": len(ss),
        **ss_fractions(ss),
        "n_helices": len(helices),
        "n_strands": len(strands),
        "helix_lengths": helix_lengths,
        "longest_helix": max(helix_lengths, default=0),
    }
