"""
sequences.py
------------
Flexible linker sequence generation and composition utilities.

The gold standard for antibody engineering flexible linkers is the (G4S)n
repeat — glycine-rich for flexibility, serine for solubility and to prevent
aggregation. This module generates these sequences and variants.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

# Physical properties per residue (approximate)
_MW_PER_RESIDUE = {
    "G": 57.05,
    "S": 87.08,
    "A": 71.08,
    "P": 97.12,
    "T": 101.10,
    "E": 129.12,
    "K": 128.17,
}


@dataclass
class LinkerSequence:
    """
    Generate and characterise flexible peptide linker sequences.

    Parameters
    ----------
    n_repeats : int
        Number of repeating units (e.g. n_repeats=3 → (G4S)3 = 15 residues).
    motif : str
        Repeat motif. Default 'GGGGS' (G4S).
    """

    n_repeats: int
    motif: str = "GGGGS"

    def __post_init__(self) -> None:
        if self.n_repeats < 1:
            raise ValueError("n_repeats must be ≥ 1")
        self.motif = self.motif.upper()

    @property
    def sequence(self) -> str:
        """Full linker amino-acid sequence (single-letter code)."""
        return self.motif * self.n_repeats

    @property
    def n_residues(self) -> int:
        """Total number of residues in the linker."""
        return len(self.sequence)

    @property
    def molecular_weight_Da(self) -> float:
        """Approximate molecular weight in Da (ignoring water loss at peptide bonds)."""
        mw = sum(_MW_PER_RESIDUE.get(aa, 110.0) for aa in self.sequence)
        # subtract water for peptide bonds
        return mw - (self.n_residues - 1) * 18.02

    @property
    def contour_length_nm(self) -> float:
        """Fully extended contour length (nm), assuming 0.38 nm per residue."""
        return self.n_residues * 0.38

    @classmethod
    def from_n_residues(
        cls,
        n_residues: int,
        motif: str = "GGGGS",
    ) -> "LinkerSequence":
        """
        Create the shortest (G4S)n linker with at least `n_residues` residues.

        Parameters
        ----------
        n_residues : int
            Minimum number of residues required.
        motif : str
            Repeat motif.

        Returns
        -------
        LinkerSequence
        """
        repeat_len = len(motif)
        n_repeats = int(np.ceil(n_residues / repeat_len))
        return cls(n_repeats=max(1, n_repeats), motif=motif)

    def as_fasta(
        self,
        header: str = "linker",
        line_width: int = 60,
    ) -> str:
        """
        Return the linker sequence in FASTA format.

        Parameters
        ----------
        header : str
            FASTA header (without the '>' prefix).
        line_width : int
            Characters per sequence line.
        """
        seq = self.sequence
        lines = [f">{header}"]
        for i in range(0, len(seq), line_width):
            lines.append(seq[i : i + line_width])
        return "\n".join(lines)

    def composition(self) -> dict[str, int]:
        """Amino-acid composition as {residue: count}."""
        comp: dict[str, int] = {}
        for aa in self.sequence:
            comp[aa] = comp.get(aa, 0) + 1
        return dict(sorted(comp.items()))

    def __repr__(self) -> str:
        return (
            f"LinkerSequence(n_repeats={self.n_repeats}, motif='{self.motif}', "
            f"n_residues={self.n_residues}, "
            f"contour_length={self.contour_length_nm:.1f} nm)"
        )


import numpy as np  # noqa: E402 — placed here to avoid circular at module level


def recommended_linkers(
    distance_nm: float,
    nanobody_kd_M: float,
    target_occupancy: float = 0.5,
    motif: str = "GGGGS",
) -> list[dict]:
    """
    Return a ranked list of recommended (G4S)n linker variants for a given
    antigen–target geometry and desired conditional occupancy.

    Parameters
    ----------
    distance_nm : float
        Antigen–target centre-to-centre distance on the membrane (nm).
    nanobody_kd_M : float
        Intrinsic Kd of the nanobody for its target (M).
    target_occupancy : float
        Desired fractional occupancy when anchored (default 0.5).
    motif : str
        Linker repeat motif.

    Returns
    -------
    list of dict, sorted by occupancy (descending), each containing:
        'n_repeats', 'n_residues', 'sequence', 'contour_length_nm',
        'c_eff_uM', 'occupancy', 'molecular_weight_Da'
    """
    from .polymer import LinkerModel

    results = []
    for n_repeats in range(1, 31):
        ls = LinkerSequence(n_repeats=n_repeats, motif=motif)
        lm = LinkerModel(n_residues=ls.n_residues)
        c_eff = lm.effective_concentration_M(distance_nm)
        occ = lm.occupancy(nanobody_kd_M, distance_nm)
        results.append(
            {
                "n_repeats": n_repeats,
                "n_residues": ls.n_residues,
                "sequence": ls.sequence,
                "contour_length_nm": ls.contour_length_nm,
                "c_eff_uM": c_eff * 1e6,
                "occupancy": occ,
                "molecular_weight_Da": ls.molecular_weight_Da,
            }
        )

    # Filter to candidates that meet the occupancy threshold
    meeting = [r for r in results if r["occupancy"] >= target_occupancy]
    if meeting:
        # Prefer the shortest linker (fewest residues) that meets the threshold
        meeting.sort(key=lambda x: x["n_residues"])
        return meeting
    # If none meets the threshold, return the top-5 by occupancy
    results.sort(key=lambda x: -x["occupancy"])
    return results[:5]
