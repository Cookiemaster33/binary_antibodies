"""
design.py
---------
End-to-end design and optimisation of a conditional proximity-gated
nanobody construct.

A ConditionalConstruct bundles:
  - the antibody fragment (Kd for antigen)
  - the flexible linker (LinkerModel)
  - the nanobody (intrinsic Kd for target)

and exposes design-space queries such as:
  - Which linker length achieves a target occupancy at a given geometry?
  - What is the selectivity ratio (anchored vs free)?
  - How does conditional binding depend on antigen expression level?
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Sequence

from .polymer import LinkerModel


@dataclass
class ConditionalConstruct:
    """
    A conditional proximity-gated nanobody construct.

    Parameters
    ----------
    antibody_kd_M : float
        Intrinsic Kd of the antibody (or antibody fragment) for the antigen (M).
    nanobody_kd_M : float
        Intrinsic Kd of the nanobody for the target (M).
    linker : LinkerModel
        Flexible linker connecting antibody C-terminus to nanobody N-terminus.
    name : str, optional
        Human-readable construct name.
    """

    antibody_kd_M: float
    nanobody_kd_M: float
    linker: LinkerModel
    name: str = "Conditional nanobody construct"

    # ------------------------------------------------------------------
    # Binding state analysis
    # ------------------------------------------------------------------

    def anchoring_fraction(self, antigen_surface_density_per_um2: float) -> float:
        """
        Fraction of constructs anchored to the antigen at equilibrium,
        estimated from the surface density of the antigen.

        Uses a simple Langmuir isotherm with the effective 2-D concentration
        converted to a bulk-equivalent for the antibody Kd comparison.

        A surface density of ~100 molecules/µm² ≈ 0.17 nM bulk equivalent
        in a 1 µm diffusion layer.

        Parameters
        ----------
        antigen_surface_density_per_um2 : float
            Number of antigen copies per µm² of cell surface.

        Returns
        -------
        float
            Fraction of constructs anchored (0–1).
        """
        # Rough conversion: N/µm² × (1e-12 m²/µm²) / (1e-7 m diff. layer) / Nₐ
        # = N × 1e-5 / (6.022e23) mol/L  (very rough estimate)
        c_antigen_M = antigen_surface_density_per_um2 * 1e-5 / 6.022e23
        return c_antigen_M / (c_antigen_M + self.antibody_kd_M)

    def nanobody_occupancy_anchored(self, distance_nm: float) -> float:
        """
        Fractional occupancy of the nanobody on target, *given* the antibody
        is anchored to the antigen (i.e. the conditional term of the AND-gate).

        Parameters
        ----------
        distance_nm : float
            Centre-to-centre distance between antigen and target on the membrane.
        """
        return self.linker.occupancy(self.nanobody_kd_M, distance_nm)

    def nanobody_occupancy_free(self, target_bulk_concentration_M: float) -> float:
        """
        Fractional occupancy of the nanobody on target when the antibody is
        **not** anchored (free diffusion in bulk).

        Parameters
        ----------
        target_bulk_concentration_M : float
            Effective bulk concentration of accessible target (M).
        """
        return target_bulk_concentration_M / (
            target_bulk_concentration_M + self.nanobody_kd_M
        )

    def selectivity_ratio(
        self,
        distance_nm: float,
        target_bulk_concentration_M: float = 1e-9,
    ) -> float:
        """
        Selectivity ratio: occupancy_anchored / occupancy_free.

        A ratio >> 1 means the nanobody preferentially binds the target
        only when the antibody is anchored.

        Parameters
        ----------
        distance_nm : float
            Antigen–target distance (nm).
        target_bulk_concentration_M : float
            Bulk target concentration used for free-state denominator.
        """
        occ_on = self.nanobody_occupancy_anchored(distance_nm)
        occ_off = self.nanobody_occupancy_free(target_bulk_concentration_M)
        if occ_off == 0:
            return float("inf")
        return occ_on / occ_off

    # ------------------------------------------------------------------
    # Optimisation
    # ------------------------------------------------------------------

    def optimal_linker_length(
        self,
        distance_nm: float,
        target_occupancy: float = 0.5,
        model: str = "fjc",
        max_residues: int = 200,
    ) -> int | None:
        """
        Find the minimum linker length (in residues) that achieves at least
        `target_occupancy` when the antibody is anchored and the target is
        `distance_nm` away.

        Parameters
        ----------
        distance_nm : float
            Antigen–target distance (nm).
        target_occupancy : float
            Desired fractional occupancy of nanobody on target (0–1).
        model : str
            Polymer model ('fjc' or 'wlc').
        max_residues : int
            Upper bound for the search.

        Returns
        -------
        int or None
            Minimum number of linker residues, or None if not achievable within
            max_residues.
        """
        for n in range(1, max_residues + 1):
            lm = LinkerModel(n_residues=n, model=model)
            if lm.contour_length_nm < distance_nm:
                continue  # physically can't reach
            occ = lm.occupancy(self.nanobody_kd_M, distance_nm)
            if occ >= target_occupancy:
                return n
        return None

    def scan_linker_lengths(
        self,
        distance_nm: float,
        lengths: Sequence[int] | None = None,
        model: str = "fjc",
    ) -> dict[str, np.ndarray]:
        """
        Scan nanobody occupancy over a range of linker lengths for a fixed
        antigen–target distance.

        Returns
        -------
        dict with keys: 'n_residues', 'occupancy', 'c_eff_M', 'selectivity_ratio'
        """
        if lengths is None:
            lengths = list(range(5, 151, 5))
        lengths_arr = np.array(lengths, dtype=int)
        occ = np.zeros(len(lengths_arr))
        c_eff = np.zeros(len(lengths_arr))

        for i, n in enumerate(lengths_arr):
            lm = LinkerModel(n_residues=int(n), model=model)
            c_eff[i] = lm.effective_concentration_M(distance_nm)
            occ[i] = lm.occupancy(self.nanobody_kd_M, distance_nm)

        return {
            "n_residues": lengths_arr,
            "occupancy": occ,
            "c_eff_M": c_eff,
        }

    def parameter_space_heatmap(
        self,
        distances_nm: np.ndarray | None = None,
        linker_lengths: np.ndarray | None = None,
        model: str = "fjc",
    ) -> dict[str, np.ndarray]:
        """
        Compute nanobody occupancy over the 2-D space of
        (antigen–target distance) × (linker length).

        Returns
        -------
        dict with keys:
            'distances_nm'    : 1-D array (shape D)
            'linker_lengths'  : 1-D array (shape L)
            'occupancy'       : 2-D array (shape D × L)
            'c_eff_M'         : 2-D array (shape D × L)
        """
        if distances_nm is None:
            distances_nm = np.linspace(1, 20, 40)
        if linker_lengths is None:
            linker_lengths = np.arange(10, 151, 10, dtype=int)

        D, L = len(distances_nm), len(linker_lengths)
        occ = np.zeros((D, L))
        c_eff = np.zeros((D, L))

        for j, n in enumerate(linker_lengths):
            lm = LinkerModel(n_residues=int(n), model=model)
            for i, d in enumerate(distances_nm):
                c_eff[i, j] = lm.effective_concentration_M(d)
                occ[i, j] = lm.occupancy(self.nanobody_kd_M, d)

        return {
            "distances_nm": distances_nm,
            "linker_lengths": linker_lengths,
            "occupancy": occ,
            "c_eff_M": c_eff,
        }

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def summary(self, distance_nm: float) -> str:
        """
        Print a human-readable design summary for a given antigen–target distance.
        """
        c_eff = self.linker.effective_concentration_M(distance_nm)
        occ = self.nanobody_occupancy_anchored(distance_nm)
        lines = [
            f"=== {self.name} ===",
            f"  Antibody Kd (antigen):          {self.antibody_kd_M*1e9:.2f} nM",
            f"  Nanobody Kd (target, intrinsic):{self.nanobody_kd_M*1e9:.2f} nM",
            f"  Linker:                         {self.linker}",
            f"",
            f"  Antigen–target distance:        {distance_nm:.1f} nm",
            f"  Effective local [nanobody]:     {c_eff*1e6:.3f} µM",
            f"  Nanobody occupancy (anchored):  {occ:.3f}  ({occ*100:.1f}%)",
            f"  Nanobody occupancy (free, 1 nM):{self.nanobody_occupancy_free(1e-9):.6f}",
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"ConditionalConstruct("
            f"antibody_kd={self.antibody_kd_M*1e9:.1f} nM, "
            f"nanobody_kd={self.nanobody_kd_M*1e9:.1f} nM, "
            f"linker={self.linker.n_residues} res)"
        )
