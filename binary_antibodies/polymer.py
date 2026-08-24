"""
polymer.py
----------
Physics of flexible peptide linkers using the freely-jointed chain (FJC)
and worm-like chain (WLC) models.

Key quantity: the effective local concentration of the nanobody end-point
relative to the antibody attachment point, given the linker length and the
antigen–target distance.
"""

from __future__ import annotations

import numpy as np
from scipy import integrate, special
from dataclasses import dataclass, field
from typing import Literal

# Physical constants
AVOGADRO = 6.022e23          # mol⁻¹
NM3_TO_MOLAR = 1e24 / AVOGADRO  # converts nm⁻³ → M  (= 1/(Nₐ × 1e-24))

# Peptide linker geometry
CA_BOND_LENGTH_NM = 0.38     # virtual Cα–Cα bond length (nm)
PERSISTENCE_LENGTH_NM = 0.4  # persistence length for unstructured peptide (nm)


@dataclass
class LinkerModel:
    """
    Model a flexible peptide linker using the freely-jointed chain (FJC)
    or worm-like chain (WLC) approximation.

    Parameters
    ----------
    n_residues : int
        Number of amino-acid residues in the linker.
    model : 'fjc' or 'wlc'
        Which polymer model to use.
        - 'fjc': freely-jointed chain (simpler, analytical)
        - 'wlc': worm-like chain (more realistic for (G4S)n linkers)
    bond_length_nm : float
        Effective monomer (Cα–Cα) bond length in nm.
    persistence_length_nm : float
        Persistence length (only used for WLC model).
    """

    n_residues: int
    model: Literal["fjc", "wlc"] = "fjc"
    bond_length_nm: float = CA_BOND_LENGTH_NM
    persistence_length_nm: float = PERSISTENCE_LENGTH_NM

    # Derived quantities cached after first access
    _contour_length_nm: float = field(init=False, repr=False)
    _rms_end_to_end_nm: float = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if self.n_residues < 1:
            raise ValueError("n_residues must be ≥ 1")
        self._contour_length_nm = self.n_residues * self.bond_length_nm
        # FJC root-mean-square end-to-end distance
        self._rms_end_to_end_nm = self.bond_length_nm * np.sqrt(self.n_residues)

    # ------------------------------------------------------------------
    # Core quantity: P(r) — probability density of end-to-end distance r
    # ------------------------------------------------------------------

    def end_to_end_pdf(self, r_nm: np.ndarray | float) -> np.ndarray:
        """
        Radial probability density P(r) [nm⁻¹] of the end-to-end distance.

        For the FJC model this is the 3-D Gaussian:
            P(r) = 4π r² × (3 / (2π <r²>))^(3/2) × exp(-3r²/(2<r²>))

        For the WLC model we use the interpolation formula by Marko & Siggia
        (Macromolecules 1995) for the force–extension, and invert numerically,
        but here we use the simpler 3-D Gaussian with WLC-corrected variance:
            <r²>_WLC = 2 Lp Lc (1 − Lp/Lc (1 − exp(−Lc/Lp)))

        Parameters
        ----------
        r_nm : array-like
            End-to-end distance(s) in nm.

        Returns
        -------
        np.ndarray
            Probability density in nm⁻¹.
        """
        r = np.asarray(r_nm, dtype=float)
        r2_mean = self._mean_squared_end_to_end()
        # 3-D Gaussian distribution (valid when N >> 1 or Lc >> Lp)
        prefactor = 4.0 * np.pi * r**2 * (3.0 / (2.0 * np.pi * r2_mean)) ** 1.5
        exponent = np.exp(-3.0 * r**2 / (2.0 * r2_mean))
        return prefactor * exponent

    def _mean_squared_end_to_end(self) -> float:
        """<r²> in nm² for the chosen model."""
        if self.model == "fjc":
            return self.n_residues * self.bond_length_nm**2
        elif self.model == "wlc":
            lp = self.persistence_length_nm
            lc = self._contour_length_nm
            # Exact WLC formula for <r²>
            return 2.0 * lp * lc * (1.0 - (lp / lc) * (1.0 - np.exp(-lc / lp)))
        else:
            raise ValueError(f"Unknown model '{self.model}'. Choose 'fjc' or 'wlc'.")

    @property
    def rms_end_to_end_nm(self) -> float:
        """Root-mean-square end-to-end distance √<r²> in nm."""
        return float(np.sqrt(self._mean_squared_end_to_end()))

    @property
    def contour_length_nm(self) -> float:
        """Fully extended contour length of the linker in nm."""
        return self._contour_length_nm

    # ------------------------------------------------------------------
    # Effective local concentration
    # ------------------------------------------------------------------

    def effective_concentration_M(self, distance_nm: float) -> float:
        """
        Effective local concentration (M) of the nanobody end-point at a
        given distance from the antibody attachment point.

        This is the probability density P(r) evaluated at r = distance_nm,
        converted to molar units (mol/L).

        C_eff(r) = P(r) / (4π r² dr) × (1 / Nₐ) integrated over a small
        shell — but practically this equals:

            C_eff = P(r=d) × (3/(4π r²)) ÷ Nₐ

        More precisely we integrate P(r) over a 1-nm shell around d:

            C_eff = ∫_{d-0.5}^{d+0.5} P(r) dr × NM3_TO_MOLAR × 1/nm³

        Parameters
        ----------
        distance_nm : float
            Antigen–target distance on the membrane (nm).

        Returns
        -------
        float
            Effective local concentration in molar (M).
        """
        if distance_nm < 0:
            raise ValueError("distance_nm must be non-negative")
        if distance_nm > self.contour_length_nm:
            return 0.0  # target is unreachable

        r2_mean = self._mean_squared_end_to_end()

        # Volume density at distance d (nm⁻³):
        # ρ(d) = P(r=d) / (4π d²)  [= number per nm³ per unit r at d]
        # But for d=0 the formula diverges; use the Gaussian volumetric form:
        #   ρ(d) = (3/(2π<r²>))^(3/2) × exp(-3d²/(2<r²>))
        rho_nm3 = (3.0 / (2.0 * np.pi * r2_mean)) ** 1.5 * np.exp(
            -3.0 * distance_nm**2 / (2.0 * r2_mean)
        )
        # Convert nm⁻³ → M
        return float(rho_nm3 * NM3_TO_MOLAR)

    def conditional_kd_M(self, nanobody_kd_M: float, distance_nm: float) -> float:
        """
        Apparent conditional Kd of the nanobody for its target, given that
        the antibody is anchored to the antigen at distance `distance_nm`.

        The effective local concentration raises the apparent on-rate, so:

            Kd_apparent = Kd_intrinsic / C_eff   (if C_eff >> Kd)

        More precisely, the conditional occupancy at equilibrium is:

            θ = C_eff / (C_eff + Kd_intrinsic)

        and the effective Kd at bulk concentration → 0 is Kd_intrinsic / C_eff.

        Parameters
        ----------
        nanobody_kd_M : float
            Intrinsic Kd of the nanobody for its target (M).
        distance_nm : float
            Antigen–target distance on the membrane (nm).

        Returns
        -------
        float
            Effective conditional Kd in M (lower = stronger conditional binding).
        """
        c_eff = self.effective_concentration_M(distance_nm)
        if c_eff == 0.0:
            return float("inf")
        return nanobody_kd_M / c_eff

    def occupancy(self, nanobody_kd_M: float, distance_nm: float) -> float:
        """
        Fractional occupancy of nanobody on target when the antibody is bound
        to the antigen.

            θ = C_eff / (C_eff + Kd)

        Parameters
        ----------
        nanobody_kd_M : float
            Intrinsic Kd of the nanobody (M).
        distance_nm : float
            Antigen–target distance on the membrane (nm).

        Returns
        -------
        float
            Fractional occupancy ∈ [0, 1].
        """
        c_eff = self.effective_concentration_M(distance_nm)
        return c_eff / (c_eff + nanobody_kd_M)

    def max_reachable_distance_nm(self, min_occupancy: float = 0.5) -> float:
        """
        Maximum antigen–target distance (nm) at which the nanobody can still
        achieve `min_occupancy` fractional occupancy, for a range of intrinsic
        Kd values.

        Returns the contour length as a hard upper bound.
        """
        return self.contour_length_nm

    # ------------------------------------------------------------------
    # Scanning helpers
    # ------------------------------------------------------------------

    def scan_distance(
        self,
        nanobody_kd_M: float,
        distances_nm: np.ndarray | None = None,
    ) -> dict[str, np.ndarray]:
        """
        Compute effective concentration, occupancy, and apparent Kd over a
        range of antigen–target distances.

        Parameters
        ----------
        nanobody_kd_M : float
            Intrinsic Kd of the nanobody (M).
        distances_nm : array-like, optional
            Distances to scan (nm). Defaults to 0–contour_length in 100 steps.

        Returns
        -------
        dict with keys:
            'distance_nm', 'c_eff_M', 'occupancy', 'kd_apparent_M'
        """
        if distances_nm is None:
            distances_nm = np.linspace(0.01, self.contour_length_nm, 200)
        distances_nm = np.asarray(distances_nm)

        c_eff = np.array([self.effective_concentration_M(d) for d in distances_nm])
        occ = c_eff / (c_eff + nanobody_kd_M)
        kd_app = np.where(c_eff > 0, nanobody_kd_M / c_eff, np.inf)

        return {
            "distance_nm": distances_nm,
            "c_eff_M": c_eff,
            "occupancy": occ,
            "kd_apparent_M": kd_app,
        }

    def __repr__(self) -> str:
        return (
            f"LinkerModel(n_residues={self.n_residues}, model='{self.model}', "
            f"contour_length={self.contour_length_nm:.1f} nm, "
            f"rms_end_to_end={self.rms_end_to_end_nm:.1f} nm)"
        )
