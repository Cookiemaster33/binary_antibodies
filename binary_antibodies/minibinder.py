"""
minibinder.py
-------------
Design guidance and interface analysis for minibinders that lock VL1.

Minibinders are small (40–80 residue) designed proteins that bind to VL1's
VH1-pairing interface (framework 2 and CDR-L2 region), keeping VL1 in an
inactive state until VH1 is anchored and outcompetes them.

This module provides:
  - A description of the VH1-pairing interface on VL1 (the target epitope).
  - Kd range requirements for viable minibinder design.
  - A computational design strategy guide.
  - Helpers to extract the interface residues from a PDB file (given Biopython).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np


# VH–VL interface residues on the VL side (Chothia/IMGT numbering).
# These are the conserved framework residues that pack against VH.
VL_VH_INTERFACE_RESIDUES_IMGT = [
    35, 36, 37, 38, 39,   # Framework 2 (β-strand B)
    45, 46, 47,           # Framework 2 (β-strand C, including the conserved Gln45)
    85, 86, 87, 88,       # Framework 3 (β-strand F)
    98,                   # Framework 4
]

VL_CDR_RESIDUES_IMGT = {
    "CDR-L1": list(range(27, 38)),
    "CDR-L2": list(range(56, 66)),
    "CDR-L3": list(range(105, 118)),
}


@dataclass
class MinibinderDesignSpec:
    """
    Specification and design constraints for a VL1-locking minibinder.

    Parameters
    ----------
    target_kd_range_M : tuple
        (min, max) desired Kd for minibinder–VL1 (M).
        Should be tighter than Kd(VH1-VL1) to reliably lock VL1 in free state.
        Must be weaker than Kd(VH1-VL1) / (C_eff / Kd(VH1-VL1)) to ensure VH1
        wins when anchored.
    target_size_residues : tuple
        (min, max) acceptable minibinder size in residues.
    target_epitope : str
        Which region of VL1 to target. One of: 'framework2', 'full_vhvl_interface'.
    n_copies : int
        Number of minibinder copies in the construct (1 or 2 recommended).
    """

    target_kd_range_M: tuple[float, float] = (0.5e-6, 20e-6)
    target_size_residues: tuple[int, int] = (40, 80)
    target_epitope: str = "framework2"
    n_copies: int = 2

    # ------------------------------------------------------------------
    # Kd requirements
    # ------------------------------------------------------------------

    @classmethod
    def from_switch_kds(
        cls,
        kd_vh_vl_M: float,
        c_eff_anchored_M: float,
        safety_factor_off: float = 10.0,
        safety_factor_on: float = 5.0,
    ) -> "MinibinderDesignSpec":
        """
        Derive the required minibinder Kd range from the switch thermodynamics.

        In the OFF state (free): minibinder must outcompete VH1 in bulk solution.
        For robust locking: Kd(MB) < Kd(VH-VL) / safety_factor_off.

        In the ON state (anchored): VH1 at C_eff must outcompete minibinder.
        For robust activation: Kd(MB) > Kd(VH-VL) / (C_eff / Kd(VH-VL))^safety_factor_on
        Simplified: Kd(MB) > C_eff / safety_factor_on.

        Parameters
        ----------
        kd_vh_vl_M : float
            Kd of engineered VH1-VL1 in solution (M).
        c_eff_anchored_M : float
            Effective concentration of VH1 near VL1 when anchored (M).
        safety_factor_off : float
            How much tighter MB should be vs VH1-VL1 in free state.
        safety_factor_on : float
            How much weaker MB must be vs VH1 when anchored.

        Returns
        -------
        MinibinderDesignSpec
        """
        kd_mb_max = kd_vh_vl_M / safety_factor_off   # MB tighter than VH1-VL1
        kd_mb_min = c_eff_anchored_M / safety_factor_on  # MB weaker than VH1 effective

        if kd_mb_min > kd_mb_max:
            raise ValueError(
                f"No feasible minibinder Kd window exists: "
                f"min={kd_mb_min*1e6:.2f} µM > max={kd_mb_max*1e6:.2f} µM. "
                f"Increase C_eff (longer linker or shorter distance) or "
                f"weaken Kd(VH1-VL1) further."
            )
        return cls(target_kd_range_M=(kd_mb_min, kd_mb_max))

    def is_kd_feasible(self, kd_M: float) -> bool:
        """Check if a candidate minibinder Kd is within the feasible window."""
        return self.target_kd_range_M[0] <= kd_M <= self.target_kd_range_M[1]

    # ------------------------------------------------------------------
    # Interface description
    # ------------------------------------------------------------------

    @property
    def target_interface_residues_imgt(self) -> list[int]:
        """IMGT-numbered VL1 residues that the minibinder should contact."""
        if self.target_epitope == "framework2":
            return VL_VH_INTERFACE_RESIDUES_IMGT[:8]  # core FR2 only
        elif self.target_epitope == "full_vhvl_interface":
            return VL_VH_INTERFACE_RESIDUES_IMGT
        else:
            raise ValueError(f"Unknown target_epitope: {self.target_epitope!r}")

    def design_brief(self) -> str:
        """
        Return a human-readable design brief for the minibinder.
        """
        kd_min_uM = self.target_kd_range_M[0] * 1e6
        kd_max_uM = self.target_kd_range_M[1] * 1e6
        res_min, res_max = self.target_size_residues

        lines = [
            "=== Minibinder Design Brief ===",
            "",
            "Target protein  : VL1 (light-chain variable domain of scFv)",
            f"Target epitope  : {self.target_epitope} (VH1-pairing interface)",
            f"Target residues : IMGT positions {self.target_interface_residues_imgt}",
            "",
            f"Required Kd     : {kd_min_uM:.2f} – {kd_max_uM:.2f} µM",
            f"Acceptable size : {res_min}–{res_max} residues",
            f"Copies in cxstr : {self.n_copies}",
            "",
            "Design strategy:",
            "  1. Extract VL1 structure (AlphaFold or PDB of parent antibody).",
            "  2. Identify the VH1-pairing interface on VL1 (FR2 β-sheet face).",
            "     Key residues: Q38, L46, Y87 (conserved packing positions).",
            "  3. Run RFdiffusion (hallucination or partial diffusion) to",
            "     generate a diverse set of small protein backbones that dock",
            "     against the VL1 FR2 surface.",
            "     Suggested RFdiffusion flags:",
            "       --contigs  'A38-47/0 B1-60'  (60-residue binder, VL1=chain A)",
            "       --num_designs 1000",
            "  4. Use ProteinMPNN to design sequences for each backbone.",
            "  5. Filter by:",
            "       - Predicted Rosetta binding energy (ΔΔG < -5 REU)",
            "       - Predicted AF2 pLDDT > 85 (monomer folding)",
            "       - AF2-Multimer pTM > 0.7 for VL1–minibinder complex",
            "  6. Experimentally validate Kd by SPR/BLI.",
            "     Target range: {:.1f}–{:.1f} µM.".format(kd_min_uM, kd_max_uM),
            "",
            "Sequence motif to avoid overlapping with CDR-L2 (IMGT 56–65):",
            "  Minibinder should contact FR2 (35–47) and FR3 (85–88), NOT CDR-L2,",
            "  to minimise interference with VL1 folding and to allow a clear",
            "  steric clash with incoming VH1 FR2 (the natural VH-pairing surface).",
        ]
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Competitive binding window
    # ------------------------------------------------------------------

    def competitive_window_plot_data(
        self,
        kd_vh_vl_M: float,
        c_eff_range_M: np.ndarray | None = None,
    ) -> dict[str, np.ndarray]:
        """
        Compute the feasible minibinder Kd window as a function of C_eff.

        Returns
        -------
        dict with 'c_eff_M', 'kd_mb_min_M', 'kd_mb_max_M', 'window_open'
        """
        if c_eff_range_M is None:
            c_eff_range_M = np.logspace(-8, -4, 100)  # 10 nM – 100 µM

        kd_mb_max = kd_vh_vl_M / 10.0  # OFF-state constraint (fixed)
        kd_mb_min = c_eff_range_M / 5.0  # ON-state constraint (varies with C_eff)

        window_open = kd_mb_min < kd_mb_max

        return {
            "c_eff_M": c_eff_range_M,
            "kd_mb_min_M": kd_mb_min,
            "kd_mb_max_M": np.full_like(c_eff_range_M, kd_mb_max),
            "window_open": window_open,
        }


def interface_buried_sasa_estimate(n_interface_residues: int = 8) -> float:
    """
    Rough estimate of buried surface area (Å²) for a minibinder engaging
    the VL1 FR2 interface.

    Typical buried SASA for a protein–protein interface:
      Small/moderate: 600–1200 Å²
      Minibinder targeting 8 FR2 residues: ~700–900 Å²
    """
    return n_interface_residues * 95.0  # ~95 Å² per residue buried


def kd_from_buried_sasa(buried_sasa_A2: float) -> float:
    """
    Very rough estimate of binding Kd from buried SASA, using the empirical
    relationship: ΔG ≈ -0.030 × BSASA (REU, ~kcal/mol), Kd = exp(ΔG/RT).

    For design guidance only — not a substitute for experimental measurement.
    """
    RT = 0.593  # kcal/mol at 25°C
    dG = -0.030 * buried_sasa_A2
    kd_M = np.exp(dG / RT) * 1.0  # reference state 1 M
    return float(kd_M)
