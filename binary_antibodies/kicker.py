"""
kicker.py
---------
Geometric and thermodynamic analysis of the CH1-kicker mechanism.

In the refined design, the CL (constant light chain) domain physically
displaces the minibinder from the nanobody's CDR face when VH1-VL1
pairing occurs. This is a STRUCTURAL switch, not a thermodynamic
competition — the CL "kicks" the minibinder off the nanobody.

Key quantities
--------------
1. Kicking distance: how far does the CL move upon VH1-VL1 pairing,
   and does it overlap with the minibinder's footprint on the nanobody?
2. Linker geometry: what (G4S)n length positions the nanobody such that
   the CL can kick the minibinder when VH1-VL1 pairs?
3. Minibinder requirements: what affinity (Kd) for the nanobody CDR face
   is needed for the OFF state to be robust while still allowing kicking?
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Optional

from .polymer import LinkerModel, NM3_TO_MOLAR

# ── Immunoglobulin domain geometry constants ─────────────────────────

# Approximate dimensions of an Ig domain (β-sandwich)
IG_DOMAIN_LENGTH_NM = 4.0     # long axis (nm)
IG_DOMAIN_WIDTH_NM  = 2.5     # short axis (nm)

# In a Fab, distance from VL-CL junction to the tip of CL: ~3.5 nm
VL_CL_JUNCTION_TO_CL_TIP_NM = 3.5

# Typical displacement of CL tip when VL goes from "free" to "VH-paired":
# VH-VL pairing rotates VL ~15–25° relative to VH, which translates to
# a ~2–4 nm shift of the CL domain tip.
CL_TIP_DISPLACEMENT_RANGE_NM = (2.0, 4.0)   # min, max displacement (nm)

# Nanobody CDR face dimensions (approx)
NANOBODY_CDR_FACE_WIDTH_NM = 2.5    # approximate footprint width
NANOBODY_CDR_FACE_HEIGHT_NM = 1.5   # approximate footprint height


@dataclass
class KickerGeometry:
    """
    Model the geometry of the CH1/CL kicker mechanism.

    Parameters
    ----------
    linker_n_residues : int
        Number of (G4S) linker residues between CL C-terminus and the
        nanobody N-terminus. This sets how far the nanobody sits from CL,
        and thus how well the CL kicker can reach the minibinder.
    minibinder_kd_M : float
        Kd of the minibinder for the nanobody CDR face (M).
        This is what the CL must overcome to displace the minibinder.
    nanobody_kd_target_M : float
        Intrinsic Kd of the free nanobody for its target (M).
    cl_displacement_nm : float
        Estimated CL tip displacement upon VH1-VL1 pairing (nm).
        Use midpoint of CL_TIP_DISPLACEMENT_RANGE_NM (~3 nm) as default.
    """

    linker_n_residues: int
    minibinder_kd_M: float
    nanobody_kd_target_M: float
    cl_displacement_nm: float = 3.0  # typical CL tip displacement on pairing

    # ------------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------------

    @property
    def linker_rms_nm(self) -> float:
        """RMS end-to-end distance of the (G4S)n linker (nm)."""
        return LinkerModel(self.linker_n_residues).rms_end_to_end_nm

    @property
    def linker_contour_nm(self) -> float:
        """Contour length of the linker (nm)."""
        return self.linker_n_residues * 0.38

    def cl_nanobody_distance_nm(self) -> float:
        """
        Expected distance between CL C-terminus and nanobody centroid
        when the linker is at its RMS end-to-end distance.
        """
        return self.linker_rms_nm

    def kicking_overlap_nm(self) -> float:
        """
        Estimated overlap between the displaced CL domain and the
        minibinder footprint on the nanobody CDR face.

        A positive value means the CL clashes with the minibinder
        (kicking is geometrically feasible).

        overlap ≈ CL_displacement - (linker_rms - NANOBODY_CDR_FACE_WIDTH/2)

        If the nanobody is held at distance d from CL (via the linker),
        and CL moves by Δ upon pairing, the overlap is:
            overlap = Δ - max(0, d - CDR_face_radius)
        """
        d = self.cl_nanobody_distance_nm()
        cdr_radius = NANOBODY_CDR_FACE_WIDTH_NM / 2.0
        gap_to_bridge = max(0.0, d - cdr_radius)
        return self.cl_displacement_nm - gap_to_bridge

    def is_kicking_feasible(self) -> bool:
        """
        Returns True if the CL displacement is large enough to physically
        reach and clash with the minibinder on the nanobody CDR face.
        """
        return self.kicking_overlap_nm() > 0

    # ------------------------------------------------------------------
    # Thermodynamic model
    # ------------------------------------------------------------------

    def displacement_free_energy_kcal(self) -> float:
        """
        Approximate free energy driving displacement of the minibinder,
        arising from the steric clash between CL and minibinder.

        Estimated from the overlap distance using a simple linear spring:
            ΔG_steric ≈ k_spring × overlap²  (Hookean approximation)
        with k_spring ≈ 50 kcal/mol/nm² (typical protein clash stiffness).

        Returns
        -------
        float
            ΔG in kcal/mol. Positive = driving displacement.
        """
        overlap = max(0.0, self.kicking_overlap_nm())
        k_spring = 50.0  # kcal/mol/nm²
        return 0.5 * k_spring * overlap ** 2

    def minibinder_binding_free_energy_kcal(self) -> float:
        """
        Free energy of minibinder binding to nanobody CDR face.

        ΔG_bind = RT × ln(Kd)  (negative = stable binding)
        """
        RT = 0.593  # kcal/mol at 25°C
        return RT * np.log(self.minibinder_kd_M)

    def kicking_is_thermodynamically_favourable(self) -> bool:
        """
        Returns True if the steric displacement free energy exceeds the
        minibinder binding free energy (kicking is thermodynamically driven).
        """
        dg_steric = self.displacement_free_energy_kcal()
        dg_bind = abs(self.minibinder_binding_free_energy_kcal())
        return dg_steric > dg_bind

    def effective_displacement_fraction(self) -> float:
        """
        Estimated fraction of minibinders displaced by the CL kick.

        Uses a simple two-state model:
            f_displaced = 1 / (1 + exp(-ΔΔG / RT))
        where ΔΔG = ΔG_steric - |ΔG_bind|

        Returns
        -------
        float
            Fraction of minibinders displaced (0–1).
        """
        RT = 0.593
        dg_steric = self.displacement_free_energy_kcal()
        dg_bind = abs(self.minibinder_binding_free_energy_kcal())
        ddg = dg_steric - dg_bind
        return float(1.0 / (1.0 + np.exp(-ddg / RT)))

    # ------------------------------------------------------------------
    # Optimal linker
    # ------------------------------------------------------------------

    @classmethod
    def optimal_linker_for_kicking(
        cls,
        minibinder_kd_M: float,
        nanobody_kd_target_M: float,
        min_displacement_fraction: float = 0.9,
        cl_displacement_nm: float = 3.0,
        max_residues: int = 100,
    ) -> Optional[int]:
        """
        Find the minimum linker length (residues) that achieves at least
        `min_displacement_fraction` displacement of the minibinder.

        Shorter linkers keep the nanobody close to CL → better kicking.
        But too short → nanobody can't reach the target.

        Parameters
        ----------
        minibinder_kd_M : float
            Kd of minibinder for nanobody CDRs (M).
        nanobody_kd_target_M : float
            Kd of nanobody for its target (M).
        min_displacement_fraction : float
            Minimum required kicking efficiency.
        cl_displacement_nm : float
            CL tip displacement upon VH1-VL1 pairing (nm).
        max_residues : int
            Maximum linker length to search.

        Returns
        -------
        int or None
            Minimum linker residues, or None if not achievable.
        """
        # Prefer shorter linkers (better kicking), but need minimum for
        # nanobody to physically separate from CL
        best = None
        for n in range(3, max_residues + 1):
            kg = cls(
                linker_n_residues=n,
                minibinder_kd_M=minibinder_kd_M,
                nanobody_kd_target_M=nanobody_kd_target_M,
                cl_displacement_nm=cl_displacement_nm,
            )
            if kg.effective_displacement_fraction() >= min_displacement_fraction:
                best = n
                break  # first (shortest) one that works
        return best

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def summary(self) -> str:
        feasible = self.is_kicking_feasible()
        thermo = self.kicking_is_thermodynamically_favourable()
        frac = self.effective_displacement_fraction()
        overlap = self.kicking_overlap_nm()

        status = "FEASIBLE ✓" if (feasible and thermo) else "REVIEW NEEDED ⚠"
        lines = [
            f"=== CH1-Kicker Geometry Analysis  [{status}] ===",
            f"",
            f"  Linker                  : {self.linker_n_residues} res "
            f"(Lc={self.linker_contour_nm:.1f} nm, r_rms={self.linker_rms_nm:.1f} nm)",
            f"  CL tip displacement     : {self.cl_displacement_nm:.1f} nm "
            f"(upon VH1-VL1 pairing)",
            f"  CL-nanobody distance    : {self.cl_nanobody_distance_nm():.1f} nm",
            f"  Kicking overlap         : {overlap:.2f} nm  "
            f"({'clash — kicking possible' if overlap > 0 else 'no overlap — kicking impossible'})",
            f"",
            f"  Minibinder Kd (CDR face): {self.minibinder_kd_M*1e9:.1f} nM",
            f"  ΔG_bind (minibinder)    : {self.minibinder_binding_free_energy_kcal():.1f} kcal/mol",
            f"  ΔG_steric (clash)       : {self.displacement_free_energy_kcal():.1f} kcal/mol",
            f"",
            f"  Displacement fraction   : {frac:.1%}",
            f"  Kicking thermodynamically favourable: {thermo}",
        ]
        return "\n".join(lines)


@dataclass
class KickerConstruct:
    """
    Full construct analysis for the CH1-kicker design.

    Chain B architecture (N→C):
        [VL1] – [CL] – (G4S)_linker – [Minibinder] – (G4S)3 – [Nanobody]

    Parameters
    ----------
    kd_antigen_M : float
        Kd of VH1 for the antigen (M).
    kd_vh_vl_M : float
        Kd of VH1 for VL1 in solution (M). Should be µM — weak enough
        that spontaneous pairing is rare, but driven on membrane.
    kd_minibinder_cdrs_M : float
        Kd of minibinder for nanobody CDR face (M). Target: nM–low µM.
    kd_nanobody_M : float
        Intrinsic Kd of nanobody for its target (M).
    cl_to_nb_linker_n : int
        (G4S)n residues between CL C-terminus and nanobody N-terminus.
    cl_displacement_nm : float
        CL tip displacement upon VH1-VL1 pairing (nm). Default 3.0.
    """

    kd_antigen_M: float
    kd_vh_vl_M: float
    kd_minibinder_cdrs_M: float
    kd_nanobody_M: float
    cl_to_nb_linker_n: int
    cl_displacement_nm: float = 3.0
    name: str = "CH1-kicker conditional nanobody"

    @property
    def kicker(self) -> KickerGeometry:
        return KickerGeometry(
            linker_n_residues=self.cl_to_nb_linker_n,
            minibinder_kd_M=self.kd_minibinder_cdrs_M,
            nanobody_kd_target_M=self.kd_nanobody_M,
            cl_displacement_nm=self.cl_displacement_nm,
        )

    def off_state_nanobody_availability(self) -> float:
        """
        Fraction of nanobody with free CDRs in the OFF state (no VH1-VL1 pairing).
        Should be close to 0 for a good OFF state.

        Approximated as: 1 - Kd_MB / (Kd_MB + C_eff_MB)
        where C_eff_MB is the intramolecular minibinder effective concentration
        (set by the (G4S)3 linker between minibinder and nanobody).
        """
        from .polymer import LinkerModel
        # Short (G4S)3 linker between minibinder and nanobody
        mb_spacer = LinkerModel(n_residues=15)  # ~(G4S)3
        c_eff_mb = mb_spacer.effective_concentration_M(0.0)
        # Fraction unblocked = 1/(1 + C_eff/Kd_MB)
        # Fraction blocked = C_eff/(C_eff + Kd_MB) = 1 - fraction unblocked
        frac_blocked = c_eff_mb / (c_eff_mb + self.kd_minibinder_cdrs_M)
        return 1.0 - frac_blocked  # fraction NOT blocked (available for target)

    def on_state_nanobody_availability(self) -> float:
        """
        Fraction of nanobody with free CDRs in the ON state (VH1-VL1 paired,
        CL kicks minibinder off CDRs).
        """
        return self.kicker.effective_displacement_fraction()

    def summary(self) -> str:
        off = self.off_state_nanobody_availability()
        on = self.on_state_nanobody_availability()
        sel = on / off if off > 0 else float("inf")

        lines = [
            f"=== {self.name} ===",
            f"",
            f"  VH1 Kd (antigen)       : {self.kd_antigen_M*1e9:.2f} nM",
            f"  VH1-VL1 Kd             : {self.kd_vh_vl_M*1e6:.1f} µM",
            f"  Minibinder Kd (CDRs)   : {self.kd_minibinder_cdrs_M*1e9:.1f} nM",
            f"  Nanobody Kd (target)   : {self.kd_nanobody_M*1e9:.1f} nM",
            f"  CL-nanobody linker     : {self.cl_to_nb_linker_n} residues",
            f"",
            self.kicker.summary(),
            f"",
            f"  Nanobody available OFF : {off:.2%}  (CDRs blocked by minibinder)",
            f"  Nanobody available ON  : {on:.2%}  (CDRs freed by CL kick)",
            f"  Selectivity ratio      : {sel:.0f}×",
        ]
        return "\n".join(lines)
