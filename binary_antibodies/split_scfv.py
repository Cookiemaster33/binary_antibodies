"""
split_scfv.py
-------------
Thermodynamic analysis of the split-scFv competitive displacement switch.

KEY DESIGN ARCHITECTURE (two-chain)
====================================
Chain A : VH1  (binds antigen on membrane, expressed separately)
Chain B : [Minibinder]–(spacer)–VL1–(G4S)n–Nanobody

In the free state, the intramolecular minibinder on chain B occupies VL1's
VH1-pairing interface.  When VH1 anchors to the antigen on the cell surface,
the high *membrane surface concentration* of VH1 outcompetes the
minibinder and drives VH1–VL1 pairing.

Why NOT a single chain (VH1–linker–VL1)?
=========================================
If VH1 and VL1 are on the same chain, VH1's intramolecular effective
concentration (tens of mM for a 60-res linker) is enormous regardless of
whether VH1 is bound to antigen.  That eliminates the conditional switch —
VH1 and VL1 would pair spontaneously.  The two-chain architecture is
required.

Membrane surface concentration
================================
When VH1 is anchored at surface density Γ (molecules/µm²), the effective
3-D concentration experienced by a solution-phase VL1 molecule near the
membrane is:

    C_eff_surf = Γ / (Nₐ × λ)

where λ ≈ 5–10 nm is the characteristic encounter distance (roughly the
radius of the interacting protein).  For Γ = 1000 molecules/µm² and λ = 5 nm:

    C_eff_surf ≈ 330 µM

This can be much larger than a typical weakened Kd(VH–VL) of 1–50 µM.

Minibinder intramolecular effective concentration
==================================================
The minibinder is covalently tethered to chain B via a flexible spacer of
N_spacer residues.  Its effective local concentration near VL1 is:

    C_eff_MB = ρ(r=0) = (3 / (2π N_spacer l²))^(3/2) × (1/Nₐ × 10^24)   [µM]

where l = 0.38 nm.  For N_spacer = 30 residues: C_eff_MB ≈ 60 mM (!).
This is fixed and independent of antigen binding — its role is to lock VL1
in the absence of membrane VH1.  For VH1 to win, C_eff_surf must compete
with C_eff_MB × Kd_VH-VL / Kd_MB.

Design feasibility window
==========================
The switch works when both conditions are satisfied simultaneously:

  OFF (free):   C_eff_MB / Kd_MB >> 1       (minibinder locks VL1)
  ON (anchored): C_eff_surf / Kd_VH-VL >> C_eff_MB / Kd_MB  (VH1 wins)

Because C_eff_MB is fixed by the spacer length, the key lever is EITHER:
  (a) increase C_eff_surf (higher antigen expression or shorter λ), OR
  (b) use allosteric VH1: antigen binding INCREASES VH1 affinity for VL1
      (free Kd >> µM; antigen-bound Kd << µM).
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Literal

from .polymer import LinkerModel, NM3_TO_MOLAR, CA_BOND_LENGTH_NM, AVOGADRO


def intramolecular_ceff_M(n_spacer_residues: int, model: str = "fjc") -> float:
    """
    Effective local concentration (M) of a domain at the end of a
    `n_spacer_residues`-residue flexible linker, evaluated at r = 0 (contact).

    This is the relevant quantity for intramolecular competitive binding.
    """
    lm = LinkerModel(n_residues=n_spacer_residues, model=model)
    return lm.effective_concentration_M(0.0)


def membrane_surface_ceff_M(
    surface_density_per_um2: float,
    encounter_distance_nm: float = 5.0,
) -> float:
    """
    Effective 3-D concentration (M) of a membrane-anchored protein experienced
    by a solution-phase binding partner near the membrane surface.

    Parameters
    ----------
    surface_density_per_um2 : float
        Number of protein copies per µm² of cell surface.
    encounter_distance_nm : float
        Characteristic encounter distance (nm); typically 3–10 nm for proteins.

    Returns
    -------
    float
        Effective concentration in M.
    """
    # Γ [molecules/m²]  = surface_density [molecules/µm²] × 10^12 [µm²/m²]
    gamma_m2 = surface_density_per_um2 * 1e12
    # C_eff [mol/m³] = Γ / (Nₐ × λ)
    c_m3 = gamma_m2 / (AVOGADRO * encounter_distance_nm * 1e-9)
    return float(c_m3 * 1e-3)  # mol/m³ → mol/L


@dataclass
class SplitScFvSwitch:
    """
    Thermodynamic model of the split-scFv competitive displacement switch.

    Two-chain architecture:
        Chain A : VH1  (binds antigen; anchors to membrane)
        Chain B : [Minibinder]–(N_spacer residues)–VL1–(G4S)n–Nanobody

    Competitive binding equilibrium for VL1:

        VH1  + VL1  ⇌  VH1:VL1     Kd_VH-VL   (driven by membrane surface conc.)
        MB   + VL1  ⇌  MB:VL1      Kd_MB       (intramolecular, always present)

    Fraction of VL1 paired with VH1 (low-[VL1] approximation):

        f_VH1 = (C_surf / Kd_VH-VL)
              / (1 + C_surf / Kd_VH-VL + C_eff_MB / Kd_MB)

    Parameters
    ----------
    kd_vh_vl_M : float
        Kd of the engineered VH1–VL1 interaction (M).
        Target: 1–50 µM (weakened from natural sub-nM).
    kd_minibinder_vl_M : float
        Kd of the designed minibinder for VL1 (M).
        Must be < Kd_VH-VL to reliably lock VL1 in the free state.
    kd_nanobody_M : float
        Intrinsic Kd of the nanobody for the Target (M).
    minibinder_spacer_residues : int
        Number of residues in the flexible spacer between the minibinder
        and VL1 on chain B.  Determines C_eff_MB.
    n_minibinders : int
        Number of minibinder copies on chain B (typically 1 or 2).
    """

    kd_vh_vl_M: float
    kd_minibinder_vl_M: float
    kd_nanobody_M: float
    minibinder_spacer_residues: int = 30
    n_minibinders: int = 2

    @property
    def c_eff_mb_M(self) -> float:
        """Total intramolecular minibinder effective concentration near VL1 (M)."""
        single = intramolecular_ceff_M(self.minibinder_spacer_residues)
        return single * self.n_minibinders

    # ------------------------------------------------------------------
    # Core equilibrium
    # ------------------------------------------------------------------

    def vl1_paired_fraction(
        self,
        c_surf_vh1_M: float,
    ) -> float:
        """
        Fraction of VL1 in the VH1:VL1 state.

        Parameters
        ----------
        c_surf_vh1_M : float
            Effective surface concentration of membrane-anchored VH1 (M).
            Use `membrane_surface_ceff_M()` to compute this from antigen density.
        """
        num = c_surf_vh1_M / self.kd_vh_vl_M
        denom = 1.0 + num + self.c_eff_mb_M / self.kd_minibinder_vl_M
        return float(num / denom)

    def vl1_locked_fraction(self, c_surf_vh1_M: float = 0.0) -> float:
        """Fraction of VL1 locked by minibinders."""
        return 1.0 - self.vl1_paired_fraction(c_surf_vh1_M)

    # ------------------------------------------------------------------
    # Full construct analysis
    # ------------------------------------------------------------------

    def nanobody_occupancy_on(
        self,
        distance_nm: float,
        linker_n_residues: int,
        surface_density_per_um2: float,
        encounter_distance_nm: float = 5.0,
        model: str = "fjc",
    ) -> float:
        """
        Fractional nanobody occupancy on the target in the ON state (VH1 anchored).

        This is: P(VH1 paired with VL1) × P(nanobody bound given VH1 anchored).

        Parameters
        ----------
        distance_nm : float
            Antigen–target distance on the membrane (nm).
        linker_n_residues : int
            (G4S)n linker length between VL1 and the nanobody on chain B.
        surface_density_per_um2 : float
            Antigen (and thus VH1) surface density (molecules/µm²).
        encounter_distance_nm : float
            Protein encounter distance used in surface concentration calculation.
        """
        c_surf = membrane_surface_ceff_M(surface_density_per_um2, encounter_distance_nm)
        f_paired = self.vl1_paired_fraction(c_surf)

        lm = LinkerModel(n_residues=linker_n_residues, model=model)
        nb_occ = lm.occupancy(self.kd_nanobody_M, distance_nm)

        return f_paired * nb_occ

    def nanobody_occupancy_off(self, bulk_construct_M: float = 1e-9) -> float:
        """
        Fractional nanobody occupancy in the OFF state (free chain B in solution,
        no membrane VH1).
        """
        f_paired_free = self.vl1_paired_fraction(c_surf_vh1_M=0.0)
        # In free solution the nanobody must find the target at bulk concentration
        nb_bulk_occ = bulk_construct_M / (bulk_construct_M + self.kd_nanobody_M)
        return f_paired_free * nb_bulk_occ

    def selectivity_ratio(
        self,
        distance_nm: float,
        linker_n_residues: int,
        surface_density_per_um2: float,
        bulk_construct_M: float = 1e-9,
        encounter_distance_nm: float = 5.0,
        model: str = "fjc",
    ) -> float:
        """ON-state occupancy / OFF-state occupancy."""
        on = self.nanobody_occupancy_on(
            distance_nm, linker_n_residues, surface_density_per_um2,
            encounter_distance_nm, model,
        )
        off = self.nanobody_occupancy_off(bulk_construct_M)
        if off == 0:
            return float("inf")
        return on / off

    # ------------------------------------------------------------------
    # Design feasibility
    # ------------------------------------------------------------------

    def feasibility(self, surface_density_per_um2: float) -> dict:
        """
        Assess design feasibility for a given antigen surface density.

        Returns a dict with computed metrics and a 'feasible' flag.
        """
        c_surf = membrane_surface_ceff_M(surface_density_per_um2)
        c_mb = self.c_eff_mb_M

        f_paired_on = self.vl1_paired_fraction(c_surf)
        f_locked_off = self.vl1_locked_fraction(0.0)

        # OFF-state quality: C_eff_MB / Kd_MB >> 1
        mb_locking_term = c_mb / self.kd_minibinder_vl_M
        # ON-state quality: C_surf / Kd_VH-VL >> C_eff_MB / Kd_MB
        on_vs_off = (c_surf / self.kd_vh_vl_M) / max(mb_locking_term, 1e-30)

        flags = []
        if f_locked_off < 0.90:
            flags.append(
                f"OFF-state locking is poor ({f_locked_off:.1%}). "
                f"Increase minibinder affinity (Kd_MB << {self.kd_minibinder_vl_M*1e6:.1f} µM) "
                f"or increase n_minibinders / reduce spacer length."
            )
        if f_paired_on < 0.20:
            flags.append(
                f"ON-state pairing is weak ({f_paired_on:.1%}). "
                f"Options: (a) higher antigen expression (current: {surface_density_per_um2:.0f}/µm²), "
                f"(b) longer MB spacer to reduce C_eff_MB, "
                f"(c) tighter VH1-VL1 affinity (current Kd = {self.kd_vh_vl_M*1e6:.1f} µM), "
                f"(d) allosteric VH1 (antigen binding increases VL1 affinity)."
            )

        return {
            "surface_density_per_um2": surface_density_per_um2,
            "c_surf_M": c_surf,
            "c_eff_mb_M": c_mb,
            "mb_locking_term": mb_locking_term,
            "on_vs_off_ratio": on_vs_off,
            "f_vl1_paired_on": f_paired_on,
            "f_vl1_locked_off": f_locked_off,
            "feasible": len(flags) == 0,
            "flags": flags,
        }

    # ------------------------------------------------------------------
    # Scanning helpers
    # ------------------------------------------------------------------

    def scan_antigen_density(
        self,
        linker_n_residues: int,
        distance_nm: float,
        densities_per_um2: np.ndarray | None = None,
        model: str = "fjc",
    ) -> dict[str, np.ndarray]:
        """Scan nanobody ON-occupancy over a range of antigen surface densities."""
        if densities_per_um2 is None:
            densities_per_um2 = np.logspace(1, 5, 80)  # 10 – 100,000 / µm²
        occ = np.array([
            self.nanobody_occupancy_on(distance_nm, linker_n_residues, rho)
            for rho in densities_per_um2
        ])
        f_paired = np.array([
            self.vl1_paired_fraction(membrane_surface_ceff_M(rho))
            for rho in densities_per_um2
        ])
        return {
            "density_per_um2": densities_per_um2,
            "c_surf_M": np.array([membrane_surface_ceff_M(r) for r in densities_per_um2]),
            "f_vl1_paired": f_paired,
            "nanobody_occupancy_on": occ,
        }

    def scan_spacer_length(
        self,
        surface_density_per_um2: float,
        spacer_lengths: np.ndarray | None = None,
    ) -> dict[str, np.ndarray]:
        """Scan f_locked (OFF) and f_paired (ON) over minibinder spacer lengths."""
        if spacer_lengths is None:
            spacer_lengths = np.arange(5, 201, 5, dtype=int)

        c_surf = membrane_surface_ceff_M(surface_density_per_um2)
        f_locked = []
        f_paired = []

        for n_sp in spacer_lengths:
            sw = SplitScFvSwitch(
                kd_vh_vl_M=self.kd_vh_vl_M,
                kd_minibinder_vl_M=self.kd_minibinder_vl_M,
                kd_nanobody_M=self.kd_nanobody_M,
                minibinder_spacer_residues=int(n_sp),
                n_minibinders=self.n_minibinders,
            )
            f_locked.append(sw.vl1_locked_fraction(0.0))
            f_paired.append(sw.vl1_paired_fraction(c_surf))

        return {
            "spacer_residues": spacer_lengths,
            "f_vl1_locked_off": np.array(f_locked),
            "f_vl1_paired_on": np.array(f_paired),
        }

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def summary(
        self,
        distance_nm: float,
        linker_n_residues: int,
        surface_density_per_um2: float = 1000.0,
        model: str = "fjc",
    ) -> str:
        c_surf = membrane_surface_ceff_M(surface_density_per_um2)
        c_mb = self.c_eff_mb_M
        feas = self.feasibility(surface_density_per_um2)
        occ_on = self.nanobody_occupancy_on(
            distance_nm, linker_n_residues, surface_density_per_um2, model=model
        )
        occ_off = self.nanobody_occupancy_off(1e-9)
        sel = occ_on / occ_off if occ_off > 0 else float("inf")
        status = "FEASIBLE" if feas["feasible"] else "REVIEW"

        lm = LinkerModel(n_residues=linker_n_residues, model=model)

        lines = [
            f"=== Split-scFv Switch Analysis  [{status}] ===",
            f"",
            f"  Architecture          : Two-chain (Chain A: VH1 | Chain B: MB–VL1–Nb)",
            f"  Linker (VL1→Nanobody) : {linker_n_residues} res "
            f"(Lc={lm.contour_length_nm:.1f} nm, r_rms={lm.rms_end_to_end_nm:.1f} nm)",
            f"  MB spacer             : {self.minibinder_spacer_residues} res "
            f"({self.n_minibinders} copy/copies)",
            f"",
            f"  Antigen–target dist.  : {distance_nm:.1f} nm",
            f"  Antigen density       : {surface_density_per_um2:.0f} /µm²",
            f"  C_surf(VH1)           : {c_surf*1e6:.1f} µM  (membrane effective conc.)",
            f"  C_eff(MB, intra)      : {c_mb*1e6:.0f} µM  (intramolecular minibinder)",
            f"",
            f"  Kd(VH1–VL1)          : {self.kd_vh_vl_M*1e6:.1f} µM",
            f"  Kd(minibinder–VL1)   : {self.kd_minibinder_vl_M*1e6:.2f} µM",
            f"  Kd(nanobody–target)  : {self.kd_nanobody_M*1e9:.1f} nM",
            f"",
            f"  OFF-state (free):     VL1 locked = {feas['f_vl1_locked_off']:.1%}",
            f"  ON-state (anchored):  VL1 paired = {feas['f_vl1_paired_on']:.1%}",
            f"",
            f"  Nanobody occ. ON  : {occ_on:.1%}",
            f"  Nanobody occ. OFF : {occ_off:.4%}  (at 1 nM bulk)",
            f"  Selectivity ratio : {sel:.0f}×",
        ]
        if feas["flags"]:
            lines.append(f"")
            lines.append(f"  ⚠ Design warnings:")
            for f in feas["flags"]:
                lines.append(f"    • {f}")
        else:
            lines.append(f"")
            lines.append(f"  ✓ All design constraints satisfied.")
        return "\n".join(lines)
