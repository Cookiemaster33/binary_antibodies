"""
structures.py
-------------
PDB structure downloading, antibody domain extraction, and VH–VL interface
analysis for the conditional nanobody design project.

Example biological system
--------------------------
Antigen  : EGFR (EGF receptor) — used by VH1 to anchor to the cell surface
Target   : HER2 (ErbB2)        — recognised by VH1+VL1 scFv and by the nanobody
Rationale: EGFR and HER2 are co-overexpressed in ~30% of breast cancers and many
           lung/colorectal cancers; a construct that requires BOTH to be present
           is more tumour-selective than either alone.

PDB entries used
-----------------
Trastuzumab Fab  PDB 1N8Z  chains A (VH) / B (VL) — anti-HER2 scFv source
Cetuximab Fab    PDB 1YY9  chains H (VH) / L (VL) — anti-EGFR VH template
Anti-HER2 Nb     PDB 5MY6  chain A                 — 2Rs15d nanobody vs HER2
"""

from __future__ import annotations

import os
import json
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import requests

# BioPython imports — soft dependency
try:
    from Bio import PDB
    from Bio.PDB import PDBParser, PDBIO, Select
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

# ── Constants ─────────────────────────────────────────────────────────

RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_DATA_API = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"

# Contact threshold for interface analysis (Cα–Cα distance, nm)
INTERFACE_CONTACT_THRESHOLD_NM = 0.8  # 8 Å

# VH–VL interface residues in Chothia numbering on the VL side
# (key positions that pack against VH framework)
VL_FRAMEWORK_CONTACT_POSITIONS = {
    35, 36, 37, 38, 39,   # FR2 core β-strand
    44, 45, 46, 47,       # FR2 / CDR-L2 junction
    85, 86, 87,           # FR3
    98,                   # FR4
}

# Approximate Cα coordinates of key VH-pairing residues (relative, for guidance)
VH_FRAMEWORK_CONTACT_POSITIONS = {
    37, 38, 39, 40, 41,   # FR2 on VH side
    44, 45, 46, 47,       # the Q/R44 conserved residue and neighbours
    89, 90, 91,           # FR3 on VH side
}


# ── PDB utilities ────────────────────────────────────────────────────

def _check_biopython() -> None:
    if not BIOPYTHON_AVAILABLE:
        raise ImportError(
            "BioPython is required for structure analysis. "
            "Install it with: pip install biopython"
        )


def download_pdb(
    pdb_id: str,
    output_dir: str | Path = "structures",
    overwrite: bool = False,
) -> Path:
    """
    Download a PDB file from RCSB.

    Parameters
    ----------
    pdb_id : str
        4-letter PDB accession code (case-insensitive).
    output_dir : str or Path
        Directory to save the file.
    overwrite : bool
        Re-download even if file exists.

    Returns
    -------
    Path
        Path to the downloaded .pdb file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pdb_id = pdb_id.upper()
    out_path = output_dir / f"{pdb_id}.pdb"

    if out_path.exists() and not overwrite:
        print(f"  {pdb_id}: already present at {out_path}")
        return out_path

    url = RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id.lower())
    print(f"  Downloading {pdb_id} from {url} ...")
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    out_path.write_text(resp.text)
    print(f"  Saved → {out_path}  ({len(resp.text)//1024} KB)")
    return out_path


def load_structure(pdb_path: str | Path, pdb_id: str | None = None) -> "Structure":
    """
    Parse a PDB file into a BioPython Structure object.
    """
    _check_biopython()
    pdb_path = Path(pdb_path)
    model_id = pdb_id or pdb_path.stem
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(model_id, str(pdb_path))
    return structure


# ── Domain extraction ────────────────────────────────────────────────

@dataclass
class AntibodyDomain:
    """
    Extracted antibody domain (VH or VL) with metadata.

    Attributes
    ----------
    domain_type : 'VH' or 'VL' or 'VHH'
    chain_id : str
    residues : list[Residue]
    pdb_id : str
    source_antibody : str
    """
    domain_type: str
    chain_id: str
    residues: list
    pdb_id: str
    source_antibody: str = ""
    sequence: str = field(default="", repr=False)

    def __post_init__(self) -> None:
        if not self.sequence and self.residues:
            from Bio.SeqUtils import seq1
            try:
                self.sequence = "".join(
                    seq1(r.get_resname()) for r in self.residues
                    if r.get_resname() not in ("HOH", "WAT")
                )
            except Exception:
                pass

    def __len__(self) -> int:
        return len(self.residues)

    def ca_coords(self) -> np.ndarray:
        """Cα coordinates of all residues as (N, 3) array in Å."""
        coords = []
        for res in self.residues:
            if "CA" in res:
                coords.append(res["CA"].get_coord())
        return np.array(coords)

    def fasta(self, header: str | None = None) -> str:
        header = header or f"{self.pdb_id}_{self.chain_id}_{self.domain_type}"
        seq = self.sequence
        lines = [f">{header}"]
        for i in range(0, len(seq), 60):
            lines.append(seq[i:i+60])
        return "\n".join(lines)


def extract_variable_domain(
    structure: "Structure",
    chain_id: str,
    domain_type: str = "VH",
    residue_range: tuple[int, int] | None = None,
    model_idx: int = 0,
) -> AntibodyDomain:
    """
    Extract a variable domain from a PDB chain.

    Parameters
    ----------
    structure : Structure
        BioPython Structure object.
    chain_id : str
        Chain ID to extract from.
    domain_type : str
        'VH', 'VL', or 'VHH'.
    residue_range : (start, end) optional
        Residue sequence numbers to extract (inclusive). If None, take all
        residues with sequence number ≤ 130 (covers most variable domains).
    model_idx : int
        Model index (0-based).
    """
    _check_biopython()
    model = list(structure.get_models())[model_idx]
    chain = model[chain_id]

    if residue_range is None:
        residue_range = (1, 130)

    residues = [
        r for r in chain.get_residues()
        if r.get_id()[1] >= residue_range[0]
        and r.get_id()[1] <= residue_range[1]
        and r.get_id()[0] == " "  # exclude HETATMs
    ]

    return AntibodyDomain(
        domain_type=domain_type,
        chain_id=chain_id,
        residues=residues,
        pdb_id=structure.get_id(),
    )


# ── VH–VL interface analysis ─────────────────────────────────────────

@dataclass
class InterfaceContact:
    """A single contact between a VH residue and a VL residue."""
    vh_resnum: int
    vh_resname: str
    vl_resnum: int
    vl_resname: str
    ca_distance_A: float
    min_heavy_distance_A: float


def analyse_vh_vl_interface(
    vh: AntibodyDomain,
    vl: AntibodyDomain,
    contact_threshold_A: float = 8.0,
) -> list[InterfaceContact]:
    """
    Identify contacting residue pairs between VH and VL.

    A contact is defined as any heavy-atom distance < contact_threshold_A.

    Parameters
    ----------
    vh, vl : AntibodyDomain
        Extracted variable domains.
    contact_threshold_A : float
        Heavy-atom distance cutoff in Å.

    Returns
    -------
    list of InterfaceContact, sorted by min_heavy_distance ascending.
    """
    _check_biopython()
    contacts = []

    for vh_res in vh.residues:
        if vh_res.get_resname() in ("HOH", "WAT"):
            continue
        vh_atoms = list(vh_res.get_atoms())

        for vl_res in vl.residues:
            if vl_res.get_resname() in ("HOH", "WAT"):
                continue
            vl_atoms = list(vl_res.get_atoms())

            # Fast screen: Cα distance
            if "CA" not in vh_res or "CA" not in vl_res:
                continue
            ca_dist = float(np.linalg.norm(
                vh_res["CA"].get_coord() - vl_res["CA"].get_coord()
            ))
            if ca_dist > contact_threshold_A + 5:
                continue  # skip remote pairs

            # Full heavy-atom scan
            min_dist = float("inf")
            for a1 in vh_atoms:
                for a2 in vl_atoms:
                    d = float(np.linalg.norm(a1.get_coord() - a2.get_coord()))
                    if d < min_dist:
                        min_dist = d

            if min_dist <= contact_threshold_A:
                contacts.append(InterfaceContact(
                    vh_resnum=vh_res.get_id()[1],
                    vh_resname=vh_res.get_resname(),
                    vl_resnum=vl_res.get_id()[1],
                    vl_resname=vl_res.get_resname(),
                    ca_distance_A=ca_dist,
                    min_heavy_distance_A=min_dist,
                ))

    contacts.sort(key=lambda c: c.min_heavy_distance_A)
    return contacts


def interface_summary(contacts: list[InterfaceContact]) -> dict:
    """
    Summarise VH–VL interface contacts.

    Returns
    -------
    dict with:
        'n_contacts', 'vl_contact_residues', 'vh_contact_residues',
        'vl_framework_contacts', 'closest_contacts'
    """
    vl_resnums = sorted(set(c.vl_resnum for c in contacts))
    vh_resnums = sorted(set(c.vh_resnum for c in contacts))
    vl_fw_contacts = [r for r in vl_resnums if r in VL_FRAMEWORK_CONTACT_POSITIONS]

    return {
        "n_contacts": len(contacts),
        "vl_contact_residues": vl_resnums,
        "vh_contact_residues": vh_resnums,
        "vl_framework_contacts": vl_fw_contacts,
        "n_vl_framework_contacts": len(vl_fw_contacts),
        "closest_contacts": contacts[:10],
    }


def buried_sasa_estimate(contacts: list[InterfaceContact]) -> float:
    """
    Rough estimate of VL-side buried SASA (Å²) based on contact count.
    Approximately 80–120 Å² per contacting VL residue.
    """
    n_vl_residues = len(set(c.vl_resnum for c in contacts))
    return n_vl_residues * 95.0  # Å² per residue (empirical average)


# ── FASTA and sequence utilities ──────────────────────────────────────

def write_fasta(
    domains: list[AntibodyDomain],
    output_path: str | Path,
    include_metadata: bool = True,
) -> None:
    """Write multiple domains to a FASTA file."""
    lines = []
    for d in domains:
        header = f"{d.pdb_id}|{d.chain_id}|{d.domain_type}"
        if d.source_antibody:
            header += f"|{d.source_antibody}"
        if include_metadata:
            header += f"  [{len(d.residues)} residues]"
        lines.append(f">{header}")
        for i in range(0, len(d.sequence), 60):
            lines.append(d.sequence[i:i+60])
        lines.append("")
    Path(output_path).write_text("\n".join(lines))
    print(f"Wrote {len(domains)} sequences → {output_path}")


# ── Pre-defined example constructs ───────────────────────────────────

EXAMPLE_STRUCTURES = {
    "trastuzumab_fab": {
        "pdb_id": "1N8Z",
        "description": "Trastuzumab Fab (anti-HER2). VH=chain A, VL=chain B.",
        "role": "Source of VH1 + VL1 (anti-HER2 scFv, binds the TARGET)",
        "vh_chain": "A",
        "vl_chain": "B",
        "target": "HER2 (ErbB2)",
        "reference": "Cho et al. (2003) Nature 421:756",
    },
    "cetuximab_fab": {
        "pdb_id": "1YY9",
        "description": "Cetuximab Fab (anti-EGFR). VH=chain D, VL=chain C.",
        "role": (
            "Anti-EGFR VH framework reference. "
            "VH1 must be engineered to bind EGFR (ANTIGEN) "
            "while preserving VL1-pairing geometry."
        ),
        "vh_chain": "D",
        "vl_chain": "C",
        "target": "EGFR (ErbB1)",
        "reference": "Li et al. (2005) Cancer Cell 7:301",
    },
    "her2_nanobody": {
        "pdb_id": "5MY6",
        "description": "2Rs15d anti-HER2 nanobody. VHH=chain B.",
        "role": "Nanobody component (binds HER2 epitope II, distinct from trastuzumab)",
        "vh_chain": "B",
        "vl_chain": None,
        "target": "HER2 (ErbB2) domain II",
        "reference": "Desmyter et al. (2002); Kijanka et al. (2016)",
    },
}


def print_design_overview() -> None:
    """Print the example biological system design overview."""
    lines = [
        "=== Example Biological System ===",
        "",
        "Cell type : HER2+/EGFR+ cancer cells (breast, gastric, lung)",
        "",
        "  EGFR (Antigen) ─── co-expressed ─── HER2 (Target)",
        "      ↑                                     ↑",
        "  VH1 anchors                      VH1+VL1 scFv binds",
        "  to EGFR                          HER2 when VH1 anchored",
        "                                   Nanobody also binds HER2",
        "",
        "Construct:",
        "  Chain A : VH1  (anti-EGFR, anchors to antigen)",
        "  Chain B : [Minibinder]–(150-res spacer)–VL1–(G4S)60–Nanobody",
        "            └─ anti-HER2 VL ─┘              └─ anti-HER2 VHH ─┘",
        "",
        "Selectivity: ONLY activates on cells co-expressing EGFR + HER2.",
        "             Cells with EGFR only or HER2 only are NOT targeted.",
        "",
    ]
    for k, v in EXAMPLE_STRUCTURES.items():
        lines.append(f"  [{v['pdb_id']}] {v['description']}")
        lines.append(f"         Role: {v['role']}")
        lines.append("")
    print("\n".join(lines))
