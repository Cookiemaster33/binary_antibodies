#!/usr/bin/env python3
"""
build_stage_a_design_target.py
------------------------------
Build the Stage A RFd3 design target for the hidden trivalent minibinder approach.

Stage A designs a connected VH–minibinder–VL hub nestled in a native Fab,
with hotspots on CH1 (chain C) and VL (chain B). A HER2 epitope stub (chain T)
is placed inside the VH–VL groove as fixed steric context for Stage B masking.

Output
------
  structures/domains/fab_hidden_minibinder_stage_a.pdb
  structures/interface/stage_a_hidden_minibinder_config.json

Usage
-----
    python scripts/build_stage_a_design_target.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

try:
    from Bio.PDB import PDBIO, PDBParser, Structure as S, Model as M, Chain as Ch, Residue as Res, Atom as At
except ImportError as e:
    raise SystemExit("BioPython required: pip install biopython") from e

FAB_PDB = ROOT / "structures" / "1N8Z.pdb"
OUT_PDB = ROOT / "structures" / "domains" / "fab_hidden_minibinder_stage_a.pdb"
OUT_JSON = ROOT / "structures" / "interface" / "stage_a_hidden_minibinder_config.json"

# Trastuzumab / HER2 epitope fragment (domain IV, PDB 1N8Z antigen context)
EPITOPE_SEQ = "CPACPAPELLGG"

# Domain boundaries in 1N8Z (Chothia-style numbering on chain)
VH_END = 113
VL_END = 107

# VL framework hotspots at the VH–VL interface (exclude CDR loops)
VL_HOTSPOTS = [35, 37, 39, 43, 44, 45, 46, 47, 95, 99, 103, 104]

# CH1 residues on heavy chain closest to VL (within ~15 Å); mapped to chain C later
CH1_HOTSPOTS_HEAVY = [139, 140, 142, 143, 159, 160, 161, 162, 163, 164, 165, 166, 173]


def _residues(chain, rmin: int, rmax: int) -> list:
    return [r for r in chain.get_residues() if r.id[0] == " " and rmin <= r.id[1] <= rmax]


def _copy_residue(res, new_chain: Ch, new_num: int) -> Res.Residue:
    new_id = (" ", new_num, " ")
    new_res = Res.Residue(new_id, res.get_resname(), res.get_segid())
    for atom in res.get_atoms():
        new_res.add(atom.copy())
    new_chain.add(new_res)
    return new_res


def _build_chain(residues: list, chain_id: str) -> Ch.Chain:
    chain = Ch.Chain(chain_id)
    for i, res in enumerate(residues, start=1):
        _copy_residue(res, chain, i)
    return chain


def _centroid(residues: list) -> np.ndarray:
    coords = [res["CA"].get_coord() for res in residues if "CA" in res]
    return np.mean(coords, axis=0)


def _place_epitope_stub(vh_res: list, vl_res: list, seq: str) -> Ch.Chain:
    """Build a short peptide at the VH–VL groove centroid, buried toward the Fab."""
    vh_cent = _centroid(vh_res)
    vl_cent = _centroid(vl_res)
    groove = 0.5 * (vh_cent + vl_cent)

    # Direction from groove toward VL (stub sits on VL side of cleft)
    toward_vl = vl_cent - vh_cent
    toward_vl /= np.linalg.norm(toward_vl) + 1e-8

    chain = Ch.Chain("T")
    # Simple extended stub along VL direction, 3.8 Å per residue
    for i, aa in enumerate(seq):
        res = Res.Residue((" ", i + 1, " "), aa, " ")
        ca = groove + toward_vl * (i * 3.8) - toward_vl * 2.0  # pull slightly into groove
        for name, offset in [
            ("N", np.array([-1.5, 0.0, 0.0])),
            ("CA", np.array([0.0, 0.0, 0.0])),
            ("C", np.array([1.5, 0.0, 0.0])),
            ("O", np.array([2.2, 1.1, 0.0])),
        ]:
            elem = name[0] if name != "CA" else "C"
            atom = At.Atom(name, ca + offset, 0.0, 1.0, " ", name, i + 1, elem)
            res.add(atom)
        chain.add(res)
    return chain


def build_design_target() -> tuple[Path, Path]:
    parser = PDBParser(QUIET=True)
    src = parser.get_structure("1N8Z", str(FAB_PDB))
    model = list(src.get_models())[0]
    heavy = model["A"]
    light = model["B"]

    vh_res = _residues(heavy, 1, VH_END)
    ch1_res = _residues(heavy, VH_END + 1, 9999)
    vl_res = _residues(light, 1, VL_END)
    cl_res = _residues(light, VL_END + 1, 9999)

    struct = S.Structure("stage_a")
    model_out = M.Model(0)
    model_out.add(_build_chain(vh_res, "A"))
    model_out.add(_build_chain(vl_res, "B"))
    model_out.add(_build_chain(ch1_res, "C"))
    model_out.add(_build_chain(cl_res, "D"))
    model_out.add(_place_epitope_stub(vh_res, vl_res, EPITOPE_SEQ))
    struct.add(model_out)

    OUT_PDB.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(OUT_PDB))

    ch1_hotspots_c = [r - VH_END for r in CH1_HOTSPOTS_HEAVY if r > VH_END]
    config = {
        "stage": "A",
        "approach": "hidden_trivalent_minibinder",
        "description": (
            "Stage A: design connected VH–minibinder–VL hub seated in native Fab. "
            "Hotspots on CH1 (chain C) and VL (chain B). Epitope stub (chain T) "
            "provides steric masking context for Stage B target arm."
        ),
        "input_pdb": str(OUT_PDB.relative_to(ROOT)),
        "source_fab": "1N8Z (trastuzumab)",
        "chains": {
            "A": f"VH 1-{len(vh_res)}",
            "B": f"VL 1-{len(vl_res)}",
            "C": f"CH1 1-{len(ch1_res)}",
            "D": f"CL 1-{len(cl_res)} (steric context only)",
            "T": f"HER2 epitope stub 1-{len(EPITOPE_SEQ)} (steric context, Stage B)",
        },
        "rfd3": {
            "contig": f"A1-{len(vh_res)},35-55,B1-{len(vl_res)}",
            "mb_length_range": "35-55",
            "select_fixed_atoms": {
                f"C1-{len(ch1_res)}": "ALL",
                f"D1-{len(cl_res)}": "ALL",
                f"T1-{len(EPITOPE_SEQ)}": "ALL",
            },
            "hotspots_ch1_chain_c": [f"C{r}" for r in ch1_hotspots_c],
            "hotspots_vl_chain_b": [f"B{r}" for r in VL_HOTSPOTS],
            "select_hotspots": ",".join(
                [f"C{r}" for r in ch1_hotspots_c]
                + [f"B{r}" for r in VL_HOTSPOTS]
            ),
        },
        "validation": {
            "metric": "global_assembly_rmsd",
            "note": "Score Boltz refold vs RFd3 for full A–MB–B chain with C,D,T as context",
        },
        "stage_b_note": (
            "Stage B will add target-arm hotspots on chain T and extend contig "
            "to include a third binding patch against the epitope stub."
        ),
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(config, indent=2) + "\n")
    return OUT_PDB, OUT_JSON


def main() -> None:
    pdb_path, json_path = build_design_target()
    print(f"Wrote design target PDB → {pdb_path}")
    print(f"Wrote RFd3 config       → {json_path}")
    cfg = json.loads(json_path.read_text())
    print(f"  contig: {cfg['rfd3']['contig']}")
    print(f"  hotspots: {len(cfg['rfd3']['select_hotspots'].split(','))} residues")


if __name__ == "__main__":
    main()
