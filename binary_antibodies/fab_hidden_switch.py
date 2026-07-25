"""
fab_hidden_switch.py
--------------------
Shared constants and PDB builders for the hidden trivalent minibinder switch.

Pipeline order
--------------
Stage 0 : weaken VH–VL framework interface (apo-frustrated, target-compatible)
Stage A : design VH–hub–VL minibinder into the Stage 0 Fab
Stage B : add target arm against epitope stub
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np

try:
    from Bio.PDB import PDBIO, PDBParser, Structure as S, Model as M, Chain as Ch, Residue as Res, Atom as At
except ImportError as e:
    raise ImportError("BioPython required: pip install biopython") from e

ROOT = Path(__file__).resolve().parents[1]
FAB_PDB = ROOT / "structures" / "1N8Z.pdb"

VH_END = 113
VL_END = 107
EPITOPE_SEQ = "CPACPAPELLGG"

# Chothia CDR ranges on isolated VH (A) / VL (B) chains
VH_CDR_RANGES = [(26, 35), (50, 65), (95, 102)]
VL_CDR_RANGES = [(24, 34), (50, 56), (89, 97)]

# Framework residues at the VH–VL interface (from vh_vl_contacts.csv, <4.5 Å, CDRs excluded)
VH_INTERFACE_FW = [34, 38, 42, 43, 44, 45, 46, 47, 49, 87, 89]
VL_INTERFACE_FW = [35, 37, 39, 43, 44, 45, 46, 47, 104, 105, 106, 107, 108, 109, 110, 111, 112]

# CH1 hotspots (heavy-chain numbering) for Stage A
CH1_HOTSPOTS_HEAVY = [139, 140, 142, 143, 159, 160, 161, 162, 163, 164, 165, 166, 173]
VL_HOTSPOTS_STAGE_A = [35, 37, 39, 43, 44, 45, 46, 47, 95, 99, 103, 104]


def in_ranges(resnum: int, ranges: Iterable[tuple[int, int]]) -> bool:
    return any(lo <= resnum <= hi for lo, hi in ranges)


def cdr_residue_numbers(chain: str, vh_len: int = VH_END, vl_len: int = VL_END) -> list[int]:
    ranges = VH_CDR_RANGES if chain == "A" else VL_CDR_RANGES if chain == "B" else None
    if ranges is None:
        raise ValueError(f"Unsupported chain: {chain}")
    length = vh_len if chain == "A" else vl_len
    return [r for r in range(1, length + 1) if in_ranges(r, ranges)]


def fv_framework_residue_numbers(chain: str, vh_len: int = VH_END, vl_len: int = VL_END) -> list[int]:
    """Non-CDR Fv residues (for structural alignment)."""
    length = vh_len if chain == "A" else vl_len
    ranges = VH_CDR_RANGES if chain == "A" else VL_CDR_RANGES
    return [r for r in range(1, length + 1) if not in_ranges(r, ranges)]


def cdr_fixed_atoms(vh_len: int = VH_END, vl_len: int = VL_END) -> dict[str, str]:
    """RFd3 select_fixed_atoms entries for CDR regions on chains A and B."""
    fixed: dict[str, str] = {}
    for lo, hi in VH_CDR_RANGES:
        fixed[f"A{lo}-{min(hi, vh_len)}"] = "ALL"
    for lo, hi in VL_CDR_RANGES:
        fixed[f"B{lo}-{min(hi, vl_len)}"] = "ALL"
    return fixed


def framework_design_residues(chain: str, vh_len: int = VH_END, vl_len: int = VL_END) -> list[str]:
    """All non-CDR Fv residues (for MPNN) on chain A or B."""
    if chain == "A":
        return [f"A{r}" for r in range(1, vh_len + 1) if not in_ranges(r, VH_CDR_RANGES)]
    if chain == "B":
        return [f"B{r}" for r in range(1, vl_len + 1) if not in_ranges(r, VL_CDR_RANGES)]
    raise ValueError(f"Unsupported chain: {chain}")


def epitope_hotspots(length: int | None = None) -> str:
    n = length or len(EPITOPE_SEQ)
    return ",".join(f"T{i}" for i in range(1, n + 1))


def _residues(chain, rmin: int, rmax: int) -> list:
    return [r for r in chain.get_residues() if r.id[0] == " " and rmin <= r.id[1] <= rmax]


def _copy_residue(res, new_chain: Ch.Chain, new_num: int) -> Res.Residue:
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


def _one_to_three(aa: str) -> str:
    from Bio.SeqUtils import seq3
    return seq3(aa).upper()


def _place_epitope_stub(vh_res: list, vl_res: list, seq: str = EPITOPE_SEQ) -> Ch.Chain:
    """Full-atom epitope stub at the VH–VL groove (RFd3 requires complete residues)."""
    from biotite.structure.info import residue as bt_residue

    vh_cent = _centroid(vh_res)
    vl_cent = _centroid(vl_res)
    groove = 0.5 * (vh_cent + vl_cent)
    toward_vl = vl_cent - vh_cent
    toward_vl /= np.linalg.norm(toward_vl) + 1e-8

    chain = Ch.Chain("T")
    for i, aa in enumerate(seq):
        bt_arr = bt_residue(_one_to_three(aa)).copy()
        ca_pos = groove + toward_vl * (i * 3.8) - toward_vl * 2.0
        bt_ca = bt_arr.coord[bt_arr.atom_name == "CA"][0]
        bt_arr.coord += ca_pos - bt_ca

        new_res = Res.Residue((" ", i + 1, " "), aa, " ")
        for j in range(bt_arr.array_length()):
            name = str(bt_arr.atom_name[j])
            coord = bt_arr.coord[j]
            elem = str(bt_arr.element[j])
            atom = At.Atom(name, coord, 0.0, 1.0, " ", name, i + 1, elem)
            new_res.add(atom)
        chain.add(new_res)
    return chain


def build_fab_context_pdb(
    out_pdb: Path,
    source_pdb: Path | None = None,
    struct_name: str = "fab_context",
) -> tuple[int, int, int, int]:
    """
    Build 5-chain Fab context PDB: A=VH, B=VL, C=CH1, D=CL, T=epitope stub.

    Returns (vh_len, vl_len, ch1_len, cl_len).
    """
    parser = PDBParser(QUIET=True)
    src_path = source_pdb or FAB_PDB
    src = parser.get_structure("src", str(src_path))
    model = list(src.get_models())[0]

    if {"A", "B"}.issubset(model.child_dict):
        heavy = model["A"]
        light = model["B"]
        vh_res = _residues(heavy, 1, VH_END)
        ch1_res = _residues(heavy, VH_END + 1, 9999) if len(_residues(heavy, VH_END + 1, 9999)) else []
        vl_res = _residues(light, 1, VL_END)
        cl_res = _residues(light, VL_END + 1, 9999)
        if not ch1_res and "C" in model:
            ch1_res = list(model["C"].get_residues())
        if not cl_res and "D" in model:
            cl_res = list(model["D"].get_residues())
    else:
        heavy = model["A"]
        light = model["B"]
        vh_res = _residues(heavy, 1, VH_END)
        ch1_res = _residues(heavy, VH_END + 1, 9999)
        vl_res = _residues(light, 1, VL_END)
        cl_res = _residues(light, VL_END + 1, 9999)

    struct = S.Structure(struct_name)
    model_out = M.Model(0)
    model_out.add(_build_chain(vh_res, "A"))
    model_out.add(_build_chain(vl_res, "B"))
    if ch1_res:
        model_out.add(_build_chain(ch1_res, "C"))
    if cl_res:
        model_out.add(_build_chain(cl_res, "D"))
    model_out.add(_place_epitope_stub(vh_res, vl_res, EPITOPE_SEQ))
    struct.add(model_out)

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(out_pdb))
    return len(vh_res), len(vl_res), len(ch1_res), len(cl_res)


def ch1_hotspots_chain_c(vh_len: int = VH_END) -> list[str]:
    return [f"C{r - vh_len}" for r in CH1_HOTSPOTS_HEAVY if r > vh_len]


def stage_a_hotspots(vh_len: int = VH_END) -> str:
    return ",".join(ch1_hotspots_chain_c(vh_len) + [f"B{r}" for r in VL_HOTSPOTS_STAGE_A])


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")
