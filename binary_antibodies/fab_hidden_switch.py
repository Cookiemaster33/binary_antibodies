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
    from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Structure as S, Model as M, Chain as Ch, Residue as Res, Atom as At
except ImportError as e:
    raise ImportError("BioPython required: pip install biopython") from e

ROOT = Path(__file__).resolve().parents[1]
FAB_PDB = ROOT / "structures" / "1N8Z.pdb"
VH_VL_CONTACTS_CSV = ROOT / "structures" / "interface" / "vh_vl_contacts.csv"

# Tightest framework-interface pairs (min heavy-atom distance, Å) kept native during degrease.
INTERFACE_CORE_MAX_HEAVY_A = 3.35  # conservative (14 rim residues) — contact distance Å
INTERFACE_CORE_MAX_HEAVY_A_AGGRESSIVE = 3.20  # aggressive (19 rim residues)
PISA_CORE_MIN_BURIED_SASA_A2 = 10.0  # keep native on deeply buried PISA interface residues
PISA_CORE_MIN_BURIED_SASA_A2_AGGRESSIVE = 5.0
DEFAULT_SPLIT_SEPARATION_A = 30.0
# Stage A: target CH1–CL centroid distance after auto-split.
# 50 Å gives CH1-VL hotspot gap ~44 Å, clearance ~16 Å, matching the user's manually-placed
# reference file (rank079_s0_native_split…_split.cif, CH1-CL centroid ~51 Å).
DEFAULT_STAGE_A_SEPARATION_A = 50.0

VH_END = 113
VL_END = 107
EPITOPE_SEQ = "CPACPAPELLGG"

# Chothia CDR ranges on isolated VH (A) / VL (B) chains
VH_CDR_RANGES = [(26, 35), (50, 65), (95, 102)]
VL_CDR_RANGES = [(24, 34), (50, 56), (89, 97)]

# Framework residues at the VH–VL interface (from vh_vl_contacts.csv, <4.5 Å, CDRs excluded)
VH_INTERFACE_FW = [34, 38, 42, 43, 44, 45, 46, 47, 49, 87, 89]
VL_INTERFACE_FW = [35, 37, 39, 43, 44, 45, 46, 47, 104, 105, 106, 107, 108, 109, 110, 111, 112]

# CH1 hotspots (heavy-chain numbering) for Stage A.
# C35 (C-strand), C60 (DE-loop), C73 (E/F-strand) in chain-C local numbering
# → heavy-chain numbers = local + VH_END (113).
CH1_HOTSPOTS_HEAVY = [148, 173, 186]
# VL hotspots in chain-B local numbering.
# B36 (post-CDR-L1 framework), B87 (pre-CDR-L3 framework).
VL_HOTSPOTS_STAGE_A = [36, 87]

# The CH1 and VL hotspot surfaces of the split Fab sit ~30-50 A apart, so a
# minibinder that touches both has to be elongated. An alpha helix rises 1.5 A
# per residue, so a ~30-residue helix spans ~45 A: a helical hairpin / bundle of
# 60-85 residues is the topology that can actually bridge the gap. The older
# 35-55 window only fits one spanning helix plus loops.
DEFAULT_MB_LENGTH_RANGE = "60-85"
COMPACT_MB_LENGTH_RANGE = "40-55"

# Low-temperature sampler settings from the upstream RFd3 protein-binder example.
# Higher step_scale and lower gamma_0 trade diversity for designability, which
# also shifts the secondary-structure distribution towards helices.
STAGE_A_SAMPLER = {"step_scale": 3.0, "gamma_0": 0.2}


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


def _contact_rows() -> list[dict[str, str | float]]:
    if not VH_VL_CONTACTS_CSV.exists():
        return []
    import csv

    rows: list[dict[str, str | float]] = []
    with VH_VL_CONTACTS_CSV.open() as fh:
        for row in csv.DictReader(fh):
            rows.append(
                {
                    "vh_resnum": int(row["vh_resnum"]),
                    "vl_resnum": int(row["vl_resnum"]),
                    "min_heavy_distance_A": float(row["min_heavy_distance_A"]),
                }
            )
    return rows


def interface_closure_core(
    vh_len: int = VH_END,
    vl_len: int = VL_END,
    max_heavy_a: float = INTERFACE_CORE_MAX_HEAVY_A,
    vh_interface_fw: list[int] | None = None,
    vl_interface_fw: list[int] | None = None,
) -> tuple[list[int], list[int]]:
    """
    Deepest-buried VH/VL framework-interface residues — keep native sequence during
    partial de-greasing so holo closure remains plausible.
    """
    vh_set = {r for r in (vh_interface_fw or VH_INTERFACE_FW) if r <= vh_len}
    vl_set = {r for r in (vl_interface_fw or VL_INTERFACE_FW) if r <= vl_len}
    vh_core: set[int] = set()
    vl_core: set[int] = set()
    for row in _contact_rows():
        if row["min_heavy_distance_A"] > max_heavy_a:
            continue
        vh, vl = int(row["vh_resnum"]), int(row["vl_resnum"])
        if vh in vh_set:
            vh_core.add(vh)
        if vl in vl_set:
            vl_core.add(vl)
    return sorted(vh_core), sorted(vl_core)


def pisa_closure_core_pair(
    residue_bsa: dict[str, dict[int, float]],
    chain_a: str,
    chain_b: str,
    iface_a: list[int],
    iface_b: list[int],
    core_min_buried_sasa_A2: float = PISA_CORE_MIN_BURIED_SASA_A2,
) -> tuple[list[int], list[int]]:
    """Keep native sequence on PISA interface residues with highest buried SASA."""
    core_a = sorted(
        r for r in iface_a if residue_bsa.get(chain_a, {}).get(r, 0.0) >= core_min_buried_sasa_A2
    )
    core_b = sorted(
        r for r in iface_b if residue_bsa.get(chain_b, {}).get(r, 0.0) >= core_min_buried_sasa_A2
    )
    return core_a, core_b


def pisa_closure_core(
    residue_bsa: dict[str, dict[int, float]],
    vh_interface_fw: list[int],
    vl_interface_fw: list[int],
    core_min_buried_sasa_A2: float = PISA_CORE_MIN_BURIED_SASA_A2,
) -> tuple[list[int], list[int]]:
    """Keep native sequence on PISA interface residues with highest buried SASA."""
    vh_core = sorted(
        r for r in vh_interface_fw if residue_bsa.get("A", {}).get(r, 0.0) >= core_min_buried_sasa_A2
    )
    vl_core = sorted(
        r for r in vl_interface_fw if residue_bsa.get("B", {}).get(r, 0.0) >= core_min_buried_sasa_A2
    )
    return vh_core, vl_core


def interface_degrease_residues(
    chain: str,
    vh_len: int = VH_END,
    vl_len: int = VL_END,
    max_heavy_a: float = INTERFACE_CORE_MAX_HEAVY_A,
    vh_interface_fw: list[int] | None = None,
    vl_interface_fw: list[int] | None = None,
    vh_core: list[int] | None = None,
    vl_core: list[int] | None = None,
) -> list[str]:
    """Rim interface framework residues to redesign on separated Fv chains."""
    if chain not in {"A", "B"}:
        raise ValueError(f"Unsupported chain: {chain}")
    iface = (vh_interface_fw or VH_INTERFACE_FW) if chain == "A" else (vl_interface_fw or VL_INTERFACE_FW)
    length = vh_len if chain == "A" else vl_len
    if vh_core is None or vl_core is None:
        vh_core, vl_core = interface_closure_core(
            vh_len, vl_len, max_heavy_a=max_heavy_a,
            vh_interface_fw=vh_interface_fw, vl_interface_fw=vl_interface_fw,
        )
    core = vh_core if chain == "A" else vl_core
    return [f"{chain}{r}" for r in iface if r <= length and r not in core]


def split_mpnn_designed_residues(
    vh_len: int = VH_END,
    vl_len: int = VL_END,
    max_heavy_a: float = INTERFACE_CORE_MAX_HEAVY_A,
    vh_interface_fw: list[int] | None = None,
    vl_interface_fw: list[int] | None = None,
    vh_core: list[int] | None = None,
    vl_core: list[int] | None = None,
) -> list[str]:
    return interface_degrease_residues(
        "A", vh_len, vl_len, max_heavy_a,
        vh_interface_fw=vh_interface_fw, vl_interface_fw=vl_interface_fw,
        vh_core=vh_core, vl_core=vl_core,
    ) + interface_degrease_residues(
        "B", vh_len, vl_len, max_heavy_a,
        vh_interface_fw=vh_interface_fw, vl_interface_fw=vl_interface_fw,
        vh_core=vh_core, vl_core=vl_core,
    )


def interface_design_residues(
    vh_len: int = VH_END,
    vl_len: int = VL_END,
    vh_interface_fw: list[int] | None = None,
    vl_interface_fw: list[int] | None = None,
    ch1_interface_fw: list[int] | None = None,
    cl_interface_fw: list[int] | None = None,
    interface_scope: str = "fv",
) -> list[str]:
    """Framework interface positions to redesign (CDRs excluded on Fv)."""
    vh = [f"A{r}" for r in (vh_interface_fw or VH_INTERFACE_FW) if r <= vh_len]
    vl = [f"B{r}" for r in (vl_interface_fw or VL_INTERFACE_FW) if r <= vl_len]
    designed = vh + vl
    if interface_scope == "full_fab" and ch1_interface_fw and cl_interface_fw:
        designed += [f"C{r}" for r in ch1_interface_fw]
        designed += [f"D{r}" for r in cl_interface_fw]
    return designed


def split_mpnn_designed_residues_fab(
    vh_len: int,
    vl_len: int,
    ch1_len: int,
    cl_len: int,
    *,
    max_heavy_a: float = INTERFACE_CORE_MAX_HEAVY_A,
    vh_interface_fw: list[int] | None = None,
    vl_interface_fw: list[int] | None = None,
    ch1_interface_fw: list[int] | None = None,
    cl_interface_fw: list[int] | None = None,
    vh_core: list[int] | None = None,
    vl_core: list[int] | None = None,
    ch1_core: list[int] | None = None,
    cl_core: list[int] | None = None,
    interface_scope: str = "fv",
) -> list[str]:
    """Rim interface residues for partial de-grease across one or both Fab interfaces."""
    designed = split_mpnn_designed_residues(
        vh_len, vl_len, max_heavy_a=max_heavy_a,
        vh_interface_fw=vh_interface_fw, vl_interface_fw=vl_interface_fw,
        vh_core=vh_core, vl_core=vl_core,
    )
    if interface_scope != "full_fab" or not ch1_interface_fw or not cl_interface_fw:
        return designed
    ch1_core = ch1_core or []
    cl_core = cl_core or []
    designed += [f"C{r}" for r in ch1_interface_fw if r <= ch1_len and r not in ch1_core]
    designed += [f"D{r}" for r in cl_interface_fw if r <= cl_len and r not in cl_core]
    return designed


def cdr_residue_set(chain: str, vh_len: int = VH_END, vl_len: int = VL_END) -> set[int]:
    return set(cdr_residue_numbers(chain, vh_len, vl_len))


def _load_biopython_structure(path: Path | str):
    """Load PDB or mmCIF into a BioPython Structure."""
    path = Path(path)
    if path.suffix.lower() in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def extract_chain_sequences(pdb_path: Path) -> dict[str, str]:
    """One-letter sequences per chain id from a PDB/mmCIF."""
    from Bio.SeqUtils import seq1

    struct = _load_biopython_structure(pdb_path)
    model = list(struct.get_models())[0]
    out: dict[str, str] = {}
    for chain in model:
        letters: list[str] = []
        for res in chain:
            if res.id[0] != " ":
                continue
            try:
                letters.append(seq1(res.get_resname()))
            except Exception:
                letters.append("X")
        if letters:
            out[chain.id] = "".join(letters)
    return out


def _shift_residues(residues: list, shift: np.ndarray) -> list:
    """Translate all atoms in every residue by a fixed 3-D vector. Returns new copies."""
    shifted: list = []
    for res in residues:
        new_res = res.copy()
        for atom in new_res:
            atom.set_coord(atom.get_coord() + shift)
        shifted.append(new_res)
    return shifted


def _auto_split_vl_cl(
    ch1_res: list,
    cl_res: list,
    vl_res: list,
    target_ch1_cl_centroid_a: float,
) -> tuple[list, list]:
    """Translate VL+CL as one rigid arm along the CH1→CL centroid axis.

    This is the geometrically correct way to open the Fab:

    * The translation direction is the vector from the CH1 centroid to the CL
      centroid — the natural interface-opening direction.
    * VL and CL receive the **same** shift vector, so the arm stays rigid and
      the VL–CL covalent geometry is preserved.
    * Because the shift is purely translational along an axis that already lies
      in the Fab plane, the two arms remain coplanar and CH1/CL continue to
      face each other in parallel.

    The old ``_shift_residues_perpendicular`` computed the shift direction as the
    cross-product of the VH–VL axis with world-Z, which (a) changes depending on
    how the molecule is oriented in global space and (b) moves VL and CL in
    *different* directions, breaking the rigid-arm requirement.

    Returns (new_vl_res, new_cl_res).
    """
    ch1_cent = _centroid(ch1_res)
    cl_cent = _centroid(cl_res)
    axis = cl_cent - ch1_cent
    current_dist = float(np.linalg.norm(axis))
    if current_dist < 1e-3:
        raise ValueError("CH1 and CL centroids are coincident; cannot define split axis.")
    axis_norm = axis / current_dist
    shift = axis_norm * (target_ch1_cl_centroid_a - current_dist)
    return _shift_residues(vl_res, shift), _shift_residues(cl_res, shift)


def _shift_residues_perpendicular(
    residues: list,
    ref_centroid: np.ndarray,
    other_centroid: np.ndarray,
    separation_a: float,
) -> list:
    """Legacy: translate residues perpendicular to the ref→other axis.

    Kept for split-chain MPNN usage only.  Stage A arm splitting should use
    ``_auto_split_vl_cl`` instead, which preserves rigid-arm geometry.
    """
    axis = other_centroid - ref_centroid
    axis /= np.linalg.norm(axis) + 1e-8
    ref = np.array([0.0, 0.0, 1.0])
    perp = np.cross(axis, ref)
    if np.linalg.norm(perp) < 1e-6:
        perp = np.cross(axis, np.array([0.0, 1.0, 0.0]))
    perp = perp / (np.linalg.norm(perp) + 1e-8)
    shift = perp * float(separation_a)
    shifted: list = []
    for res in residues:
        new_res = res.copy()
        for atom in new_res:
            atom.set_coord(atom.get_coord() + shift)
        shifted.append(new_res)
    return shifted


def build_split_fab_mpnn_pdb(
    out_pdb: Path,
    source_pdb: Path | None = None,
    separation_a: float = DEFAULT_SPLIT_SEPARATION_A,
    struct_name: str = "split_fab_mpnn",
) -> tuple[int, int, int, int]:
    """
    Build Fab PDB with VH/VL and CH1/CL interfaces separated for split-chain MPNN.

    Chains A (VH) and C (CH1) stay fixed; B (VL) and D (CL) are translated apart.
    """
    parser = PDBParser(QUIET=True)
    src_path = source_pdb or FAB_PDB
    src = parser.get_structure("src", str(src_path))
    model = list(src.get_models())[0]
    for chain_id in ("A", "B", "C", "D"):
        if chain_id not in model.child_dict:
            raise ValueError(f"Source PDB must contain chains A–D; missing {chain_id}")

    vh_res = list(model["A"].get_residues())
    vl_res = list(model["B"].get_residues())
    ch1_res = list(model["C"].get_residues())
    cl_res = list(model["D"].get_residues())

    vl_shifted = _shift_residues_perpendicular(
        vl_res, _centroid(vh_res), _centroid(vl_res), separation_a
    )
    cl_shifted = _shift_residues_perpendicular(
        cl_res, _centroid(ch1_res), _centroid(cl_res), separation_a
    )

    struct = S.Structure(struct_name)
    model_out = M.Model(0)
    model_out.add(_build_chain(vh_res, "A"))
    model_out.add(_build_chain(vl_shifted, "B"))
    model_out.add(_build_chain(ch1_res, "C"))
    model_out.add(_build_chain(cl_shifted, "D"))
    struct.add(model_out)

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(out_pdb))
    return len(vh_res), len(vl_res), len(ch1_res), len(cl_res)


def build_split_fv_mpnn_pdb(
    out_pdb: Path,
    source_pdb: Path | None = None,
    separation_a: float = DEFAULT_SPLIT_SEPARATION_A,
    struct_name: str = "split_fv_mpnn",
) -> tuple[int, int]:
    """
    Build Fv-only PDB with VH (A) and VL (B) translated apart for split-chain MPNN.

    VL is shifted along the axis perpendicular to the VH→VL vector so interface
    residues become solvent-exposed on both chains.
    """
    parser = PDBParser(QUIET=True)
    src_path = source_pdb or FAB_PDB
    src = parser.get_structure("src", str(src_path))
    model = list(src.get_models())[0]

    if {"A", "B"}.issubset(model.child_dict):
        heavy = model["A"]
        light = model["B"]
        vh_res = _residues(heavy, 1, VH_END)
        vl_res = _residues(light, 1, VL_END)
    else:
        raise ValueError("Source PDB must contain chains A (VH) and B (VL)")

    vh_cent = _centroid(vh_res)
    vl_cent = _centroid(vl_res)
    axis = vl_cent - vh_cent
    axis /= np.linalg.norm(axis) + 1e-8
    # Perpendicular shift (arbitrary but stable): cross with global Z unless parallel.
    ref = np.array([0.0, 0.0, 1.0])
    perp = np.cross(axis, ref)
    if np.linalg.norm(perp) < 1e-6:
        perp = np.cross(axis, np.array([0.0, 1.0, 0.0]))
    perp = perp / (np.linalg.norm(perp) + 1e-8)
    shift = perp * float(separation_a)

    vl_shifted: list = []
    for res in vl_res:
        new_res = res.copy()
        for atom in new_res:
            atom.set_coord(atom.get_coord() + shift)
        vl_shifted.append(new_res)

    struct = S.Structure(struct_name)
    model_out = M.Model(0)
    model_out.add(_build_chain(vh_res, "A"))
    model_out.add(_build_chain(vl_shifted, "B"))
    struct.add(model_out)

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(out_pdb))
    return len(vh_res), len(vl_res)


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

        resname = _one_to_three(aa)
        new_res = Res.Residue((" ", i + 1, " "), resname, " ")
        for j in range(bt_arr.array_length()):
            name = str(bt_arr.atom_name[j])
            elem = str(bt_arr.element[j])
            if name == "OXT" or elem == "H":
                continue
            coord = bt_arr.coord[j]
            atom = At.Atom(name, coord, 0.0, 1.0, " ", name, i + 1, elem)
            new_res.add(atom)
        chain.add(new_res)
    return chain


def _chain_nterm_sequence(chain, n: int = 12) -> str:
    from Bio.SeqUtils import seq1

    letters: list[str] = []
    for res in chain:
        if res.id[0] != " ":
            continue
        if len(letters) >= n:
            break
        try:
            letters.append(seq1(res.get_resname()))
        except Exception:
            letters.append("X")
    return "".join(letters)


def _looks_like_vh_nterm(seq: str) -> bool:
    """Heavy-chain Fv N-terminus (trastuzumab / human IgG1: EVQLVES...)."""
    return seq.startswith(("EVQL", "QVQL", "QMQL"))


def _identify_boltz_ig_chain_ids(model) -> tuple[str, str]:
    """
    Return (immunoglobulin_heavy_chain_id, immunoglobulin_light_chain_id) for Boltz H/L.

    Boltz YAML uses id H/L, but some holo outputs assign H=light (VL+CL) and L=heavy
    (VH+CH1). Detect from Fv N-terminal sequence rather than chain letter.
    """
    if "H" not in model.child_dict or "L" not in model.child_dict:
        raise ValueError("Expected Boltz chains H and L")
    h_n = _chain_nterm_sequence(model["H"])
    l_n = _chain_nterm_sequence(model["L"])
    h_is_vh = _looks_like_vh_nterm(h_n)
    l_is_vh = _looks_like_vh_nterm(l_n)
    if l_is_vh and not h_is_vh:
        return "L", "H"
    if h_is_vh and not l_is_vh:
        return "H", "L"
    # Fallback: canonical fuse convention H=VH+CH1, L=VL+CL
    return "H", "L"


def _load_fab_chain_residues(
    source_pdb: Path | None = None,
) -> tuple[list, list, list, list, Ch.Chain | None]:
    """Return (vh, vl, ch1, cl, epitope_chain_or_none) from Fab PDB/CIF."""
    src_path = source_pdb or FAB_PDB
    src = _load_biopython_structure(src_path)
    model = list(src.get_models())[0]
    chains = set(model.child_dict)
    epitope_from_source: Ch.Chain | None = None

    if "H" in chains and "L" in chains:
        ig_heavy_id, ig_light_id = _identify_boltz_ig_chain_ids(model)
        heavy = model[ig_heavy_id]
        light = model[ig_light_id]
        vh_res = _residues(heavy, 1, VH_END)
        ch1_res = _residues(heavy, VH_END + 1, 9999)
        vl_res = _residues(light, 1, VL_END)
        cl_res = _residues(light, VL_END + 1, 9999)
        if "T" in chains:
            epitope_from_source = model["T"]
    elif {"A", "B"}.issubset(chains):
        heavy = model["A"]
        light = model["B"]
        vh_res = _residues(heavy, 1, VH_END)
        ch1_res = _residues(heavy, VH_END + 1, 9999) if len(_residues(heavy, VH_END + 1, 9999)) else []
        vl_res = _residues(light, 1, VL_END)
        cl_res = _residues(light, VL_END + 1, 9999)
        if not ch1_res and "C" in chains:
            ch1_res = [r for r in model["C"].get_residues() if r.id[0] == " "]
        if not cl_res and "D" in chains:
            cl_res = [r for r in model["D"].get_residues() if r.id[0] == " "]
        if "T" in chains:
            epitope_from_source = model["T"]
    else:
        raise ValueError("Source must contain fused chains H/L or split chains A–D")

    return vh_res, vl_res, ch1_res, cl_res, epitope_from_source


def _write_structure(struct: S.Structure, out_path: Path) -> None:
    """Write BioPython structure to PDB or mmCIF based on suffix."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.suffix.lower() in {".cif", ".mmcif"}:
        from Bio.PDB import MMCIFIO

        mmcif_io = MMCIFIO()
        mmcif_io.set_structure(struct)
        mmcif_io.save(str(out_path))
    else:
        io = PDBIO()
        io.set_structure(struct)
        io.save(str(out_path))


def fab_has_split_chains(source_pdb: Path | str) -> bool:
    """True when structure already uses logical Fab chains A–D (not fused H/L)."""
    src = _load_biopython_structure(source_pdb)
    model = list(src.get_models())[0]
    chains = set(model.child_dict)
    return {"A", "B"}.issubset(chains) and "H" not in chains and "L" not in chains


def _is_hl_format(model) -> bool:
    """True when the structure uses the two-arm H/L chain naming.

    H = VH+CH1 (one rigid arm), L = VL+CL (the other rigid arm).  This is the
    format the user produces when manually repositioning the Fab arms in a
    molecular-graphics tool before Stage A.  It is distinct from the fused Boltz
    H/L output (which has the arms packed together) and from the legacy A/B/C/D
    split (four separate domain chains).
    """
    chains = set(model.child_dict)
    return "H" in chains and "L" in chains and not {"A", "B", "C", "D"}.intersection(chains)


def fab_arms_are_separated(source_pdb: Path | str, min_gap_a: float = 8.0) -> bool:
    """True when the VH/CH1 and VL/CL arms are already pulled apart in the input.

    Works with both formats:
    - A/B/C/D split: arm 1 = chains A+C, arm 2 = chains B+D
    - H/L two-arm: arm 1 = chain H (VH+CH1), arm 2 = chain L (VL+CL)

    An assembled Fab has packed VH-VL and CH1-CL interfaces (~4-5 Å contacts), so
    any gap beyond ``min_gap_a`` means the arms have been separated deliberately.
    """
    model = list(_load_biopython_structure(source_pdb).get_models())[0]

    def arm_ca(chain_ids: tuple[str, ...]) -> np.ndarray:
        coords = [
            res["CA"].get_coord()
            for cid in chain_ids
            if cid in model.child_dict
            for res in model[cid]
            if res.id[0] == " " and "CA" in res
        ]
        return np.array(coords)

    if _is_hl_format(model):
        arm1, arm2 = arm_ca(("H",)), arm_ca(("L",))
    else:
        arm1, arm2 = arm_ca(("A", "C")), arm_ca(("B", "D"))

    if len(arm1) == 0 or len(arm2) == 0:
        return False
    gap = float(np.linalg.norm(arm1[:, None, :] - arm2[None, :, :], axis=2).min())
    return gap >= min_gap_a


def infer_already_split(source_pdb: Path | str | None) -> bool:
    """True when Fab coordinates should be used as-is (no VL/CL re-translation).

    Handles both A/B/C/D split files and H/L two-arm files.  Detection is based
    on the structure's own geometry, not its filename: any file where the two Fab
    arms are already separated is a deliberate placement (e.g. arms repositioned
    by hand in PyMOL or ChimeraX) and must be passed through verbatim.
    """
    if source_pdb is None:
        return False
    model = list(_load_biopython_structure(source_pdb).get_models())[0]
    if _is_hl_format(model):
        # H/L with arms pulled apart is always a hand-placed file.
        return fab_arms_are_separated(source_pdb)
    if not fab_has_split_chains(source_pdb):
        return False
    return fab_arms_are_separated(source_pdb)


def verify_design_target_matches_source(
    design_pdb: Path | str, source_pdb: Path | str, tol_a: float = 1e-3
) -> dict[str, float]:
    """Assert the built design target reproduces the source coordinates exactly.

    Guards the promise that a hand-placed Fab is passed through untouched.  Raises
    if any CA moves by more than ``tol_a``.  Returns the per-chain maximum deviation.

    Handles both source formats:
    - A/B/C/D source: chains are compared directly by ID.
    - H/L source: H is split at VH_END into virtual chains A (VH) and C (CH1);
      L is split at VL_END into virtual chains B (VL) and D (CL).  The design
      target (always A/B/C/D) is compared against these virtual chains.
    """
    def _ca_by_chain(path: Path | str) -> dict[str, np.ndarray]:
        model = list(_load_biopython_structure(path).get_models())[0]
        return {
            chain.id: np.array(
                [res["CA"].get_coord() for res in chain if res.id[0] == " " and "CA" in res]
            )
            for chain in model.get_chains()
        }

    design = _ca_by_chain(design_pdb)
    source_raw = _ca_by_chain(source_pdb)

    # For H/L source, synthesise virtual A/B/C/D CAs from the two-arm chains.
    # Use _identify_boltz_ig_chain_ids to correctly handle the case where the
    # user's CIF has H=light-arm and L=heavy-arm (e.g. the chain labelled "L"
    # actually carries EVQL/VH at its N-terminus).
    if "H" in source_raw and "L" in source_raw and "A" not in source_raw:
        src_model = list(_load_biopython_structure(source_pdb).get_models())[0]
        ig_heavy_id, ig_light_id = _identify_boltz_ig_chain_ids(src_model)
        heavy_ca = source_raw[ig_heavy_id]   # VH+CH1 arm
        light_ca = source_raw[ig_light_id]   # VL+CL arm
        source: dict[str, np.ndarray] = {
            "A": heavy_ca[:VH_END],           # VH → design chain A
            "C": heavy_ca[VH_END:],           # CH1 → design chain C
            "B": light_ca[:VL_END],           # VL → design chain B
            "D": light_ca[VL_END:],           # CL → design chain D
        }
    else:
        source = source_raw

    deviations: dict[str, float] = {}
    problems: list[str] = []
    for cid, coords in design.items():
        if cid == "T":
            # Epitope stub is generated, not taken from the source.
            continue
        ref = source.get(cid)
        if ref is None or len(ref) != len(coords):
            if ref is not None:
                problems.append(
                    f"chain {cid}: length mismatch (design {len(coords)} vs source {len(ref)})"
                )
            continue
        dev = float(np.abs(coords - ref).max()) if len(coords) else 0.0
        deviations[cid] = round(dev, 6)
        if dev > tol_a:
            problems.append(f"chain {cid}: moved by {dev:.3f} Å")
    if problems:
        raise ValueError(
            "Design target does not reproduce the source Fab placement: "
            + "; ".join(problems)
        )
    return deviations


def _assemble_fab_structure(
    vh_res: list,
    vl_res: list,
    ch1_res: list,
    cl_res: list,
    epitope_from_source: Ch.Chain | None,
    struct_name: str,
    *,
    place_stub_if_missing: bool = True,
) -> S.Structure:
    struct = S.Structure(struct_name)
    model_out = M.Model(0)
    model_out.add(_build_chain(vh_res, "A"))
    model_out.add(_build_chain(vl_res, "B"))
    if ch1_res:
        model_out.add(_build_chain(ch1_res, "C"))
    if cl_res:
        model_out.add(_build_chain(cl_res, "D"))
    if epitope_from_source is not None:
        t_res = [r for r in epitope_from_source.get_residues() if r.id[0] == " "]
        model_out.add(_build_chain(t_res, "T"))
    elif place_stub_if_missing:
        model_out.add(_place_epitope_stub(vh_res, vl_res, EPITOPE_SEQ))
    struct.add(model_out)
    return struct


def build_holo_split_cif(
    out_path: Path,
    holo_cif: Path,
    separation_a: float = DEFAULT_STAGE_A_SEPARATION_A,
    struct_name: str = "holo_split",
) -> tuple[int, int, int, int]:
    """
    Expand fused Boltz holo (H/L[/T]) to A/B/C/D[/T] and separate VL/CL from VH/CH1.

    ``separation_a`` is the target CH1–CL centroid distance in Å after splitting.
    VL and CL are translated together as one rigid arm along the CH1→CL axis so
    the arms remain coplanar and CH1/CL continue to face each other.

    Writes a *_split.cif (or .pdb) suitable for PyMOL inspection and Stage A --fab-pdb.
    """
    vh_res, vl_res, ch1_res, cl_res, epitope_from_source = _load_fab_chain_residues(holo_cif)
    vl_res, cl_res = _auto_split_vl_cl(ch1_res, cl_res, vl_res, separation_a)
    struct = _assemble_fab_structure(
        vh_res, vl_res, ch1_res, cl_res, epitope_from_source, struct_name, place_stub_if_missing=True
    )
    _write_structure(struct, out_path)
    return len(vh_res), len(vl_res), len(ch1_res), len(cl_res)


DEFAULT_AUTO_SPLIT_DISTANCES: tuple[float, ...] = (30.0, 40.0, 50.0, 60.0)
"""CH1–CL centroid separations (Å) tested by default in auto_split_fab_arms.

The native assembled Fab has ~17 Å centroid separation. The user's manual
placement (rank079 v2 CIF) sits at ~51 Å centroid / 21 Å closest approach.
The range 30–60 Å spans from a narrow gap (short minibinder ≈ 40 aa helical
hairpin) to a wide gap (long minibinder ≈ 75 aa three-helix bundle).
"""


def auto_split_fab_arms(
    source_pdb: Path | str,
    out_dir: Path,
    distances_a: Iterable[float] = DEFAULT_AUTO_SPLIT_DISTANCES,
    prefix: str = "fab_split",
    struct_name_prefix: str = "split",
) -> list[dict]:
    """Generate A/B/C/D split PDBs at multiple CH1–CL centroid distances.

    For each distance ``d`` in ``distances_a``, writes
    ``<out_dir>/<prefix>_<d>a.pdb`` and returns a list of dicts with keys:

    * ``distance_a``: the requested CH1–CL centroid target (Å)
    * ``actual_centroid_a``: measured CH1–CL centroid after translation
    * ``closest_approach_a``: minimum CA–CA distance between CH1 and CL
    * ``ch1_vl_centroid_a``: gap between CH1 and VL hotspot centroids
    * ``path``: Path of the written PDB

    The VH+CH1 arm is fixed; VL+CL translate together along the CH1→CL axis,
    so arms stay coplanar and CH1/CL remain facing each other.

    Accepts both fused H/L and already-split A/B/C/D sources (the latter are
    re-split from their current geometry to the requested distances).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    vh_res, vl_res, ch1_res, cl_res, epitope_from_source = _load_fab_chain_residues(source_pdb)

    results = []
    for d in distances_a:
        d = float(d)
        new_vl, new_cl = _auto_split_vl_cl(ch1_res, cl_res, vl_res, d)

        ch1_cent = _centroid(ch1_res)
        cl_cent = _centroid(new_cl)
        actual_centroid = float(np.linalg.norm(ch1_cent - cl_cent))

        ch1_ca = np.array([r["CA"].get_coord() for r in ch1_res if "CA" in r])
        cl_ca = np.array([r["CA"].get_coord() for r in new_cl if "CA" in r])
        closest = float(np.linalg.norm(ch1_ca[:, None, :] - cl_ca[None, :, :], axis=2).min())

        vl_ca = np.array([r["CA"].get_coord() for r in new_vl if "CA" in r])
        ch1_vl_cent = float(np.linalg.norm(ch1_cent - vl_ca.mean(0)))

        struct = _assemble_fab_structure(
            vh_res, new_vl, ch1_res, new_cl, None,
            f"{struct_name_prefix}_{int(d)}a",
            place_stub_if_missing=False,
        )
        out_path = out_dir / f"{prefix}_{int(d)}a.pdb"
        _write_structure(struct, out_path)

        results.append({
            "distance_a": d,
            "actual_centroid_a": round(actual_centroid, 2),
            "closest_approach_a": round(closest, 2),
            "ch1_vl_centroid_a": round(ch1_vl_cent, 2),
            "path": out_path,
        })

    return results


def build_fab_context_pdb(
    out_pdb: Path,
    source_pdb: Path | None = None,
    struct_name: str = "fab_context",
) -> tuple[int, int, int, int]:
    """
    Build 5-chain Fab context PDB: A=VH, B=VL, C=CH1, D=CL, T=epitope stub.

    Returns (vh_len, vl_len, ch1_len, cl_len).
    """
    vh_res, vl_res, ch1_res, cl_res, epitope_from_source = _load_fab_chain_residues(source_pdb)
    struct = _assemble_fab_structure(
        vh_res, vl_res, ch1_res, cl_res, epitope_from_source, struct_name
    )
    _write_structure(struct, out_pdb)
    return len(vh_res), len(vl_res), len(ch1_res), len(cl_res)


def build_stage_a_design_target_pdb(
    out_pdb: Path,
    source_pdb: Path | None = None,
    separation_a: float = DEFAULT_STAGE_A_SEPARATION_A,
    struct_name: str = "stage_a",
    *,
    already_split: bool = False,
    include_epitope: bool = False,
) -> tuple[int, int, int, int]:
    """
    Build Stage A RFd3 target: Fab chains A–D with VL/CL translated apart.

    The epitope stub (chain T) is **excluded by default**.  In the split
    configuration the stub is geometrically displaced far from the inter-arm
    gap (it was synthesised at the VH–VL groove of the assembled Fab), so it
    contributes no useful context to RFd3 and can mislead the origin placement.
    Pass ``include_epitope=True`` only for debugging or legacy compatibility.

    ``separation_a`` is the **target CH1–CL centroid distance** in Å.  VL and CL
    are moved together as a rigid arm along the CH1→CL axis so the arms stay
    coplanar and the CH1/CL interface faces remain parallel.

    Pass a pre-built *_split.cif (or H/L CIF) with ``already_split=True`` to use
    coordinates verbatim without any re-translation.
    """
    vh_res, vl_res, ch1_res, cl_res, epitope_from_source = _load_fab_chain_residues(source_pdb)

    if not already_split and separation_a > 0:
        vl_res, cl_res = _auto_split_vl_cl(ch1_res, cl_res, vl_res, separation_a)

    epitope = epitope_from_source if include_epitope else None
    struct = _assemble_fab_structure(
        vh_res, vl_res, ch1_res, cl_res, epitope, struct_name,
        place_stub_if_missing=False,
    )
    _write_structure(struct, out_pdb)
    return len(vh_res), len(vl_res), len(ch1_res), len(cl_res)


def stage_a_contig(
    vh_len: int,
    vl_len: int,
    ch1_len: int,
    cl_len: int,
    mb_length_range: str = DEFAULT_MB_LENGTH_RANGE,
    epitope_len: int = 0,
) -> str:
    """RFd3 contig in the canonical binder-design layout.

    The designed minibinder comes first as its own chain, then a chain break,
    then every Fab chain as a separate fixed target chain:

        ``60-85,/0,A1-113,/0,B1-107,/0,C1-107,/0,D1-107``

    Chain T (epitope stub) is omitted by default: in the split Fab geometry it
    sits far from the inter-arm gap and adds no useful RFd3 context.

    This matters. A designed segment that sits *between* two motif segments in a
    contig (e.g. ``B1-107/0,35-55,C1-107``) is covalently bonded to both of them,
    so RFd3 has to translate and rotate the Fab domains until a single polymer can
    physically connect them - that is what moved the manually placed arms in
    earlier runs - and it forces the "minibinder" to be an extended linker rather
    than a folded binder. Keeping the minibinder on its own chain removes both
    failure modes.
    """
    fab_chains = [f"A1-{vh_len}", f"B1-{vl_len}", f"C1-{ch1_len}", f"D1-{cl_len}"]
    if epitope_len > 0:
        fab_chains.append(f"T1-{epitope_len}")
    return f"{mb_length_range},/0," + ",/0,".join(fab_chains)


def stage_a_rfd3_config(
    vh_len: int,
    vl_len: int,
    ch1_len: int,
    cl_len: int,
    mb_length_range: str = DEFAULT_MB_LENGTH_RANGE,
    epitope_len: int = 0,
    *,
    ori_token: list[float] | None = None,
    is_non_loopy: bool = True,
) -> dict[str, object]:
    """RFd3 inputs that keep the full Fab fixed in 3D while designing a helical MB.

    ``is_non_loopy`` is RFd3's only secondary-structure lever (there is no
    ``select_ss``); upstream reports it yields "a lot more helices and fewer loops
    (and less sheets)". Pair it with the low-temperature sampler settings in
    :data:`STAGE_A_SAMPLER` for helix-rich, designable backbones.
    """
    cfg: dict[str, object] = {
        "dialect": 2,
        "contig": stage_a_contig(vh_len, vl_len, ch1_len, cl_len, mb_length_range, epitope_len),
        "select_fixed_atoms": stage_a_fixed_atoms(vh_len, vl_len, ch1_len, cl_len, epitope_len),
        "select_hotspots": stage_a_hotspots(vh_len),
        "is_non_loopy": is_non_loopy,
        "mb_length_range": mb_length_range,
        "sampler": dict(STAGE_A_SAMPLER),
    }
    if ori_token is not None:
        # Place the designed chain's centre of mass in the inter-arm gap. The
        # documented "hotspots" strategy offsets the origin 10 A *outward* from the
        # hotspot centroid, which would push the minibinder out of the gap, so pin
        # the origin explicitly instead.
        cfg["ori_token"] = [round(float(v), 3) for v in ori_token]
    else:
        cfg["infer_ori_strategy"] = "hotspots"
    return cfg


def stage_a_fixed_atoms(
    vh_len: int,
    vl_len: int,
    ch1_len: int,
    cl_len: int,
    epitope_len: int = 0,
) -> dict[str, str]:
    """All Fab chains fixed; only the unlinked minibinder is designed.

    Chain T (epitope stub) is excluded by default for Stage A.
    """
    fixed = {
        f"A1-{vh_len}": "ALL",
        f"B1-{vl_len}": "ALL",
        f"C1-{ch1_len}": "ALL",
        f"D1-{cl_len}": "ALL",
    }
    if epitope_len > 0:
        fixed[f"T1-{epitope_len}"] = "ALL"
    return fixed


def ch1_hotspots_chain_c(vh_len: int = VH_END) -> list[str]:
    return [f"C{r - vh_len}" for r in CH1_HOTSPOTS_HEAVY if r > vh_len]


def stage_a_hotspots(vh_len: int = VH_END) -> str:
    return ",".join(ch1_hotspots_chain_c(vh_len) + [f"B{r}" for r in VL_HOTSPOTS_STAGE_A])


def stage_a_hotspot_geometry(design_pdb: Path, vh_len: int = VH_END) -> dict[str, object]:
    """Measure the CH1 <-> VL gap the minibinder has to bridge.

    Returns the two hotspot centroids, their separation, the midpoint (used as the
    RFd3 ``ori_token``) and the closest approach between the two hotspot surfaces.
    """
    model = list(_load_biopython_structure(design_pdb).get_models())[0]

    def hotspot_ca(chain_id: str, resnums: Iterable[int]) -> np.ndarray:
        by_num = {
            res.id[1]: res["CA"].get_coord()
            for res in model[chain_id]
            if res.id[0] == " " and "CA" in res
        }
        found = [by_num[n] for n in resnums if n in by_num]
        if not found:
            raise ValueError(f"No hotspot CA atoms found in chain {chain_id} of {design_pdb}")
        return np.array(found)

    ch1_ca = hotspot_ca("C", [r - vh_len for r in CH1_HOTSPOTS_HEAVY if r > vh_len])
    vl_ca = hotspot_ca("B", VL_HOTSPOTS_STAGE_A)
    ch1_centroid, vl_centroid = ch1_ca.mean(0), vl_ca.mean(0)
    pair_dists = np.linalg.norm(ch1_ca[:, None, :] - vl_ca[None, :, :], axis=2)

    # Full CH1-CL centroid separation (arm-opening metric, independent of hotspot selection).
    ch1_all_ca = np.array([
        res["CA"].get_coord() for res in model["C"] if res.id[0] == " " and "CA" in res
    ]) if "C" in model.child_dict else np.empty((0, 3))
    cl_all_ca = np.array([
        res["CA"].get_coord() for res in model["D"] if res.id[0] == " " and "CA" in res
    ]) if "D" in model.child_dict else np.empty((0, 3))
    ch1_cl_centroid_a = (
        round(float(np.linalg.norm(ch1_all_ca.mean(0) - cl_all_ca.mean(0))), 2)
        if len(ch1_all_ca) and len(cl_all_ca) else None
    )

    return {
        "ch1_hotspot_centroid": [round(float(v), 3) for v in ch1_centroid],
        "vl_hotspot_centroid": [round(float(v), 3) for v in vl_centroid],
        "centroid_separation_a": round(float(np.linalg.norm(ch1_centroid - vl_centroid)), 2),
        "midpoint": [round(float(v), 3) for v in (ch1_centroid + vl_centroid) / 2.0],
        "closest_hotspot_approach_a": round(float(pair_dists.min()), 2),
        "n_ch1_hotspots": int(len(ch1_ca)),
        "n_vl_hotspots": int(len(vl_ca)),
        "ch1_cl_centroid_a": ch1_cl_centroid_a,
    }


def stage_a_gap_centre(
    design_pdb: Path,
    vh_len: int = VH_END,
    *,
    axis_range: tuple[float, float] = (0.2, 0.8),
    max_perpendicular_a: float = 10.0,
) -> dict[str, object]:
    """Find the roomiest point in the CH1 <-> VL gap for the designed chain's origin.

    The plain midpoint of the two hotspot centroids lies on the line joining them,
    which for this Fab passes within ~4.6 A of Fab atoms - effectively on the
    protein surface. Search along that axis and a bounded perpendicular disc for
    the point with the largest clearance to any Fab heavy atom, so the minibinder
    starts in open space between the arms rather than buried against one of them.
    """
    geometry = stage_a_hotspot_geometry(design_pdb, vh_len)
    start = np.array(geometry["ch1_hotspot_centroid"], dtype=float)
    end = np.array(geometry["vl_hotspot_centroid"], dtype=float)
    axis = end - start
    axis_unit = axis / np.linalg.norm(axis)

    helper = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(helper, axis_unit))) > 0.9:
        helper = np.array([1.0, 0.0, 0.0])
    perp1 = np.cross(axis_unit, helper)
    perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(axis_unit, perp1)

    model = list(_load_biopython_structure(design_pdb).get_models())[0]
    fab = np.array(
        [
            atom.get_coord()
            for chain in model.get_chains()
            for res in chain
            if res.id[0] == " "
            for atom in res
            if atom.element != "H"
        ]
    )

    def clearance(point: np.ndarray) -> float:
        return float(np.linalg.norm(fab - point, axis=1).min())

    best_point = start + 0.5 * axis
    best_clearance = clearance(best_point)
    best_offset = 0.0
    for t in np.linspace(axis_range[0], axis_range[1], 25):
        base = start + t * axis
        for radius in np.linspace(0.0, max_perpendicular_a, 6):
            angles = [0.0] if radius == 0.0 else np.linspace(0, 2 * np.pi, 12, endpoint=False)
            for angle in angles:
                point = base + radius * (np.cos(angle) * perp1 + np.sin(angle) * perp2)
                value = clearance(point)
                # Prefer clearance, but break ties towards the axis so the origin
                # stays between the two hotspot surfaces.
                if value > best_clearance + 1e-6 or (
                    abs(value - best_clearance) <= 1e-6 and radius < best_offset
                ):
                    best_point, best_clearance, best_offset = point, value, float(radius)

    return {
        "ori_token": [round(float(v), 3) for v in best_point],
        "clearance_a": round(best_clearance, 2),
        "perpendicular_offset_a": round(best_offset, 2),
        "midpoint_clearance_a": round(clearance(start + 0.5 * axis), 2),
        **geometry,
    }


def stage_a_ori_token(design_pdb: Path, vh_len: int = VH_END) -> list[float]:
    """RFd3 origin token: the roomiest point in the gap between CH1 and VL."""
    return list(stage_a_gap_centre(design_pdb, vh_len)["ori_token"])  # type: ignore[arg-type]


def _kabsch_transform(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (R, mobile_centroid, target_centroid) mapping mobile → target frame."""
    mc = mobile - mobile.mean(0)
    tc = target - target.mean(0)
    v, _, wt = np.linalg.svd(mc.T @ tc)
    d = np.sign(np.linalg.det(v @ wt))
    r = v @ np.diag([1.0, 1.0, d]) @ wt
    return r, mobile.mean(0), target.mean(0)


def _read_output_chain_residues(cif_path: Path) -> tuple[list, str]:
    """Return (residues, one-letter sequence) from RFd3 single-chain CIF output."""
    from Bio.SeqUtils import seq1

    struct = _load_biopython_structure(cif_path)
    chain = next(list(struct.get_models())[0].get_chains())
    residues = [r for r in chain if r.id[0] == " "]
    seq = "".join(seq1(r.get_resname()) for r in residues)
    return residues, seq


def _infer_mb_segment(
    seq_len: int,
    vh_len: int,
    vl_len: int,
    ch1_len: int,
    cl_len: int,
    mb_min: int = 20,
    mb_max: int = 120,
) -> tuple[int, int, int]:
    """Return (mb_start, mb_end, mb_len) for RFd3 output (0-based residue indices)."""
    fixed = vh_len + vl_len + ch1_len + cl_len
    if seq_len >= fixed + mb_min:
        mb_len = seq_len - fixed
        if mb_min <= mb_len <= mb_max:
            return vh_len + vl_len, vh_len + vl_len + mb_len, mb_len
    compact = vl_len + ch1_len
    if seq_len >= compact + mb_min:
        mb_len = seq_len - compact
        if mb_min <= mb_len <= mb_max:
            return vl_len, vl_len + mb_len, mb_len
    raise ValueError(f"Cannot infer minibinder segment from output length {seq_len}")


def _ordered_output_chains(cif_path: Path) -> list[tuple[str, list]]:
    """Output chains in file order, as (chain_id, standard residues)."""
    model = list(_load_biopython_structure(cif_path).get_models())[0]
    return [
        (chain.id, [res for res in chain if res.id[0] == " "])
        for chain in model.get_chains()
        if any(res.id[0] == " " for res in chain)
    ]


def _split_minibinder_chain(
    chains: list[tuple[str, list]],
    fab_lengths: list[int],
    fab_centroids: list[np.ndarray] | None = None,
) -> tuple[list, list[list]]:
    """Separate the designed chain from the fixed Fab chains in a multi-chain output.

    The Fab chains normally come back in contig order (A, B, C, D, T), so the
    designed chain is the one whose removal leaves exactly the expected Fab residue
    counts in order. If RFd3 ever reorders chains, fall back to matching each
    output chain to the nearest input chain of equal length by centroid - the Fab
    is returned in the input frame, so centroids identify chains unambiguously even
    though VL, CH1 and CL all have 107 residues.
    """
    for idx in range(len(chains)):
        rest = [residues for i, (_cid, residues) in enumerate(chains) if i != idx]
        if [len(r) for r in rest] == fab_lengths:
            return chains[idx][1], rest

    if fab_centroids is not None and len(chains) == len(fab_lengths) + 1:
        def centroid(residues: list) -> np.ndarray:
            return np.array(
                [res["CA"].get_coord() for res in residues if "CA" in res]
            ).mean(0)

        for idx in range(len(chains)):
            candidates = [residues for i, (_cid, residues) in enumerate(chains) if i != idx]
            remaining = list(candidates)
            ordered: list[list] = []
            ok = True
            for want_len, want_centroid in zip(fab_lengths, fab_centroids):
                matches = [r for r in remaining if len(r) == want_len]
                if not matches:
                    ok = False
                    break
                best = min(matches, key=lambda r: float(np.linalg.norm(centroid(r) - want_centroid)))
                ordered.append(best)
                remaining.remove(best)
            if ok and not remaining:
                return chains[idx][1], ordered

    observed = {cid: len(residues) for cid, residues in chains}
    raise ValueError(
        f"Cannot identify designed chain; output chains {observed} do not leave "
        f"Fab counts {fab_lengths} after removing one chain"
    )


def graft_stage_a_minibinder(
    out_path: Path,
    input_pdb: Path,
    rfd3_cif: Path,
    vh_len: int = VH_END,
    vl_len: int = VL_END,
    ch1_len: int | None = None,
    cl_len: int | None = None,
    mb_chain_id: str = "M",
) -> dict[str, object]:
    """
    Build output = exact input Fab (A,B,C,D,T) + designed minibinder (chain M).

    The Fab chains are copied verbatim from ``input_pdb``, so the final file is
    guaranteed to carry the user's placement of the arms bit-for-bit. Only the
    minibinder comes from RFd3.

    Two RFd3 output layouts are supported:

    * current - the designed minibinder is its own chain alongside the fixed Fab
      chains. RFd3 returns the Fab in the input frame, so ``fab_drift_rmsd_a``
      reports how far it actually moved and no superposition is applied unless
      that drift exceeds ``max_drift_a``.
    * legacy - the whole assembly is one fused polymer (VH|VL|MB|CH1|CL). The
      minibinder is superimposed via VL.

    Returns a dict with ``mb_length``, ``layout`` and ``fab_drift_rmsd_a``.
    """
    max_drift_a = 1.0
    src = _load_biopython_structure(input_pdb)
    in_model = list(src.get_models())[0]
    for cid in ("A", "B", "C", "D"):
        if cid not in in_model.child_dict:
            raise ValueError(f"Input Fab must contain chain {cid}")

    def std_residues(chain_id: str) -> list:
        return [res for res in in_model[chain_id].get_residues() if res.id[0] == " "]

    if ch1_len is None:
        ch1_len = len(std_residues("C"))
    if cl_len is None:
        cl_len = len(std_residues("D"))

    def ca_coords(residues: list) -> np.ndarray:
        return np.array([res["CA"].get_coord() for res in residues if "CA" in res])

    out_chains = _ordered_output_chains(rfd3_cif)
    in_chain_ids = [cid for cid in ("A", "B", "C", "D", "T") if cid in in_model.child_dict]

    if len(out_chains) > 1:
        layout = "separate_chain"
        fab_lengths = [len(std_residues(cid)) for cid in in_chain_ids]
        fab_centroids = [ca_coords(std_residues(cid)).mean(0) for cid in in_chain_ids]
        mb_res, out_fab = _split_minibinder_chain(out_chains, fab_lengths, fab_centroids)
        out_fab_ca = np.vstack([ca_coords(residues) for residues in out_fab])
        in_fab_ca = np.vstack([ca_coords(std_residues(cid)) for cid in in_chain_ids])
        drift = float(np.sqrt(((out_fab_ca - in_fab_ca) ** 2).sum(1).mean()))
        if drift <= max_drift_a:
            rot = np.eye(3)
            out_cent = in_cent = np.zeros(3)
        else:
            rot, out_cent, in_cent = _kabsch_transform(out_fab_ca, in_fab_ca)
    else:
        layout = "fused_polymer"
        out_residues, out_seq = _read_output_chain_residues(rfd3_cif)
        mb_start, mb_end, _mb_len = _infer_mb_segment(
            len(out_seq), vh_len, vl_len, ch1_len, cl_len
        )
        mb_res = out_residues[mb_start:mb_end]
        out_vl = (
            out_residues[:vl_len]
            if mb_start == vl_len
            else out_residues[vh_len : vh_len + vl_len]
        )
        rot, out_cent, in_cent = _kabsch_transform(ca_coords(out_vl), ca_coords(std_residues("B")))
        if mb_start == vl_len:
            segments = [("B", 0, vl_len), ("C", mb_end, mb_end + ch1_len)]
        else:
            segments = [
                ("A", 0, vh_len),
                ("B", vh_len, vh_len + vl_len),
                ("C", mb_end, mb_end + ch1_len),
                ("D", mb_end + ch1_len, mb_end + ch1_len + cl_len),
            ]
        # Drift after superposing on VL: how far the rest of the Fab was moved by
        # RFd3 relative to the input placement.
        sq_dev: list[np.ndarray] = []
        for cid, lo, hi in segments:
            seg = ca_coords(out_residues[lo:hi])
            ref = ca_coords(std_residues(cid))
            if len(seg) != len(ref):
                continue
            sq_dev.append((((seg - out_cent) @ rot + in_cent - ref) ** 2).sum(1))
        drift = float(np.sqrt(np.concatenate(sq_dev).mean())) if sq_dev else float("nan")

    struct = S.Structure("grafted")
    out_model = M.Model(0)
    for cid in ("A", "B", "C", "D"):
        chain = Ch.Chain(cid)
        for i, res in enumerate(std_residues(cid), start=1):
            _copy_residue(res, chain, i)
        out_model.add(chain)

    mb_chain = Ch.Chain(mb_chain_id)
    for i, res in enumerate(mb_res, start=1):
        new_res = Res.Residue((" ", i, " "), res.get_resname(), " ")
        for atom in res.get_atoms():
            xformed = (atom.get_coord() - out_cent) @ rot + in_cent
            new_res.add(
                At.Atom(
                    atom.name,
                    xformed,
                    atom.bfactor,
                    atom.occupancy,
                    atom.altloc,
                    atom.fullname,
                    atom.serial_number,
                    atom.element,
                )
            )
        mb_chain.add(new_res)
    out_model.add(mb_chain)

    if "T" in in_model.child_dict:
        t_chain = Ch.Chain("T")
        for i, res in enumerate(std_residues("T"), start=1):
            _copy_residue(res, t_chain, i)
        out_model.add(t_chain)

    struct.add(out_model)
    _write_structure(struct, out_path)
    return {
        "mb_length": len(mb_res),
        "layout": layout,
        "fab_drift_rmsd_a": round(drift, 3),
        "realigned": bool(drift > max_drift_a),
    }


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")
