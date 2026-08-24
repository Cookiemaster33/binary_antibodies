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
# Stage A: spread VL/CL away from VH/CH1 so CH1 and VL hotspot surfaces are farther apart.
DEFAULT_STAGE_A_SEPARATION_A = 45.0

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


def _load_biopython_structure(path: Path):
    """Load PDB or mmCIF into a BioPython Structure."""
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


def _shift_residues_perpendicular(
    residues: list,
    ref_centroid: np.ndarray,
    other_centroid: np.ndarray,
    separation_a: float,
) -> list:
    """Translate residues perpendicular to the ref→other axis by separation_a."""
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


def fab_has_split_chains(source_pdb: Path) -> bool:
    """True when structure already uses logical Fab chains A–D (not fused H/L)."""
    src = _load_biopython_structure(source_pdb)
    model = list(src.get_models())[0]
    chains = set(model.child_dict)
    return {"A", "B"}.issubset(chains) and "H" not in chains and "L" not in chains


def infer_already_split(source_pdb: Path | None) -> bool:
    """True when Fab coordinates should be used as-is (no VL/CL re-translation)."""
    if source_pdb is None:
        return False
    name = source_pdb.name.lower()
    if name.endswith("_split.cif") or name.endswith("_split.pdb") or name.endswith("_split.mmcif"):
        return True
    return "_split" in source_pdb.stem and fab_has_split_chains(source_pdb)


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

    Writes a *_split.cif (or .pdb) suitable for PyMOL inspection and Stage A --fab-pdb.
    """
    vh_res, vl_res, ch1_res, cl_res, epitope_from_source = _load_fab_chain_residues(holo_cif)
    vl_res = _shift_residues_perpendicular(
        vl_res, _centroid(vh_res), _centroid(vl_res), separation_a
    )
    cl_res = _shift_residues_perpendicular(
        cl_res, _centroid(ch1_res), _centroid(cl_res), separation_a
    )
    struct = _assemble_fab_structure(
        vh_res, vl_res, ch1_res, cl_res, epitope_from_source, struct_name, place_stub_if_missing=True
    )
    _write_structure(struct, out_path)
    return len(vh_res), len(vl_res), len(ch1_res), len(cl_res)


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
) -> tuple[int, int, int, int]:
    """
    Build Stage A RFd3 target: full Fab (A–D + epitope T) with VL/CL translated apart.

    Pass a pre-built *_split.cif with already_split=True (or separation_a=0) to use
    coordinates as-is without re-translating VL/CL.
    """
    vh_res, vl_res, ch1_res, cl_res, epitope_from_source = _load_fab_chain_residues(source_pdb)

    if not already_split and separation_a > 0:
        vl_res = _shift_residues_perpendicular(
            vl_res, _centroid(vh_res), _centroid(vl_res), separation_a
        )
        cl_res = _shift_residues_perpendicular(
            cl_res, _centroid(ch1_res), _centroid(cl_res), separation_a
        )

    struct = _assemble_fab_structure(
        vh_res, vl_res, ch1_res, cl_res, epitope_from_source, struct_name
    )
    _write_structure(struct, out_pdb)
    return len(vh_res), len(vl_res), len(ch1_res), len(cl_res)


def stage_a_contig(
    vh_len: int,
    vl_len: int,
    ch1_len: int,
    cl_len: int,
    mb_length_range: str = "35-55",
) -> str:
    """RFd3 contig: VL and CH1 in output polymer with designed MB between them.

    VH, CL, and epitope T are provided as unindexed fixed context (see stage_a_unindex).
    """
    _ = vh_len, cl_len  # retained for call-site compatibility
    return f"B1-{vl_len}/0,{mb_length_range},C1-{ch1_len}"


def stage_a_unindex(
    vh_len: int,
    vl_len: int,
    cl_len: int,
    epitope_len: int = len(EPITOPE_SEQ),
) -> str:
    """Unindexed fixed Fab context: VH, CL, epitope stay in 3D space but not in output polymer."""
    _ = vl_len
    return f"A1-{vh_len},D1-{cl_len},T1-{epitope_len}"


def stage_a_rfd3_config(
    vh_len: int,
    vl_len: int,
    ch1_len: int,
    cl_len: int,
    mb_length_range: str = "35-55",
    epitope_len: int = len(EPITOPE_SEQ),
) -> dict[str, str | dict[str, str]]:
    """RFd3 inputs that keep the full Fab fixed in 3D while designing MB between VL and CH1."""
    return {
        "contig": stage_a_contig(vh_len, vl_len, ch1_len, cl_len, mb_length_range),
        "unindex": stage_a_unindex(vh_len, vl_len, cl_len, epitope_len),
        "select_fixed_atoms": stage_a_fixed_atoms(vh_len, vl_len, ch1_len, cl_len, epitope_len),
        "select_hotspots": stage_a_hotspots(vh_len),
        "mb_length_range": mb_length_range,
    }


def stage_a_fixed_atoms(
    vh_len: int,
    vl_len: int,
    ch1_len: int,
    cl_len: int,
    epitope_len: int = len(EPITOPE_SEQ),
) -> dict[str, str]:
    """All Fab + epitope chains fixed; only the unlinked minibinder is designed."""
    return {
        f"A1-{vh_len}": "ALL",
        f"B1-{vl_len}": "ALL",
        f"C1-{ch1_len}": "ALL",
        f"D1-{cl_len}": "ALL",
        f"T1-{epitope_len}": "ALL",
    }


def ch1_hotspots_chain_c(vh_len: int = VH_END) -> list[str]:
    return [f"C{r - vh_len}" for r in CH1_HOTSPOTS_HEAVY if r > vh_len]


def stage_a_hotspots(vh_len: int = VH_END) -> str:
    return ",".join(ch1_hotspots_chain_c(vh_len) + [f"B{r}" for r in VL_HOTSPOTS_STAGE_A])


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
    mb_min: int = 35,
    mb_max: int = 55,
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


def graft_stage_a_minibinder(
    out_path: Path,
    input_pdb: Path,
    rfd3_cif: Path,
    vh_len: int = VH_END,
    vl_len: int = VL_END,
    ch1_len: int | None = None,
    cl_len: int | None = None,
    mb_chain_id: str = "M",
) -> int:
    """
    Build output = exact input Fab (A,B,C,D,T) + designed minibinder (chain M).

    RFd3 only supplies minibinder coordinates; Fab chains are copied verbatim from
    input_pdb. The minibinder is superimposed via VL (chain B) alignment.
    Returns minibinder length.
    """
    src = _load_biopython_structure(input_pdb)
    in_model = list(src.get_models())[0]
    for cid in ("A", "B", "C", "D"):
        if cid not in in_model.child_dict:
            raise ValueError(f"Input Fab must contain chain {cid}")

    if ch1_len is None:
        ch1_len = sum(1 for r in in_model["C"] if r.id[0] == " ")
    if cl_len is None:
        cl_len = sum(1 for r in in_model["D"] if r.id[0] == " ")

    out_residues, out_seq = _read_output_chain_residues(rfd3_cif)
    mb_start, mb_end, mb_len = _infer_mb_segment(
        len(out_seq), vh_len, vl_len, ch1_len, cl_len
    )

    def ca_coords(residues: list) -> np.ndarray:
        return np.array([r["CA"].get_coord() for r in residues if "CA" in r])

    out_vl = out_residues[:vl_len] if mb_start == vl_len else out_residues[vh_len:vh_len + vl_len]
    in_vl = [r for r in in_model["B"].get_residues() if r.id[0] == " "]
    rot, out_cent, in_cent = _kabsch_transform(ca_coords(out_vl), ca_coords(in_vl))

    mb_res = out_residues[mb_start:mb_end]
    struct = S.Structure("grafted")
    out_model = M.Model(0)

    for cid in ("A", "B", "C", "D"):
        chain = Ch.Chain(cid)
        for i, res in enumerate([r for r in in_model[cid].get_residues() if r.id[0] == " "], start=1):
            _copy_residue(res, chain, i)
        out_model.add(chain)

    mb_chain = Ch.Chain(mb_chain_id)
    for i, res in enumerate(mb_res, start=1):
        new_res = Res.Residue((" ", i, " "), res.get_resname(), " ")
        for atom in res.get_atoms():
            coord = atom.get_coord()
            xformed = (coord - out_cent) @ rot + in_cent
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
        for i, res in enumerate([r for r in in_model["T"].get_residues() if r.id[0] == " "], start=1):
            _copy_residue(res, t_chain, i)
        out_model.add(t_chain)

    struct.add(out_model)
    _write_structure(struct, out_path)
    return mb_len


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")
