"""
stage_0_scoring.py
------------------
Score Stage 0 VH–VL interface weakening designs (split MPNN de-grease).

Selection criteria
------------------
(a) Static VH–VL interface contacts on the native Fab backbone, weighted by the
    designed sequence, should be much lower than WT (~113 heavy-atom contacts).
(b) Holo-only Boltz (A+B+C+D+T) must preserve epitope engagement and Fab-like
    geometry — the weakened interface must not block target binding.

Apo Boltz folds are not required.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Inlined from fab_hidden_switch (no BioPython — runs in foundry Docker)
VH_INTERFACE_FW = [34, 38, 42, 43, 44, 45, 46, 47, 49, 87, 89]
VL_INTERFACE_FW = [35, 37, 39, 43, 44, 45, 46, 47, 104, 105, 106, 107, 108, 109, 110, 111, 112]
VH_CDR_RANGES = [(26, 35), (50, 65), (95, 102)]
VL_CDR_RANGES = [(24, 34), (50, 56), (89, 97)]

HYDROPHOBIC = frozenset("AILMFWV")
POSITIVE = frozenset("KRH")
NEGATIVE = frozenset("DE")
POLAR = frozenset("STNQ")


def in_ranges(resnum: int, ranges: list[tuple[int, int]]) -> bool:
    return any(lo <= resnum <= hi for lo, hi in ranges)


def cdr_residue_numbers(chain: str) -> list[int]:
    if chain == "A":
        return [r for r in range(1, 114) if in_ranges(r, VH_CDR_RANGES)]
    return [r for r in range(1, 108) if in_ranges(r, VL_CDR_RANGES)]


def fv_framework_residue_numbers(chain: str) -> list[int]:
    ranges = VH_CDR_RANGES if chain == "A" else VL_CDR_RANGES
    length = 113 if chain == "A" else 107
    return [r for r in range(1, length + 1) if not in_ranges(r, ranges)]


CONTACT_CUTOFF_A = 6.0
HEAVY_CONTACT_CUTOFF_A = 5.0
CLASH_CUTOFF_A = 4.0
INTERFACE_CLASH_CUTOFF_A = 2.8
INTER_CHAIN_CLASH_CUTOFF_A = 3.0
CDR_CONTACT_CUTOFF_A = 6.0

DEFAULT_FILTERS = {
    "max_static_fraction_of_native_contacts": 0.45,
    "max_holo_clashes_4A": 50,
    "max_holo_fv_framework_rmsd_A": 3.5,
    "max_holo_fv_cdr_rmsd_A": 6.0,
    "min_holo_cdr_epitope_contacts": 6,
    "max_holo_vh_vl_interface_clashes": 0,
}


def load_structure(path: Path):
    import biotite.structure.io.pdbx as pdbx
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if path.suffix.lower() in {".cif", ".mmcif"}:
            cif = pdbx.CIFFile.read(str(path))
            return pdbx.get_structure(cif, model=1, include_bonds=False)
        from biotite.structure.io.pdb import PDBFile

        pdb = PDBFile.read(str(path))
        return pdb.get_structure(model=1)


def load_filters(config_path: Path | None) -> dict:
    filters = dict(DEFAULT_FILTERS)
    if config_path and config_path.exists():
        cfg = json.loads(config_path.read_text())
        filters.update(cfg.get("validation", {}).get("filters", {}))
    return filters


def chain_ca(aa, chain_id: str, resnums: list[int] | None = None) -> dict[int, np.ndarray]:
    mask = (aa.chain_id == chain_id) & (aa.atom_name == "CA")
    out: dict[int, np.ndarray] = {}
    for i in np.where(mask)[0]:
        r = int(aa.res_id[i])
        if resnums is None or r in resnums:
            out[r] = aa.coord[i]
    return out


def _interface_residue_set(chain: str) -> set[int]:
    return set(VH_INTERFACE_FW if chain == "A" else VL_INTERFACE_FW)


def interface_contacts_heavy(aa) -> int:
    """Heavy-atom contacts between VH/VL framework interface residues (< 5 Å)."""
    vh_iface = _interface_residue_set("A")
    vl_iface = _interface_residue_set("B")
    vh_atoms = aa[(aa.chain_id == "A") & np.isin(aa.res_id, list(vh_iface)) & (aa.element != "H")]
    vl_atoms = aa[(aa.chain_id == "B") & np.isin(aa.res_id, list(vl_iface)) & (aa.element != "H")]
    if len(vh_atoms) == 0 or len(vl_atoms) == 0:
        return 0
    count = 0
    for coord in vh_atoms.coord:
        d = np.linalg.norm(vl_atoms.coord - coord, axis=1)
        count += int(np.sum(d < HEAVY_CONTACT_CUTOFF_A))
    return count


def interface_centroid_distance(vh_ca: dict[int, np.ndarray], vl_ca: dict[int, np.ndarray]) -> float:
    vh_pts = [vh_ca[r] for r in VH_INTERFACE_FW if r in vh_ca]
    vl_pts = [vl_ca[r] for r in VL_INTERFACE_FW if r in vl_ca]
    if not vh_pts or not vl_pts:
        return float("nan")
    return float(np.linalg.norm(np.mean(vh_pts, axis=0) - np.mean(vl_pts, axis=0)))


def fr4_ca_distance(vh_ca: dict[int, np.ndarray], vl_ca: dict[int, np.ndarray]) -> float:
    vh_fr4 = [r for r in range(100, 114) if r in vh_ca]
    vl_fr4 = [r for r in range(100, 114) if r in vl_ca]
    if not vh_fr4 or not vl_fr4:
        return float("nan")
    vh_cent = np.mean([vh_ca[r] for r in vh_fr4], axis=0)
    vl_cent = np.mean([vl_ca[r] for r in vl_fr4], axis=0)
    return float(np.linalg.norm(vh_cent - vl_cent))


def count_inter_chain_clashes(
    aa,
    cutoff: float = INTER_CHAIN_CLASH_CUTOFF_A,
    exclude_chains: tuple[str, ...] = ("T",),
) -> int:
    from scipy.spatial import cKDTree

    skip = set(exclude_chains)
    heavy = aa[aa.element != "H"]
    chains = heavy.chain_id
    coords = heavy.coord
    unique = [c for c in sorted(set(chains)) if c not in skip]
    clashes = 0
    for i, c1 in enumerate(unique):
        m1 = chains == c1
        pts1 = coords[m1]
        if len(pts1) == 0:
            continue
        for c2 in unique[i + 1 :]:
            m2 = chains == c2
            pts2 = coords[m2]
            if len(pts2) == 0:
                continue
            tree = cKDTree(pts2)
            for neighbors in tree.query_ball_point(pts1, cutoff):
                clashes += len(neighbors)
    return clashes


def count_vh_vl_interface_clashes(aa, cutoff: float = INTERFACE_CLASH_CUTOFF_A) -> int:
    from scipy.spatial import cKDTree

    heavy = aa[aa.element != "H"]
    vh = heavy[(heavy.chain_id == "A") & np.isin(heavy.res_id, list(VH_INTERFACE_FW))]
    vl = heavy[(heavy.chain_id == "B") & np.isin(heavy.res_id, list(VL_INTERFACE_FW))]
    if len(vh) == 0 or len(vl) == 0:
        return 0
    tree = cKDTree(vl.coord)
    clashes = 0
    for pt in vh.coord:
        clashes += len(tree.query_ball_point(pt, cutoff))
    return clashes


def count_cdr_epitope_contacts(aa, cutoff: float = CDR_CONTACT_CUTOFF_A) -> int:
    cdr_nums_a = set(cdr_residue_numbers("A"))
    cdr_nums_b = set(cdr_residue_numbers("B"))
    cdr_mask = (
        ((aa.chain_id == "A") & np.isin(aa.res_id, list(cdr_nums_a)))
        | ((aa.chain_id == "B") & np.isin(aa.res_id, list(cdr_nums_b)))
    )
    t_mask = aa.chain_id == "T"
    cdr_atoms = aa[cdr_mask & (aa.element != "H")]
    t_atoms = aa[t_mask & (aa.element != "H")]
    if len(cdr_atoms) == 0 or len(t_atoms) == 0:
        return 0
    count = 0
    for c in cdr_atoms.coord:
        d = np.linalg.norm(t_atoms.coord - c, axis=1)
        count += int(np.sum(d < cutoff))
    return count


def _collect_fv_framework_ca(aa, chain_id: str) -> dict[int, np.ndarray]:
    fw = fv_framework_residue_numbers(chain_id)
    return chain_ca(aa, chain_id, fw)


def sequence_residue(chain_seq: str, resnum: int) -> str:
    idx = resnum - 1
    if 0 <= idx < len(chain_seq):
        return chain_seq[idx]
    return "X"


def pair_contact_weight(aa_vh: str, aa_vl: str) -> float:
    """Compatibility weight for a VH/VL residue pair maintaining a native geometry contact."""
    if aa_vh == "X" or aa_vl == "X":
        return 0.5
    if aa_vh in POSITIVE and aa_vl in POSITIVE:
        return 0.0
    if aa_vh in NEGATIVE and aa_vl in NEGATIVE:
        return 0.0
    if (aa_vh in POSITIVE and aa_vl in NEGATIVE) or (aa_vh in NEGATIVE and aa_vl in POSITIVE):
        return 0.55
    if "G" in (aa_vh, aa_vl):
        return 0.15
    if "P" in (aa_vh, aa_vl):
        return 0.25
    if "A" in (aa_vh, aa_vl):
        return 0.35
    if aa_vh in HYDROPHOBIC and aa_vl in HYDROPHOBIC:
        return 1.0
    if aa_vh in POLAR and aa_vl in POLAR:
        return 0.6
    if (aa_vh in HYDROPHOBIC and aa_vl in POLAR) or (aa_vh in POLAR and aa_vl in HYDROPHOBIC):
        return 0.45
    return 0.5


def contact_pair_weight(
    seq_a: str,
    seq_b: str,
    wt_a: str,
    wt_b: str,
    vh_res: int,
    vl_res: int,
) -> float:
    """
    Weight for one native heavy-atom contact pair under a designed sequence.
    Unchanged WT pairs count as 1.0; mutations scale by designed compatibility.
    """
    des_vh = sequence_residue(seq_a, vh_res)
    des_vl = sequence_residue(seq_b, vl_res)
    wt_vh = sequence_residue(wt_a, vh_res)
    wt_vl = sequence_residue(wt_b, vl_res)
    if des_vh == wt_vh and des_vl == wt_vl:
        return 1.0
    return pair_contact_weight(des_vh, des_vl)


def static_interface_contacts_heavy(
    aa_template,
    seq_a: str,
    seq_b: str,
    wt_a: str | None = None,
    wt_b: str | None = None,
) -> float:
    """
    Predict VH–VL interface contacts on the native Fab geometry, weighting each
    native heavy-atom pair by designed-sequence compatibility at those residues.
    """
    if wt_a is None or wt_b is None:
        wt_a = seq_a
        wt_b = seq_b
    vh_iface = _interface_residue_set("A")
    vl_iface = _interface_residue_set("B")
    vh_atoms = aa_template[
        (aa_template.chain_id == "A")
        & np.isin(aa_template.res_id, list(vh_iface))
        & (aa_template.element != "H")
    ]
    vl_atoms = aa_template[
        (aa_template.chain_id == "B")
        & np.isin(aa_template.res_id, list(vl_iface))
        & (aa_template.element != "H")
    ]
    if len(vh_atoms) == 0 or len(vl_atoms) == 0:
        return 0.0

    total = 0.0
    for vh_idx in range(len(vh_atoms)):
        vh_res = int(vh_atoms.res_id[vh_idx])
        aa_vh = sequence_residue(seq_a, vh_res)
        coord = vh_atoms.coord[vh_idx]
        d = np.linalg.norm(vl_atoms.coord - coord, axis=1)
        close = np.where(d < HEAVY_CONTACT_CUTOFF_A)[0]
        for j in close:
            vl_res = int(vl_atoms.res_id[j])
            total += contact_pair_weight(seq_a, seq_b, wt_a, wt_b, vh_res, vl_res)
    return total


def _collect_fv_cdr_ca(aa, chain_id: str) -> dict[int, np.ndarray]:
    return chain_ca(aa, chain_id, cdr_residue_numbers(chain_id))


def _ca_rmsd(
    aligned: dict[tuple[str, int], np.ndarray],
    native_ca: dict[int, np.ndarray],
    holo_ca: dict[int, np.ndarray],
    resnums: list[int],
    chain: str,
) -> float:
    common = sorted(set(holo_ca) & set(native_ca) & set(resnums))
    if len(common) < 3:
        return float("nan")
    pts_p = np.vstack([aligned[(chain, r)] for r in common])
    pts_q = np.vstack([native_ca[r] for r in common])
    return float(np.sqrt(np.mean(np.sum((pts_p - pts_q) ** 2, axis=1))))


class NativeFvReference:
    """Native VH+VL framework from the Stage 0 design target PDB."""

    def __init__(self, path: Path) -> None:
        self.path = path
        aa = load_structure(path)
        self.aa = aa
        self.vh_ca = _collect_fv_framework_ca(aa, "A")
        self.vl_ca = _collect_fv_framework_ca(aa, "B")
        self.vh_cdr_ca = _collect_fv_cdr_ca(aa, "A")
        self.vl_cdr_ca = _collect_fv_cdr_ca(aa, "B")
        self.vh_iface_ca = chain_ca(aa, "A", VH_INTERFACE_FW)
        self.vl_iface_ca = chain_ca(aa, "B", VL_INTERFACE_FW)
        self.vh_fw_nums = fv_framework_residue_numbers("A")
        self.vl_fw_nums = fv_framework_residue_numbers("B")
        self.vh_cdr_nums = cdr_residue_numbers("A")
        self.vl_cdr_nums = cdr_residue_numbers("B")
        self.native_wt_contacts = interface_contacts_heavy(aa)
        self.native_seq_a = _chain_sequence(aa, "A")
        self.native_seq_b = _chain_sequence(aa, "B")
        self.native_static_contacts = static_interface_contacts_heavy(
            aa, self.native_seq_a, self.native_seq_b, self.native_seq_a, self.native_seq_b
        )

    def holo_fv_rmsd_metrics(self, aa) -> dict[str, float]:
        """
        Align holo VH+VL to native on all framework Cα (CDRs excluded), then report:
          - holo_fv_framework_rmsd_A: full Fv framework
          - holo_fv_interface_framework_rmsd_A: VH–VL interface rim only
          - holo_fv_cdr_rmsd_A: CDR loops (expected to move with epitope stub T)
        """
        vh = _collect_fv_framework_ca(aa, "A")
        vl = _collect_fv_framework_ca(aa, "B")
        vh_cdr = _collect_fv_cdr_ca(aa, "A")
        vl_cdr = _collect_fv_cdr_ca(aa, "B")
        common_vh = sorted(set(vh) & set(self.vh_ca))
        common_vl = sorted(set(vl) & set(self.vl_ca))
        if len(common_vh) < 8 or len(common_vl) < 8:
            nan = float("nan")
            return {
                "holo_fv_framework_rmsd_A": nan,
                "holo_fv_interface_framework_rmsd_A": nan,
                "holo_fv_cdr_rmsd_A": nan,
            }

        P = np.vstack([vh[r] for r in common_vh] + [vl[r] for r in common_vl])
        Q = np.vstack([self.vh_ca[r] for r in common_vh] + [self.vl_ca[r] for r in common_vl])
        Pc = P - P.mean(0)
        Qc = Q - Q.mean(0)
        U, _S, Vt = np.linalg.svd(Pc.T @ Qc)
        d = np.linalg.det(Vt.T @ U.T)
        R = Vt.T @ np.diag([1, 1, d]) @ U.T
        P_aligned = Pc @ R.T + Q.mean(0)

        aligned: dict[tuple[str, int], np.ndarray] = {}
        idx = 0
        for r in common_vh:
            aligned[("A", r)] = P_aligned[idx]
            idx += 1
        for r in common_vl:
            aligned[("B", r)] = P_aligned[idx]
            idx += 1

        # CDR coordinates after the same rigid-body transform (not used for alignment).
        for r in sorted(set(vh_cdr) & set(self.vh_cdr_ca)):
            aligned[("A", r)] = (vh_cdr[r] - P.mean(0)) @ R.T + Q.mean(0)
        for r in sorted(set(vl_cdr) & set(self.vl_cdr_ca)):
            aligned[("B", r)] = (vl_cdr[r] - P.mean(0)) @ R.T + Q.mean(0)

        fw_rmsd = _ca_rmsd(aligned, self.vh_ca, vh, self.vh_fw_nums, "A")
        vl_fw = _ca_rmsd(aligned, self.vl_ca, vl, self.vl_fw_nums, "B")
        if not np.isnan(fw_rmsd) and not np.isnan(vl_fw):
            fw_rmsd = float(np.sqrt((fw_rmsd**2 + vl_fw**2) / 2))
        elif np.isnan(fw_rmsd):
            fw_rmsd = vl_fw

        iface_vh = sorted(set(vh) & set(self.vh_iface_ca))
        iface_vl = sorted(set(vl) & set(self.vl_iface_ca))
        if iface_vh and iface_vl:
            pts_p = np.vstack(
                [aligned[("A", r)] for r in iface_vh] + [aligned[("B", r)] for r in iface_vl]
            )
            pts_q = np.vstack(
                [self.vh_iface_ca[r] for r in iface_vh] + [self.vl_iface_ca[r] for r in iface_vl]
            )
            iface_rmsd = float(np.sqrt(np.mean(np.sum((pts_p - pts_q) ** 2, axis=1))))
        else:
            iface_rmsd = float("nan")

        vh_cdr_r = _ca_rmsd(aligned, self.vh_cdr_ca, vh_cdr, self.vh_cdr_nums, "A")
        vl_cdr_r = _ca_rmsd(aligned, self.vl_cdr_ca, vl_cdr, self.vl_cdr_nums, "B")
        if not np.isnan(vh_cdr_r) and not np.isnan(vl_cdr_r):
            cdr_rmsd = float(np.sqrt((vh_cdr_r**2 + vl_cdr_r**2) / 2))
        elif not np.isnan(vh_cdr_r):
            cdr_rmsd = vh_cdr_r
        else:
            cdr_rmsd = vl_cdr_r

        return {
            "holo_fv_framework_rmsd_A": fw_rmsd,
            "holo_fv_interface_framework_rmsd_A": iface_rmsd,
            "holo_fv_cdr_rmsd_A": cdr_rmsd,
        }

    def holo_fv_framework_rmsd(self, aa) -> float:
        """Backward-compatible alias: full Fv framework RMSD."""
        return self.holo_fv_rmsd_metrics(aa)["holo_fv_framework_rmsd_A"]


def _chain_sequence(aa, chain_id: str) -> str:
    three_to_one = {
        "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
        "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
        "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    }
    mask = aa.chain_id == chain_id
    residues = sorted({int(r) for r in aa.res_id[mask]})
    chars: list[str] = []
    for resnum in residues:
        res_mask = mask & (aa.res_id == resnum)
        names = aa.res_name[res_mask]
        res_name = str(names[0]).upper()
        chars.append(three_to_one.get(res_name, "X"))
    return "".join(chars)


def holo_fv_framework_rmsd(aa, ref: NativeFvReference) -> float:
    return ref.holo_fv_framework_rmsd(aa)


def score_structure(
    path: Path,
    *,
    has_epitope: bool,
    ref: NativeFvReference | None = None,
) -> dict:
    aa = load_structure(path)
    vh_ca = chain_ca(aa, "A")
    vl_ca = chain_ca(aa, "B")
    metrics: dict = {
        "structure": path.name,
        "vh_vl_interface_contacts": interface_contacts_heavy(aa),
        "vh_vl_interface_centroid_distance": interface_centroid_distance(
            chain_ca(aa, "A", VH_INTERFACE_FW),
            chain_ca(aa, "B", VL_INTERFACE_FW),
        ),
        "vh_vl_fr4_ca_distance": fr4_ca_distance(vh_ca, vl_ca),
        "inter_chain_clashes_4A": count_inter_chain_clashes(aa),
    }
    if has_epitope:
        metrics["cdr_epitope_contacts"] = count_cdr_epitope_contacts(aa)
        metrics["vh_vl_interface_clashes"] = count_vh_vl_interface_clashes(aa)
        if ref is not None:
            metrics.update(ref.holo_fv_rmsd_metrics(aa))
    return metrics


def apply_filters(rec: dict, filters: dict, native_contacts: int) -> dict:
    holo = rec.get("holo", {})
    static_c = rec.get("static_vh_vl_interface_contacts", 999)
    static_frac = static_c / native_contacts if native_contacts > 0 else 1.0

    checks: dict[str, bool] = {}

    if "max_static_fraction_of_native_contacts" in filters:
        checks["passes_interface_weakened"] = (
            static_frac <= filters["max_static_fraction_of_native_contacts"]
        )
    elif "max_static_interface_contacts" in filters:
        checks["passes_interface_weakened"] = static_c <= filters["max_static_interface_contacts"]

    checks["passes_holo_global_clashes"] = holo.get("inter_chain_clashes_4A", 999) <= filters[
        "max_holo_clashes_4A"
    ]
    checks["passes_holo_framework_geometry"] = holo.get("holo_fv_framework_rmsd_A", 999) <= filters[
        "max_holo_fv_framework_rmsd_A"
    ]
    if "max_holo_fv_cdr_rmsd_A" in filters:
        checks["passes_holo_cdr_geometry"] = holo.get("holo_fv_cdr_rmsd_A", 999) <= filters[
            "max_holo_fv_cdr_rmsd_A"
        ]
    checks["passes_cdr_engagement"] = holo.get("cdr_epitope_contacts", 0) >= filters[
        "min_holo_cdr_epitope_contacts"
    ]
    checks["passes_holo_interface_clean"] = holo.get("vh_vl_interface_clashes", 999) <= filters[
        "max_holo_vh_vl_interface_clashes"
    ]

    rec["static_fraction_of_native_contacts"] = round(static_frac, 4)
    rec["filter_checks"] = checks
    rec["passes_stage_0"] = all(checks.values())
    return rec


def score_design(
    holo_cif: Path | None,
    design_id: str,
    ref: NativeFvReference,
    filters: dict,
    *,
    seq_a: str | None = None,
    seq_b: str | None = None,
) -> dict:
    rec: dict = {"design_id": design_id}
    if seq_a and seq_b:
        static_c = static_interface_contacts_heavy(
            ref.aa, seq_a, seq_b, ref.native_seq_a, ref.native_seq_b
        )
        rec["static_vh_vl_interface_contacts"] = round(static_c, 2)
        rec["chains"] = {"A": seq_a, "B": seq_b}
    if holo_cif and holo_cif.exists():
        rec["holo"] = score_structure(holo_cif, has_epitope=True, ref=ref)
    if "static_vh_vl_interface_contacts" in rec and "holo" in rec:
        apply_filters(rec, filters, ref.native_wt_contacts)
    return rec


def rank_key(rec: dict) -> tuple:
    """Prefer filter passes, then weaker static interface, then holo binding quality."""
    checks = rec.get("filter_checks", {})
    n_pass = sum(1 for v in checks.values() if v)
    static_c = rec.get("static_vh_vl_interface_contacts", 999)
    return (
        n_pass,
        -static_c,
        rec.get("holo", {}).get("cdr_epitope_contacts", 0),
        -(rec.get("holo", {}).get("holo_fv_framework_rmsd_A", 999)),
        -(rec.get("holo", {}).get("holo_fv_cdr_rmsd_A", 999)),
    )


def load_sequence_lookup(pipeline: Path) -> dict[str, tuple[str, str]]:
    """Map design_id → (seq_a, seq_b) from MPNN output."""
    mpnn_path = pipeline / "outputs" / "mpnn_stage_0" / "all_sequences.json"
    top_path = pipeline / "final" / "top_designs_stage_0.json"
    lookup: dict[str, tuple[str, str]] = {}

    if top_path.exists():
        top = json.loads(top_path.read_text())
        for d in top.get("designs", []):
            chains = d.get("chains", {})
            if chains.get("A") and chains.get("B"):
                lookup[d["name"]] = (chains["A"], chains["B"])

    if mpnn_path.exists() and not lookup:
        seqs = json.loads(mpnn_path.read_text())
        if isinstance(seqs, dict):
            seqs = seqs.get("sequences", [])
        for d in seqs:
            chains = d.get("chains", {})
            if chains.get("A") and chains.get("B"):
                key = f"s{d['seq_idx']}"
                lookup[key] = (chains["A"], chains["B"])
    return lookup


def design_id_to_sequences(design_id: str, lookup: dict[str, tuple[str, str]]) -> tuple[str, str] | tuple[None, None]:
    if design_id in lookup:
        return lookup[design_id]
    if "_s" in design_id:
        suffix = design_id.rsplit("_s", 1)[-1]
        if suffix.isdigit():
            key = f"s{suffix}"
            for k, val in lookup.items():
                if k.endswith(suffix) or k == key:
                    return val
    return None, None


def main() -> None:
    p = argparse.ArgumentParser(description="Score Stage 0 VH–VL interface designs (holo-only).")
    p.add_argument("--pipeline-dir", type=Path, default=Path("/workspace"))
    p.add_argument("--boltz-holo-subdir", default="boltz_outputs_holo")
    p.add_argument("--reference-pdb", type=Path, default=None)
    p.add_argument("--config-json", type=Path, default=None)
    p.add_argument("--results-file", default="final/stage_0_results.json")
    p.add_argument("--top-n", type=int, default=50)
    args = p.parse_args()

    pipeline = args.pipeline_dir
    ref_path = args.reference_pdb or (ROOT / "structures/domains/fab_stage_0_vhvL_interface.pdb")
    config_path = args.config_json or (ROOT / "structures/interface/stage_0_vhvL_interface_config.json")
    filters = load_filters(config_path)
    ref = NativeFvReference(ref_path)
    seq_lookup = load_sequence_lookup(pipeline)

    holo_base = pipeline / args.boltz_holo_subdir
    results: list[dict] = []

    for holo_pred in sorted(holo_base.glob("boltz_results_*/predictions/*/*_model_0.cif")):
        name = holo_pred.parent.name
        seq_a, seq_b = design_id_to_sequences(name, seq_lookup)
        results.append(score_design(holo_pred, name, ref, filters, seq_a=seq_a, seq_b=seq_b))

    results.sort(key=rank_key, reverse=True)
    top = results[: args.top_n]
    out = {
        "scoring_mode": "holo_only_static_interface",
        "n_scored": len(results),
        "n_passing": sum(1 for r in results if r.get("passes_stage_0")),
        "filters": filters,
        "reference_pdb": str(ref_path),
        "native_wt_interface_contacts": ref.native_wt_contacts,
        "native_static_interface_contacts": round(ref.native_static_contacts, 2),
        "designs": top,
    }
    out_path = pipeline / args.results_file
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Stage 0 scoring: {out['n_passing']}/{out['n_scored']} pass all filters → {out_path}")
    for r in top[:5]:
        static_c = r.get("static_vh_vl_interface_contacts", "?")
        static_f = r.get("static_fraction_of_native_contacts", "?")
        fw_rmsd = r.get("holo", {}).get("holo_fv_framework_rmsd_A", "?")
        iface_rmsd = r.get("holo", {}).get("holo_fv_interface_framework_rmsd_A", "?")
        cdr_rmsd = r.get("holo", {}).get("holo_fv_cdr_rmsd_A", "?")
        cdr_t = r.get("holo", {}).get("cdr_epitope_contacts", "?")
        holo_c = r.get("holo", {}).get("vh_vl_interface_contacts", "?")
        print(
            f"  {r['design_id']}: static={static_c} ({static_f}×WT) holo_iface={holo_c} "
            f"fw_rmsd={fw_rmsd} iface_rmsd={iface_rmsd} cdr_rmsd={cdr_rmsd} cdr-T={cdr_t} "
            f"pass={r.get('passes_stage_0')}"
        )


if __name__ == "__main__":
    main()
