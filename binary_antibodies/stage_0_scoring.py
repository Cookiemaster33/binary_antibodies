"""
stage_0_scoring.py
------------------
Score Stage 0 VH–VL interface weakening designs.

Core idea: select sequences with a large apo→holo pairing gap where holo still
looks like a real closed Fab engaging the epitope — not steric wedges.

Metrics (per design)
--------------------
  vh_vl_interface_contacts     — framework Cα pairs < 6 Å (want apo low, holo high)
  vh_vl_fr4_ca_distance      — FR4 centroid separation (want apo large, holo small)
  holo_minus_apo_contact_delta — holo − apo contacts (want positive)
  holo_fv_framework_rmsd_A     — Kabsch RMSD vs native Fab on VH+VL framework Cα
  cdr_epitope_contacts         — holo only: CDR (A/B) ↔ epitope T heavy-atom contacts
  holo_vh_vl_interface_clashes — holo only: cross-chain clashes at interface FW
  wedge_suspect                — apo very open but holo geometry/clashes fail (plug-like)

Run inside foundry Docker or any env with biotite + numpy.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from binary_antibodies.fab_hidden_switch import (  # noqa: E402
    VH_INTERFACE_FW,
    VL_INTERFACE_FW,
    cdr_residue_numbers,
    fv_framework_residue_numbers,
)

CONTACT_CUTOFF_A = 6.0
HEAVY_CONTACT_CUTOFF_A = 5.0
CLASH_CUTOFF_A = 4.0
INTERFACE_CLASH_CUTOFF_A = 2.8
INTER_CHAIN_CLASH_CUTOFF_A = 3.0
CDR_CONTACT_CUTOFF_A = 6.0

DEFAULT_FILTERS = {
    "max_apo_interface_contacts": 12,
    "min_holo_interface_contacts": 50,
    "max_holo_clashes_4A": 50,
    "max_holo_fv_framework_rmsd_A": 3.5,
    "min_holo_cdr_epitope_contacts": 6,
    "max_holo_vh_vl_interface_clashes": 0,
    "min_holo_minus_apo_contact_delta": 15,
    "min_apo_interface_centroid_distance_A": 14.0,
    "max_holo_interface_centroid_distance_A": 12.0,
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


def interface_contacts(vh_ca: dict[int, np.ndarray], vl_ca: dict[int, np.ndarray]) -> int:
    """Legacy Cα contact count (kept for logging)."""
    vh_keys = [r for r in VH_INTERFACE_FW if r in vh_ca]
    vl_keys = [r for r in VL_INTERFACE_FW if r in vl_ca]
    count = 0
    for rv in vh_keys:
        for rl in vl_keys:
            if np.linalg.norm(vh_ca[rv] - vl_ca[rl]) < CONTACT_CUTOFF_A:
                count += 1
    return count


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
    """Distance between VH/VL interface-framework Cα centroids (open vs closed)."""
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
    """Inter-chain heavy-atom overlaps; exclude epitope stub (rough placeholder geometry)."""
    skip = set(exclude_chains)
    heavy = aa[aa.element != "H"]
    n = len(heavy)
    clashes = 0
    chain = heavy.chain_id
    coords = heavy.coord
    for i in range(n):
        if chain[i] in skip:
            continue
        for j in range(i + 1, n):
            if chain[j] in skip:
                continue
            if chain[i] == chain[j]:
                continue
            if np.linalg.norm(coords[i] - coords[j]) < cutoff:
                clashes += 1
    return clashes


def _interface_residue_set(chain: str) -> set[int]:
    return set(VH_INTERFACE_FW if chain == "A" else VL_INTERFACE_FW)


def count_vh_vl_interface_clashes(aa, cutoff: float = INTERFACE_CLASH_CUTOFF_A) -> int:
    """Cross-chain A↔B clashes involving at least one interface-framework residue."""
    vh_iface = _interface_residue_set("A")
    vl_iface = _interface_residue_set("B")
    heavy = aa[aa.element != "H"]
    clashes = 0
    for i in range(len(heavy)):
        if heavy.chain_id[i] not in ("A", "B"):
            continue
        ri = int(heavy.res_id[i])
        for j in range(i + 1, len(heavy)):
            cj = heavy.chain_id[j]
            if (heavy.chain_id[i], cj) not in (("A", "B"), ("B", "A")):
                continue
            rj = int(heavy.res_id[j])
            iface_i = (
                (heavy.chain_id[i] == "A" and ri in vh_iface)
                or (heavy.chain_id[i] == "B" and ri in vl_iface)
            )
            iface_j = (
                (cj == "A" and rj in vh_iface)
                or (cj == "B" and rj in vl_iface)
            )
            if not (iface_i or iface_j):
                continue
            if np.linalg.norm(heavy.coord[i] - heavy.coord[j]) < cutoff:
                clashes += 1
    return clashes


def count_cdr_epitope_contacts(aa, cutoff: float = CDR_CONTACT_CUTOFF_A) -> int:
    """Heavy-atom contacts between CDR residues (A/B) and epitope stub chain T."""
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


class NativeFvReference:
    """Native VH+VL framework Cα from the Stage 0 design target PDB."""

    def __init__(self, path: Path) -> None:
        aa = load_structure(path)
        self.vh_ca = _collect_fv_framework_ca(aa, "A")
        self.vl_ca = _collect_fv_framework_ca(aa, "B")
        self.vh_iface_ca = chain_ca(aa, "A", VH_INTERFACE_FW)
        self.vl_iface_ca = chain_ca(aa, "B", VL_INTERFACE_FW)
        self.native_apo_contacts = interface_contacts_heavy(aa)

    def holo_fv_framework_rmsd(self, aa) -> float:
        """Align holo A+B framework to native, return interface-framework Cα RMSD."""
        vh = _collect_fv_framework_ca(aa, "A")
        vl = _collect_fv_framework_ca(aa, "B")
        common_vh = sorted(set(vh) & set(self.vh_ca))
        common_vl = sorted(set(vl) & set(self.vl_ca))
        if len(common_vh) < 8 or len(common_vl) < 8:
            return float("nan")

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

        iface_vh = sorted(set(vh) & set(self.vh_iface_ca))
        iface_vl = sorted(set(vl) & set(self.vl_iface_ca))
        if not iface_vh or not iface_vl:
            return float("nan")

        pts_p = np.vstack(
            [aligned[("A", r)] for r in iface_vh] + [aligned[("B", r)] for r in iface_vl]
        )
        pts_q = np.vstack(
            [self.vh_iface_ca[r] for r in iface_vh] + [self.vl_iface_ca[r] for r in iface_vl]
        )
        return float(np.sqrt(np.mean(np.sum((pts_p - pts_q) ** 2, axis=1))))


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
            metrics["holo_fv_framework_rmsd_A"] = holo_fv_framework_rmsd(aa, ref)
    return metrics


def wedge_suspect(apo: dict, holo: dict) -> bool:
    """Flag plug-like designs: very open apo but holo cannot close properly."""
    apo_open = apo.get("vh_vl_interface_contacts", 99) <= 5
    rmsd = holo.get("holo_fv_framework_rmsd_A", float("nan"))
    holo_bad_rmsd = not np.isnan(rmsd) and rmsd > 4.0
    holo_iface_clash = holo.get("vh_vl_interface_clashes", 0) > 0
    holo_low_cdr = holo.get("cdr_epitope_contacts", 0) < 3
    return apo_open and (holo_bad_rmsd or holo_iface_clash or holo_low_cdr)


def apply_filters(rec: dict, filters: dict) -> dict:
    apo = rec.get("apo", {})
    holo = rec.get("holo", {})
    checks = {
        "passes_apo_contacts": apo.get("vh_vl_interface_contacts", 999)
        <= filters["max_apo_interface_contacts"],
        "passes_holo_contacts": holo.get("vh_vl_interface_contacts", 0)
        >= filters["min_holo_interface_contacts"],
        "passes_holo_global_clashes": holo.get("inter_chain_clashes_4A", 999)
        <= filters["max_holo_clashes_4A"],
        "passes_contact_delta": rec.get("holo_minus_apo_contact_delta", -999)
        >= filters["min_holo_minus_apo_contact_delta"],
        "passes_apo_open": apo.get("vh_vl_interface_centroid_distance", 0)
        >= filters["min_apo_interface_centroid_distance_A"],
        "passes_holo_closed": holo.get("vh_vl_interface_centroid_distance", 999)
        <= filters["max_holo_interface_centroid_distance_A"],
        "passes_holo_geometry": holo.get("holo_fv_framework_rmsd_A", 999)
        <= filters["max_holo_fv_framework_rmsd_A"],
        "passes_cdr_engagement": holo.get("cdr_epitope_contacts", 0)
        >= filters["min_holo_cdr_epitope_contacts"],
        "passes_anti_wedge": holo.get("vh_vl_interface_clashes", 999)
        <= filters["max_holo_vh_vl_interface_clashes"]
        and not rec.get("wedge_suspect", True),
    }
    rec["filter_checks"] = checks
    rec["passes_stage_0"] = all(checks.values())
    return rec


def score_design(
    apo_cif: Path | None,
    holo_cif: Path | None,
    design_id: str,
    ref: NativeFvReference,
    filters: dict,
) -> dict:
    rec: dict = {"design_id": design_id}
    if apo_cif and apo_cif.exists():
        rec["apo"] = score_structure(apo_cif, has_epitope=False, ref=None)
    if holo_cif and holo_cif.exists():
        rec["holo"] = score_structure(holo_cif, has_epitope=True, ref=ref)
    if "apo" in rec and "holo" in rec:
        rec["holo_minus_apo_contact_delta"] = (
            rec["holo"]["vh_vl_interface_contacts"] - rec["apo"]["vh_vl_interface_contacts"]
        )
        rec["fr4_distance_delta"] = (
            rec["apo"]["vh_vl_fr4_ca_distance"] - rec["holo"]["vh_vl_fr4_ca_distance"]
        )
        rec["wedge_suspect"] = wedge_suspect(rec["apo"], rec["holo"])
        apply_filters(rec, filters)
    return rec


def rank_key(rec: dict) -> tuple:
    """Prefer large conditional gap with good holo geometry and CDR engagement."""
    checks = rec.get("filter_checks", {})
    n_pass = sum(1 for v in checks.values() if v)
    return (
        n_pass,
        rec.get("holo_minus_apo_contact_delta") or -999,
        rec.get("holo", {}).get("cdr_epitope_contacts", 0),
        -(rec.get("holo", {}).get("holo_fv_framework_rmsd_A", 999)),
        -(rec.get("apo", {}).get("vh_vl_interface_contacts", 999)),
    )


def main() -> None:
    p = argparse.ArgumentParser(description="Score Stage 0 VH–VL interface designs.")
    p.add_argument("--pipeline-dir", type=Path, default=Path("/workspace"))
    p.add_argument("--boltz-apo-subdir", default="boltz_outputs_apo")
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

    apo_base = pipeline / args.boltz_apo_subdir
    holo_base = pipeline / args.boltz_holo_subdir
    results: list[dict] = []

    for holo_pred in sorted(holo_base.glob("boltz_results_*/predictions/*/*_model_0.cif")):
        name = holo_pred.parent.name
        apo_pred = None
        for candidate in apo_base.glob(f"boltz_results_*/predictions/{name}/*_model_0.cif"):
            apo_pred = candidate
            break
        results.append(score_design(apo_pred, holo_pred, name, ref, filters))

    results.sort(key=rank_key, reverse=True)
    top = results[: args.top_n]
    out = {
        "n_scored": len(results),
        "n_passing": sum(1 for r in results if r.get("passes_stage_0")),
        "filters": filters,
        "reference_pdb": str(ref_path),
        "native_apo_interface_contacts": ref.native_apo_contacts,
        "designs": top,
    }
    out_path = pipeline / args.results_file
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Stage 0 scoring: {out['n_passing']}/{out['n_scored']} pass all filters → {out_path}")
    for r in top[:5]:
        apo_c = r.get("apo", {}).get("vh_vl_interface_contacts", "?")
        holo_c = r.get("holo", {}).get("vh_vl_interface_contacts", "?")
        delta = r.get("holo_minus_apo_contact_delta", "?")
        rmsd = r.get("holo", {}).get("holo_fv_framework_rmsd_A", "?")
        cdr_t = r.get("holo", {}).get("cdr_epitope_contacts", "?")
        wedge = r.get("wedge_suspect", "?")
        print(
            f"  {r['design_id']}: apo={apo_c} holo={holo_c} Δ={delta} "
            f"rmsd={rmsd} cdr-T={cdr_t} wedge={wedge} pass={r.get('passes_stage_0')}"
        )


if __name__ == "__main__":
    main()
