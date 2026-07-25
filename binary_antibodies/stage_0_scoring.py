"""
stage_0_scoring.py
------------------
Score Stage 0 VH–VL interface weakening designs.

Metrics (per design):
  - vh_vl_interface_contacts: Cβ/Cα pairs < 6 Å across framework interface residues
  - vh_vl_fr4_ca_distance: mean FR4 Cα separation (lower = more paired)
  - apo_minus_holo_contact_delta: apo_contacts - holo_contacts (want positive)

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

from binary_antibodies.fab_hidden_switch import VH_INTERFACE_FW, VL_INTERFACE_FW  # noqa: E402

CONTACT_CUTOFF_A = 6.0
CLASH_CUTOFF_A = 4.0


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


def chain_ca(aa, chain_id: str, resnums: list[int] | None = None) -> dict[int, np.ndarray]:
    mask = (aa.chain_id == chain_id) & (aa.atom_name == "CA")
    out: dict[int, np.ndarray] = {}
    for i in np.where(mask)[0]:
        r = int(aa.res_id[i])
        if resnums is None or r in resnums:
            out[r] = aa.coord[i]
    return out


def interface_contacts(vh_ca: dict[int, np.ndarray], vl_ca: dict[int, np.ndarray]) -> int:
    vh_keys = [r for r in VH_INTERFACE_FW if r in vh_ca]
    vl_keys = [r for r in VL_INTERFACE_FW if r in vl_ca]
    count = 0
    for rv in vh_keys:
        for rl in vl_keys:
            if np.linalg.norm(vh_ca[rv] - vl_ca[rl]) < CONTACT_CUTOFF_A:
                count += 1
    return count


def fr4_ca_distance(vh_ca: dict[int, np.ndarray], vl_ca: dict[int, np.ndarray]) -> float:
    vh_fr4 = [r for r in range(100, 114) if r in vh_ca]
    vl_fr4 = [r for r in range(100, 114) if r in vl_ca]
    if not vh_fr4 or not vl_fr4:
        return float("nan")
    vh_cent = np.mean([vh_ca[r] for r in vh_fr4], axis=0)
    vl_cent = np.mean([vl_ca[r] for r in vl_fr4], axis=0)
    return float(np.linalg.norm(vh_cent - vl_cent))


def count_clashes(aa, cutoff: float = CLASH_CUTOFF_A) -> int:
    heavy = aa[aa.element != "H"]
    clashes = 0
    for i in range(len(heavy)):
        d = np.linalg.norm(heavy.coord[i + 1 :] - heavy.coord[i], axis=1)
        clashes += int(np.sum(d < cutoff))
    return clashes


def score_structure(path: Path) -> dict:
    aa = load_structure(path)
    vh_ca = chain_ca(aa, "A")
    vl_ca = chain_ca(aa, "B")
    return {
        "structure": path.name,
        "vh_vl_interface_contacts": interface_contacts(vh_ca, vl_ca),
        "vh_vl_fr4_ca_distance": fr4_ca_distance(vh_ca, vl_ca),
        "clashes_4A": count_clashes(aa),
    }


def score_design(
    apo_cif: Path | None,
    holo_cif: Path | None,
    design_id: str,
) -> dict:
    rec: dict = {"design_id": design_id}
    if apo_cif and apo_cif.exists():
        apo = score_structure(apo_cif)
        rec["apo"] = apo
    if holo_cif and holo_cif.exists():
        holo = score_structure(holo_cif)
        rec["holo"] = holo
    if "apo" in rec and "holo" in rec:
        rec["apo_minus_holo_contact_delta"] = (
            rec["apo"]["vh_vl_interface_contacts"] - rec["holo"]["vh_vl_interface_contacts"]
        )
        rec["passes_apo_filter"] = rec["apo"]["vh_vl_interface_contacts"] <= 12
        rec["passes_holo_filter"] = (
            rec["holo"]["vh_vl_interface_contacts"] >= 20 and rec["holo"]["clashes_4A"] == 0
        )
        rec["passes_stage_0"] = rec["passes_apo_filter"] and rec["passes_holo_filter"]
    return rec


def main() -> None:
    p = argparse.ArgumentParser(description="Score Stage 0 VH–VL interface designs.")
    p.add_argument("--pipeline-dir", type=Path, default=Path("/workspace"))
    p.add_argument("--boltz-apo-subdir", default="boltz_outputs_apo")
    p.add_argument("--boltz-holo-subdir", default="boltz_outputs_holo")
    p.add_argument("--results-file", default="final/stage_0_results.json")
    p.add_argument("--top-n", type=int, default=50)
    args = p.parse_args()

    pipeline = args.pipeline_dir
    apo_base = pipeline / args.boltz_apo_subdir
    holo_base = pipeline / args.boltz_holo_subdir
    results: list[dict] = []

    for holo_pred in sorted(holo_base.glob("boltz_results_*/predictions/*/*_model_0.cif")):
        name = holo_pred.parent.name
        apo_pred = None
        for candidate in apo_base.glob(f"boltz_results_*/predictions/{name}/*_model_0.cif"):
            apo_pred = candidate
            break
        results.append(score_design(apo_pred, holo_pred, name))

    results.sort(
        key=lambda r: (
            -(r.get("apo_minus_holo_contact_delta") or -999),
            r.get("apo", {}).get("vh_vl_interface_contacts", 999),
        )
    )
    top = results[: args.top_n]
    out = {
        "n_scored": len(results),
        "n_passing": sum(1 for r in results if r.get("passes_stage_0")),
        "designs": top,
    }
    out_path = pipeline / args.results_file
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Stage 0 scoring: {out['n_passing']}/{out['n_scored']} pass filters → {out_path}")
    for r in top[:5]:
        apo_c = r.get("apo", {}).get("vh_vl_interface_contacts", "?")
        holo_c = r.get("holo", {}).get("vh_vl_interface_contacts", "?")
        delta = r.get("apo_minus_holo_contact_delta", "?")
        print(f"  {r['design_id']}: apo={apo_c} holo={holo_c} Δ={delta} pass={r.get('passes_stage_0')}")


if __name__ == "__main__":
    main()
