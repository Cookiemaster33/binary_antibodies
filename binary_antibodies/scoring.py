"""
scoring.py
----------
Boltz-2 output parsing and global assembly RMSD scoring for the bispecific
bridging minibinder pipeline. Runs inside the foundry Docker container
(has atomworks + biotite).

Validation RMSD: Kabsch-aligned Cα RMSD over the full VH1 + minibinder +
nanobody assembly (chain A). This checks that Boltz-2 preserves the designed
bispecific geometry, not just local minibinder backbone shape.

Key facts about Boltz-2 output format:
  - Predictions are at: <out_dir>/boltz_results_<input_dir_name>/predictions/<name>/
  - Confidence JSON: predictions/<name>/confidence_<name>_model_0.json
  - Structure CIF:   predictions/<name>/<name>_model_0.cif
  - complex_plddt:   0–1 scale (NOT 0–100)
  - pair_chains_iptm: uses NUMERIC keys "0","1","2" (not "A","B","C")
  - For single-chain: use ptm (not iptm, which is 0 for 1 chain)
  - For 3-chain:      ipTM(MB↔VH1) = pair["0"]["1"], ipTM(MB↔Nb) = pair["1"]["2"]
"""

from __future__ import annotations

import json
import glob
import os
import shutil
from pathlib import Path
from typing import Literal

import numpy as np


# ── Kabsch-aligned RMSD ───────────────────────────────────────────────────────

def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Compute Kabsch-aligned RMSD between two (N×3) Cα coordinate arrays."""
    P = P - P.mean(0)
    Q = Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(P.T @ Q)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return float(np.sqrt(np.mean(np.sum((P @ R.T - Q) ** 2, axis=1))))


# ── Structure loading ─────────────────────────────────────────────────────────

def design_chain_id(aa) -> str:
    """Return the connected VH1-MB-Nb design chain (always A when present)."""
    chains = sorted(set(aa.chain_id))
    return "A" if "A" in chains else chains[0]


def load_ca_rfd3(cif_path: Path, res_min: int, res_max: int) -> np.ndarray:
    """
    Load Cα coordinates from an RFdiffusion3 CIF (single chain A,
    residues 1-285 in the connected chain design).
    Uses atomworks (available in foundry Docker).
    """
    from atomworks.io.utils.io_utils import load_any
    raw = load_any(str(cif_path))
    aa = raw[0] if hasattr(raw, "__getitem__") else raw
    ch = design_chain_id(aa)
    mask = (
        (aa.chain_id == ch)
        & np.isin(aa.atom_name, ["CA"])
        & (aa.res_id >= res_min)
        & (aa.res_id <= res_max)
    )
    return aa.coord[mask]


def load_ca_boltz(cif_path: Path, chain_id: str, res_min: int = 1, res_max: int = 9999) -> np.ndarray:
    """
    Load Cα coordinates from a Boltz-2 CIF.
    Uses biotite (available everywhere).
    """
    import biotite.structure.io.pdbx as pdbx
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cif = pdbx.CIFFile.read(str(cif_path))
        aa = pdbx.get_structure(cif, model=1, include_bonds=False)
    mask = (
        (aa.atom_name == "CA")
        & (aa.chain_id == chain_id)
        & (aa.res_id >= res_min)
        & (aa.res_id <= res_max)
    )
    return aa.coord[mask]


# ── Boltz-2 path discovery ────────────────────────────────────────────────────

def find_boltz_prediction_dir(boltz_out_dir: Path, design_name: str) -> Path | None:
    """
    Find the Boltz-2 prediction directory for a design.

    Boltz-2 output layout:
        <boltz_out_dir>/
            boltz_results_<input_dir_name>/
                predictions/
                    <design_name>/
                        confidence_<design_name>_model_0.json
                        <design_name>_model_0.cif
    """
    # The subdir is named after the input directory: "boltz_results_<inputs_dir>"
    for pred_base in boltz_out_dir.glob("boltz_results_*/predictions"):
        candidate = pred_base / design_name
        if candidate.exists():
            return candidate

    # Fallback: exhaustive search
    matches = list(boltz_out_dir.glob(f"**/predictions/{design_name}"))
    return matches[0] if matches else None


def parse_boltz_confidence(
    pred_dir: Path,
    design_name: str,
    mode: Literal["single_chain", "three_chain"] = "single_chain",
) -> dict:
    """
    Parse Boltz-2 confidence JSON for a design.

    Returns dict with:
        plddt_01     : complex_plddt on 0–1 scale
        plddt_pct    : complex_plddt as percentage (0–100)
        ptm          : overall pTM (use for single_chain)
        iptm_overall : overall ipTM (0 for single_chain)
        iptm_MB_VH1  : ipTM between chain 1 (MB) and chain 0 (VH1), 3-chain only
        iptm_MB_Nb   : ipTM between chain 1 (MB) and chain 2 (Nb), 3-chain only
    """
    conf_path = pred_dir / f"confidence_{design_name}_model_0.json"
    if not conf_path.exists():
        return {"plddt_01": 0.0, "plddt_pct": 0.0, "ptm": 0.0,
                "iptm_overall": 0.0, "iptm_MB_VH1": 0.0, "iptm_MB_Nb": 0.0}

    c = json.load(open(conf_path))
    plddt_01 = float(c.get("complex_plddt", 0) or 0)
    ptm      = float(c.get("ptm",           0) or 0)
    iptm_all = float(c.get("iptm",          0) or 0)

    # Chain-pair ipTM uses NUMERIC keys "0","1","2" (not "A","B","C")
    pair = c.get("pair_chains_iptm", {})
    iptm_MB_VH1 = float(
        pair.get("0", {}).get("1", 0) or pair.get("1", {}).get("0", 0)
    )
    iptm_MB_Nb = float(
        pair.get("1", {}).get("2", 0) or pair.get("2", {}).get("1", 0)
    )

    return {
        "plddt_01":     round(plddt_01, 4),
        "plddt_pct":    round(plddt_01 * 100, 1),  # for human readability
        "ptm":          round(ptm,      4),
        "iptm_overall": round(iptm_all, 4),
        "iptm_MB_VH1":  round(iptm_MB_VH1, 4),
        "iptm_MB_Nb":   round(iptm_MB_Nb,  4),
    }


# ── Global assembly RMSD ───────────────────────────────────────────────────────

def assembly_length(rfd3_cif: Path, vh1_end: int = 115, nb_len: int = 115) -> int:
    """Return total residue count on the connected VH1-MB-Nb design chain."""
    from atomworks.io.utils.io_utils import load_any
    raw = load_any(str(rfd3_cif))
    aa = raw[0] if hasattr(raw, "__getitem__") else raw
    ch = design_chain_id(aa)
    return len(set(aa.res_id[aa.chain_id == ch]))


def compute_global_rmsd_single_chain(
    rfd3_cif: Path,
    boltz_cif: Path,
    vh1_end: int = 115,
    nb_len: int = 115,
) -> float:
    """
    Global Kabsch RMSD over all Cα in the VH1 + minibinder + nanobody assembly.

    Both RFd3 and Boltz-2 use a single connected chain A with residues:
        VH1 (1–vh1_end) | MB (vh1_end+1 …) | Nb (last nb_len residues)
    """
    try:
        total = assembly_length(rfd3_cif, vh1_end, nb_len)
        rv = load_ca_rfd3(rfd3_cif, 1, total)
        bv = load_ca_boltz(boltz_cif, "A", 1, total)
        n = min(len(rv), len(bv))
        if n < vh1_end + nb_len + 5:
            return 999.0
        return kabsch_rmsd(bv[:n], rv[:n])

    except Exception as e:
        print(f"  global RMSD error: {e}")
        return 999.0


def compute_global_rmsd_three_chain(
    rfd3_cif: Path,
    boltz_cif: Path,
    vh1_end: int = 115,
    nb_len: int = 115,
) -> float:
    """
    Global Kabsch RMSD for a 3-chain Boltz-2 prediction vs single-chain RFd3.

    RFd3 is one connected chain; Boltz-2 has chains A (VH1), B (MB), C (Nb).
    All domain Cα atoms are stacked and aligned together.
    """
    try:
        total = assembly_length(rfd3_cif, vh1_end, nb_len)
        mb_len = total - vh1_end - nb_len
        mb_start = vh1_end + 1
        mb_end = vh1_end + mb_len
        nb_start = mb_end + 1

        rv = load_ca_rfd3(rfd3_cif, 1, vh1_end)
        rm = load_ca_rfd3(rfd3_cif, mb_start, mb_end)
        rn = load_ca_rfd3(rfd3_cif, nb_start, total)
        Q = np.vstack([rv, rm, rn])

        bv = load_ca_boltz(boltz_cif, "A", 1, vh1_end)
        bm = load_ca_boltz(boltz_cif, "B", mb_start, mb_end)
        bn = load_ca_boltz(boltz_cif, "C", nb_start, total)
        P = np.vstack([bv, bm, bn])

        n = min(len(P), len(Q))
        if n < vh1_end + nb_len + 5:
            return 999.0
        return kabsch_rmsd(P[:n], Q[:n])

    except Exception as e:
        print(f"  global RMSD error (3-chain): {e}")
        return 999.0


# Backwards-compatible aliases (old metric was anchor-aligned minibinder-only).
compute_scrmsdsingle_chain = compute_global_rmsd_single_chain
compute_scrmsd_three_chain = compute_global_rmsd_three_chain


# ── Main scoring function ─────────────────────────────────────────────────────

def score_all_designs(
    pipeline_dir: Path,
    boltz_mode: Literal["single_chain", "three_chain"] = "single_chain",
    plddt_threshold: float = 0.60,   # 0-1 scale — Boltz-2 reports 0→1, NOT 0→100
    global_rmsd_threshold: float = 20.0,
    n_top_cifs: int = 10,
    *,
    rfd3_subdir: str | None = None,
    mpnn_subdir: str = "mpnn",
    boltz_subdir: str = "boltz_outputs",
    results_filename: str = "final_results.json",
    copy_top_cifs: bool = True,
    structs_subdir: str = "top_structures",
) -> list[dict]:
    """
    Score all designs in the pipeline directory.

    Parameters
    ----------
    pipeline_dir      Root pipeline directory (contains outputs/, boltz_outputs/, final/).
    boltz_mode        'single_chain' or 'three_chain'.
    plddt_threshold   Minimum pLDDT to pass (0-1 scale).
    global_rmsd_threshold  Maximum global assembly RMSD to pass (Å).
    n_top_cifs        Number of top designs to copy CIF files for.
    """
    MPNN_OUT   = pipeline_dir / "outputs" / mpnn_subdir
    rfd3_sub = rfd3_subdir or os.environ.get("RFD3_ACTIVE_SUBDIR", "rfd3")
    RFD3_OUT   = pipeline_dir / "outputs" / rfd3_sub
    BOLTZ_OUT  = pipeline_dir / boltz_subdir
    FINAL      = pipeline_dir / "final"
    STRUCTS    = pipeline_dir / structs_subdir
    for d in [FINAL, STRUCTS]:
        d.mkdir(parents=True, exist_ok=True)

    # Load MPNN sequences and build top-N design list
    mpnn_data = json.load(open(MPNN_OUT / "all_sequences.json"))
    mpnn_data.sort(key=lambda x: x["sequence_recovery"])
    seen, top_designs = set(), []
    for r in mpnn_data:
        if r["backbone"] not in seen:
            seen.add(r["backbone"])
            top_designs.append(r)

    # Map design name → Boltz-2 prediction dir
    pred_dirs: dict[str, Path] = {}
    for pred_base in BOLTZ_OUT.glob("boltz_results_*/predictions"):
        for d in pred_base.iterdir():
            if d.is_dir():
                pred_dirs[d.name] = d

    print(f"Scoring {len(top_designs)} designs (Boltz-2 mode: {boltz_mode})")
    print(f"RFd3 structures: {RFD3_OUT}")
    print(f"pLDDT threshold: >{plddt_threshold*100:.0f}%  |  global RMSD threshold: <{global_rmsd_threshold}Å")
    print(f"Found {len(pred_dirs)} Boltz-2 prediction directories")

    results = []
    for idx, r in enumerate(top_designs):
        bb, mb, rec = r["backbone"], r["minibinder_sequence"], r["sequence_recovery"]
        rank = idx + 1
        name = f"rank{rank:02d}_{bb}_rec{rec:.3f}"
        pred_dir = pred_dirs.get(name)

        # Parse confidence
        conf = parse_boltz_confidence(pred_dir, name, mode=boltz_mode) if pred_dir else {}
        plddt_01  = conf.get("plddt_01", 0.0)
        ptm       = conf.get("ptm", 0.0)
        iptm_vh1  = conf.get("iptm_MB_VH1", 0.0)
        iptm_nb   = conf.get("iptm_MB_Nb",  0.0)

        # Compute global assembly RMSD (VH1 + MB + Nb)
        rfd3_cif = RFD3_OUT / f"{bb}.cif"
        global_rmsd = 999.0
        if pred_dir and rfd3_cif.exists():
            boltz_cif = pred_dir / f"{name}_model_0.cif"
            if boltz_cif.exists():
                if boltz_mode == "single_chain":
                    global_rmsd = compute_global_rmsd_single_chain(rfd3_cif, boltz_cif)
                else:
                    global_rmsd = compute_global_rmsd_three_chain(rfd3_cif, boltz_cif)

        # Pass/fail
        if boltz_mode == "single_chain":
            passed = global_rmsd < global_rmsd_threshold and plddt_01 > plddt_threshold
        else:
            passed = (global_rmsd < global_rmsd_threshold and plddt_01 > plddt_threshold
                      and iptm_vh1 > 0.3 and iptm_nb > 0.3)

        results.append({
            "rank": rank,
            "backbone": bb,
            "minibinder_sequence": mb,
            "sequence_recovery": round(rec, 3),
            "boltz_plddt_pct": round(plddt_01 * 100, 1),
            "boltz_ptm": round(ptm, 4),
            "iptm_MB_VH1": round(iptm_vh1, 4),
            "iptm_MB_Nb":  round(iptm_nb,  4),
            "global_rmsd_A": round(global_rmsd, 3),
            "pass_filter": passed,
        })

    results.sort(key=lambda x: (x["global_rmsd_A"] if x["global_rmsd_A"] < 900 else 999,
                                 -x["boltz_plddt_pct"]))
    for i, r in enumerate(results):
        r["rank"] = i + 1

    results_path = FINAL / results_filename
    json.dump(results, open(results_path, "w"), indent=2)
    passing = [r for r in results if r["pass_filter"]]
    fasta_name = "validated_minibinders.fasta" if results_filename == "final_results.json" else f"validated_{results_filename.replace('.json', '')}.fasta"
    with open(FINAL / fasta_name, "w") as f:
        for r in passing:
            f.write(f">{r['backbone']}__globalRMSD{r['global_rmsd_A']:.2f}"
                    f"__pLDDT{r['boltz_plddt_pct']:.0f}\n"
                    f"{r['minibinder_sequence']}\n")

    if copy_top_cifs:
        for i, r in enumerate(results[:n_top_cifs]):
            bb = r["backbone"]
            rec = r["sequence_recovery"]
            rank = r["rank"]
            name = f"rank{rank:02d}_{bb}_rec{rec:.3f}"
            pred_dir = pred_dirs.get(name)
            rfd3_src = RFD3_OUT / f"{bb}.cif"
            if rfd3_src.exists():
                shutil.copy(rfd3_src, STRUCTS / f"top{i+1:02d}_{bb}_rfd3.cif")
            if pred_dir:
                bs = pred_dir / f"{name}_model_0.cif"
                if bs.exists():
                    shutil.copy(bs, STRUCTS / f"top{i+1:02d}_{bb}_boltz2.cif")

        summary_name = "top_summary.tsv" if results_filename == "final_results.json" else f"top_{results_filename.replace('.json', '')}.tsv"
        with open(STRUCTS / summary_name, "w") as f:
            f.write("rank\tbackbone\tglobal_rmsd_A\tboltz_plddt_pct\tboltz_ptm\t"
                    "iptm_MB_VH1\tiptm_MB_Nb\tsequence_recovery\tminibinder_sequence\n")
            for i, r in enumerate(results[:n_top_cifs]):
                f.write(f"{i+1}\t{r['backbone']}\t{r['global_rmsd_A']}\t{r['boltz_plddt_pct']}\t"
                        f"{r['boltz_ptm']}\t{r['iptm_MB_VH1']}\t{r['iptm_MB_Nb']}\t"
                        f"{r['sequence_recovery']}\t{r['minibinder_sequence']}\n")

    # Print table
    print(f"\n{'='*76}")
    col = "pTM" if boltz_mode == "single_chain" else "↔VH1/↔Nb"
    print(f"{'R':>3}  {'Backbone':>28}  {'globalRMSD':>10}  {'pLDDT%':>6}  {col:>10}  Pass")
    print("  " + "-"*74)
    for r in results[:20]:
        if boltz_mode == "single_chain":
            extra = f"{r['boltz_ptm']:>10.4f}"
        else:
            extra = f"{r['iptm_MB_VH1']:>5.3f}/{r['iptm_MB_Nb']:>5.3f}"
        print(f"{r['rank']:>3}.  {r['backbone']:>28}  {r['global_rmsd_A']:>10.3f}  "
              f"{r['boltz_plddt_pct']:>6.1f}  {extra}{'  ✓' if r['pass_filter'] else ''}")

    print(f"\nPassed (global RMSD<{global_rmsd_threshold}Å & pLDDT>{plddt_threshold*100:.0f}%): "
          f"{len(passing)}/{len(results)}")
    print(f"Results written to {results_path}")
    if passing:
        p = passing[0]
        print(f"\nBest: {p['backbone']}  global RMSD={p['global_rmsd_A']:.2f}Å  pLDDT={p['boltz_plddt_pct']:.0f}%")
        print(f"  Sequence: {p['minibinder_sequence']}")

    return results


# ── CLI ────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse, os
    p = argparse.ArgumentParser(description="Score Boltz-2 + RFd3 pipeline outputs.")
    p.add_argument("--pipeline-dir", default=os.path.expanduser("~/pipeline"),
                   help="Root pipeline directory on the instance.")
    p.add_argument("--mode", choices=["single_chain", "three_chain"], default="single_chain",
                   help="Boltz-2 prediction mode (single_chain = connected polypeptide).")
    p.add_argument("--plddt-threshold", type=float, default=0.60,
                   help="Min pLDDT to pass (0-1 scale; default 0.60 = 60%%).")
    p.add_argument("--global-rmsd-threshold", type=float, default=20.0,
                   help="Max global assembly RMSD to pass (Å).")
    p.add_argument("--scrmsd-threshold", type=float, default=None,
                   help="Deprecated alias for --global-rmsd-threshold.")
    p.add_argument("--n-top-cifs", type=int, default=10,
                   help="Number of top designs to copy CIF files for.")
    p.add_argument("--rfd3-subdir", default=None,
                   help="Subdirectory under outputs/ for RFd3 CIFs (default: rfd3).")
    p.add_argument("--mpnn-subdir", default="mpnn",
                   help="Subdirectory under outputs/ for MPNN JSON.")
    p.add_argument("--boltz-subdir", default="boltz_outputs",
                   help="Boltz-2 output directory name under pipeline root.")
    p.add_argument("--results-file", default="final_results.json",
                   help="Filename for JSON results under final/.")
    p.add_argument("--no-copy-cifs", action="store_true",
                   help="Skip copying top CIF files (use for intermediate round-1 scoring).")
    args = p.parse_args()

    rmsd_threshold = (args.scrmsd_threshold if args.scrmsd_threshold is not None
                      else args.global_rmsd_threshold)
    score_all_designs(
        pipeline_dir=Path(args.pipeline_dir),
        boltz_mode=args.mode,
        plddt_threshold=args.plddt_threshold,
        global_rmsd_threshold=rmsd_threshold,
        n_top_cifs=args.n_top_cifs,
        rfd3_subdir=args.rfd3_subdir,
        mpnn_subdir=args.mpnn_subdir,
        boltz_subdir=args.boltz_subdir,
        results_filename=args.results_file,
        copy_top_cifs=not args.no_copy_cifs,
    )
