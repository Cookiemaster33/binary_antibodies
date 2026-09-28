#!/usr/bin/env python3
"""
validate_stage_a_fab_fixed.py
-----------------------------
Check whether Stage A RFd3 outputs preserved the input Fab domain coordinates.

For legacy runs (full Fab in contig with /0 breaks), domains keep internal
structure but reorient ~15–20 Å relative to each other when aligned globally.

Usage
-----
    python scripts/validate_stage_a_fab_fixed.py \\
        --input pipeline_results/stage_a_rank079_split_cif_v4/inputs/fab_hidden_minibinder_stage_a.pdb \\
        --cif pipeline_results/stage_a_rank079_split_cif_v4/structures/sample_cifs/sa_b000_000.cif
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def read_cif_seq(cif: Path) -> str:
    text = cif.read_text()
    m = re.search(r"_entity_poly\.pdbx_seq_one_letter_code\s*\n;([^\n;]+(?:\n[^;]+)*)", text, re.S)
    if not m:
        raise ValueError(f"No sequence in {cif}")
    return re.sub(r"\s+", "", m.group(1))


def chain_ca(struct, chain_id: str) -> np.ndarray:
    from Bio.PDB import MMCIFParser, PDBParser

    model = list(struct.get_models())[0]
    chain = model[chain_id]
    return np.array([r["CA"].get_coord() for r in chain if r.id[0] == " " and "CA" in r])


def kabsch_rmsd(mobile: np.ndarray, target: np.ndarray) -> float:
    mc = mobile - mobile.mean(0)
    tc = target - target.mean(0)
    v, _, wt = np.linalg.svd(mc.T @ tc)
    d = np.sign(np.linalg.det(v @ wt))
    r = v @ np.diag([1.0, 1.0, d]) @ wt
    return float(np.sqrt(((mc @ r - tc) ** 2).sum(1).mean()))


def main() -> None:
    p = argparse.ArgumentParser(description="Validate Fab fixed coordinates in Stage A RFd3 output.")
    p.add_argument("--input", type=Path, required=True, help="Stage A input PDB (chains A–D,T)")
    p.add_argument("--cif", type=Path, required=True, help="RFd3 output CIF")
    p.add_argument("--layout", choices=["legacy", "unindex", "auto"], default="auto")
    args = p.parse_args()

    from Bio.PDB import MMCIFParser, PDBParser

    inp = PDBParser(QUIET=True).get_structure("in", str(args.input))
    out = MMCIFParser(QUIET=True).get_structure("out", str(args.cif))

    out_model = list(out.get_models())[0]
    out_chains = list(out_model.child_dict)
    seq = read_cif_seq(args.cif)

    vh_len, vl_len, ch1_len, cl_len = 113, 107, 107, 107
    layout = args.layout
    if layout == "auto":
        layout = "unindex" if len(seq) < vh_len + vl_len + ch1_len + cl_len + 35 else "legacy"

    print(f"Input: {args.input}")
    print(f"Output: {args.cif}")
    print(f"Layout: {layout}  (output len={len(seq)})")

    if layout == "unindex":
        mb_len = len(seq) - vl_len - ch1_len
        out_chain = list(out_model.get_chains())[0]
        out_ca = np.array([r["CA"].get_coord() for r in out_chain if r.id[0] == " " and "CA" in r])
        vl_out, ch1_out = out_ca[:vl_len], out_ca[vl_len + mb_len : vl_len + mb_len + ch1_len]
        vl_in, ch1_in = chain_ca(inp, "B"), chain_ca(inp, "C")
        print(f"  MB length: {mb_len}")
        print(f"  VL RMSD (per-domain align): {kabsch_rmsd(vl_out, vl_in):.2f} Å")
        print(f"  CH1 RMSD (per-domain align): {kabsch_rmsd(ch1_out, ch1_in):.2f} Å")
        print("  VH/CL/T not in output polymer (unindexed context only).")
        return

    # Legacy: VH|VL|MB|CH1|CL on one chain
    fixed = vh_len + vl_len + ch1_len + cl_len
    mb_len = len(seq) - fixed
    out_chain = list(out_model.get_chains())[0]
    out_ca = np.array([r["CA"].get_coord() for r in out_chain if r.id[0] == " " and "CA" in r])
    segs = [
        ("VH", "A", 0, vh_len),
        ("VL", "B", vh_len, vh_len + vl_len),
        ("CH1", "C", vh_len + vl_len + mb_len, vh_len + vl_len + mb_len + ch1_len),
        ("CL", "D", vh_len + vl_len + mb_len + ch1_len, vh_len + vl_len + mb_len + ch1_len + cl_len),
    ]
    print(f"  MB length: {mb_len}")
    print("  Per-domain RMSD (each domain aligned separately):")
    for name, cid, lo, hi in segs:
        rmsd = kabsch_rmsd(out_ca[lo:hi], chain_ca(inp, cid))
        print(f"    {name}: {rmsd:.2f} Å")

    vh_in = chain_ca(inp, "A")
    vh_out = out_ca[:vh_len]
    mc = vh_out - vh_out.mean(0)
    tc = vh_in - vh_in.mean(0)
    v, _, wt = np.linalg.svd(mc.T @ tc)
    d = np.sign(np.linalg.det(v @ wt))
    r = v @ np.diag([1.0, 1.0, d]) @ wt
    out_aligned = (out_ca - vh_out.mean(0)) @ r

    print("  After global VH alignment (legacy runs often fail here):")
    for name, cid, lo, hi in segs:
        if name == "VH":
            continue
        rmsd = float(np.sqrt(((out_aligned[lo:hi] - chain_ca(inp, cid)) ** 2).sum(1).mean()))
        print(f"    {name}: {rmsd:.2f} Å")


if __name__ == "__main__":
    main()
