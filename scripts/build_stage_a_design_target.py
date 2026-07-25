#!/usr/bin/env python3
"""
build_stage_a_design_target.py
------------------------------
Build the Stage A RFd3 design target for the hidden trivalent minibinder approach.

Stage A designs a connected VH–minibinder–VL hub nestled in a native Fab,
with hotspots on CH1 (chain C) and VL (chain B). A HER2 epitope stub (chain T)
is placed inside the VH–VL groove as fixed steric context for Stage B masking.

After Stage 0, pass the refined Fab with --fab-pdb (chains A–D from top hit).

Output
------
  structures/domains/fab_hidden_minibinder_stage_a.pdb
  structures/interface/stage_a_hidden_minibinder_config.json

Usage
-----
    python scripts/build_stage_a_design_target.py
    python scripts/build_stage_a_design_target.py --fab-pdb pipeline_results/stage_0/structures/top01_apo.cif
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from binary_antibodies.fab_hidden_switch import (  # noqa: E402
    EPITOPE_SEQ,
    build_fab_context_pdb,
    stage_a_hotspots,
    write_json,
)

OUT_PDB = ROOT / "structures" / "domains" / "fab_hidden_minibinder_stage_a.pdb"
OUT_JSON = ROOT / "structures" / "interface" / "stage_a_hidden_minibinder_config.json"


def build_config(vh_len: int, vl_len: int, ch1_len: int, cl_len: int, source: str) -> dict:
    return {
        "stage": "A",
        "approach": "hidden_trivalent_minibinder",
        "description": (
            "Stage A: design connected VH–minibinder–VL hub seated in Stage 0 Fab. "
            "Hotspots on CH1 (chain C) and VL (chain B). Epitope stub (chain T) "
            "provides steric masking context for Stage B target arm."
        ),
        "input_pdb": str(OUT_PDB.relative_to(ROOT)),
        "source_fab": source,
        "chains": {
            "A": f"VH 1-{vh_len}",
            "B": f"VL 1-{vl_len}",
            "C": f"CH1 1-{ch1_len}",
            "D": f"CL 1-{cl_len} (steric context only)",
            "T": f"HER2 epitope stub 1-{len(EPITOPE_SEQ)} (steric context, Stage B)",
        },
        "rfd3": {
            "contig": f"A1-{vh_len},35-55,B1-{vl_len}",
            "mb_length_range": "35-55",
            "select_fixed_atoms": {
                f"C1-{ch1_len}": "ALL",
                f"D1-{cl_len}": "ALL",
                f"T1-{len(EPITOPE_SEQ)}": "ALL",
            },
            "select_hotspots": stage_a_hotspots(vh_len),
        },
        "validation": {
            "metric": "global_assembly_rmsd",
            "note": "Score Boltz refold vs RFd3 for full A–MB–B chain with C,D,T as context",
        },
        "stage_b_note": (
            "Stage B will add target-arm hotspots on chain T and extend contig "
            "to include a third binding patch against the epitope stub."
        ),
        "requires_stage_0": True,
    }


def main() -> None:
    p = argparse.ArgumentParser(description="Build Stage A hidden minibinder design target.")
    p.add_argument(
        "--fab-pdb",
        type=Path,
        default=None,
        help="Stage 0 refined Fab (PDB/CIF with chains A–D). Default: native 1N8Z.",
    )
    p.add_argument("--out-pdb", type=Path, default=OUT_PDB)
    p.add_argument("--out-json", type=Path, default=OUT_JSON)
    args = p.parse_args()

    source = str(args.fab_pdb) if args.fab_pdb else "1N8Z (trastuzumab)"
    vh_len, vl_len, ch1_len, cl_len = build_fab_context_pdb(
        args.out_pdb,
        source_pdb=args.fab_pdb,
        struct_name="stage_a",
    )
    config = build_config(vh_len, vl_len, ch1_len, cl_len, source)
    write_json(args.out_json, config)

    print(f"Wrote design target PDB → {args.out_pdb}")
    print(f"Wrote RFd3 config       → {args.out_json}")
    print(f"  contig: {config['rfd3']['contig']}")
    print(f"  hotspots: {len(config['rfd3']['select_hotspots'].split(','))} residues")
    if args.fab_pdb:
        print(f"  Fab source: {args.fab_pdb}")
    else:
        print("  WARNING: using native Fab — run Stage 0 first for production designs.")


if __name__ == "__main__":
    main()
