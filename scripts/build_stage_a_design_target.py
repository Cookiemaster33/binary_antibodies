#!/usr/bin/env python3
"""
build_stage_a_design_target.py
------------------------------
Build the Stage A RFd3 design target for the hidden trivalent minibinder approach.

Stage A designs a **separate** minibinder that bridges CH1 (chain C) and VL (chain B)
hotspots while the **full Fab** (A–D) and epitope stub (T) provide fixed steric context.
VL/CL are translated apart from VH/CH1 to open space between the two binding surfaces.

After Stage 0, pass the refined Fab with --fab-pdb (chains A–D or fused H/L from top hit).

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
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from binary_antibodies.fab_hidden_switch import (  # noqa: E402
    DEFAULT_STAGE_A_SEPARATION_A,
    EPITOPE_SEQ,
    build_stage_a_design_target_pdb,
    stage_a_contig,
    stage_a_fixed_atoms,
    stage_a_hotspots,
    write_json,
)

OUT_PDB = ROOT / "structures" / "domains" / "fab_hidden_minibinder_stage_a.pdb"
OUT_JSON = ROOT / "structures" / "interface" / "stage_a_hidden_minibinder_config.json"
MB_LENGTH_RANGE = "35-55"


def build_config(
    vh_len: int,
    vl_len: int,
    ch1_len: int,
    cl_len: int,
    source: str,
    separation_a: float,
    epitope_len: int = len(EPITOPE_SEQ),
) -> dict:
    return {
        "stage": "A",
        "approach": "hidden_trivalent_minibinder",
        "description": (
            "Stage A: design unlinked minibinder bridging CH1 (C) and VL (B) hotspots "
            "with full Fab (A–D) + epitope stub (T) as fixed context. "
            "VL/CL translated apart from VH/CH1 to open space between binding surfaces."
        ),
        "input_pdb": str(OUT_PDB.relative_to(ROOT)),
        "source_fab": source,
        "layout": {
            "separation_a": separation_a,
            "minibinder_linked": False,
            "full_fab_context": True,
        },
        "chains": {
            "A": f"VH 1-{vh_len} (fixed)",
            "B": f"VL 1-{vl_len} (fixed, translated apart)",
            "C": f"CH1 1-{ch1_len} (fixed, hotspot)",
            "D": f"CL 1-{cl_len} (fixed, translated with VL)",
            "T": f"HER2 epitope stub 1-{epitope_len} (fixed, Stage B)",
        },
        "rfd3": {
            "contig": stage_a_contig(vh_len, vl_len, ch1_len, cl_len, MB_LENGTH_RANGE),
            "mb_length_range": MB_LENGTH_RANGE,
            "select_fixed_atoms": stage_a_fixed_atoms(vh_len, vl_len, ch1_len, cl_len, epitope_len),
            "select_hotspots": stage_a_hotspots(vh_len),
        },
        "validation": {
            "metric": "global_assembly_rmsd",
            "note": "Score Boltz refold of Fab + unlinked minibinder vs RFd3 with A–D,T fixed",
        },
        "stage_b_note": (
            "Stage B will add target-arm hotspots on chain T. A flexible linker between "
            "minibinder and Fab can be added after selecting a backbone."
        ),
        "requires_stage_0": True,
    }


def main() -> None:
    p = argparse.ArgumentParser(description="Build Stage A hidden minibinder design target.")
    p.add_argument(
        "--fab-pdb",
        type=Path,
        default=None,
        help="Stage 0 refined Fab (PDB/CIF with chains A–D or fused H/L). Default: native 1N8Z.",
    )
    p.add_argument(
        "--separation-a",
        type=float,
        default=DEFAULT_STAGE_A_SEPARATION_A,
        help=f"Å translation of VL/CL away from VH/CH1 (default {DEFAULT_STAGE_A_SEPARATION_A})",
    )
    p.add_argument("--out-pdb", type=Path, default=OUT_PDB)
    p.add_argument("--out-json", type=Path, default=OUT_JSON)
    args = p.parse_args()

    source = str(args.fab_pdb) if args.fab_pdb else "1N8Z (trastuzumab)"
    vh_len, vl_len, ch1_len, cl_len = build_stage_a_design_target_pdb(
        args.out_pdb,
        source_pdb=args.fab_pdb,
        separation_a=args.separation_a,
    )
    config = build_config(vh_len, vl_len, ch1_len, cl_len, source, args.separation_a)
    write_json(args.out_json, config)

    print(f"Wrote design target PDB → {args.out_pdb}")
    print(f"Wrote RFd3 config       → {args.out_json}")
    print(f"  contig: {config['rfd3']['contig']}")
    print(f"  separation: {args.separation_a} Å (VL/CL away from VH/CH1)")
    print(f"  hotspots: {len(config['rfd3']['select_hotspots'].split(','))} residues")
    if args.fab_pdb:
        print(f"  Fab source: {args.fab_pdb}")
    else:
        print("  WARNING: using native Fab — run Stage 0 first for production designs.")


if __name__ == "__main__":
    main()
