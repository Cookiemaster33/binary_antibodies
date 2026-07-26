#!/usr/bin/env python3
"""
build_stage_0_design_target.py
------------------------------
Stage 0: split-chain ProteinMPNN on the VH–VL framework interface (no RFd3).

Modes
-----
  --whole-interface   redesign all ~23 interface FW residues (default for brute-force)
  --conservative      partial de-grease (3.35 Å core, 14 rim residues)

Outputs
-------
  structures/domains/fab_stage_0_vhvL_interface.pdb
  structures/domains/fab_stage_0_split_mpnn.pdb
  structures/interface/stage_0_vhvL_interface_config.json
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from binary_antibodies.fab_hidden_switch import (  # noqa: E402
    DEFAULT_SPLIT_SEPARATION_A,
    EPITOPE_SEQ,
    INTERFACE_CORE_MAX_HEAVY_A,
    INTERFACE_CORE_MAX_HEAVY_A_AGGRESSIVE,
    build_fab_context_pdb,
    build_split_fv_mpnn_pdb,
    extract_chain_sequences,
    interface_closure_core,
    interface_degrease_residues,
    interface_design_residues,
    split_mpnn_designed_residues,
    write_json,
)

OUT_PDB = ROOT / "structures" / "domains" / "fab_stage_0_vhvL_interface.pdb"
OUT_SPLIT_PDB = ROOT / "structures" / "domains" / "fab_stage_0_split_mpnn.pdb"
OUT_JSON = ROOT / "structures" / "interface" / "stage_0_vhvL_interface_config.json"


def build_config(
    vh_len: int,
    vl_len: int,
    ch1_len: int,
    cl_len: int,
    source: str,
    separation_a: float,
    designed: list[str],
    approach: str,
    description: str,
    n_mpnn_seqs: int,
    top_n_boltz: int,
    core_max_heavy_a: float | None,
) -> dict:
    vh_core, vl_core = (
        interface_closure_core(vh_len, vl_len, max_heavy_a=core_max_heavy_a)
        if core_max_heavy_a is not None
        else ([], [])
    )

    return {
        "stage": "0",
        "approach": approach,
        "description": description,
        "input_pdb": str(OUT_PDB.relative_to(ROOT)),
        "split_mpnn_pdb": str(OUT_SPLIT_PDB.relative_to(ROOT)),
        "source_fab": source,
        "chains": {
            "A": f"VH 1-{vh_len}",
            "B": f"VL 1-{vl_len}",
            "C": f"CH1 1-{ch1_len}",
            "D": f"CL 1-{cl_len} (steric context only)",
            "T": f"HER2 epitope stub 1-{len(EPITOPE_SEQ)} (holo closure context)",
        },
        "split_mpnn": {
            "separation_A": separation_a,
            "core_max_heavy_A": core_max_heavy_a,
            "designed_residues": designed,
            "vh_closure_core": [f"A{r}" for r in vh_core],
            "vl_closure_core": [f"B{r}" for r in vl_core],
            "n_sequences": n_mpnn_seqs,
            "top_n_boltz": top_n_boltz,
            "mpnn_rank_by": "static_vh_vl_interface_contacts_asc",
        },
        "native_chain_sequences": {},
        "validation": {
            "mode": "holo_only_static_interface",
            "static_interface": "Weighted VH–VL contacts on native Fab geometry from designed A/B sequence",
            "holo_fold": "Boltz A+B+C+D+T — epitope engagement and Fab-like geometry",
            "metrics": [
                "static_vh_vl_interface_contacts",
                "static_fraction_of_native_contacts",
                "holo_fv_framework_rmsd_A",
                "cdr_epitope_contacts",
                "vh_vl_interface_clashes",
            ],
            "filters": {
                "max_static_fraction_of_native_contacts": 0.45,
                "max_holo_clashes_4A": 50,
                "max_holo_fv_framework_rmsd_A": 3.5,
                "min_holo_cdr_epitope_contacts": 6,
                "max_holo_vh_vl_interface_clashes": 0,
            },
        },
        "next_stage": (
            "Feed top Stage 0 Fab (chains A–D) into build_stage_a_design_target.py "
            "via --fab-pdb, then run Stage A minibinder design."
        ),
    }


def main() -> None:
    p = argparse.ArgumentParser(description="Build Stage 0 split-MPNN design targets.")
    p.add_argument("--source-pdb", type=Path, default=None, help="Parent Fab PDB")
    p.add_argument("--out-pdb", type=Path, default=OUT_PDB)
    p.add_argument("--out-split-pdb", type=Path, default=OUT_SPLIT_PDB)
    p.add_argument("--out-json", type=Path, default=OUT_JSON)
    p.add_argument("--separation-A", type=float, default=DEFAULT_SPLIT_SEPARATION_A)
    p.add_argument(
        "--conservative",
        action="store_true",
        help="Partial de-grease with 3.35 Å core (14 rim residues)",
    )
    p.add_argument(
        "--aggressive",
        action="store_true",
        help="Aggressive de-grease with 3.20 Å core (19 rim residues)",
    )
    p.add_argument("--core-max-heavy-A", type=float, default=None)
    p.add_argument("--n-mpnn-seqs", type=int, default=1000)
    p.add_argument("--top-n-boltz", type=int, default=100)
    args = p.parse_args()

    source = str(args.source_pdb) if args.source_pdb else "1N8Z (trastuzumab)"
    vh_len, vl_len, ch1_len, cl_len = build_fab_context_pdb(
        args.out_pdb, source_pdb=args.source_pdb, struct_name="stage_0"
    )
    build_split_fv_mpnn_pdb(args.out_split_pdb, source_pdb=args.out_pdb, separation_a=args.separation_A)

    if args.conservative:
        core = INTERFACE_CORE_MAX_HEAVY_A
        designed = split_mpnn_designed_residues(vh_len, vl_len, max_heavy_a=core)
        approach = "split_mpnn_partial_degrease"
        desc = f"Partial de-grease: {len(designed)} rim residues; core ≤ {core} Å fixed."
    elif args.aggressive:
        core = INTERFACE_CORE_MAX_HEAVY_A_AGGRESSIVE
        designed = split_mpnn_designed_residues(vh_len, vl_len, max_heavy_a=core)
        approach = "split_mpnn_aggressive_degrease"
        desc = f"Aggressive de-grease: {len(designed)} rim residues; core ≤ {core} Å fixed."
    else:
        core = None
        designed = interface_design_residues(vh_len, vl_len)
        approach = "split_mpnn_whole_interface"
        desc = (
            f"Whole-interface: all {len(designed)} VH/VL framework interface residues redesigned; "
            f"top {args.top_n_boltz} by static VH–VL weakening → holo Boltz."
        )

    config = build_config(
        vh_len, vl_len, ch1_len, cl_len, source, args.separation_A,
        designed, approach, desc, args.n_mpnn_seqs, args.top_n_boltz, core,
    )
    config["native_chain_sequences"] = extract_chain_sequences(args.out_pdb)
    write_json(args.out_json, config)

    sm = config["split_mpnn"]
    print(f"Wrote Fab context PDB  → {args.out_pdb}")
    print(f"Wrote split MPNN PDB   → {args.out_split_pdb}")
    print(f"Wrote config           → {args.out_json}")
    print(f"  approach:     {config['approach']}")
    print(f"  designed:     {len(sm['designed_residues'])} residues")
    print(f"  MPNN / Boltz: {sm['n_sequences']} → top {sm['top_n_boltz']}")


if __name__ == "__main__":
    main()
