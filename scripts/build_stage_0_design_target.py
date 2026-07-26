#!/usr/bin/env python3
"""
build_stage_0_design_target.py
------------------------------
Stage 0: partial de-greasing of the VH–VL framework interface via split-chain
ProteinMPNN (no RFd3). CDRs, CH1, CL, and epitope stub stay native.

Outputs
-------
  structures/domains/fab_stage_0_vhvL_interface.pdb   — full Fab context (Boltz / scoring)
  structures/domains/fab_stage_0_split_mpnn.pdb     — separated Fv for MPNN
  structures/interface/stage_0_vhvL_interface_config.json

Usage
-----
    python scripts/build_stage_0_design_target.py
    python scripts/build_stage_0_design_target.py --source-pdb path/to/parent.pdb
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
    build_fab_context_pdb,
    build_split_fv_mpnn_pdb,
    extract_chain_sequences,
    interface_closure_core,
    interface_degrease_residues,
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
) -> dict:
    vh_core, vl_core = interface_closure_core(vh_len, vl_len)
    designed = split_mpnn_designed_residues(vh_len, vl_len)

    return {
        "stage": "0",
        "approach": "split_mpnn_partial_degrease",
        "description": (
            "Stage 0: split-chain ProteinMPNN partial de-greasing of the VH–VL "
            "framework interface (rim residues only; closure core fixed). Skips RFd3. "
            "Goal: weaken apo VH–VL coupling while preserving holo closure."
        ),
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
            "designed_residues": designed,
            "vh_closure_core": [f"A{r}" for r in vh_core],
            "vl_closure_core": [f"B{r}" for r in vl_core],
            "vh_degrease": interface_degrease_residues("A", vh_len, vl_len),
            "vl_degrease": interface_degrease_residues("B", vh_len, vl_len),
            "n_sequences": 32,
        },
        "native_chain_sequences": {},
        "validation": {
            "apo_fold": "Boltz A+B+C+D — prefer weaker VH–VL contacts vs native",
            "holo_fold": "Boltz A+B+C+D+T — prefer closed Fv engaging epitope",
            "metrics": [
                "vh_vl_interface_contacts",
                "holo_minus_apo_contact_delta",
                "apo_fraction_of_native_contacts",
                "holo_fraction_of_native_contacts",
                "holo_fv_framework_rmsd_A",
                "cdr_epitope_contacts",
                "vh_vl_interface_clashes",
                "wedge_suspect",
            ],
            "filters": {
                "max_apo_fraction_of_native_contacts": 0.45,
                "min_holo_fraction_of_native_contacts": 0.45,
                "min_holo_minus_apo_contact_delta": 15,
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
    p.add_argument(
        "--separation-A",
        type=float,
        default=DEFAULT_SPLIT_SEPARATION_A,
        help="VH–VL translation for split MPNN input (Å)",
    )
    args = p.parse_args()

    source = str(args.source_pdb) if args.source_pdb else "1N8Z (trastuzumab)"
    vh_len, vl_len, ch1_len, cl_len = build_fab_context_pdb(
        args.out_pdb,
        source_pdb=args.source_pdb,
        struct_name="stage_0",
    )
    build_split_fv_mpnn_pdb(
        args.out_split_pdb,
        source_pdb=args.out_pdb,
        separation_a=args.separation_A,
    )
    config = build_config(vh_len, vl_len, ch1_len, cl_len, source, args.separation_A)
    config["native_chain_sequences"] = extract_chain_sequences(args.out_pdb)
    write_json(args.out_json, config)

    sm = config["split_mpnn"]
    print(f"Wrote Fab context PDB  → {args.out_pdb}")
    print(f"Wrote split MPNN PDB   → {args.out_split_pdb}")
    print(f"Wrote config           → {args.out_json}")
    print(f"  separation:   {sm['separation_A']} Å")
    print(f"  degrease:     {len(sm['designed_residues'])} residues "
          f"(VH core {len(sm['vh_closure_core'])}, VL core {len(sm['vl_closure_core'])})")
    print(f"  designed:     {', '.join(sm['designed_residues'][:8])}"
          f"{'...' if len(sm['designed_residues']) > 8 else ''}")


if __name__ == "__main__":
    main()
