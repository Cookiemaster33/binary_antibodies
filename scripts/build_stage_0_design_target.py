#!/usr/bin/env python3
"""
build_stage_0_design_target.py
------------------------------
Stage 0: weaken the VH–VL framework interface while keeping CDRs, CH1, CL,
and the HER2 epitope stub fixed.

RFd3 uses partial diffusion on chains A+B with epitope hotspots on chain T so
designs remain compatible with target-bound closure geometry.

Output
------
  structures/domains/fab_stage_0_vhvL_interface.pdb
  structures/interface/stage_0_vhvL_interface_config.json

Usage
-----
    python scripts/build_stage_0_design_target.py
    python scripts/build_stage_0_design_target.py --source-pdb path/to/stage0_parent.pdb
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
    VH_INTERFACE_FW,
    VL_INTERFACE_FW,
    build_fab_context_pdb,
    cdr_fixed_atoms,
    epitope_hotspots,
    write_json,
)

OUT_PDB = ROOT / "structures" / "domains" / "fab_stage_0_vhvL_interface.pdb"
OUT_JSON = ROOT / "structures" / "interface" / "stage_0_vhvL_interface_config.json"


def build_config(vh_len: int, vl_len: int, ch1_len: int, cl_len: int, source: str) -> dict:
    fixed = cdr_fixed_atoms(vh_len, vl_len)
    fixed[f"C1-{ch1_len}"] = "ALL"
    fixed[f"D1-{cl_len}"] = "ALL"
    fixed[f"T1-{len(EPITOPE_SEQ)}"] = "ALL"

    return {
        "stage": "0",
        "approach": "vh_vl_interface_weakening",
        "description": (
            "Stage 0: partial-diffusion redesign of VH–VL framework interface "
            "(CDRs + constant domains + epitope stub fixed). Goal: reduce apo "
            "VH–VL affinity while preserving target-bound closure geometry."
        ),
        "input_pdb": str(OUT_PDB.relative_to(ROOT)),
        "source_fab": source,
        "chains": {
            "A": f"VH 1-{vh_len}",
            "B": f"VL 1-{vl_len}",
            "C": f"CH1 1-{ch1_len}",
            "D": f"CL 1-{cl_len} (steric context only)",
            "T": f"HER2 epitope stub 1-{len(EPITOPE_SEQ)} (hotspot / closure context)",
        },
        "rfd3": {
            "contig": f"A1-{vh_len}/0,B1-{vl_len}",
            "partial_t": 12.0,
            "select_hotspots": epitope_hotspots(),
            "select_fixed_atoms": fixed,
            "vh_interface_framework": [f"A{r}" for r in VH_INTERFACE_FW if r <= vh_len],
            "vl_interface_framework": [f"B{r}" for r in VL_INTERFACE_FW if r <= vl_len],
        },
        "validation": {
            "apo_fold": "Boltz chains A+B+C+D (no T) — prefer fewer VH–VL interface contacts",
            "holo_fold": "Boltz chains A+B+C+D+T — prefer closed Fv around T without clashes",
            "metrics": [
                "vh_vl_interface_contacts",
                "vh_vl_fr4_ca_distance",
                "apo_minus_holo_contact_delta",
            ],
            "filters": {
                "max_apo_interface_contacts": 12,
                "min_holo_interface_contacts": 20,
                "max_holo_clashes_4A": 0,
            },
        },
        "next_stage": (
            "Feed top Stage 0 Fab (chains A–D) into build_stage_a_design_target.py "
            "via --fab-pdb, then run Stage A minibinder design."
        ),
    }


def main() -> None:
    p = argparse.ArgumentParser(description="Build Stage 0 VH–VL interface design target.")
    p.add_argument(
        "--source-pdb",
        type=Path,
        default=None,
        help="Parent Fab PDB (default: structures/1N8Z.pdb)",
    )
    p.add_argument("--out-pdb", type=Path, default=OUT_PDB)
    p.add_argument("--out-json", type=Path, default=OUT_JSON)
    args = p.parse_args()

    source = str(args.source_pdb) if args.source_pdb else "1N8Z (trastuzumab)"
    vh_len, vl_len, ch1_len, cl_len = build_fab_context_pdb(
        args.out_pdb,
        source_pdb=args.source_pdb,
        struct_name="stage_0",
    )
    config = build_config(vh_len, vl_len, ch1_len, cl_len, source)
    write_json(args.out_json, config)

    print(f"Wrote design target PDB → {args.out_pdb}")
    print(f"Wrote RFd3 config       → {args.out_json}")
    print(f"  contig:     {config['rfd3']['contig']}")
    print(f"  partial_t:  {config['rfd3']['partial_t']} Å")
    print(f"  hotspots:   {len(config['rfd3']['select_hotspots'].split(','))} epitope residues")
    print(f"  fixed CDRs: {len(config['rfd3']['select_fixed_atoms'])} atom groups")


if __name__ == "__main__":
    main()
