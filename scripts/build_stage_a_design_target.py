#!/usr/bin/env python3
"""
build_stage_a_design_target.py
------------------------------
Build the Stage A RFd3 design target for the hidden trivalent minibinder approach.

Stage A designs a **separate, helix-biased** minibinder that bridges CH1 (chain C)
and VL (chain B) hotspots while the **full Fab** (A–D) and epitope stub (T) act as a
rigid target. Every Fab chain is its own fixed target chain in the contig and the
minibinder is its own designed chain, so RFd3 has no reason to move the arms.

After Stage 0, pass the refined Fab with --fab-pdb (fused H/L holo, or a structure
whose arms are already positioned — those coordinates are used verbatim).

Output
------
  structures/domains/fab_hidden_minibinder_stage_a.pdb
  structures/interface/stage_a_hidden_minibinder_config.json

Usage
-----
    python scripts/build_stage_a_design_target.py
    python scripts/build_stage0_holo_split_cif.py path/to/rank079_holo.cif
    python scripts/build_stage_a_design_target.py \\
        --fab-pdb path/to/rank079_holo_split.cif
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from binary_antibodies.fab_hidden_switch import (  # noqa: E402
    DEFAULT_MB_LENGTH_RANGE,
    DEFAULT_STAGE_A_SEPARATION_A,
    EPITOPE_SEQ,
    build_stage_a_design_target_pdb,
    infer_already_split,
    stage_a_hotspot_geometry,
    stage_a_ori_token,
    stage_a_rfd3_config,
    verify_design_target_matches_source,
    write_json,
)

OUT_PDB = ROOT / "structures" / "domains" / "fab_hidden_minibinder_stage_a.pdb"
OUT_JSON = ROOT / "structures" / "interface" / "stage_a_hidden_minibinder_config.json"


def build_config(
    vh_len: int,
    vl_len: int,
    ch1_len: int,
    cl_len: int,
    source: str,
    separation_a: float,
    mb_length_range: str,
    geometry: dict,
    ori_token: list[float] | None,
    epitope_len: int = len(EPITOPE_SEQ),
) -> dict:
    return {
        "stage": "A",
        "approach": "hidden_trivalent_minibinder",
        "description": (
            "Stage A: design an unlinked, helix-biased minibinder that bridges the CH1 (C) "
            "and VL (B) hotspot surfaces of the split Fab. The full Fab (A–D) plus epitope "
            "stub (T) is a rigid target: every Fab chain is a separate fixed target chain in "
            "the contig and the minibinder is its own designed chain, so RFd3 never has to "
            "move the arms to satisfy chain connectivity."
        ),
        "input_pdb": str(OUT_PDB.relative_to(ROOT)),
        "source_fab": source,
        "layout": {
            "separation_a": separation_a,
            "minibinder_linked": False,
            "full_fab_context": True,
            "minibinder_own_chain": True,
        },
        "chains": {
            "A": f"VH 1-{vh_len} (fixed target)",
            "B": f"VL 1-{vl_len} (fixed target, hotspot)",
            "C": f"CH1 1-{ch1_len} (fixed target, hotspot)",
            "D": f"CL 1-{cl_len} (fixed target)",
            "T": f"HER2 epitope stub 1-{epitope_len} (fixed target, Stage B)",
            "M": "designed minibinder (added by graft_stage_a_outputs.py)",
        },
        "hotspot_geometry": geometry,
        "secondary_structure": {
            "goal": "alpha-helical minibinder (helical hairpin / bundle)",
            "levers": [
                "is_non_loopy=true — RFd3's only SS lever; more helices, fewer loops and sheets",
                "low-temperature sampling (step_scale=3.0, gamma_0=0.2)",
                f"length window {mb_length_range} aa — long enough for helices that span the gap",
                "minibinder on its own chain so it folds as a binder, not a linker",
            ],
            "qc": "scripts/analyze_stage_a_minibinders.py reports helix fraction per design",
        },
        "rfd3": {
            **stage_a_rfd3_config(
                vh_len,
                vl_len,
                ch1_len,
                cl_len,
                mb_length_range,
                epitope_len,
                ori_token=ori_token,
            ),
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
        help="Stage 0 holo (fused H/L) or pre-built *_split.cif. Default: native 1N8Z.",
    )
    p.add_argument(
        "--separation-a",
        type=float,
        default=DEFAULT_STAGE_A_SEPARATION_A,
        help=f"Å translation of VL/CL away from VH/CH1 (default {DEFAULT_STAGE_A_SEPARATION_A}; "
        "ignored for *_split.cif inputs)",
    )
    p.add_argument(
        "--already-split",
        action="store_true",
        help="Fab coordinates already separated; do not re-translate VL/CL (auto for *_split.cif)",
    )
    p.add_argument(
        "--mb-length-range",
        default=DEFAULT_MB_LENGTH_RANGE,
        help=f"Minibinder length window, 'min-max' (default {DEFAULT_MB_LENGTH_RANGE}; "
        "long enough for helices that span the CH1↔VL gap)",
    )
    p.add_argument(
        "--infer-ori-strategy",
        action="store_true",
        help="Let RFd3 infer the origin from hotspots instead of pinning it to the gap midpoint",
    )
    p.add_argument("--out-pdb", type=Path, default=OUT_PDB)
    p.add_argument("--out-json", type=Path, default=OUT_JSON)
    args = p.parse_args()

    source = str(args.fab_pdb) if args.fab_pdb else "1N8Z (trastuzumab)"
    already_split = args.already_split or infer_already_split(args.fab_pdb)
    effective_sep = 0.0 if already_split else args.separation_a
    vh_len, vl_len, ch1_len, cl_len = build_stage_a_design_target_pdb(
        args.out_pdb,
        source_pdb=args.fab_pdb,
        separation_a=effective_sep,
        already_split=already_split,
    )

    passthrough: dict[str, float] | None = None
    if already_split and args.fab_pdb is not None:
        passthrough = verify_design_target_matches_source(args.out_pdb, args.fab_pdb)

    geometry = stage_a_hotspot_geometry(args.out_pdb)
    ori_token = None if args.infer_ori_strategy else stage_a_ori_token(args.out_pdb)
    config = build_config(
        vh_len,
        vl_len,
        ch1_len,
        cl_len,
        source,
        effective_sep,
        args.mb_length_range,
        geometry,
        ori_token,
    )
    if already_split:
        config["layout"]["input_already_split"] = True
    if passthrough is not None:
        config["layout"]["source_passthrough_max_dev_a"] = passthrough
    write_json(args.out_json, config)

    print(f"Wrote design target PDB → {args.out_pdb}")
    print(f"Wrote RFd3 config       → {args.out_json}")
    print(f"  contig: {config['rfd3']['contig']}")
    if already_split:
        print("  input: arms already separated — coordinates used verbatim (no re-translation)")
        if passthrough is not None:
            worst = max(passthrough.values())
            print(f"  verified passthrough: max CA deviation vs source = {worst:.6f} Å")
    else:
        print(f"  separation: {args.separation_a} Å (VL/CL away from VH/CH1)")
    print(f"  hotspots: {len(config['rfd3']['select_hotspots'].split(','))} residues")
    print(
        f"  CH1↔VL hotspot gap: {geometry['centroid_separation_a']} Å centroid-to-centroid, "
        f"{geometry['closest_hotspot_approach_a']} Å closest approach"
    )
    print(f"  minibinder length: {args.mb_length_range} aa")
    print(f"  helical conditioning: is_non_loopy={config['rfd3']['is_non_loopy']}, "
          f"sampler={config['rfd3']['sampler']}")
    if ori_token is not None:
        print(f"  ori_token (gap midpoint): {ori_token}")
    else:
        print("  ori_token: inferred from hotspots")
    if args.fab_pdb:
        print(f"  Fab source: {args.fab_pdb}")
    else:
        print("  WARNING: using native Fab — run Stage 0 first for production designs.")


if __name__ == "__main__":
    main()
