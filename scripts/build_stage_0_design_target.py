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
    PISA_CORE_MIN_BURIED_SASA_A2,
    PISA_CORE_MIN_BURIED_SASA_A2_AGGRESSIVE,
    VH_INTERFACE_FW,
    VL_INTERFACE_FW,
    build_fab_context_pdb,
    build_split_fv_mpnn_pdb,
    cdr_residue_set,
    extract_chain_sequences,
    interface_closure_core,
    interface_design_residues,
    pisa_closure_core,
    split_mpnn_designed_residues,
    write_json,
)
from binary_antibodies.pisa_scoring import identify_fv_interface_residues  # noqa: E402

OUT_PDB = ROOT / "structures" / "domains" / "fab_stage_0_vhvL_interface.pdb"
OUT_SPLIT_PDB = ROOT / "structures" / "domains" / "fab_stage_0_split_mpnn.pdb"
OUT_JSON = ROOT / "structures" / "interface" / "stage_0_vhvL_interface_config.json"
PISA_WORK = ROOT / "structures" / "interface" / "pisa_wt_fv"
PISA_XML = PISA_WORK / "wt_fv_interfaces.xml"


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
    interface_definition: dict | None = None,
    vh_core: list[int] | None = None,
    vl_core: list[int] | None = None,
) -> dict:
    if vh_core is None or vl_core is None:
        vh_core, vl_core = (
            interface_closure_core(vh_len, vl_len, max_heavy_a=core_max_heavy_a)
            if core_max_heavy_a is not None
            else ([], [])
        )

    cfg = {
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
                "holo_fv_interface_framework_rmsd_A",
                "holo_fv_cdr_rmsd_A",
                "cdr_epitope_contacts",
                "vh_vl_interface_clashes",
            ],
            "filters": {
                "max_static_fraction_of_native_contacts": 0.45,
                "max_holo_clashes_4A": 50,
                "max_holo_fv_framework_rmsd_A": 3.5,
                "max_holo_fv_cdr_rmsd_A": 6.0,
                "min_holo_cdr_epitope_contacts": 6,
                "max_holo_vh_vl_interface_clashes": 0,
            },
        },
        "next_stage": (
            "Feed top Stage 0 Fab (chains A–D) into build_stage_a_design_target.py "
            "via --fab-pdb, then run Stage A minibinder design."
        ),
    }
    if interface_definition:
        cfg["interface_definition"] = interface_definition
    return cfg


def resolve_interface_residues(
    fab_pdb: Path,
    vh_len: int,
    vl_len: int,
    *,
    use_pisa: bool,
    min_buried_sasa_A2: float,
) -> tuple[list[int], list[int], dict | None, dict[str, dict[int, float]] | None]:
    """Return (vh_iface_fw, vl_iface_fw, interface_definition, residue_bsa)."""
    legacy_vh = [r for r in VH_INTERFACE_FW if r <= vh_len]
    legacy_vl = [r for r in VL_INTERFACE_FW if r <= vl_len]
    if not use_pisa:
        return legacy_vh, legacy_vl, None, None

    exclude_cdr = {
        "A": cdr_residue_set("A", vh_len, vl_len),
        "B": cdr_residue_set("B", vh_len, vl_len),
    }
    pisa_result = identify_fv_interface_residues(
        fab_pdb,
        work_dir=PISA_WORK,
        vh_len=vh_len,
        vl_len=vl_len,
        min_buried_sasa_A2=min_buried_sasa_A2,
        exclude_cdr_resnums=exclude_cdr,
        session_name="wt_fv_interface",
    )
    if pisa_result.get("pisa_status") != "ok":
        print(
            f"  WARNING: PISA interface detection failed ({pisa_result.get('pisa_error', pisa_result.get('pisa_status'))}); "
            "using legacy contact-based interface lists."
        )
        return legacy_vh, legacy_vl, None, None

    fw = pisa_result["framework_interface_residues"]
    vh_iface = fw.get("A", [])
    vl_iface = fw.get("B", [])
    if not vh_iface or not vl_iface:
        print("  WARNING: PISA returned empty framework interface; using legacy lists.")
        return legacy_vh, legacy_vl, None, None

    PISA_WORK.mkdir(parents=True, exist_ok=True)
    PISA_XML.write_text(pisa_result.pop("pisa_xml_text", ""))

    interface_definition = {
        "source": "pisa",
        "structure": str(fab_pdb.relative_to(ROOT)),
        "fv_chains": ["A", "B"],
        "min_buried_sasa_A2": min_buried_sasa_A2,
        "exclude_cdrs": True,
        "vh_framework_interface": vh_iface,
        "vl_framework_interface": vl_iface,
        "legacy_vh_framework_interface": legacy_vh,
        "legacy_vl_framework_interface": legacy_vl,
        "pisa_metrics": {
            k: pisa_result[k]
            for k in (
                "pisa_int_area_A2",
                "pisa_int_solv_en_kcal",
                "pisa_interface_id",
                "buried_interface_residues",
            )
            if k in pisa_result
        },
        "pisa_xml": str(PISA_XML.relative_to(ROOT)),
    }
    residue_bsa = pisa_result.get("residue_buried_sasa_A2")
    print(f"  PISA WT interface: A={len(vh_iface)} B={len(vl_iface)} framework residues")
    print(f"    A: {', '.join(f'A{r}' for r in vh_iface)}")
    print(f"    B: {', '.join(f'B{r}' for r in vl_iface)}")
    return vh_iface, vl_iface, interface_definition, residue_bsa


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
    p.add_argument(
        "--pisa-interface",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use PISA buried SASA on WT Fv to define VH–VL interface residues (default: on)",
    )
    p.add_argument(
        "--pisa-min-bsa",
        type=float,
        default=0.0,
        help="Minimum PISA buried SASA (Å²) for a residue to count as interface",
    )
    args = p.parse_args()

    source = str(args.source_pdb) if args.source_pdb else "1N8Z (trastuzumab)"
    vh_len, vl_len, ch1_len, cl_len = build_fab_context_pdb(
        args.out_pdb, source_pdb=args.source_pdb, struct_name="stage_0"
    )
    build_split_fv_mpnn_pdb(args.out_split_pdb, source_pdb=args.out_pdb, separation_a=args.separation_A)

    vh_iface, vl_iface, interface_definition, residue_bsa = resolve_interface_residues(
        args.out_pdb,
        vh_len,
        vl_len,
        use_pisa=args.pisa_interface,
        min_buried_sasa_A2=args.pisa_min_bsa,
    )

    if args.conservative:
        core = INTERFACE_CORE_MAX_HEAVY_A
        if interface_definition and residue_bsa:
            vh_core, vl_core = pisa_closure_core(
                residue_bsa, vh_iface, vl_iface, core_min_buried_sasa_A2=PISA_CORE_MIN_BURIED_SASA_A2
            )
        else:
            vh_core, vl_core = interface_closure_core(
                vh_len, vl_len, max_heavy_a=core,
                vh_interface_fw=vh_iface, vl_interface_fw=vl_iface,
            )
        designed = split_mpnn_designed_residues(
            vh_len, vl_len, max_heavy_a=core,
            vh_interface_fw=vh_iface, vl_interface_fw=vl_iface,
            vh_core=vh_core, vl_core=vl_core,
        )
        approach = "split_mpnn_partial_degrease"
        desc = f"Partial de-grease: {len(designed)} rim residues; core fixed."
    elif args.aggressive:
        core = INTERFACE_CORE_MAX_HEAVY_A_AGGRESSIVE
        if interface_definition and residue_bsa:
            vh_core, vl_core = pisa_closure_core(
                residue_bsa, vh_iface, vl_iface,
                core_min_buried_sasa_A2=PISA_CORE_MIN_BURIED_SASA_A2_AGGRESSIVE,
            )
        else:
            vh_core, vl_core = interface_closure_core(
                vh_len, vl_len, max_heavy_a=core,
                vh_interface_fw=vh_iface, vl_interface_fw=vl_iface,
            )
        designed = split_mpnn_designed_residues(
            vh_len, vl_len, max_heavy_a=core,
            vh_interface_fw=vh_iface, vl_interface_fw=vl_iface,
            vh_core=vh_core, vl_core=vl_core,
        )
        approach = "split_mpnn_aggressive_degrease"
        desc = f"Aggressive de-grease: {len(designed)} rim residues; core fixed."
    else:
        core = None
        vh_core, vl_core = [], []
        designed = interface_design_residues(vh_len, vl_len, vh_iface, vl_iface)
        approach = "split_mpnn_whole_interface"
        desc = (
            f"Whole-interface: all {len(designed)} PISA/framework interface residues redesigned; "
            f"top {args.top_n_boltz} by static VH–VL weakening → holo Boltz."
        )

    config = build_config(
        vh_len, vl_len, ch1_len, cl_len, source, args.separation_A,
        designed, approach, desc, args.n_mpnn_seqs, args.top_n_boltz, core,
        interface_definition=interface_definition,
        vh_core=vh_core, vl_core=vl_core,
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
