#!/usr/bin/env python3
"""
build_stage_0_design_target.py
------------------------------
Stage 0: split-chain ProteinMPNN on Fab interface framework residues (no RFd3).

Interface scope
---------------
  fv (default)     — VH–VL only (chains A+B)
  full_fab         — VH–VL (A+B) and CH1–CL (C+D)

Modes
-----
  --whole-interface   redesign all PISA interface framework residues (default)
  --conservative      partial de-grease (PISA or distance core + rim)
  --aggressive        aggressive partial de-grease

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
    build_split_fab_mpnn_pdb,
    build_split_fv_mpnn_pdb,
    cdr_residue_set,
    extract_chain_sequences,
    interface_closure_core,
    interface_design_residues,
    pisa_closure_core,
    pisa_closure_core_pair,
    split_mpnn_designed_residues,
    split_mpnn_designed_residues_fab,
    write_json,
)
from binary_antibodies.pisa_scoring import (  # noqa: E402
    INTERFACE_SCOPE_FULL_FAB,
    INTERFACE_SCOPE_FV,
    identify_fab_interface_residues,
)

OUT_PDB = ROOT / "structures" / "domains" / "fab_stage_0_vhvL_interface.pdb"
OUT_SPLIT_PDB = ROOT / "structures" / "domains" / "fab_stage_0_split_mpnn.pdb"
OUT_JSON = ROOT / "structures" / "interface" / "stage_0_vhvL_interface_config.json"
PISA_WORK = ROOT / "structures" / "interface" / "pisa_wt_fv"
PISA_XML = PISA_WORK / "wt_fab_interfaces.xml"


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
    interface_scope: str,
    interface_definition: dict | None = None,
    vh_core: list[int] | None = None,
    vl_core: list[int] | None = None,
    ch1_core: list[int] | None = None,
    cl_core: list[int] | None = None,
) -> dict:
    if vh_core is None or vl_core is None:
        vh_core, vl_core = (
            interface_closure_core(vh_len, vl_len, max_heavy_a=core_max_heavy_a)
            if core_max_heavy_a is not None
            else ([], [])
        )
    ch1_core = ch1_core or []
    cl_core = cl_core or []

    rank_by = (
        "static_total_interface_contacts_asc"
        if interface_scope == INTERFACE_SCOPE_FULL_FAB
        else "static_vh_vl_interface_contacts_asc"
    )
    metrics = [
        "static_vh_vl_interface_contacts",
        "static_fraction_of_native_contacts",
        "holo_fv_framework_rmsd_A",
        "holo_fv_interface_framework_rmsd_A",
        "holo_fv_cdr_rmsd_A",
        "cdr_epitope_contacts",
        "vh_vl_interface_clashes",
    ]
    filters = {
        "max_static_fraction_of_native_contacts": 0.45,
        "max_holo_clashes_4A": 50,
        "max_holo_fv_framework_rmsd_A": 3.5,
        "max_holo_fv_cdr_rmsd_A": 6.0,
        "min_holo_cdr_epitope_contacts": 6,
        "max_holo_vh_vl_interface_clashes": 0,
    }
    static_desc = "Weighted VH–VL contacts on native Fab geometry from designed A/B sequence"
    if interface_scope == INTERFACE_SCOPE_FULL_FAB:
        metrics.extend(
            [
                "static_ch1_cl_interface_contacts",
                "static_ch1_cl_fraction_of_native_contacts",
                "ch1_cl_interface_clashes",
            ]
        )
        filters["max_static_ch1_cl_fraction_of_native_contacts"] = 0.45
        filters["max_holo_ch1_cl_interface_clashes"] = 0
        static_desc += "; CH1–CL contacts from designed C/D sequence"

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
            "D": f"CL 1-{cl_len} (steric context; redesigned in full_fab scope)",
            "T": f"HER2 epitope stub 1-{len(EPITOPE_SEQ)} (holo closure context)",
        },
        "split_mpnn": {
            "separation_A": separation_a,
            "core_max_heavy_A": core_max_heavy_a,
            "designed_residues": designed,
            "vh_closure_core": [f"A{r}" for r in vh_core],
            "vl_closure_core": [f"B{r}" for r in vl_core],
            "ch1_closure_core": [f"C{r}" for r in ch1_core],
            "cl_closure_core": [f"D{r}" for r in cl_core],
            "n_sequences": n_mpnn_seqs,
            "top_n_boltz": top_n_boltz,
            "mpnn_rank_by": rank_by,
        },
        "native_chain_sequences": {},
        "validation": {
            "mode": "holo_only_static_interface",
            "static_interface": static_desc,
            "holo_fold": "Boltz A+B+C+D+T — epitope engagement and Fab-like geometry",
            "metrics": metrics,
            "filters": filters,
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
    ch1_len: int,
    cl_len: int,
    *,
    interface_scope: str,
    use_pisa: bool,
    min_buried_sasa_A2: float,
) -> tuple[
    list[int],
    list[int],
    list[int],
    list[int],
    dict | None,
    dict[str, dict[int, float]] | None,
]:
    """Return (vh, vl, ch1, cl iface lists, interface_definition, residue_bsa)."""
    legacy_vh = [r for r in VH_INTERFACE_FW if r <= vh_len]
    legacy_vl = [r for r in VL_INTERFACE_FW if r <= vl_len]
    legacy_ch1: list[int] = []
    legacy_cl: list[int] = []

    if not use_pisa:
        return legacy_vh, legacy_vl, legacy_ch1, legacy_cl, None, None

    exclude_cdr = {
        "A": cdr_residue_set("A", vh_len, vl_len),
        "B": cdr_residue_set("B", vh_len, vl_len),
    }
    pisa_result = identify_fab_interface_residues(
        fab_pdb,
        work_dir=PISA_WORK,
        interface_scope=interface_scope,
        chain_lengths={"A": vh_len, "B": vl_len, "C": ch1_len, "D": cl_len},
        min_buried_sasa_A2=min_buried_sasa_A2,
        exclude_cdr_resnums=exclude_cdr,
        session_name=f"wt_{interface_scope}_interface",
    )
    if pisa_result.get("pisa_status") != "ok":
        print(
            f"  WARNING: PISA interface detection failed ({pisa_result.get('pisa_error', pisa_result.get('pisa_status'))}); "
            "using legacy VH/VL lists only."
        )
        return legacy_vh, legacy_vl, legacy_ch1, legacy_cl, None, None

    fw = pisa_result["framework_interface_residues"]
    vh_iface = fw.get("A", [])
    vl_iface = fw.get("B", [])
    ch1_iface = fw.get("C", []) if interface_scope == INTERFACE_SCOPE_FULL_FAB else []
    cl_iface = fw.get("D", []) if interface_scope == INTERFACE_SCOPE_FULL_FAB else []
    if not vh_iface or not vl_iface:
        print("  WARNING: PISA returned empty VH/VL framework interface; using legacy lists.")
        return legacy_vh, legacy_vl, legacy_ch1, legacy_cl, None, None
    if interface_scope == INTERFACE_SCOPE_FULL_FAB and (not ch1_iface or not cl_iface):
        print("  WARNING: PISA returned empty CH1/CL interface; falling back to VH/VL scope only.")
        interface_scope = INTERFACE_SCOPE_FV
        ch1_iface, cl_iface = [], []

    PISA_WORK.mkdir(parents=True, exist_ok=True)
    PISA_XML.write_text(pisa_result.pop("pisa_xml_text", ""))

    interface_definition = {
        "source": "pisa",
        "scope": interface_scope,
        "structure": str(fab_pdb.relative_to(ROOT)),
        "min_buried_sasa_A2": min_buried_sasa_A2,
        "exclude_cdrs": True,
        "vh_framework_interface": vh_iface,
        "vl_framework_interface": vl_iface,
        "legacy_vh_framework_interface": legacy_vh,
        "legacy_vl_framework_interface": legacy_vl,
        "pisa_metrics": pisa_result.get("interfaces", {}),
        "pisa_xml": str(PISA_XML.relative_to(ROOT)),
    }
    if interface_scope == INTERFACE_SCOPE_FULL_FAB:
        interface_definition["ch1_framework_interface"] = ch1_iface
        interface_definition["cl_framework_interface"] = cl_iface

    residue_bsa = pisa_result.get("residue_buried_sasa_A2")
    print(f"  PISA WT interface ({interface_scope}): A={len(vh_iface)} B={len(vl_iface)}", end="")
    if interface_scope == INTERFACE_SCOPE_FULL_FAB:
        print(f" C={len(ch1_iface)} D={len(cl_iface)}", end="")
    print(" framework residues")
    print(f"    A: {', '.join(f'A{r}' for r in vh_iface)}")
    print(f"    B: {', '.join(f'B{r}' for r in vl_iface)}")
    if ch1_iface:
        print(f"    C: {', '.join(f'C{r}' for r in ch1_iface)}")
    if cl_iface:
        print(f"    D: {', '.join(f'D{r}' for r in cl_iface)}")
    return vh_iface, vl_iface, ch1_iface, cl_iface, interface_definition, residue_bsa


def main() -> None:
    p = argparse.ArgumentParser(description="Build Stage 0 split-MPNN design targets.")
    p.add_argument("--source-pdb", type=Path, default=None, help="Parent Fab PDB")
    p.add_argument("--out-pdb", type=Path, default=OUT_PDB)
    p.add_argument("--out-split-pdb", type=Path, default=OUT_SPLIT_PDB)
    p.add_argument("--out-json", type=Path, default=OUT_JSON)
    p.add_argument("--separation-A", type=float, default=DEFAULT_SPLIT_SEPARATION_A)
    p.add_argument(
        "--interface-scope",
        choices=[INTERFACE_SCOPE_FV, INTERFACE_SCOPE_FULL_FAB],
        default=INTERFACE_SCOPE_FV,
        help="fv = VH/VL only; full_fab = VH/VL + CH1/CL interfaces",
    )
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
        help="Use PISA buried SASA on WT Fab to define interface residues (default: on)",
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
    if args.interface_scope == INTERFACE_SCOPE_FULL_FAB:
        build_split_fab_mpnn_pdb(args.out_split_pdb, source_pdb=args.out_pdb, separation_a=args.separation_A)
    else:
        build_split_fv_mpnn_pdb(args.out_split_pdb, source_pdb=args.out_pdb, separation_a=args.separation_A)

    vh_iface, vl_iface, ch1_iface, cl_iface, interface_definition, residue_bsa = resolve_interface_residues(
        args.out_pdb,
        vh_len,
        vl_len,
        ch1_len,
        cl_len,
        interface_scope=args.interface_scope,
        use_pisa=args.pisa_interface,
        min_buried_sasa_A2=args.pisa_min_bsa,
    )
    effective_scope = (
        interface_definition.get("scope", args.interface_scope)
        if interface_definition
        else args.interface_scope
    )

    ch1_core: list[int] = []
    cl_core: list[int] = []
    if args.conservative:
        core = INTERFACE_CORE_MAX_HEAVY_A
        if interface_definition and residue_bsa:
            vh_core, vl_core = pisa_closure_core(
                residue_bsa, vh_iface, vl_iface, core_min_buried_sasa_A2=PISA_CORE_MIN_BURIED_SASA_A2
            )
            if effective_scope == INTERFACE_SCOPE_FULL_FAB:
                ch1_core, cl_core = pisa_closure_core_pair(
                    residue_bsa, "C", "D", ch1_iface, cl_iface,
                    core_min_buried_sasa_A2=PISA_CORE_MIN_BURIED_SASA_A2,
                )
        else:
            vh_core, vl_core = interface_closure_core(
                vh_len, vl_len, max_heavy_a=core,
                vh_interface_fw=vh_iface, vl_interface_fw=vl_iface,
            )
        designed = split_mpnn_designed_residues_fab(
            vh_len, vl_len, ch1_len, cl_len, max_heavy_a=core,
            vh_interface_fw=vh_iface, vl_interface_fw=vl_iface,
            ch1_interface_fw=ch1_iface, cl_interface_fw=cl_iface,
            vh_core=vh_core, vl_core=vl_core, ch1_core=ch1_core, cl_core=cl_core,
            interface_scope=effective_scope,
        )
        approach = "split_mpnn_partial_degrease"
        desc = f"Partial de-grease ({effective_scope}): {len(designed)} rim residues; core fixed."
    elif args.aggressive:
        core = INTERFACE_CORE_MAX_HEAVY_A_AGGRESSIVE
        if interface_definition and residue_bsa:
            vh_core, vl_core = pisa_closure_core(
                residue_bsa, vh_iface, vl_iface,
                core_min_buried_sasa_A2=PISA_CORE_MIN_BURIED_SASA_A2_AGGRESSIVE,
            )
            if effective_scope == INTERFACE_SCOPE_FULL_FAB:
                ch1_core, cl_core = pisa_closure_core_pair(
                    residue_bsa, "C", "D", ch1_iface, cl_iface,
                    core_min_buried_sasa_A2=PISA_CORE_MIN_BURIED_SASA_A2_AGGRESSIVE,
                )
        else:
            vh_core, vl_core = interface_closure_core(
                vh_len, vl_len, max_heavy_a=core,
                vh_interface_fw=vh_iface, vl_interface_fw=vl_iface,
            )
        designed = split_mpnn_designed_residues_fab(
            vh_len, vl_len, ch1_len, cl_len, max_heavy_a=core,
            vh_interface_fw=vh_iface, vl_interface_fw=vl_iface,
            ch1_interface_fw=ch1_iface, cl_interface_fw=cl_iface,
            vh_core=vh_core, vl_core=vl_core, ch1_core=ch1_core, cl_core=cl_core,
            interface_scope=effective_scope,
        )
        approach = "split_mpnn_aggressive_degrease"
        desc = f"Aggressive de-grease ({effective_scope}): {len(designed)} rim residues; core fixed."
    else:
        core = None
        vh_core, vl_core = [], []
        designed = interface_design_residues(
            vh_len, vl_len, vh_iface, vl_iface, ch1_iface, cl_iface, effective_scope
        )
        approach = "split_mpnn_whole_interface"
        desc = (
            f"Whole-interface ({effective_scope}): all {len(designed)} PISA framework interface residues redesigned; "
            f"top {args.top_n_boltz} by static weakening → holo Boltz."
        )

    config = build_config(
        vh_len, vl_len, ch1_len, cl_len, source, args.separation_A,
        designed, approach, desc, args.n_mpnn_seqs, args.top_n_boltz, core,
        effective_scope,
        interface_definition=interface_definition,
        vh_core=vh_core, vl_core=vl_core, ch1_core=ch1_core, cl_core=cl_core,
    )
    config["native_chain_sequences"] = extract_chain_sequences(args.out_pdb)
    write_json(args.out_json, config)

    sm = config["split_mpnn"]
    print(f"Wrote Fab context PDB  → {args.out_pdb}")
    print(f"Wrote split MPNN PDB   → {args.out_split_pdb}")
    print(f"  interface scope: {effective_scope}")
    print(f"Wrote config           → {args.out_json}")
    print(f"  approach:     {config['approach']}")
    print(f"  designed:     {len(sm['designed_residues'])} residues")
    print(f"  MPNN / Boltz: {sm['n_sequences']} → top {sm['top_n_boltz']}")


if __name__ == "__main__":
    main()
