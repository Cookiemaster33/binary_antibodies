#!/usr/bin/env python3
"""
build_boltz_stage_0_inputs.py
-----------------------------
Rank MPNN sequences by predicted static VH–VL weakening, write holo-only Boltz YAMLs.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from binary_antibodies.boltz_fab_chains import (  # noqa: E402
    boltz_holo_chain_sequences,
    format_boltz_yaml,
    load_boltz_chain_mode,
)
from binary_antibodies.stage_0_scoring import (  # noqa: E402
    NativeFvReference,
    load_interface_fw_lists,
    static_chain_pair_interface_contacts_heavy,
    static_interface_contacts_heavy,
)


def main() -> None:
    p = argparse.ArgumentParser(description="Build holo-only Boltz inputs for Stage 0.")
    p.add_argument("--pipeline-dir", type=Path, default=Path("/workspace"))
    p.add_argument("--reference-pdb", type=Path, default=None)
    p.add_argument("--config-json", type=Path, default=None)
    p.add_argument("--top-n", type=int, default=None)
    args = p.parse_args()

    pipeline = args.pipeline_dir
    config_path = args.config_json or (
        pipeline / "inputs" / os.environ.get("CONFIG_JSON", "stage_0_vhvL_interface_config.json")
    )
    ref_path = args.reference_pdb or (pipeline / "inputs" / "fab_stage_0_vhvL_interface.pdb")
    mpnn_path = pipeline / "outputs" / "mpnn_stage_0" / "all_sequences.json"
    holo_dir = pipeline / "boltz_inputs_holo"
    holo_dir.mkdir(parents=True, exist_ok=True)

    cfg = json.loads(config_path.read_text())
    native = cfg.get("native_chain_sequences", {})
    epitope = native.get("T", "")
    top_n = args.top_n or int(os.environ.get("TOP_N", cfg.get("split_mpnn", {}).get("top_n_boltz", 100)))

    vh_iface, vl_iface, ch1_iface, cl_iface, interface_scope = load_interface_fw_lists(config_path)
    ref = NativeFvReference(
        ref_path, vh_iface, vl_iface, ch1_iface, cl_iface, interface_scope=interface_scope
    )
    rank_key = cfg.get("split_mpnn", {}).get("mpnn_rank_by", "static_vh_vl_interface_contacts_asc")
    boltz_mode = load_boltz_chain_mode(cfg, interface_scope)
    data = json.loads(mpnn_path.read_text())
    if isinstance(data, dict):
        data = data.get("sequences", [])

    for r in data:
        chains = r.get("chains", {})
        seq_a = chains.get("A", "")
        seq_b = chains.get("B", "")
        seq_c = chains.get("C", "")
        seq_d = chains.get("D", "")
        static_c = static_interface_contacts_heavy(
            ref.aa,
            seq_a,
            seq_b,
            ref.native_seq_a,
            ref.native_seq_b,
            vh_iface=ref.vh_iface_fw,
            vl_iface=ref.vl_iface_fw,
        )
        r["static_vh_vl_interface_contacts"] = round(static_c, 2)
        r["static_fraction_of_native"] = round(static_c / ref.native_static_contacts, 4)
        if interface_scope == "full_fab" and ref.ch1_iface_fw and ref.cl_iface_fw and seq_c and seq_d:
            static_cc = static_chain_pair_interface_contacts_heavy(
                ref.aa,
                "C",
                "D",
                seq_c,
                seq_d,
                ref.ch1_iface_fw,
                ref.cl_iface_fw,
                ref.native_seq_c,
                ref.native_seq_d,
            )
            r["static_ch1_cl_interface_contacts"] = round(static_cc, 2)
            r["static_ch1_cl_fraction_of_native"] = round(
                static_cc / ref.native_static_ch1_cl_contacts, 4
            )
            r["static_total_interface_contacts"] = round(static_c + static_cc, 2)
        else:
            r["static_total_interface_contacts"] = round(static_c, 2)

    sort_field = (
        "static_total_interface_contacts"
        if rank_key == "static_total_interface_contacts_asc"
        else "static_vh_vl_interface_contacts"
    )
    data.sort(key=lambda x: x[sort_field])
    picks = data[:top_n]

    design_list = []
    for i, r in enumerate(picks):
        chains = dict(r["chains"])
        for ch, seq in native.items():
            if ch not in chains and seq:
                chains[ch] = seq
        rank = i + 1
        name = f"rank{rank:03d}_s0_{r['backbone']}_s{r['seq_idx']}"
        holo_entries = boltz_holo_chain_sequences(
            chains, mode=boltz_mode, epitope=epitope
        )
        (holo_dir / f"{name}.yaml").write_text(format_boltz_yaml(holo_entries))
        design_list.append(
            {
                "rank": rank,
                "name": name,
                "backbone": r["backbone"],
                "seq_idx": r["seq_idx"],
                "mpnn_score": r.get("mpnn_score"),
                "static_vh_vl_interface_contacts": r["static_vh_vl_interface_contacts"],
                "static_fraction_of_native": r["static_fraction_of_native"],
                "static_total_interface_contacts": r.get("static_total_interface_contacts"),
                "static_ch1_cl_interface_contacts": r.get("static_ch1_cl_interface_contacts"),
                "chains": {
                    ch: chains.get(ch, "")
                    for ch in ("A", "B", "C", "D")
                    if chains.get(ch)
                },
            }
        )

    out = {
        "designs": design_list,
        "n_mpnn_total": len(data),
        "top_n": top_n,
        "rank_by": rank_key,
        "boltz_chain_mode": boltz_mode,
        "native_static_interface_contacts": round(ref.native_static_contacts, 2),
        "native_static_ch1_cl_contacts": round(ref.native_static_ch1_cl_contacts, 2),
    }
    final_dir = pipeline / "final"
    final_dir.mkdir(parents=True, exist_ok=True)
    (final_dir / "top_designs_stage_0.json").write_text(json.dumps(out, indent=2) + "\n")

    if picks:
        lo, hi = picks[0][sort_field], picks[-1][sort_field]
        print(
            f"Boltz holo inputs: {len(design_list)} designs "
            f"(weakest {sort_field} {lo:.1f}–{hi:.1f} of {len(data)} MPNN)"
        )
    else:
        print("No MPNN sequences to process.")


if __name__ == "__main__":
    main()
