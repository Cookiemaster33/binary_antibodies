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

from binary_antibodies.stage_0_scoring import (  # noqa: E402
    NativeFvReference,
    load_interface_fw_lists,
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

    vh_iface, vl_iface = load_interface_fw_lists(config_path)
    ref = NativeFvReference(ref_path, vh_iface, vl_iface)
    data = json.loads(mpnn_path.read_text())
    if isinstance(data, dict):
        data = data.get("sequences", [])

    for r in data:
        chains = r.get("chains", {})
        seq_a = chains.get("A", "")
        seq_b = chains.get("B", "")
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

    data.sort(key=lambda x: x["static_vh_vl_interface_contacts"])
    picks = data[:top_n]

    design_list = []
    for i, r in enumerate(picks):
        chains = dict(r["chains"])
        for ch, seq in native.items():
            if ch not in chains and seq:
                chains[ch] = seq
        rank = i + 1
        name = f"rank{rank:03d}_s0_{r['backbone']}_s{r['seq_idx']}"
        holo_yaml = "sequences:\n"
        for ch in ("A", "B", "C", "D"):
            if chains.get(ch):
                holo_yaml += (
                    f"  - protein:\n      id: {ch}\n      sequence: \"{chains[ch]}\"\n      msa: empty\n"
                )
        if epitope:
            holo_yaml += f"  - protein:\n      id: T\n      sequence: \"{epitope}\"\n      msa: empty\n"
        (holo_dir / f"{name}.yaml").write_text(holo_yaml)
        design_list.append(
            {
                "rank": rank,
                "name": name,
                "backbone": r["backbone"],
                "seq_idx": r["seq_idx"],
                "mpnn_score": r.get("mpnn_score"),
                "static_vh_vl_interface_contacts": r["static_vh_vl_interface_contacts"],
                "static_fraction_of_native": r["static_fraction_of_native"],
                "chains": {"A": chains.get("A", ""), "B": chains.get("B", "")},
            }
        )

    out = {
        "designs": design_list,
        "n_mpnn_total": len(data),
        "top_n": top_n,
        "rank_by": "static_vh_vl_interface_contacts_asc",
        "native_static_interface_contacts": round(ref.native_static_contacts, 2),
    }
    final_dir = pipeline / "final"
    final_dir.mkdir(parents=True, exist_ok=True)
    (final_dir / "top_designs_stage_0.json").write_text(json.dumps(out, indent=2) + "\n")

    if picks:
        lo, hi = picks[0]["static_vh_vl_interface_contacts"], picks[-1]["static_vh_vl_interface_contacts"]
        print(
            f"Boltz holo inputs: {len(design_list)} designs "
            f"(weakest static contacts {lo:.1f}–{hi:.1f} of {len(data)} MPNN)"
        )
    else:
        print("No MPNN sequences to process.")


if __name__ == "__main__":
    main()
