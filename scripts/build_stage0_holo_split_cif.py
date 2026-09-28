#!/usr/bin/env python3
"""
build_stage0_holo_split_cif.py
------------------------------
Expand a fused Boltz holo (chains H/L[/T]) to A/B/C/D[/T] and separate VL/CL
from VH/CH1 for PyMOL inspection and Stage A input.

Usage
-----
    python scripts/build_stage0_holo_split_cif.py \\
        pipeline_results/stage_0_full_fab_fused_t025/structures/holo/rank079_s0_native_split_s296_model_0.cif

    # Custom output path (default: <input_stem>_split.cif next to input)
    python scripts/build_stage0_holo_split_cif.py holo.cif -o structures/top5_holo/rank079_split.cif
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from binary_antibodies.fab_hidden_switch import (  # noqa: E402
    DEFAULT_STAGE_A_SEPARATION_A,
    build_holo_split_cif,
    extract_chain_sequences,
    fab_has_split_chains,
)


def default_split_path(holo_cif: Path) -> Path:
    stem = holo_cif.stem
    if stem.endswith("_split"):
        return holo_cif.with_name(f"{stem}.cif" if holo_cif.suffix else f"{stem}_out.cif")
    return holo_cif.with_name(f"{stem}_split{holo_cif.suffix or '.cif'}")


def main() -> None:
    p = argparse.ArgumentParser(description="Build *_split.cif from fused Stage 0 holo.")
    p.add_argument("holo_cif", type=Path, help="Fused Boltz holo (H/L[/T]) or existing split CIF")
    p.add_argument(
        "-o",
        "--out",
        type=Path,
        default=None,
        help="Output path (default: <stem>_split.cif beside input)",
    )
    p.add_argument(
        "--separation-a",
        type=float,
        default=DEFAULT_STAGE_A_SEPARATION_A,
        help=f"Å translation of VL/CL away from VH/CH1 (default {DEFAULT_STAGE_A_SEPARATION_A})",
    )
    args = p.parse_args()

    holo = args.holo_cif.resolve()
    if not holo.exists():
        sys.exit(f"ERROR: holo not found: {holo}")

    out = (args.out or default_split_path(holo)).resolve()
    if fab_has_split_chains(holo) and holo == out:
        sys.exit(f"ERROR: input already has split chains A–D: {holo}")

    vh_len, vl_len, ch1_len, cl_len = build_holo_split_cif(
        out,
        holo,
        separation_a=args.separation_a,
        struct_name=out.stem,
    )
    seqs = extract_chain_sequences(out)
    chains = sorted(seqs)

    print(f"Wrote split Fab → {out}")
    print(f"  chains: {', '.join(chains)}")
    print(f"  lengths: VH={vh_len} VL={vl_len} CH1={ch1_len} CL={cl_len}")
    print(f"  separation: {args.separation_a} Å (VL/CL away from VH/CH1)")
    print()
    print("Stage A:")
    print(f"  python scripts/build_stage_a_design_target.py --fab-pdb {out}")


if __name__ == "__main__":
    main()
