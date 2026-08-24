#!/usr/bin/env python3
"""
graft_stage_a_outputs.py
------------------------
Graft RFd3 minibinder designs onto the exact input Fab (A,B,C,D,T unchanged).

Final structures: input Fab + chain M (minibinder).

Usage
-----
    python scripts/graft_stage_a_outputs.py \\
        --input pipeline_results/stage_a_rank079_split_cif_v4/inputs/fab_hidden_minibinder_stage_a.pdb \\
        --rfd3-dir pipeline_results/stage_a_rank079_split_cif_v4/structures/rfd3_stage_a \\
        --out-dir pipeline_results/stage_a_rank079_split_cif_v4/structures/grafted
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from binary_antibodies.fab_hidden_switch import graft_stage_a_minibinder  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser(description="Graft RFd3 MB onto fixed input Fab.")
    p.add_argument("--input", type=Path, required=True, help="Stage A input Fab PDB")
    p.add_argument("--rfd3-dir", type=Path, required=True, help="Directory of raw RFd3 CIFs")
    p.add_argument("--out-dir", type=Path, required=True, help="Output directory for grafted PDBs")
    args = p.parse_args()

    cifs = sorted(args.rfd3_dir.glob("sa_*.cif"))
    if not cifs:
        sys.exit(f"No RFd3 CIFs in {args.rfd3_dir}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    ok = 0
    for cif in cifs:
        out = args.out_dir / f"{cif.stem}_grafted.pdb"
        try:
            mb_len = graft_stage_a_minibinder(out, args.input, cif)
            ok += 1
            if ok <= 3:
                print(f"  {cif.name} → {out.name}  (MB={mb_len} aa)")
        except Exception as exc:
            print(f"  SKIP {cif.name}: {exc}")

    print(f"Grafted {ok}/{len(cifs)} → {args.out_dir}")


if __name__ == "__main__":
    main()
