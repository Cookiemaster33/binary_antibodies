#!/usr/bin/env python3
"""
optimise_linker.py
------------------
CLI tool: find the optimal flexible linker length for a conditional
proximity-gated nanobody construct, given the antigen–target geometry.

Usage
-----
    python scripts/optimise_linker.py \
        --distance 8.0 \
        --kd-nanobody 10e-9 \
        --kd-antibody 1e-9 \
        --target-occupancy 0.5 \
        --model fjc
"""

import argparse
import sys
import os

# Allow running from repo root without installing the package
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np
from binary_antibodies.polymer import LinkerModel
from binary_antibodies.design import ConditionalConstruct
from binary_antibodies.sequences import LinkerSequence, recommended_linkers


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Optimise flexible linker length for a conditional nanobody construct.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--distance",
        type=float,
        required=True,
        help="Antigen–target centre-to-centre distance on the membrane (nm).",
    )
    p.add_argument(
        "--kd-nanobody",
        type=float,
        default=10e-9,
        help="Intrinsic Kd of the nanobody for its target (M), e.g. 10e-9 for 10 nM.",
    )
    p.add_argument(
        "--kd-antibody",
        type=float,
        default=1e-9,
        help="Intrinsic Kd of the antibody for the antigen (M).",
    )
    p.add_argument(
        "--target-occupancy",
        type=float,
        default=0.5,
        help="Minimum desired nanobody occupancy when the antibody is anchored (0–1).",
    )
    p.add_argument(
        "--model",
        choices=["fjc", "wlc"],
        default="fjc",
        help="Polymer model for the linker (fjc = freely-jointed chain, wlc = worm-like chain).",
    )
    p.add_argument(
        "--motif",
        type=str,
        default="GGGGS",
        help="Linker repeat motif (e.g. GGGGS for G4S).",
    )
    p.add_argument(
        "--scan",
        action="store_true",
        help="Print a table scanning occupancy vs linker length.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    print(
        f"\nConditional Nanobody Linker Optimiser"
        f"\n{'='*45}"
        f"\n  Antigen–target distance : {args.distance:.1f} nm"
        f"\n  Nanobody Kd (intrinsic) : {args.kd_nanobody*1e9:.2f} nM"
        f"\n  Antibody Kd (antigen)   : {args.kd_antibody*1e9:.2f} nM"
        f"\n  Target occupancy        : {args.target_occupancy:.0%}"
        f"\n  Polymer model           : {args.model.upper()}"
        f"\n  Linker motif            : ({args.motif})n"
    )

    # ── Optimisation ────────────────────────────────────────────────────
    construct = ConditionalConstruct(
        antibody_kd_M=args.kd_antibody,
        nanobody_kd_M=args.kd_nanobody,
        linker=LinkerModel(n_residues=1, model=args.model),  # placeholder
    )

    opt_n = construct.optimal_linker_length(
        distance_nm=args.distance,
        target_occupancy=args.target_occupancy,
        model=args.model,
    )

    print(f"\n{'─'*45}")
    if opt_n is not None:
        lm_opt = LinkerModel(n_residues=opt_n, model=args.model)
        c_eff = lm_opt.effective_concentration_M(args.distance)
        occ = lm_opt.occupancy(args.kd_nanobody, args.distance)
        print(f"  Optimal linker length   : {opt_n} residues")
        print(f"  Contour length          : {lm_opt.contour_length_nm:.1f} nm")
        print(f"  RMS end-to-end          : {lm_opt.rms_end_to_end_nm:.1f} nm")
        print(f"  Effective [nanobody]    : {c_eff*1e6:.3f} µM")
        print(f"  Nanobody occupancy      : {occ:.3f}  ({occ*100:.1f}%)")

        # Recommend the nearest (G4S)n linker
        recs = recommended_linkers(
            distance_nm=args.distance,
            nanobody_kd_M=args.kd_nanobody,
            target_occupancy=args.target_occupancy,
            motif=args.motif,
        )
        if recs:
            best = recs[0]
            ls = LinkerSequence(n_repeats=best["n_repeats"], motif=args.motif)
            print(f"\n  Recommended ({args.motif})n linker:")
            print(f"    n_repeats  : {best['n_repeats']}  →  {best['n_residues']} residues")
            print(f"    Sequence   : {best['sequence']}")
            print(f"    Contour Lc : {best['contour_length_nm']:.1f} nm")
            print(f"    Occupancy  : {best['occupancy']:.3f}  ({best['occupancy']*100:.1f}%)")
            print(f"    Molar mass : {best['molecular_weight_Da']:.0f} Da")
            print(f"\n  FASTA:\n{ls.as_fasta(header='flexible_linker')}")
    else:
        print(
            f"  ✗ Could not achieve {args.target_occupancy:.0%} occupancy within 200 residues"
            f" at {args.distance:.1f} nm distance."
        )
        print(
            f"    Consider: (1) using a nanobody with lower Kd, "
            f"(2) engineering a shorter antigen–target distance, "
            f"or (3) using a rigid/structured scaffold."
        )

    # ── Optional scan table ─────────────────────────────────────────────
    if args.scan:
        print(f"\n{'─'*45}")
        print(f"  Scan: occupancy vs linker length  (distance = {args.distance:.1f} nm)")
        print(f"  {'Residues':>10}  {'Contour(nm)':>12}  {'Ceff(µM)':>10}  {'Occupancy':>10}")
        print(f"  {'─'*10}  {'─'*12}  {'─'*10}  {'─'*10}")
        for n in range(5, 155, 5):
            lm = LinkerModel(n_residues=n, model=args.model)
            c = lm.effective_concentration_M(args.distance)
            o = lm.occupancy(args.kd_nanobody, args.distance)
            marker = " ◀" if opt_n and n == opt_n else ""
            print(
                f"  {n:>10}  {lm.contour_length_nm:>12.1f}  "
                f"{c*1e6:>10.4f}  {o:>10.4f}{marker}"
            )

    print()


if __name__ == "__main__":
    main()
