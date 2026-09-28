#!/usr/bin/env python3
"""
analyze_stage_a_minibinders.py
------------------------------
QC and ranking for Stage A grafted minibinder designs.

Answers the two questions that matter for the hidden-switch minibinder:

1. Is the Fab still exactly where it was placed?
   ``fab_max_dev_a`` is the largest CA deviation of chains A–D/T from the input.
   Grafted files must report 0.000.

2. Is the minibinder alpha-helical, and does it actually cross-link the two arms?
   ``helix_frac`` / ``strand_frac`` / ``loop_frac`` come from a CA-trace P-SEA
   assignment; ``ch1_contacts`` and ``vl_contacts`` count minibinder residues
   within ``--contact-cutoff`` of the CH1 and VL hotspot surfaces.

Designs are ranked by (bridges both arms, helix fraction, total hotspot contacts).

Usage
-----
    python scripts/analyze_stage_a_minibinders.py \\
        --input pipeline_results/<run>/inputs/fab_hidden_minibinder_stage_a.pdb \\
        --grafted-dir pipeline_results/<run>/structures/grafted \\
        --out-csv pipeline_results/<run>/final/stage_a_minibinder_qc.csv
"""

from __future__ import annotations

import argparse
import csv
import statistics as stats
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from binary_antibodies.backbone_ss import summarize  # noqa: E402
from binary_antibodies.fab_hidden_switch import (  # noqa: E402
    CH1_HOTSPOTS_HEAVY,
    VH_END,
    VL_HOTSPOTS_STAGE_A,
    _load_biopython_structure,
)

FAB_CHAINS = ("A", "B", "C", "D", "T")
MB_CHAIN = "M"


def chain_ca(model, chain_id: str) -> np.ndarray:
    if chain_id not in model.child_dict:
        return np.empty((0, 3))
    return np.array(
        [res["CA"].get_coord() for res in model[chain_id] if res.id[0] == " " and "CA" in res]
    )


def hotspot_ca(model, chain_id: str, resnums: list[int]) -> np.ndarray:
    by_num = {
        res.id[1]: res["CA"].get_coord()
        for res in model[chain_id]
        if res.id[0] == " " and "CA" in res
    }
    return np.array([by_num[n] for n in resnums if n in by_num])


def heavy_atoms(model, chain_ids: tuple[str, ...]) -> np.ndarray:
    coords = [
        atom.get_coord()
        for cid in chain_ids
        if cid in model.child_dict
        for res in model[cid]
        if res.id[0] == " "
        for atom in res
        if atom.element != "H"
    ]
    return np.array(coords)


def analyze(path: Path, ref_model, ref_ca, ch1_hs, vl_hs, cutoff: float) -> dict[str, object]:
    model = list(_load_biopython_structure(path).get_models())[0]

    fab_dev = 0.0
    for cid in FAB_CHAINS:
        out, ref = chain_ca(model, cid), ref_ca[cid]
        if len(ref) == 0:
            continue
        if len(out) != len(ref):
            fab_dev = float("nan")
            break
        fab_dev = max(fab_dev, float(np.abs(out - ref).max()))

    mb_ca = chain_ca(model, MB_CHAIN)
    if len(mb_ca) < 5:
        raise ValueError(f"{path.name}: minibinder chain {MB_CHAIN} missing or too short")
    ss = summarize(mb_ca)

    def contacts(hotspots: np.ndarray) -> tuple[int, float]:
        if len(hotspots) == 0:
            return 0, float("nan")
        d = np.linalg.norm(mb_ca[:, None, :] - hotspots[None, :, :], axis=2)
        return int((d.min(1) <= cutoff).sum()), float(d.min())

    ch1_contacts, ch1_min = contacts(ch1_hs)
    vl_contacts, vl_min = contacts(vl_hs)

    mb_heavy = heavy_atoms(model, (MB_CHAIN,))
    fab_heavy = heavy_atoms(ref_model, FAB_CHAINS)
    clash = float(np.linalg.norm(mb_heavy[:, None, :] - fab_heavy[None, :, :], axis=2).min())

    centred = mb_ca - mb_ca.mean(0)
    rg = float(np.sqrt((centred**2).sum(1).mean()))
    span = float(np.linalg.norm(mb_ca[:, None, :] - mb_ca[None, :, :], axis=2).max())

    return {
        "design": path.stem.replace("_grafted", ""),
        "mb_length": ss["length"],
        "helix_frac": round(float(ss["helix_frac"]), 3),
        "strand_frac": round(float(ss["strand_frac"]), 3),
        "loop_frac": round(float(ss["loop_frac"]), 3),
        "n_helices": ss["n_helices"],
        "longest_helix": ss["longest_helix"],
        "ch1_contacts": ch1_contacts,
        "vl_contacts": vl_contacts,
        "ch1_min_dist_a": round(ch1_min, 2),
        "vl_min_dist_a": round(vl_min, 2),
        "bridges_both": int(ch1_contacts > 0 and vl_contacts > 0),
        "min_fab_clash_a": round(clash, 2),
        "radius_gyration_a": round(rg, 2),
        "max_span_a": round(span, 2),
        "fab_max_dev_a": round(fab_dev, 4) if fab_dev == fab_dev else "NA",
        "ss_string": ss["ss"],
    }


def main() -> None:
    p = argparse.ArgumentParser(description="QC Stage A grafted minibinder designs.")
    p.add_argument("--input", type=Path, required=True, help="Stage A input Fab PDB (A–D,T)")
    p.add_argument("--grafted-dir", type=Path, required=True, help="Directory of *_grafted.pdb")
    p.add_argument("--out-csv", type=Path, default=None)
    p.add_argument(
        "--contact-cutoff",
        type=float,
        default=10.0,
        help="CA–CA cutoff (Å) counting a minibinder residue as touching a hotspot surface",
    )
    p.add_argument("--top", type=int, default=15, help="How many top designs to print")
    args = p.parse_args()

    pdbs = sorted(args.grafted_dir.glob("*_grafted.pdb"))
    if not pdbs:
        sys.exit(f"No *_grafted.pdb in {args.grafted_dir}")

    ref_model = list(_load_biopython_structure(args.input).get_models())[0]
    ref_ca = {cid: chain_ca(ref_model, cid) for cid in FAB_CHAINS}
    ch1_hs = hotspot_ca(ref_model, "C", [r - VH_END for r in CH1_HOTSPOTS_HEAVY if r > VH_END])
    vl_hs = hotspot_ca(ref_model, "B", VL_HOTSPOTS_STAGE_A)

    rows: list[dict[str, object]] = []
    for pdb in pdbs:
        try:
            rows.append(analyze(pdb, ref_model, ref_ca, ch1_hs, vl_hs, args.contact_cutoff))
        except Exception as exc:
            print(f"  SKIP {pdb.name}: {exc}")

    if not rows:
        sys.exit("No designs could be analysed")

    rows.sort(
        key=lambda r: (
            -int(r["bridges_both"]),
            -float(r["helix_frac"]),
            -(int(r["ch1_contacts"]) + int(r["vl_contacts"])),
        )
    )

    helix = [float(r["helix_frac"]) for r in rows]
    strand = [float(r["strand_frac"]) for r in rows]
    loop = [float(r["loop_frac"]) for r in rows]
    devs = [float(r["fab_max_dev_a"]) for r in rows if r["fab_max_dev_a"] != "NA"]

    print(f"Analysed {len(rows)} designs from {args.grafted_dir}")
    print("\nFab fidelity (chains A–D,T vs input):")
    if devs:
        print(f"  max CA deviation across all designs: {max(devs):.4f} Å")
        print(f"  designs with a perfectly preserved Fab: {sum(1 for d in devs if d < 1e-3)}/{len(rows)}")
    else:
        print("  could not be evaluated (chain length mismatch)")

    print("\nMinibinder secondary structure:")
    print(f"  helix  frac: mean {stats.mean(helix):.3f}  median {stats.median(helix):.3f}  max {max(helix):.3f}")
    print(f"  strand frac: mean {stats.mean(strand):.3f}  median {stats.median(strand):.3f}")
    print(f"  loop   frac: mean {stats.mean(loop):.3f}  median {stats.median(loop):.3f}")
    for threshold in (0.3, 0.5, 0.7):
        n = sum(1 for h in helix if h >= threshold)
        print(f"  designs with helix_frac >= {threshold:.1f}: {n}/{len(rows)} ({100 * n / len(rows):.0f}%)")

    bridging = sum(int(r["bridges_both"]) for r in rows)
    print(f"\nCross-linking (CA–CA <= {args.contact_cutoff:.0f} Å to both hotspot surfaces):")
    print(f"  designs bridging CH1 and VL: {bridging}/{len(rows)}")
    helical_bridges = [r for r in rows if r["bridges_both"] and float(r["helix_frac"]) >= 0.5]
    print(f"  of which >= 50% helical: {len(helical_bridges)}")

    by_length: dict[int, list[float]] = {}
    for r in rows:
        by_length.setdefault(int(r["mb_length"]) // 10 * 10, []).append(float(r["helix_frac"]))
    print("\nHelix fraction by minibinder length decade:")
    for decade in sorted(by_length):
        vals = by_length[decade]
        print(f"  {decade}-{decade + 9} aa: n={len(vals):3d}  mean helix {stats.mean(vals):.3f}")

    print(f"\nTop {min(args.top, len(rows))} designs (bridging, then helicity):")
    header = f"  {'design':28s} {'len':>4s} {'H':>5s} {'E':>5s} {'nH':>3s} {'maxH':>5s} {'CH1':>4s} {'VL':>4s} {'clash':>6s}"
    print(header)
    for r in rows[: args.top]:
        print(
            f"  {str(r['design']):28s} {r['mb_length']:4d} {float(r['helix_frac']):5.2f} "
            f"{float(r['strand_frac']):5.2f} {r['n_helices']:3d} {r['longest_helix']:5d} "
            f"{r['ch1_contacts']:4d} {r['vl_contacts']:4d} {float(r['min_fab_clash_a']):6.2f}"
        )

    if args.out_csv:
        args.out_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.out_csv.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nWrote {args.out_csv}")


if __name__ == "__main__":
    main()
