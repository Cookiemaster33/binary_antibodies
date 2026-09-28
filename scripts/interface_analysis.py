"""
interface_analysis.py
---------------------
Visualise and rank the VH–VL interface residues of the Trastuzumab scFv.
Produces plots used to guide minibinder design.

Usage
-----
    python scripts/interface_analysis.py --structures-dir structures
"""

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from binary_antibodies.structures import (
    load_structure, extract_variable_domain,
    analyse_vh_vl_interface, interface_summary,
    VL_FRAMEWORK_CONTACT_POSITIONS,
)


plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150,
})


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Visualise Trastuzumab VH-VL interface for minibinder design.",
    )
    p.add_argument("--structures-dir", default="structures")
    p.add_argument("--output-dir", default="figures")
    return p.parse_args()


def contact_heatmap(
    contacts,
    summary: dict,
    output_path: str,
) -> None:
    """Plot a contact map between VH and VL residues."""
    vh_res = sorted(summary["vh_contact_residues"])
    vl_res = sorted(summary["vl_contact_residues"])

    # Build distance matrix
    dist_matrix = np.full((len(vl_res), len(vh_res)), np.nan)
    vl_idx = {r: i for i, r in enumerate(vl_res)}
    vh_idx = {r: i for i, r in enumerate(vh_res)}

    for c in contacts:
        i = vl_idx.get(c.vl_resnum)
        j = vh_idx.get(c.vh_resnum)
        if i is not None and j is not None:
            curr = dist_matrix[i, j]
            if np.isnan(curr) or c.min_heavy_distance_A < curr:
                dist_matrix[i, j] = c.min_heavy_distance_A

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # ── Left: contact map ─────────────────────────────────────────────
    ax = axes[0]
    masked = np.ma.masked_where(np.isnan(dist_matrix), dist_matrix)
    im = ax.imshow(masked, cmap="YlOrRd_r", aspect="auto", vmin=2, vmax=8)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Min heavy-atom distance (Å)")

    # Label VL framework contact residues
    fw_rows = [vl_idx[r] for r in summary["vl_framework_contacts"] if r in vl_idx]
    for row in fw_rows:
        ax.axhline(row, color="blue", alpha=0.25, lw=4)

    ax.set_yticks(range(len(vl_res)))
    ax.set_yticklabels(vl_res, fontsize=6)
    ax.set_xticks(range(0, len(vh_res), 5))
    ax.set_xticklabels(vh_res[::5], fontsize=6, rotation=90)
    ax.set_xlabel("VH residue number")
    ax.set_ylabel("VL residue number")
    ax.set_title("Trastuzumab VH–VL contact map\n(blue bars = VL framework 2 residues → minibinder target)")

    fw_patch = mpatches.Patch(color="blue", alpha=0.3, label="VL framework 2 (minibinder target)")
    ax.legend(handles=[fw_patch], fontsize=8, loc="upper right")

    # ── Right: VL residue contact frequency ──────────────────────────
    ax = axes[1]
    vl_contact_count = {r: 0 for r in vl_res}
    for c in contacts:
        if c.vl_resnum in vl_contact_count:
            vl_contact_count[c.vl_resnum] += 1

    colors = [
        "#d62728" if r in VL_FRAMEWORK_CONTACT_POSITIONS else "#1f77b4"
        for r in vl_res
    ]
    ax.bar(range(len(vl_res)), [vl_contact_count[r] for r in vl_res], color=colors)
    ax.set_xticks(range(0, len(vl_res), 3))
    ax.set_xticklabels(vl_res[::3], fontsize=7, rotation=90)
    ax.set_xlabel("VL residue number")
    ax.set_ylabel("Number of VH contacts")
    ax.set_title("VL contact frequency\n(red = framework 2 — primary minibinder target)")

    fw_patch = mpatches.Patch(color="#d62728", label="Framework 2 residues")
    other_patch = mpatches.Patch(color="#1f77b4", label="Other contact residues")
    ax.legend(handles=[fw_patch, other_patch], fontsize=8)

    fig.suptitle(
        "Trastuzumab VH–VL Interface Analysis\n"
        "Identifying VL1 residues for minibinder design",
        fontsize=12, y=1.01,
    )
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def vl_sequence_annotation(
    vl_domain,
    contact_resnums: list[int],
    fw_resnums: list[int],
    output_path: str,
) -> None:
    """Plot VL sequence with interface residues annotated."""
    from Bio.SeqUtils import seq1

    seq_data = []
    for res in vl_domain.residues:
        try:
            aa = seq1(res.get_resname())
        except Exception:
            aa = "X"
        resnum = res.get_id()[1]
        seq_data.append((resnum, aa))

    resnums = [r for r, _ in seq_data]
    aas = [a for _, a in seq_data]

    fig, ax = plt.subplots(figsize=(14, 3))
    ax.axis("off")

    cols = 40  # residues per row
    for idx, (resnum, aa) in enumerate(seq_data):
        row = idx // cols
        col = idx % cols
        x = col / cols
        y = 1.0 - row * 0.35

        if resnum in fw_resnums:
            bg = "#ff4444"
            color = "white"
            weight = "bold"
        elif resnum in contact_resnums:
            bg = "#4488ff"
            color = "white"
            weight = "bold"
        else:
            bg = "white"
            color = "black"
            weight = "normal"

        ax.text(x, y, aa,
                ha="center", va="center", fontsize=8,
                fontfamily="monospace", fontweight=weight,
                color=color,
                bbox=dict(boxstyle="square,pad=0.15", fc=bg, ec="gray", lw=0.3))
        ax.text(x, y - 0.12, str(resnum), ha="center", va="center",
                fontsize=4, color="gray")

    fw_patch = mpatches.Patch(color="#ff4444", label="Framework 2 (minibinder target)")
    other_patch = mpatches.Patch(color="#4488ff", label="Other VH contact residues")
    ax.legend(handles=[fw_patch, other_patch], fontsize=8,
              loc="lower right", bbox_to_anchor=(1, -0.1))
    ax.set_title(
        "Trastuzumab VL1 sequence — interface residues annotated\n"
        "Red = minibinder must bind; blue = additional VH contacts",
        fontsize=10,
    )

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    pdb_dir = Path(args.structures_dir)
    print("\nLoading structures ...")
    trast = load_structure(pdb_dir / "1N8Z.pdb", "1N8Z")
    vh = extract_variable_domain(trast, "A", "VH")
    vl = extract_variable_domain(trast, "B", "VL")
    print(f"  VH: {len(vh)} residues | VL: {len(vl)} residues")

    print("\nAnalysing VH–VL interface ...")
    contacts = analyse_vh_vl_interface(vh, vl, contact_threshold_A=8.0)
    summary = interface_summary(contacts)
    print(f"  {summary['n_contacts']} contacts found")
    print(f"  VL framework contacts (minibinder target): {summary['vl_framework_contacts']}")

    print("\nGenerating figures ...")
    contact_heatmap(
        contacts, summary,
        os.path.join(args.output_dir, "vh_vl_contact_map.png"),
    )
    vl_sequence_annotation(
        vl,
        summary["vl_contact_residues"],
        summary["vl_framework_contacts"],
        os.path.join(args.output_dir, "vl_sequence_annotation.png"),
    )
    print()


if __name__ == "__main__":
    main()
