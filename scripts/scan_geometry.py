#!/usr/bin/env python3
"""
scan_geometry.py
----------------
Generate and save a heatmap of nanobody occupancy over the 2-D parameter
space of (antigen–target distance) × (linker length), and a set of
1-D profile plots.

Output files are saved to ./figures/.

Usage
-----
    python scripts/scan_geometry.py \
        --kd-nanobody 10e-9 \
        --model fjc \
        --output-dir figures
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from binary_antibodies.design import ConditionalConstruct
from binary_antibodies.polymer import LinkerModel


# ── Styling ────────────────────────────────────────────────────────────
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 150,
    }
)

PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Scan geometry for conditional nanobody construct design.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--kd-nanobody", type=float, default=10e-9,
                   help="Nanobody intrinsic Kd (M).")
    p.add_argument("--kd-antibody", type=float, default=1e-9,
                   help="Antibody Kd for antigen (M).")
    p.add_argument("--model", choices=["fjc", "wlc"], default="fjc")
    p.add_argument("--output-dir", type=str, default="figures")
    return p.parse_args()


def heatmap_plot(data: dict, kd_nM: float, output_path: str) -> None:
    """Heatmap of nanobody occupancy over (distance × linker length)."""
    occ = data["occupancy"]
    distances = data["distances_nm"]
    lengths = data["linker_lengths"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # ── Left: occupancy heatmap ────────────────────────────────────────
    ax = axes[0]
    im = ax.imshow(
        occ,
        origin="lower",
        aspect="auto",
        cmap="RdYlGn",
        vmin=0,
        vmax=1,
        extent=[lengths[0], lengths[-1], distances[0], distances[-1]],
    )
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Nanobody occupancy (anchored)", fontsize=10)
    # 50% contour line
    ax.contour(
        lengths,
        distances,
        occ,
        levels=[0.5],
        colors="white",
        linewidths=2,
        linestyles="--",
    )
    ax.set_xlabel("Linker length (residues)")
    ax.set_ylabel("Antigen–target distance (nm)")
    ax.set_title(
        f"Conditional nanobody occupancy\n"
        f"(Kd = {kd_nM:.0f} nM, model = FJC)",
        fontsize=11,
    )
    ax.text(
        0.98, 0.98,
        "-- 50% occupancy contour",
        transform=ax.transAxes,
        ha="right", va="top",
        color="white", fontsize=9,
    )

    # ── Right: effective concentration heatmap ─────────────────────────
    ax = axes[1]
    c_eff_uM = data["c_eff_M"] * 1e6
    log_c = np.log10(np.where(c_eff_uM > 0, c_eff_uM, 1e-12))
    im2 = ax.imshow(
        log_c,
        origin="lower",
        aspect="auto",
        cmap="viridis",
        extent=[lengths[0], lengths[-1], distances[0], distances[-1]],
    )
    cb2 = fig.colorbar(im2, ax=ax, fraction=0.046, pad=0.04)
    cb2.set_label("log₁₀ C_eff (µM)", fontsize=10)
    ax.set_xlabel("Linker length (residues)")
    ax.set_ylabel("Antigen–target distance (nm)")
    ax.set_title(
        "Effective local concentration of nanobody\n(when antibody is anchored)",
        fontsize=11,
    )

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def profile_plots(
    construct: ConditionalConstruct,
    output_path: str,
) -> None:
    """1-D profiles: occupancy vs distance for several linker lengths."""
    distances = np.linspace(0.5, 25, 300)
    linker_lengths = [20, 40, 60, 80, 120]
    models = ["fjc", "wlc"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    for ax, model in zip(axes, models):
        for color, n in zip(PALETTE, linker_lengths):
            lm = LinkerModel(n_residues=n, model=model)
            occ = np.array([lm.occupancy(construct.nanobody_kd_M, d) for d in distances])
            ax.plot(
                distances,
                occ,
                color=color,
                lw=2,
                label=f"{n} res  (Lc={lm.contour_length_nm:.0f} nm)",
            )
        ax.axhline(0.5, color="black", linestyle=":", lw=1, alpha=0.5)
        ax.set_xlabel("Antigen–target distance (nm)")
        ax.set_ylabel("Nanobody occupancy (anchored)" if ax == axes[0] else "")
        ax.set_title(f"Model: {model.upper()}")
        ax.set_ylim(0, 1.05)
        ax.legend(fontsize=9, title="Linker length")
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    fig.suptitle(
        f"Conditional nanobody occupancy vs antigen–target distance\n"
        f"Nanobody Kd = {construct.nanobody_kd_M*1e9:.0f} nM",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def ceff_vs_distance_plot(
    construct: ConditionalConstruct,
    output_path: str,
) -> None:
    """Effective local concentration vs distance for several linker lengths."""
    distances = np.linspace(0.1, 25, 300)
    linker_lengths = [20, 40, 60, 80, 120]

    fig, ax = plt.subplots(figsize=(8, 5))

    for color, n in zip(PALETTE, linker_lengths):
        lm = LinkerModel(n_residues=n)
        c_eff = np.array([lm.effective_concentration_M(d) * 1e6 for d in distances])
        ax.semilogy(
            distances,
            np.where(c_eff > 0, c_eff, 1e-12),
            color=color,
            lw=2,
            label=f"{n} res",
        )

    # Mark Kd line
    kd_uM = construct.nanobody_kd_M * 1e6
    ax.axhline(kd_uM, color="black", linestyle="--", lw=1.5, label=f"Kd = {kd_uM*1000:.1f} nM")
    ax.fill_between(distances, kd_uM, 1e3, alpha=0.07, color="green", label="C_eff > Kd (bound)")

    ax.set_xlabel("Antigen–target distance (nm)")
    ax.set_ylabel("Effective [nanobody] (µM)")
    ax.set_title("Effective local nanobody concentration vs distance (FJC model)")
    ax.legend(fontsize=9)
    ax.set_ylim(bottom=1e-8)

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    construct = ConditionalConstruct(
        antibody_kd_M=args.kd_antibody,
        nanobody_kd_M=args.kd_nanobody,
        linker=LinkerModel(n_residues=60),  # placeholder for parameter space scan
    )

    distances_nm = np.linspace(1.0, 22.0, 50)
    linker_lengths = np.arange(10, 151, 5, dtype=int)

    print(f"\nScanning {len(distances_nm)} distances × {len(linker_lengths)} linker lengths ...")
    data = construct.parameter_space_heatmap(
        distances_nm=distances_nm,
        linker_lengths=linker_lengths,
        model=args.model,
    )

    heatmap_plot(
        data,
        kd_nM=args.kd_nanobody * 1e9,
        output_path=os.path.join(args.output_dir, "occupancy_heatmap.png"),
    )
    profile_plots(
        construct,
        output_path=os.path.join(args.output_dir, "occupancy_profiles.png"),
    )
    ceff_vs_distance_plot(
        construct,
        output_path=os.path.join(args.output_dir, "ceff_vs_distance.png"),
    )

    print(f"\nAll figures saved to ./{args.output_dir}/\n")


if __name__ == "__main__":
    main()
