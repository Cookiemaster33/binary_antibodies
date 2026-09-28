#!/usr/bin/env python3
"""
setup_example.py
----------------
Download example PDB structures, extract variable domains, analyse the
VH–VL interface, and produce the inputs needed for the minibinder design
pipeline.

Usage
-----
    python scripts/setup_example.py --output-dir structures

Output
------
  structures/
    1N8Z.pdb         Trastuzumab Fab  (anti-HER2) — VH1 + VL1 source
    1YY9.pdb         Cetuximab Fab    (anti-EGFR) — VH template for antigen binding
    5MY6.pdb         2Rs15d nanobody  (anti-HER2) — conditional nanobody
  structures/domains/
    trastuzumab_VH.pdb
    trastuzumab_VL.pdb
    cetuximab_VH.pdb
    her2_nanobody_VHH.pdb
    sequences.fasta
  structures/interface/
    vh_vl_contacts.csv         All VH–VL contacts with distances
    vl_minibinder_target.txt   VL residues to target for minibinder design
    interface_summary.json
"""

import argparse
import csv
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np

from binary_antibodies.structures import (
    download_pdb,
    load_structure,
    extract_variable_domain,
    analyse_vh_vl_interface,
    interface_summary,
    buried_sasa_estimate,
    write_fasta,
    EXAMPLE_STRUCTURES,
    print_design_overview,
    AntibodyDomain,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Download and process example structures for conditional nanobody design.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--output-dir", type=str, default="structures")
    p.add_argument("--skip-download", action="store_true",
                   help="Skip PDB download if files already exist.")
    return p.parse_args()


def save_domain_pdb(domain: AntibodyDomain, output_path: Path) -> None:
    """Save extracted domain residues as a minimal PDB file."""
    try:
        from Bio.PDB import PDBIO, Structure as S, Model as M, Chain as Ch
        from Bio import PDB

        struct = PDB.Structure.Structure(domain.pdb_id)
        model = PDB.Model.Model(0)
        chain = PDB.Chain.Chain(domain.chain_id)
        for res in domain.residues:
            chain.add(res.copy())
        model.add(chain)
        struct.add(model)

        io = PDBIO()
        io.set_structure(struct)
        io.save(str(output_path))
        print(f"    Saved domain PDB → {output_path}")
    except Exception as e:
        print(f"    Warning: could not save domain PDB: {e}")


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    domains_dir = output_dir / "domains"
    interface_dir = output_dir / "interface"
    domains_dir.mkdir(parents=True, exist_ok=True)
    interface_dir.mkdir(parents=True, exist_ok=True)

    print_design_overview()
    print("=" * 60)
    print("Step 1: Download PDB structures")
    print("=" * 60)

    # ── Download ──────────────────────────────────────────────────────
    pdb_paths = {}
    for name, info in EXAMPLE_STRUCTURES.items():
        pdb_paths[name] = download_pdb(
            info["pdb_id"],
            output_dir=output_dir,
            overwrite=not args.skip_download,
        )

    print()
    print("=" * 60)
    print("Step 2: Extract variable domains")
    print("=" * 60)

    domains = {}

    # Trastuzumab VH + VL
    trast_struct = load_structure(pdb_paths["trastuzumab_fab"], "1N8Z")
    trast_vh = extract_variable_domain(trast_struct, chain_id="A", domain_type="VH")
    trast_vh.source_antibody = "Trastuzumab (anti-HER2)"
    trast_vl = extract_variable_domain(trast_struct, chain_id="B", domain_type="VL")
    trast_vl.source_antibody = "Trastuzumab (anti-HER2)"
    domains["trastuzumab_VH"] = trast_vh
    domains["trastuzumab_VL"] = trast_vl
    print(f"  Trastuzumab VH: {len(trast_vh)} residues, chain A")
    print(f"  Trastuzumab VL: {len(trast_vl)} residues, chain B")

    # Cetuximab VH (template for antigen-binding engineering)
    cetux_struct = load_structure(pdb_paths["cetuximab_fab"], "1YY9")
    cetux_vh = extract_variable_domain(cetux_struct, chain_id="D", domain_type="VH")
    cetux_vh.source_antibody = "Cetuximab (anti-EGFR)"
    cetux_vl = extract_variable_domain(cetux_struct, chain_id="C", domain_type="VL")
    cetux_vl.source_antibody = "Cetuximab (anti-EGFR)"
    domains["cetuximab_VH"] = cetux_vh
    domains["cetuximab_VL"] = cetux_vl
    print(f"  Cetuximab VH:   {len(cetux_vh)} residues, chain D")
    print(f"  Cetuximab VL:   {len(cetux_vl)} residues, chain C")

    # Anti-HER2 nanobody
    nb_struct = load_structure(pdb_paths["her2_nanobody"], "5MY6")
    nb_vh = extract_variable_domain(nb_struct, chain_id="B", domain_type="VHH",
                                     residue_range=(1, 140))
    nb_vh.source_antibody = "2Rs15d anti-HER2 nanobody"
    domains["her2_nanobody_VHH"] = nb_vh
    print(f"  HER2 nanobody:  {len(nb_vh)} residues, chain B")

    # Save domains as PDB files
    for name, dom in domains.items():
        save_domain_pdb(dom, domains_dir / f"{name}.pdb")

    # Write FASTA
    write_fasta(list(domains.values()), domains_dir / "sequences.fasta")

    print()
    print("=" * 60)
    print("Step 3: Analyse Trastuzumab VH–VL interface")
    print("=" * 60)
    print("  (This identifies the VL1 residues that minibinders must target)")
    print()

    contacts = analyse_vh_vl_interface(trast_vh, trast_vl, contact_threshold_A=8.0)
    summary = interface_summary(contacts)
    bsasa = buried_sasa_estimate(contacts)

    print(f"  Total VH–VL contacts (< 8 Å): {summary['n_contacts']}")
    print(f"  VL contact residues:  {summary['vl_contact_residues']}")
    print(f"  VH contact residues:  {summary['vh_contact_residues']}")
    print(f"  VL framework contacts (minibinder target): {summary['vl_framework_contacts']}")
    print(f"  Estimated buried SASA (VL side): {bsasa:.0f} Å²")
    print()
    print(f"  Top 5 closest VH–VL contacts:")
    for c in summary["closest_contacts"][:5]:
        from Bio.SeqUtils import seq1
        try:
            vh_aa = seq1(c.vh_resname)
            vl_aa = seq1(c.vl_resname)
        except Exception:
            vh_aa = c.vh_resname[:3]
            vl_aa = c.vl_resname[:3]
        print(f"    VH {vh_aa}{c.vh_resnum} ↔ VL {vl_aa}{c.vl_resnum}  "
              f"min dist = {c.min_heavy_distance_A:.2f} Å  "
              f"Cα dist = {c.ca_distance_A:.2f} Å")

    # Save contacts CSV
    csv_path = interface_dir / "vh_vl_contacts.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "vh_resnum", "vh_resname", "vl_resnum", "vl_resname",
            "ca_distance_A", "min_heavy_distance_A",
        ])
        writer.writeheader()
        for c in contacts:
            writer.writerow({
                "vh_resnum": c.vh_resnum,
                "vh_resname": c.vh_resname,
                "vl_resnum": c.vl_resnum,
                "vl_resname": c.vl_resname,
                "ca_distance_A": round(c.ca_distance_A, 3),
                "min_heavy_distance_A": round(c.min_heavy_distance_A, 3),
            })
    print(f"\n  Saved contacts → {csv_path}")

    # Save summary JSON
    summary_out = {
        "trastuzumab_vh_vl_interface": {
            "n_contacts": summary["n_contacts"],
            "vl_contact_residues": summary["vl_contact_residues"],
            "vh_contact_residues": summary["vh_contact_residues"],
            "vl_framework_contacts": summary["vl_framework_contacts"],
            "estimated_buried_sasa_A2": round(bsasa, 1),
        }
    }
    json_path = interface_dir / "interface_summary.json"
    json_path.write_text(json.dumps(summary_out, indent=2))
    print(f"  Saved summary → {json_path}")

    # Save minibinder target residue list
    mb_target_path = interface_dir / "vl_minibinder_target.txt"
    mb_target_path.write_text(
        "# VL1 residues to target for minibinder design\n"
        "# Source: Trastuzumab VL chain B, residues contacting VH (< 8 Å)\n"
        "# These residues should be within the minibinder binding footprint.\n\n"
        "VL_framework_contacts:\n"
        + "\n".join(f"  {r}" for r in summary["vl_framework_contacts"])
        + "\n\nAll_VL_contact_residues:\n"
        + "\n".join(f"  {r}" for r in summary["vl_contact_residues"])
    )
    print(f"  Saved minibinder target → {mb_target_path}")

    print()
    print("=" * 60)
    print("Step 4: Design summary")
    print("=" * 60)
    print()
    print("  Next steps for wet-lab / computational pipeline:")
    print()
    print("  [1] VH1 engineering (anti-EGFR + VL1-pairing)")
    print("      • Start from Trastuzumab VH (structures/domains/trastuzumab_VH.pdb)")
    print("      • Graft CDRs of anti-EGFR antibody (Cetuximab: structures/domains/cetuximab_VH.pdb)")
    print("        onto the Trastuzumab VH framework to maintain VL1-pairing geometry")
    print("      • OR: use RFdiffusion to de novo design a VH that binds EGFR while")
    print("        preserving the canonical VH–VL interface geometry")
    print("      • Weaken VH–VL affinity via FR2 mutations (target Kd_VH-VL ≈ 1–10 µM)")
    print()
    print("  [2] Minibinder design (locks VL1 in free state)")
    print("      • Target: VL1 framework 2 residues")
    print(f"        Residues: {summary['vl_framework_contacts']}")
    print("      • Input PDB: structures/domains/trastuzumab_VL.pdb")
    print("      • Tool: RFdiffusion (partial diffusion mode) + ProteinMPNN")
    print("      • Target Kd: 10–200 µM (weak intrinsic; intramolecular C_eff locks VL1)")
    print("      • Validation: SPR or BLI against VL1 alone and VH1–VL1 complex")
    print()
    print("  [3] Nanobody (conditional arm)")
    print("      • Ready: structures/domains/her2_nanobody_VHH.pdb")
    print("      • 2Rs15d (PDB 5MY6) binds HER2 domain II")
    print("      • Different epitope from Trastuzumab → no steric clash expected")
    print()
    print("  [4] Chain B assembly sequence (N→C)")
    print("      SP – [Minibinder] – (G4S)30 – [Minibinder] – (G4S)30 – [VL1] – (G4S)12 – [VHH]")
    print()
    print("All structure files saved to:", str(output_dir))


if __name__ == "__main__":
    main()
