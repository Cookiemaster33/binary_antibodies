#!/usr/bin/env python3
"""
nanobody_cdr_analysis.py
------------------------
Identify the CDR face of the 2Rs15d anti-HER2 nanobody (PDB 5MY6).
These residues are the TARGET for minibinder design in the CH1-kicker construct.

Usage
-----
    python scripts/nanobody_cdr_analysis.py --structures-dir structures
"""

import argparse, os, sys, json, csv
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np

# IMGT CDR definitions for VHH nanobodies
# CDR1: 27-32, CDR2: 52-56, CDR3: 99-111
NANOBODY_CDR_IMGT = {
    "CDR1": list(range(27, 33)),
    "CDR2": list(range(52, 57)),
    "CDR3": list(range(99, 112)),
}

# Framework residues flanking CDRs (also relevant for binding)
NANOBODY_PARATOPE_ADJACENT = [26, 33, 50, 51, 57, 97, 98, 112]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--structures-dir", default="structures")
    p.add_argument("--output-dir", default="figures")
    return p.parse_args()


def identify_cdr_residues(nanobody_domain):
    """
    Identify CDR residues present in the nanobody domain structure.
    Uses IMGT numbering (approximate based on sequence position for VHH).
    """
    all_res = sorted(set(
        r.get_id()[1] for r in nanobody_domain.residues
        if r.get_id()[0] == " "
    ))

    cdr_residues = {}
    for cdr_name, imgt_range in NANOBODY_CDR_IMGT.items():
        found = [r for r in imgt_range if r in all_res]
        cdr_residues[cdr_name] = found

    return cdr_residues


def compute_cdr_face_center(nanobody_domain) -> np.ndarray:
    """Compute the geometric center (centroid) of all CDR Cα atoms."""
    all_cdr_res = sum(NANOBODY_CDR_IMGT.values(), [])
    ca_coords = []
    for res in nanobody_domain.residues:
        if res.get_id()[1] in all_cdr_res and "CA" in res:
            ca_coords.append(res["CA"].get_coord())
    if not ca_coords:
        return np.zeros(3)
    return np.mean(ca_coords, axis=0)


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    from binary_antibodies.structures import load_structure, extract_variable_domain

    pdb_path = Path(args.structures_dir) / "5MY6.pdb"
    print(f"\nLoading 2Rs15d nanobody from {pdb_path}")
    struct = load_structure(pdb_path, "5MY6")
    nb = extract_variable_domain(struct, chain_id="B", domain_type="VHH",
                                  residue_range=(1, 140))
    print(f"  {len(nb)} residues loaded")

    # Identify CDR residues
    cdr_residues = identify_cdr_residues(nb)
    all_cdr_res = sum(cdr_residues.values(), [])

    print(f"\n=== CDR residues (IMGT numbering) ===")
    for cdr, res in cdr_residues.items():
        print(f"  {cdr}: {res}")

    print(f"\n=== CDR face for minibinder design ===")
    print(f"  All CDR residues: {all_cdr_res}")
    print(f"  CDR + flanking:   {sorted(set(all_cdr_res + NANOBODY_PARATOPE_ADJACENT))}")

    # Save for use in RFdiffusion
    target_spec = {
        "pdb_file": "structures/domains/her2_nanobody_VHH.pdb",
        "chain": "A",  # chain in the domain PDB
        "cdr_residues": {k: v for k, v in cdr_residues.items()},
        "all_cdr_residues": all_cdr_res,
        "paratope_residues": sorted(set(all_cdr_res + NANOBODY_PARATOPE_ADJACENT)),
        "rfd3_hotspot_string": ",".join(
            f"A{r}" for r in sorted(set(all_cdr_res + NANOBODY_PARATOPE_ADJACENT))
            if r <= 130
        ),
        "design_note": (
            "Minibinder must bind the NANOBODY CDR face (anti-idiotypic design). "
            "This blocks the nanobody from binding HER2 in the OFF state. "
            "The CH1 kicker displaces this minibinder upon VH1-VL1 pairing."
        )
    }

    out_path = Path(args.structures_dir) / "interface" / "nanobody_cdr_target.json"
    out_path.parent.mkdir(exist_ok=True)
    out_path.write_text(json.dumps(target_spec, indent=2))
    print(f"\nSaved minibinder design target → {out_path}")

    print(f"\n=== RFdiffusion3 hotspot string ===")
    print(f"  select_hotspots: \"{target_spec['rfd3_hotspot_string']}\"")

    print(f"\n=== Kicker geometry analysis ===")
    from binary_antibodies.kicker import KickerConstruct

    # Example construct with reasonable parameters
    construct = KickerConstruct(
        kd_antigen_M=1e-9,
        kd_vh_vl_M=5e-6,
        kd_minibinder_cdrs_M=100e-9,  # 100 nM minibinder on nanobody CDRs
        kd_nanobody_M=10e-9,
        cl_to_nb_linker_n=30,          # short linker — keep nanobody close to CL
        name="2Rs15d anti-HER2 CH1-kicker construct",
    )
    print(construct.summary())

    # Scan linker lengths
    print(f"\n=== Linker length scan ===")
    print(f"{'Linker (res)':>12}  {'Overlap (nm)':>12}  {'Disp. frac.':>11}  {'OFF avail.':>10}")
    print("  " + "-"*52)
    for n in [10, 15, 20, 25, 30, 40, 50]:
        c = KickerConstruct(
            kd_antigen_M=1e-9, kd_vh_vl_M=5e-6,
            kd_minibinder_cdrs_M=100e-9, kd_nanobody_M=10e-9,
            cl_to_nb_linker_n=n,
        )
        kg = c.kicker
        print(f"  {n:>10}  {kg.kicking_overlap_nm():>12.2f}  "
              f"{kg.effective_displacement_fraction():>11.1%}  "
              f"{c.off_state_nanobody_availability():>10.2%}")


if __name__ == "__main__":
    main()
