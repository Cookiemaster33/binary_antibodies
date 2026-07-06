"""
RFdiffusion3 round-2 partial diffusion refinement.

Takes top round-1 backbones, adds controlled noise (partial_t), and re-diffuses
only the minibinder region while pinning VH1 and nanobody anchors.

Designed to run inside the foundry Docker container.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import torch
from atomworks.io.utils.io_utils import load_any, to_cif_file
from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine
from rfd3.inference.input_parsing import DesignInputSpecification

VH1_END = 115
NB_LEN = 115

# VH1 CH1-face (excluding VH-VL overlap 39, 41, 85)
VH1_HOTSPOT_RES = [
    11, 12, 13, 14, 15, 40,
    79, 80, 81, 82, 83, 84,
    103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115,
]

# Nanobody CDR residues on original chain C
NB_CDR_RES = (
    list(range(27, 34))
    + list(range(52, 58))
    + list(range(99, 113))
)


def detect_domains(total_len: int) -> tuple[int, int, int]:
    """Return (mb_start, mb_end, nb_start) for a connected VH1-MB-Nb chain."""
    mb_len = total_len - VH1_END - NB_LEN
    if mb_len < 1:
        raise ValueError(f"Invalid chain length {total_len}: minibinder length {mb_len}")
    mb_start = VH1_END + 1
    mb_end = VH1_END + mb_len
    nb_start = mb_end + 1
    return mb_start, mb_end, nb_start


def remap_hotspots(mb_len: int) -> str:
    """Map round-1 A/C hotspots onto the connected chain A numbering."""
    _, _, nb_start = detect_domains(VH1_END + mb_len + NB_LEN)
    parts = [f"A{r}" for r in VH1_HOTSPOT_RES]
    parts += [f"A{nb_start + r - 1}" for r in NB_CDR_RES]
    return ",".join(parts)


def select_templates(round1_dir: Path) -> list[Path]:
    """Pick round-1 RFd3 CIF templates for partial diffusion."""
    explicit = os.environ.get("RFD3_ROUND2_TEMPLATE_BACKBONES", "").strip()
    if explicit:
        templates: list[Path] = []
        for stem in explicit.split(","):
            stem = stem.strip()
            if not stem:
                continue
            path = round1_dir / (stem if stem.endswith(".cif") else f"{stem}.cif")
            if path.exists():
                templates.append(path)
            else:
                print(f"  WARNING: template not found, skipping: {path}")
        return templates

    # Prefer lowest global assembly RMSD backbones from round-1 Boltz validation.
    results_path = Path(
        os.environ.get(
            "RFD3_ROUND1_RESULTS",
            "/workspace/final/round1_final_results.json",
        )
    )
    n = int(os.environ.get("RFD3_ROUND2_TEMPLATES", 5))
    if results_path.exists():
        results = json.load(open(results_path))
        results.sort(
            key=lambda r: (
                r["global_rmsd_A"] if r.get("global_rmsd_A", 999) < 900 else 999,
                -r.get("boltz_plddt_pct", 0),
            )
        )
        templates = []
        seen: set[str] = set()
        for row in results:
            bb = row.get("backbone", "")
            if not bb or bb in seen:
                continue
            if "global_rmsd_A" not in row:
                print(f"  WARNING: {bb} missing global_rmsd_A — re-score round 1 first")
                continue
            rmsd = row["global_rmsd_A"]
            if rmsd >= 900:
                continue
            path = round1_dir / f"{bb}.cif"
            if not path.exists():
                print(f"  WARNING: round-1 RFd3 CIF missing for {bb}")
                continue
            templates.append(path)
            seen.add(bb)
            print(f"  template {len(templates)}: {bb} (global RMSD={rmsd:.2f} Å)")
            if len(templates) >= n:
                break
        if templates:
            print(f"Selected {len(templates)} templates from round-1 global RMSD ranking")
            return templates
        print("WARNING: round-1 results found but no valid templates — falling back")

    # Fallback: evenly sample round-1 CIFs (no scoring available).
    cifs = sorted(round1_dir.glob("mb_*.cif"))
    if not cifs:
        return []
    if len(cifs) <= n:
        return cifs
    step = len(cifs) / n
    return [cifs[int(i * step)] for i in range(n)]


def build_round2_input(
    template_cif: Path,
    original_pdb: Path,
    out_cif: Path,
    *,
    vl_context: bool,
) -> tuple[int, int, int]:
    """
    Build RFd3 input for partial diffusion.

    Returns (mb_start, mb_end, nb_start).
 
    When vl_context=True, chain B (VL) from the original input PDB is appended
    as fixed steric context (same as round 1).
    """
    raw = load_any(str(template_cif))
    aa = raw[0] if hasattr(raw, "__getitem__") else raw

    ch = list(set(aa.chain_id))[0]
    total = len(set(aa.res_id[aa.chain_id == ch]))
    mb_start, mb_end, nb_start = detect_domains(total)

    if vl_context and original_pdb.exists():
        raw_orig = load_any(str(original_pdb))
        aa_orig = raw_orig[0] if hasattr(raw_orig, "__getitem__") else raw_orig
        vl_mask = aa_orig.chain_id == "B"
        if np.any(vl_mask):
            aa_vl = aa_orig[vl_mask]
            # atomworks/biotite-style concatenation
            combined = aa + aa_vl
            to_cif_file(combined, str(out_cif))
            return mb_start, mb_end, nb_start

    to_cif_file(aa, str(out_cif))
    return mb_start, mb_end, nb_start


def run_round2() -> int:
    torch.set_float32_matmul_precision("high")

    round1_dir = Path(os.environ.get("RFD3_ROUND1_DIR", "/workspace/outputs/rfd3"))
    round2_dir = Path(os.environ.get("RFD3_ROUND2_DIR", "/workspace/outputs/rfd3_round2"))
    scratch_dir = Path(os.environ.get("RFD3_ROUND2_SCRATCH", "/workspace/outputs/rfd3_round2_inputs"))
    original_pdb = Path(
        os.environ.get(
            "RFD3_INPUT_PDB",
            "/workspace/inputs/vh1_vl_nanobody_design_target.pdb",
        )
    )

    partial_t = float(os.environ.get("RFD3_PARTIAL_T", "2.0"))
    designs_per = int(os.environ.get("RFD3_ROUND2_DESIGNS_PER_TEMPLATE", 8))
    batch_size = int(os.environ.get("RFD3_ROUND2_BATCH_SIZE", designs_per))
    vl_context = os.environ.get("RFD3_ROUND2_VL_CONTEXT", "1") not in ("0", "false", "False")

    round2_dir.mkdir(parents=True, exist_ok=True)
    scratch_dir.mkdir(parents=True, exist_ok=True)

    templates = select_templates(round1_dir)
    if not templates:
        print("Round 2: no templates found — skipping.")
        return 0

    print(f"Round 2 partial diffusion: {len(templates)} templates × {designs_per} designs")
    print(f"  partial_t={partial_t} Å | VL context={vl_context}")

    model = RFD3InferenceEngine(
        **RFD3InferenceConfig(diffusion_batch_size=batch_size)
    )
    saved = 0

    for template in templates:
        stem = template.stem
        print(f"\n  Template: {stem}", flush=True)

        input_cif = scratch_dir / f"{stem}_input.cif"
        try:
            mb_start, mb_end, nb_start = build_round2_input(
                template, original_pdb, input_cif, vl_context=vl_context
            )
        except Exception as exc:
            print(f"    ERROR building input: {exc}")
            continue

        nb_end = nb_start + NB_LEN - 1
        fixed_atoms = {
            f"A1-{VH1_END}": "ALL",
            f"A{nb_start}-{nb_end}": "ALL",
        }
        if vl_context:
            fixed_atoms["B1-115"] = "ALL"

        hotspots = remap_hotspots(mb_end - VH1_END)
        print(f"    MB region A{mb_start}-{mb_end} | fixed VH1 + Nb"
              f"{' + VL' if vl_context else ''}")

        spec = DesignInputSpecification.safe_init(
            input=str(input_cif),
            partial_t=partial_t,
            select_hotspots=hotspots,
            select_fixed_atoms=fixed_atoms,
        )

        n_batches = max(1, (designs_per + batch_size - 1) // batch_size)
        produced = 0
        for bi in range(n_batches):
            for _key, designs in model.run(inputs=spec, out_dir=None, n_batches=1).items():
                for i, design in enumerate(designs):
                    if produced >= designs_per:
                        break
                    out_name = f"r2_{stem}_{produced:03d}.cif"
                    to_cif_file(design.atom_array, str(round2_dir / out_name))
                    produced += 1
                    saved += 1
                if produced >= designs_per:
                    break

        print(f"    Saved {produced} refined designs")

    print(f"\nRound 2 done: {saved} designs in {round2_dir}")
    return saved


if __name__ == "__main__":
    run_round2()
