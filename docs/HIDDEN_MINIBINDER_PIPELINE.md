# Hidden Trivalent Minibinder Switch — Pipeline

## Mechanism

A native Fab is re-engineered so that a **trivalent minibinder hub** (CH1-arm + VL-arm + target-arm) sits hidden in the VH–VL groove. VH–VL pairing is the **only** trigger: when antigen drives closure, the hub is ejected and the target face becomes available.

Critical requirement: **designed VH–VL interface must be much weaker than WT**, while **target-bound holo** must still engage the epitope normally.

## Pipeline order

```
Stage 0  →  Stage A  →  Stage B
VH–VL       minibinder    target arm
weaken      hub design    (de novo)
```

| Stage | What is designed | What stays fixed |
|-------|------------------|------------------|
| **0** | VH + VL **interface rim** (partial de-grease) | CDRs, closure-core interface, CH1, CL, epitope stub T |
| **A** | VH–hub–VL **minibinder** (35–55 aa) | Stage 0 Fv, CH1, CL, epitope stub |
| **B** | Third arm against epitope / target | Stage A assembly |

## Stage 0 — Split-chain MPNN de-greasing

**Goal:** Reduce intrinsic VH–VL coupling relative to WT while preserving epitope engagement when stub `T` is present. Apo Boltz folds are **not** required.

**Design (no RFd3):**
1. Build native Fab context PDB (`fab_stage_0_vhvL_interface.pdb`)
2. **PISA on WT Fab** defines buried interface residues (CDRs excluded on Fv):
   - `fv` (default): VH–VL only (chains A+B)
   - `full_fab`: VH–VL (A+B) **and** CH1–CL (C+D)
3. Build **split MPNN PDB** — Fv-only (`fv`) or full Fab with both interfaces separated (`full_fab`)
4. **ProteinMPNN** redesigns PISA-defined interface framework residues
   - **Whole interface (default):** all PISA framework interface residues
   - **Aggressive:** PISA core (BSA ≥ 5 Å²) fixed → rim redesigned
   - **Conservative:** PISA core (BSA ≥ 10 Å²) fixed → rim redesigned
5. Rank MPNN output by **predicted static VH–VL contacts** on the native Fab backbone (weakest first)
6. **Boltz holo only:** A+B+C+D+T

**Validation:**

| Criterion | Metric | Target |
|-----------|--------|--------|
| (a) Weaker interface than WT | `static_vh_vl_interface_contacts` / native | ≤ 45% of WT |
| (b) Target binding preserved | CDR ↔ epitope contacts (holo) | ≥ 6 |
| (b) Fab-like holo geometry | Fv **framework** Cα RMSD vs native (holo) | ≤ 3.5 Å |
| (b) CDR geometry vs native | Fv **CDR** Cα RMSD vs native (holo, same alignment) | ≤ 6.0 Å |
| (b) No interface steric block | VH–VL interface cross-clashes (holo) | 0 |
| (c) Weaker VH–VL energetics (ranking) | PISA `pisa_delta_int_solv_en_vs_native_kcal` on holo Fv | higher = weaker |

We select on **lowest static interface contacts** among designs passing holo filters, with **PISA solvation-energy delta** as the primary tie-breaker when static scores are similar.

**PISA:** Post-Boltz, `pdbegroup/pisa` Docker scores VH–VL (chains A+B) interface area and solvation energy vs native WT reference. Positive `pisa_delta_int_solv_en_vs_native_kcal` means a less favorable (weaker) interface than native.

**Build & launch:**
```bash
python scripts/build_stage_0_design_target.py                        # fv scope (default)
python scripts/build_stage_0_design_target.py --interface-scope full_fab
python scripts/launch_stage_0.py --interface-scope fv --no-wait --no-terminate
python scripts/launch_stage_0.py --interface-scope full_fab --no-wait --no-terminate
```

**Outputs:** `pipeline_results/stage_0_vhvL_interface/`

## Stage A — Hidden minibinder hub

**Prerequisite:** Top Stage 0 Fab (chains A–D).

```bash
python scripts/build_stage_a_design_target.py \
  --fab-pdb path/to/stage0_top_holo.cif
python scripts/launch_stage_a.py --no-wait --no-terminate
```

**RFd3 contig:** `A1-113,35-55,B1-107`  
**Hotspots:** CH1 (C) + VL framework (B)

## Stage B — Target arm (future)

Add hotspots on chain `T` and extend contig with a third binding patch against the HER2 epitope stub.

## Key files

| File | Purpose |
|------|---------|
| `binary_antibodies/fab_hidden_switch.py` | Shared constants, split PDB builder, degrease residue lists |
| `scripts/build_stage_0_design_target.py` | Stage 0 PDBs + JSON config |
| `scripts/gpu_setup/run_stage_0_vhvL_interface.sh` | Split MPNN → Boltz → scoring |
| `scripts/launch_stage_0.py` | Lambda launcher |
| `binary_antibodies/stage_0_scoring.py` | Static + holo + PISA interface metrics |
| `binary_antibodies/pisa_scoring.py` | PISA Docker wrapper and XML parser |
| `scripts/build_stage_a_design_target.py` | Stage A (accepts `--fab-pdb`) |
