# Hidden Trivalent Minibinder Switch — Pipeline

## Mechanism

A native Fab is re-engineered so that a **trivalent minibinder hub** (CH1-arm + VL-arm + target-arm) sits hidden in the VH–VL groove. VH–VL pairing is the **only** trigger: when antigen drives closure, the hub is ejected and the target face becomes available.

Critical requirement: **apo VH–VL pairing must be weak**; target-bound closure must still be geometrically allowed.

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

## Stage 0 — Split-chain MPNN partial de-greasing

**Goal:** Reduce intrinsic VH–VL coupling in apo while preserving holo closure when epitope `T` is present. We do **not** expect Boltz apo to show fully dissociated Fv — success is **relative weakening** plus a large holo−apo contact delta.

**Design (no RFd3):**
1. Build native Fab context PDB (`fab_stage_0_vhvL_interface.pdb`)
2. Build **split Fv** PDB with VH/VL translated 30 Å apart (`fab_stage_0_split_mpnn.pdb`)
3. **ProteinMPNN** redesigns **rim** interface framework residues only; **closure core** (tightest buried pairs from `vh_vl_contacts.csv`) stays native
4. **Boltz** apo: A+B+C+D; holo: A+B+C+D+T

**Validation (apo vs holo differential, relative to native ~113 interface contacts):**

| Metric | Target | Purpose |
|--------|--------|---------|
| Apo contacts / native | ≤ 45% | Weakened apo coupling |
| Holo contacts / native | ≥ 45% | Holo still pairs |
| Holo − apo contact delta | ≥ 15 | Conditional switch gap |
| Fv framework RMSD vs native (holo) | ≤ 3.5 Å | Holo looks like real Fab |
| CDR ↔ epitope contacts (holo) | ≥ 6 | Co-binding with target |
| VH–VL interface cross-clashes (holo) | 0 | Anti-wedge |
| `wedge_suspect` | false | Reject plug-like false positives |

We select on **maximum holo−apo contact delta** among designs passing filters.

**Build & launch:**
```bash
python scripts/build_stage_0_design_target.py
python scripts/launch_stage_0.py --no-wait --no-terminate
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
| `binary_antibodies/stage_0_scoring.py` | Apo/holo interface metrics |
| `scripts/build_stage_a_design_target.py` | Stage A (accepts `--fab-pdb`) |
