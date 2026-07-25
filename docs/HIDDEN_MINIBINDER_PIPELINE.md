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
| **0** | VH + VL **framework** at the Fv interface | CDRs, CH1, CL, epitope stub T |
| **A** | VH–hub–VL **minibinder** (35–55 aa) | Stage 0 Fv, CH1, CL, epitope stub |
| **B** | Third arm against epitope / target | Stage A assembly |

## Stage 0 — VH–VL interface weakening

**Goal:** Reduce spontaneous VH–VL pairing in apo while preserving a clash-free closed pose when epitope `T` is present.

**RFd3 setup:**
- Input: `structures/domains/fab_stage_0_vhvL_interface.pdb`
- Contig: `A1-113/0,B1-107` (both Fv chains, no inserted domain)
- `partial_t = 12 Å` — larger than round-2 refinement; explores frustrated interfaces
- Hotspots: epitope stub `T1–T12` (closure geometry)
- Fixed: all CDRs, CH1, CL, epitope coordinates

**Validation (apo vs holo differential):**

| Metric | Apo | Holo | Purpose |
|--------|-----|------|---------|
| VH–VL interface contacts (heavy-atom) | ≤ 12 | ≥ 50 | Conditional pairing gap |
| Interface centroid distance | ≥ 14 Å | ≤ 12 Å | Open vs closed |
| Fv framework RMSD vs native | — | ≤ 3.5 Å | Holo looks like real Fab |
| CDR ↔ epitope contacts | — | ≥ 6 | Co-binding, not just VH–VL slam |
| VH–VL interface cross-clashes | — | 0 | Anti-wedge (no steric plugs) |
| `wedge_suspect` flag | — | false | Reject plug-like false positives |

We select on **maximum holo−apo contact delta** among designs passing all filters — not minimum apo affinity alone.

**Build & launch:**
```bash
python scripts/build_stage_0_design_target.py
python scripts/launch_stage_0.py --no-wait --no-terminate
```

**Outputs:** `pipeline_results/stage_0_vhvL_interface/` (after push script, TBD)

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

## Current PoC run (minibinder-first)

Instance `0926bc09188c4cebb0cd8872131e29cb` is running **Stage A on native Fab** as a structural feasibility check. Results are exploratory only.

**After PoC completes:** switch to Stage 0 → Stage A order using this document.

## Key files

| File | Purpose |
|------|---------|
| `binary_antibodies/fab_hidden_switch.py` | Shared constants, PDB builder |
| `scripts/build_stage_0_design_target.py` | Stage 0 PDB + JSON |
| `scripts/gpu_setup/run_stage_0_vhvL_interface.sh` | Full Stage 0 GPU pipeline |
| `scripts/launch_stage_0.py` | Lambda launcher |
| `binary_antibodies/stage_0_scoring.py` | Apo/holo interface metrics |
| `scripts/build_stage_a_design_target.py` | Stage A (accepts `--fab-pdb`) |
