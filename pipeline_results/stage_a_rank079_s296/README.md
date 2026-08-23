# Stage A — rank079 hidden minibinder (RFd3)

RFd3 backbone generation for VH–minibinder–VL hub seated in Stage 0 Fab **rank079**
(fused H/L Boltz holo, top PISA Δ from `stage_0_full_fab_fused_t025`).

## Sequences

| File | Description |
|------|-------------|
| `final/stage_a_sequences.json` | Full per-design data: VH, minibinder, VL |
| `final/stage_a_minibinder_sequences.fasta` | Minibinder segment only (200 designs) |
| `final/stage_a_minibinder_summary.csv` | Spreadsheet-friendly summary |

Minibinder lengths sampled by RFd3: **35–55 aa** (contig `35-55`).

## Structures

- `structures/sample_cifs/` — 5 spread samples for PyMOL/ChimeraX
- `structures/rfd3_stage_a_cifs.tar.gz` — all **200** RFd3 output CIFs

## Not included

MPNN sequence design and Boltz global RMSD validation were not run (RFd3-only).
