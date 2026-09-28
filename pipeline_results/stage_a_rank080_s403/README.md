# Stage A — rank080 hidden minibinder (RFd3)

RFd3 backbone generation for VH–minibinder–VL hub seated in Stage 0 Fab **rank080**.

## Sequences (pulled from RFd3 CIFs)

| File | Description |
|------|-------------|
| `final/stage_a_sequences.json` | Full per-design data: VH, minibinder, VL, context chains |
| `final/stage_a_minibinder_sequences.fasta` | Minibinder segment only (200 designs) |
| `final/stage_a_minibinder_summary.csv` | Spreadsheet-friendly summary |

Minibinder lengths sampled by RFd3: **36–55 aa** (contig `35-55`).

## Structures

- `structures/sample_cifs/` — 5 spread samples for PyMOL/ChimeraX
- All **200 CIFs** remain on Lambda instance (see `run_config.json`)

## Not included

MPNN sequence design and Boltz global RMSD validation were not run.
