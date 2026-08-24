# Stage A PoC — RFd3 only (native Fab)

Exploratory proof-of-concept: hidden minibinder hub (`A1-113,35-55,B1-107`) designed
into **native** trastuzumab VH–VL without prior Stage 0 interface weakening.

## Contents

| File | Description |
|------|-------------|
| `structures/rfd3_stage_a_cifs.tar.gz` | All 200 RFd3 output CIFs |
| `manifest.json` | Design filenames |
| `stage_a_pipeline.log` | Full pipeline log |
| `inputs/` | Design target PDB + RFd3 config |
| `run_config.json` | Run metadata |

## Not included

MPNN, Boltz, and global RMSD validation were not run (RFd3-only PoC).

## Next steps

Production pipeline: **Stage 0** (VH–VL weakening) → **Stage A** with `--fab-pdb`.
