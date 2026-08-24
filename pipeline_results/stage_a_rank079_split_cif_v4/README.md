# Stage A v4 — rank079 corrected split CIF (RFd3 + graft)

- **Input:** `rank079_s0_native_split_s296_model_0_split.cif` (A=VH, B=VL, C=CH1, D=CL, T=epitope)
- **Design:** minibinder (chain M) cross-linking VL (B) and CH1 (C) hotspots

## Files for PyMOL

| File | Description |
|------|-------------|
| **`structures/grafted/*_grafted.pdb`** | **Use these.** Input Fab copied exactly (A,B,C,D,T) + chain M |
| `structures/rfd3_stage_a_cifs.tar.gz` | Raw RFd3 output (Fab domains reorient — do not use for Fab inspection) |
| `final/stage_a_minibinder_summary.csv` | Minibinder sequences |

Grafted PDBs: chains A,B,C,D,T are **identical coordinates** to the input Fab; only chain M is new.
