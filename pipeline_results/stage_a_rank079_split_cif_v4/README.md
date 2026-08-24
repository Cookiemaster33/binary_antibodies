# Stage A v4 — rank079 corrected split CIF (RFd3)

RFd3 backbone generation with **corrected H/L orientation** in the split Fab input:

- **Input:** `rank079_s0_native_split_s296_model_0_split.cif` (A=VH, B=VL, C=CH1, D=CL)
- Boltz H/L swap fixed: H was light chain, L was heavy chain in holo output
- **VL (B) and CH1 (C) face each other** (~46 Å) for minibinder cross-link design
- **Unlinked minibinder** 35–55 aa between VL and CH1 hotspots

## Files

| File | Description |
|------|-------------|
| `final/stage_a_sequences.json` | Per-design VH/VL/MB/CH1/CL segments |
| `final/stage_a_minibinder_sequences.fasta` | Minibinder sequences only |
| `final/stage_a_minibinder_summary.csv` | Spreadsheet summary |
| `structures/rfd3_stage_a_cifs.tar.gz` | All 200 RFd3 CIFs |
| `structures/sample_cifs/` | 5 samples for PyMOL |

Contig: `A1-113,B1-107/0,35-55,C1-107,D1-107`
