# Stage A v2 — rank079 unlinked minibinder (RFd3)

RFd3 backbone generation with redesigned Stage A layout:

- **Full Fab** fixed (A=VH, B=VL, C=CH1, D=CL, T=epitope)
- **VL/CL translated 45 Å** from VH/CH1 (hotspot surfaces farther apart)
- **Unlinked minibinder** between VL and CH1 (`B/0,35-55,C` contig — no fusion to VH/VL)

Source Fab: Stage 0 **rank079** from `stage_0_full_fab_fused_t025`.

## Files

| File | Description |
|------|-------------|
| `final/stage_a_sequences.json` | Per-design VH/VL/MB/CH1/CL segments |
| `final/stage_a_minibinder_sequences.fasta` | Minibinder sequences only |
| `final/stage_a_minibinder_summary.csv` | Spreadsheet summary |
| `structures/rfd3_stage_a_cifs.tar.gz` | All 200 RFd3 CIFs |
| `structures/sample_cifs/` | 5 spread samples |

## RFd3 contig

`A1-113,B1-107/0,35-55,C1-101,D1-113`

Output CIFs concatenate segments on a single chain in contig order: VH | VL | MB | CH1 | CL.
