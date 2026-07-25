# Stage 0 — VH–VL Interface Weakening

Partial-diffusion redesign of trastuzumab VH/VL framework (CDRs + CH1/CL/epitope fixed).

## Summary

| Metric | Value |
|--------|-------|
| RFd3 designs | 200 |
| MPNN sequences | 1600 |
| Boltz validated (top 50) | 50 apo + 50 holo |
| Passed all filters | **0/50** |

**Best candidate:** `rank32_s0_b014_003_s0`
- Holo − apo contact delta: **41**
- Holo Fv framework RMSD vs native: **1.86 Å**
- Apo interface contacts: **22** (filter ≤12 — still too paired in apo Boltz fold)
- CDR–epitope contacts (holo): 245

## Files

| File | Description |
|------|-------------|
| `stage_0_results.json` | Full scoring results (ranked top 50) |
| `top_designs_stage_0.json` | Boltz input metadata |
| `mpnn_all_sequences.json` | MPNN sequences for all backbones |
| `stage_0_pipeline.log` | Full pipeline log |
| `inputs/` | Design target PDB + RFd3 config |
| `structures/stage_0_rfd3_cifs.tar.gz` | All 200 RFd3 CIFs |
| `structures/stage_0_boltz_apo.tar.gz` | Boltz apo predictions (A+B+C+D) |
| `structures/stage_0_boltz_holo.tar.gz` | Boltz holo predictions (A+B+C+D+T) |
| `structures/stage_0_top_structures.tar.gz` | Top 10 apo+holo Boltz CIFs |

## Next steps

Review `stage_0_results.json` and decide:
- Relax apo filters (Boltz may fold Fv as paired even without epitope)
- Rerun with more designs / higher `partial_t`
- Proceed to Stage A with best partial hit for exploratory minibinder design
