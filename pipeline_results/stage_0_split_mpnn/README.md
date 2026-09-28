# Stage 0 — Split-chain MPNN Partial De-grease

Split Fv ProteinMPNN on 14 rim interface residues (closure core fixed); no RFd3.

## Summary

| Metric | Value |
|--------|-------|
| MPNN sequences | 32 |
| Boltz validated | 32 apo + 32 holo |
| Passed all filters | **0/32** |
| Runtime | ~21 min (A100, us-west-2) |
| Instance | `feda9ebe3c464be69007fe9eec23ac2e` (terminated) |

**Best candidate:** `rank23_s0_native_split_s0`
- Holo − apo contact delta: **51**
- Holo Fv framework RMSD vs native: **1.02 Å**
- Apo interface contacts: **68** (60% of native ~113 — filter ≤45%)
- Holo interface contacts: **119** (105% of native)
- CDR–epitope contacts (holo): 216
- Fails: apo not weakened enough; 1 VH–VL interface clash

**Second best:** `rank22_s0_native_split_s0` — Δ=40, holo=115, RMSD=4.4 Å

## vs previous RFd3 run (`stage_0_vhvL_interface/`)

| | Split MPNN | RFd3 + MPNN |
|---|---|---|
| Best apo contacts | 62 (55% native) | **22** (19% native) |
| Best holo−apo Δ | **51** | 41 |
| Best holo RMSD | **1.02 Å** | 1.86 Å |
| Designs with Δ ≥ 15 | 2/32 | several |

Split MPNN improved holo fidelity and conditional gap for the top hit, but **did not weaken apo** as much as RFd3. 23/32 designs showed inverted coupling (holo contacts < apo).

## Files

| File | Description |
|------|-------------|
| `stage_0_results.json` | Full scoring results (all 32 designs, ranked) |
| `top_designs_stage_0.json` | Boltz input metadata |
| `mpnn_all_sequences.json` | All MPNN sequences |
| `stage_0_pipeline.log` | Full pipeline log |
| `run_config.json` | Run parameters |
| `inputs/` | Fab PDB, split MPNN PDB, config |
| `structures/stage_0_boltz_apo.tar.gz` | Boltz apo predictions (A+B+C+D) |
| `structures/stage_0_boltz_holo.tar.gz` | Boltz holo predictions (A+B+C+D+T) |
| `structures/stage_0_top_structures.tar.gz` | Top 10 apo+holo Boltz CIFs |

## Next steps

- Review `rank23` and `rank22` for Stage A (strong holo closure, large Δ)
- Consider hybrid: RFd3 backbone + split MPNN de-greasing
- Expand degrease residue set or relax apo filter for selection
