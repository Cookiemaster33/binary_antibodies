# Stage 0 — holo-only + static interface scoring

**Instance:** `bfffb352fc8c450684b1d9cb39d301e9` (A10, us-east-1) — **terminated**  
**Completed:** 2026-07-26  
**Pipeline:** 1000 MPNN → static rank → holo Boltz (100, `--no_kernels`) → scoring

## Result

**0/100** pass all filters (re-scored with split framework/CDR RMSD).

## Files

| File | Contents |
|------|----------|
| `stage_0_results.json` | Top 100 scored designs |
| `top_designs_stage_0.json` | Boltz input list + static ranks |
| `all_sequences.json` | Full 1000 MPNN sequences |
| `stage_0_pipeline.log` | Full pipeline log |
| `boltz_holo_run.log` | Boltz holo run (retry with --no_kernels) |
| `run_config.json` | Run metadata |

## Best candidate: rank009 (`s682`)

| Metric | Value |
|--------|-------|
| Static contacts | 53.8 (47.6% WT) |
| Holo VH–VL iface | 85 |
| Framework RMSD | ~2.5 Å |
| Interface rim RMSD | ~2.87 Å |
| CDR RMSD | ~0.9 Å |
| CDR–epitope | 209 |

Fails only: `interface_weakened` (static > 45% WT).
