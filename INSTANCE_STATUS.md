# Lambda Instances

No active instances.

## Last run — Stage A PoC (terminated)

| Field | Value |
|---|---|
| Instance ID | `0926bc09188c4cebb0cd8872131e29cb` (terminated) |
| Results | `pipeline_results/stage_a_poc_rfd3_only/` |
| Designs | 200 RFd3 CIFs (tarball) |
| Note | Exploratory minibinder-first PoC on native Fab |

## Next run (when ready)

Stage 0 → Stage A — see `docs/HIDDEN_MINIBINDER_PIPELINE.md`

```bash
python scripts/build_stage_0_design_target.py
python scripts/launch_stage_0.py --no-wait --no-terminate
```
