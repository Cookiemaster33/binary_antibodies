# Pipeline Run — v5 Two-Round Partial Diffusion

## Status: **Ready to launch** (waiting for `LAMBDA_API_KEY`)

Results will be saved to a **new folder** (does not overwrite previous runs):
`pipeline_results/v5_two_round_refine/`

## Launch command

Add `LAMBDA_API_KEY` to your Cursor Cloud Agent secrets, then re-run the agent or:

```bash
export LAMBDA_API_KEY=<your-key>
python scripts/launch_full_pipeline.py \
    --ssh-key ~/.ssh/lambda_agent_key \
    --rfd3-rounds 2 \
    --results-dir pipeline_results/v5_two_round_refine
```

## Run configuration

| Parameter | Value |
|-----------|-------|
| `RFD3_ROUNDS` | 2 |
| Round 1 designs | 200 |
| MB length range | 35–70 |
| Round 2 templates | 5 (lowest round-1 scRMSD RFd3 backbones) |
| Designs per template | 8 |
| `partial_t` | 2.0 Å |

## Workflow

```
Round 1 RFd3 → MPNN → Boltz → scRMSD
         ↓ pick top 5 by scRMSD (RFd3 CIFs)
Round 2 partial diffusion → MPNN → Boltz → final scRMSD
         ↓
GitHub: pipeline_results/v5_two_round_refine/
```

## Monitor / terminate

```bash
python scripts/launch_full_pipeline.py --status
python scripts/launch_full_pipeline.py --terminate --instance-id <id>
```
