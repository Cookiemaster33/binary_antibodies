# Lambda Instances

## Active — Stage 0 holo-only (static rank)

| Field | Value |
|---|---|
| Instance ID | `bfffb352fc8c450684b1d9cb39d301e9` |
| IP | `157.151.241.144` |
| Region | `us-east-1` |
| Type | `gpu_1x_a10` (A100 unavailable) |
| Mode | Whole interface (23 res), 1000 MPNN → static rank → holo Boltz top 100 |
| SSH key | `~/.ssh/cursor_lambda_ephemeral` |

**Monitor:**
```bash
ssh -i ~/.ssh/cursor_lambda_ephemeral ubuntu@157.151.241.144 'tail -f /home/ubuntu/pipeline/stage_0_pipeline.log'
```

**Terminate when done:**
```bash
python3 scripts/launch_stage_0.py --terminate --instance-id bfffb352fc8c450684b1d9cb39d301e9
```

**Expected runtime:** ~30 min MPNN + ~2 h holo Boltz (100 designs, no apo)
