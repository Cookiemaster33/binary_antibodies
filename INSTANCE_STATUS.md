# Active Lambda Instance — v5 Two-Round Partial Diffusion

| Field | Value |
|---|---|
| Instance ID | `5314d4de9afc487f8611143bbfc70ffe` |
| IP | `132.145.134.114` |
| Status | **Booting → pipeline starting** |
| Results folder | `pipeline_results/v5_two_round_refine/` |

## Configuration

- `RFD3_ROUNDS=2` — validate round 1, then partial diffusion on top scRMSD RFd3 backbones
- Round 1: 200 designs, MB length 35–70
- Round 2: top 5 templates × 8 designs, `partial_t=2.0`

## Monitor

```bash
tail -f /workspace/pipeline_launch.log

# Or on the instance once SSH is up:
ssh -i ~/.ssh/lambda_agent_key ubuntu@132.145.134.114 \
  "tmux -f /exec-daemon/tmux.portal.conf attach -t pipeline"
```

## Terminate when done

```bash
python3 scripts/launch_full_pipeline.py --terminate --instance-id 5314d4de9afc487f8611143bbfc70ffe
```
