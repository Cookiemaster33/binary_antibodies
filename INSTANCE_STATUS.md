# Lambda Instances

## Active — Stage 0 aggressive split MPNN (running)

| Field | Value |
|---|---|
| Instance ID | `6195d0acc0ae4c1d9b038db8fc60ce15` |
| IP | `150.136.35.109` |
| Region | `us-east-1` |
| Type | `gpu_1x_a100_sxm4` |
| Mode | Aggressive de-grease (19 rim residues, core ≤ 3.20 Å) |
| MPNN seqs | 64 |
| SSH key | `~/.ssh/cursor_lambda_ephemeral` |

**Monitor:**
```bash
ssh -i ~/.ssh/cursor_lambda_ephemeral ubuntu@150.136.35.109 'tail -f /home/ubuntu/pipeline/stage_0_pipeline.log'
```

**Terminate when done:**
```bash
python3 scripts/launch_stage_0.py --terminate --instance-id 6195d0acc0ae4c1d9b038db8fc60ce15
```

## Previous — Stage 0 partial split MPNN (terminated)

| Field | Value |
|---|---|
| Instance ID | `feda9ebe3c464be69007fe9eec23ac2e` |
| Results | `pipeline_results/stage_0_split_mpnn/` |
| Outcome | 0/32 passed; best Δ=51, apo=62 (55% native) |
