# Active Lambda Instance — Stage 0 VH–VL Interface

| Field | Value |
|---|---|
| Instance ID | `5400b4d26ae3488e9f418f7666d9268e` |
| IP | `129.213.129.28` |
| Status | **Running — Stage 0 RFd3 in progress** |
| Pipeline | RFd3 → MPNN → Boltz apo/holo → differential scoring |
| Designs | 200 (partial_t=12 Å) |

**Note:** First two launches failed due to epitope stub PDB issues (fixed: full-atom, three-letter residue names).

## Monitor

```bash
ssh -i ~/.ssh/cursor_lambda_ephemeral ubuntu@129.213.129.28 \
  "tail -f /home/ubuntu/pipeline/stage_0_pipeline.log"

ssh -i ~/.ssh/cursor_lambda_ephemeral ubuntu@129.213.129.28 \
  "tmux attach -t stage_0"
```

## Terminate

```bash
python3 scripts/launch_stage_0.py --terminate --instance-id 5400b4d26ae3488e9f418f7666d9268e
```
