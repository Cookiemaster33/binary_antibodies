# Active Lambda Instance — Stage 0 VH–VL Interface

| Field | Value |
|---|---|
| Instance ID | `bbf9df5599e849fda5e5fc7e82608ed7` |
| IP | `150.136.216.155` |
| Status | **Running — Stage 0 full pipeline** |
| Pipeline | RFd3 (partial_t=12) → MPNN → Boltz apo/holo → scoring |
| Designs | 200 |

## Monitor

```bash
ssh -i ~/.ssh/cursor_lambda_ephemeral ubuntu@150.136.216.155 \
  "tail -f /home/ubuntu/pipeline/stage_0_pipeline.log"

ssh -i ~/.ssh/cursor_lambda_ephemeral ubuntu@150.136.216.155 \
  "tmux attach -t stage_0"
```

## Terminate

```bash
python3 scripts/launch_stage_0.py --terminate --instance-id bbf9df5599e849fda5e5fc7e82608ed7
```

## After completion

Top hits from `final/stage_0_results.json` → feed into Stage A:

```bash
python scripts/build_stage_a_design_target.py --fab-pdb <top_holo.cif>
python scripts/launch_stage_a.py
```
