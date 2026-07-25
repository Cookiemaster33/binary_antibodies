# Active Lambda Instance — Stage A PoC (minibinder-first)

| Field | Value |
|---|---|
| Instance ID | `0926bc09188c4cebb0cd8872131e29cb` |
| IP | `150.136.64.212` |
| Status | **Running — Stage A PoC RFd3 (native Fab, exploratory)** |
| Production order | **Stage 0 → Stage A** (see `docs/HIDDEN_MINIBINDER_PIPELINE.md`) |

> This run uses native trastuzumab VH–VL as a feasibility check. After it finishes,
> launch Stage 0 first, then rebuild Stage A with `--fab-pdb` from the top Stage 0 hit.

## What happened to the previous instance?

Instance `76af49c4...` @ `150.136.119.26` was terminated. Its cloud-init bootstrap never started RFd3 (likely failed on `build_stage_a_design_target.py` — BioPython not installed on a fresh Lambda image). Relaunched with SSH-based deploy (same pattern as the v6 full pipeline).

## Monitor (cloud agent ephemeral key)

```bash
ssh -i ~/.ssh/cursor_lambda_ephemeral ubuntu@150.136.64.212 \
  "tail -f /home/ubuntu/pipeline/stage_a_pipeline.log"

ssh -i ~/.ssh/cursor_lambda_ephemeral ubuntu@150.136.64.212 \
  "tmux attach -t stage_a"

ssh -i ~/.ssh/cursor_lambda_ephemeral ubuntu@150.136.64.212 \
  "nvidia-smi"
```

If you use your local `cursor-agent` key, it will **not** work on this instance — only the ephemeral key registered at launch is authorized.

## Terminate when done

```bash
python3 scripts/launch_stage_a.py --terminate --instance-id 0926bc09188c4cebb0cd8872131e29cb
```
