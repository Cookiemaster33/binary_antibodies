# Active Lambda Instance — Stage A Hidden Minibinder

| Field | Value |
|---|---|
| Instance ID | `76af49c4b04046ea94ae0d97d84b1402` |
| IP | `150.136.119.26` |
| Status | **Running — cloud-init bootstrap → RFd3 Stage A** |
| Branch | `cursor/conditional-nanobody-design-992c` |
| Designs | 200 (VH–hub–VL, CH1+VL hotspots) |

## Bootstrap (automatic)

Cloud-init clones repo, builds design target, runs `setup_pipeline_rfd3.sh`, then Stage A RFd3 in tmux session `stage_a`.

## Monitor (from machine with Lambda SSH key)

```bash
ssh -i ~/.ssh/lambda_agent_key ubuntu@150.136.119.26 \
  "tail -f /home/ubuntu/stage_a_bootstrap.log"

ssh -i ~/.ssh/lambda_agent_key ubuntu@150.136.119.26 \
  "tail -f /home/ubuntu/pipeline/stage_a_pipeline.log"

ssh -i ~/.ssh/lambda_agent_key ubuntu@150.136.119.26 \
  "tmux attach -t stage_a"
```

## Terminate when done

```bash
python3 scripts/launch_full_pipeline.py --terminate --instance-id 76af49c4b04046ea94ae0d97d84b1402
```
