# Active Lambda Instance — v5 Two-Round (relaunch)

| Field | Value |
|---|---|
| Instance ID | `d83660305eea427ba4664752954b13f8` |
| IP | `150.136.116.198` |
| Status | **Running — RFd3 round 1 in progress** |
| Results folder | `pipeline_results/v5_two_round_refine/` |

## What happened to the first instance?

Instance `5314d4de` terminated because the pipeline **crashed within seconds** of starting:
- `cp` tried to copy `rfd3_round2.py` onto itself → exit code 1
- `set -e` in the shell script treated that as fatal
- Launch script saw tmux session die → auto-terminated the instance

**Fixed:** skip self-copy, keep instance alive on failure, relaunched.

## Monitor

```bash
tail -f /workspace/pipeline_launch.log

ssh -i ~/.ssh/lambda_agent_key ubuntu@150.136.116.198 \
  "tail -f ~/pipeline/full_pipeline.log"
```
