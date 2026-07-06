# Active Lambda Instance — v6 Global RMSD Rerun

| Field | Value |
|---|---|
| Instance ID | `142cd0e053a947b093ef23dc2b18e051` |
| IP | `150.136.146.146` |
| Status | **Running — RFd3 round 1 in progress** |
| Results folder | `pipeline_results/v6_global_rmsd_rerun` |
| Metric | Global assembly RMSD (VH1 + MB + Nb) for round 1 & 2 ranking |

## Monitor

```bash
tail -f /workspace/pipeline_launch_v6.log

ssh -i ~/.ssh/lambda_agent_key ubuntu@150.136.146.146 \
  "tail -f ~/pipeline/full_pipeline.log"
```

Launch script will auto-terminate the instance on successful completion.
