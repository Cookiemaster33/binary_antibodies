# Active Lambda Cloud Instance

| Field | Value |
|---|---|
| Instance ID | `28b31e6dafc845d48ebfd92931f1ee5b` |
| IP | `132.145.135.76` |
| Type | `gpu_1x_a100_sxm4` (40 GB VRAM) |
| Region | `us-east-1` |
| SSH key name | `cursor-agent` |
| Started | 2026-07-04 ~02:10 UTC |
| Status | **Setting up pipeline** (step 1/6: system packages) |

## Monitor progress

```bash
# From an agent with LAMBDA_SSH_KEY secret set:
ssh -i $LAMBDA_SSH_KEY ubuntu@132.145.135.76 "tail -f ~/pipeline/setup_pipeline.log"

# Or check if setup is done:
ssh -i $LAMBDA_SSH_KEY ubuntu@132.145.135.76 "grep -c 'Setup complete' ~/pipeline/setup_pipeline.log"
```

## What runs after setup (~20 min)

The `run_minibinder_design.sh` script will be launched automatically in a
`design` tmux session and will:
1. RFdiffusion — 200 designs against VL1 FR2 hotspot (B35-39, B44-47, B98)
2. ProteinMPNN — 8 sequences × 200 backbones
3. AF2-Multimer — score top 20 by ipTM

Results will be in `~/pipeline/outputs/` on the instance.

## To add LAMBDA_SSH_KEY secret

Go to cursor.com/settings → Cloud Agents → Secrets → add `LAMBDA_SSH_KEY`
with the value of the private key printed during setup.

## To terminate (DO NOT FORGET — costs $1.99/hr)

```bash
curl -X POST https://cloud.lambda.ai/api/v1/instance-operations/terminate \
  -H "Authorization: Bearer $LAMBDA_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"instance_ids": ["28b31e6dafc845d48ebfd92931f1ee5b"]}'
```
