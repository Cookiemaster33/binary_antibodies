# Active Lambda Cloud Instance — PIPELINE RUNNING

| Field | Value |
|---|---|
| Instance ID | `28b31e6dafc845d48ebfd92931f1ee5b` |
| IP | `132.145.135.76` |
| Type | `gpu_1x_a100_sxm4` (40 GB VRAM) |
| Region | `us-east-1` |
| SSH key name | `cursor-agent` |
| Status | **RFdiffusion3 running** — 3/20 batches done, ~22 CIF files generated |

## Pipeline

RFdiffusion3 (rc-foundry Docker) + ProteinMPNN, running in tmux session `design`.

```
Spec:  contig="B1-130/0 60"   (fix VL1, design 60-res minibinder)
       select_hotspots="B35-39,B44-47,B98"  (VL1 FR2 face)
N_DESIGNS=200  BATCH_SIZE=10  N_MPNN_SEQS=8
```

## Monitor

```bash
ssh -i $LAMBDA_SSH_KEY ubuntu@132.145.135.76 "tail -f ~/pipeline/run.log"
ssh -i $LAMBDA_SSH_KEY ubuntu@132.145.135.76 "ls ~/pipeline/outputs/rfd3/*.cif | wc -l"
```

## Download results when done

```bash
scp -r -i $LAMBDA_SSH_KEY ubuntu@132.145.135.76:~/pipeline/outputs ./pipeline_results
```

## Terminate instance (DO NOT FORGET — $1.99/hr)

```bash
curl -X POST https://cloud.lambda.ai/api/v1/instance-operations/terminate \
  -H "Authorization: Bearer $LAMBDA_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"instance_ids": ["28b31e6dafc845d48ebfd92931f1ee5b"]}'
```
