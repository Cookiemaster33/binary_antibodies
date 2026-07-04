# Active Lambda Cloud Instance — MPNN RUNNING

| Field | Value |
|---|---|
| Instance ID | `28b31e6dafc845d48ebfd92931f1ee5b` |
| IP | `132.145.135.76` |
| SSH key name | `cursor-agent` |
| Status | **ProteinMPNN running** — 200 RFD3 designs complete, MPNN sequencing |

## What was fixed

Original run used wrong contig `"B1-130/0 60"` → output was only the 130-res VL1 target.
Fixed contig: `"B1-130/0,60"` → output is 190 residues (1-130 = VL1, 131-190 = designed binder).

## Pipeline

- RFD3: 200 × 190-residue designs (VL1 + 60-res binder) ✓ DONE
- MPNN: 8 sequences per backbone, binder portion only (res 131-190), VL1 fixed ← RUNNING
- Results: `~/pipeline/outputs/mpnn/binder_sequences.fasta`

## Monitor

```bash
ssh -i $LAMBDA_SSH_KEY ubuntu@132.145.135.76 "tail -5 ~/pipeline/run.log"
```

## Terminate when done

```bash
curl -X POST https://cloud.lambda.ai/api/v1/instance-operations/terminate \
  -H "Authorization: Bearer $LAMBDA_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"instance_ids": ["28b31e6dafc845d48ebfd92931f1ee5b"]}'
```
