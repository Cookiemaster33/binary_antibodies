# Active Lambda Cloud Instance — INTEGRATED PIPELINE RUNNING

| Field | Value |
|---|---|
| Instance ID | `b8d850dd4e184f1d829754f1413cc4c6` |
| IP | `129.146.33.224` |
| Region | `us-west-2` |
| Type | `gpu_1x_a100_sxm4` |
| SSH key | `cursor-agent` |
| Status | **Step 1/4: RFdiffusion3 running** |

## What's running

Integrated pipeline (RFd3 → MPNN → AF2 self-consistency), all on one instance.

**Design context: 3-chain (new)**
- Chain B: 2Rs15d nanobody (115 res) — hotspot: CDR1+CDR2+CDR3
- Chain C: Trastuzumab CL domain (113 res) — kicker context
- Contig: `B1-115,50,C1-113` → 278-res output (115+50+113)
- Minibinder (residues 116-165) bridges nanobody CDR face AND CL domain

**Validation: self-consistency (Step 4)**
- AF2 monomer on 50-res minibinder sequence alone
- scRMSD vs RFd3 backbone — filter: scRMSD < 2 Å, pLDDT > 60

## Monitor

```bash
ssh -i $LAMBDA_SSH_KEY ubuntu@129.146.33.224 \
  "grep -v 'WARNING\|Cached\|MACE\|not set\|bashrc' ~/pipeline/pipeline.log | tail -10"
```

## Terminate

```bash
curl -X POST https://cloud.lambda.ai/api/v1/instance-operations/terminate \
  -H "Authorization: Bearer $LAMBDA_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"instance_ids": ["b8d850dd4e184f1d829754f1413cc4c6"]}'
```
