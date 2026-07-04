# Active Lambda Cloud Instance — PIPELINE RUNNING

| Field | Value |
|---|---|
| Instance ID | `790ebff8788b4ac8814822b5f0fc419b` |
| IP | `150.136.66.45` |
| Type | `gpu_1x_a100_sxm4` (40 GB) |
| Region | `us-east-1` |
| SSH key | `cursor-agent` |
| Status | **RFdiffusion3 running** — ~20/200 designs done |

## What's running

Anti-idiotypic minibinder design against **2Rs15d nanobody CDR face** (CH1-kicker design).

```
Input:     her2_nanobody_VHH.pdb (chain B, 115 residues)
Contig:    B1-115/0,50  →  fix nanobody, design 50-res minibinder
Hotspots:  CDR1 (B26-33) + CDR2 (B50-57) + CDR3 (B97-112)
Output:    CIF files with 165 res (1-115=nanobody, 116-165=minibinder)
```

Pipeline: RFdiffusion3 (200 designs) → ProteinMPNN (8 seqs/backbone) → `minibinder_sequences.fasta`

## Monitor

```bash
ssh -i $LAMBDA_SSH_KEY ubuntu@150.136.66.45 \
  "grep -v WARNING ~/pipeline/run.log | tail -5; \
   echo 'CIF count:' \$(ls ~/pipeline/outputs/rfd3/*.cif 2>/dev/null | wc -l)"
```

## Terminate when done

```bash
curl -X POST https://cloud.lambda.ai/api/v1/instance-operations/terminate \
  -H "Authorization: Bearer $LAMBDA_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"instance_ids": ["790ebff8788b4ac8814822b5f0fc419b"]}'
```
