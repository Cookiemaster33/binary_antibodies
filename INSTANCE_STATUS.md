# Active Lambda Cloud Instance — 500-DESIGN BRUTE FORCE RUN

| Field | Value |
|---|---|
| Instance ID | `71d4883603a9499983f6c10a89675df5` |
| IP | `129.146.97.5` |
| Region | `us-west-2` |
| Status | **Step 1/4: RFdiffusion3 — batch ~3/50** |

## Run parameters

| Parameter | Value |
|---|---|
| N_DESIGNS | **500** (was 200) |
| TOP_N for Boltz-2 | **100** (was 50) |
| Everything else | same as v2 (4-chain RFd3, single-chain Boltz-2, fixed scoring) |

## Timeline (~165 min total)

| Step | Time |
|---|---|
| RFdiffusion3 (500 designs) | ~50 min |
| ProteinMPNN (4000 sequences) | ~25 min |
| Boltz-2 (100 complexes) | ~90 min |
| scRMSD scoring | ~3 min |

## Monitor

```bash
ssh -i $LAMBDA_SSH_KEY ubuntu@129.146.97.5 \
  "grep -v 'WARNING\|Cached\|MACE\|not set\|bashrc\|AMP\|Tensor\|networkx' \
   ~/pipeline/full_pipeline.log | tail -5; \
   echo CIFs: \$(ls ~/pipeline/outputs/rfd3/*.cif 2>/dev/null | wc -l)/500"
```

## Terminate

```bash
curl -X POST https://cloud.lambda.ai/api/v1/instance-operations/terminate \
  -H "Authorization: Bearer $LAMBDA_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"instance_ids": ["71d4883603a9499983f6c10a89675df5"]}'
```
