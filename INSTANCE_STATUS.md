# Active Lambda Cloud Instance — VARIABLE-LENGTH + FIXED-ATOMS RUN

| Field | Value |
|---|---|
| Instance ID | `4d240eaac41b400aa781b393583fca80` |
| IP | `129.146.177.146` |
| Status | **Step 1/4: RFdiffusion3 batch ~2/20** |

## New in this run

1. **Variable-length minibinder** (`MB_LENGTH_RANGE=35-70`)
   RFd3 samples a length from 35–70 residues per design.
   Shorter designs → fewer DOF → expected lower scRMSD.

2. **`select_fixed_atoms` on chains A, B, C**
   VH1, VL, Nanobody frozen at exact input coordinates (hard pin).
   Previously `select_hotspots` was soft — chains could drift.

## Monitor

```bash
ssh -i $LAMBDA_SSH_KEY ubuntu@129.146.177.146 \
  "grep -v 'WARNING\|Cached\|MACE\|not set\|bashrc\|networkx' \
   ~/pipeline/full_pipeline.log | tail -5; \
   echo CIFs: \$(ls ~/pipeline/outputs/rfd3/*.cif 2>/dev/null | wc -l)/200"
```

## Terminate when done

```bash
curl -X POST https://cloud.lambda.ai/api/v1/instance-operations/terminate \
  -H "Authorization: Bearer $LAMBDA_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"instance_ids": ["4d240eaac41b400aa781b393583fca80"]}'
```
