# Active Lambda Cloud Instance — PIPELINE v2 RUNNING

| Field | Value |
|---|---|
| Instance ID | `95d1d0703932446c82dff3d34648532b` |
| IP | `129.146.79.157` |
| Region | `us-west-2` |
| Status | **Step 1/4: RFdiffusion3 running — batch ~2/20** |

## Pipeline v2 improvements

1. **4-chain RFd3 design**: VH1(A) + VL(B, steric context) + Nanobody(C)
   - VL blocks the VH-VL pairing face during design — minibinder CANNOT clash with VL
   - Removed 3 overlapping hotspots (39,41,85): they contact BOTH CH1 face AND VL

2. **Single-chain Boltz-2 validation**: folds VH1+Minibinder+Nanobody as ONE polypeptide
   - Same topology as the RFd3 output (no chain separation)
   - scRMSD < 2 Å is now meaningful — the chain is connected so topology is preserved

## ETA

~90 min total. Results push automatically to GitHub when done.

👉 https://github.com/Cookiemaster33/binary_antibodies/pull/1

## Terminate when done

```bash
curl -X POST https://cloud.lambda.ai/api/v1/instance-operations/terminate \
  -H "Authorization: Bearer $LAMBDA_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"instance_ids": ["95d1d0703932446c82dff3d34648532b"]}'
```
