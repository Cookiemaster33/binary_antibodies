# Active Lambda Cloud Instance — FULL PIPELINE RUNNING (with structures)

| Field | Value |
|---|---|
| Instance ID | `59a43a5ff4c647b99baeb366a74876ac` |
| IP | `129.146.52.56` |
| Region | `us-west-2` |
| Status | **Running — Step 0: Docker pull + Boltz-2 install** |

## This run adds structure files

In addition to JSON scores, this run will push to GitHub:
- `top01_*_rfd3_backbone.cif` through `top10_*_rfd3_backbone.cif`
- `top01_*_boltz2_complex.cif` through `top10_*_boltz2_complex.cif`
- `top10_summary.tsv` — table with all metrics

Ranked by √(ipTM_VH1 × ipTM_Nb) — designs where Boltz-2 is most confident
about BOTH interfaces simultaneously.

## Expected timeline (~90 min total)

| Step | Time |
|---|---|
| Docker + Boltz-2 install | ~5 min |
| RFdiffusion3 (200 designs) | ~20 min |
| ProteinMPNN | ~10 min |
| Boltz-2 (50 complexes) | ~45 min |
| scRMSD + top-10 CIF collection | ~3 min |
| GitHub push | ~1 min |

## Results will appear in the PR automatically when done

👉 https://github.com/Cookiemaster33/binary_antibodies/pull/1

## Terminate when done

```bash
curl -X POST https://cloud.lambda.ai/api/v1/instance-operations/terminate \
  -H "Authorization: Bearer $LAMBDA_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"instance_ids": ["59a43a5ff4c647b99baeb366a74876ac"]}'
```
