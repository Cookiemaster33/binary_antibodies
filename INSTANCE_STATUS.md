# Active Lambda Cloud Instance — FULL INTEGRATED PIPELINE RUNNING

| Field | Value |
|---|---|
| Instance ID | `52d978c10b4543a7b54ba965a3ea09b9` |
| IP | `129.146.164.146` |
| Region | `us-west-2` |
| Status | **Step 0: Docker pull + Boltz-2 install (parallel)** |

## What will run (all on this instance)

```
Step 1: RFdiffusion3  — 200 bispecific bridging minibinder designs
                         Contig: A1-115,55,B1-115 (VH1+MB+Nanobody)
                         CIFs saved to ~/pipeline/outputs/rfd3/

Step 2: ProteinMPNN   — 8 sequences per backbone (1600 total)

Step 3: Boltz-2       — top-50 unique backbones as 3-chain complexes
                         Pocket constraints on VH1-CH1-face + Nanobody CDRs

Step 4: scRMSD        — superimpose Boltz-2 on VH1+Nb (fixed), measure
                         RMSD of Boltz-2 minibinder vs RFd3 backbone
                         CIFs are available because same instance!

Filter: scRMSD < 2Å AND pLDDT > 60 AND ipTM(MB↔VH1) > 0.3 AND ipTM(MB↔Nb) > 0.3
```

## Estimated time

| Step | Time |
|---|---|
| Docker + Boltz-2 install | ~5 min (parallel) |
| RFdiffusion3 (200 designs) | ~20 min |
| ProteinMPNN | ~10 min |
| Boltz-2 (50 complexes) | ~45 min |
| scRMSD scoring | ~2 min |
| **Total** | **~80 min, ~$2.65** |

## Monitor

```bash
ssh -i $LAMBDA_SSH_KEY ubuntu@129.146.164.146 \
  "grep -v 'WARNING\|Cached\|MACE\|not set\|bashrc\|AMP\|Tensor\|DEBUG' \
   ~/pipeline/full_pipeline.log | tail -10; \
   echo CIFs: \$(ls ~/pipeline/outputs/rfd3/*.cif 2>/dev/null | wc -l)"
```

## Terminate when done

```bash
curl -X POST https://cloud.lambda.ai/api/v1/instance-operations/terminate \
  -H "Authorization: Bearer $LAMBDA_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"instance_ids": ["52d978c10b4543a7b54ba965a3ea09b9"]}'
```
