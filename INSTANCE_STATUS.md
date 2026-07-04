# Active Lambda Cloud Instance — BISPECIFIC MINIBINDER RUNNING

| Field | Value |
|---|---|
| Instance ID | `663bf9f3a1764ad0beecb4eec9ad51a8` |
| IP | `161.153.122.32` |
| Region | `us-west-2` |
| Status | **RFdiffusion3 running** — batch ~14/20 |

## Design: Bispecific bridging minibinder

Input: `vh1_nanobody_design_target.pdb`
- Chain A: VH1 (115 res) — hotspot: VH1-CH1-contact face (FR1+FR2+FR3)
- Chain B: Nanobody positioned in CH1 slot (115 res) — hotspot: CDR1+CDR2+CDR3
- Contig: `A1-115,55,B1-115` → 285 res total

Minibinder (residues 116-170) bridges:
- VH1-CH1-face (where CL/VL1 will dock upon activation) 
- Nanobody CDR face (CDRs it must block in OFF state)

## Monitor

```bash
ssh -i $LAMBDA_SSH_KEY ubuntu@161.153.122.32 \
  "grep -v 'WARNING\|Cached\|MACE\|not set\|bashrc' ~/pipeline/pipeline.log | tail -5; \
   echo CIFs: \$(ls ~/pipeline/outputs/rfd3/*.cif 2>/dev/null | wc -l)"
```

## Terminate

```bash
curl -X POST https://cloud.lambda.ai/api/v1/instance-operations/terminate \
  -H "Authorization: Bearer $LAMBDA_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"instance_ids": ["663bf9f3a1764ad0beecb4eec9ad51a8"]}'
```
