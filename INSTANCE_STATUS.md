# Active Lambda Cloud Instance — BOLTZ-2 RUNNING

| Field | Value |
|---|---|
| Instance ID | `94c63717df0a4432b814247acbd82dd2` |
| IP | `129.146.164.146` |
| Region | `us-west-2` |
| Status | **Boltz-2 downloading data, then predicting 50 complexes** |

## What's running

Boltz-2 v2.2.1 predicting 50 three-chain complexes:
- Chain A: VH1 (Trastuzumab VH, 115 res)
- Chain B: Designed minibinder (56 res, from MPNN)
- Chain C: 2Rs15d anti-HER2 nanobody (115 res)

**Pocket constraints** guide Boltz-2 to the designed binding surfaces:
- VH1-CH1-face residues (28 contacts on chain A)
- Nanobody CDR residues (27 contacts on chain C)

**Scoring**: pLDDT + chain-pair ipTM(A,B) + ipTM(B,C)
Filter: pLDDT > 60 AND ipTM(MB↔VH1) > 0.3 AND ipTM(MB↔Nb) > 0.3

## Note on scRMSD

Full scRMSD (Boltz-2 vs RFd3 backbone) requires CIFs from the design run
on the same instance. Next time: use integrated pipeline that does
RFd3 → MPNN → Boltz-2 all on one instance (CIFs never lost).

## Monitor

```bash
ssh -i $LAMBDA_SSH_KEY ubuntu@129.146.164.146 \
  "tail -10 ~/pipeline/boltz2.log | grep -E 'rank|pLDDT|Passed|Error'"
```

## Terminate

```bash
curl -X POST https://cloud.lambda.ai/api/v1/instance-operations/terminate \
  -H "Authorization: Bearer $LAMBDA_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"instance_ids": ["94c63717df0a4432b814247acbd82dd2"]}'
```
