# Active Lambda Cloud Instance — BOLTZ-2 RUNNING

| Field | Value |
|---|---|
| Instance ID | `52d978c10b4543a7b54ba965a3ea09b9` |
| IP | `129.146.164.146` |
| Region | `us-west-2` |
| tmux session | `boltz2` |
| Status | **Boltz-2 predicting 50 complexes (~45 min remaining)** |

## What's already done on this instance

- ✅ RFdiffusion3: 200 CIFs in `~/pipeline/outputs/rfd3/`
- ✅ ProteinMPNN: 1600 sequences in `~/pipeline/outputs/mpnn/`
- ✅ Boltz-2 YAMLs: 50 inputs in `~/pipeline/boltz_inputs/`
- 🔄 Boltz-2 prediction: running now in tmux session `boltz2`
- ⏳ scRMSD scoring: will run automatically after Boltz-2 finishes

## What to do when done

```bash
# Download results
scp -r -i $LAMBDA_SSH_KEY ubuntu@129.146.164.146:~/pipeline/final ./results_final

# Terminate instance
curl -X POST https://cloud.lambda.ai/api/v1/instance-operations/terminate \
  -H "Authorization: Bearer $LAMBDA_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"instance_ids": ["52d978c10b4543a7b54ba965a3ea09b9"]}'
```

## Monitor

```bash
ssh -i $LAMBDA_SSH_KEY ubuntu@129.146.164.146 \
  "grep -v 'WARNING\|Cached\|MACE\|not set\|bashrc\|AMP\|Tensor\|networkx' \
   ~/pipeline/full_pipeline.log | tail -8 ; \
   echo Boltz CIFs: \$(find ~/pipeline/boltz_outputs -name '*.cif' 2>/dev/null | wc -l)/50"
```

## Results location (after completion)

| File | Content |
|---|---|
| `~/pipeline/final/final_results.json` | scRMSD + pLDDT + ipTM for all 50 designs |
| `~/pipeline/final/validated_minibinders.fasta` | Passing designs (scRMSD<2Å, pLDDT>60, both ipTM>0.3) |
