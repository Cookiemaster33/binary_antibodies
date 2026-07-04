# Instance Running — Results Will Auto-Push to GitHub

| Field | Value |
|---|---|
| Instance ID | `52d978c10b4543a7b54ba965a3ea09b9` |
| IP | `129.146.164.146` |
| Status | **Boltz-2 predicting 50 complexes (~45 min)** |

## Auto-push configured ✓

When Boltz-2 + scRMSD scoring finish, results will be **pushed automatically
to this GitHub PR** — no action needed from you.

**Where to see results:**
👉 https://github.com/Cookiemaster33/binary_antibodies/pull/1

Look for a new commit titled "results: Boltz-2 + scRMSD validation complete"
in the PR. Files will appear at:
- `pipeline_results/v3_bispecific_validated/final_results.json`
- `pipeline_results/v3_bispecific_validated/validated_minibinders.fasta`

## Terminate the instance when done

After results appear in the PR, terminate to avoid charges:

```bash
curl -X POST https://cloud.lambda.ai/api/v1/instance-operations/terminate \
  -H "Authorization: Bearer $LAMBDA_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"instance_ids": ["52d978c10b4543a7b54ba965a3ea09b9"]}'
```

Or just ask a new Cursor agent: "terminate the Lambda instance 52d978c10b4543a7b54ba965a3ea09b9"
