#!/usr/bin/env bash
# ============================================================
# run_peptide_cftr_boltz2.sh
# --------------------------
# Boltz-2 templated binding prediction: peptide LGQQQPFPPQQPY
# vs CFTR structures 5UAK and 6O2P.
#
# Input:
#   $PIPELINE/yaml/            — Boltz-2 YAML inputs (2 files)
#   $PIPELINE/templates/       — CIF template files (5UAK.cif, 6O2P.cif)
# Output:
#   $PIPELINE/boltz_results/   — Boltz-2 predictions
#   $PIPELINE/results/         — ipTM summary
# ============================================================
set -eo pipefail

PIPELINE=/home/ubuntu/pipeline
YAML_DIR=$PIPELINE/yaml
TEMPLATE_DIR=$PIPELINE/templates
RESULTS_DIR=$PIPELINE/boltz_results
SUMMARY_DIR=$PIPELINE/results

mkdir -p "$YAML_DIR" "$TEMPLATE_DIR" "$RESULTS_DIR" "$SUMMARY_DIR"

echo "===== Peptide-CFTR Boltz-2 screen: $(date) ====="
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo n/a)"

N_YAML=$(ls "$YAML_DIR"/*.yaml 2>/dev/null | wc -l)
echo "YAML inputs: $N_YAML"
if [ "$N_YAML" -eq 0 ]; then echo "ERROR: no YAML files in $YAML_DIR"; exit 1; fi


# ── 1. Install Boltz-2 in an isolated venv ───────────────────
echo ""
echo "=== [1/3] Boltz-2 setup ==="
BOLTZ_VENV="$HOME/boltz_venv"
BOLTZ="$BOLTZ_VENV/bin/boltz"

if [ ! -f "$BOLTZ" ]; then
    echo "  Creating venv..."
    python3 -m venv "$BOLTZ_VENV"
    PIP="$BOLTZ_VENV/bin/pip"
    "$PIP" install -q --upgrade pip 2>&1 | tail -2
    # Install PyTorch for CUDA 12.4 (compatible with Lambda driver)
    "$PIP" install -q torch torchvision \
        --index-url https://download.pytorch.org/whl/cu124 2>&1 | tail -4
    # Install boltz base (no cuda extras to avoid cuequivariance_ops_torch issue)
    "$PIP" install -q boltz 2>&1 | tail -4
    # Remove cuequivariance-torch if boltz pulled it in — the compiled ops
    # extension is not available for this driver; force pure-PyTorch fallback
    "$PIP" uninstall -y cuequivariance-torch cuequivariance cuequivariance-ops 2>/dev/null || true
    echo "  Boltz-2 venv ready."
fi

"$BOLTZ_VENV/bin/python3" -c "
import torch
print(f'  torch {torch.__version__}  CUDA={torch.cuda.is_available()}'
      f'  device={torch.cuda.get_device_name(0) if torch.cuda.is_available() else \"n/a\"}')
"


# ── 2. Run Boltz-2 templated predictions ─────────────────────
echo ""
echo "=== [2/3] Boltz-2 predictions ==="

"$BOLTZ" predict "$YAML_DIR" \
    --out_dir "$RESULTS_DIR" \
    --devices 1 \
    --num_workers 2 \
    --override \
    2>&1 | tee "$PIPELINE/boltz_run.log"

echo "  Predictions complete."


# ── 3. Extract and report ipTM scores ────────────────────────
echo ""
echo "=== [3/3] ipTM summary ==="

python3 - << 'PYEOF'
import json
from pathlib import Path

results = Path("/home/ubuntu/pipeline/boltz_results")
summary = Path("/home/ubuntu/pipeline/results")
summary.mkdir(exist_ok=True)

scores = []
for f in sorted(results.rglob("confidence_*.json")):
    try:
        data = json.loads(f.read_text())
        name = f.parent.parent.name  # boltz_results/<name>/predictions/

        raw_iptm = data.get("iptm", {})
        if isinstance(raw_iptm, dict):
            # Chain-pair ipTM: AB = peptide vs receptor
            iptm_AB = float(raw_iptm.get("AB", 0.0))
            iptm_all = {k: round(float(v), 3) for k, v in raw_iptm.items()}
        else:
            iptm_AB = float(raw_iptm)
            iptm_all = {"overall": iptm_AB}

        ptm   = float(data.get("ptm",          0.0))
        plddt = float(data.get("mean_plddt",   data.get("complex_plddt", 0.0)))

        scores.append({
            "name":     name,
            "iptm_AB":  round(iptm_AB, 3),
            "iptm_all": iptm_all,
            "ptm":      round(ptm,   3),
            "plddt":    round(plddt, 1),
        })
    except Exception as e:
        print(f"  skip {f}: {e}")

scores.sort(key=lambda x: x["iptm_AB"], reverse=True)

print(f"\n{'Name':<30} {'ipTM(A↔B)':>10} {'pTM':>6} {'pLDDT':>7}")
print("-" * 58)
for s in scores:
    flag = " ★" if s["iptm_AB"] > 0.5 else ""
    print(f"{s['name']:<30} {s['iptm_AB']:>10.3f} {s['ptm']:>6.3f} {s['plddt']:>7.1f}{flag}")
    for pair, val in s["iptm_all"].items():
        if pair != "AB":
            print(f"  {'ipTM(' + pair + ')':>28}: {val:.3f}")

import json as _json
(_json.loads("{}"))  # noop
(summary / "iptm_scores.json").write_text(json.dumps(scores, indent=2))
print(f"\nSaved: {summary}/iptm_scores.json")
PYEOF

echo ""
echo "===== Done: $(date) ====="
