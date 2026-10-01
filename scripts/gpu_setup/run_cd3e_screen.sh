#!/usr/bin/env bash
# Run Boltz-1 multimer screen: top Stage A minibinders vs CD3e ECD
set -euo pipefail

PIPELINE=/home/ubuntu/pipeline
SCREEN_DIR=$PIPELINE/cd3e_screen
RESULTS_DIR=$PIPELINE/cd3e_screen_results

mkdir -p "$SCREEN_DIR" "$RESULTS_DIR"

N_YAML=$(ls "$SCREEN_DIR"/*.yaml 2>/dev/null | wc -l)
echo "=== CD3e screen: $N_YAML complexes ==="
echo "$(date)"

if [ "$N_YAML" -eq 0 ]; then
    echo "ERROR: no YAML files in $SCREEN_DIR"; exit 1
fi

# Run Boltz-1 (no MSA server — both chains have msa: empty)
$HOME/.local/bin/boltz predict "$SCREEN_DIR" \
    --out_dir "$RESULTS_DIR" \
    --accelerator gpu \
    --devices 1 \
    --override \
    2>&1 | tee "$PIPELINE/cd3e_screen.log"

echo ""
echo "=== CD3e screen complete: $(date) ==="

# Extract ipTM scores
echo ""
echo "=== ipTM summary ==="
python3 - << 'PYEOF'
import json
from pathlib import Path

results = Path("/home/ubuntu/pipeline/cd3e_screen_results")
scores = []
for f in sorted(results.glob("**/confidence_*.json")):
    try:
        data = json.loads(f.read_text())
        name = f.parent.parent.name   # boltz_results/<name>/predictions/
        iptm = data.get("iptm", data.get("complex_plddt", 0.0))
        ptm  = data.get("ptm", 0.0)
        scores.append((iptm, ptm, name))
    except Exception as e:
        print(f"  skip {f}: {e}")

scores.sort(reverse=True)
print(f"{'rank':<6} {'ipTM':>6} {'pTM':>6}  design")
for i, (iptm, ptm, name) in enumerate(scores, 1):
    print(f"{i:<6} {iptm:>6.3f} {ptm:>6.3f}  {name}")
PYEOF
