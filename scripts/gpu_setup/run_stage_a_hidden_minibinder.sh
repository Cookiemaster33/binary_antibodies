#!/bin/bash
# ============================================================
# run_stage_a_hidden_minibinder.sh
# Stage A: RFd3 design of unlinked minibinder bridging CH1 + VL in full Fab context.
#
# Hotspots: CH1 (chain C) + VL (chain B)
# Fixed context: full Fab A–D + HER2 epitope stub (T)
# Minibinder: separate chain (not linked to VH/VL); linker added later
#
# Prerequisite (run locally before GPU launch):
#   python scripts/build_stage_a_design_target.py
# ============================================================
set -eo pipefail

PIPELINE=${PIPELINE_DIR:-/home/ubuntu/pipeline}
LOG=$PIPELINE/stage_a_pipeline.log
exec > >(tee -a "$LOG") 2>&1

N_DESIGNS=${N_DESIGNS:-200}
BATCH_SIZE=${BATCH_SIZE:-10}
N_MPNN_SEQS=${N_MPNN_SEQS:-8}
TOP_N=${TOP_N:-50}
MB_LENGTH_RANGE=${MB_LENGTH_RANGE:-35-55}
N_BATCHES=$(( (N_DESIGNS + BATCH_SIZE - 1) / BATCH_SIZE ))

INPUT_PDB_NAME=${INPUT_PDB:-fab_hidden_minibinder_stage_a.pdb}
CONFIG_JSON=${CONFIG_JSON:-stage_a_hidden_minibinder_config.json}

echo "===== Stage A hidden minibinder pipeline: $(date) ====="
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "Designs: $N_DESIGNS | MB length: $MB_LENGTH_RANGE"

mkdir -p "$PIPELINE/inputs" "$PIPELINE/outputs/rfd3_stage_a" \
         "$PIPELINE/outputs/mpnn_stage_a" "$PIPELINE/boltz_inputs_stage_a" \
         "$PIPELINE/boltz_outputs_stage_a" "$PIPELINE/final"

sudo chmod 666 /var/run/docker.sock 2>/dev/null || true
docker pull rosettacommons/foundry:latest > "$PIPELINE/docker_pull.log" 2>&1 || true

# Copy design target + config into pipeline inputs
cp "$PIPELINE/inputs/$INPUT_PDB_NAME" "$PIPELINE/inputs/$INPUT_PDB_NAME" 2>/dev/null || true

cat > "$PIPELINE/run_rfd3_stage_a.py" << 'PYEOF'
import json, os, torch
from pathlib import Path
from atomworks.io.utils.io_utils import to_cif_file
from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine
from rfd3.inference.input_parsing import DesignInputSpecification

torch.set_float32_matmul_precision("high")

PIPELINE = Path("/workspace")
INPUT  = PIPELINE / "inputs" / os.environ.get("INPUT_PDB", "fab_hidden_minibinder_stage_a.pdb")
CONFIG = PIPELINE / "inputs" / os.environ.get("CONFIG_JSON", "stage_a_hidden_minibinder_config.json")
OUT    = Path("/workspace/outputs/rfd3_stage_a")
N      = int(os.environ.get("N_DESIGNS", 200))
BATCH  = int(os.environ.get("BATCH_SIZE", 10))
NBATCH = max(1, N // BATCH)
MB_LEN = os.environ.get("MB_LENGTH_RANGE", "35-55")
OUT.mkdir(parents=True, exist_ok=True)

cfg = json.load(open(CONFIG))
rfd3 = cfg["rfd3"]
contig = rfd3["contig"]
mb_range = rfd3.get("mb_length_range", os.environ.get("MB_LENGTH_RANGE", "35-55"))

print(f"Stage A RFd3: {N} designs")
print(f"  input: {INPUT}")
print(f"  contig: {contig}")
print(f"  unindex: {rfd3.get('unindex')}")
print(f"  layout: VL+CH1 in output; VH/CL/T unindexed fixed context; MB {mb_range} aa")
print(f"  hotspots: {rfd3['select_hotspots'][:80]}...")

spec = DesignInputSpecification.safe_init(
    input=str(INPUT),
    contig=contig,
    unindex=rfd3.get("unindex"),
    select_hotspots=rfd3["select_hotspots"],
    select_fixed_atoms=rfd3["select_fixed_atoms"],
)

model = RFD3InferenceEngine(**RFD3InferenceConfig(diffusion_batch_size=BATCH))
saved = 0
for bi in range(NBATCH):
    print(f"  Batch {bi+1}/{NBATCH} ...", flush=True)
    for _key, designs in model.run(inputs=spec, out_dir=None, n_batches=1).items():
        for i, d in enumerate(designs):
            to_cif_file(d.atom_array, str(OUT / f"sa_b{bi:03d}_{i:03d}.cif"))
            saved += 1
print(f"Stage A RFd3 done: {saved} designs → {OUT}")
PYEOF

echo ""
echo "=== Step 1: RFdiffusion3 Stage A ==="
docker run --rm --gpus all \
    -v "$PIPELINE:/workspace" \
    -e N_DESIGNS=$N_DESIGNS \
    -e BATCH_SIZE=$BATCH_SIZE \
    -e MB_LENGTH_RANGE=$MB_LENGTH_RANGE \
    -e INPUT_PDB=$INPUT_PDB_NAME \
    -e CONFIG_JSON=$CONFIG_JSON \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_rfd3_stage_a.py

echo "RFd3 Stage A: $(ls $PIPELINE/outputs/rfd3_stage_a/sa_*.cif 2>/dev/null | wc -l) designs"
echo ""
echo "===== Stage A RFd3 complete: $(date) ====="
echo "Next: MPNN → Boltz → global RMSD validation (wire in follow-up)"
