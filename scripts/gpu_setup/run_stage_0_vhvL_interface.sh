#!/bin/bash
# ============================================================
# run_stage_0_vhvL_interface.sh
# Stage 0: split-chain ProteinMPNN on VH–VL interface → holo Boltz validation
#
# Split MPNN (N sequences) → rank by static VH–VL weakening → holo Boltz → scoring
# ============================================================
set -eo pipefail

PIPELINE=${PIPELINE_DIR:-/home/ubuntu/pipeline}
LOG=$PIPELINE/stage_0_pipeline.log
exec > >(tee -a "$LOG") 2>&1

N_MPNN_SEQS=${N_MPNN_SEQS:-1000}
MPNN_BATCH=${MPNN_BATCH:-100}
TOP_N=${TOP_N:-100}
INPUT_PDB_NAME=${INPUT_PDB:-fab_stage_0_vhvL_interface.pdb}
SPLIT_PDB_NAME=${SPLIT_PDB:-fab_stage_0_split_mpnn.pdb}
INTERFACE_SCOPE=${INTERFACE_SCOPE:-fv}

echo "===== Stage 0 split-MPNN pipeline: $(date) ====="
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "MPNN: $N_MPNN_SEQS (batch $MPNN_BATCH) | Holo Boltz top: $TOP_N (static rank)"

mkdir -p "$PIPELINE/inputs" "$PIPELINE/outputs/mpnn_stage_0" \
         "$PIPELINE/boltz_inputs_holo" "$PIPELINE/boltz_outputs_holo" \
         "$PIPELINE/final" "$PIPELINE/binary_antibodies" "$PIPELINE/scripts" \
         "$PIPELINE/structures/interface/pisa_wt_fv"

sudo chmod 666 /var/run/docker.sock 2>/dev/null || true
docker pull rosettacommons/foundry:latest > "$PIPELINE/docker_pull.log" 2>&1 || true
docker pull pdbegroup/pisa:latest >> "$PIPELINE/docker_pull.log" 2>&1 || true

cp -f "$PIPELINE/../binary_antibodies/stage_0_scoring.py" "$PIPELINE/binary_antibodies/" 2>/dev/null || \
  cp -f /workspace/binary_antibodies/stage_0_scoring.py "$PIPELINE/binary_antibodies/" 2>/dev/null || true
cp -f /workspace/binary_antibodies/pisa_scoring.py "$PIPELINE/binary_antibodies/" 2>/dev/null || true
cp -f /workspace/binary_antibodies/fab_hidden_switch.py "$PIPELINE/binary_antibodies/" 2>/dev/null || true
cp -f "$PIPELINE/../scripts/build_boltz_stage_0_inputs.py" "$PIPELINE/scripts/" 2>/dev/null || \
  cp -f /workspace/scripts/build_boltz_stage_0_inputs.py "$PIPELINE/scripts/" 2>/dev/null || true

echo ""
echo "=== Step 0: Build Stage 0 targets (PISA WT interface → MPNN residue set) ==="
pip install -q "numpy<2" "networkx>=3.0" biopython biotite 2>/dev/null || true
export PYTHONPATH="$PIPELINE:$PYTHONPATH"
python3 "$PIPELINE/scripts/build_stage_0_design_target.py" \
    --source-pdb "$PIPELINE/structures/1N8Z.pdb" \
    --out-pdb "$PIPELINE/inputs/$INPUT_PDB_NAME" \
    --out-split-pdb "$PIPELINE/inputs/$SPLIT_PDB_NAME" \
    --out-json "$PIPELINE/inputs/$CONFIG_JSON" \
    --interface-scope "$INTERFACE_SCOPE" \
    --n-mpnn-seqs "$N_MPNN_SEQS" \
    --top-n-boltz "$TOP_N"

python3 - <<'PYEOF'
import json, os, sys
from pathlib import Path

pipeline = Path(os.environ.get("PIPELINE_DIR", "/home/ubuntu/pipeline"))
config_path = pipeline / "inputs" / os.environ.get("CONFIG_JSON", "stage_0_vhvL_interface_config.json")
cfg = json.loads(config_path.read_text())
idef = cfg.get("interface_definition", {})
source = idef.get("source")
scope = idef.get("scope", "fv")
n_designed = len(cfg.get("split_mpnn", {}).get("designed_residues", []))
print(f"  interface_definition.source = {source!r}")
print(f"  interface_definition.scope = {scope!r}")
print(f"  designed_residues = {n_designed}")
if source != "pisa":
    print("ERROR: PISA interface definition missing — aborting test run.")
    sys.exit(1)
if scope == "full_fab":
    ch1 = idef.get("ch1_framework_interface", [])
    cl = idef.get("cl_framework_interface", [])
    print(f"  CH1/CL interface residues: C={len(ch1)} D={len(cl)}")
    if not ch1 or not cl:
        print("ERROR: full_fab scope missing CH1/CL interface residues.")
        sys.exit(1)
PYEOF

# ── Split-chain ProteinMPNN ─────────────────────────────────────
cat > "$PIPELINE/run_mpnn_stage_0.py" << 'PYEOF'
import json, os
from pathlib import Path
from atomworks.io.utils.io_utils import load_any
from mpnn.inference_engines.mpnn import MPNNInferenceEngine

PIPELINE = Path("/workspace")
SPLIT_PDB = PIPELINE / "inputs" / os.environ.get("SPLIT_PDB", "fab_stage_0_split_mpnn.pdb")
CONFIG = PIPELINE / "inputs" / os.environ.get("CONFIG_JSON", "stage_0_vhvL_interface_config.json")
OUT  = Path("/workspace/outputs/mpnn_stage_0")
NSEQS = int(os.environ.get("N_MPNN_SEQS", 1000))
BATCH = int(os.environ.get("MPNN_BATCH", 100))
OUT.mkdir(parents=True, exist_ok=True)

cfg = json.load(open(CONFIG))
designed = cfg.get("split_mpnn", {}).get("designed_residues", [])
if not designed:
    raise SystemExit("Config missing split_mpnn.designed_residues")

print(f"Stage 0 split MPNN: {SPLIT_PDB}")
print(f"  approach: {cfg.get('approach')}")
print(f"  designed residues ({len(designed)}): {', '.join(designed[:10])}{'...' if len(designed) > 10 else ''}")
print(f"  target sequences: {NSEQS} (batch size {BATCH})")

raw = load_any(str(SPLIT_PDB))
aa = raw[0] if hasattr(raw, "__getitem__") else raw
chains = sorted(set(aa.chain_id))

engine = MPNNInferenceEngine(
    model_type="protein_mpnn", is_legacy_weights=True,
    out_directory=None, write_structures=False, write_fasta=False,
)

def per_seq_scores(out_dict, n_seqs):
    for key in ("scores", "global_scores", "score", "mpnn_scores"):
        val = out_dict.get(key)
        if isinstance(val, (list, tuple)) and len(val) >= n_seqs:
            return [float(v) for v in val[:n_seqs]]
    return None

results = []
seq_counter = 0
remaining = NSEQS

while remaining > 0:
    bs = min(BATCH, remaining)
    print(f"  MPNN batch: {bs} sequences ({seq_counter} done)", flush=True)
    result = engine.run(atom_arrays=[aa], input_dicts=[{
        "batch_size": bs, "remove_waters": True, "designed_residues": designed,
    }])
    batch_results = result if isinstance(result, list) else [result]
    for r in batch_results:
        od = r.output_dict
        seqs = od.get("sequences", [])
        if not seqs:
            seq = od.get("designed_sequence", "")
            seqs = [seq] if seq else []
        scores = per_seq_scores(od, len(seqs))
        batch_recovery = float(od.get("sequence_recovery", 0))
        for si, seq in enumerate(seqs):
            chain_seqs = {}
            pos = 0
            for ch in chains:
                n = len(sorted(set(aa[aa.chain_id == ch].res_id)))
                chain_seqs[ch] = seq[pos:pos + n] if len(seq) >= pos + n else ""
                pos += n
            mpnn_score = scores[si] if scores else batch_recovery
            results.append({
                "backbone": "native_split",
                "seq_idx": seq_counter,
                "mpnn_score": float(mpnn_score),
                "sequence_recovery": float(batch_recovery),
                "chains": chain_seqs,
            })
            seq_counter += 1
    remaining -= bs

json.dump(results, open(OUT / "all_sequences.json", "w"), indent=2)
print(f"MPNN Stage 0: {len(results)} sequences → {OUT}")
if results:
  scores = [r["mpnn_score"] for r in results]
  print(f"  MPNN score range: {min(scores):.4f} .. {max(scores):.4f}")
PYEOF

echo ""
echo "=== Step 1: Split-chain ProteinMPNN ==="
docker run --rm --gpus all \
    -v "$PIPELINE:/workspace" \
    -e N_MPNN_SEQS=$N_MPNN_SEQS \
    -e MPNN_BATCH=$MPNN_BATCH \
    -e SPLIT_PDB=$SPLIT_PDB_NAME \
    -e CONFIG_JSON=$CONFIG_JSON \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_mpnn_stage_0.py

echo ""
echo "=== Step 2: Build holo Boltz inputs (rank by static VH–VL weakening) ==="
docker run --rm --gpus all \
    -v "$PIPELINE:/workspace" \
    -e TOP_N=$TOP_N \
    -e CONFIG_JSON=$CONFIG_JSON \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    -e PYTHONPATH=/workspace \
    rosettacommons/foundry:latest \
    python3 /workspace/scripts/build_boltz_stage_0_inputs.py \
        --pipeline-dir /workspace \
        --reference-pdb /workspace/inputs/$INPUT_PDB_NAME \
        --config-json /workspace/inputs/$CONFIG_JSON \
        --top-n $TOP_N

pip install boltz[cuda] -U -q > "$PIPELINE/boltz_install.log" 2>&1 || true
pip install -q 'networkx>=3.0' 'platformdirs>=3.0' 2>/dev/null || true

# cuequivariance CUDA kernels often mismatch pip-installed boltz (kv_lengths crash).
# Use PyTorch attention for all GPUs — slower but reliable.
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo "unknown")
BOLTZ_EXTRA_ARGS=(--no_kernels)
echo "GPU=$GPU_NAME — boltz --no_kernels (avoids cuequivariance version skew)"

N_HOLO=$(ls $PIPELINE/boltz_inputs_holo/*.yaml 2>/dev/null | wc -l)
echo ""
echo "=== Step 3: Boltz holo (A+B+C+D+T) — $N_HOLO designs ==="
set +e
$HOME/.local/bin/boltz predict $PIPELINE/boltz_inputs_holo \
    --out_dir $PIPELINE/boltz_outputs_holo \
    --devices 1 --num_workers 2 --override \
    "${BOLTZ_EXTRA_ARGS[@]}" 2>&1 | tee "$PIPELINE/boltz_holo_run.log" | \
    grep -E "Predicting|Saving|Done|Error|failed|Traceback" | tail -30
BOLTZ_RC=${PIPESTATUS[0]}
set -e
N_CIF=$(find "$PIPELINE/boltz_outputs_holo" -name "*_model_0.cif" 2>/dev/null | wc -l)
echo "Boltz holo exit=$BOLTZ_RC structures=$N_CIF / $N_HOLO"
if [[ $BOLTZ_RC -ne 0 || "$N_CIF" -eq 0 ]]; then
  echo "ERROR: Boltz holo failed or produced no structures. See $PIPELINE/boltz_holo_run.log"
  tail -40 "$PIPELINE/boltz_holo_run.log" || true
  exit 1
fi

echo ""
echo "=== Step 4: Stage 0 scoring (static + holo + PISA VH–VL) ==="
docker pull pdbegroup/pisa:latest >> "$PIPELINE/docker_pull.log" 2>&1 || true
pip install -q biotite numpy scipy 2>/dev/null || true
python3 "$PIPELINE/binary_antibodies/stage_0_scoring.py" \
    --pipeline-dir "$PIPELINE" \
    --reference-pdb "$PIPELINE/inputs/$INPUT_PDB_NAME" \
    --config-json "$PIPELINE/inputs/$CONFIG_JSON" \
    --results-file final/stage_0_results.json \
    --top-n "$TOP_N"

echo ""
echo "===== Stage 0 complete: $(date) ====="
