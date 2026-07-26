#!/bin/bash
# ============================================================
# run_stage_0_vhvL_interface.sh
# Stage 0: split-chain ProteinMPNN partial de-greasing of VH–VL
# framework interface (no RFd3).
#
# Split MPNN (rim degrease) → Boltz apo (A+B+C+D) + holo (+T) → scoring
# ============================================================
set -eo pipefail

PIPELINE=${PIPELINE_DIR:-/home/ubuntu/pipeline}
LOG=$PIPELINE/stage_0_pipeline.log
exec > >(tee -a "$LOG") 2>&1

N_MPNN_SEQS=${N_MPNN_SEQS:-32}
TOP_N=${TOP_N:-50}
INPUT_PDB_NAME=${INPUT_PDB:-fab_stage_0_vhvL_interface.pdb}
SPLIT_PDB_NAME=${SPLIT_PDB:-fab_stage_0_split_mpnn.pdb}
CONFIG_JSON=${CONFIG_JSON:-stage_0_vhvL_interface_config.json}

echo "===== Stage 0 split-MPNN pipeline: $(date) ====="
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "MPNN sequences: $N_MPNN_SEQS | top Boltz: $TOP_N"

mkdir -p "$PIPELINE/inputs" "$PIPELINE/outputs/mpnn_stage_0" \
         "$PIPELINE/boltz_inputs_apo" "$PIPELINE/boltz_inputs_holo" \
         "$PIPELINE/boltz_outputs_apo" "$PIPELINE/boltz_outputs_holo" \
         "$PIPELINE/final" "$PIPELINE/binary_antibodies"

sudo chmod 666 /var/run/docker.sock 2>/dev/null || true
docker pull rosettacommons/foundry:latest > "$PIPELINE/docker_pull.log" 2>&1 || true

cp -f "$PIPELINE/../binary_antibodies/stage_0_scoring.py" "$PIPELINE/binary_antibodies/" 2>/dev/null || \
  cp -f /workspace/binary_antibodies/stage_0_scoring.py "$PIPELINE/binary_antibodies/" 2>/dev/null || true
cp -f /workspace/binary_antibodies/fab_hidden_switch.py "$PIPELINE/binary_antibodies/" 2>/dev/null || true

# ── Split-chain ProteinMPNN (partial de-grease) ─────────────────
cat > "$PIPELINE/run_mpnn_stage_0.py" << 'PYEOF'
import json, os
from pathlib import Path
from atomworks.io.utils.io_utils import load_any
from mpnn.inference_engines.mpnn import MPNNInferenceEngine

PIPELINE = Path("/workspace")
SPLIT_PDB = PIPELINE / "inputs" / os.environ.get("SPLIT_PDB", "fab_stage_0_split_mpnn.pdb")
CONFIG = PIPELINE / "inputs" / os.environ.get("CONFIG_JSON", "stage_0_vhvL_interface_config.json")
OUT  = Path("/workspace/outputs/mpnn_stage_0")
NSEQS = int(os.environ.get("N_MPNN_SEQS", 32))
OUT.mkdir(parents=True, exist_ok=True)

cfg = json.load(open(CONFIG))
designed = cfg.get("split_mpnn", {}).get("designed_residues", [])
if not designed:
    raise SystemExit("Config missing split_mpnn.designed_residues")

print(f"Stage 0 split MPNN: {SPLIT_PDB}")
print(f"  designed residues ({len(designed)}): {', '.join(designed[:12])}{'...' if len(designed) > 12 else ''}")

raw = load_any(str(SPLIT_PDB))
aa = raw[0] if hasattr(raw, "__getitem__") else raw
chains = sorted(set(aa.chain_id))

engine = MPNNInferenceEngine(
    model_type="protein_mpnn", is_legacy_weights=True,
    out_directory=None, write_structures=False, write_fasta=False,
)

result = engine.run(atom_arrays=[aa], input_dicts=[{
    "batch_size": NSEQS, "remove_waters": True, "designed_residues": designed,
}])

results = []
for r in (result if isinstance(result, list) else [result]):
    seqs = r.output_dict.get("sequences", [])
    if not seqs:
        seq = r.output_dict.get("designed_sequence", "")
        seqs = [seq] if seq else []
    for si, seq in enumerate(seqs):
        chain_seqs = {}
        pos = 0
        for ch in chains:
            n = len(sorted(set(aa[aa.chain_id == ch].res_id)))
            chain_seqs[ch] = seq[pos:pos + n] if len(seq) >= pos + n else ""
            pos += n
        results.append({
            "backbone": "native_split",
            "seq_idx": si,
            "sequence_recovery": float(r.output_dict.get("sequence_recovery", 0)),
            "chains": chain_seqs,
        })

json.dump(results, open(OUT / "all_sequences.json", "w"), indent=2)
print(f"MPNN Stage 0: {len(results)} sequences → {OUT}")
PYEOF

echo ""
echo "=== Step 1: Split-chain ProteinMPNN (partial de-grease) ==="
docker run --rm --gpus all \
    -v "$PIPELINE:/workspace" \
    -e N_MPNN_SEQS=$N_MPNN_SEQS \
    -e SPLIT_PDB=$SPLIT_PDB_NAME \
    -e CONFIG_JSON=$CONFIG_JSON \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_mpnn_stage_0.py

# ── Boltz inputs (apo + holo) ────────────────────────────────────
cat > "$PIPELINE/build_boltz_stage_0.py" << 'PYEOF'
import json, os
from pathlib import Path

PIPELINE = Path("/workspace")
MPNN = PIPELINE / "outputs" / "mpnn_stage_0" / "all_sequences.json"
CONFIG = PIPELINE / "inputs" / os.environ.get("CONFIG_JSON", "stage_0_vhvL_interface_config.json")
APO  = PIPELINE / "boltz_inputs_apo"
HOLO = PIPELINE / "boltz_inputs_holo"
TOP_N = int(os.environ.get("TOP_N", 50))
APO.mkdir(parents=True, exist_ok=True)
HOLO.mkdir(parents=True, exist_ok=True)

cfg = json.load(open(CONFIG))
native = cfg.get("native_chain_sequences", {})
epitope = native.get("T", "")

data = json.load(open(MPNN))
data.sort(key=lambda x: x["sequence_recovery"])
picks = data[:TOP_N]

design_list = []
for i, r in enumerate(picks):
    chains = dict(r["chains"])
    for ch, seq in native.items():
        if ch not in chains and seq:
            chains[ch] = seq
    rank = i + 1
    name = f"rank{rank:02d}_s0_{r['backbone']}_s{r['seq_idx']}"
    apo_yaml = "sequences:\n"
    holo_yaml = "sequences:\n"
    for ch in ("A", "B", "C", "D"):
        if chains.get(ch):
            block = f"  - protein:\n      id: {ch}\n      sequence: \"{chains[ch]}\"\n      msa: empty\n"
            apo_yaml += block
            holo_yaml += block
    if epitope:
        holo_yaml += (
            f"  - protein:\n      id: T\n      sequence: \"{epitope}\"\n      msa: empty\n"
        )
    (APO / f"{name}.yaml").write_text(apo_yaml)
    (HOLO / f"{name}.yaml").write_text(holo_yaml)
    design_list.append({"rank": rank, "name": name, "backbone": r["backbone"], "seq_idx": r["seq_idx"]})

json.dump({"designs": design_list}, open(PIPELINE / "final" / "top_designs_stage_0.json", "w"), indent=2)
print(f"Boltz inputs: {len(design_list)} apo + holo YAML pairs")
PYEOF

echo ""
echo "=== Step 2: Build Boltz inputs ==="
docker run --rm --gpus all \
    -v "$PIPELINE:/workspace" \
    -e TOP_N=$TOP_N \
    -e CONFIG_JSON=$CONFIG_JSON \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/build_boltz_stage_0.py

pip install boltz[cuda] -U -q > "$PIPELINE/boltz_install.log" 2>&1 || true
pip install -q 'networkx>=3.0' 'platformdirs>=3.0' 2>/dev/null || true

N_APO=$(ls $PIPELINE/boltz_inputs_apo/*.yaml 2>/dev/null | wc -l)
echo ""
echo "=== Step 3a: Boltz apo (A+B+C+D) ==="
$HOME/.local/bin/boltz predict $PIPELINE/boltz_inputs_apo \
    --out_dir $PIPELINE/boltz_outputs_apo \
    --devices 1 --num_workers 2 --override 2>&1 | \
    grep -E "Predicting|Saving|Done|Error" | head -40 || true

echo ""
echo "=== Step 3b: Boltz holo (A+B+C+D+T) ==="
$HOME/.local/bin/boltz predict $PIPELINE/boltz_inputs_holo \
    --out_dir $PIPELINE/boltz_outputs_holo \
    --devices 1 --num_workers 2 --override 2>&1 | \
    grep -E "Predicting|Saving|Done|Error" | head -40 || true

echo ""
echo "=== Step 4: Stage 0 scoring ==="
docker run --rm --gpus all \
    -v "$PIPELINE:/workspace" \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/binary_antibodies/stage_0_scoring.py \
        --pipeline-dir /workspace \
        --reference-pdb /workspace/inputs/fab_stage_0_vhvL_interface.pdb \
        --config-json /workspace/inputs/stage_0_vhvL_interface_config.json \
        --results-file final/stage_0_results.json

echo ""
echo "===== Stage 0 complete: $(date) ====="
