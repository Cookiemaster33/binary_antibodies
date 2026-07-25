#!/bin/bash
# ============================================================
# run_stage_0_vhvL_interface.sh
# Stage 0: partial-diffusion redesign of VH–VL framework interface.
#
# RFd3 → ProteinMPNN (Fv framework) → Boltz apo (A+B+C+D) + holo (A+B+C+D+T)
# → stage_0_scoring.py
# ============================================================
set -eo pipefail

PIPELINE=${PIPELINE_DIR:-/home/ubuntu/pipeline}
LOG=$PIPELINE/stage_0_pipeline.log
exec > >(tee -a "$LOG") 2>&1

N_DESIGNS=${N_DESIGNS:-200}
BATCH_SIZE=${BATCH_SIZE:-10}
N_MPNN_SEQS=${N_MPNN_SEQS:-8}
TOP_N=${TOP_N:-50}
PARTIAL_T=${PARTIAL_T:-12.0}
N_BATCHES=$(( (N_DESIGNS + BATCH_SIZE - 1) / BATCH_SIZE ))

INPUT_PDB_NAME=${INPUT_PDB:-fab_stage_0_vhvL_interface.pdb}
CONFIG_JSON=${CONFIG_JSON:-stage_0_vhvL_interface_config.json}

echo "===== Stage 0 VH–VL interface pipeline: $(date) ====="
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "Designs: $N_DESIGNS | partial_t: ${PARTIAL_T} Å"

mkdir -p "$PIPELINE/inputs" "$PIPELINE/outputs/rfd3_stage_0" \
         "$PIPELINE/outputs/mpnn_stage_0" "$PIPELINE/boltz_inputs_apo" \
         "$PIPELINE/boltz_inputs_holo" "$PIPELINE/boltz_outputs_apo" \
         "$PIPELINE/boltz_outputs_holo" "$PIPELINE/final" "$PIPELINE/binary_antibodies"

sudo chmod 666 /var/run/docker.sock 2>/dev/null || true
docker pull rosettacommons/foundry:latest > "$PIPELINE/docker_pull.log" 2>&1 || true

cp -f "$PIPELINE/../binary_antibodies/stage_0_scoring.py" "$PIPELINE/binary_antibodies/" 2>/dev/null || \
  cp -f /workspace/binary_antibodies/stage_0_scoring.py "$PIPELINE/binary_antibodies/" 2>/dev/null || true
cp -f /workspace/binary_antibodies/fab_hidden_switch.py "$PIPELINE/binary_antibodies/" 2>/dev/null || true

# ── RFd3 partial diffusion ───────────────────────────────────────
cat > "$PIPELINE/run_rfd3_stage_0.py" << 'PYEOF'
import json, os, torch
from pathlib import Path
from atomworks.io.utils.io_utils import to_cif_file
from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine
from rfd3.inference.input_parsing import DesignInputSpecification

torch.set_float32_matmul_precision("high")

PIPELINE = Path("/workspace")
INPUT  = PIPELINE / "inputs" / os.environ.get("INPUT_PDB", "fab_stage_0_vhvL_interface.pdb")
CONFIG = PIPELINE / "inputs" / os.environ.get("CONFIG_JSON", "stage_0_vhvL_interface_config.json")
OUT    = Path("/workspace/outputs/rfd3_stage_0")
N      = int(os.environ.get("N_DESIGNS", 200))
BATCH  = int(os.environ.get("BATCH_SIZE", 10))
NBATCH = max(1, N // BATCH)
PARTIAL_T = float(os.environ.get("PARTIAL_T", "12.0"))
OUT.mkdir(parents=True, exist_ok=True)

cfg = json.load(open(CONFIG))
rfd3 = cfg["rfd3"]
vh_end = int(cfg["chains"]["A"].split("-")[-1].split()[-1])
vl_end = int(cfg["chains"]["B"].split("-")[-1].split()[-1])
contig = rfd3.get("contig", f"A1-{vh_end}/0,B1-{vl_end}")

print(f"Stage 0 RFd3: {N} designs | partial_t={PARTIAL_T}")
print(f"  input: {INPUT}")
print(f"  contig: {contig}")
print(f"  hotspots: {rfd3['select_hotspots'][:60]}...")

spec = DesignInputSpecification.safe_init(
    input=str(INPUT),
    contig=contig,
    partial_t=PARTIAL_T,
    select_hotspots=rfd3["select_hotspots"],
    select_fixed_atoms=rfd3["select_fixed_atoms"],
)

model = RFD3InferenceEngine(**RFD3InferenceConfig(diffusion_batch_size=BATCH))
saved = 0
for bi in range(NBATCH):
    print(f"  Batch {bi+1}/{NBATCH} ...", flush=True)
    for _key, designs in model.run(inputs=spec, out_dir=None, n_batches=1).items():
        for i, d in enumerate(designs):
            to_cif_file(d.atom_array, str(OUT / f"s0_b{bi:03d}_{i:03d}.cif"))
            saved += 1
print(f"Stage 0 RFd3 done: {saved} designs → {OUT}")
PYEOF

echo ""
echo "=== Step 1: RFdiffusion3 Stage 0 (partial diffusion) ==="
docker run --rm --gpus all \
    -v "$PIPELINE:/workspace" \
    -e N_DESIGNS=$N_DESIGNS \
    -e BATCH_SIZE=$BATCH_SIZE \
    -e PARTIAL_T=$PARTIAL_T \
    -e INPUT_PDB=$INPUT_PDB_NAME \
    -e CONFIG_JSON=$CONFIG_JSON \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_rfd3_stage_0.py

echo "RFd3 Stage 0: $(ls $PIPELINE/outputs/rfd3_stage_0/s0_*.cif 2>/dev/null | wc -l) designs"

# ── ProteinMPNN on Fv framework ──────────────────────────────────
cat > "$PIPELINE/run_mpnn_stage_0.py" << 'PYEOF'
import json, os
from pathlib import Path
from atomworks.io.utils.io_utils import load_any
from mpnn.inference_engines.mpnn import MPNNInferenceEngine

RFD3 = Path("/workspace/outputs/rfd3_stage_0")
OUT  = Path("/workspace/outputs/mpnn_stage_0")
NSEQS = int(os.environ.get("N_MPNN_SEQS", 8))
OUT.mkdir(parents=True, exist_ok=True)

VH_CDR = [(26,35),(50,65),(95,102)]
VL_CDR = [(24,34),(50,56),(89,97)]

def in_ranges(r, ranges):
    return any(lo <= r <= hi for lo, hi in ranges)

def designed_framework(chain, res_ids, cdr_ranges):
    return [f"{chain}{r}" for r in sorted(res_ids) if not in_ranges(r, cdr_ranges)]

engine = MPNNInferenceEngine(
    model_type="protein_mpnn", is_legacy_weights=True,
    out_directory=None, write_structures=False, write_fasta=False,
)

results = []
for cif in sorted(RFD3.glob("s0_*.cif")):
    raw = load_any(str(cif)); aa = raw[0] if hasattr(raw, "__getitem__") else raw
    chains = sorted(set(aa.chain_id))
    designed = []
    if "A" in chains:
        res_a = sorted(set(aa[aa.chain_id == "A"].res_id))
        designed += designed_framework("A", res_a, VH_CDR)
    if "B" in chains:
        res_b = sorted(set(aa[aa.chain_id == "B"].res_id))
        designed += designed_framework("B", res_b, VL_CDR)
    # C, D, T stay fixed (not in designed_residues)

    aa_all = aa
    result = engine.run(atom_arrays=[aa_all], input_dicts=[{
        "batch_size": NSEQS, "remove_waters": True, "designed_residues": designed,
    }])
    for r in (result if isinstance(result, list) else [result]):
        seqs = r.output_dict.get("sequences", [])
        if not seqs:
            seq = r.output_dict.get("designed_sequence", "")
            seqs = [seq] if seq else []
        for si, seq in enumerate(seqs):
            chain_seqs = {}
            pos = 0
            for ch in chains:
                n = len(sorted(set(aa_all[aa_all.chain_id == ch].res_id)))
                chain_seqs[ch] = seq[pos:pos+n] if len(seq) >= pos+n else ""
                pos += n
            results.append({
                "backbone": cif.stem, "seq_idx": si,
                "sequence_recovery": float(r.output_dict.get("sequence_recovery", 0)),
                "chains": chain_seqs,
            })

json.dump(results, open(OUT / "all_sequences.json", "w"), indent=2)
print(f"MPNN Stage 0: {len(results)} sequences")
PYEOF

echo ""
echo "=== Step 2: ProteinMPNN (Fv framework) ==="
docker run --rm --gpus all \
    -v "$PIPELINE:/workspace" \
    -e N_MPNN_SEQS=$N_MPNN_SEQS \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_mpnn_stage_0.py

# ── Boltz inputs (apo + holo) ────────────────────────────────────
cat > "$PIPELINE/build_boltz_stage_0.py" << 'PYEOF'
import json, os
from pathlib import Path

PIPELINE = Path("/workspace")
MPNN = PIPELINE / "outputs" / "mpnn_stage_0" / "all_sequences.json"
APO  = PIPELINE / "boltz_inputs_apo"
HOLO = PIPELINE / "boltz_inputs_holo"
TOP_N = int(os.environ.get("TOP_N", 50))
APO.mkdir(parents=True, exist_ok=True)
HOLO.mkdir(parents=True, exist_ok=True)

data = json.load(open(MPNN))
data.sort(key=lambda x: x["sequence_recovery"])
seen, picks = set(), []
for r in data:
    if r["backbone"] not in seen:
        seen.add(r["backbone"]); picks.append(r)
    if len(picks) >= TOP_N: break

design_list = []
for i, r in enumerate(picks):
    chains = r["chains"]
    rank = i + 1
    name = f"rank{rank:02d}_{r['backbone']}_s{r['seq_idx']}"
    apo_yaml = "sequences:\n"
    holo_yaml = "sequences:\n"
    for ch in ("A", "B", "C", "D"):
        if ch in chains and chains[ch]:
            block = f"  - protein:\n      id: {ch}\n      sequence: \"{chains[ch]}\"\n      msa: empty\n"
            apo_yaml += block
            holo_yaml += block
    if "T" in chains and chains["T"]:
        holo_yaml += (
            f"  - protein:\n      id: T\n      sequence: \"{chains['T']}\"\n      msa: empty\n"
        )
    (APO / f"{name}.yaml").write_text(apo_yaml)
    (HOLO / f"{name}.yaml").write_text(holo_yaml)
    design_list.append({"rank": rank, "name": name, "backbone": r["backbone"]})

json.dump({"designs": design_list}, open(PIPELINE / "final" / "top_designs_stage_0.json", "w"), indent=2)
print(f"Boltz inputs: {len(design_list)} apo + holo YAML pairs")
PYEOF

echo ""
echo "=== Step 3: Build Boltz inputs ==="
docker run --rm --gpus all \
    -v "$PIPELINE:/workspace" \
    -e TOP_N=$TOP_N \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/build_boltz_stage_0.py

pip install boltz[cuda] -U -q > "$PIPELINE/boltz_install.log" 2>&1 || true
pip install -q 'networkx>=3.0' 'platformdirs>=3.0' 2>/dev/null || true

N_APO=$(ls $PIPELINE/boltz_inputs_apo/*.yaml 2>/dev/null | wc -l)
echo ""
echo "=== Step 4a: Boltz apo (A+B+C+D) ==="
$HOME/.local/bin/boltz predict $PIPELINE/boltz_inputs_apo \
    --out_dir $PIPELINE/boltz_outputs_apo \
    --devices 1 --num_workers 2 --override 2>&1 | \
    grep -E "Predicting|Saving|Done|Error" | head -30 || true

echo ""
echo "=== Step 4b: Boltz holo (A+B+C+D+T) ==="
$HOME/.local/bin/boltz predict $PIPELINE/boltz_inputs_holo \
    --out_dir $PIPELINE/boltz_outputs_holo \
    --devices 1 --num_workers 2 --override 2>&1 | \
    grep -E "Predicting|Saving|Done|Error" | head -30 || true

echo ""
echo "=== Step 5: Stage 0 scoring ==="
docker run --rm --gpus all \
    -v "$PIPELINE:/workspace" \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/binary_antibodies/stage_0_scoring.py \
        --pipeline-dir /workspace \
        --results-file final/stage_0_results.json

echo ""
echo "===== Stage 0 complete: $(date) ====="
