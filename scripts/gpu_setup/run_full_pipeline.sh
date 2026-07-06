#!/bin/bash
# ============================================================
# run_full_pipeline.sh  v4
# RFdiffusion3 → ProteinMPNN → Boltz-2 → global assembly RMSD
# Optional round 2: score round 1, then partial-diffusion on
#   lowest global-RMSD round-1 RFd3 backbones (not Boltz structures)
#
# RFD3_ROUNDS=1 (default): single de novo pass + validation
# RFD3_ROUNDS=2: round 1 → validate → top global RMSD RFd3 CIFs → round 2 → validate
# ============================================================
set -eo pipefail
PIPELINE=/home/ubuntu/pipeline
LOG=$PIPELINE/full_pipeline.log
exec > >(tee -a "$LOG") 2>&1

# ── User-configurable parameters ─────────────────────────────────
N_DESIGNS=${N_DESIGNS:-200}
BATCH_SIZE=${BATCH_SIZE:-10}
N_MPNN_SEQS=${N_MPNN_SEQS:-8}
N_BATCHES=$(( (N_DESIGNS + BATCH_SIZE - 1) / BATCH_SIZE ))
TOP_N=${TOP_N:-50}
MB_LENGTH_RANGE=${MB_LENGTH_RANGE:-35-70}

RFD3_ROUNDS=${RFD3_ROUNDS:-1}
RFD3_ROUND2_TEMPLATES=${RFD3_ROUND2_TEMPLATES:-5}
RFD3_ROUND2_DESIGNS_PER_TEMPLATE=${RFD3_ROUND2_DESIGNS_PER_TEMPLATE:-8}
RFD3_ROUND2_BATCH_SIZE=${RFD3_ROUND2_BATCH_SIZE:-8}
RFD3_PARTIAL_T=${RFD3_PARTIAL_T:-2.0}
RFD3_ROUND2_VL_CONTEXT=${RFD3_ROUND2_VL_CONTEXT:-1}
RFD3_ROUND2_TEMPLATE_BACKBONES=${RFD3_ROUND2_TEMPLATE_BACKBONES:-}

PIPELINE_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "===== Full integrated pipeline v4: $(date) ====="
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "RFd3 rounds: $RFD3_ROUNDS | MB length range: $MB_LENGTH_RANGE"
if [ "$RFD3_ROUNDS" = "2" ]; then
    echo "Two-round mode: validate round 1 → top $RFD3_ROUND2_TEMPLATES by global RMSD → partial diffusion"
    echo "  round 2: $RFD3_ROUND2_DESIGNS_PER_TEMPLATE designs/template | partial_t=${RFD3_PARTIAL_T} Å"
fi
mkdir -p $PIPELINE/outputs/rfd3 $PIPELINE/outputs/rfd3_round2 \
         $PIPELINE/outputs/rfd3_round2_inputs \
         $PIPELINE/outputs/mpnn $PIPELINE/outputs/mpnn_round1 \
         $PIPELINE/boltz_inputs $PIPELINE/boltz_inputs_round1 \
         $PIPELINE/boltz_outputs $PIPELINE/boltz_outputs_round1 \
         $PIPELINE/final $PIPELINE/top_structures

# Copy helper scripts when launched from a different directory (skip if already in place).
mkdir -p "$PIPELINE/binary_antibodies"
for script in rfd3_round2.py push_results.sh; do
    src="$PIPELINE_SCRIPT_DIR/$script"
    dst="$PIPELINE/$script"
    if [ "$src" != "$dst" ] && [ -f "$src" ]; then
        cp "$src" "$dst"
    fi
done
scoring_src="$PIPELINE_SCRIPT_DIR/../../binary_antibodies/scoring.py"
if [ -f "$scoring_src" ]; then
    cp "$scoring_src" "$PIPELINE/binary_antibodies/scoring.py"
fi
chmod +x "$PIPELINE/push_results.sh" 2>/dev/null || true

# ── Step 0: Install Boltz-2 + pull Docker in parallel ────────────
pip install boltz[cuda] -U -q > $PIPELINE/boltz_install.log 2>&1 &
BOLTZ_PID=$!
sudo chmod 666 /var/run/docker.sock
docker pull rosettacommons/foundry:latest > $PIPELINE/docker_pull.log 2>&1 &
DOCKER_PID=$!
echo "Installing Boltz-2 (PID $BOLTZ_PID) + pulling Docker (PID $DOCKER_PID) in parallel..."
wait $DOCKER_PID && echo "Docker ready."

# ── Embedded Python helpers (written once, reused per phase) ─────
cat > $PIPELINE/run_rfd3.py << 'PYEOF'
import os, torch
from pathlib import Path
from atomworks.io.utils.io_utils import to_cif_file
from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine
from rfd3.inference.input_parsing import DesignInputSpecification
torch.set_float32_matmul_precision('high')

INPUT  = "/workspace/inputs/vh1_vl_nanobody_design_target.pdb"
OUT    = Path("/workspace/outputs/rfd3")
N      = int(os.environ.get("N_DESIGNS", 200))
BATCH  = int(os.environ.get("BATCH_SIZE", 10))
NBATCH = max(1, N // BATCH)
OUT.mkdir(parents=True, exist_ok=True)

HOTSPOTS = (
    "A11,A12,A13,A14,A15,A40,"
    "A79,A80,A81,A82,A83,A84,"
    "A103,A104,A105,A106,A107,A108,A109,A110,A111,A112,A113,A114,A115,"
    "C27,C28,C29,C30,C31,C32,C33,"
    "C52,C53,C54,C55,C56,C57,"
    "C99,C100,C101,C102,C103,C104,C105,C106,C107,C108,C109,C110,C111,C112"
)
FIXED_ATOMS = {
    "A1-115": "ALL",
    "B1-115": "ALL",
    "C1-115": "ALL",
}
MB_LENGTH_RANGE = os.environ.get("MB_LENGTH_RANGE", "35-70")
CONTIG = f"A1-115,{MB_LENGTH_RANGE},C1-115"
print(f"RFd3 round 1: {N} designs | minibinder length {MB_LENGTH_RANGE}")
spec = DesignInputSpecification.safe_init(
    input=INPUT, contig=CONTIG,
    select_hotspots=HOTSPOTS, select_fixed_atoms=FIXED_ATOMS,
)
model = RFD3InferenceEngine(**RFD3InferenceConfig(diffusion_batch_size=BATCH))
saved = 0
for bi in range(NBATCH):
    print(f"  Batch {bi+1}/{NBATCH} ...", flush=True)
    for key, designs in model.run(inputs=spec, out_dir=None, n_batches=1).items():
        for i, d in enumerate(designs):
            to_cif_file(d.atom_array, str(OUT / f"mb_b{bi:03d}_{i:03d}.cif"))
            saved += 1
print(f"RFd3 round 1 done: {saved} designs")
PYEOF

cat > $PIPELINE/run_mpnn.py << 'PYEOF'
import os, json
from pathlib import Path
from atomworks.io.utils.io_utils import load_any
from mpnn.inference_engines.mpnn import MPNNInferenceEngine

RFD3_SUBDIR = os.environ.get("RFD3_ACTIVE_SUBDIR", "rfd3")
MPNN_SUBDIR = os.environ.get("MPNN_OUTPUT_SUBDIR", "mpnn")
RFD3  = Path("/workspace/outputs") / RFD3_SUBDIR
OUT   = Path("/workspace/outputs") / MPNN_SUBDIR
NSEQS = int(os.environ.get("N_MPNN_SEQS", 8))
OUT.mkdir(parents=True, exist_ok=True)

engine = MPNNInferenceEngine(model_type="protein_mpnn", is_legacy_weights=True,
    out_directory=None, write_structures=False, write_fasta=False)

glob_pat = "r2_*.cif" if RFD3_SUBDIR == "rfd3_round2" else "mb_*.cif"
cifs = sorted(RFD3.glob(glob_pat))
print(f"MPNN ({RFD3_SUBDIR} → {MPNN_SUBDIR}): {len(cifs)} × {NSEQS}")

results = []
for idx, cif in enumerate(cifs):
    raw = load_any(str(cif)); aa = raw[0] if hasattr(raw,"__getitem__") else raw
    chains = sorted(set(aa.chain_id))
    # Connected VH1-MB-Nb is always chain A; chain B is optional VL steric context (round 2).
    ch = "A" if "A" in chains else chains[0]
    aa_chain = aa[aa.chain_id == ch]
    total = len(sorted(set(aa_chain.res_id)))
    mb_len = total - 230
    mb_start, mb_end = 116, 115 + mb_len
    chain_res = set(aa_chain.res_id)
    designed = [f"{ch}{r}" for r in range(mb_start, mb_end+1) if r in chain_res]
    result = engine.run(atom_arrays=[aa_chain], input_dicts=[{
        "batch_size": NSEQS, "remove_waters": True, "designed_residues": designed}])
    for r in (result if isinstance(result,list) else [result]):
        seq = r.output_dict.get("designed_sequence","")
        rec = float(r.output_dict.get("sequence_recovery",0))
        mb_seq = seq[mb_start-1:mb_end] if len(seq)>=mb_end else seq[mb_start-1:]
        results.append({"backbone":cif.stem,"full_sequence":seq,
                        "minibinder_sequence":mb_seq,"sequence_recovery":rec,
                        "mb_start":mb_start,"mb_end":mb_end,"mb_len":mb_len})
    if (idx+1)%20==0: print(f"  {idx+1}/{len(cifs)} done...", flush=True)

results.sort(key=lambda x: x["sequence_recovery"])
json.dump(results, open(OUT/"all_sequences.json","w"), indent=2)
with open(OUT/"minibinder_sequences.fasta","w") as f:
    for i,r in enumerate(results):
        f.write(f">rank_{i+1}_{r['backbone']}_rec{r['sequence_recovery']:.3f}\n{r['minibinder_sequence']}\n")
print(f"MPNN done: {len(results)} sequences → {OUT}")
PYEOF

cat > $PIPELINE/build_boltz_inputs.py << 'PYEOF'
import json, os
from pathlib import Path

PIPELINE = Path("/workspace")
MPNN_SUBDIR = os.environ.get("MPNN_OUTPUT_SUBDIR", "mpnn")
BOLTZ_IN_NAME = os.environ.get("BOLTZ_INPUT_DIR", "boltz_inputs")
MPNN_OUT = PIPELINE / "outputs" / MPNN_SUBDIR
BOLTZ_IN = PIPELINE / BOLTZ_IN_NAME
BOLTZ_IN.mkdir(parents=True, exist_ok=True)

mpnn_data = json.load(open(MPNN_OUT/"all_sequences.json"))
mpnn_data.sort(key=lambda x: x["sequence_recovery"])
seen, top_designs = set(), []
for r in mpnn_data:
    if r["backbone"] not in seen:
        seen.add(r["backbone"]); top_designs.append(r)
    if len(top_designs) >= int(os.environ.get("TOP_N", 50)): break

print(f"Building {len(top_designs)} Boltz-2 inputs in {BOLTZ_IN_NAME}...")
design_list = []
for i, r in enumerate(top_designs):
    bb, full_seq = r["backbone"], r["full_sequence"]
    if len(full_seq) < 171:
        print(f"  Skip {bb}: short sequence ({len(full_seq)})"); continue
    rank, rec = i + 1, r["sequence_recovery"]
    name = f"rank{rank:02d}_{bb}_rec{rec:.3f}"
    yaml = f"""sequences:
  - protein:
      id: A
      sequence: "{full_seq}"
      msa: empty
"""
    (BOLTZ_IN / f"{name}.yaml").write_text(yaml)
    design_list.append({"rank":rank,"backbone":bb,"minibinder_sequence":r["minibinder_sequence"],
                        "sequence_recovery":rec,
                        "mb_start":r.get("mb_start",116),
                        "mb_end":r.get("mb_end",170),
                        "mb_len":r.get("mb_len",55)})

meta_name = "top_designs_round1.json" if "round1" in BOLTZ_IN_NAME else "top_designs.json"
json.dump({"designs":design_list}, open(PIPELINE/meta_name,"w"), indent=2)
print(f"Wrote {len(design_list)} YAML files")
PYEOF

# ── Helper: MPNN → Boltz → score for one phase ───────────────────
run_validation_phase() {
    local PHASE_LABEL=$1
    local RFD3_SUBDIR=$2
    local MPNN_SUBDIR=$3
    local BOLTZ_IN_DIR=$4
    local BOLTZ_OUT_DIR=$5
    local RESULTS_FILE=$6
    local COPY_CIFS=${7:-yes}

    echo ""
    echo "=== $PHASE_LABEL: ProteinMPNN ==="
    docker run --rm --gpus all \
        -v $PIPELINE:/workspace \
        -e N_MPNN_SEQS=$N_MPNN_SEQS \
        -e RFD3_ACTIVE_SUBDIR=$RFD3_SUBDIR \
        -e MPNN_OUTPUT_SUBDIR=$MPNN_SUBDIR \
        -e FOUNDRY_CHECKPOINT_DIRS=/weights \
        rosettacommons/foundry:latest \
        python3 /workspace/run_mpnn.py

    echo ""
    echo "=== $PHASE_LABEL: Boltz-2 ==="
    wait $BOLTZ_PID 2>/dev/null || pip install boltz[cuda] -U -q
    pip install -q 'networkx>=3.0' 'platformdirs>=3.0' 2>/dev/null || true

    docker run --rm --gpus all \
        -v $PIPELINE:/workspace \
        -e TOP_N=$TOP_N \
        -e MPNN_OUTPUT_SUBDIR=$MPNN_SUBDIR \
        -e BOLTZ_INPUT_DIR=$BOLTZ_IN_DIR \
        -e FOUNDRY_CHECKPOINT_DIRS=/weights \
        rosettacommons/foundry:latest \
        python3 /workspace/build_boltz_inputs.py

    N_YAML=$(ls $PIPELINE/$BOLTZ_IN_DIR/*.yaml 2>/dev/null | wc -l)
    echo "Running Boltz-2 on $N_YAML sequences → $BOLTZ_OUT_DIR ..."
    $HOME/.local/bin/boltz predict $PIPELINE/$BOLTZ_IN_DIR \
        --out_dir $PIPELINE/$BOLTZ_OUT_DIR \
        --devices 1 --num_workers 2 --override 2>&1 | \
        grep -E "Predicting|Saving|Done|Error" | head -20

    echo ""
    echo "=== $PHASE_LABEL: global assembly RMSD scoring ==="
    local SCORE_ARGS=(
        --pipeline-dir /workspace
        --mode single_chain
        --plddt-threshold 0.60
        --global-rmsd-threshold 20.0
        --rfd3-subdir "$RFD3_SUBDIR"
        --mpnn-subdir "$MPNN_SUBDIR"
        --boltz-subdir "$BOLTZ_OUT_DIR"
        --results-file "$RESULTS_FILE"
    )
    if [ "$COPY_CIFS" = "no" ]; then
        SCORE_ARGS+=(--no-copy-cifs)
    else
        SCORE_ARGS+=(--n-top-cifs 10)
    fi

    docker run --rm --gpus all \
        -v $PIPELINE:/workspace \
        -e FOUNDRY_CHECKPOINT_DIRS=/weights \
        rosettacommons/foundry:latest \
        python3 /workspace/binary_antibodies/scoring.py "${SCORE_ARGS[@]}"
}

# ── Step 1: RFdiffusion3 round 1 ─────────────────────────────────
echo ""
echo "=== Step 1: RFdiffusion3 round 1 ==="
docker run --rm --gpus all \
    -v $PIPELINE:/workspace \
    -e N_DESIGNS=$N_DESIGNS -e BATCH_SIZE=$BATCH_SIZE \
    -e MB_LENGTH_RANGE=$MB_LENGTH_RANGE \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_rfd3.py
echo "RFd3 round 1: $(ls $PIPELINE/outputs/rfd3/mb_*.cif 2>/dev/null | wc -l) designs"

# ── Round 1 validation (always; feeds round-2 template selection) ─
if [ "$RFD3_ROUNDS" = "2" ]; then
    run_validation_phase \
        "Round 1 validation" \
        rfd3 mpnn_round1 \
        boltz_inputs_round1 boltz_outputs_round1 \
        round1_final_results.json no

    # ── Step 2: RFdiffusion3 round 2 (partial diffusion on top global RMSD) ─
    echo ""
    echo "=== Step 2: RFdiffusion3 round 2 (top global RMSD round-1 RFd3 CIFs) ==="
    docker run --rm --gpus all \
        -v $PIPELINE:/workspace \
        -e RFD3_ROUND2_TEMPLATES=$RFD3_ROUND2_TEMPLATES \
        -e RFD3_ROUND2_DESIGNS_PER_TEMPLATE=$RFD3_ROUND2_DESIGNS_PER_TEMPLATE \
        -e RFD3_ROUND2_BATCH_SIZE=$RFD3_ROUND2_BATCH_SIZE \
        -e RFD3_PARTIAL_T=$RFD3_PARTIAL_T \
        -e RFD3_ROUND2_VL_CONTEXT=$RFD3_ROUND2_VL_CONTEXT \
        -e RFD3_ROUND2_TEMPLATE_BACKBONES="$RFD3_ROUND2_TEMPLATE_BACKBONES" \
        -e RFD3_ROUND1_RESULTS=/workspace/final/round1_final_results.json \
        -e FOUNDRY_CHECKPOINT_DIRS=/weights \
        rosettacommons/foundry:latest \
        python3 /workspace/rfd3_round2.py
    echo "RFd3 round 2: $(ls $PIPELINE/outputs/rfd3_round2/r2_*.cif 2>/dev/null | wc -l) designs"

    # ── Round 2 validation (final results) ───────────────────────
    run_validation_phase \
        "Round 2 validation" \
        rfd3_round2 mpnn \
        boltz_inputs boltz_outputs \
        final_results.json yes
else
    # ── Single-round validation ──────────────────────────────────
    run_validation_phase \
        "Validation" \
        rfd3 mpnn \
        boltz_inputs boltz_outputs \
        final_results.json yes
fi

# ── Push to GitHub ────────────────────────────────────────────────
echo ""
echo "=== Push to GitHub ==="
export RESULTS_GITHUB_DIR=${RESULTS_GITHUB_DIR:-pipeline_results/v5_two_round_refine}
export GITHUB_BRANCH=${GITHUB_BRANCH:-cursor/conditional-nanobody-design-992c}
bash $PIPELINE/push_results.sh
echo ""
echo "===== Complete: $(date) ====="
