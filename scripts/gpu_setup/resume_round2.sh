#!/bin/bash
# Resume round-2 from MPNN after fixing multi-chain MPNN bug.
set -eo pipefail
PIPELINE=/home/ubuntu/pipeline
export RFD3_ACTIVE_SUBDIR=rfd3_round2
export MPNN_OUTPUT_SUBDIR=mpnn
export N_MPNN_SEQS=${N_MPNN_SEQS:-8}
export TOP_N=${TOP_N:-50}
export RESULTS_GITHUB_DIR=${RESULTS_GITHUB_DIR:-pipeline_results/v5_two_round_refine}
export GITHUB_BRANCH=${GITHUB_BRANCH:-cursor/conditional-nanobody-design-992c}
export GITHUB_TOKEN=${GITHUB_TOKEN:-}

echo "=== Resume round 2: MPNN → Boltz → score → push ==="

docker run --rm --gpus all \
    -v $PIPELINE:/workspace \
    -e N_MPNN_SEQS=$N_MPNN_SEQS \
    -e RFD3_ACTIVE_SUBDIR=$RFD3_ACTIVE_SUBDIR \
    -e MPNN_OUTPUT_SUBDIR=$MPNN_OUTPUT_SUBDIR \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_mpnn.py

docker run --rm --gpus all \
    -v $PIPELINE:/workspace \
    -e TOP_N=$TOP_N \
    -e MPNN_OUTPUT_SUBDIR=mpnn \
    -e BOLTZ_INPUT_DIR=boltz_inputs \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/build_boltz_inputs.py

N_YAML=$(ls $PIPELINE/boltz_inputs/*.yaml 2>/dev/null | wc -l)
echo "Running Boltz-2 on $N_YAML sequences..."
$HOME/.local/bin/boltz predict $PIPELINE/boltz_inputs \
    --out_dir $PIPELINE/boltz_outputs \
    --devices 1 --num_workers 2 --override

docker run --rm --gpus all \
    -v $PIPELINE:/workspace \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/binary_antibodies/scoring.py \
        --pipeline-dir /workspace \
        --mode single_chain \
        --rfd3-subdir rfd3_round2 \
        --mpnn-subdir mpnn \
        --boltz-subdir boltz_outputs \
        --results-file final_results.json \
        --n-top-cifs 10

bash $PIPELINE/push_results.sh
echo "===== Resume complete: $(date) ====="
