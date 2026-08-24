#!/bin/bash
# ============================================================
# setup_pipeline_rfd3.sh
# ----------------------
# Setup for RFdiffusion3 (rc-foundry) via Docker + ProteinMPNN
# on a Lambda Cloud GPU instance.
#
# Much simpler than rfd1 — all dependencies are in the Docker image.
# ============================================================
set -eo pipefail

LOG="$HOME/setup_pipeline.log"
exec > >(tee -a "$LOG") 2>&1

echo "===== Pipeline setup started: $(date) ====="
echo "Host: $(hostname)  GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null)"

# ── 1. Docker permissions ─────────────────────────────────────
echo ""
echo "--- [1/4] Docker ---"
sudo chmod 666 /var/run/docker.sock
docker --version

# ── 2. Pull foundry image (if not already present) ───────────
echo ""
echo "--- [2/4] Foundry Docker image ---"
if docker images rosettacommons/foundry | grep -q "foundry"; then
    echo "  Image already present."
    docker images rosettacommons/foundry
else
    echo "  Pulling rosettacommons/foundry:latest ..."
    docker pull rosettacommons/foundry:latest
fi

# ── 3. Download rfd3 model weights ───────────────────────────
echo ""
echo "--- [3/4] Downloading rfd3 model weights ---"
mkdir -p "$HOME/.foundry/checkpoints"

docker run --rm \
    --gpus all \
    -v "$HOME/.foundry:/root/.foundry" \
    rosettacommons/foundry:latest \
    foundry install base-models 2>&1 | tail -20

echo "  Weights:"
ls "$HOME/.foundry/checkpoints/" 2>/dev/null || echo "  (check ~/.foundry/checkpoints)"

# ── 4. Workspace ──────────────────────────────────────────────
echo ""
echo "--- [4/4] Workspace ---"
mkdir -p "$HOME/pipeline/inputs" \
         "$HOME/pipeline/outputs/rfd3" \
         "$HOME/pipeline/outputs/mpnn" \
         "$HOME/pipeline/outputs/af2"

echo ""
echo "===== Setup complete: $(date) ====="
echo ""
echo "Test with:"
echo "  docker run --rm --gpus all \\"
echo "    -v \$HOME/.foundry:/root/.foundry \\"
echo "    -v \$HOME/pipeline:/workspace \\"
echo "    rosettacommons/foundry:latest \\"
echo "    foundry list-installed"
