#!/usr/bin/env bash
# ============================================================
# setup_boltz.sh
# ---------------
# Install Boltz-1 and pre-fetch model weights on a Lambda GPU
# instance.  This script is run once before the CD3e screen.
# ============================================================
set -euo pipefail

LOG="$HOME/setup_boltz.log"
exec > >(tee -a "$LOG") 2>&1

echo "===== Boltz setup started: $(date) ====="
echo "Host: $(hostname)  GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null)"

# ── 1. Install Boltz ────────────────────────────────────────
echo ""
echo "--- [1/3] Install Boltz ---"
if command -v boltz &>/dev/null; then
    echo "  Boltz already installed: $(boltz --version 2>/dev/null || echo ok)"
else
    pip install boltz -U -q 2>&1 | tail -5 \
      || pip install --break-system-packages boltz -U -q 2>&1 | tail -5
    echo "  Installed: $(boltz --version 2>/dev/null || $HOME/.local/bin/boltz --version 2>/dev/null || echo ok)"
fi

# Ensure the boltz binary is on PATH
BOLTZ_BIN=$(command -v boltz 2>/dev/null || echo "$HOME/.local/bin/boltz")
export PATH="$HOME/.local/bin:$PATH"

# ── 2. Pre-fetch model weights ───────────────────────────────
echo ""
echo "--- [2/3] Pre-fetch Boltz weights ---"
WEIGHTS_DIR="$HOME/.boltz"
if [ -d "$WEIGHTS_DIR" ] && [ "$(ls -A "$WEIGHTS_DIR" 2>/dev/null)" ]; then
    echo "  Weights already present in $WEIGHTS_DIR"
    ls "$WEIGHTS_DIR"
else
    echo "  Triggering weight download via a minimal predict run..."
    TMPDIR_PRED=$(mktemp -d)
    cat > "$TMPDIR_PRED/dummy.yaml" << 'YAML'
sequences:
  - protein:
      id: A
      sequence: MGSSHHHHHHSQDPMSSYQHFMKLNLNPVVAALN
      msa: empty
YAML
    # Run predict – will fail on single-chain "complex" but weights get cached first
    $BOLTZ_BIN predict "$TMPDIR_PRED" \
        --out_dir "$TMPDIR_PRED/out" \
        --accelerator gpu \
        --devices 1 \
        --override 2>&1 | tail -10 || true
    rm -rf "$TMPDIR_PRED"
    echo "  Weight pre-fetch done."
fi

# ── 3. Workspace ────────────────────────────────────────────
echo ""
echo "--- [3/3] Workspace ---"
mkdir -p "$HOME/pipeline/cd3e_screen" \
         "$HOME/pipeline/cd3e_screen_results"

echo ""
echo "===== Setup complete: $(date) ====="
