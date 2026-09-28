#!/bin/bash
# ============================================================
# setup_pipeline.sh
# -----------------
# One-shot setup script for an RFdiffusion + ProteinMPNN +
# ColabFold (AF2) pipeline on a fresh Lambda Cloud GPU instance.
#
# Run as: ubuntu@<instance-ip>
# Estimated time: 15-25 min (mostly model downloads)
# ============================================================
set -euo pipefail

LOG="$HOME/setup_pipeline.log"
exec > >(tee -a "$LOG") 2>&1

echo "===== Pipeline setup started: $(date) ====="
echo "Host: $(hostname)  GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'no GPU detected')"

# ── 1. System packages ────────────────────────────────────────
echo ""
echo "--- [1/6] System packages ---"
sudo apt-get update -q
sudo apt-get install -y -q git curl wget rsync tree tmux pigz

# ── 2. Mamba / conda ─────────────────────────────────────────
echo ""
echo "--- [2/6] Mamba ---"
if ! command -v mamba &>/dev/null; then
    wget -q "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" \
        -O /tmp/miniforge.sh
    bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
    rm /tmp/miniforge.sh
fi
source "$HOME/miniforge3/etc/profile.d/conda.sh" || source "$HOME/anaconda3/etc/profile.d/conda.sh"
echo "conda: $(conda --version)"

# ── 3. RFdiffusion ───────────────────────────────────────────
echo ""
echo "--- [3/6] RFdiffusion ---"
RF_DIR="$HOME/RFdiffusion"
if [ ! -d "$RF_DIR" ]; then
    git clone https://github.com/RosettaCommons/RFdiffusion.git "$RF_DIR"
fi
cd "$RF_DIR"

if ! conda env list | grep -q "SE3nv"; then
    mamba env create -f env/SE3nv.yml -n SE3nv -y
fi
conda activate SE3nv

# Install RFdiffusion as package
pip install -e . -q

# Download models
echo "  Downloading RFdiffusion models (~2.5 GB) ..."
mkdir -p "$HOME/rf_models"
bash "$RF_DIR/scripts/download_models.sh" "$HOME/rf_models"
echo "  Models:"
ls "$HOME/rf_models"

# ── 4. ProteinMPNN ───────────────────────────────────────────
echo ""
echo "--- [4/6] ProteinMPNN ---"
MPNN_DIR="$HOME/ProteinMPNN"
if [ ! -d "$MPNN_DIR" ]; then
    git clone https://github.com/dauparas/ProteinMPNN.git "$MPNN_DIR"
fi
# ProteinMPNN runs in the same SE3nv env (torch already installed)
pip install -q biotite

# ── 5. ColabFold / AF2-Multimer for scoring ──────────────────
echo ""
echo "--- [5/6] ColabFold (AF2 scoring) ---"
if ! conda env list | grep -q "colabfold"; then
    mamba create -n colabfold python=3.10 -y
fi
conda activate colabfold
pip install -q "colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold"
# Download AF2 weights
python -c "from colabfold.download import download_alphafold_params; download_alphafold_params('multimer', '$HOME/af2_params')" || \
    python -m colabfold.download --help 2>&1 | head -5
conda activate SE3nv

# ── 6. Workspace ──────────────────────────────────────────────
echo ""
echo "--- [6/6] Workspace ---"
mkdir -p "$HOME/pipeline/inputs" "$HOME/pipeline/outputs/rfdiffusion" \
         "$HOME/pipeline/outputs/mpnn" "$HOME/pipeline/outputs/af2"
echo "  Directory tree:"
tree -d "$HOME/pipeline"

echo ""
echo "===== Setup complete: $(date) ====="
echo "Conda envs: $(conda env list)"
echo ""
echo "To run RFdiffusion:"
echo "  conda activate SE3nv"
echo "  cd $RF_DIR"
echo "  python scripts/run_inference.py <args>"
