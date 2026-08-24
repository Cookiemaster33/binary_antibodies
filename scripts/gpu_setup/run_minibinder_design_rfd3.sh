#!/bin/bash
# ============================================================
# run_minibinder_design_rfd3.sh
# -----------------------------
# Minibinder design pipeline using RFdiffusion3 (rc-foundry Docker)
#
# Inputs (expected in ~/pipeline/inputs/):
#   trastuzumab_VL.pdb   — VL1 target structure
#   design_config.py     — rfd3 Python config (written by this script)
#
# Outputs in ~/pipeline/outputs/:
#   rfd3/                — designed backbones (CIF format)
#   mpnn/                — designed sequences
#   summary.csv          — ranked results
# ============================================================
set -eo pipefail

LOG="$HOME/pipeline/run.log"
exec > >(tee -a "$LOG") 2>&1

PIPELINE="$HOME/pipeline"
INPUT_PDB="$PIPELINE/inputs/trastuzumab_VL.pdb"
RFD3_OUT="$PIPELINE/outputs/rfd3"
MPNN_OUT="$PIPELINE/outputs/mpnn"
AF2_OUT="$PIPELINE/outputs/af2"
N_DESIGNS="${N_DESIGNS:-200}"
N_MPNN_SEQS="${N_MPNN_SEQS:-8}"
BATCH_SIZE="${BATCH_SIZE:-10}"   # designs per GPU batch

sudo chmod 666 /var/run/docker.sock 2>/dev/null || true

echo "===== Minibinder design run (rfd3): $(date) ====="
echo "N_DESIGNS=$N_DESIGNS  N_MPNN_SEQS=$N_MPNN_SEQS"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"

mkdir -p "$RFD3_OUT" "$MPNN_OUT" "$AF2_OUT"

# ── Step 1: Write rfd3 design Python script ──────────────────
cat > "$PIPELINE/run_rfd3.py" << 'PYEOF'
"""
RFdiffusion3 binder design: generate minibinders targeting VL1 FR2 hotspot.
"""
import os, json
from pathlib import Path
import numpy as np
import biotite.structure.io.pdb as pdb_io
from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine
from atomworks.io.utils.io_utils import to_cif_file

INPUT_PDB  = "/workspace/inputs/trastuzumab_VL.pdb"
OUT_DIR    = Path("/workspace/outputs/rfd3")
N_DESIGNS  = int(os.environ.get("N_DESIGNS", 200))
BATCH_SIZE = int(os.environ.get("BATCH_SIZE", 10))
N_BATCHES  = max(1, N_DESIGNS // BATCH_SIZE)

# VL1 FR2 hotspot residues (Chothia numbering, chain B in the domain PDB)
# Residues: 35-39 (FR2 core), 44-47 (FR2/CDR-L2 junction), 98 (FR4)
HOTSPOT_RESIDUES = [35, 36, 37, 38, 39, 44, 45, 46, 47, 98]

print(f"Generating {N_DESIGNS} minibinder designs ({N_BATCHES} batches × {BATCH_SIZE})")
print(f"Target: VL1 FR2 hotspot residues {HOTSPOT_RESIDUES}")

# Load target structure
from atomworks.io.utils.io_utils import from_pdb_file
target = from_pdb_file(INPUT_PDB)
print(f"Target loaded: {len(target)} atoms")

config = RFD3InferenceConfig(
    specification={
        # Binder design: design a protein that binds to the target
        'target': target,
        'length': 60,          # 60-residue minibinder
        'select_hotspots': {
            'chain': 'B',
            'res_ids': HOTSPOT_RESIDUES,
        },
        'extra': {},
    },
    diffusion_batch_size=BATCH_SIZE,
)

model = RFD3InferenceEngine(**config)

all_outputs = {}
for batch_idx in range(N_BATCHES):
    print(f"  Batch {batch_idx+1}/{N_BATCHES} ...", flush=True)
    outputs = model.run(inputs=None, out_dir=None, n_batches=1)
    all_outputs.update({f"batch{batch_idx}_{k}": v for k, v in outputs.items()})

# Save outputs as CIF files
saved = 0
for key, designs in all_outputs.items():
    for i, design in enumerate(designs):
        out_path = OUT_DIR / f"minibinder_{key}_{i:03d}.cif"
        to_cif_file(design.atom_array, str(out_path))
        saved += 1

print(f"Saved {saved} designs to {OUT_DIR}")
PYEOF

# ── Step 2: Run rfd3 in Docker ───────────────────────────────
echo ""
echo "=== Step 1/3: RFdiffusion3 binder design ==="

docker run --rm \
    --gpus all \
    -v "$HOME/.foundry:/root/.foundry" \
    -v "$PIPELINE:/workspace" \
    -e N_DESIGNS="$N_DESIGNS" \
    -e BATCH_SIZE="$BATCH_SIZE" \
    rosettacommons/foundry:latest \
    python3 /workspace/run_rfd3.py

N_GENERATED=$(ls "$RFD3_OUT"/*.cif 2>/dev/null | wc -l)
echo "  RFdiffusion3 done. $N_GENERATED backbones generated."

# ── Step 3: ProteinMPNN via foundry Docker ───────────────────
echo ""
echo "=== Step 2/3: ProteinMPNN sequence design ==="

cat > "$PIPELINE/run_mpnn.py" << 'PYEOF'
"""Run ProteinMPNN on all rfd3-generated backbones."""
import os, json
from pathlib import Path
from mpnn.inference_engines.mpnn import MPNNInferenceEngine
from atomworks.io.utils.io_utils import from_cif_file

RFD3_OUT   = Path("/workspace/outputs/rfd3")
MPNN_OUT   = Path("/workspace/outputs/mpnn")
N_SEQS     = int(os.environ.get("N_MPNN_SEQS", 8))

MPNN_OUT.mkdir(parents=True, exist_ok=True)

engine = MPNNInferenceEngine(
    model_type="protein_mpnn",
    is_legacy_weights=True,
    out_directory=None,
    write_structures=False,
    write_fasta=False,
)

all_results = []
cif_files = sorted(RFD3_OUT.glob("*.cif"))
print(f"Running MPNN on {len(cif_files)} backbones, {N_SEQS} seqs each...")

for cif_path in cif_files:
    atom_array = from_cif_file(str(cif_path))
    results = engine.run(
        inputs=atom_array,
        input_configs=[{"batch_size": N_SEQS, "remove_waters": True}],
    )
    for result in results:
        all_results.append({
            "backbone": cif_path.stem,
            "sequence": result.sequence,
            "mpnn_score": float(result.score) if hasattr(result, "score") else 0.0,
        })

# Save as FASTA + JSON
all_results.sort(key=lambda x: -x["mpnn_score"])
fasta_path = MPNN_OUT / "all_sequences.fasta"
with open(fasta_path, "w") as f:
    for rank, r in enumerate(all_results):
        f.write(f">rank_{rank+1}__backbone_{r['backbone']}__score_{r['mpnn_score']:.4f}\n")
        f.write(r["sequence"] + "\n")

json_path = MPNN_OUT / "all_sequences.json"
with open(json_path, "w") as f:
    json.dump(all_results, f, indent=2)

print(f"Saved {len(all_results)} sequences → {fasta_path}")
print(f"Top sequence: {all_results[0]['sequence'][:40]}...  score={all_results[0]['mpnn_score']:.4f}")
PYEOF

docker run --rm \
    --gpus all \
    -v "$HOME/.foundry:/root/.foundry" \
    -v "$PIPELINE:/workspace" \
    -e N_MPNN_SEQS="$N_MPNN_SEQS" \
    rosettacommons/foundry:latest \
    python3 /workspace/run_mpnn.py

echo "  ProteinMPNN done."

# ── Step 4: Quick MPNN-score ranking ─────────────────────────
echo ""
echo "=== Step 3/3: Ranking by MPNN score ==="

python3 - << 'PYEOF'
import json, csv, os
from pathlib import Path

results_file = os.path.expanduser("~/pipeline/outputs/mpnn/all_sequences.json")
csv_out = os.path.expanduser("~/pipeline/outputs/summary.csv")

if not os.path.exists(results_file):
    print("No MPNN results yet.")
    exit(0)

with open(results_file) as f:
    results = json.load(f)

results.sort(key=lambda x: -x.get("mpnn_score", 0))
with open(csv_out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["rank", "backbone", "sequence", "mpnn_score"])
    w.writeheader()
    for i, r in enumerate(results):
        w.writerow({"rank": i+1, **r})

print(f"Ranked {len(results)} designs → {csv_out}")
print("Top 5:")
for r in results[:5]:
    print(f"  {r['backbone']:40s}  score={r.get('mpnn_score',0):.4f}")
PYEOF

echo ""
echo "===== Pipeline complete: $(date) ====="
ls -lh "$PIPELINE/outputs/"
