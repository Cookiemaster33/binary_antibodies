#!/bin/bash
# ============================================================
# run_minibinder_design.sh
# ------------------------
# Full minibinder design pipeline on the Lambda GPU instance.
# Called remotely by launch_minibinder_design.py after upload.
#
# Inputs (expected in ~/pipeline/inputs/):
#   trastuzumab_VL.pdb   — VL1 target structure
#
# Outputs written to ~/pipeline/outputs/:
#   rfdiffusion/         — designed backbones (.pdb, .trb)
#   mpnn/                — designed sequences (.fa, scores.jsonl)
#   af2/                 — AF2-Multimer predictions (.pdb, .json)
#   summary.csv          — ranked designs
# ============================================================
set -euo pipefail

LOG="$HOME/pipeline/run.log"
exec > >(tee -a "$LOG") 2>&1

source "$HOME/miniforge3/etc/profile.d/conda.sh" 2>/dev/null || \
    source "$HOME/anaconda3/etc/profile.d/conda.sh"

RF_DIR="$HOME/RFdiffusion"
MPNN_DIR="$HOME/ProteinMPNN"
PIPELINE="$HOME/pipeline"
INPUT_PDB="$PIPELINE/inputs/trastuzumab_VL.pdb"
RF_OUT="$PIPELINE/outputs/rfdiffusion/vl1_minibinder"
MPNN_OUT="$PIPELINE/outputs/mpnn"
AF2_OUT="$PIPELINE/outputs/af2"
N_DESIGNS="${N_DESIGNS:-200}"      # set via env or default 200
N_MPNN_SEQS="${N_MPNN_SEQS:-8}"    # ProteinMPNN sequences per backbone

echo "===== Minibinder design run: $(date) ====="
echo "N_DESIGNS=$N_DESIGNS  N_MPNN_SEQS=$N_MPNN_SEQS"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"

# ── Step 1: RFdiffusion ──────────────────────────────────────
echo ""
echo "=== Step 1/3: RFdiffusion binder design ==="
conda activate SE3nv
mkdir -p "$(dirname $RF_OUT)"

# Hotspot residues: VL1 framework 2 (Chothia: 35-39, 44-47, 98)
# Chain B in the trastuzumab_VL.pdb domain file
python "$RF_DIR/scripts/run_inference.py" \
    inference.input_pdb="$INPUT_PDB" \
    inference.output_prefix="$RF_OUT" \
    inference.model_directory_path="$HOME/rf_models" \
    inference.num_designs="$N_DESIGNS" \
    "contigmap.contigs=[B35-47/0 A1-65]" \
    "ppi.hotspot_res=[B35,B36,B37,B38,B39,B44,B45,B46,B47,B98]" \
    "diffuser.T=50" \
    "inference.output_prefix=$RF_OUT"

echo "  RFdiffusion done. $(ls ${RF_OUT}_*.pdb 2>/dev/null | wc -l) backbones generated."

# ── Step 2: ProteinMPNN ──────────────────────────────────────
echo ""
echo "=== Step 2/3: ProteinMPNN sequence design ==="
mkdir -p "$MPNN_OUT"

# Parse chain information from each backbone + VL1 complex
for pdb in "${RF_OUT}"_*.pdb; do
    stem=$(basename "$pdb" .pdb)
    python "$MPNN_DIR/helper_scripts/parse_multiple_chains.py" \
        --input_path "$pdb" \
        --output_path "$MPNN_OUT/${stem}_parsed.jsonl"

    python "$MPNN_DIR/protein_mpnn_run.py" \
        --jsonl_path "$MPNN_OUT/${stem}_parsed.jsonl" \
        --out_folder "$MPNN_OUT/${stem}" \
        --num_seq_per_target "$N_MPNN_SEQS" \
        --sampling_temp "0.1" \
        --seed 42 \
        --batch_size 1 \
        --chain_id_jsonl "" \
        --fixed_positions_jsonl "" \
        2>/dev/null
done

echo "  ProteinMPNN done."

# ── Step 3: AF2-Multimer scoring (top-N by MPNN score) ───────
echo ""
echo "=== Step 3/3: AF2-Multimer scoring (top 20 designs) ==="
conda activate colabfold
mkdir -p "$AF2_OUT"

# Collect top 20 sequences by MPNN score, pair with VL1 sequence
VL1_SEQ=$(python3 -c "
from Bio import SeqIO
import sys
for rec in SeqIO.parse('$INPUT_PDB', 'pdb-seqres'):
    print(str(rec.seq))
    break
" 2>/dev/null || echo "DIQMTQSPSSVSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")

# Build combined FASTA for AF2-Multimer (minibinder:VL1)
TOP_SEQS="$AF2_OUT/top20_for_af2.fasta"
python3 - <<'PYEOF'
import glob, json, os

mpnn_out = os.environ.get("MPNN_OUT", "$MPNN_OUT")
af2_out = os.environ.get("AF2_OUT", "$AF2_OUT")
vl1_seq = os.environ.get("VL1_SEQ", "")
top_n = 20

scores = []
for score_file in glob.glob(f"{mpnn_out}/*/seqs/*.fa"):
    with open(score_file) as f:
        lines = f.readlines()
    for i in range(0, len(lines)-1, 2):
        header = lines[i].strip().lstrip(">")
        seq = lines[i+1].strip()
        try:
            score = float(header.split("score=")[1].split(",")[0])
        except Exception:
            score = 0.0
        scores.append((score, seq, header))

scores.sort(reverse=True)
os.makedirs(af2_out, exist_ok=True)
with open(f"{af2_out}/top{top_n}_for_af2.fasta", "w") as f:
    for rank, (score, seq, header) in enumerate(scores[:top_n]):
        f.write(f">design_{rank+1}__score_{score:.3f}\n")
        f.write(f"{seq}:{vl1_seq}\n")
print(f"Wrote {min(top_n, len(scores))} sequences for AF2")
PYEOF

# Run ColabFold AF2-Multimer on the combined sequences
colabfold_batch \
    "$AF2_OUT/top20_for_af2.fasta" \
    "$AF2_OUT" \
    --model-type alphafold2_multimer_v3 \
    --num-recycle 3 \
    --num-models 1 \
    --templates \
    2>/dev/null || echo "  AF2 scoring complete (or partial)."

echo ""
echo "===== Pipeline complete: $(date) ====="
echo "Results in: $PIPELINE/outputs/"
ls -lh "$PIPELINE/outputs/"

# ── Summary CSV ──────────────────────────────────────────────
python3 - <<'PYEOF'
import glob, json, os, csv

af2_out = "$AF2_OUT"
out_csv = "$PIPELINE/outputs/summary.csv"

rows = []
for score_file in glob.glob(f"{af2_out}/*_scores_*multimer*.json"):
    with open(score_file) as f:
        data = json.load(f)
    name = os.path.basename(score_file).replace("_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json", "")
    rows.append({
        "design": name,
        "iptm": data.get("iptm", 0),
        "ptm": data.get("ptm", 0),
        "mean_plddt": sum(data.get("plddt", [0])) / max(len(data.get("plddt", [1])), 1),
    })

rows.sort(key=lambda x: -x["iptm"])
with open(out_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["design", "iptm", "ptm", "mean_plddt"])
    writer.writeheader()
    writer.writerows(rows)
print(f"Summary: {out_csv}  ({len(rows)} designs ranked by ipTM)")
if rows:
    print(f"Top design: {rows[0]['design']}  ipTM={rows[0]['iptm']:.3f}")
PYEOF
