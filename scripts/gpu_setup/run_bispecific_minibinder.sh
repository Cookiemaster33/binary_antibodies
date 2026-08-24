#!/bin/bash
# ============================================================
# run_bispecific_minibinder.sh
# Bispecific bridging minibinder design:
#   VH1-CH1-face ↔ Nanobody CDRs
#
# Input: vh1_nanobody_design_target.pdb
#   Chain A: VH1 (115 res) — hotspot: CH1-contact residues
#   Chain B: Nanobody in CH1 slot (115 res) — hotspot: CDR1+2+3
# Contig: A1-115,55,B1-115 → 285 res (VH1 + 55-res MB + Nanobody)
# ============================================================
set -eo pipefail
sudo chmod 666 /var/run/docker.sock 2>/dev/null || true
PIPELINE=/home/ubuntu/pipeline
LOG=$PIPELINE/pipeline.log
exec > >(tee -a "$LOG") 2>&1

N_DESIGNS=${N_DESIGNS:-200}
BATCH_SIZE=${BATCH_SIZE:-10}
N_MPNN_SEQS=${N_MPNN_SEQS:-8}
N_BATCHES=$(( (N_DESIGNS + BATCH_SIZE - 1) / BATCH_SIZE ))

echo "===== Bispecific minibinder pipeline: $(date) ====="
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "Contig: A1-115,55,B1-115  (VH1 + 55-res minibinder + Nanobody)"
mkdir -p $PIPELINE/outputs/rfd3 $PIPELINE/outputs/mpnn $PIPELINE/outputs/final

# ── RFdiffusion3 ──────────────────────────────────────────────────
cat > $PIPELINE/run_rfd3.py << 'PYEOF'
import os, torch
from pathlib import Path
from atomworks.io.utils.io_utils import to_cif_file
from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine
from rfd3.inference.input_parsing import DesignInputSpecification
torch.set_float32_matmul_precision('high')

INPUT  = "/workspace/inputs/vh1_nanobody_design_target.pdb"
OUT    = Path("/workspace/outputs/rfd3")
N      = int(os.environ.get("N_DESIGNS", 200))
BATCH  = int(os.environ.get("BATCH_SIZE", 10))
NBATCH = max(1, N // BATCH)
OUT.mkdir(parents=True, exist_ok=True)

# VH1-CH1-face hotspots + Nanobody CDR hotspots
HOTSPOTS = (
    # VH1 residues that contact CH1 (the surface the minibinder must bind on VH1 side)
    "A11,A12,A13,A14,A15,A39,A40,A41,"
    "A79,A80,A81,A82,A83,A84,A85,"
    "A103,A104,A105,A106,A107,A108,A109,A110,A111,A112,A113,A114,A115,"
    # Nanobody CDR1 + CDR2 + CDR3 (the surface the minibinder must bind on nanobody side)
    "B27,B28,B29,B30,B31,B32,B33,"
    "B52,B53,B54,B55,B56,B57,"
    "B99,B100,B101,B102,B103,B104,B105,B106,B107,B108,B109,B110,B111,B112"
)

print(f"=== RFd3 Bispecific Bridging Minibinder ===")
print(f"N={N}  batches={NBATCH}×{BATCH}")
print(f"Design: 55-res minibinder bridging VH1-CH1-face ↔ Nanobody CDRs")

spec = DesignInputSpecification.safe_init(
    input=INPUT,
    contig="A1-115,55,B1-115",
    select_hotspots=HOTSPOTS,
)
model = RFD3InferenceEngine(**RFD3InferenceConfig(diffusion_batch_size=BATCH))

saved = 0
for bi in range(NBATCH):
    print(f"  Batch {bi+1}/{NBATCH} ...", flush=True)
    for key, designs in model.run(inputs=spec, out_dir=None, n_batches=1).items():
        for i, d in enumerate(designs):
            to_cif_file(d.atom_array, str(OUT / f"mb_b{bi:03d}_{i:03d}.cif"))
            saved += 1
print(f"RFd3 done: {saved} designs")
print(f"Structure: A=VH1(1-115), MB=minibinder(116-170), B=Nanobody(171-285)")
PYEOF

docker run --rm --gpus all \
    -v $PIPELINE:/workspace \
    -e N_DESIGNS=$N_DESIGNS -e BATCH_SIZE=$BATCH_SIZE \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_rfd3.py

N_CIF=$(ls $PIPELINE/outputs/rfd3/mb_*.cif 2>/dev/null | wc -l)
echo "Generated: $N_CIF designs"

# ── ProteinMPNN ───────────────────────────────────────────────────
cat > $PIPELINE/run_mpnn.py << 'PYEOF'
import os, json
from pathlib import Path
from atomworks.io.utils.io_utils import load_any
from mpnn.inference_engines.mpnn import MPNNInferenceEngine

RFD3  = Path("/workspace/outputs/rfd3")
OUT   = Path("/workspace/outputs/mpnn")
NSEQS = int(os.environ.get("N_MPNN_SEQS", 8))
OUT.mkdir(parents=True, exist_ok=True)

engine = MPNNInferenceEngine(model_type="protein_mpnn", is_legacy_weights=True,
    out_directory=None, write_structures=False, write_fasta=False)

results = []
cifs = sorted(RFD3.glob("mb_*.cif"))
print(f"MPNN: {len(cifs)} × {NSEQS}  |  Designing: residues 116-170 (minibinder)")

for idx, cif in enumerate(cifs):
    raw = load_any(str(cif)); aa = raw[0] if hasattr(raw,"__getitem__") else raw
    chain = list(set(aa.chain_id))[0]
    # Design only minibinder (116-170); fix VH1 (1-115) and nanobody (171-285)
    designed = [f"{chain}{r}" for r in range(116,171) if r in set(aa.res_id)]
    result = engine.run(atom_arrays=[aa], input_dicts=[{
        "batch_size": NSEQS, "remove_waters": True, "designed_residues": designed}])
    for r in (result if isinstance(result,list) else [result]):
        seq = r.output_dict.get('designed_sequence','')
        rec = float(r.output_dict.get('sequence_recovery',0))
        mb_seq = seq[115:171] if len(seq)>=171 else seq[115:]
        results.append({"backbone":cif.stem,"full_sequence":seq,
                        "minibinder_sequence":mb_seq,"sequence_recovery":rec})
    if (idx+1)%20==0:
        print(f"  {idx+1}/{len(cifs)} done...", flush=True)

results.sort(key=lambda x: x["sequence_recovery"])
json.dump(results, open(OUT/"all_sequences.json","w"), indent=2)
fasta = OUT/"minibinder_sequences.fasta"
with open(fasta,"w") as f:
    for i,r in enumerate(results):
        f.write(f">rank_{i+1}_{r['backbone']}_rec{r['sequence_recovery']:.3f}\n{r['minibinder_sequence']}\n")
print(f"MPNN done: {len(results)} sequences → {fasta}")
if results:
    print(f"Top: {results[0]['minibinder_sequence'][:50]}  rec={results[0]['sequence_recovery']:.3f}")
PYEOF

docker run --rm --gpus all \
    -v $PIPELINE:/workspace \
    -e N_MPNN_SEQS=$N_MPNN_SEQS \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_mpnn.py

echo ""
echo "===== Pipeline complete: $(date) ====="
ls -lh $PIPELINE/outputs/mpnn/
