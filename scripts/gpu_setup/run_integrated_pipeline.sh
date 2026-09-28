#!/bin/bash
# ============================================================
# run_integrated_pipeline.sh
# --------------------------
# Full pipeline on one instance:
#   1. RFdiffusion3  — design minibinder bridging nanobody CDRs + CL domain
#   2. ProteinMPNN   — design sequences for the minibinder
#   3. ColabFold/AF2 — self-consistency validation (scRMSD + pLDDT)
#   4. Rank and save final results
#
# Input:
#   /workspace/inputs/nanobody_cl_design_target.pdb
#     chain B (1-115): 2Rs15d anti-HER2 nanobody
#     chain C (1-113): Trastuzumab CL domain (kicker context)
#
# Design target:
#   50-residue minibinder that bridges nanobody CDR face + CL domain
#   Contig: B1-115/0,50/C1-113
#   Hotspots: CDR1(B26-33) + CDR2(B50-57) + CDR3(B97-112)
#
# Validation:
#   AF2 monomer on the 50-residue minibinder sequence alone.
#   scRMSD = RMSD of AF2 predicted backbone vs RFd3 designed backbone.
#   Filter: scRMSD < 2.0 Å AND pLDDT > 60.
# ============================================================
set -eo pipefail

PIPELINE=/home/ubuntu/pipeline
LOG=$PIPELINE/pipeline.log
exec > >(tee -a "$LOG") 2>&1

N_DESIGNS=${N_DESIGNS:-200}
BATCH_SIZE=${BATCH_SIZE:-10}
N_MPNN_SEQS=${N_MPNN_SEQS:-8}
N_BATCHES=$(( (N_DESIGNS + BATCH_SIZE - 1) / BATCH_SIZE ))

sudo chmod 666 /var/run/docker.sock 2>/dev/null || true

echo "===== Integrated pipeline start: $(date) ====="
echo "N_DESIGNS=$N_DESIGNS  BATCH_SIZE=$BATCH_SIZE  N_MPNN_SEQS=$N_MPNN_SEQS"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"

mkdir -p $PIPELINE/outputs/rfd3 $PIPELINE/outputs/mpnn \
         $PIPELINE/outputs/af2  $PIPELINE/outputs/final

# ── 1. RFdiffusion3 ──────────────────────────────────────────────────
echo ""
echo "=== Step 1/4: RFdiffusion3 ==="

cat > $PIPELINE/run_rfd3.py << 'PYEOF'
"""
RFdiffusion3 — design minibinder bridging nanobody CDRs and CL domain.
Input: two-chain PDB (nanobody chain B + CL domain chain C)
Design: 50-residue minibinder inserted between chains
Hotspots: nanobody CDR face (CDR1+CDR2+CDR3)
"""
import os, torch
from pathlib import Path
from atomworks.io.utils.io_utils import to_cif_file
from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine
from rfd3.inference.input_parsing import DesignInputSpecification
torch.set_float32_matmul_precision('high')

INPUT  = "/workspace/inputs/nanobody_cl_design_target.pdb"
OUT    = Path("/workspace/outputs/rfd3")
N      = int(os.environ.get("N_DESIGNS", 200))
BATCH  = int(os.environ.get("BATCH_SIZE", 10))
NBATCH = max(1, N // BATCH)
OUT.mkdir(parents=True, exist_ok=True)

print(f"RFd3: {N} designs ({NBATCH} batches × {BATCH})")
print(f"Input: nanobody (chain B, 115 res) + CL domain (chain C, 113 res)")
print(f"Design: 50-res minibinder between them")

HOTSPOTS = (
    "B26,B27,B28,B29,B30,B31,B32,B33,"    # CDR1
    "B50,B51,B52,B53,B54,B55,B56,B57,"    # CDR2
    "B97,B98,B99,B100,B101,B102,B103,"    # CDR3
    "B104,B105,B106,B107,B108,B109,B110,B111,B112"
)

spec = DesignInputSpecification.safe_init(
    input=INPUT,
    contig="B1-115/0,50/C1-113",    # nanobody | minibinder (50 res) | CL domain
    select_hotspots=HOTSPOTS,
)

config = RFD3InferenceConfig(diffusion_batch_size=BATCH)
model  = RFD3InferenceEngine(**config)

saved = 0
for bi in range(NBATCH):
    print(f"  Batch {bi+1}/{NBATCH} ...", flush=True)
    for key, designs in model.run(inputs=spec, out_dir=None, n_batches=1).items():
        for i, d in enumerate(designs):
            to_cif_file(d.atom_array, str(OUT / f"mb_b{bi:03d}_{i:03d}.cif"))
            saved += 1

print(f"RFd3 done: {saved} designs")
print(f"CIF structure: chain B(1-115 nanobody) + minibinder(116-165) + chain C(166-278 CL)")
PYEOF

docker run --rm --gpus all \
    -v $PIPELINE:/workspace \
    -e N_DESIGNS=$N_DESIGNS -e BATCH_SIZE=$BATCH_SIZE \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_rfd3.py

N_CIF=$(ls $PIPELINE/outputs/rfd3/mb_*.cif 2>/dev/null | wc -l)
echo "  Generated: $N_CIF backbones"

# ── 2. ProteinMPNN ───────────────────────────────────────────────────
echo ""
echo "=== Step 2/4: ProteinMPNN ==="

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
print(f"MPNN: {len(cifs)} backbones × {NSEQS} seqs")
print(f"Designing: residues 116-165 (minibinder)")
print(f"Fixed:     residues 1-115 (nanobody) + 166-end (CL domain)")

for idx, cif_path in enumerate(cifs):
    raw = load_any(str(cif_path))
    aa  = raw[0] if hasattr(raw, '__getitem__') else raw
    chain = list(set(aa.chain_id))[0]
    all_res = sorted(set(aa.res_id))
    
    # Fix nanobody (1-115) and CL domain (166+), design minibinder (116-165)
    designed = [f"{chain}{r}" for r in range(116, 166) if r in set(aa.res_id)]
    fixed    = [f"{chain}{r}" for r in all_res if r not in range(116, 166)]
    
    result = engine.run(atom_arrays=[aa], input_dicts=[{
        "batch_size": NSEQS, "remove_waters": True,
        "fixed_residues": fixed,
        "designed_residues": designed,
    }])
    
    for r in (result if isinstance(result, list) else [result]):
        seq = r.output_dict.get('designed_sequence', '')
        rec = float(r.output_dict.get('sequence_recovery', 0))
        mb_seq = seq[115:165] if len(seq) >= 165 else seq[115:]
        results.append({"backbone": cif_path.stem, "full_sequence": seq,
                        "minibinder_sequence": mb_seq, "sequence_recovery": rec})
    
    if (idx+1) % 20 == 0:
        print(f"  {idx+1}/{len(cifs)} done...", flush=True)

results.sort(key=lambda x: x["sequence_recovery"])
json.dump(results, open(OUT / "all_sequences.json", "w"), indent=2)
fasta = OUT / "minibinder_sequences.fasta"
with open(fasta, "w") as f:
    for i, r in enumerate(results):
        f.write(f">rank_{i+1}_{r['backbone']}_rec{r['sequence_recovery']:.3f}\n{r['minibinder_sequence']}\n")
print(f"MPNN done: {len(results)} sequences → {fasta}")
PYEOF

docker run --rm --gpus all \
    -v $PIPELINE:/workspace \
    -e N_MPNN_SEQS=$N_MPNN_SEQS \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_mpnn.py

# ── 3. Install ColabFold for AF2 self-consistency ────────────────────
echo ""
echo "=== Step 3/4: Install ColabFold (if needed) ==="

if ! command -v ~/.local/bin/colabfold_batch &>/dev/null; then
    echo "  Installing ColabFold..."
    pip install -q "colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold"
    pip install -q --upgrade "jax[cuda12]"
    python3 -m colabfold.download 2>&1 | tail -3
    echo "  ColabFold installed."
else
    echo "  ColabFold already installed."
fi

# ── 4. AF2 self-consistency validation ──────────────────────────────
echo ""
echo "=== Step 4/4: AF2 Self-Consistency Validation ==="

cat > $PIPELINE/run_selfconsistency.py << 'PYEOF'
"""
AF2 self-consistency validation.

For each minibinder sequence (top 50 unique backbones × 1 best MPNN seq):
1. Run ColabFold AF2 MONOMER on the 50-residue minibinder sequence alone
2. Compare AF2-predicted backbone to the RFd3-designed backbone
3. Compute scRMSD (Cα RMSD after superposition)
4. Filter: scRMSD < 2.0 Å AND pLDDT > 60

This is the standard self-consistency protocol from the RFdiffusion paper.
"""
import json, subprocess, os, glob
from pathlib import Path
import numpy as np

MPNN_OUT = Path("/workspace/outputs/mpnn")
RFD3_OUT = Path("/workspace/outputs/rfd3")
AF2_OUT  = Path("/workspace/outputs/af2_sc")
FINAL    = Path("/workspace/outputs/final")
AF2_OUT.mkdir(parents=True, exist_ok=True)
FINAL.mkdir(parents=True, exist_ok=True)

results = json.load(open(MPNN_OUT / "all_sequences.json"))

# Take top 50 unique backbone designs (lowest recovery)
seen = {}
top_designs = []
for r in results:
    bb = r["backbone"]
    if bb not in seen:
        seen[bb] = r
        top_designs.append(r)
    if len(top_designs) >= 50:
        break

print(f"Self-consistency validation on {len(top_designs)} designs")

# Write FASTA for AF2 (monomer: just the 50-res minibinder)
sc_fasta = AF2_OUT / "sc_sequences.fasta"
with open(sc_fasta, "w") as f:
    for i, r in enumerate(top_designs):
        mb = r["minibinder_sequence"]
        f.write(f">{r['backbone']}\n{mb}\n")

print(f"Running AF2 monomer on {len(top_designs)} sequences...")
result = subprocess.run([
    os.path.expanduser("~/.local/bin/colabfold_batch"),
    str(sc_fasta), str(AF2_OUT),
    "--model-type", "alphafold2",          # monomer
    "--msa-mode", "single_sequence",
    "--num-recycle", "3",
    "--num-models", "1",
    "--use-gpu-relax",
], capture_output=False, text=True)

# Parse AF2 results and compute scRMSD
from atomworks.io.utils.io_utils import load_any
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import gzip, io

sc_results = []
for design in top_designs:
    bb = design["backbone"]
    mb_seq = design["minibinder_sequence"]
    
    # Find the AF2 output PDB for this design
    af2_pdbs = glob.glob(str(AF2_OUT / f"{bb}_relaxed_rank_001*.pdb"))
    if not af2_pdbs:
        af2_pdbs = glob.glob(str(AF2_OUT / f"{bb}*.pdb"))
    
    # Find AF2 score
    score_files = glob.glob(str(AF2_OUT / f"{bb}_scores_rank_001*.json"))
    plddt = 0.0
    ptm = 0.0
    if score_files:
        sc = json.load(open(score_files[0]))
        plddt = float(np.mean(sc.get("plddt", [0])))
        ptm   = float(sc.get("ptm", 0))
    
    # Compute scRMSD if AF2 structure available
    sc_rmsd = 999.0
    if af2_pdbs:
        try:
            # Load AF2 predicted structure (monomer minibinder)
            from Bio.PDB import PDBParser
            import warnings
            warnings.filterwarnings('ignore')
            parser = PDBParser(QUIET=True)
            af2_struct = parser.get_structure("af2", af2_pdbs[0])
            af2_chain = list(list(af2_struct.get_models())[0].get_chains())[0]
            af2_ca = np.array([r['CA'].get_coord() for r in af2_chain.get_residues()
                               if r.get_id()[0] == ' ' and 'CA' in r])
            
            # Load RFd3 designed backbone (minibinder portion: residues 116-165)
            cif_path = RFD3_OUT / f"{bb}.cif"
            raw = load_any(str(cif_path))
            aa = raw[0] if hasattr(raw, '__getitem__') else raw
            chain_id = list(set(aa.chain_id))[0]
            rfd3_mask = (aa.chain_id == chain_id) & \
                        np.isin(aa.atom_name, ["CA"]) & \
                        (aa.res_id >= 116) & (aa.res_id <= 165)
            rfd3_ca = aa.coord[rfd3_mask]
            
            # Align and compute RMSD (Kabsch algorithm)
            if len(af2_ca) == len(rfd3_ca) and len(af2_ca) > 0:
                # Center
                af2_centered = af2_ca - af2_ca.mean(0)
                rfd3_centered = rfd3_ca - rfd3_ca.mean(0)
                # SVD for optimal rotation
                H = af2_centered.T @ rfd3_centered
                U, S, Vt = np.linalg.svd(H)
                R = Vt.T @ U.T
                af2_rotated = af2_centered @ R.T
                sc_rmsd = float(np.sqrt(np.mean(np.sum((af2_rotated - rfd3_centered)**2, axis=1))))
        except Exception as e:
            print(f"  RMSD error for {bb}: {e}")
    
    sc_results.append({
        "backbone": bb,
        "minibinder_sequence": mb_seq,
        "sequence_recovery": design["sequence_recovery"],
        "af2_mean_plddt": round(plddt, 1),
        "af2_ptm": round(ptm, 3),
        "sc_rmsd": round(sc_rmsd, 3),
        "pass_filter": sc_rmsd < 2.0 and plddt > 60,
    })

sc_results.sort(key=lambda x: (x["sc_rmsd"] if x["sc_rmsd"] < 900 else 999))

# Save results
json.dump(sc_results, open(FINAL / "sc_results.json", "w"), indent=2)

# Write passing designs FASTA
passing = [r for r in sc_results if r["pass_filter"]]
with open(FINAL / "validated_minibinders.fasta", "w") as f:
    for r in passing:
        f.write(f">{r['backbone']}__scRMSD{r['sc_rmsd']:.2f}__pLDDT{r['af2_mean_plddt']:.0f}\n"
                f"{r['minibinder_sequence']}\n")

print(f"\n=== Self-Consistency Results ===")
print(f"{'Rank':>4}  {'Backbone':>30}  {'scRMSD':>7}  {'pLDDT':>6}  {'pTM':>5}  Pass")
print("  " + "-"*65)
for i, r in enumerate(sc_results[:20], 1):
    mark = " ✓" if r["pass_filter"] else ""
    print(f"  {i:>3}.  {r['backbone']:>30}  {r['sc_rmsd']:>7.3f}  "
          f"{r['af2_mean_plddt']:>6.1f}  {r['af2_ptm']:>5.3f}{mark}")

print(f"\nPassed filter (scRMSD < 2Å AND pLDDT > 60): {len(passing)}/{len(sc_results)}")
if passing:
    print(f"Best: {passing[0]['backbone']}  scRMSD={passing[0]['sc_rmsd']:.2f}  "
          f"pLDDT={passing[0]['af2_mean_plddt']:.0f}")
print(f"Saved: {FINAL}/sc_results.json  and  {FINAL}/validated_minibinders.fasta")
PYEOF

python3 $PIPELINE/run_selfconsistency.py

echo ""
echo "===== Pipeline complete: $(date) ====="
ls -lh $PIPELINE/outputs/final/
