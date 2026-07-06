#!/bin/bash
# ============================================================
# run_full_pipeline.sh  v3
# RFdiffusion3 (+ optional round-2 partial diffusion)
#   → ProteinMPNN → Boltz-2 (single-chain) → scRMSD
#
# Design changes:
#   - 4-chain input: VH1(A) + VL(B, steric context) + Nanobody(C)
#   - VH1 hotspots exclude VH-VL interface overlap (39,41,85 removed)
#   - Boltz-2 folds the CONNECTED single chain (not 3 separate chains)
#     → scRMSD is meaningful because topology is preserved
#   - Optional round 2: partial diffusion on round-1 templates (RFD3_ROUNDS=2)
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

# RFd3 rounds: 1 = de novo only (default), 2 = add partial-diffusion refinement
RFD3_ROUNDS=${RFD3_ROUNDS:-1}
RFD3_ROUND2_TEMPLATES=${RFD3_ROUND2_TEMPLATES:-5}
RFD3_ROUND2_DESIGNS_PER_TEMPLATE=${RFD3_ROUND2_DESIGNS_PER_TEMPLATE:-8}
RFD3_ROUND2_BATCH_SIZE=${RFD3_ROUND2_BATCH_SIZE:-8}
RFD3_PARTIAL_T=${RFD3_PARTIAL_T:-2.0}
RFD3_ROUND2_VL_CONTEXT=${RFD3_ROUND2_VL_CONTEXT:-1}
# Optional comma-separated backbone stems, e.g. mb_b003_000,mb_b006_005
RFD3_ROUND2_TEMPLATE_BACKBONES=${RFD3_ROUND2_TEMPLATE_BACKBONES:-}

if [ "$RFD3_ROUNDS" = "2" ]; then
    RFD3_ACTIVE_SUBDIR=rfd3_round2
else
    RFD3_ACTIVE_SUBDIR=rfd3
fi

PIPELINE_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "===== Full integrated pipeline v3: $(date) ====="
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "RFd3 rounds: $RFD3_ROUNDS | MB length range: $MB_LENGTH_RANGE"
if [ "$RFD3_ROUNDS" = "2" ]; then
    echo "Round 2: $RFD3_ROUND2_TEMPLATES templates × $RFD3_ROUND2_DESIGNS_PER_TEMPLATE designs"
    echo "         partial_t=${RFD3_PARTIAL_T} Å | VL context=$RFD3_ROUND2_VL_CONTEXT"
fi
mkdir -p $PIPELINE/outputs/rfd3 $PIPELINE/outputs/rfd3_round2 \
         $PIPELINE/outputs/rfd3_round2_inputs \
         $PIPELINE/outputs/mpnn \
         $PIPELINE/boltz_inputs $PIPELINE/boltz_outputs \
         $PIPELINE/final $PIPELINE/top_structures

cp "$PIPELINE_SCRIPT_DIR/rfd3_round2.py" "$PIPELINE/rfd3_round2.py"

# ── Step 0: Install Boltz-2 + pull Docker in parallel ────────────
pip install boltz[cuda] -U -q > $PIPELINE/boltz_install.log 2>&1 &
BOLTZ_PID=$!
sudo chmod 666 /var/run/docker.sock
docker pull rosettacommons/foundry:latest > $PIPELINE/docker_pull.log 2>&1 &
DOCKER_PID=$!
echo "Installing Boltz-2 (PID $BOLTZ_PID) + pulling Docker (PID $DOCKER_PID) in parallel..."
wait $DOCKER_PID && echo "Docker ready."

# ── Step 1: RFdiffusion3 ─────────────────────────────────────────
echo ""
echo "=== Step 1/5: RFdiffusion3 round 1 (4-chain context) ==="
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

# VH1 CH1-face hotspots (EXCLUDING VH-VL overlap residues 39,41,85)
# VL is chain B — in input PDB but NOT in contig → steric clash avoidance
HOTSPOTS = (
    "A11,A12,A13,A14,A15,A40,"
    "A79,A80,A81,A82,A83,A84,"
    "A103,A104,A105,A106,A107,A108,A109,A110,A111,A112,A113,A114,A115,"
    "C27,C28,C29,C30,C31,C32,C33,"
    "C52,C53,C54,C55,C56,C57,"
    "C99,C100,C101,C102,C103,C104,C105,C106,C107,C108,C109,C110,C111,C112"
)
# select_fixed_atoms: HARD PIN — VH1, VL, Nanobody frozen at input coordinates
# select_hotspots:    SOFT ATTRACTOR — minibinder drawn toward these residues
# contig length range: RFd3 samples MB length from 35-70 residues per design
FIXED_ATOMS = {
    "A1-115": "ALL",   # VH1
    "B1-115": "ALL",   # VL (steric context, frozen)
    "C1-115": "ALL",   # Nanobody
}
MB_LENGTH_RANGE = os.environ.get("MB_LENGTH_RANGE", "35-70")
CONTIG = f"A1-115,{MB_LENGTH_RANGE},C1-115"
print(f"RFd3: {N} designs | variable-length minibinder ({MB_LENGTH_RANGE} residues)")
print(f"  select_fixed_atoms: chains A,B,C pinned at input coordinates")
print(f"  select_hotspots: soft attractors on VH1-CH1-face + Nanobody CDRs")
spec = DesignInputSpecification.safe_init(
    input=INPUT,
    contig=CONTIG,
    select_hotspots=HOTSPOTS,
    select_fixed_atoms=FIXED_ATOMS,
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
PYEOF

docker run --rm --gpus all \
    -v $PIPELINE:/workspace \
    -e N_DESIGNS=$N_DESIGNS -e BATCH_SIZE=$BATCH_SIZE \
    -e MB_LENGTH_RANGE=$MB_LENGTH_RANGE \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_rfd3.py

echo "RFd3 round 1: $(ls $PIPELINE/outputs/rfd3/mb_*.cif 2>/dev/null | wc -l) designs done"

# ── Step 1b: RFdiffusion3 round 2 (optional partial diffusion) ──
if [ "$RFD3_ROUNDS" = "2" ]; then
    echo ""
    echo "=== Step 1b/5: RFdiffusion3 round 2 (partial diffusion) ==="
    docker run --rm --gpus all \
        -v $PIPELINE:/workspace \
        -e RFD3_ROUND2_TEMPLATES=$RFD3_ROUND2_TEMPLATES \
        -e RFD3_ROUND2_DESIGNS_PER_TEMPLATE=$RFD3_ROUND2_DESIGNS_PER_TEMPLATE \
        -e RFD3_ROUND2_BATCH_SIZE=$RFD3_ROUND2_BATCH_SIZE \
        -e RFD3_PARTIAL_T=$RFD3_PARTIAL_T \
        -e RFD3_ROUND2_VL_CONTEXT=$RFD3_ROUND2_VL_CONTEXT \
        -e RFD3_ROUND2_TEMPLATE_BACKBONES="$RFD3_ROUND2_TEMPLATE_BACKBONES" \
        -e FOUNDRY_CHECKPOINT_DIRS=/weights \
        rosettacommons/foundry:latest \
        python3 /workspace/rfd3_round2.py
    echo "RFd3 round 2: $(ls $PIPELINE/outputs/rfd3_round2/r2_*.cif 2>/dev/null | wc -l) designs done"
fi

# ── Step 2: ProteinMPNN ──────────────────────────────────────────
echo ""
echo "=== Step 2/5: ProteinMPNN ==="
cat > $PIPELINE/run_mpnn.py << 'PYEOF'
import os, json
from pathlib import Path
from atomworks.io.utils.io_utils import load_any
from mpnn.inference_engines.mpnn import MPNNInferenceEngine

RFD3_ROUNDS = os.environ.get("RFD3_ROUNDS", "1")
RFD3_SUBDIR = os.environ.get("RFD3_ACTIVE_SUBDIR", "rfd3_round2" if RFD3_ROUNDS == "2" else "rfd3")
RFD3  = Path("/workspace/outputs") / RFD3_SUBDIR
OUT   = Path("/workspace/outputs/mpnn")
NSEQS = int(os.environ.get("N_MPNN_SEQS", 8))
OUT.mkdir(parents=True, exist_ok=True)

engine = MPNNInferenceEngine(model_type="protein_mpnn", is_legacy_weights=True,
    out_directory=None, write_structures=False, write_fasta=False)

results = []
glob_pat = "r2_*.cif" if RFD3_ROUNDS == "2" else "mb_*.cif"
cifs = sorted(RFD3.glob(glob_pat))
print(f"MPNN ({RFD3_SUBDIR}): {len(cifs)} × {NSEQS}")

for idx, cif in enumerate(cifs):
    raw = load_any(str(cif)); aa = raw[0] if hasattr(raw,"__getitem__") else raw
    ch = list(set(aa.chain_id))[0]
    all_res = sorted(set(aa.res_id))
    total = len(all_res)
    # Detect variable minibinder length: total = 115 (VH1) + MB + 115 (Nb)
    mb_len = total - 230
    mb_start, mb_end = 116, 115 + mb_len
    nb_start = mb_end + 1
    designed = [f"{ch}{r}" for r in range(mb_start, mb_end+1) if r in set(aa.res_id)]
    result = engine.run(atom_arrays=[aa], input_dicts=[{
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
print(f"MPNN done: {len(results)} sequences")
PYEOF

docker run --rm --gpus all \
    -v $PIPELINE:/workspace \
    -e N_MPNN_SEQS=$N_MPNN_SEQS \
    -e RFD3_ROUNDS=$RFD3_ROUNDS \
    -e RFD3_ACTIVE_SUBDIR=$RFD3_ACTIVE_SUBDIR \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_mpnn.py

# ── Step 3: Boltz-2 SINGLE-CHAIN validation ──────────────────────
echo ""
echo "=== Step 3/5: Boltz-2 (single connected chain) ==="

wait $BOLTZ_PID 2>/dev/null || pip install boltz[cuda] -U -q
pip install -q 'networkx>=3.0' 'platformdirs>=3.0' 2>/dev/null || true

cat > $PIPELINE/build_boltz_inputs.py << 'PYEOF'
"""
Build Boltz-2 inputs for SINGLE-CHAIN validation.

For each top-N design, build ONE sequence:
  VH1(115) + Minibinder(55) + Nanobody(115) = 285 residues
  (same as the RFd3 CIF single chain — no explicit linkers)

Boltz-2 predicts the folded structure of this single chain.
scRMSD: align Boltz-2 on VH1(1-115)+Nb(171-285), measure RMSD of MB(116-170).

This is the correct self-consistency check for a CONNECTED bispecific binder.
"""
import json
from pathlib import Path
from atomworks.io.utils.io_utils import load_any
import biotite.sequence as bseq

PIPELINE = Path("/workspace")
MPNN_OUT = PIPELINE / "outputs/mpnn"
RFD3_SUBDIR = __import__('os').environ.get("RFD3_ACTIVE_SUBDIR", "rfd3")
RFD3_OUT = PIPELINE / "outputs" / RFD3_SUBDIR
BOLTZ_IN = PIPELINE / "boltz_inputs"
BOLTZ_IN.mkdir(exist_ok=True)

mpnn_data = json.load(open(MPNN_OUT/"all_sequences.json"))
mpnn_data.sort(key=lambda x: x["sequence_recovery"])
seen, top_designs = set(), []
for r in mpnn_data:
    if r["backbone"] not in seen:
        seen.add(r["backbone"]); top_designs.append(r)
    if len(top_designs) >= int(__import__('os').environ.get("TOP_N", 50)): break

print(f"Building {len(top_designs)} single-chain Boltz-2 inputs...")
design_list = []
for i, r in enumerate(top_designs):
    bb = r["backbone"]
    # Extract full 285-residue sequence from the RFd3 MPNN output
    full_seq = r["full_sequence"]  # VH1(1-115) + Minibinder(116-170) + Nanobody(171-285)
    if len(full_seq) < 171:
        print(f"  Skip {bb}: short sequence ({len(full_seq)})"); continue

    rank = i + 1
    rec  = r["sequence_recovery"]
    name = f"rank{rank:02d}_{bb}_rec{rec:.3f}"

    # Single-chain YAML — Boltz-2 folds it as one polypeptide
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

json.dump({"designs":design_list,
           "validation_note":"Single connected chain: VH1(1-115)+MB(116-170)+Nb(171-285). "
                             "scRMSD aligns on VH1+Nb portions, measures MB RMSD."},
          open(PIPELINE/"top_designs.json","w"), indent=2)
print(f"Wrote {len(design_list)} YAML files")
PYEOF

docker run --rm --gpus all \
    -v $PIPELINE:/workspace \
    -e TOP_N=$TOP_N \
    -e RFD3_ACTIVE_SUBDIR=$RFD3_ACTIVE_SUBDIR \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/build_boltz_inputs.py

N_YAML=$(ls $PIPELINE/boltz_inputs/*.yaml | wc -l)
echo "Running Boltz-2 on $N_YAML single-chain sequences..."

$HOME/.local/bin/boltz predict $PIPELINE/boltz_inputs \
    --out_dir $PIPELINE/boltz_outputs \
    --devices 1 --num_workers 2 --override 2>&1 | \
    grep -E "Predicting|Saving|Done|Error" | head -20

# ── Step 4: scRMSD scoring ────────────────────────────────────────
echo ""
echo "=== Step 4/5: scRMSD + ipTM scoring ==="

cat > $PIPELINE/final_score.py << 'PYEOF'
"""
scRMSD: compare Boltz-2 single-chain prediction to RFd3 backbone.
Alignment: VH1(A,1-115) + Nanobody(A,171-285) as rigid anchors.
Metric:    RMSD of Minibinder(A,116-170) after alignment.
"""
import json, glob, shutil
from pathlib import Path
import numpy as np

PIPELINE   = Path("/workspace")
BOLTZ_BASE = PIPELINE / "boltz_outputs/boltz_results_boltz_inputs/predictions"
RFD3_OUT   = PIPELINE / "outputs/rfd3"
MPNN_OUT   = PIPELINE / "outputs/mpnn"
FINAL      = PIPELINE / "final"
STRUCTS    = PIPELINE / "top_structures"
for d in [FINAL, STRUCTS]: d.mkdir(exist_ok=True)

from atomworks.io.utils.io_utils import load_any
import biotite.structure.io.pdbx as pdbx

def kabsch_rmsd(P, Q):
    P=P-P.mean(0); Q=Q-Q.mean(0)
    U,S,Vt=np.linalg.svd(P.T@Q); d=np.linalg.det(Vt.T@U.T)
    R=Vt.T@np.diag([1,1,d])@U.T
    return float(np.sqrt(np.mean(np.sum((P@R.T-Q)**2,axis=1))))

def rfd3_ca(cif,r1,r2):
    raw=load_any(str(cif)); aa=raw[0] if hasattr(raw,"__getitem__") else raw
    ch=list(set(aa.chain_id))[0]
    m=(aa.chain_id==ch)&np.isin(aa.atom_name,["CA"])&(aa.res_id>=r1)&(aa.res_id<=r2)
    return aa.coord[m]

def boltz_single_ca(cif, r1, r2):
    """Boltz-2 single-chain prediction — all residues on chain A."""
    aa=pdbx.get_structure(pdbx.CIFFile.read(str(cif)),model=1,include_bonds=False)
    m=(aa.atom_name=="CA")&(aa.chain_id=="A")&(aa.res_id>=r1)&(aa.res_id<=r2)
    return aa.coord[m]

top_data = json.load(open(PIPELINE/"top_designs.json"))
designs  = top_data["designs"]
pred_dirs = {d.name:d for d in BOLTZ_BASE.iterdir() if d.is_dir()}

results = []
for r in designs:
    rank,bb,mb,rec = r["rank"],r["backbone"],r["minibinder_sequence"],r["sequence_recovery"]
    name = f"rank{rank:02d}_{bb}_rec{rec:.3f}"
    pred_dir = pred_dirs.get(name)

    plddt = iptm = 0.0
    sc_rmsd = 999.0

    if pred_dir:
        conf = pred_dir / f"confidence_{name}_model_0.json"
        cif  = pred_dir / f"{name}_model_0.cif"
        if conf.exists():
            c = json.load(open(conf))
            plddt = float(c.get("complex_plddt", 0))
            iptm  = float(c.get("iptm", 0))  # single chain — use overall pTM/pLDDT

        if cif.exists():
            try:
                rc = RFD3_OUT / f"{bb}.cif"
                if rc.exists():
                    # RFd3 backbone: single chain, VH1=1-115, MB=116-170, Nb=171-285
                    rv = rfd3_ca(rc, 1,   115)
                    rm = rfd3_ca(rc, 116, 170)
                    rn = rfd3_ca(rc, 171, 285)
                    # Boltz-2: single chain A, same residue numbering
                    bv = boltz_single_ca(cif, 1,   115)
                    bm = boltz_single_ca(cif, 116, 170)
                    bn = boltz_single_ca(cif, 171, 285)
                    nv=min(len(bv),len(rv)); nn=min(len(bn),len(rn)); nm=min(len(bm),len(rm))
                    if nv>10 and nn>10 and nm>5:
                        P=np.vstack([bv[:nv],bn[:nn]]); Q=np.vstack([rv[:nv],rn[:nn]])
                        Pc=P-P.mean(0); Qc=Q-Q.mean(0)
                        U,S,Vt=np.linalg.svd(Pc.T@Qc); d=np.linalg.det(Vt.T@U.T)
                        Rmat=Vt.T@np.diag([1,1,d])@U.T
                        sc_rmsd=kabsch_rmsd((bm[:nm]-P.mean(0))@Rmat.T+Q.mean(0),rm[:nm])
            except Exception as e: print(f"  scRMSD err {bb}: {e}")

    # For single-chain: pLDDT>0.60 AND scRMSD<2Å is a passing design
    passed = sc_rmsd < 2.0 and plddt > 0.60
    results.append({"rank":rank,"backbone":bb,"minibinder_sequence":mb,
                    "sequence_recovery":rec,"boltz_plddt":round(plddt,3),
                    "boltz_plddt_pct":round(plddt*100,1),
                    "boltz_iptm":round(iptm,3),
                    "sc_rmsd_A":round(sc_rmsd,3),"pass_filter":passed})

results.sort(key=lambda x:(x["sc_rmsd_A"] if x["sc_rmsd_A"]<900 else 999, -x["boltz_plddt"]))
json.dump(results, open(FINAL/"final_results.json","w"), indent=2)
passing=[r for r in results if r["pass_filter"]]
with open(FINAL/"validated_minibinders.fasta","w") as f:
    for r in passing:
        f.write(f">{r['backbone']}__scRMSD{r['sc_rmsd_A']:.2f}__pLDDT{r['boltz_plddt_pct']:.0f}\n{r['minibinder_sequence']}\n")

# Top-10 by pLDDT x low_scRMSD
scored=[(r,r["boltz_plddt"]/(r["sc_rmsd_A"] if r["sc_rmsd_A"]<900 else 999+1)) for r in results]
scored.sort(key=lambda x:-x[1])
top10=[r for r,_ in scored[:10]]

for i,r in enumerate(top10):
    bb=r["backbone"]; rec=r["sequence_recovery"]; rank=r["rank"]
    name=f"rank{rank:02d}_{bb}_rec{rec:.3f}"
    pred_dir=pred_dirs.get(name)
    rfd3_src=RFD3_OUT/f"{bb}.cif"
    if rfd3_src.exists(): shutil.copy(rfd3_src, STRUCTS/f"top{i+1:02d}_{bb}_rfd3.cif")
    if pred_dir:
        bs=pred_dir/f"{name}_model_0.cif"
        if bs.exists(): shutil.copy(bs, STRUCTS/f"top{i+1:02d}_{bb}_boltz2.cif")

with open(STRUCTS/"top10_summary.tsv","w") as f:
    f.write("rank\tbackbone\tsc_rmsd_A\tboltz_plddt_pct\tboltz_iptm\tsequence_recovery\tminibinder_sequence\n")
    for i,r in enumerate(top10):
        f.write(f"{i+1}\t{r['backbone']}\t{r['sc_rmsd_A']}\t{r['boltz_plddt_pct']}\t"
                f"{r['boltz_iptm']}\t{r['sequence_recovery']}\t{r['minibinder_sequence']}\n")

print(f"\n{'='*72}")
print(f"{'R':>3}  {'Backbone':>28}  {'scRMSD':>7}  {'pLDDT%':>6}  {'pTM':>5}  Pass")
print("  "+"-"*70)
for r in results[:20]:
    print(f"{r['rank']:>3}.  {r['backbone']:>28}  {r['sc_rmsd_A']:>7.3f}  "
          f"{r['boltz_plddt_pct']:>6.1f}  {r['boltz_iptm']:>5.3f}{'  ✓' if r['pass_filter'] else ''}")
print(f"\nPassed (scRMSD<2Å & pLDDT>60%): {len(passing)}/{len(results)}")
if passing:
    p=passing[0]
    print(f"\nBest: {p['backbone']}")
    print(f"  Sequence: {p['minibinder_sequence']}")
    print(f"  scRMSD={p['sc_rmsd_A']:.2f}Å  pLDDT={p['boltz_plddt_pct']:.0f}%")
print(f"Saved {len(list(STRUCTS.glob('*.cif')))} CIF files")
PYEOF

docker run --rm --gpus all \
    -v $PIPELINE:/workspace \
    -e RFD3_ACTIVE_SUBDIR=$RFD3_ACTIVE_SUBDIR \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/binary_antibodies/scoring.py \
        --pipeline-dir /workspace \
        --mode single_chain \
        --plddt-threshold 0.60 \
        --scrmsd-threshold 2.0 \
        --n-top-cifs 10

# ── Step 5: Push to GitHub ────────────────────────────────────────
echo ""
echo "=== Step 5/5: Push to GitHub ==="
bash $PIPELINE/push_results.sh
echo ""
echo "===== Complete: $(date) ====="
