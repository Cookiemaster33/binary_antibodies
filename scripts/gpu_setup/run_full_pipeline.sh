#!/bin/bash
# ============================================================
# run_full_pipeline.sh
# --------------------
# Fully integrated pipeline on ONE instance:
#   1. RFdiffusion3 — design bispecific bridging minibinder
#   2. ProteinMPNN  — design sequences
#   3. Boltz-2      — predict full 3-chain complex for top-50
#   4. scRMSD       — compare Boltz-2 minibinder to RFd3 backbone
#                     (possible because CIFs are on the same instance)
#
# Design: VH1(A, 115 res) + Minibinder(55 res) + Nanobody(B, 115 res)
# Contig: A1-115,55,B1-115
# ============================================================
set -eo pipefail
PIPELINE=/home/ubuntu/pipeline
LOG=$PIPELINE/full_pipeline.log
exec > >(tee -a "$LOG") 2>&1

N_DESIGNS=${N_DESIGNS:-200}
BATCH_SIZE=${BATCH_SIZE:-10}
N_MPNN_SEQS=${N_MPNN_SEQS:-8}
N_BATCHES=$(( (N_DESIGNS + BATCH_SIZE - 1) / BATCH_SIZE ))
TOP_N=${TOP_N:-50}   # how many designs to run through Boltz-2

echo "===== Full integrated pipeline: $(date) ====="
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "N_DESIGNS=$N_DESIGNS  N_MPNN_SEQS=$N_MPNN_SEQS  TOP_N=$TOP_N"
mkdir -p $PIPELINE/outputs/rfd3 $PIPELINE/outputs/mpnn \
         $PIPELINE/boltz_inputs $PIPELINE/boltz_outputs $PIPELINE/final

# ── 1. Install Boltz-2 ───────────────────────────────────────────
echo ""
echo "=== Step 0: Install Boltz-2 (background) ==="
pip install boltz[cuda] -U -q > $PIPELINE/boltz_install.log 2>&1 &
BOLTZ_INSTALL_PID=$!
echo "  Boltz-2 installing in background (PID $BOLTZ_INSTALL_PID)..."

# ── 2. Pull foundry Docker (parallel with Boltz install) ─────────
echo "=== Step 0b: Pull foundry Docker (parallel) ==="
sudo chmod 666 /var/run/docker.sock
docker pull rosettacommons/foundry:latest > $PIPELINE/docker_pull.log 2>&1 &
DOCKER_PID=$!
echo "  Docker pull in background (PID $DOCKER_PID)..."

# Wait for Docker (needed first for RFd3)
echo "  Waiting for Docker pull..."
wait $DOCKER_PID && echo "  Docker ready." || echo "  Docker pull error."

# ── 3. RFdiffusion3 ──────────────────────────────────────────────
echo ""
echo "=== Step 1/4: RFdiffusion3 (bispecific bridging minibinder) ==="

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

HOTSPOTS = (
    "A11,A12,A13,A14,A15,A39,A40,A41,"
    "A79,A80,A81,A82,A83,A84,A85,"
    "A103,A104,A105,A106,A107,A108,A109,A110,A111,A112,A113,A114,A115,"
    "B27,B28,B29,B30,B31,B32,B33,"
    "B52,B53,B54,B55,B56,B57,"
    "B99,B100,B101,B102,B103,B104,B105,B106,B107,B108,B109,B110,B111,B112"
)

print(f"RFd3: {N} designs  contig=A1-115,55,B1-115")
spec   = DesignInputSpecification.safe_init(
    input=INPUT, contig="A1-115,55,B1-115", select_hotspots=HOTSPOTS)
model  = RFD3InferenceEngine(**RFD3InferenceConfig(diffusion_batch_size=BATCH))

saved = 0
for bi in range(NBATCH):
    print(f"  Batch {bi+1}/{NBATCH} ...", flush=True)
    for key, designs in model.run(inputs=spec, out_dir=None, n_batches=1).items():
        for i, d in enumerate(designs):
            to_cif_file(d.atom_array, str(OUT / f"mb_b{bi:03d}_{i:03d}.cif"))
            saved += 1
print(f"RFd3 done: {saved} designs  (A=VH1 1-115, MB=116-170, B=Nanobody 171-285)")
PYEOF

docker run --rm --gpus all \
    -v $PIPELINE:/workspace \
    -e N_DESIGNS=$N_DESIGNS -e BATCH_SIZE=$BATCH_SIZE \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_rfd3.py

echo "RFd3 done: $(ls $PIPELINE/outputs/rfd3/mb_*.cif | wc -l) designs"

# ── 4. ProteinMPNN ───────────────────────────────────────────────
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
print(f"MPNN: {len(cifs)} × {NSEQS}  (designing residues 116-170)")

for idx, cif in enumerate(cifs):
    raw = load_any(str(cif)); aa = raw[0] if hasattr(raw,"__getitem__") else raw
    ch = list(set(aa.chain_id))[0]
    designed = [f"{ch}{r}" for r in range(116,171) if r in set(aa.res_id)]
    result = engine.run(atom_arrays=[aa], input_dicts=[{
        "batch_size": NSEQS, "remove_waters": True, "designed_residues": designed}])
    for r in (result if isinstance(result,list) else [result]):
        seq = r.output_dict.get("designed_sequence","")
        rec = float(r.output_dict.get("sequence_recovery",0))
        results.append({"backbone": cif.stem, "full_sequence": seq,
                        "minibinder_sequence": seq[115:171] if len(seq)>=171 else seq[115:],
                        "sequence_recovery": rec})
    if (idx+1)%20==0: print(f"  {idx+1}/{len(cifs)} done...", flush=True)

results.sort(key=lambda x: x["sequence_recovery"])
json.dump(results, open(OUT/"all_sequences.json","w"), indent=2)
fasta = OUT/"minibinder_sequences.fasta"
with open(fasta,"w") as f:
    for i,r in enumerate(results):
        f.write(f">rank_{i+1}_{r['backbone']}_rec{r['sequence_recovery']:.3f}\n{r['minibinder_sequence']}\n")
print(f"MPNN done: {len(results)} sequences")
if results: print(f"Top: {results[0]['minibinder_sequence'][:50]}  rec={results[0]['sequence_recovery']:.3f}")
PYEOF

docker run --rm --gpus all \
    -v $PIPELINE:/workspace \
    -e N_MPNN_SEQS=$N_MPNN_SEQS \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_mpnn.py

# ── 5. Wait for Boltz-2 install to finish ────────────────────────
echo ""
echo "=== Waiting for Boltz-2 install to complete ==="
wait $BOLTZ_INSTALL_PID 2>/dev/null || true
if ! $HOME/.local/bin/boltz --help &>/dev/null; then
    echo "  Re-installing Boltz-2..."
    pip install boltz[cuda] -U -q
fi
echo "  Boltz-2 ready: $(python3 -c 'import boltz; print(boltz.__version__)')"

# Also download Boltz CCD data now (first-run download)
echo "  Pre-downloading Boltz CCD data..."
python3 -c "from boltz.main import download; download('$HOME/.boltz')" 2>/dev/null || true

# ── 6. Build Boltz-2 inputs + run ────────────────────────────────
echo ""
echo "=== Step 3/4: Boltz-2 complex prediction (top $TOP_N designs) ==="

python3 << PYEOF
import json
from pathlib import Path

mpnn_data = json.load(open("$PIPELINE/outputs/mpnn/all_sequences.json"))
validation = json.load(open("$PIPELINE/validation_inputs.json"))
VH1_SEQ = validation["VH1_seq"]
NB_SEQ  = validation["NB_seq"]

# Top-N unique backbones
seen, top_designs = set(), []
for r in mpnn_data:
    if r["backbone"] not in seen:
        seen.add(r["backbone"]); top_designs.append(r)
    if len(top_designs) >= $TOP_N: break

VH1_CONTACTS = [11,12,13,14,15,39,40,41,79,80,81,82,83,84,85,103,104,105,106,107,108,109,110,111,112,113,114,115]
NB_CDR       = list(range(27,34)) + list(range(52,58)) + list(range(99,113))

outdir = Path("$PIPELINE/boltz_inputs")
outdir.mkdir(exist_ok=True)

# Save design list for scRMSD step
json.dump({"designs": [{"rank":i+1,"backbone":r["backbone"],
                         "minibinder_sequence":r["minibinder_sequence"],
                         "sequence_recovery":r["sequence_recovery"]}
                        for i,r in enumerate(top_designs)],
            "VH1_seq": VH1_SEQ, "NB_seq": NB_SEQ},
          open("$PIPELINE/top_designs.json","w"), indent=2)

contacts_str = (
    "[" +
    ", ".join(f'["A", {res}]' for res in VH1_CONTACTS) + ", " +
    ", ".join(f'["C", {res}]' for res in NB_CDR) +
    "]"
)
for i, r in enumerate(top_designs):
    name = f"rank{i+1:02d}_{r['backbone']}_rec{r['sequence_recovery']:.3f}"
    yaml = f"""sequences:
  - protein:
      id: A
      sequence: "{VH1_SEQ}"
      msa: empty
  - protein:
      id: B
      sequence: "{r['minibinder_sequence']}"
      msa: empty
  - protein:
      id: C
      sequence: "{NB_SEQ}"
      msa: empty
constraints:
  - pocket:
      binder: B
      contacts: {contacts_str}
      max_distance: 10.0
"""
    (outdir / f"{name}.yaml").write_text(yaml)

print(f"Wrote {len(top_designs)} YAML files for Boltz-2")
PYEOF

$HOME/.local/bin/boltz predict $PIPELINE/boltz_inputs \
    --out_dir $PIPELINE/boltz_outputs \
    --devices 1 \
    --num_workers 2 \
    --override 2>&1 | tee $PIPELINE/boltz_run.log

echo "Boltz-2 predictions done."

# ── 7. scRMSD + ipTM scoring ─────────────────────────────────────
echo ""
echo "=== Step 4/4: scRMSD (Boltz-2 vs RFd3 backbone) + ipTM ==="

python3 << 'PYEOF'
import json, glob
from pathlib import Path
import numpy as np

PIPELINE   = Path("/home/ubuntu/pipeline")
RFD3_OUT   = PIPELINE / "outputs/rfd3"
BOLTZ_OUT  = PIPELINE / "boltz_outputs"
FINAL      = PIPELINE / "final"
FINAL.mkdir(exist_ok=True)

top_data = json.load(open(PIPELINE / "top_designs.json"))
designs  = top_data["designs"]

# ── Helpers ──────────────────────────────────────────────────────
def kabsch_rmsd(P, Q):
    P = P - P.mean(0); Q = Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(P.T @ Q)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1,1,d]) @ U.T
    return float(np.sqrt(np.mean(np.sum((P @ R.T - Q)**2, axis=1))))

def load_cif_ca(cif_path, chain_id=None, res_min=None, res_max=None):
    """Load Cα coords from a CIF file using biotite."""
    import biotite.structure.io.pdbx as pdbx
    cif = pdbx.CIFFile.read(str(cif_path))
    aa  = pdbx.get_structure(cif, model=1, include_bonds=False)
    mask = aa.atom_name == "CA"
    if chain_id:  mask &= (aa.chain_id == chain_id)
    if res_min:   mask &= (aa.res_id >= res_min)
    if res_max:   mask &= (aa.res_id <= res_max)
    return aa.coord[mask], aa.res_id[mask]

def load_rfd3_ca(cif_path, res_min, res_max):
    """Load Cα from RFd3 single-chain CIF (all chain A)."""
    from atomworks.io.utils.io_utils import load_any
    raw = load_any(str(cif_path))
    aa  = raw[0] if hasattr(raw, "__getitem__") else raw
    ch  = list(set(aa.chain_id))[0]
    mask = (aa.chain_id == ch) & np.isin(aa.atom_name, ["CA"]) & \
           (aa.res_id >= res_min) & (aa.res_id <= res_max)
    return aa.coord[mask]

results = []
for r in designs:
    rank, bb, mb = r["rank"], r["backbone"], r["minibinder_sequence"]
    rec  = r["sequence_recovery"]
    name = f"rank{rank:02d}_{bb}_rec{rec:.3f}"

    # Boltz-2 output CIF
    boltz_cifs = sorted(glob.glob(
        str(BOLTZ_OUT / f"**/{name}/**/*.cif"), recursive=True))
    # Confidence JSON
    conf_files = sorted(glob.glob(
        str(BOLTZ_OUT / f"**/{name}*/**/*confidence*.json"), recursive=True))

    plddt = iptm_ab = iptm_bc = 0.0
    if conf_files:
        try:
            c = json.load(open(conf_files[0]))
            plddt   = float(c.get("mean_plddt", 0) or 0)
            pair    = c.get("pair_chains_iptm", c.get("chain_pair_iptm", {}))
            if isinstance(pair, dict):
                iptm_ab = float(pair.get("A-B", pair.get("AB", 0)) or 0)
                iptm_bc = float(pair.get("B-C", pair.get("BC", 0)) or 0)
            if not iptm_ab:
                iptm_ab = iptm_bc = float(c.get("iptm", 0) or 0)
        except Exception as e:
            print(f"  Conf parse error {name}: {e}")

    # scRMSD: superimpose Boltz VH1(A)+Nb(C) onto RFd3, then measure MB RMSD
    sc_rmsd = 999.0
    if boltz_cifs:
        try:
            rfd3_cif = RFD3_OUT / f"{bb}.cif"
            if rfd3_cif.exists():
                # RFd3 backbone (single chain A): VH1=1-115, MB=116-170, Nb=171-285
                rfd3_vh1 = load_rfd3_ca(rfd3_cif, 1,   115)
                rfd3_mb  = load_rfd3_ca(rfd3_cif, 116, 170)
                rfd3_nb  = load_rfd3_ca(rfd3_cif, 171, 285)

                # Boltz-2 output (3 chains: A=VH1, B=MB, C=Nb)
                boltz_vh1, _ = load_cif_ca(boltz_cifs[0], chain_id="A")
                boltz_mb,  _ = load_cif_ca(boltz_cifs[0], chain_id="B")
                boltz_nb,  _ = load_cif_ca(boltz_cifs[0], chain_id="C")

                # Build alignment anchor: VH1 + Nanobody (the fixed parts)
                n_vh1 = min(len(boltz_vh1), len(rfd3_vh1))
                n_nb  = min(len(boltz_nb),  len(rfd3_nb))
                n_mb  = min(len(boltz_mb),  len(rfd3_mb))

                if n_vh1 > 10 and n_nb > 10 and n_mb > 5:
                    P_anch = np.vstack([boltz_vh1[:n_vh1], boltz_nb[:n_nb]])
                    Q_anch = np.vstack([rfd3_vh1[:n_vh1],  rfd3_nb[:n_nb]])

                    P_c = P_anch - P_anch.mean(0)
                    Q_c = Q_anch - Q_anch.mean(0)
                    U, S, Vt = np.linalg.svd(P_c.T @ Q_c)
                    d = np.linalg.det(Vt.T @ U.T)
                    Rmat = Vt.T @ np.diag([1,1,d]) @ U.T

                    # Apply rotation+translation to Boltz MB
                    boltz_mb_aligned = (boltz_mb[:n_mb] - P_anch.mean(0)) @ Rmat.T + Q_anch.mean(0)
                    sc_rmsd = kabsch_rmsd(boltz_mb_aligned, rfd3_mb[:n_mb])
                else:
                    print(f"  Length issue {name}: VH1={n_vh1} Nb={n_nb} MB={n_mb}")
        except Exception as e:
            print(f"  scRMSD error {name}: {e}")

    passed = sc_rmsd < 2.0 and plddt > 60 and iptm_ab > 0.3 and iptm_bc > 0.3
    results.append({
        "rank": rank, "backbone": bb,
        "minibinder_sequence": mb,
        "sequence_recovery": rec,
        "boltz_mean_plddt": round(plddt, 1),
        "iptm_MB_VH1": round(iptm_ab, 3),
        "iptm_MB_Nb":  round(iptm_bc, 3),
        "sc_rmsd_A": round(sc_rmsd, 3),
        "pass_filter": passed,
    })
    print(f"  rank{rank:02d} {bb}  scRMSD={sc_rmsd:.2f}Å  "
          f"pLDDT={plddt:.1f}  ipTM↔VH1={iptm_ab:.3f}  ipTM↔Nb={iptm_bc:.3f}"
          f"{'  ✓' if passed else ''}")

# Sort by scRMSD (ascending), then ipTMs
results.sort(key=lambda x: (x["sc_rmsd_A"] if x["sc_rmsd_A"]<900 else 999,
                             -(x["iptm_MB_VH1"]*x["iptm_MB_Nb"])**0.5))

json.dump(results, open(FINAL/"final_results.json","w"), indent=2)
passing = [r for r in results if r["pass_filter"]]

with open(FINAL/"validated_minibinders.fasta","w") as f:
    for r in passing:
        f.write(f">{r['backbone']}"
                f"__scRMSD{r['sc_rmsd_A']:.2f}"
                f"__pLDDT{r['boltz_mean_plddt']:.0f}"
                f"__ipTM_VH1_{r['iptm_MB_VH1']:.2f}"
                f"__ipTM_Nb_{r['iptm_MB_Nb']:.2f}\n"
                f"{r['minibinder_sequence']}\n")

print(f"\n{'='*70}")
print(f"{'R':>3}  {'Backbone':>28}  {'scRMSD':>7}  {'pLDDT':>6}  {'↔VH1':>5}  {'↔Nb':>5}  Pass")
print("  " + "-"*68)
for r in results[:20]:
    print(f"{r['rank']:>3}.  {r['backbone']:>28}  {r['sc_rmsd_A']:>7.3f}  "
          f"{r['boltz_mean_plddt']:>6.1f}  {r['iptm_MB_VH1']:>5.3f}  "
          f"{r['iptm_MB_Nb']:>5.3f}{'  ✓' if r['pass_filter'] else ''}")

print(f"\nPassed filter (scRMSD<2Å, pLDDT>60, both ipTM>0.3): {len(passing)}/{len(results)}")
if passing:
    p = passing[0]
    print(f"\nBest validated design:")
    print(f"  Backbone:   {p['backbone']}")
    print(f"  Sequence:   {p['minibinder_sequence']}")
    print(f"  scRMSD:     {p['sc_rmsd_A']:.2f} Å  (< 2.0 Å = correct position ✓)")
    print(f"  pLDDT:      {p['boltz_mean_plddt']:.0f}")
    print(f"  ipTM↔VH1:   {p['iptm_MB_VH1']:.3f}  (> 0.3 = confident interface ✓)")
    print(f"  ipTM↔Nb:    {p['iptm_MB_Nb']:.3f}  (> 0.3 = confident interface ✓)")
PYEOF

echo ""
echo "===== Full pipeline complete: $(date) ====="
echo "Results in: $PIPELINE/final/"
ls -lh $PIPELINE/final/

# ── 8. Collect top-10 structure files ────────────────────────────
echo ""
echo "=== Collecting top-10 CIF structures ==="

python3 << 'PYEOF'
import json, shutil, os
from pathlib import Path
import numpy as np

HOME      = Path.home()
PIPELINE  = HOME / "pipeline"
BOLTZ_BASE = PIPELINE / "boltz_outputs/boltz_results_boltz_inputs/predictions"
RFD3_OUT  = PIPELINE / "outputs/rfd3"
MPNN_OUT  = PIPELINE / "outputs/mpnn"
STRUCTS   = PIPELINE / "top_structures"
STRUCTS.mkdir(exist_ok=True)

# Load final results
results = json.load(open(PIPELINE / "final/final_results.json"))

# Rank by geometric mean of ipTMs (both interfaces must be good)
results_with_score = []
for r in results:
    if r["iptm_MB_VH1"] > 0 and r["iptm_MB_Nb"] > 0:
        combined = (r["iptm_MB_VH1"] * r["iptm_MB_Nb"]) ** 0.5
    else:
        combined = 0
    results_with_score.append((combined, r))

results_with_score.sort(reverse=True, key=lambda x: x[0])
top10 = [r for _, r in results_with_score[:10]]

print(f"Top 10 designs by √(ipTM_VH1 × ipTM_Nb):")
for i, r in enumerate(top10):
    bb = r["backbone"]
    name = f"rank{r['rank']:02d}_{bb}_rec{r['sequence_recovery']:.3f}"
    combined = (r["iptm_MB_VH1"] * r["iptm_MB_Nb"]) ** 0.5

    # Copy RFd3 CIF (the designed backbone)
    rfd3_src = RFD3_OUT / f"{bb}.cif"
    if rfd3_src.exists():
        shutil.copy(rfd3_src, STRUCTS / f"top{i+1:02d}_{bb}_rfd3_backbone.cif")
        print(f"  top{i+1:02d} {bb}  √ipTM={combined:.3f}  pLDDT={r['boltz_plddt_pct']:.0f}%  ✓ RFd3 CIF")
    else:
        print(f"  top{i+1:02d} {bb}  ✗ RFd3 CIF not found")

    # Copy Boltz-2 CIF (the predicted complex)
    boltz_cif = BOLTZ_BASE / name / f"{name}_model_0.cif"
    if boltz_cif.exists():
        shutil.copy(boltz_cif, STRUCTS / f"top{i+1:02d}_{bb}_boltz2_complex.cif")
        print(f"           ✓ Boltz-2 complex CIF")
    else:
        print(f"           ✗ Boltz-2 CIF not found ({boltz_cif})")

# Save a summary TSV for easy reading
with open(STRUCTS / "top10_summary.tsv", "w") as f:
    f.write("rank\tbackbone\tcombined_iptm\tiptm_VH1\tiptm_Nb\tpLDDT_pct\tsc_rmsd_A\tminibinder_sequence\n")
    for i, r in enumerate(top10):
        combined = (r["iptm_MB_VH1"] * r["iptm_MB_Nb"]) ** 0.5
        f.write(f"{i+1}\t{r['backbone']}\t{combined:.3f}\t{r['iptm_MB_VH1']:.3f}\t"
                f"{r['iptm_MB_Nb']:.3f}\t{r['boltz_plddt_pct']:.1f}\t{r['sc_rmsd_A']:.2f}\t"
                f"{r['minibinder_sequence']}\n")

print(f"\nSaved {len(list(STRUCTS.glob('*.cif')))} CIF files + summary TSV to {STRUCTS}")
PYEOF

# Push structures to GitHub
echo ""
echo "=== Pushing structures to GitHub ==="
STRUCTS_DIR=$HOME/pipeline/top_structures
REPO_STRUCTS=$HOME/repo/pipeline_results/v3_bispecific_validated/structures
mkdir -p $REPO_STRUCTS

cp -f $STRUCTS_DIR/*.cif $REPO_STRUCTS/ 2>/dev/null || true
cp -f $STRUCTS_DIR/*.tsv $REPO_STRUCTS/ 2>/dev/null || true

cd $HOME/repo
git add pipeline_results/v3_bispecific_validated/
git commit -m "results: top-10 structure CIFs (RFd3 backbone + Boltz-2 complex)

For each top design (ranked by √(ipTM_VH1 × ipTM_Nb)):
  *_rfd3_backbone.cif   — RFdiffusion3 designed backbone
  *_boltz2_complex.cif  — Boltz-2 predicted 3-chain complex

Open in PyMOL/ChimeraX to visualize the bispecific bridging minibinder
between VH1-CH1-face (chain A) and nanobody CDRs (chain C)." 2>&1

GH_TOKEN_VAR=$(cat /tmp/gh_token 2>/dev/null)
git push "https://x-access-token:${GH_TOKEN_VAR}@github.com/Cookiemaster33/binary_antibodies.git" \
    cursor/conditional-nanobody-design-992c 2>&1 | tail -3

echo "Structures pushed to GitHub ✓"
echo "View at: https://github.com/Cookiemaster33/binary_antibodies/pull/1"
