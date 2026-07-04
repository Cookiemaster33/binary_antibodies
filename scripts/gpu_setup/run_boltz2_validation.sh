#!/bin/bash
# ============================================================
# run_boltz2_validation.sh
# ------------------------
# Boltz-2 validation of bispecific bridging minibinder designs.
#
# For each top-50 design:
#   1. Build 3-chain YAML (VH1 | Minibinder | Nanobody)
#      with pocket constraints on BOTH binding surfaces
#   2. Run boltz predict (--use_msa_server for VH1+Nb, empty for MB)
#   3. Parse output: scRMSD + ipTM(MB↔VH1) + ipTM(MB↔Nb)
#   4. Rank and save validated_minibinders.fasta
#
# scRMSD: superimpose Boltz-2 structure on VH1+Nb (fixed parts),
#         then measure RMSD of minibinder vs RFd3 designed backbone.
# ============================================================
set -eo pipefail
PIPELINE=/home/ubuntu/pipeline
LOG=$PIPELINE/boltz_validation.log
exec > >(tee -a "$LOG") 2>&1

echo "===== Boltz-2 validation: $(date) ====="
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"

mkdir -p $PIPELINE/boltz_inputs $PIPELINE/boltz_outputs $PIPELINE/final

# ── 1. Install Boltz-2 ───────────────────────────────────────────
echo ""
echo "=== Installing Boltz-2 ==="
if ! command -v boltz &>/dev/null; then
    pip install boltz[cuda] -U -q 2>&1 | tail -3
    echo "Boltz-2 installed."
else
    echo "Boltz-2 already installed: $(boltz --version 2>/dev/null || echo 'ok')"
fi

# ── 2. Generate YAML inputs ──────────────────────────────────────
echo ""
echo "=== Building YAML inputs for top-50 designs ==="

python3 << 'PYEOF'
import json, os
from pathlib import Path

data = json.load(open("/home/ubuntu/pipeline/validation_inputs.json"))
VH1_SEQ = data['VH1_seq']
NB_SEQ  = data['NB_seq']
designs = data['designs']

# VH1-CH1-face contact residues (hotspots on chain A)
VH1_CONTACTS = [11,12,13,14,15,39,40,41,79,80,81,82,83,84,85,
                103,104,105,106,107,108,109,110,111,112,113,114,115]
# Nanobody CDR residues (hotspots on chain C)
NB_CDR = list(range(27,34)) + list(range(52,58)) + list(range(99,113))

out_dir = Path("/home/ubuntu/pipeline/boltz_inputs")
out_dir.mkdir(exist_ok=True)

for r in designs:
    rank = r['rank']
    bb   = r['backbone']
    mb   = r['minibinder_sequence']
    rec  = r['sequence_recovery']
    name = f"rank{rank:02d}_{bb}_rec{rec:.3f}"

    # Pocket constraint contacts: VH1-CH1-face (A) + Nanobody CDRs (C)
    vh1_contacts = [[f"A", res] for res in VH1_CONTACTS]
    nb_contacts  = [[f"C", res] for res in NB_CDR]
    all_contacts = vh1_contacts + nb_contacts

    yaml_content = f"""sequences:
  - protein:
      id: A
      sequence: "{VH1_SEQ}"
      msa: empty
  - protein:
      id: B
      sequence: "{mb}"
      msa: empty
  - protein:
      id: C
      sequence: "{NB_SEQ}"
      msa: empty
constraints:
  - pocket:
      binder: B
      contacts: {str([[c[0], c[1]] for c in all_contacts]).replace("'", "")}
      max_distance: 10.0
"""
    (out_dir / f"{name}.yaml").write_text(yaml_content)

print(f"Wrote {len(designs)} YAML files to {out_dir}")
PYEOF

N_YAML=$(ls $PIPELINE/boltz_inputs/*.yaml | wc -l)
echo "Generated $N_YAML YAML files."

# ── 3. Run Boltz-2 ───────────────────────────────────────────────
echo ""
echo "=== Running Boltz-2 predictions ==="
echo "Predicting $N_YAML 3-chain complexes..."

boltz predict $PIPELINE/boltz_inputs \
    --out_dir $PIPELINE/boltz_outputs \
    --num_workers 2 \
    --devices 1 \
    --override 2>&1 | grep -v "^$" | tee $PIPELINE/boltz_run.log

echo "Boltz-2 predictions complete."

# ── 4. Compute scRMSD + parse ipTM ──────────────────────────────
echo ""
echo "=== Computing scRMSD and parsing ipTM scores ==="

python3 << 'PYEOF'
import json, glob, os
from pathlib import Path
import numpy as np

try:
    import biotite.structure as struc
    import biotite.structure.io.pdbx as pdbx
    HAS_BIOTITE = True
except ImportError:
    HAS_BIOTITE = False

PIPELINE = Path("/home/ubuntu/pipeline")
RFD3_OUT = PIPELINE / "outputs/rfd3"
BOLTZ_OUT = PIPELINE / "boltz_outputs"
FINAL    = PIPELINE / "final"
FINAL.mkdir(exist_ok=True)

data = json.load(open(PIPELINE / "validation_inputs.json"))
designs = data['designs']

def kabsch_rmsd(P, Q):
    """Kabsch-aligned RMSD."""
    P = P - P.mean(0); Q = Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(P.T @ Q)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1,1,d]) @ U.T
    return float(np.sqrt(np.mean(np.sum((P @ R.T - Q)**2, axis=1))))

def get_ca_coords(aa, chain_id=None, res_start=None, res_end=None):
    """Get Cα coordinates from AtomArray with optional filters."""
    mask = aa.atom_name == "CA"
    if chain_id:
        mask &= (aa.chain_id == chain_id)
    if res_start:
        mask &= (aa.res_id >= res_start)
    if res_end:
        mask &= (aa.res_id <= res_end)
    return aa.coord[mask], aa.res_id[mask]

results = []
for r in designs:
    rank = r['rank']
    bb   = r['backbone']
    mb   = r['minibinder_sequence']
    rec  = r['sequence_recovery']
    name = f"rank{rank:02d}_{bb}_rec{rec:.3f}"

    # ── Find Boltz-2 output CIF ──────────────────────────────────
    boltz_cifs = list(BOLTZ_OUT.glob(f"**/{name}/**/*.cif"))
    if not boltz_cifs:
        boltz_cifs = list(BOLTZ_OUT.glob(f"**/{name}*/**/*.cif"))

    boltz_conf = 0.0
    iptm_mb_vh1 = 0.0
    iptm_mb_nb  = 0.0
    sc_rmsd = 999.0

    if boltz_cifs:
        try:
            cif_file = pdbx.CIFFile.read(str(boltz_cifs[0]))
            boltz_aa = pdbx.get_structure(cif_file, model=1, include_bonds=False)
            boltz_conf = float(np.mean(boltz_aa.b_factor[boltz_aa.atom_name=="CA"]))

            # Extract chain Cα coords from Boltz-2 output
            boltz_A_ca, _ = get_ca_coords(boltz_aa, chain_id="A")  # VH1
            boltz_B_ca, _ = get_ca_coords(boltz_aa, chain_id="B")  # Minibinder
            boltz_C_ca, _ = get_ca_coords(boltz_aa, chain_id="C")  # Nanobody

            # ── RFd3 backbone ────────────────────────────────────
            from atomworks.io.utils.io_utils import load_any
            cif_rfd3 = RFD3_OUT / f"{bb}.cif"
            if cif_rfd3.exists():
                raw = load_any(str(cif_rfd3))
                rfd3_aa = raw[0] if hasattr(raw, "__getitem__") else raw
                ch = list(set(rfd3_aa.chain_id))[0]

                rfd3_vh1_ca = rfd3_aa.coord[(rfd3_aa.atom_name=="CA") &
                                             (rfd3_aa.res_id >= 1) & (rfd3_aa.res_id <= 115)]
                rfd3_mb_ca  = rfd3_aa.coord[(rfd3_aa.atom_name=="CA") &
                                             (rfd3_aa.res_id >= 116) & (rfd3_aa.res_id <= 170)]
                rfd3_nb_ca  = rfd3_aa.coord[(rfd3_aa.atom_name=="CA") &
                                             (rfd3_aa.res_id >= 171) & (rfd3_aa.res_id <= 285)]

                # ── Superimpose Boltz-2 on RFd3 fixed parts ──────
                # Concatenate VH1 + Nanobody as alignment anchors
                n_vh1 = min(len(boltz_A_ca), len(rfd3_vh1_ca))
                n_nb  = min(len(boltz_C_ca), len(rfd3_nb_ca))

                if n_vh1 > 10 and n_nb > 10:
                    P_anchor = np.vstack([boltz_A_ca[:n_vh1], boltz_C_ca[:n_nb]])
                    Q_anchor = np.vstack([rfd3_vh1_ca[:n_vh1], rfd3_nb_ca[:n_nb]])

                    # Kabsch rotation for anchor
                    P_c = P_anchor - P_anchor.mean(0)
                    Q_c = Q_anchor - Q_anchor.mean(0)
                    U2, S2, Vt2 = np.linalg.svd(P_c.T @ Q_c)
                    d2 = np.linalg.det(Vt2.T @ U2.T)
                    R_mat = Vt2.T @ np.diag([1,1,d2]) @ U2.T
                    t_P = P_anchor.mean(0)
                    t_Q = Q_anchor.mean(0)

                    # Apply to minibinder
                    n_mb = min(len(boltz_B_ca), len(rfd3_mb_ca))
                    if n_mb > 5:
                        boltz_B_aligned = (boltz_B_ca[:n_mb] - t_P) @ R_mat.T + t_Q
                        sc_rmsd = kabsch_rmsd(boltz_B_aligned, rfd3_mb_ca[:n_mb])

        except Exception as e:
            print(f"  Error for {name}: {e}")

    # ── Parse Boltz confidence JSON for ipTM ─────────────────────
    conf_jsons = list(BOLTZ_OUT.glob(f"**/{name}*/**/*confidence*.json"))
    if conf_jsons:
        try:
            conf = json.load(open(conf_jsons[0]))
            iptm_mb_vh1 = float(conf.get("iptm", {}).get("AB", 0) or
                                conf.get("chain_pair_iptm", {}).get("AB", 0) or
                                conf.get("iptm", 0))
            iptm_mb_nb  = float(conf.get("iptm", {}).get("BC", 0) or
                                conf.get("chain_pair_iptm", {}).get("BC", 0) or
                                conf.get("iptm", 0))
            boltz_conf  = float(conf.get("mean_plddt", boltz_conf) or boltz_conf)
        except Exception as e:
            pass

    passed = sc_rmsd < 2.0 and boltz_conf > 60
    results.append({
        "rank": rank, "backbone": bb,
        "minibinder_sequence": mb,
        "sequence_recovery": rec,
        "sc_rmsd_A": round(sc_rmsd, 3),
        "boltz_mean_plddt": round(boltz_conf, 1),
        "iptm_MB_VH1": round(iptm_mb_vh1, 3),
        "iptm_MB_Nb": round(iptm_mb_nb, 3),
        "pass_filter": passed,
    })
    status = " ✓" if passed else ""
    print(f"  rank{rank:02d} {bb}  scRMSD={sc_rmsd:.2f}Å  pLDDT={boltz_conf:.1f}  "
          f"ipTM(MB↔VH1)={iptm_mb_vh1:.3f}  ipTM(MB↔Nb)={iptm_mb_nb:.3f}{status}")

results.sort(key=lambda x: x["sc_rmsd_A"] if x["sc_rmsd_A"]<900 else 999)

json.dump(results, open(FINAL/"boltz2_sc_results.json","w"), indent=2)
passing = [r for r in results if r["pass_filter"]]
with open(FINAL/"validated_minibinders_boltz2.fasta","w") as f:
    for r in passing:
        f.write(f">{r['backbone']}__scRMSD{r['sc_rmsd_A']:.2f}__pLDDT{r['boltz_mean_plddt']:.0f}"
                f"__ipTM_VH1_{r['iptm_MB_VH1']:.2f}__ipTM_Nb_{r['iptm_MB_Nb']:.2f}\n"
                f"{r['minibinder_sequence']}\n")

print(f"\n=== Boltz-2 Self-Consistency Results ===")
print(f"{'R':>3}  {'Backbone':>30}  {'scRMSD':>7}  {'pLDDT':>6}  {'ipTM↔VH1':>9}  {'ipTM↔Nb':>8}  Pass")
print("  " + "-"*80)
for r in results[:20]:
    print(f"{r['rank']:>3}.  {r['backbone']:>30}  {r['sc_rmsd_A']:>7.3f}  "
          f"{r['boltz_mean_plddt']:>6.1f}  {r['iptm_MB_VH1']:>9.3f}  "
          f"{r['iptm_MB_Nb']:>8.3f}{'  ✓' if r['pass_filter'] else ''}")

print(f"\nPassed: {len(passing)}/{len(results)}")
if passing:
    p = passing[0]
    print(f"Best: {p['backbone']}  scRMSD={p['sc_rmsd_A']:.2f}Å  pLDDT={p['boltz_mean_plddt']:.0f}")
    print(f"Sequence: {p['minibinder_sequence']}")
PYEOF

echo ""
echo "===== Boltz-2 validation complete: $(date) ====="
ls -lh $PIPELINE/final/
