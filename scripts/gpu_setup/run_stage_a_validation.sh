#!/usr/bin/env bash
# ============================================================
# run_stage_a_validation.sh
# -------------------------
# Stage A minibinder validation: ProteinMPNN redesign + Boltz-2
# self-consistency on a Lambda GPU instance.
#
# Input:  $PIPELINE/top_structures/*.pdb
#         (chains A=VH, B=VL, C=CH1, D=CL, M=minibinder)
# Output: $PIPELINE/stage_a_validation/
#           stage_a_boltz2_sc.csv    — all results ranked by scRMSD
#           stage_a_validated.fasta  — passing candidates (scRMSD < 2.5 Å)
#
# Tools:
#   ProteinMPNN → run inside rosettacommons/foundry Docker (atomworks MPNN)
#   Boltz-2     → pip install boltz[cuda], run on host
# ============================================================
set -eo pipefail
export PATH="$HOME/.local/bin:$PATH"

PIPELINE=/home/ubuntu/pipeline
PDB_DIR=$PIPELINE/top_structures
MPNN_OUT=$PIPELINE/mpnn_outputs
BOLTZ_IN=$PIPELINE/boltz_sc_inputs
BOLTZ_OUT=$PIPELINE/boltz_sc_outputs
RESULTS=$PIPELINE/stage_a_validation

mkdir -p "$MPNN_OUT" "$BOLTZ_IN" "$BOLTZ_OUT" "$RESULTS"

N_MPNN_SEQS=${N_MPNN_SEQS:-8}

echo "===== Stage A validation: $(date) ====="
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo n/a)"
echo "MPNN seqs per design: $N_MPNN_SEQS"

N_PDBS=$(ls "$PDB_DIR"/*.pdb 2>/dev/null | wc -l)
echo "Found $N_PDBS PDBs"
if [ "$N_PDBS" -eq 0 ]; then echo "ERROR: no PDB files in $PDB_DIR"; exit 1; fi


# ── 1. Docker + Boltz-2 setup (parallel) ─────────────────────
echo ""
echo "=== [1/5] Setup: Docker + Boltz-2 ==="

# Docker access
sudo chmod 666 /var/run/docker.sock 2>/dev/null || true
docker --version

# Pull foundry image (for MPNN)
if docker images rosettacommons/foundry | grep -q "foundry"; then
    echo "  Foundry image already present."
else
    echo "  Pulling rosettacommons/foundry:latest..."
    docker pull rosettacommons/foundry:latest 2>&1 | tail -5
fi

# Install Boltz-2 in a clean Python venv to avoid conflicts with Lambda's
# system packages (numpy, matplotlib, torchvision, etc. all have version
# mismatches that break Boltz when pip-installed into the system Python).
BOLTZ_VENV="$HOME/boltz_venv"
BOLTZ="$BOLTZ_VENV/bin/boltz"

if [ ! -f "$BOLTZ" ]; then
    echo "  Creating Boltz-2 venv at $BOLTZ_VENV..."
    python3 -m venv "$BOLTZ_VENV"
    PIP="$BOLTZ_VENV/bin/pip"
    echo "  Installing torch (cu124)..."
    "$PIP" install -q --upgrade pip 2>&1 | tail -2
    "$PIP" install -q \
        torch torchvision \
        --index-url https://download.pytorch.org/whl/cu124 2>&1 | tail -5
    echo "  Installing cuequivariance (Boltz-2 triangular-mult kernel)..."
    "$PIP" install -q \
        cuequivariance-torch cuequivariance \
        --extra-index-url https://pypi.nvidia.com 2>&1 | tail -5 || \
    "$PIP" install -q cuequivariance-torch cuequivariance 2>&1 | tail -5 || \
        echo "  WARNING: cuequivariance not installed — may use fallback"
    echo "  Installing boltz..."
    "$PIP" install -q boltz 2>&1 | tail -5
fi

# Verify GPU access
"$BOLTZ_VENV/bin/python3" -c "
import torch
print(f'  torch {torch.__version__}  CUDA={torch.cuda.is_available()}  '
      f'device={torch.cuda.get_device_name(0) if torch.cuda.is_available() else \"n/a\"}')
"
echo "  Boltz-2: $BOLTZ"

# Install biopython (host-side, for sequence extraction)
pip install -q biopython 2>&1 | tail -2 || true


# ── 2. ProteinMPNN via foundry Docker ────────────────────────
echo ""
echo "=== [2/5] ProteinMPNN redesign (foundry Docker) ==="

# Write the MPNN python script that runs inside the container
cat > "$PIPELINE/run_mpnn_multi.py" << 'PYEOF'
"""
Run ProteinMPNN on multi-chain Stage A PDB files.
  - Design chain: M (minibinder)
  - Fixed chains: A (VH), B (VL), C (CH1), D (CL)
  - N_SEQS sequences per design
Runs inside rosettacommons/foundry Docker (has atomworks + mpnn packages).
"""
import json, os, sys
from pathlib import Path

PDB_DIR  = Path("/workspace/top_structures")
MPNN_OUT = Path("/workspace/mpnn_outputs")
N_SEQS   = int(os.environ.get("N_MPNN_SEQS", 8))
MPNN_OUT.mkdir(parents=True, exist_ok=True)

from atomworks.io.utils.io_utils import load_any
from mpnn.inference_engines.mpnn import MPNNInferenceEngine

engine = MPNNInferenceEngine(
    model_type="protein_mpnn",
    is_legacy_weights=True,
    out_directory=None,
    write_structures=False,
    write_fasta=False,
)

all_results = []
pdbs = sorted(PDB_DIR.glob("*.pdb"))
print(f"Processing {len(pdbs)} PDBs × {N_SEQS} MPNN sequences...")

for idx, pdb_file in enumerate(pdbs):
    name = pdb_file.stem
    try:
        raw = load_any(str(pdb_file))
        aa  = raw[0] if hasattr(raw, "__getitem__") else raw

        # Identify chain M residue numbers to design
        chain_m_mask = aa.chain_id == "M"
        if chain_m_mask.sum() == 0:
            print(f"  WARNING: no chain M in {name}")
            continue
        chain_m_res = sorted(set(aa.res_id[chain_m_mask].tolist()))
        designed    = [f"M{r}" for r in chain_m_res]

        # Run MPNN: full multi-chain array, only M designed
        result = engine.run(
            atom_arrays=[aa],
            input_dicts=[{
                "batch_size": N_SEQS,
                "designed_residues": designed,
                "remove_waters": True,
            }],
        )

        results = result if isinstance(result, list) else [result]
        for sample_idx, r in enumerate(results):
            full_seq = r.output_dict.get("designed_sequence", "")
            seq_rec  = float(r.output_dict.get("sequence_recovery", 0.0))

            # Extract chain M sequence from the designed output.
            # MPNNInferenceEngine with multi-chain: the output sequence
            # contains ALL chains concatenated. Chain M is last (A,B,C,D,M order).
            # Compute fixed-chain total length to slice.
            fixed_len = int(chain_m_mask.__invert__().sum() // 1)  # non-M atoms → use res count
            # Use unique residue count per chain for accurate slicing
            chain_lens = {}
            for cid in ["A", "B", "C", "D", "M"]:
                mask = aa.chain_id == cid
                chain_lens[cid] = len(set(aa.res_id[mask].tolist())) if mask.sum() > 0 else 0
            fixed_total = sum(chain_lens[c] for c in ["A", "B", "C", "D"])
            m_len = chain_lens["M"]

            if len(full_seq) == m_len:
                mb_seq = full_seq
            elif len(full_seq) >= fixed_total + m_len:
                mb_seq = full_seq[fixed_total: fixed_total + m_len]
            else:
                mb_seq = full_seq[-m_len:] if m_len > 0 else full_seq

            all_results.append({
                "backbone": name,
                "sample": sample_idx + 1,
                "mb_sequence": mb_seq,
                "mb_len": len(mb_seq),
                "sequence_recovery": round(seq_rec, 4),
                "chain_lens": chain_lens,
            })

        if (idx + 1) % 5 == 0 or (idx + 1) == len(pdbs):
            print(f"  {idx+1}/{len(pdbs)} done", flush=True)

    except Exception as exc:
        print(f"  ERROR on {name}: {exc}", file=sys.stderr)

json.dump(all_results, open(MPNN_OUT / "mpnn_sequences.json", "w"), indent=2)

# Write per-backbone FASTA files (for reference)
from itertools import groupby
by_bb = {}
for r in all_results:
    by_bb.setdefault(r["backbone"], []).append(r)
(MPNN_OUT / "seqs").mkdir(exist_ok=True)
for bb, entries in by_bb.items():
    with open(MPNN_OUT / "seqs" / f"{bb}.fa", "w") as f:
        for e in entries:
            f.write(f">sample={e['sample']} rec={e['sequence_recovery']:.4f}\n"
                    f"{e['mb_sequence']}\n")

print(f"\nMPNN done: {len(all_results)} sequences from {len(pdbs)} designs → {MPNN_OUT}")
PYEOF

# Run MPNN inside foundry Docker
docker run --rm --gpus all \
    -e N_MPNN_SEQS="$N_MPNN_SEQS" \
    -v "$PIPELINE:/workspace" \
    rosettacommons/foundry:latest \
    python3 /workspace/run_mpnn_multi.py \
    2>&1 | tee "$PIPELINE/mpnn.log"

N_SEQS_DONE=$(python3 -c "
import json
from pathlib import Path
try:
    d = json.load(open('$PIPELINE/mpnn_outputs/mpnn_sequences.json'))
    print(len(d))
except:
    print(0)
")
echo "  MPNN complete: $N_SEQS_DONE sequences generated"


# ── 3. Build Boltz-2 YAML inputs ──────────────────────────────
echo ""
echo "=== [3/5] Building Boltz-2 YAML inputs ==="

python3 << PYEOF
import json
from pathlib import Path

PDB_DIR    = Path("$PDB_DIR")
MPNN_OUT   = Path("$MPNN_OUT")
BOLTZ_IN   = Path("$BOLTZ_IN")
RESULTS    = Path("$RESULTS")
BOLTZ_IN.mkdir(exist_ok=True)

def chain_seq_from_pdb(pdb_path, chain_id):
    AA3 = {
        "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q",
        "GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K",
        "MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",
        "TYR":"Y","VAL":"V",
    }
    seen, seq = set(), []
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            chain     = line[21]
            resname   = line[17:20].strip()
            resnum    = int(line[22:26])
            if chain == chain_id and atom_name == "CA" and resnum not in seen:
                seen.add(resnum)
                seq.append(AA3.get(resname, "X"))
    return "".join(seq)

mpnn_data = json.loads((MPNN_OUT / "mpnn_sequences.json").read_text())

# Build per-backbone VL and CH1 sequences (same for all samples of a backbone)
fab_seqs = {}
for pdb in sorted(PDB_DIR.glob("*.pdb")):
    fab_seqs[pdb.stem] = {
        "vl":  chain_seq_from_pdb(str(pdb), "B"),
        "ch1": chain_seq_from_pdb(str(pdb), "C"),
    }

input_record = {}
yaml_count = 0

for entry in mpnn_data:
    bb     = entry["backbone"]
    sample = entry["sample"]
    mb_seq = entry["mb_sequence"]
    if bb not in fab_seqs:
        continue
    vl_seq  = fab_seqs[bb]["vl"]
    ch1_seq = fab_seqs[bb]["ch1"]

    yaml_stem = f"{bb}__s{sample:02d}"
    yaml_path = BOLTZ_IN / f"{yaml_stem}.yaml"
    yaml_path.write_text(
        f'sequences:\n'
        f'  - protein:\n'
        f'      id: A\n'
        f'      sequence: "{mb_seq}"\n'
        f'      msa: empty\n'
        f'  - protein:\n'
        f'      id: B\n'
        f'      sequence: "{vl_seq}"\n'
        f'      msa: empty\n'
        f'  - protein:\n'
        f'      id: C\n'
        f'      sequence: "{ch1_seq}"\n'
        f'      msa: empty\n'
    )
    input_record[yaml_stem] = {
        "pdb": bb,
        "sample": sample,
        "mb_seq": mb_seq,
        "vl_seq": vl_seq,
        "ch1_seq": ch1_seq,
        "seq_rec": entry.get("sequence_recovery", 0.0),
    }
    yaml_count += 1

(RESULTS / "input_record.json").write_text(json.dumps(input_record, indent=2))
print(f"  Built {yaml_count} YAML inputs → {BOLTZ_IN}")
PYEOF

N_YAML=$(ls "$BOLTZ_IN"/*.yaml 2>/dev/null | wc -l)
echo "  Total YAML inputs: $N_YAML"
if [ "$N_YAML" -eq 0 ]; then
    echo "ERROR: 0 YAML inputs built — check mpnn.log"; exit 1
fi


# ── 4. Boltz-2 self-consistency predictions ───────────────────
echo ""
echo "=== [4/5] Boltz-2 self-consistency ($N_YAML predictions) ==="

"$BOLTZ" predict "$BOLTZ_IN" \
    --out_dir "$BOLTZ_OUT" \
    --devices 1 \
    --num_workers 4 \
    --override \
    2>&1 | tee "$PIPELINE/boltz_sc.log"

echo "  Boltz-2 complete."


# ── 5. Compute scRMSD + ipTM and rank ─────────────────────────
echo ""
echo "=== [5/5] scRMSD + ipTM ranking ==="

python3 << PYEOF
import json, math, csv
import numpy as np
from pathlib import Path

PIPELINE  = Path("$PIPELINE")
PDB_DIR   = Path("$PDB_DIR")
BOLTZ_OUT = Path("$BOLTZ_OUT")
RESULTS   = Path("$RESULTS")
INPUT_REC = json.loads((RESULTS / "input_record.json").read_text())


def parse_pdb_ca(path, chain_id):
    coords = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() == "CA" and line[21] == chain_id:
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords) if coords else np.zeros((0, 3))


def parse_cif_ca(cif_path, chain_id):
    """Parse Cα coords from Boltz-2 mmCIF output (auth_asym_id = chain label)."""
    coords = []
    col_map = {}
    header_done = False
    with open(cif_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("_atom_site."):
                field = line.split(".")[1].strip()
                col_map[field] = len(col_map)
                continue
            if col_map and not header_done and not line.startswith("_"):
                header_done = True
            if not header_done:
                continue
            parts = line.split()
            if not parts or parts[0] not in ("ATOM", "HETATM"):
                continue
            try:
                atom_id   = parts[col_map.get("label_atom_id",   2)]
                auth_asym = parts[col_map.get("auth_asym_id",    22 if len(parts) > 22 else 6)]
                label_asym = parts[col_map.get("label_asym_id",  6)]
                x = float(parts[col_map.get("Cartn_x", 10)])
                y = float(parts[col_map.get("Cartn_y", 11)])
                z = float(parts[col_map.get("Cartn_z", 12)])
                use_chain = auth_asym if auth_asym not in (".", "?") else label_asym
                if atom_id == "CA" and use_chain == chain_id:
                    coords.append([x, y, z])
            except (IndexError, ValueError):
                pass
    return np.array(coords) if coords else np.zeros((0, 3))


def kabsch_superpose(P, Q):
    """Rotate/translate P to best fit Q. Return (R, t)."""
    P_c = P - P.mean(0); Q_c = Q - Q.mean(0)
    H   = P_c.T @ Q_c
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = Q.mean(0) - P.mean(0) @ R.T
    return R, t


rows = []
conf_files = sorted(BOLTZ_OUT.rglob("confidence_*.json"))
print(f"  Confidence files found: {len(conf_files)}")

for conf_f in conf_files:
    yaml_stem = conf_f.parent.parent.name
    if yaml_stem not in INPUT_REC:
        continue
    rec      = INPUT_REC[yaml_stem]
    pdb_name = rec["pdb"]

    # Parse confidence
    try:
        conf      = json.loads(conf_f.read_text())
        raw_iptm  = conf.get("iptm", {})
        if isinstance(raw_iptm, dict):
            iptm_mb_vl  = float(raw_iptm.get("AB", 0.0))
            iptm_mb_ch1 = float(raw_iptm.get("AC", 0.0))
        else:
            iptm_mb_vl = iptm_mb_ch1 = float(raw_iptm)
        plddt = float(conf.get("mean_plddt", conf.get("complex_plddt", 0.0)))
        ptm   = float(conf.get("ptm", 0.0))
    except Exception as e:
        print(f"  conf error {yaml_stem}: {e}")
        continue

    # scRMSD
    sc_rmsd = float("nan")
    cif_files = sorted(conf_f.parent.glob("*.cif"))
    if not cif_files:
        cif_files = sorted(conf_f.parent.parent.rglob("*.cif"))

    if cif_files:
        cif_f     = cif_files[0]
        input_pdb = PDB_DIR / f"{pdb_name}.pdb"
        if input_pdb.exists():
            ref_vl  = parse_pdb_ca(str(input_pdb), "B")
            ref_ch1 = parse_pdb_ca(str(input_pdb), "C")
            ref_mb  = parse_pdb_ca(str(input_pdb), "M")
            pred_vl  = parse_cif_ca(str(cif_f), "B")
            pred_ch1 = parse_cif_ca(str(cif_f), "C")
            pred_mb  = parse_cif_ca(str(cif_f), "A")  # minibinder is chain A in YAML
            n_vl  = min(len(ref_vl),  len(pred_vl))
            n_ch1 = min(len(ref_ch1), len(pred_ch1))
            n_mb  = min(len(ref_mb),  len(pred_mb))
            if n_vl + n_ch1 >= 20 and n_mb >= 5:
                pred_anc = np.vstack([pred_vl[:n_vl],  pred_ch1[:n_ch1]])
                ref_anc  = np.vstack([ref_vl[:n_vl],   ref_ch1[:n_ch1]])
                R, t     = kabsch_superpose(pred_anc, ref_anc)
                pred_mb_aln = pred_mb[:n_mb] @ R.T + t
                diff    = pred_mb_aln - ref_mb[:n_mb]
                sc_rmsd = float(np.sqrt((diff ** 2).sum() / n_mb))

    rows.append({
        "yaml_stem":   yaml_stem,
        "pdb":         pdb_name,
        "sample":      rec["sample"],
        "seq_rec":     round(rec.get("seq_rec", 0.0), 4),
        "mb_seq":      rec["mb_seq"],
        "mb_len":      len(rec["mb_seq"]),
        "sc_rmsd_A":   round(sc_rmsd, 3) if not math.isnan(sc_rmsd) else float("nan"),
        "iptm_mb_vl":  round(iptm_mb_vl,  3),
        "iptm_mb_ch1": round(iptm_mb_ch1, 3),
        "mean_iptm":   round((iptm_mb_vl + iptm_mb_ch1) / 2, 3),
        "plddt":       round(plddt, 1),
        "ptm":         round(ptm, 3),
    })

valid   = [r for r in rows if not math.isnan(r["sc_rmsd_A"])]
invalid = [r for r in rows if math.isnan(r["sc_rmsd_A"])]
valid.sort(key=lambda r: (r["sc_rmsd_A"], -r["mean_iptm"]))

print(f"\n  Total: {len(rows)}  Valid: {len(valid)}  No-CIF: {len(invalid)}")
print(f"\n{'Rank':<5} {'scRMSD':>7} {'ipTM_VL':>8} {'ipTM_CH1':>9} {'pLDDT':>7} {'seqrec':>7}  Design")
print("  " + "-" * 78)
for i, r in enumerate(valid[:30], 1):
    flag = " ★" if r["sc_rmsd_A"] < 2.0 and r["mean_iptm"] > 0.5 else ""
    print(f"{i:<5} {r['sc_rmsd_A']:>7.3f} {r['iptm_mb_vl']:>8.3f}"
          f" {r['iptm_mb_ch1']:>9.3f} {r['plddt']:>7.1f} {r['seq_rec']:>7.4f}  {r['yaml_stem']}{flag}")

csv_path = RESULTS / "stage_a_boltz2_sc.csv"
fieldnames = ["rank", "yaml_stem", "pdb", "sample", "mb_len", "seq_rec",
              "sc_rmsd_A", "iptm_mb_vl", "iptm_mb_ch1", "mean_iptm",
              "plddt", "ptm", "mb_seq"]
all_out = valid + invalid
for i, r in enumerate(all_out, 1):
    r["rank"] = i
with open(csv_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(all_out)
print(f"\n  CSV: {csv_path}")

passing = [r for r in valid if r["sc_rmsd_A"] < 2.5]
fasta_path = RESULTS / "stage_a_validated.fasta"
with open(fasta_path, "w") as f:
    for i, r in enumerate(passing, 1):
        f.write(f">rank{i:03d}_{r['pdb']}__s{r['sample']:02d}"
                f"__scRMSD{r['sc_rmsd_A']:.2f}__ipTM{r['mean_iptm']:.2f}\n"
                f"{r['mb_seq']}\n")
print(f"  FASTA: {fasta_path}  ({len(passing)} with scRMSD < 2.5 Å)")

if valid:
    rmsds = [r["sc_rmsd_A"] for r in valid]
    iptms = [r["mean_iptm"]  for r in valid]
    print(f"\n  scRMSD — median: {sorted(rmsds)[len(rmsds)//2]:.2f} Å"
          f"  best: {min(rmsds):.2f} Å  >2.0 Å: {sum(x>2.0 for x in rmsds)}/{len(rmsds)}")
    print(f"  ipTM   — median: {sorted(iptms)[len(iptms)//2]:.3f}"
          f"  best: {max(iptms):.3f}  >0.5: {sum(x>0.5 for x in iptms)}/{len(iptms)}")
PYEOF

echo ""
echo "===== Stage A validation complete: $(date) ====="
