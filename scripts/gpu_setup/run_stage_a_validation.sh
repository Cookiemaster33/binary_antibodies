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
#           stage_a_boltz2_sc.csv    — all results
#           stage_a_validated.fasta  — passing candidates
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


# ── 1. Install ProteinMPNN ────────────────────────────────────
echo ""
echo "=== [1/5] ProteinMPNN setup ==="
MPNN_REPO="$HOME/ProteinMPNN"
if [ ! -d "$MPNN_REPO" ]; then
    echo "  Cloning ProteinMPNN..."
    git clone --depth 1 https://github.com/dauparas/ProteinMPNN "$MPNN_REPO" 2>&1 | tail -3
fi
pip install -q biopython 2>&1 | tail -2 || true
echo "  ProteinMPNN ready."


# ── 2. Install Boltz-2 ────────────────────────────────────────
echo ""
echo "=== [2/5] Boltz-2 setup ==="
if ! command -v boltz &>/dev/null; then
    pip install -q "boltz[cuda]" 2>&1 | tail -5 \
      || pip install --break-system-packages -q "boltz[cuda]" 2>&1 | tail -5
fi
BOLTZ=$(command -v boltz 2>/dev/null || echo "$HOME/.local/bin/boltz")
echo "  Boltz: $BOLTZ"


# ── 3. Run ProteinMPNN ────────────────────────────────────────
echo ""
echo "=== [3/5] ProteinMPNN redesign ==="

python3 << PYEOF
import json, subprocess, sys
from pathlib import Path

pdb_dir    = Path("$PDB_DIR")
mpnn_out   = Path("$MPNN_OUT")
mpnn_repo  = Path("$MPNN_REPO")
n_seqs     = int("$N_MPNN_SEQS")
chain_jsonl = mpnn_out / "chain_ids.jsonl"

# Build chain_id_jsonl
entries = {}
for pdb in sorted(pdb_dir.glob("*.pdb")):
    entries[pdb.stem] = {
        "designed_chain_list": ["M"],
        "fixed_chain_list": ["A", "B", "C", "D"],
    }
with open(chain_jsonl, "w") as f:
    for name, val in entries.items():
        f.write(json.dumps({name: val}) + "\n")
print(f"  chain_ids.jsonl: {len(entries)} entries")

# Run ProteinMPNN
cmd = [
    sys.executable, str(mpnn_repo / "protein_mpnn_run.py"),
    "--pdb_path",         str(pdb_dir),
    "--chain_id_jsonl",   str(chain_jsonl),
    "--out_folder",       str(mpnn_out),
    "--num_seq_per_target", str(n_seqs),
    "--sampling_temp",    "0.1",
    "--seed",             "37",
    "--batch_size",       "1",
]
print("  Running:", " ".join(cmd))
r = subprocess.run(cmd, capture_output=False)
if r.returncode != 0:
    print(f"  WARNING: ProteinMPNN exited {r.returncode}")
else:
    seqs_dir = mpnn_out / "seqs"
    n_fa = len(list(seqs_dir.glob("*.fa"))) if seqs_dir.exists() else 0
    print(f"  Done: {n_fa} FASTA files in {seqs_dir}")
PYEOF


# ── 4. Build Boltz-2 YAML inputs ──────────────────────────────
echo ""
echo "=== [4/5] Building Boltz-2 YAML inputs ==="

python3 << PYEOF
import json, re
from pathlib import Path

try:
    from Bio import PDB as biopdb, SeqUtils
    def chain_seq(struct, cid):
        chain = struct[0][cid]
        res = [r for r in chain.get_residues() if biopdb.is_aa(r)]
        return SeqUtils.seq1("".join(r.resname for r in res))
    parser = biopdb.PDBParser(QUIET=True)
    USE_BIO = True
except ImportError:
    USE_BIO = False

def chain_seq_raw(pdb_path, chain_id):
    """Extract Cα residue sequence directly from PDB ATOM records."""
    AA3 = {
        "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q",
        "GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K",
        "MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",
        "TYR":"Y","VAL":"V",
    }
    seen = set()
    seq = []
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

pdb_dir    = Path("$PDB_DIR")
mpnn_out   = Path("$MPNN_OUT")
boltz_in   = Path("$BOLTZ_IN")
results    = Path("$RESULTS")
boltz_in.mkdir(exist_ok=True)

# Parse MPNN FASTA outputs
# Layout: mpnn_out/seqs/<pdb_stem>.fa
# Header: >score=..., global_score=..., sample=N, T=0.1, seq_rec=..., pdb=...
# Sequence: full concatenated sequence of all chains

def parse_mpnn_fa(fa_path, chain_lengths):
    """
    Parse MPNN output FASTA and extract the designed-chain (M) sequence.
    chain_lengths: dict {chain_id: length} for chains in PDB order
    Returns list of (sample_idx, mb_sequence)
    """
    entries = []
    lines = fa_path.read_text().splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith(">"):
            header = line
            i += 1
            seq = ""
            while i < len(lines) and not lines[i].startswith(">"):
                seq += lines[i].strip()
                i += 1
            m = re.search(r"sample=(\d+)", header)
            sample = int(m.group(1)) if m else -1
            if sample == 0:
                continue  # sample=0 is the native reference; skip

            # If sequence length == M chain length: just the designed chain
            m_len = chain_lengths.get("M", 0)
            total_fixed = sum(v for k, v in chain_lengths.items() if k != "M")

            if len(seq) == m_len:
                mb_seq = seq
            elif len(seq) == total_fixed + m_len:
                # Full concatenated sequence — M chain is at the end
                # (PDB chain order: A, B, C, D, M)
                mb_seq = seq[total_fixed:]
            else:
                # Unexpected length; try to take last m_len chars
                mb_seq = seq[-m_len:] if m_len > 0 else seq
            entries.append((sample, mb_seq))
        else:
            i += 1
    return entries

input_record = {}
yaml_count = 0

for pdb in sorted(pdb_dir.glob("*.pdb")):
    pdb_stem = pdb.stem
    fa_path = mpnn_out / "seqs" / f"{pdb_stem}.fa"
    if not fa_path.exists():
        print(f"  WARNING: no FASTA for {pdb_stem}")
        continue

    # Get chain lengths from input PDB
    chain_lengths = {}
    for cid in ["A", "B", "C", "D", "M"]:
        seq = chain_seq_raw(str(pdb), cid)
        chain_lengths[cid] = len(seq)

    # Get VL and CH1 sequences for Boltz-2 YAML
    vl_seq  = chain_seq_raw(str(pdb), "B")
    ch1_seq = chain_seq_raw(str(pdb), "C")

    seqs = parse_mpnn_fa(fa_path, chain_lengths)
    if not seqs:
        print(f"  WARNING: no sequences parsed from {fa_path}")
        continue

    for sample, mb_seq in seqs:
        yaml_stem = f"{pdb_stem}__s{sample:02d}"
        yaml_path = boltz_in / f"{yaml_stem}.yaml"
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
            "pdb": pdb_stem,
            "sample": sample,
            "mb_seq": mb_seq,
            "vl_seq": vl_seq,
            "ch1_seq": ch1_seq,
            "chain_lengths": chain_lengths,
        }
        yaml_count += 1

(results / "input_record.json").write_text(json.dumps(input_record, indent=2))
print(f"  Built {yaml_count} YAML inputs → {boltz_in}")
PYEOF

N_YAML=$(ls "$BOLTZ_IN"/*.yaml 2>/dev/null | wc -l)
echo "  Total YAML inputs: $N_YAML"


# ── 5. Boltz-2 self-consistency predictions ───────────────────
echo ""
echo "=== [5/5] Boltz-2 self-consistency ($N_YAML predictions) ==="

"$BOLTZ" predict "$BOLTZ_IN" \
    --out_dir "$BOLTZ_OUT" \
    --devices 1 \
    --num_workers 4 \
    --override \
    2>&1 | tee "$PIPELINE/boltz_sc.log"

echo "  Boltz-2 complete."


# ── 6. Compute scRMSD + ipTM and rank ─────────────────────────
echo ""
echo "=== Results: scRMSD + ipTM ranking ==="

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
    """Parse Cα coordinates from Boltz-2 mmCIF output."""
    coords = []
    in_atom_site = False
    col_map = {}
    with open(cif_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("_atom_site."):
                col_name = line.strip().split(".")[1]
                col_map[col_name] = len(col_map)
                in_atom_site = True
                continue
            if in_atom_site and (line.startswith("#") or line.startswith("_")):
                in_atom_site = False
            if not line or line.startswith("#") or line.startswith("_"):
                continue
            parts = line.split()
            if not parts or parts[0] not in ("ATOM", "HETATM"):
                continue
            # Fallback: try positional parsing if col_map is available
            try:
                if col_map:
                    label_atom = parts[col_map.get("label_atom_id", 2)]
                    label_chain = parts[col_map.get("label_asym_id", 6)]
                    auth_chain  = parts[col_map.get("auth_asym_id", 22 if len(parts) > 22 else 6)]
                    x = float(parts[col_map.get("Cartn_x", 10)])
                    y = float(parts[col_map.get("Cartn_y", 11)])
                    z = float(parts[col_map.get("Cartn_z", 12)])
                    use_chain = auth_chain if auth_chain != "." else label_chain
                else:
                    # Positional fallback
                    label_atom = parts[2]
                    use_chain  = parts[6]
                    x, y, z = float(parts[10]), float(parts[11]), float(parts[12])
                if label_atom == "CA" and use_chain == chain_id:
                    coords.append([x, y, z])
            except (IndexError, ValueError):
                pass
    return np.array(coords) if coords else np.zeros((0, 3))


def kabsch_rmsd_and_rot(P, Q):
    """Return (R, t, rmsd) that maps P → Q (minimum RMSD)."""
    P_c = P - P.mean(0)
    Q_c = Q - Q.mean(0)
    H = P_c.T @ Q_c
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = Q.mean(0) - P.mean(0) @ R.T
    P_rot = P @ R.T + t
    rmsd = float(np.sqrt(((P_rot - Q) ** 2).sum() / len(P)))
    return R, t, rmsd


rows = []
conf_files = list(BOLTZ_OUT.rglob("confidence_*.json"))
print(f"  Confidence files found: {len(conf_files)}")

for conf_f in sorted(conf_files):
    yaml_stem = conf_f.parent.parent.name
    if yaml_stem not in INPUT_REC:
        continue
    rec = INPUT_REC[yaml_stem]
    pdb_name = rec["pdb"]

    # ipTM
    try:
        conf = json.loads(conf_f.read_text())
        raw_iptm = conf.get("iptm", {})
        if isinstance(raw_iptm, dict):
            iptm_mb_vl  = float(raw_iptm.get("AB", 0.0))
            iptm_mb_ch1 = float(raw_iptm.get("AC", 0.0))
        else:
            iptm_mb_vl = iptm_mb_ch1 = float(raw_iptm)
        plddt = float(conf.get("mean_plddt", conf.get("complex_plddt", 0.0)))
        ptm   = float(conf.get("ptm", 0.0))
    except Exception as e:
        print(f"  conf parse error {yaml_stem}: {e}")
        continue

    # scRMSD
    sc_rmsd = float("nan")
    cif_files = sorted(conf_f.parent.glob("*.cif"))
    if not cif_files:
        cif_files = sorted(conf_f.parent.parent.rglob("*.cif"))

    if cif_files:
        cif_f = cif_files[0]
        input_pdb = PDB_DIR / f"{pdb_name}.pdb"
        if input_pdb.exists():
            ref_vl  = parse_pdb_ca(str(input_pdb), "B")
            ref_ch1 = parse_pdb_ca(str(input_pdb), "C")
            ref_mb  = parse_pdb_ca(str(input_pdb), "M")

            pred_vl  = parse_cif_ca(str(cif_f), "B")
            pred_ch1 = parse_cif_ca(str(cif_f), "C")
            pred_mb  = parse_cif_ca(str(cif_f), "A")

            n_vl  = min(len(ref_vl),  len(pred_vl))
            n_ch1 = min(len(ref_ch1), len(pred_ch1))
            n_mb  = min(len(ref_mb),  len(pred_mb))

            if n_vl + n_ch1 >= 20 and n_mb >= 5:
                pred_anc = np.vstack([pred_vl[:n_vl],  pred_ch1[:n_ch1]])
                ref_anc  = np.vstack([ref_vl[:n_vl],   ref_ch1[:n_ch1]])
                R, t, _ = kabsch_rmsd_and_rot(pred_anc, ref_anc)
                pred_mb_aln = pred_mb[:n_mb] @ R.T + t
                diff = pred_mb_aln - ref_mb[:n_mb]
                sc_rmsd = float(np.sqrt((diff ** 2).sum() / n_mb))

    mean_iptm = (iptm_mb_vl + iptm_mb_ch1) / 2
    rows.append({
        "yaml_stem":   yaml_stem,
        "pdb":         pdb_name,
        "sample":      rec["sample"],
        "mb_seq":      rec["mb_seq"],
        "mb_len":      len(rec["mb_seq"]),
        "sc_rmsd_A":   round(sc_rmsd, 3) if not math.isnan(sc_rmsd) else float("nan"),
        "iptm_mb_vl":  round(iptm_mb_vl,  3),
        "iptm_mb_ch1": round(iptm_mb_ch1, 3),
        "mean_iptm":   round(mean_iptm, 3),
        "plddt":       round(plddt, 1),
        "ptm":         round(ptm, 3),
    })

valid   = [r for r in rows if not math.isnan(r["sc_rmsd_A"])]
invalid = [r for r in rows if math.isnan(r["sc_rmsd_A"])]
valid.sort(key=lambda r: (r["sc_rmsd_A"], -r["mean_iptm"]))

print(f"\n  Total: {len(rows)}  Valid: {len(valid)}  No-structure: {len(invalid)}")
print(f"\n{'Rank':<5} {'scRMSD':>7} {'ipTM_VL':>8} {'ipTM_CH1':>9} {'pLDDT':>7}  Design")
print("  " + "-" * 72)
for i, r in enumerate(valid[:30], 1):
    flag = " ★" if r["sc_rmsd_A"] < 2.0 and r["mean_iptm"] > 0.5 else ""
    print(f"{i:<5} {r['sc_rmsd_A']:>7.3f} {r['iptm_mb_vl']:>8.3f}"
          f" {r['iptm_mb_ch1']:>9.3f} {r['plddt']:>7.1f}  {r['yaml_stem']}{flag}")

csv_path = RESULTS / "stage_a_boltz2_sc.csv"
fieldnames = ["rank", "yaml_stem", "pdb", "sample", "mb_len",
              "sc_rmsd_A", "iptm_mb_vl", "iptm_mb_ch1", "mean_iptm",
              "plddt", "ptm", "mb_seq"]
all_rows = valid + invalid
for i, r in enumerate(all_rows, 1):
    r["rank"] = i
with open(csv_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(all_rows)
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
