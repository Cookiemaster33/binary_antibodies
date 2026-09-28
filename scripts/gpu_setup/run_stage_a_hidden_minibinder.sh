#!/bin/bash
# ============================================================
# run_stage_a_hidden_minibinder.sh
# Stage A: RFd3 design of an unlinked, alpha-helical minibinder that cross-links
# the CH1 and VL surfaces of the split Fab.
#
# Target:     full Fab A-D + HER2 epitope stub T, each a separate FIXED chain
# Designed:   the minibinder only, on its own chain (chain M after grafting)
# Helicity:   is_non_loopy=true + low-temperature sampling + helical length window
#
# Prerequisite (run locally before GPU launch):
#   python scripts/build_stage_a_design_target.py --fab-pdb <placed Fab>
# ============================================================
set -eo pipefail

PIPELINE=${PIPELINE_DIR:-/home/ubuntu/pipeline}
LOG=$PIPELINE/stage_a_pipeline.log
exec > >(tee -a "$LOG") 2>&1

N_DESIGNS=${N_DESIGNS:-200}
BATCH_SIZE=${BATCH_SIZE:-10}
# Comma-separated length windows. The first is the helical-bridge regime; the
# second is the original compact window, kept so the two can be compared.
MB_LENGTH_RANGES=${MB_LENGTH_RANGES:-60-85,40-55}
IS_NON_LOOPY=${IS_NON_LOOPY:-1}
STEP_SCALE=${STEP_SCALE:-3.0}
GAMMA_0=${GAMMA_0:-0.2}
LOW_MEMORY_MODE=${LOW_MEMORY_MODE:-0}

INPUT_PDB_NAME=${INPUT_PDB:-fab_hidden_minibinder_stage_a.pdb}
CONFIG_JSON=${CONFIG_JSON:-stage_a_hidden_minibinder_config.json}

echo "===== Stage A helical minibinder pipeline: $(date) ====="
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "Designs: $N_DESIGNS | length windows: $MB_LENGTH_RANGES"
echo "Helical conditioning: is_non_loopy=$IS_NON_LOOPY step_scale=$STEP_SCALE gamma_0=$GAMMA_0"

mkdir -p "$PIPELINE/inputs" "$PIPELINE/outputs/rfd3_stage_a" \
         "$PIPELINE/outputs/stage_a_grafted" "$PIPELINE/final"

sudo chmod 666 /var/run/docker.sock 2>/dev/null || true
docker pull rosettacommons/foundry:latest > "$PIPELINE/docker_pull.log" 2>&1 || true

cat > "$PIPELINE/run_rfd3_stage_a.py" << 'PYEOF'
import json, os, torch
from pathlib import Path
from atomworks.io.utils.io_utils import to_cif_file
from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine
from rfd3.inference.input_parsing import DesignInputSpecification

torch.set_float32_matmul_precision("high")

PIPELINE = Path("/workspace")
INPUT  = PIPELINE / "inputs" / os.environ.get("INPUT_PDB", "fab_hidden_minibinder_stage_a.pdb")
CONFIG = PIPELINE / "inputs" / os.environ.get("CONFIG_JSON", "stage_a_hidden_minibinder_config.json")
OUT    = Path("/workspace/outputs/rfd3_stage_a")
N      = int(os.environ.get("N_DESIGNS", 200))
BATCH  = int(os.environ.get("BATCH_SIZE", 10))
RANGES = [r.strip() for r in os.environ.get("MB_LENGTH_RANGES", "60-85").split(",") if r.strip()]
NON_LOOPY = os.environ.get("IS_NON_LOOPY", "1") == "1"
STEP_SCALE = float(os.environ.get("STEP_SCALE", 3.0))
GAMMA_0 = float(os.environ.get("GAMMA_0", 0.2))
OUT.mkdir(parents=True, exist_ok=True)

cfg = json.load(open(CONFIG))
rfd3 = cfg["rfd3"]
base_contig = rfd3["contig"]
default_range = rfd3.get("mb_length_range", RANGES[0])


def contig_for(length_range):
    """Swap the designed-length window at the head of the contig."""
    head, _, tail = base_contig.partition(",/0,")
    assert head == default_range, f"unexpected contig head {head!r} != {default_range!r}"
    return f"{length_range},/0,{tail}"


def build_spec(length_range):
    kwargs = dict(
        input=str(INPUT),
        contig=contig_for(length_range),
        select_hotspots=rfd3["select_hotspots"],
        select_fixed_atoms=rfd3["select_fixed_atoms"],
    )
    if "dialect" in rfd3:
        kwargs["dialect"] = rfd3["dialect"]
    if NON_LOOPY:
        kwargs["is_non_loopy"] = True
    if rfd3.get("ori_token") is not None:
        kwargs["ori_token"] = rfd3["ori_token"]
    elif rfd3.get("infer_ori_strategy"):
        kwargs["infer_ori_strategy"] = rfd3["infer_ori_strategy"]
    try:
        return DesignInputSpecification.safe_init(**kwargs)
    except TypeError as exc:
        # Drop optional keys this build of RFd3 does not accept rather than
        # losing the whole run.
        for key in ("ori_token", "infer_ori_strategy", "is_non_loopy", "dialect"):
            if key in kwargs and key in str(exc):
                print(f"  WARNING: RFd3 rejected {key!r} ({exc}); retrying without it")
                kwargs.pop(key)
                return DesignInputSpecification.safe_init(**kwargs)
        raise


def build_engine():
    config_kwargs = {"diffusion_batch_size": BATCH}
    if os.environ.get("LOW_MEMORY_MODE", "0") == "1":
        config_kwargs["low_memory_mode"] = True
    try:
        engine_kwargs = dict(RFD3InferenceConfig(**config_kwargs))
    except TypeError as exc:
        print(f"  WARNING: RFd3 rejected {sorted(config_kwargs)} ({exc}); using batch size only")
        engine_kwargs = dict(RFD3InferenceConfig(diffusion_batch_size=BATCH))
    sampler = engine_kwargs.get("inference_sampler")
    applied = {}
    for key, value in (("step_scale", STEP_SCALE), ("gamma_0", GAMMA_0)):
        if sampler is None:
            continue
        if hasattr(sampler, key):
            setattr(sampler, key, value)
            applied[key] = value
        elif isinstance(sampler, dict) and key in sampler:
            sampler[key] = value
            applied[key] = value
    print(f"  sampler overrides applied: {applied or 'none (using RFd3 defaults)'}")
    return RFD3InferenceEngine(**engine_kwargs)


print(f"Stage A RFd3: {N} designs across {len(RANGES)} length window(s)")
print(f"  input:    {INPUT}")
print(f"  target:   full Fab A-D + T, all chains fixed; minibinder on its own chain")
print(f"  hotspots: {str(rfd3['select_hotspots'])[:90]}...")
print(f"  ori_token: {rfd3.get('ori_token')}")
print(f"  is_non_loopy: {NON_LOOPY}")

model = build_engine()
per_range = max(1, N // max(1, len(RANGES)))
nbatch = max(1, per_range // BATCH)
saved = 0
for length_range in RANGES:
    tag = length_range.replace("-", "to")
    spec = build_spec(length_range)
    print(f"\n  === length window {length_range} ({nbatch} x {BATCH}) ===", flush=True)
    print(f"  contig: {contig_for(length_range)}")
    for bi in range(nbatch):
        print(f"  Batch {bi+1}/{nbatch} ...", flush=True)
        for _key, designs in model.run(inputs=spec, out_dir=None, n_batches=1).items():
            for i, d in enumerate(designs):
                to_cif_file(d.atom_array, str(OUT / f"sa_{tag}_b{bi:03d}_{i:03d}.cif"))
                saved += 1
print(f"\nStage A RFd3 done: {saved} designs → {OUT}")
PYEOF

echo ""
echo "=== Step 1: RFdiffusion3 Stage A ==="
docker run --rm --gpus all \
    -v "$PIPELINE:/workspace" \
    -e N_DESIGNS=$N_DESIGNS \
    -e BATCH_SIZE=$BATCH_SIZE \
    -e MB_LENGTH_RANGES=$MB_LENGTH_RANGES \
    -e IS_NON_LOOPY=$IS_NON_LOOPY \
    -e STEP_SCALE=$STEP_SCALE \
    -e GAMMA_0=$GAMMA_0 \
    -e LOW_MEMORY_MODE=$LOW_MEMORY_MODE \
    -e INPUT_PDB=$INPUT_PDB_NAME \
    -e CONFIG_JSON=$CONFIG_JSON \
    -e FOUNDRY_CHECKPOINT_DIRS=/weights \
    rosettacommons/foundry:latest \
    python3 /workspace/run_rfd3_stage_a.py

echo "RFd3 Stage A: $(ls $PIPELINE/outputs/rfd3_stage_a/sa_*.cif 2>/dev/null | wc -l) designs"

echo ""
echo "=== Step 2: Graft minibinder onto exact input Fab ==="
PYTHONPATH="$PIPELINE:${PYTHONPATH:-}" python3 "$PIPELINE/graft_stage_a_outputs.py" \
    --input "$PIPELINE/inputs/$INPUT_PDB_NAME" \
    --rfd3-dir "$PIPELINE/outputs/rfd3_stage_a" \
    --out-dir "$PIPELINE/outputs/stage_a_grafted"

echo "Grafted: $(ls $PIPELINE/outputs/stage_a_grafted/*_grafted.pdb 2>/dev/null | wc -l) structures"

echo ""
echo "=== Step 3: QC — helicity, bridging, Fab fidelity ==="
PYTHONPATH="$PIPELINE:${PYTHONPATH:-}" python3 "$PIPELINE/analyze_stage_a_minibinders.py" \
    --input "$PIPELINE/inputs/$INPUT_PDB_NAME" \
    --grafted-dir "$PIPELINE/outputs/stage_a_grafted" \
    --out-csv "$PIPELINE/final/stage_a_minibinder_qc.csv" \
    || echo "  QC step failed (non-fatal)"

echo ""
echo "===== Stage A RFd3 complete: $(date) ====="
echo "Outputs: outputs/stage_a_grafted/*_grafted.pdb (input Fab unchanged + chain M)"
