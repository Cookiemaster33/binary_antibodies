#!/usr/bin/env python3
"""
launch_stage_a.py
-----------------
Launch Lambda GPU for Stage A hidden minibinder RFd3 design.

Uses SSH deploy (same pattern as launch_full_pipeline.py) because cloud-init
bootstrap is unreliable on Lambda (missing BioPython, runcmd timing).

Usage
-----
    python scripts/build_stage_a_design_target.py   # build PDB + config first
    python scripts/launch_stage_a.py
    python scripts/launch_stage_a.py --status
    python scripts/launch_stage_a.py --terminate --instance-id <id>
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from binary_antibodies.fab_hidden_switch import (
    COMPACT_MB_LENGTH_RANGE,
    DEFAULT_MB_LENGTH_RANGE,
    STAGE_A_SAMPLER,
)
from binary_antibodies.lambda_client import LambdaClient, resolve_lambda_api_key

REMOTE_USER = "ubuntu"
REMOTE_PIPELINE = "/home/ubuntu/pipeline"
REMOTE_TMUX = "tmux"
GITHUB_BRANCH = "cursor/helical-minibinder-fixed-fab-992c"
INPUT_PDB = "fab_hidden_minibinder_stage_a.pdb"
CONFIG_JSON = "stage_a_hidden_minibinder_config.json"
DEFAULT_STAGE0_FAB = (
    ROOT
    / "pipeline_results/stage_0_full_fab_fused_t025/structures/top5_holo"
    / "rank079_s0_native_split_s296_model_0_split.cif"
)
DEFAULT_STAGE0_HOLO = (
    ROOT
    / "pipeline_results/stage_0_full_fab_fused_t025/structures/holo"
    / "rank079_s0_native_split_s296_model_0.cif"
)
EPHEMERAL_KEY_NAME = "cursor-cloud-ephemeral-992c"
EPHEMERAL_KEY_PATH = Path.home() / ".ssh" / "cursor_lambda_ephemeral"


def ssh_opts(key: str) -> list[str]:
    return [
        "-i", os.path.expanduser(key),
        "-o", "StrictHostKeyChecking=no",
        "-o", "ConnectTimeout=15",
        "-o", "ServerAliveInterval=30",
    ]


def ssh_key_usable(key: str) -> bool:
    path = os.path.expanduser(key)
    return os.path.isfile(path) and os.access(path, os.R_OK)


def ssh(ip: str, key: str, cmd: str, check: bool = True) -> int:
    r = subprocess.run(
        ["ssh"] + ssh_opts(key) + [f"{REMOTE_USER}@{ip}", cmd],
        check=False,
    )
    if check and r.returncode != 0:
        raise RuntimeError(f"SSH failed ({r.returncode}): {cmd[:120]}")
    return r.returncode


def scp(local: str, ip: str, remote: str, key: str) -> None:
    subprocess.run(
        ["scp"] + ssh_opts(key) + [local, f"{REMOTE_USER}@{ip}:{remote}"],
        check=True,
    )


def scp_from(ip: str, remote: str, local: str, key: str, recursive: bool = False) -> int:
    cmd = ["scp"] + ssh_opts(key)
    if recursive:
        cmd.append("-r")
    cmd += [f"{REMOTE_USER}@{ip}:{remote}", local]
    return subprocess.run(cmd, check=False).returncode


def wait_ssh(ip: str, key: str, timeout_s: int = 600) -> None:
    print(f"  Waiting for SSH on {ip}", end="", flush=True)
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        if ssh(ip, key, "echo ok", check=False) == 0:
            print(" ready!")
            return
        print(".", end="", flush=True)
        time.sleep(10)
    raise TimeoutError(f"SSH not ready on {ip}")


def ensure_ephemeral_ssh_key(client: LambdaClient) -> str:
    """Create a local ephemeral key and register it with Lambda if needed."""
    EPHEMERAL_KEY_PATH.parent.mkdir(parents=True, exist_ok=True)
    if not EPHEMERAL_KEY_PATH.exists():
        subprocess.run(
            ["ssh-keygen", "-t", "ed25519", "-f", str(EPHEMERAL_KEY_PATH), "-N", "", "-C", EPHEMERAL_KEY_NAME],
            check=True,
        )
    pub = EPHEMERAL_KEY_PATH.with_suffix(".pub").read_text().strip()
    existing = {k.get("name"): k for k in client.list_ssh_keys()}
    if EPHEMERAL_KEY_NAME not in existing:
        print(f"  Registering ephemeral SSH key: {EPHEMERAL_KEY_NAME}")
        client.add_ssh_key(EPHEMERAL_KEY_NAME, pub)
    elif existing[EPHEMERAL_KEY_NAME].get("public_key", "").strip() != pub:
        print(f"  Updating ephemeral SSH key: {EPHEMERAL_KEY_NAME}")
        client.delete_ssh_key(EPHEMERAL_KEY_NAME)
        client.add_ssh_key(EPHEMERAL_KEY_NAME, pub)
    return str(EPHEMERAL_KEY_PATH)


MATERIALISED_KEY_PATH = Path.home() / ".ssh" / "lambda_agent_key_from_env"


def looks_like_private_key(value: str) -> bool:
    return "PRIVATE KEY" in value


def materialise_private_key(material: str) -> str:
    """Write inline private-key material to disk and return its path.

    ``LAMBDA_SSH_KEY`` sometimes carries the key itself rather than a path (this
    is how secrets are injected into Cloud Agent VMs).
    """
    MATERIALISED_KEY_PATH.parent.mkdir(parents=True, exist_ok=True)
    text = material if material.endswith("\n") else material + "\n"
    MATERIALISED_KEY_PATH.write_text(text)
    MATERIALISED_KEY_PATH.chmod(0o600)
    return str(MATERIALISED_KEY_PATH)


def resolve_ssh_key(client: LambdaClient, preferred: str, ssh_key_name: str) -> tuple[str, str]:
    """Return (private_key_path, ssh_key_name_for_launch). Lambda allows one key only.

    Never print the key itself: the value may be the private key material.
    """
    if looks_like_private_key(preferred):
        registered = {k.get("name") for k in client.list_ssh_keys()}
        if ssh_key_name in registered:
            print(f"  Using inline SSH key material for registered key '{ssh_key_name}'.")
            return materialise_private_key(preferred), ssh_key_name
        print(
            f"  Inline SSH key material provided but '{ssh_key_name}' is not registered "
            "with Lambda — using ephemeral key."
        )
        return ensure_ephemeral_ssh_key(client), EPHEMERAL_KEY_NAME

    preferred_path = os.path.expanduser(preferred)
    if os.path.abspath(preferred_path) == os.path.abspath(str(EPHEMERAL_KEY_PATH)):
        return str(EPHEMERAL_KEY_PATH), EPHEMERAL_KEY_NAME
    if ssh_key_usable(preferred):
        return preferred_path, ssh_key_name
    print(f"  SSH key not found at {preferred_path} — using ephemeral key.")
    key_path = ensure_ephemeral_ssh_key(client)
    return key_path, EPHEMERAL_KEY_NAME


def upload_stage_a(ip: str, key: str) -> None:
    print("  Uploading Stage A pipeline files ...")
    ssh(ip, key, f"mkdir -p {REMOTE_PIPELINE}/inputs {REMOTE_PIPELINE}/binary_antibodies")

    files = [
        (ROOT / "structures/domains" / INPUT_PDB, f"{REMOTE_PIPELINE}/inputs/{INPUT_PDB}"),
        (ROOT / "structures/interface" / CONFIG_JSON, f"{REMOTE_PIPELINE}/inputs/{CONFIG_JSON}"),
        (ROOT / "scripts/gpu_setup/run_stage_a_hidden_minibinder.sh", f"{REMOTE_PIPELINE}/run_stage_a.sh"),
        (ROOT / "scripts/gpu_setup/setup_pipeline_rfd3.sh", f"{REMOTE_PIPELINE}/setup_pipeline_rfd3.sh"),
        (ROOT / "scripts/graft_stage_a_outputs.py", f"{REMOTE_PIPELINE}/graft_stage_a_outputs.py"),
        (
            ROOT / "scripts/analyze_stage_a_minibinders.py",
            f"{REMOTE_PIPELINE}/analyze_stage_a_minibinders.py",
        ),
        (
            ROOT / "scripts/gpu_setup/remote_package_init.py",
            f"{REMOTE_PIPELINE}/binary_antibodies/__init__.py",
        ),
        (ROOT / "binary_antibodies/fab_hidden_switch.py", f"{REMOTE_PIPELINE}/binary_antibodies/fab_hidden_switch.py"),
        (ROOT / "binary_antibodies/backbone_ss.py", f"{REMOTE_PIPELINE}/binary_antibodies/backbone_ss.py"),
    ]
    for local, remote in files:
        if not local.exists():
            raise FileNotFoundError(f"Missing local file: {local}")
        scp(str(local), ip, remote, key)

    ssh(ip, key, f"chmod +x {REMOTE_PIPELINE}/run_stage_a.sh {REMOTE_PIPELINE}/setup_pipeline_rfd3.sh")


def run_setup(ip: str, key: str) -> None:
    print("  Running foundry Docker setup ...")
    ssh(
        ip,
        key,
        f"{REMOTE_TMUX} new-session -d -s setup -c {REMOTE_PIPELINE} "
        f"'bash {REMOTE_PIPELINE}/setup_pipeline_rfd3.sh 2>&1 | tee {REMOTE_PIPELINE}/setup.log'",
    )
    while True:
        if ssh(ip, key, f"grep -q 'Setup complete' {REMOTE_PIPELINE}/setup.log", check=False) == 0:
            print("  Setup complete.")
            return
        if ssh(ip, key, f"{REMOTE_TMUX} has-session -t setup", check=False) != 0:
            ssh(ip, key, f"tail -30 {REMOTE_PIPELINE}/setup.log", check=False)
            raise RuntimeError("Setup session ended unexpectedly")
        time.sleep(20)


def run_stage_a(
    ip: str,
    key: str,
    n_designs: int,
    mb_length_ranges: str,
    is_non_loopy: bool,
    step_scale: float,
    gamma_0: float,
    batch_size: int,
    low_memory_mode: bool,
) -> None:
    env = " ".join([
        f"N_DESIGNS={n_designs}",
        f"BATCH_SIZE={batch_size}",
        f"LOW_MEMORY_MODE={1 if low_memory_mode else 0}",
        f"INPUT_PDB={INPUT_PDB}",
        f"CONFIG_JSON={CONFIG_JSON}",
        f"MB_LENGTH_RANGES={mb_length_ranges}",
        f"IS_NON_LOOPY={1 if is_non_loopy else 0}",
        f"STEP_SCALE={step_scale}",
        f"GAMMA_0={gamma_0}",
    ])
    print("  Starting Stage A RFd3 in tmux session 'stage_a' ...")
    ssh(
        ip,
        key,
        f"{REMOTE_TMUX} new-session -d -s stage_a -c {REMOTE_PIPELINE} "
        f"\"{env} bash {REMOTE_PIPELINE}/run_stage_a.sh\"",
    )


def wait_for_pipeline(ip: str, key: str, timeout_s: int = 14400) -> bool:
    print(f"  Polling pipeline log on {ip} (up to {timeout_s // 3600}h)...")
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        if ssh(ip, key, f"grep -q 'Stage A RFd3 complete' {REMOTE_PIPELINE}/stage_a_pipeline.log", check=False) == 0:
            print("  Stage A RFd3 complete!")
            return True
        if ssh(ip, key, f"{REMOTE_TMUX} has-session -t stage_a", check=False) != 0:
            ssh(ip, key, f"tail -40 {REMOTE_PIPELINE}/stage_a_pipeline.log", check=False)
            raise RuntimeError("Stage A tmux session ended before completion")
        time.sleep(120)
    return False


def fetch_results(ip: str, key: str, results_dir: Path, fab_pdb: Path) -> None:
    """Pull grafted structures, raw RFd3 CIFs, QC table and logs into results_dir."""
    structures = results_dir / "structures"
    inputs = results_dir / "inputs"
    final = results_dir / "final"
    for d in (structures, inputs, final):
        d.mkdir(parents=True, exist_ok=True)

    print(f"  Fetching results → {results_dir}")
    ssh(
        ip,
        key,
        f"cd {REMOTE_PIPELINE}/outputs && tar czf {REMOTE_PIPELINE}/grafted_pdbs.tar.gz "
        f"-C stage_a_grafted . && tar czf {REMOTE_PIPELINE}/rfd3_stage_a_cifs.tar.gz "
        f"-C rfd3_stage_a .",
        check=False,
    )
    scp_from(ip, f"{REMOTE_PIPELINE}/grafted_pdbs.tar.gz", str(structures), key)
    scp_from(ip, f"{REMOTE_PIPELINE}/rfd3_stage_a_cifs.tar.gz", str(structures), key)
    scp_from(ip, f"{REMOTE_PIPELINE}/final/stage_a_minibinder_qc.csv", str(final), key)
    scp_from(ip, f"{REMOTE_PIPELINE}/stage_a_pipeline.log", str(final), key)
    scp_from(ip, f"{REMOTE_PIPELINE}/inputs/{INPUT_PDB}", str(inputs), key)
    scp_from(ip, f"{REMOTE_PIPELINE}/inputs/{CONFIG_JSON}", str(inputs), key)

    samples = structures / "grafted"
    samples.mkdir(parents=True, exist_ok=True)
    ssh(
        ip,
        key,
        f"mkdir -p {REMOTE_PIPELINE}/sample_out && "
        f"ls {REMOTE_PIPELINE}/outputs/stage_a_grafted/*_grafted.pdb | head -5 | "
        f"xargs -I{{}} cp {{}} {REMOTE_PIPELINE}/sample_out/",
        check=False,
    )
    scp_from(ip, f"{REMOTE_PIPELINE}/sample_out/*", str(samples), key)
    print(f"  Source Fab was: {fab_pdb}")


# Preference order for RFd3: single-GPU 40-80 GB cards first, then smaller cards,
# then multi-GPU boxes (only one GPU is used, so they are a last resort on cost).
INSTANCE_PREFERENCE = [
    "gpu_1x_a100_sxm4",
    "gpu_1x_a100",
    "gpu_1x_a100_80gb_sxm4",
    "gpu_1x_h100_pcie",
    "gpu_1x_h100_sxm5",
    "gpu_1x_gh200",
    "gpu_1x_a6000",
    "gpu_1x_a10",
    "gpu_2x_a100",
    "gpu_4x_a100",
    "gpu_8x_a100_80gb_sxm4",
    "gpu_8x_a100",
]

# RFd3 on the ~530-residue Fab + minibinder system. Smaller cards need a smaller
# diffusion batch to stay inside VRAM.
BATCH_SIZE_BY_TYPE = {
    "gpu_1x_a10": 2,
    "gpu_1x_a6000": 4,
}
DEFAULT_BATCH_SIZE = 10


def select_instance_type(
    client: LambdaClient, requested: str | None, max_price: float | None
) -> tuple[str, str, float]:
    """Pick (instance_type, region, price_per_hour) from what currently has capacity."""
    available = {t["name"]: t for t in client.available_instance_types()}
    if not available:
        raise RuntimeError("No Lambda instance types currently have capacity")

    if requested:
        if requested not in available:
            raise RuntimeError(
                f"{requested} has no capacity. Available: "
                + ", ".join(f"{n} (${available[n]['price_per_hour']:.2f}/h)" for n in available)
            )
        chosen = requested
    else:
        ordered = [n for n in INSTANCE_PREFERENCE if n in available]
        ordered += sorted(n for n in available if n not in INSTANCE_PREFERENCE)
        if max_price is not None:
            affordable = [n for n in ordered if available[n]["price_per_hour"] <= max_price]
            if not affordable:
                raise RuntimeError(
                    f"No instance type with capacity is within ${max_price:.2f}/h. Available: "
                    + ", ".join(
                        f"{n} (${available[n]['price_per_hour']:.2f}/h)" for n in ordered
                    )
                )
            ordered = affordable
        chosen = ordered[0]

    info = available[chosen]
    return chosen, info["available_regions"][0], info["price_per_hour"]


def print_monitor_commands(ip: str, key: str) -> None:
    print("\nMonitor:")
    print(f"  ssh -i {key} {REMOTE_USER}@{ip} 'tail -f {REMOTE_PIPELINE}/stage_a_pipeline.log'")
    print(f"  ssh -i {key} {REMOTE_USER}@{ip} '{REMOTE_TMUX} attach -t stage_a'")
    print(f"  ssh -i {key} {REMOTE_USER}@{ip} 'nvidia-smi'")


def main() -> None:
    p = argparse.ArgumentParser(description="Launch Stage A hidden minibinder RFd3 on Lambda.")
    p.add_argument("--ssh-key", default=os.environ.get("LAMBDA_SSH_KEY", "~/.ssh/lambda_agent_key"))
    p.add_argument("--ssh-key-name", default="cursor-agent")
    p.add_argument("--n-designs", type=int, default=200)
    p.add_argument(
        "--fab-pdb",
        type=Path,
        default=DEFAULT_STAGE0_FAB,
        help="Placed/split Fab (chains A–D[,T]) or fused holo for build_stage_a_design_target.py. "
        "If the arms are already separated the coordinates are used verbatim.",
    )
    p.add_argument(
        "--mb-length-ranges",
        default=f"{DEFAULT_MB_LENGTH_RANGE},{COMPACT_MB_LENGTH_RANGE}",
        help="Comma-separated minibinder length windows to sample "
        f"(default '{DEFAULT_MB_LENGTH_RANGE},{COMPACT_MB_LENGTH_RANGE}')",
    )
    p.add_argument(
        "--no-non-loopy",
        action="store_true",
        help="Disable is_non_loopy (the helix bias); off-target designs become loopier",
    )
    p.add_argument("--step-scale", type=float, default=STAGE_A_SAMPLER["step_scale"])
    p.add_argument("--gamma-0", type=float, default=STAGE_A_SAMPLER["gamma_0"])
    p.add_argument(
        "--instance-type",
        default="",
        help="Force a Lambda instance type; default picks the best one with capacity",
    )
    p.add_argument(
        "--max-price",
        type=float,
        default=None,
        help="Skip instance types above this $/hour when auto-selecting",
    )
    p.add_argument(
        "--wait-for-capacity-min",
        type=int,
        default=0,
        help="Poll for a suitable instance type for this many minutes before giving up",
    )
    p.add_argument(
        "--attach-instance-id",
        default="",
        help="Reuse an instance that is already launched instead of starting a new one",
    )
    p.add_argument(
        "--boot-timeout-s",
        type=int,
        default=1800,
        help="How long to wait for the instance to become active (default 1800)",
    )
    p.add_argument(
        "--batch-size",
        type=int,
        default=None,
        help="RFd3 diffusion batch size (default depends on the GPU)",
    )
    p.add_argument(
        "--results-dir",
        type=Path,
        default=None,
        help="Local directory to pull results into when the run completes",
    )
    p.add_argument("--no-terminate", action="store_true")
    p.add_argument("--no-wait", action="store_true", help="Start pipeline and exit without waiting for completion.")
    p.add_argument("--status", action="store_true")
    p.add_argument("--terminate", action="store_true")
    p.add_argument("--instance-id", default="")
    args = p.parse_args()

    api_key = resolve_lambda_api_key()
    if not api_key:
        sys.exit("ERROR: Lambda API key not found")

    client = LambdaClient(api_key=api_key)

    if args.status:
        for inst in client.list_instances():
            print(
                f"  {inst['id']}  {inst.get('name', '?')}  "
                f"status={inst.get('status')}  ip={inst.get('ip')}"
            )
        return

    if args.terminate:
        if not args.instance_id:
            sys.exit("ERROR: --terminate requires --instance-id")
        client.terminate(args.instance_id)
        return

    fab_pdb = args.fab_pdb
    if not fab_pdb.exists() and fab_pdb == DEFAULT_STAGE0_FAB and DEFAULT_STAGE0_HOLO.exists():
        print(f"  Split CIF missing; building from holo → {fab_pdb}")
        fab_pdb.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/build_stage0_holo_split_cif.py"),
                str(DEFAULT_STAGE0_HOLO),
                "-o",
                str(fab_pdb),
            ],
            check=True,
        )
    if not fab_pdb.exists():
        sys.exit(f"ERROR: Stage 0 Fab not found: {fab_pdb}")
    print(f"  Stage 0 Fab source: {fab_pdb}")
    primary_range = args.mb_length_ranges.split(",")[0].strip()
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/build_stage_a_design_target.py"),
            "--fab-pdb",
            str(fab_pdb),
            "--mb-length-range",
            primary_range,
        ],
        check=True,
    )

    ssh_key, launch_key_name = resolve_ssh_key(client, args.ssh_key, args.ssh_key_name)

    if args.attach_instance_id:
        iid = args.attach_instance_id
        existing = client.get_instance(iid)
        instance_type = existing.get("instance_type", {}).get("name", "unknown")
        print(f"  Attaching to existing instance {iid} ({instance_type})")
        # The instance only accepts the key it was launched with.
        instance_keys = existing.get("ssh_key_names") or []
        if instance_keys and launch_key_name not in instance_keys:
            if EPHEMERAL_KEY_NAME in instance_keys and EPHEMERAL_KEY_PATH.exists():
                print(f"  Switching to '{EPHEMERAL_KEY_NAME}' to match the instance's key.")
                ssh_key, launch_key_name = str(EPHEMERAL_KEY_PATH), EPHEMERAL_KEY_NAME
            else:
                sys.exit(
                    f"ERROR: instance {iid} accepts SSH keys {instance_keys}, but the "
                    f"resolved key is '{launch_key_name}' and no matching private key "
                    "is available locally."
                )
    else:
        deadline = time.time() + args.wait_for_capacity_min * 60
        while True:
            try:
                instance_type, region, price = select_instance_type(
                    client, args.instance_type or None, args.max_price
                )
                break
            except RuntimeError as exc:
                if time.time() >= deadline:
                    sys.exit(f"ERROR: {exc}")
                print(f"  {exc}\n  Waiting for capacity ...")
                time.sleep(120)
        print(f"  Instance: {instance_type} in {region} (${price:.2f}/h)")
        inst = client.launch(
            instance_type,
            region,
            ssh_key_names=[launch_key_name],
            name="stage-a-hidden-minibinder",
        )
        iid = inst["id"]

    batch_size = args.batch_size or BATCH_SIZE_BY_TYPE.get(instance_type, DEFAULT_BATCH_SIZE)
    low_memory = instance_type in BATCH_SIZE_BY_TYPE
    print(
        f"  Batch size {batch_size}" + (", low-memory mode" if low_memory else "")
    )
    ip = ""
    try:
        active = client.wait_until_active(iid, timeout_s=args.boot_timeout_s)
        ip = active["ip"]
        print(f"\nInstance {iid} @ {ip}")
        print(f"  Branch: {GITHUB_BRANCH}")
        print(f"  Designs: {args.n_designs}  length windows: {args.mb_length_ranges}")
        print(
            f"  Helical conditioning: is_non_loopy={not args.no_non_loopy}, "
            f"step_scale={args.step_scale}, gamma_0={args.gamma_0}"
        )

        wait_ssh(ip, ssh_key)
        upload_stage_a(ip, ssh_key)
        run_setup(ip, ssh_key)
        run_stage_a(
            ip,
            ssh_key,
            args.n_designs,
            args.mb_length_ranges,
            not args.no_non_loopy,
            args.step_scale,
            args.gamma_0,
            batch_size,
            low_memory,
        )
        print_monitor_commands(ip, ssh_key)

        if not args.no_wait:
            if wait_for_pipeline(ip, ssh_key):
                ssh(ip, ssh_key, f"ls {REMOTE_PIPELINE}/outputs/rfd3_stage_a/*.cif 2>/dev/null | wc -l", check=False)
                if args.results_dir:
                    fetch_results(ip, ssh_key, args.results_dir, fab_pdb)
                if not args.no_terminate:
                    print(f"\nTerminating {iid} ...")
                    client.terminate(iid)
            elif not args.no_terminate:
                print("  Pipeline still running — leaving instance up.")
        else:
            print(f"\nInstance left running: {iid}")

    except Exception:
        print(f"\nFailed — instance left running for debug: {iid} @ {ip}")
        if ip:
            print_monitor_commands(ip, ssh_key)
        raise


if __name__ == "__main__":
    main()
