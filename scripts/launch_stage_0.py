#!/usr/bin/env python3
"""
launch_stage_0.py
-----------------
Launch Lambda GPU for Stage 0 VH–VL interface weakening.

    python scripts/build_stage_0_design_target.py
    python scripts/launch_stage_0.py --no-wait --no-terminate
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

from binary_antibodies.lambda_client import LambdaClient, resolve_lambda_api_key

REMOTE_USER = "ubuntu"
REMOTE_PIPELINE = "/home/ubuntu/pipeline"
REMOTE_TMUX = "tmux"
INPUT_PDB = "fab_stage_0_vhvL_interface.pdb"
SPLIT_PDB = "fab_stage_0_split_mpnn.pdb"
CONFIG_JSON = "stage_0_vhvL_interface_config.json"
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
    return str(EPHEMERAL_KEY_PATH)


def resolve_ssh_key(client: LambdaClient, preferred: str, ssh_key_name: str) -> tuple[str, str]:
    if ssh_key_usable(preferred):
        return os.path.expanduser(preferred), ssh_key_name
    print(f"  SSH key not found at {preferred} — using ephemeral key.")
    return ensure_ephemeral_ssh_key(client), EPHEMERAL_KEY_NAME


def upload_stage_0(ip: str, key: str) -> None:
    print("  Uploading Stage 0 pipeline files ...")
    ssh(ip, key, f"mkdir -p {REMOTE_PIPELINE}/inputs {REMOTE_PIPELINE}/binary_antibodies")

    files = [
        (ROOT / "structures/domains" / INPUT_PDB, f"{REMOTE_PIPELINE}/inputs/{INPUT_PDB}"),
        (ROOT / "structures/domains" / SPLIT_PDB, f"{REMOTE_PIPELINE}/inputs/{SPLIT_PDB}"),
        (ROOT / "structures/interface" / CONFIG_JSON, f"{REMOTE_PIPELINE}/inputs/{CONFIG_JSON}"),
        (ROOT / "scripts/gpu_setup/run_stage_0_vhvL_interface.sh", f"{REMOTE_PIPELINE}/run_stage_0.sh"),
        (ROOT / "scripts/gpu_setup/setup_pipeline_rfd3.sh", f"{REMOTE_PIPELINE}/setup_pipeline_rfd3.sh"),
        (ROOT / "binary_antibodies/stage_0_scoring.py", f"{REMOTE_PIPELINE}/binary_antibodies/stage_0_scoring.py"),
        (ROOT / "binary_antibodies/fab_hidden_switch.py", f"{REMOTE_PIPELINE}/binary_antibodies/fab_hidden_switch.py"),
    ]
    for local, remote in files:
        if not local.exists():
            raise FileNotFoundError(f"Missing local file: {local}")
        scp(str(local), ip, remote, key)

    ssh(ip, key, f"chmod +x {REMOTE_PIPELINE}/run_stage_0.sh {REMOTE_PIPELINE}/setup_pipeline_rfd3.sh")


def run_setup(ip: str, key: str) -> None:
    print("  Running foundry Docker setup ...")
    ssh(
        ip, key,
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


def run_stage_0(ip: str, key: str, n_mpnn_seqs: int) -> None:
    env = " ".join([
        f"N_MPNN_SEQS={n_mpnn_seqs}",
        f"INPUT_PDB={INPUT_PDB}",
        f"SPLIT_PDB={SPLIT_PDB}",
        f"CONFIG_JSON={CONFIG_JSON}",
    ])
    print("  Starting Stage 0 pipeline in tmux session 'stage_0' ...")
    ssh(
        ip, key,
        f"{REMOTE_TMUX} new-session -d -s stage_0 -c {REMOTE_PIPELINE} "
        f"\"{env} bash {REMOTE_PIPELINE}/run_stage_0.sh\"",
    )


def pick_region(client: LambdaClient, instance_type: str, preferred: str = "us-east-1") -> str:
    """Choose a region with capacity for the requested instance type."""
    types = client.list_instance_types()
    info = types.get(instance_type, {})
    regions = [r["name"] for r in info.get("regions_with_capacity_available", [])]
    if preferred in regions:
        return preferred
    if regions:
        print(f"  No {instance_type} capacity in {preferred}; using {regions[0]}")
        return regions[0]
    available = client.available_instance_types()
    names = ", ".join(t["name"] for t in available[:5])
    raise RuntimeError(
        f"No Lambda capacity for {instance_type}. Available types: {names}"
    )


def main() -> None:
    p = argparse.ArgumentParser(description="Launch Stage 0 VH–VL interface design on Lambda.")
    p.add_argument("--ssh-key", default=os.environ.get("LAMBDA_SSH_KEY", "~/.ssh/lambda_agent_key"))
    p.add_argument("--ssh-key-name", default="cursor-agent")
    p.add_argument("--n-mpnn-seqs", type=int, default=32, help="ProteinMPNN sequences per backbone")
    p.add_argument("--region", default="us-east-1", help="Preferred Lambda region")
    p.add_argument("--instance-type", default="gpu_1x_a100_sxm4", help="Lambda instance type")
    p.add_argument("--no-terminate", action="store_true")
    p.add_argument("--no-wait", action="store_true")
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
            print(f"  {inst['id']}  {inst.get('name','?')}  status={inst.get('status')}  ip={inst.get('ip')}")
        return

    if args.terminate:
        if not args.instance_id:
            sys.exit("ERROR: --terminate requires --instance-id")
        client.terminate(args.instance_id)
        return

    subprocess.run([sys.executable, str(ROOT / "scripts/build_stage_0_design_target.py")], check=True)

    ssh_key, launch_key_name = resolve_ssh_key(client, args.ssh_key, args.ssh_key_name)
    region = pick_region(client, args.instance_type, args.region)
    inst = client.launch(
        args.instance_type,
        region,
        ssh_key_names=[launch_key_name],
        name="stage-0-vhvL-interface",
    )
    iid = inst["id"]
    ip = ""
    try:
        active = client.wait_until_active(iid)
        ip = active["ip"]
        print(f"\nInstance {iid} @ {ip}")
        print(f"  Region: {region} | type: {args.instance_type}")
        print(f"  MPNN sequences: {args.n_mpnn_seqs}")

        wait_ssh(ip, ssh_key)
        upload_stage_0(ip, ssh_key)
        run_setup(ip, ssh_key)
        run_stage_0(ip, ssh_key, args.n_mpnn_seqs)

        print("\nMonitor:")
        print(f"  ssh -i {ssh_key} {REMOTE_USER}@{ip} 'tail -f {REMOTE_PIPELINE}/stage_0_pipeline.log'")
        print(f"  ssh -i {ssh_key} {REMOTE_USER}@{ip} '{REMOTE_TMUX} attach -t stage_0'")
        print(f"\nInstance left running: {iid}")
    except Exception:
        print(f"\nFailed — instance left running: {iid} @ {ip}")
        raise


if __name__ == "__main__":
    main()
