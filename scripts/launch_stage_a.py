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

from binary_antibodies.lambda_client import LambdaClient, resolve_lambda_api_key

REMOTE_USER = "ubuntu"
REMOTE_PIPELINE = "/home/ubuntu/pipeline"
REMOTE_TMUX = "tmux"
GITHUB_BRANCH = "cursor/conditional-nanobody-design-992c"
INPUT_PDB = "fab_hidden_minibinder_stage_a.pdb"
CONFIG_JSON = "stage_a_hidden_minibinder_config.json"
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
    return str(EPHEMERAL_KEY_PATH)


def resolve_ssh_key(client: LambdaClient, preferred: str) -> tuple[str, list[str]]:
    """Return (private_key_path, ssh_key_names_for_launch)."""
    if ssh_key_usable(preferred):
        return os.path.expanduser(preferred), [os.environ.get("LAMBDA_SSH_KEY_NAME", "cursor-agent")]
    print(f"  SSH key not found at {preferred} — using ephemeral key.")
    key_path = ensure_ephemeral_ssh_key(client)
    return key_path, ["cursor-agent", EPHEMERAL_KEY_NAME]


def upload_stage_a(ip: str, key: str) -> None:
    print("  Uploading Stage A pipeline files ...")
    ssh(ip, key, f"mkdir -p {REMOTE_PIPELINE}/inputs")

    files = [
        (ROOT / "structures/domains" / INPUT_PDB, f"{REMOTE_PIPELINE}/inputs/{INPUT_PDB}"),
        (ROOT / "structures/interface" / CONFIG_JSON, f"{REMOTE_PIPELINE}/inputs/{CONFIG_JSON}"),
        (ROOT / "scripts/gpu_setup/run_stage_a_hidden_minibinder.sh", f"{REMOTE_PIPELINE}/run_stage_a.sh"),
        (ROOT / "scripts/gpu_setup/setup_pipeline_rfd3.sh", f"{REMOTE_PIPELINE}/setup_pipeline_rfd3.sh"),
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


def run_stage_a(ip: str, key: str, n_designs: int) -> None:
    env = " ".join([
        f"N_DESIGNS={n_designs}",
        f"INPUT_PDB={INPUT_PDB}",
        f"CONFIG_JSON={CONFIG_JSON}",
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

    subprocess.run([sys.executable, str(ROOT / "scripts/build_stage_a_design_target.py")], check=True)

    ssh_key, launch_key_names = resolve_ssh_key(client, args.ssh_key)
    if ssh_key_usable(args.ssh_key):
        launch_key_names = [args.ssh_key_name]

    inst = client.launch(
        "gpu_1x_a100_sxm4",
        "us-east-1",
        ssh_key_names=launch_key_names,
        name="stage-a-hidden-minibinder",
    )
    iid = inst["id"]
    ip = ""
    try:
        active = client.wait_until_active(iid)
        ip = active["ip"]
        print(f"\nInstance {iid} @ {ip}")
        print(f"  Branch: {GITHUB_BRANCH}")
        print(f"  Designs: {args.n_designs}")

        wait_ssh(ip, ssh_key)
        upload_stage_a(ip, ssh_key)
        run_setup(ip, ssh_key)
        run_stage_a(ip, ssh_key, args.n_designs)
        print_monitor_commands(ip, ssh_key)

        if not args.no_wait:
            if wait_for_pipeline(ip, ssh_key):
                ssh(ip, ssh_key, f"ls {REMOTE_PIPELINE}/outputs/rfd3_stage_a/*.cif 2>/dev/null | wc -l", check=False)
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
