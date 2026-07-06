#!/usr/bin/env python3
"""
launch_full_pipeline.py
-----------------------
Launch Lambda Cloud GPU, run run_full_pipeline.sh (RFd3 → MPNN → Boltz → scRMSD),
push results to a new pipeline_results/ subfolder on GitHub.

Usage
-----
    export LAMBDA_API_KEY=...
    python scripts/launch_full_pipeline.py \\
        --ssh-key ~/.ssh/lambda_agent_key \\
        --rfd3-rounds 2 \\
        --results-dir pipeline_results/v5_two_round_refine

    # Monitor a running job
    python scripts/launch_full_pipeline.py --status

    # Terminate a stuck instance
    python scripts/launch_full_pipeline.py --terminate --instance-id <id>
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

from binary_antibodies.lambda_client import LambdaClient

DEFAULT_INSTANCE_TYPE = "gpu_1x_a100_sxm4"
DEFAULT_REGION = "us-east-1"
DEFAULT_SSH_KEY_NAME = os.environ.get("LAMBDA_SSH_KEY_NAME", "cursor-agent-binary-antibodies")
REMOTE_USER = "ubuntu"
REMOTE_PIPELINE = "/home/ubuntu/pipeline"
INPUT_PDB = "structures/domains/vh1_vl_nanobody_design_target.pdb"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Launch full RFd3/MPNN/Boltz pipeline on Lambda.")
    p.add_argument("--ssh-key", default=os.environ.get("LAMBDA_SSH_KEY", "~/.ssh/lambda_agent_key"))
    p.add_argument("--ssh-key-name", default=DEFAULT_SSH_KEY_NAME)
    p.add_argument("--instance-type", default=DEFAULT_INSTANCE_TYPE)
    p.add_argument("--region", default=DEFAULT_REGION)
    p.add_argument("--api-key", default=os.environ.get("LAMBDA_API_KEY", ""))
    p.add_argument("--n-designs", type=int, default=200)
    p.add_argument("--n-mpnn-seqs", type=int, default=8)
    p.add_argument("--top-n", type=int, default=50)
    p.add_argument("--mb-length-range", default="35-70")
    p.add_argument("--rfd3-rounds", type=int, default=2, choices=[1, 2])
    p.add_argument("--rfd3-round2-templates", type=int, default=5)
    p.add_argument("--rfd3-round2-designs-per-template", type=int, default=8)
    p.add_argument("--rfd3-partial-t", type=float, default=2.0)
    p.add_argument(
        "--results-dir",
        default="pipeline_results/v5_two_round_refine",
        help="GitHub subfolder for results (under repo root).",
    )
    p.add_argument("--github-branch", default="cursor/conditional-nanobody-design-992c")
    p.add_argument("--no-terminate", action="store_true")
    p.add_argument("--status", action="store_true")
    p.add_argument("--terminate", action="store_true")
    p.add_argument("--instance-id", default="")
    return p.parse_args()


def ssh_opts(key: str) -> list[str]:
    return [
        "-i", os.path.expanduser(key),
        "-o", "StrictHostKeyChecking=no",
        "-o", "ConnectTimeout=15",
        "-o", "ServerAliveInterval=30",
    ]


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


def github_token() -> str:
    tok = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if tok:
        return tok
    try:
        return subprocess.check_output(["gh", "auth", "token"], text=True).strip()
    except Exception:
        return ""


def upload_pipeline(ip: str, key: str) -> None:
    print("  Uploading pipeline files ...")
    ssh(ip, key, f"mkdir -p {REMOTE_PIPELINE}/inputs {REMOTE_PIPELINE}/binary_antibodies")

    files = [
        (ROOT / INPUT_PDB, f"{REMOTE_PIPELINE}/inputs/vh1_vl_nanobody_design_target.pdb"),
        (ROOT / "scripts/gpu_setup/run_full_pipeline.sh", f"{REMOTE_PIPELINE}/run_full_pipeline.sh"),
        (ROOT / "scripts/gpu_setup/rfd3_round2.py", f"{REMOTE_PIPELINE}/rfd3_round2.py"),
        (ROOT / "scripts/gpu_setup/push_results.sh", f"{REMOTE_PIPELINE}/push_results.sh"),
        (ROOT / "scripts/gpu_setup/setup_pipeline_rfd3.sh", f"{REMOTE_PIPELINE}/setup_pipeline_rfd3.sh"),
        (ROOT / "binary_antibodies/scoring.py", f"{REMOTE_PIPELINE}/binary_antibodies/scoring.py"),
    ]
    for local, remote in files:
        scp(str(local), ip, remote, key)

    ssh(ip, key, f"chmod +x {REMOTE_PIPELINE}/run_full_pipeline.sh {REMOTE_PIPELINE}/push_results.sh {REMOTE_PIPELINE}/setup_pipeline_rfd3.sh")


def run_setup(ip: str, key: str) -> None:
    print("  Running foundry Docker setup ...")
    ssh(ip, key,
        f"tmux -f /exec-daemon/tmux.portal.conf new-session -d -s setup -c {REMOTE_PIPELINE} "
        f"'bash {REMOTE_PIPELINE}/setup_pipeline_rfd3.sh 2>&1 | tee setup.log'")
    while True:
        if ssh(ip, key, f"grep -q 'Setup complete' {REMOTE_PIPELINE}/setup.log", check=False) == 0:
            print("  Setup complete.")
            return
        if ssh(ip, key, "tmux -f /exec-daemon/tmux.portal.conf has-session -t setup", check=False) != 0:
            ssh(ip, key, f"tail -30 {REMOTE_PIPELINE}/setup.log", check=False)
            raise RuntimeError("Setup session ended unexpectedly")
        time.sleep(20)


def run_pipeline(ip: str, key: str, args: argparse.Namespace, gh_tok: str) -> None:
    env = " ".join([
        f"N_DESIGNS={args.n_designs}",
        f"N_MPNN_SEQS={args.n_mpnn_seqs}",
        f"TOP_N={args.top_n}",
        f"MB_LENGTH_RANGE={args.mb_length_range}",
        f"RFD3_ROUNDS={args.rfd3_rounds}",
        f"RFD3_ROUND2_TEMPLATES={args.rfd3_round2_templates}",
        f"RFD3_ROUND2_DESIGNS_PER_TEMPLATE={args.rfd3_round2_designs_per_template}",
        f"RFD3_PARTIAL_T={args.rfd3_partial_t}",
        f"RESULTS_GITHUB_DIR={args.results_dir}",
        f"GITHUB_BRANCH={args.github_branch}",
        f"GITHUB_TOKEN={gh_tok}",
    ])
    print(f"  Starting pipeline (RFD3_ROUNDS={args.rfd3_rounds}) ...")
    ssh(ip, key,
        f"tmux -f /exec-daemon/tmux.portal.conf new-session -d -s pipeline -c {REMOTE_PIPELINE} "
        f"\"{env} bash {REMOTE_PIPELINE}/run_full_pipeline.sh\"")
    print(f"  Monitor: ssh -i {args.ssh_key} {REMOTE_USER}@{ip} "
          f"'tmux -f /exec-daemon/tmux.portal.conf attach -t pipeline'")

    while True:
        if ssh(ip, key, f"grep -q 'Complete:' {REMOTE_PIPELINE}/full_pipeline.log", check=False) == 0:
            print("  Pipeline complete!")
            return
        if ssh(ip, key, "tmux -f /exec-daemon/tmux.portal.conf has-session -t pipeline", check=False) != 0:
            ssh(ip, key, f"tail -40 {REMOTE_PIPELINE}/full_pipeline.log", check=False)
            raise RuntimeError("Pipeline session ended before completion")
        rc = ssh(ip, key,
                 f"grep -iE 'error|traceback|failed' {REMOTE_PIPELINE}/full_pipeline.log | tail -3",
                 check=False)
        time.sleep(60)


def main() -> None:
    args = parse_args()
    if not args.api_key:
        print("ERROR: Set LAMBDA_API_KEY in Cursor secrets or pass --api-key.")
        sys.exit(1)

    client = LambdaClient(api_key=args.api_key)

    if args.status:
        for inst in client.list_instances():
            print(f"  {inst['id']}  {inst.get('name','?')}  status={inst.get('status')}  ip={inst.get('ip')}")
        return

    if args.terminate:
        if not args.instance_id:
            print("ERROR: --terminate requires --instance-id")
            sys.exit(1)
        client.terminate(args.instance_id)
        print(f"Terminated {args.instance_id}")
        return

    gh_tok = github_token()
    if not gh_tok:
        print("WARNING: No GitHub token — results will not be pushed from instance.")

    avail = {t["name"]: t for t in client.available_instance_types()}
    region = args.region
    if args.instance_type not in avail:
        print(f"ERROR: {args.instance_type} not available.")
        sys.exit(1)
    if region not in avail[args.instance_type]["available_regions"]:
        region = avail[args.instance_type]["available_regions"][0]

    print(f"\n{'='*60}")
    print("Full pipeline launch")
    print(f"  RFD3_ROUNDS     : {args.rfd3_rounds}")
    print(f"  N designs (r1)  : {args.n_designs}")
    print(f"  Results folder  : {args.results_dir}")
    print(f"  Region          : {region}")
    print(f"{'='*60}\n")

    inst = client.launch(
        instance_type=args.instance_type,
        region=region,
        ssh_key_names=[args.ssh_key_name],
        name="binary-antibodies-pipeline",
    )
    instance_id = inst["id"]
    print(f"Instance ID: {instance_id}")
    ip = ""

    try:
        active = client.wait_until_active(instance_id)
        ip = active["ip"]
        wait_ssh(ip, args.ssh_key)
        upload_pipeline(ip, args.ssh_key)
        run_setup(ip, args.ssh_key)
        run_pipeline(ip, args.ssh_key, args, gh_tok)
        print(f"\nResults pushed to {args.results_dir} on branch {args.github_branch}")
    finally:
        if not args.no_terminate:
            print(f"\nTerminating {instance_id} ...")
            client.terminate(instance_id)
        else:
            print(f"\nInstance still running: {instance_id} @ {ip}")
            print(f"  ssh -i {args.ssh_key} {REMOTE_USER}@{ip}")


if __name__ == "__main__":
    main()
