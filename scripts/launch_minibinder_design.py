#!/usr/bin/env python3
"""
launch_minibinder_design.py
---------------------------
Launch a Lambda Cloud GPU instance, run the full minibinder design
pipeline (RFdiffusion → ProteinMPNN → AF2-Multimer), download results.

Requirements
------------
  - LAMBDA_API_KEY environment variable set
  - SSH private key path set via --ssh-key or SSH_PRIVATE_KEY_PATH env var
    (must correspond to the 'mactoby' public key registered in Lambda)

Usage
-----
    # Full run (launch instance, run pipeline, download, terminate)
    python scripts/launch_minibinder_design.py \\
        --ssh-key ~/.ssh/id_ed25519 \\
        --n-designs 200 \\
        --output-dir pipeline_results

    # Check status of a running job
    python scripts/launch_minibinder_design.py --status --instance-id <id>

    # Download results from a finished job without terminating
    python scripts/launch_minibinder_design.py \\
        --download-only --instance-id <id> --output-dir pipeline_results
"""

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from binary_antibodies.lambda_client import LambdaClient

# ── Defaults ────────────────────────────────────────────────────────
DEFAULT_INSTANCE_TYPE = "gpu_1x_a100_sxm4"
DEFAULT_REGION = "us-east-1"
SSH_KEY_NAME = "mactoby"
REMOTE_USER = "ubuntu"
INPUT_PDB_LOCAL = "structures/domains/trastuzumab_VL.pdb"
REMOTE_PIPELINE_DIR = "/home/ubuntu/pipeline"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Launch Lambda Cloud GPU instance and run minibinder design.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--ssh-key", type=str,
                   default=os.environ.get("SSH_PRIVATE_KEY_PATH", "~/.ssh/id_ed25519"),
                   help="Path to SSH private key (must match 'mactoby' key in Lambda).")
    p.add_argument("--instance-type", default=DEFAULT_INSTANCE_TYPE)
    p.add_argument("--region", default=DEFAULT_REGION)
    p.add_argument("--n-designs", type=int, default=200,
                   help="Number of RFdiffusion designs to generate.")
    p.add_argument("--n-mpnn-seqs", type=int, default=8,
                   help="ProteinMPNN sequences per backbone.")
    p.add_argument("--output-dir", default="pipeline_results",
                   help="Local directory for downloaded results.")
    p.add_argument("--no-terminate", action="store_true",
                   help="Do NOT terminate the instance after the run.")
    p.add_argument("--status", action="store_true",
                   help="Check status of a running instance.")
    p.add_argument("--download-only", action="store_true",
                   help="Download results from an already-finished instance.")
    p.add_argument("--instance-id", type=str,
                   help="Instance ID (for --status or --download-only).")
    p.add_argument("--api-key", type=str,
                   default=os.environ.get("LAMBDA_API_KEY", ""),
                   help="Lambda Cloud API key (or set LAMBDA_API_KEY env var).")
    p.add_argument("--setup-only", action="store_true",
                   help="Launch instance and run setup, but do not start the design job.")
    return p.parse_args()


# ── SSH helpers ──────────────────────────────────────────────────────

def ssh_opts(ssh_key: str) -> list[str]:
    return [
        "-i", os.path.expanduser(ssh_key),
        "-o", "StrictHostKeyChecking=no",
        "-o", "ConnectTimeout=10",
        "-o", "ServerAliveInterval=30",
    ]


def wait_for_ssh(ip: str, ssh_key: str, timeout_s: int = 300) -> None:
    """Poll until SSH is accepting connections."""
    print(f"  Waiting for SSH on {ip}", end="", flush=True)
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        try:
            result = subprocess.run(
                ["ssh"] + ssh_opts(ssh_key) + [f"{REMOTE_USER}@{ip}", "echo ok"],
                capture_output=True, timeout=15, text=True,
            )
            if result.returncode == 0:
                print(" ready!")
                return
        except Exception:
            pass
        print(".", end="", flush=True)
        time.sleep(10)
    raise TimeoutError(f"SSH not available on {ip} after {timeout_s}s")


def ssh_run(ip: str, ssh_key: str, command: str, check: bool = True) -> int:
    """Run a command on the remote instance via SSH."""
    result = subprocess.run(
        ["ssh"] + ssh_opts(ssh_key) + [f"{REMOTE_USER}@{ip}", command],
        check=False,
    )
    if check and result.returncode != 0:
        raise RuntimeError(f"Remote command failed (exit {result.returncode}): {command}")
    return result.returncode


def scp_to(local: str, ip: str, remote: str, ssh_key: str) -> None:
    """Upload a file to the instance."""
    subprocess.run(
        ["scp"] + ssh_opts(ssh_key) + [local, f"{REMOTE_USER}@{ip}:{remote}"],
        check=True,
    )


def scp_from(ip: str, remote: str, local: str, ssh_key: str) -> None:
    """Download a file/dir from the instance (rsync for directories)."""
    if remote.endswith("/"):
        subprocess.run(
            ["rsync", "-avz", "-e", "ssh " + " ".join(ssh_opts(ssh_key)),
             f"{REMOTE_USER}@{ip}:{remote}", local],
            check=True,
        )
    else:
        subprocess.run(
            ["scp"] + ssh_opts(ssh_key) + [f"{REMOTE_USER}@{ip}:{remote}", local],
            check=True,
        )


# ── Pipeline steps ───────────────────────────────────────────────────

def upload_inputs(ip: str, ssh_key: str) -> None:
    """Upload the VL1 PDB and pipeline scripts."""
    print("  Uploading inputs ...")
    ssh_run(ip, ssh_key, f"mkdir -p {REMOTE_PIPELINE_DIR}/inputs")
    scp_to(INPUT_PDB_LOCAL, ip, f"{REMOTE_PIPELINE_DIR}/inputs/", ssh_key)
    scp_to("scripts/gpu_setup/setup_pipeline.sh", ip,
           f"{REMOTE_PIPELINE_DIR}/setup_pipeline.sh", ssh_key)
    scp_to("scripts/gpu_setup/run_minibinder_design.sh", ip,
           f"{REMOTE_PIPELINE_DIR}/run_minibinder_design.sh", ssh_key)
    ssh_run(ip, ssh_key,
            f"chmod +x {REMOTE_PIPELINE_DIR}/setup_pipeline.sh "
            f"{REMOTE_PIPELINE_DIR}/run_minibinder_design.sh")
    print("  Upload done.")


def run_setup(ip: str, ssh_key: str) -> None:
    """Run the setup script on the instance (15–25 min)."""
    print("  Running pipeline setup (RFdiffusion + ProteinMPNN + AF2 install)...")
    print("  This takes ~20 minutes. Tailing log in background ...")
    # Start setup in tmux so it survives SSH disconnect
    ssh_run(ip, ssh_key,
            f"tmux new-session -d -s setup -c {REMOTE_PIPELINE_DIR} "
            f"'bash {REMOTE_PIPELINE_DIR}/setup_pipeline.sh 2>&1 | tee setup_pipeline.log'")
    # Poll log until "Setup complete"
    print("  Waiting for setup to finish", end="", flush=True)
    while True:
        rc = ssh_run(ip, ssh_key,
                     f"grep -q 'Setup complete' {REMOTE_PIPELINE_DIR}/setup_pipeline.log",
                     check=False)
        if rc == 0:
            print(" done!")
            break
        print(".", end="", flush=True)
        time.sleep(30)


def run_design(ip: str, ssh_key: str, n_designs: int, n_mpnn_seqs: int) -> None:
    """Launch the design pipeline in a tmux session."""
    print(f"  Starting minibinder design ({n_designs} RFdiffusion + {n_mpnn_seqs}×MPNN + AF2) ...")
    env = f"N_DESIGNS={n_designs} N_MPNN_SEQS={n_mpnn_seqs}"
    ssh_run(ip, ssh_key,
            f"tmux new-session -d -s design -c {REMOTE_PIPELINE_DIR} "
            f"'{env} bash {REMOTE_PIPELINE_DIR}/run_minibinder_design.sh'")
    print(f"  Design job running in tmux session 'design'.")
    print(f"  Monitor: ssh -i {ssh_key} {REMOTE_USER}@{ip} 'tmux attach -t design'")

    # Poll until completion
    print("  Waiting for design pipeline to finish", end="", flush=True)
    while True:
        rc = ssh_run(ip, ssh_key,
                     f"grep -q 'Pipeline complete' {REMOTE_PIPELINE_DIR}/run.log",
                     check=False)
        if rc == 0:
            print(" done!")
            break
        # Also check if tmux session is still running
        tmux_rc = ssh_run(ip, ssh_key,
                          "tmux has-session -t design", check=False)
        if tmux_rc != 0:
            # Session ended — check for error
            ssh_run(ip, ssh_key,
                    f"tail -20 {REMOTE_PIPELINE_DIR}/run.log", check=False)
            break
        print(".", end="", flush=True)
        time.sleep(60)


def download_results(ip: str, ssh_key: str, output_dir: str) -> None:
    """Download pipeline outputs to local directory."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    print(f"  Downloading results to {output_dir} ...")
    scp_from(ip, f"{REMOTE_PIPELINE_DIR}/outputs/", f"{output_dir}/", ssh_key)
    scp_from(ip, f"{REMOTE_PIPELINE_DIR}/run.log", f"{output_dir}/run.log", ssh_key)
    print(f"  Downloaded.")

    # Print top results if summary exists
    summary = Path(output_dir) / "summary.csv"
    if summary.exists():
        import csv
        with open(summary) as f:
            rows = list(csv.DictReader(f))
        print(f"\n  Top 5 designs by ipTM:")
        for row in rows[:5]:
            print(f"    {row['design']:50s}  ipTM={float(row.get('iptm',0)):.3f}  "
                  f"pTM={float(row.get('ptm',0)):.3f}")


# ── Main ─────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    if not args.api_key:
        print("ERROR: Set LAMBDA_API_KEY environment variable or pass --api-key.")
        sys.exit(1)

    client = LambdaClient(api_key=args.api_key)

    # ── Status check ──────────────────────────────────────────
    if args.status:
        if not args.instance_id:
            instances = client.list_instances()
            if not instances:
                print("No instances running.")
            for inst in instances:
                print(f"  {inst['id']}  {inst.get('name','?')}  "
                      f"status={inst.get('status')}  ip={inst.get('ip')}")
        else:
            inst = client.get_instance(args.instance_id)
            import json
            print(json.dumps(inst, indent=2))
        return

    # ── Download only ─────────────────────────────────────────
    if args.download_only:
        if not args.instance_id:
            print("ERROR: --download-only requires --instance-id.")
            sys.exit(1)
        inst = client.get_instance(args.instance_id)
        ip = inst["ip"]
        print(f"Downloading from {ip} ...")
        download_results(ip, args.ssh_key, args.output_dir)
        return

    # ── Full launch and run ───────────────────────────────────
    print(f"\n{'='*60}")
    print(f"Minibinder Design Pipeline on Lambda Cloud")
    print(f"{'='*60}")
    print(f"  Instance type : {args.instance_type}")
    print(f"  Region        : {args.region}")
    print(f"  N designs     : {args.n_designs}")
    print(f"  N MPNN seqs   : {args.n_mpnn_seqs}")
    print(f"  Input PDB     : {INPUT_PDB_LOCAL}")
    print(f"  Output dir    : {args.output_dir}")
    print(f"  SSH key       : {args.ssh_key}")
    print()

    # Check available capacity first
    avail = {t["name"]: t for t in client.available_instance_types()}
    if args.instance_type not in avail:
        print(f"  WARNING: {args.instance_type} not immediately available.")
        print("  Available GPU types:")
        for name, info in avail.items():
            if info["gpus"] > 0:
                print(f"    {name}  ${info['price_per_hour']:.2f}/hr  "
                      f"regions: {info['available_regions']}")
        sys.exit(1)

    region = args.region
    if region not in avail[args.instance_type]["available_regions"]:
        region = avail[args.instance_type]["available_regions"][0]
        print(f"  Switching region to {region} (requested {args.region} unavailable)")

    print(f"  Estimated cost: ~${avail[args.instance_type]['price_per_hour'] * 1.5:.2f} "
          f"(~90 min for {args.n_designs} designs + AF2)")

    # Launch
    print("\n[1/6] Launching instance ...")
    inst_info = client.launch(
        instance_type=args.instance_type,
        region=region,
        ssh_key_names=[SSH_KEY_NAME],
        name="rfdiffusion-minibinder",
    )
    instance_id = inst_info["id"]
    print(f"  Instance ID: {instance_id}")

    try:
        # Wait for active
        print("\n[2/6] Waiting for instance to become active ...")
        inst = client.wait_until_active(instance_id)
        ip = inst["ip"]

        # Wait for SSH
        print(f"\n[3/6] Waiting for SSH on {ip} ...")
        wait_for_ssh(ip, args.ssh_key)

        # Upload inputs
        print("\n[4/6] Uploading inputs and scripts ...")
        upload_inputs(ip, args.ssh_key)

        # Install pipeline
        if not args.setup_only:
            print("\n[5/6] Installing pipeline (RFdiffusion + ProteinMPNN + AF2) ...")
        run_setup(ip, args.ssh_key)

        if args.setup_only:
            print(f"\n  Setup complete. Instance {instance_id} ({ip}) is ready.")
            print(f"  Connect: ssh -i {args.ssh_key} {REMOTE_USER}@{ip}")
            return

        # Run design
        print("\n[6/6] Running design pipeline ...")
        run_design(ip, args.ssh_key, args.n_designs, args.n_mpnn_seqs)

        # Download results
        print("\n[7/7] Downloading results ...")
        download_results(ip, args.ssh_key, args.output_dir)

    finally:
        if not args.no_terminate:
            print(f"\nTerminating instance {instance_id} ...")
            client.terminate(instance_id)
        else:
            inst = client.get_instance(instance_id)
            print(f"\nInstance {instance_id} still running at {inst.get('ip')}.")
            print(f"Connect: ssh -i {args.ssh_key} {REMOTE_USER}@{inst.get('ip')}")
            print(f"Terminate manually: python scripts/launch_minibinder_design.py "
                  f"--terminate --instance-id {instance_id}")

    print(f"\nDone! Results in: {args.output_dir}/")


if __name__ == "__main__":
    main()
