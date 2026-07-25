#!/usr/bin/env python3
"""
launch_stage_a.py
-----------------
Launch Lambda GPU for Stage A hidden minibinder RFd3 design.

Usage
-----
    python scripts/build_stage_a_design_target.py   # build PDB + config first
    python scripts/launch_stage_a.py --ssh-key ~/.ssh/lambda_agent_key
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

REMOTE_PIPELINE = "/home/ubuntu/pipeline"
INPUT_PDB = "fab_hidden_minibinder_stage_a.pdb"
CONFIG_JSON = "stage_a_hidden_minibinder_config.json"


def ssh_opts(key: str) -> list[str]:
    return ["-i", os.path.expanduser(key), "-o", "StrictHostKeyChecking=no", "-o", "ConnectTimeout=15"]


def main() -> None:
    p = argparse.ArgumentParser(description="Launch Stage A hidden minibinder RFd3 on Lambda.")
    p.add_argument("--ssh-key", default=os.environ.get("LAMBDA_SSH_KEY", "~/.ssh/lambda_agent_key"))
    p.add_argument("--ssh-key-name", default="cursor-agent")
    p.add_argument("--n-designs", type=int, default=200)
    p.add_argument("--no-terminate", action="store_true")
    args = p.parse_args()

    api_key = resolve_lambda_api_key()
    if not api_key:
        sys.exit("ERROR: Lambda API key not found")

    # Ensure design target exists
    subprocess.run([sys.executable, str(ROOT / "scripts/build_stage_a_design_target.py")], check=True)

    client = LambdaClient(api_key=api_key)
    inst = client.launch("gpu_1x_a100_sxm4", "us-east-1", ssh_key_names=[args.ssh_key_name])
    iid = inst["id"]
    try:
        active = client.wait_until_active(iid)
        ip = active["ip"]
        print(f"Instance {iid} @ {ip}")

        subprocess.run(["ssh"] + ssh_opts(args.ssh_key) + [f"ubuntu@{ip}", f"mkdir -p {REMOTE_PIPELINE}/inputs"], check=True)
        for local, remote in [
            (ROOT / "structures/domains" / INPUT_PDB, f"{REMOTE_PIPELINE}/inputs/{INPUT_PDB}"),
            (ROOT / "structures/interface" / CONFIG_JSON, f"{REMOTE_PIPELINE}/inputs/{CONFIG_JSON}"),
            (ROOT / "scripts/gpu_setup/run_stage_a_hidden_minibinder.sh", f"{REMOTE_PIPELINE}/run_stage_a.sh"),
            (ROOT / "scripts/gpu_setup/setup_pipeline_rfd3.sh", f"{REMOTE_PIPELINE}/setup_pipeline_rfd3.sh"),
        ]:
            subprocess.run(["scp"] + ssh_opts(args.ssh_key) + [str(local), f"ubuntu@{ip}:{remote}"], check=True)

        subprocess.run(
            ["ssh"] + ssh_opts(args.ssh_key)
            + [f"ubuntu@{ip}", f"chmod +x {REMOTE_PIPELINE}/run_stage_a.sh && bash {REMOTE_PIPELINE}/setup_pipeline_rfd3.sh"],
            check=True,
        )
        env = f"N_DESIGNS={args.n_designs} INPUT_PDB={INPUT_PDB} CONFIG_JSON={CONFIG_JSON}"
        subprocess.run(
            ["ssh"] + ssh_opts(args.ssh_key)
            + [f"ubuntu@{ip}", f"tmux new-session -d -s stage_a '{env} bash {REMOTE_PIPELINE}/run_stage_a.sh'"],
            check=True,
        )
        print(f"Monitor: ssh -i {args.ssh_key} ubuntu@{ip} 'tail -f {REMOTE_PIPELINE}/stage_a_pipeline.log'")

        if not args.no_terminate:
            print("Leaving instance running (use --no-terminate). Poll log manually.")
    except Exception:
        print(f"Failed — instance {iid} left running for debug")
        raise


if __name__ == "__main__":
    main()
