#!/usr/bin/env python3
"""
launch_stage_a.py
-----------------
Launch Lambda GPU for Stage A hidden minibinder RFd3 design.

Bootstraps the instance via cloud-init (git clone) so SSH is only needed
for monitoring, not for file upload.

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
import textwrap
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from binary_antibodies.lambda_client import LambdaClient, resolve_lambda_api_key

REMOTE_PIPELINE = "/home/ubuntu/pipeline"
GITHUB_REPO = "Cookiemaster33/binary_antibodies"
GITHUB_BRANCH = "cursor/conditional-nanobody-design-992c"
INPUT_PDB = "fab_hidden_minibinder_stage_a.pdb"
CONFIG_JSON = "stage_a_hidden_minibinder_config.json"


def ssh_opts(key: str) -> list[str]:
    return ["-i", os.path.expanduser(key), "-o", "StrictHostKeyChecking=no", "-o", "ConnectTimeout=15"]


def bootstrap_user_data(n_designs: int) -> str:
    """Cloud-init script: clone repo, build target, run setup + Stage A pipeline."""
    return textwrap.dedent(f"""\
        #cloud-config
        runcmd:
          - |
            set -eo pipefail
            exec > /home/ubuntu/stage_a_bootstrap.log 2>&1
            export HOME=/home/ubuntu
            cd /home/ubuntu
            rm -rf binary_antibodies
            git clone --depth 1 --branch {GITHUB_BRANCH} \\
              https://github.com/{GITHUB_REPO}.git binary_antibodies
            cd binary_antibodies
            python3 scripts/build_stage_a_design_target.py
            mkdir -p {REMOTE_PIPELINE}/inputs
            cp structures/domains/{INPUT_PDB} {REMOTE_PIPELINE}/inputs/
            cp structures/interface/{CONFIG_JSON} {REMOTE_PIPELINE}/inputs/
            cp scripts/gpu_setup/run_stage_a_hidden_minibinder.sh {REMOTE_PIPELINE}/run_stage_a.sh
            cp scripts/gpu_setup/setup_pipeline_rfd3.sh {REMOTE_PIPELINE}/
            chmod +x {REMOTE_PIPELINE}/*.sh
            chown -R ubuntu:ubuntu /home/ubuntu/binary_antibodies {REMOTE_PIPELINE}
            sudo -u ubuntu bash {REMOTE_PIPELINE}/setup_pipeline_rfd3.sh
            sudo -u ubuntu tmux new-session -d -s stage_a \\
              "N_DESIGNS={n_designs} INPUT_PDB={INPUT_PDB} CONFIG_JSON={CONFIG_JSON} \\
               bash {REMOTE_PIPELINE}/run_stage_a.sh"
            echo "Stage A bootstrap complete: $(date)" >> /home/ubuntu/stage_a_bootstrap.log
    """)


def ssh_key_usable(key: str) -> bool:
    path = os.path.expanduser(key)
    return os.path.isfile(path) and os.access(path, os.R_OK)


def wait_for_pipeline(ip: str, key: str, timeout_s: int = 14400) -> bool:
    """Poll remote log for RFd3 completion."""
    if not ssh_key_usable(key):
        return False
    print(f"  Polling pipeline log on {ip} (up to {timeout_s // 3600}h)...")
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        r = subprocess.run(
            ["ssh"] + ssh_opts(key)
            + [f"ubuntu@{ip}", f"grep -q 'Stage A RFd3 done' {REMOTE_PIPELINE}/stage_a_pipeline.log 2>/dev/null"],
            capture_output=True,
        )
        if r.returncode == 0:
            print("  Stage A RFd3 complete!")
            return True
        time.sleep(120)
    return False


def main() -> None:
    p = argparse.ArgumentParser(description="Launch Stage A hidden minibinder RFd3 on Lambda.")
    p.add_argument("--ssh-key", default=os.environ.get("LAMBDA_SSH_KEY", "~/.ssh/lambda_agent_key"))
    p.add_argument("--ssh-key-name", default="cursor-agent")
    p.add_argument("--n-designs", type=int, default=200)
    p.add_argument("--no-terminate", action="store_true")
    p.add_argument("--branch", default=GITHUB_BRANCH)
    args = p.parse_args()

    api_key = resolve_lambda_api_key()
    if not api_key:
        sys.exit("ERROR: Lambda API key not found")

    subprocess.run([sys.executable, str(ROOT / "scripts/build_stage_a_design_target.py")], check=True)

    client = LambdaClient(api_key=api_key)
    user_data = bootstrap_user_data(args.n_designs)
    inst = client.launch(
        "gpu_1x_a100_sxm4",
        "us-east-1",
        ssh_key_names=[args.ssh_key_name],
        name="stage-a-hidden-minibinder",
        user_data=user_data,
    )
    iid = inst["id"]
    ip = ""
    try:
        active = client.wait_until_active(iid)
        ip = active["ip"]
        print(f"\nInstance {iid} @ {ip}")
        print(f"  Branch: {args.branch}")
        print(f"  Designs: {args.n_designs}")
        print("\nBootstrap runs automatically via cloud-init (git clone → setup → RFd3).")

        if ssh_key_usable(args.ssh_key):
            print(f"\nMonitor:")
            print(f"  ssh -i {args.ssh_key} ubuntu@{ip} 'tail -f {REMOTE_PIPELINE}/stage_a_pipeline.log'")
            print(f"  ssh -i {args.ssh_key} ubuntu@{ip} 'tail -f /home/ubuntu/stage_a_bootstrap.log'")
            if wait_for_pipeline(ip, args.ssh_key):
                subprocess.run(
                    ["ssh"] + ssh_opts(args.ssh_key)
                    + [f"ubuntu@{ip}", f"ls {REMOTE_PIPELINE}/outputs/rfd3_stage_a/*.cif 2>/dev/null | wc -l"],
                )
                if not args.no_terminate:
                    print(f"\nTerminating {iid} ...")
                    client.terminate(iid)
            elif not args.no_terminate:
                print("  Pipeline still running — leaving instance up (or use --no-terminate).")
        else:
            print(f"\nSSH key not found at {args.ssh_key} — cannot monitor from here.")
            print(f"Instance is bootstrapping. Check Lambda console or add SSH key to monitor:")
            print(f"  ssh -i <your-key> ubuntu@{ip} 'tail -f {REMOTE_PIPELINE}/stage_a_pipeline.log'")
            print(f"\nInstance left running: {iid}")

    except Exception:
        print(f"\nFailed — instance left running for debug: {iid} @ {ip}")
        raise


if __name__ == "__main__":
    main()
