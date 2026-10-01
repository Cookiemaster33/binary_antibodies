#!/usr/bin/env python3
"""Launch CD3e binding screen on a Lambda GPU instance.

Uploads the YAML inputs, sets up Boltz-2, runs the screen,
fetches confidence JSON files and prints the ipTM ranking.
"""
import argparse, json, os, subprocess, sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.launch_stage_a import (
    LambdaClient, materialise_private_key, resolve_ssh_key,
    ssh, scp, scp_from, select_instance_type,
    wait_until_reachable, REMOTE_USER,
)

ROOT = Path(__file__).resolve().parents[1]
REMOTE_PIPELINE = "/home/ubuntu/pipeline"
REMOTE_TMUX = "tmux"


def upload_files(ip: str, key: str) -> None:
    screen_dir = ROOT / "structures" / "cd3e_screen"
    run_script = ROOT / "scripts" / "gpu_setup" / "run_cd3e_screen.sh"
    setup_script = ROOT / "scripts" / "gpu_setup" / "setup_boltz.sh"

    # Create directories
    ssh(ip, key, f"mkdir -p {REMOTE_PIPELINE}/cd3e_screen {REMOTE_PIPELINE}/cd3e_screen_results")

    # Upload YAML inputs
    for yaml_file in sorted(screen_dir.glob("*.yaml")):
        scp(str(yaml_file), ip, f"{REMOTE_PIPELINE}/cd3e_screen/{yaml_file.name}", key)

    # Upload run scripts
    scp(str(run_script), ip, f"{REMOTE_PIPELINE}/run_cd3e_screen.sh", key)
    scp(str(setup_script), ip, f"{REMOTE_PIPELINE}/setup_boltz.sh", key)
    ssh(ip, key, f"chmod +x {REMOTE_PIPELINE}/run_cd3e_screen.sh {REMOTE_PIPELINE}/setup_boltz.sh")
    print(f"  Uploaded {len(list(screen_dir.glob('*.yaml')))} YAML files + run script")


def run_setup(ip: str, key: str) -> None:
    print("  Running Boltz setup (pip install + weight download)...")
    ssh(ip, key,
        f"{REMOTE_TMUX} new-session -d -s setup -c {REMOTE_PIPELINE} "
        f"'bash {REMOTE_PIPELINE}/setup_boltz.sh 2>&1 | tee {REMOTE_PIPELINE}/setup_boltz.log'")
    while True:
        if ssh(ip, key, f"grep -q 'Setup complete' {REMOTE_PIPELINE}/setup_boltz.log", check=False) == 0:
            print("  Setup complete.")
            return
        if ssh(ip, key, f"{REMOTE_TMUX} has-session -t setup", check=False) != 0:
            ssh(ip, key, f"tail -20 {REMOTE_PIPELINE}/setup_boltz.log", check=False)
            raise RuntimeError("Setup session ended unexpectedly")
        time.sleep(20)


def run_screen(ip: str, key: str) -> None:
    print("  Launching CD3e screen in tmux session 'cd3e_screen'...")
    ssh(ip, key,
        f"{REMOTE_TMUX} kill-session -t cd3e_screen 2>/dev/null || true", check=False, quiet=True)
    ssh(ip, key,
        f"{REMOTE_TMUX} new-session -d -s cd3e_screen -c {REMOTE_PIPELINE} "
        f"'bash {REMOTE_PIPELINE}/run_cd3e_screen.sh 2>&1 | tee {REMOTE_PIPELINE}/cd3e_screen.log'")


def wait_screen(ip: str, key: str, n_yaml: int, timeout_s: int = 7200) -> bool:
    print(f"  Polling screen progress ({n_yaml} complexes)...")
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        done = int(subprocess.run(
            ["ssh", "-i", key, "-o", "StrictHostKeyChecking=no",
             "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
             f"{REMOTE_USER}@{ip}",
             f"ls {REMOTE_PIPELINE}/cd3e_screen_results/*/predictions/confidence_*.json "
             f"2>/dev/null | wc -l"],
            capture_output=True, text=True).stdout.strip() or "0")
        print(f"  {time.strftime('%H:%M:%S')}  {done}/{n_yaml} predictions done", flush=True)
        if done >= n_yaml:
            return True
        # Check session still running
        if ssh(ip, key, f"{REMOTE_TMUX} has-session -t cd3e_screen", check=False, quiet=True) != 0:
            if done >= n_yaml:
                return True
            raise RuntimeError("Screen session ended before all predictions completed")
        time.sleep(120)
    return False


def fetch_results(ip: str, key: str, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"  Fetching results → {out_dir}")
    ssh(ip, key,
        f"cd {REMOTE_PIPELINE}/cd3e_screen_results && "
        f"tar czf {REMOTE_PIPELINE}/cd3e_screen_results.tar.gz .")
    scp_from(ip, f"{REMOTE_PIPELINE}/cd3e_screen_results.tar.gz", str(out_dir), key)
    subprocess.run(["tar", "xzf", str(out_dir / "cd3e_screen_results.tar.gz"),
                    "-C", str(out_dir)], check=False)
    scp_from(ip, f"{REMOTE_PIPELINE}/cd3e_screen.log",
             str(out_dir / "cd3e_screen.log"), key)


def summarise(out_dir: Path) -> None:
    scores = []
    for f in sorted(out_dir.rglob("confidence_*.json")):
        try:
            data = json.loads(f.read_text())
            name = f.parent.parent.name

            # Boltz-2: iptm is a dict {"AB": 0.7, ...}; extract chain-pair AB
            raw_iptm = data.get("iptm", 0.0)
            if isinstance(raw_iptm, dict):
                iptm = float(raw_iptm.get("AB", next(iter(raw_iptm.values()), 0.0)))
            else:
                iptm = float(raw_iptm)

            ptm   = float(data.get("ptm", 0.0))
            plddt = float(data.get("mean_plddt", data.get("complex_plddt", 0.0)))
            scores.append({"design": name, "iptm": iptm, "ptm": ptm, "plddt": plddt})
        except Exception:
            pass
    if not scores:
        print("  No confidence files found.")
        return
    scores.sort(key=lambda x: x["iptm"], reverse=True)
    print(f"\n{'Rank':<5} {'ipTM':>6} {'pTM':>6} {'pLDDT':>7}  Design")
    print("-" * 68)
    for i, s in enumerate(scores, 1):
        flag = " ★" if s["iptm"] > 0.5 else ""
        print(f"{i:<5} {s['iptm']:>6.3f} {s['ptm']:>6.3f} {s['plddt']:>7.1f}  {s['design']}{flag}")

    import csv
    csv_path = out_dir / "cd3e_screen_iptm.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["rank", "design", "iptm", "ptm", "plddt"])
        w.writeheader()
        for i, s in enumerate(scores, 1):
            w.writerow({"rank": i, **s})
    print(f"\n  Saved: {csv_path}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--results-dir", default="pipeline_results/cd3e_screen")
    p.add_argument("--no-terminate", action="store_true")
    p.add_argument("--ssh-key", default=os.environ.get("LAMBDA_SSH_KEY", ""))
    args = p.parse_args()

    client = LambdaClient(os.environ["LAMBDA_API_KEY"])

    n_yaml = len(list((ROOT / "structures" / "cd3e_screen").glob("*.yaml")))
    print(f"  CD3e screen: {n_yaml} complexes")

    itype, region, price = select_instance_type(client, None, None)
    print(f"  Launching {itype} in {region} (${price:.2f}/h)...")
    ssh_key, key_name = resolve_ssh_key(client, args.ssh_key, "cursor-agent")
    iid = client.launch(itype, region, ssh_key_names=[key_name], name="cd3e-screen")["id"]
    print(f"  Instance: {iid}")

    ip = wait_until_reachable(client, iid, ssh_key, 1800)

    try:
        upload_files(ip, ssh_key)
        run_setup(ip, ssh_key)
        run_screen(ip, ssh_key)
        wait_screen(ip, ssh_key, n_yaml)
        out_dir = ROOT / args.results_dir
        fetch_results(ip, ssh_key, out_dir)
        summarise(out_dir)
    finally:
        if not args.no_terminate:
            client.terminate(iid)
            print(f"  Terminated: {iid}")


if __name__ == "__main__":
    main()
