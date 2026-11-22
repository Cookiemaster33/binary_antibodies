#!/usr/bin/env python3
"""
One-off: Boltz-2 templated prediction of peptide LGQQQPFPPQQPY vs
CFTR structures 5UAK and 6O2P on a Lambda GPU instance.

Usage:
    python3 scripts/launch_peptide_cftr_boltz2.py

Outputs fetched to: pipeline_results/peptide_cftr_boltz2/
"""
import json, os, subprocess, sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.launch_stage_a import (
    LambdaClient, resolve_ssh_key,
    ssh, scp, scp_from,
    wait_until_reachable, REMOTE_USER,
    INSTANCE_PREFERENCE, ARM_INSTANCE_TYPES,
)

ROOT            = Path(__file__).resolve().parents[1]
REMOTE_PIPELINE = "/home/ubuntu/pipeline"
SESSION         = "peptide_cftr"


def launch_with_retry(client: LambdaClient, key_name: str) -> str:
    from requests.exceptions import HTTPError
    avail_map  = {t["name"]: t for t in client.available_instance_types()
                  if t["name"] not in ARM_INSTANCE_TYPES}
    ordered    = [n for n in INSTANCE_PREFERENCE if n in avail_map]
    ordered   += sorted(n for n in avail_map if n not in INSTANCE_PREFERENCE)
    candidates = [(itype, region, avail_map[itype]["price_per_hour"])
                  for itype in ordered
                  for region in avail_map[itype]["available_regions"]]
    for itype, region, price in candidates:
        print(f"  Trying {itype} in {region} (${price:.2f}/h)...")
        try:
            iid = client.launch(itype, region,
                                ssh_key_names=[key_name],
                                name="peptide-cftr-boltz2")["id"]
            print(f"  Launched: {iid}")
            return iid
        except (HTTPError, Exception) as e:
            print(f"  Failed ({e}), next...")
            time.sleep(2)
    raise RuntimeError("No instance launched — all regions returned errors")


def upload_files(ip: str, key: str) -> int:
    yaml_dir     = ROOT / "structures" / "peptide_cftr_boltz2"
    template_dir = ROOT / "structures" / "peptide_binding_screen"
    run_script   = ROOT / "scripts" / "gpu_setup" / "run_peptide_cftr_boltz2.sh"

    yamls     = sorted(yaml_dir.glob("*.yaml"))
    templates = sorted(template_dir.glob("*.cif"))

    if not yamls:
        raise FileNotFoundError(f"No YAML files in {yaml_dir}")
    if not templates:
        raise FileNotFoundError(f"No CIF files in {template_dir}")

    ssh(ip, key,
        f"mkdir -p {REMOTE_PIPELINE}/yaml "
        f"{REMOTE_PIPELINE}/templates "
        f"{REMOTE_PIPELINE}/results")

    print(f"  Uploading {len(yamls)} YAML + {len(templates)} CIF + run script...")
    for f in yamls:
        scp(str(f), ip, f"{REMOTE_PIPELINE}/yaml/{f.name}", key)
    for f in templates:
        scp(str(f), ip, f"{REMOTE_PIPELINE}/templates/{f.name}", key)
    scp(str(run_script), ip, f"{REMOTE_PIPELINE}/run_peptide_cftr_boltz2.sh", key)
    ssh(ip, key, f"chmod +x {REMOTE_PIPELINE}/run_peptide_cftr_boltz2.sh")

    print("  Upload complete.")
    return len(yamls)


def run_pipeline(ip: str, key: str) -> None:
    print(f"  Starting predictions in tmux session '{SESSION}'...")
    ssh(ip, key,
        f"tmux kill-session -t {SESSION} 2>/dev/null || true",
        check=False, quiet=True)
    ssh(ip, key,
        f"tmux new-session -d -s {SESSION} -c {REMOTE_PIPELINE} "
        f"'bash {REMOTE_PIPELINE}/run_peptide_cftr_boltz2.sh "
        f"2>&1 | tee {REMOTE_PIPELINE}/run.log'")


def wait_pipeline(ip: str, key: str, n_yaml: int, timeout_s: int = 7200) -> bool:
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        done = int(subprocess.run(
            ["ssh", "-i", key,
             "-o", "StrictHostKeyChecking=no", "-o", "BatchMode=yes",
             "-o", "ConnectTimeout=15",
             f"{REMOTE_USER}@{ip}",
             f"ls {REMOTE_PIPELINE}/boltz_results/*/predictions/confidence_*.json "
             f"2>/dev/null | wc -l"],
            capture_output=True, text=True,
        ).stdout.strip() or "0")

        print(f"  {time.strftime('%H:%M:%S')}  {done}/{n_yaml} predictions done",
              flush=True)
        if done >= n_yaml:
            return True

        session_up = ssh(ip, key, f"tmux has-session -t {SESSION}",
                         check=False, quiet=True) == 0
        if not session_up:
            done_check = ssh(ip, key,
                             f"test -f {REMOTE_PIPELINE}/results/iptm_scores.json",
                             check=False, quiet=True) == 0
            if done_check:
                print("  Session ended; results JSON present — success.")
                return True
            ssh(ip, key, f"tail -40 {REMOTE_PIPELINE}/run.log", check=False)
            raise RuntimeError("Pipeline session ended unexpectedly")
        time.sleep(120)
    return False


def fetch_results(ip: str, key: str, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"  Fetching results → {out_dir}")
    ssh(ip, key,
        f"tar czf {REMOTE_PIPELINE}/peptide_cftr_results.tar.gz "
        f"-C {REMOTE_PIPELINE} results/ run.log boltz_run.log 2>/dev/null || true")
    scp_from(ip, f"{REMOTE_PIPELINE}/peptide_cftr_results.tar.gz",
             str(out_dir), key)
    subprocess.run(["tar", "xzf", str(out_dir / "peptide_cftr_results.tar.gz"),
                    "-C", str(out_dir)], check=False)

    # Also fetch structure predictions (CIF files)
    ssh(ip, key,
        f"tar czf {REMOTE_PIPELINE}/boltz_structures.tar.gz "
        f"-C {REMOTE_PIPELINE} boltz_results/ 2>/dev/null || true")
    scp_from(ip, f"{REMOTE_PIPELINE}/boltz_structures.tar.gz",
             str(out_dir), key)
    subprocess.run(["tar", "xzf", str(out_dir / "boltz_structures.tar.gz"),
                    "-C", str(out_dir)], check=False)
    print("  Fetch complete.")


def summarise(out_dir: Path) -> None:
    score_file = out_dir / "results" / "iptm_scores.json"
    if not score_file.exists():
        # Try direct path
        for p in out_dir.rglob("iptm_scores.json"):
            score_file = p
            break
    if not score_file.exists():
        print("  No iptm_scores.json found.")
        return

    scores = json.loads(score_file.read_text())
    print(f"\n{'Name':<30} {'ipTM(A↔B)':>10} {'pTM':>6} {'pLDDT':>7}")
    print("-" * 58)
    for s in scores:
        flag = " ★" if s["iptm_AB"] > 0.5 else ""
        print(f"{s['name']:<30} {s['iptm_AB']:>10.3f} {s['ptm']:>6.3f} {s['plddt']:>7.1f}{flag}")


def main() -> None:
    client   = LambdaClient(os.environ["LAMBDA_API_KEY"])
    ssh_key, key_name = resolve_ssh_key(client, os.environ.get("LAMBDA_SSH_KEY", ""),
                                        "cursor-agent")

    n_yaml = len(list((ROOT / "structures" / "peptide_cftr_boltz2").glob("*.yaml")))
    print(f"  Peptide-CFTR Boltz-2 screen: {n_yaml} predictions (LGQQQPFPPQQPY vs 5UAK, 6O2P)")

    iid = launch_with_retry(client, key_name)
    ip  = wait_until_reachable(client, iid, ssh_key, 1800)

    try:
        upload_files(ip, ssh_key)
        run_pipeline(ip, ssh_key)
        wait_pipeline(ip, ssh_key, n_yaml)
        out_dir = ROOT / "pipeline_results" / "peptide_cftr_boltz2"
        fetch_results(ip, ssh_key, out_dir)
        summarise(out_dir)
    finally:
        try:
            client.terminate(iid)
            print(f"  Terminated: {iid}")
        except Exception as e:
            print(f"  WARNING: terminate failed ({e}) — {iid} may still be running")


if __name__ == "__main__":
    main()
