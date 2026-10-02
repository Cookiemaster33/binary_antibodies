#!/usr/bin/env python3
"""Launch Stage A validation (ProteinMPNN + Boltz-2 SC) on a Lambda GPU instance.

Pipeline:
  1. Upload top-25 Stage A grafted PDBs + run script
  2. Run ProteinMPNN (8 seqs per design, chain M, fixed A/B/C/D)
  3. Run Boltz-2 self-consistency on all 200 minibinder+VL+CH1 complexes
  4. Compute scRMSD + chain-pair ipTM
  5. Fetch results → pipeline_results/stage_a_validation/
"""
import argparse, json, os, subprocess, sys, time
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
REMOTE_TMUX     = "tmux"
SESSION         = "stage_a_val"


# ── launch helpers ─────────────────────────────────────────────────────────────

def launch_with_retry(client: LambdaClient, key_name: str) -> str:
    from requests.exceptions import HTTPError

    avail_raw  = client.available_instance_types()
    avail_map  = {t["name"]: t for t in avail_raw if t["name"] not in ARM_INSTANCE_TYPES}
    ordered    = [n for n in INSTANCE_PREFERENCE if n in avail_map]
    ordered   += sorted(n for n in avail_map if n not in INSTANCE_PREFERENCE)

    candidates = [
        (itype, region, avail_map[itype]["price_per_hour"])
        for itype in ordered
        for region in avail_map[itype]["available_regions"]
    ]

    for itype, region, price in candidates:
        print(f"  Trying {itype} in {region} (${price:.2f}/h)...")
        try:
            iid = client.launch(
                itype, region,
                ssh_key_names=[key_name],
                name="stage-a-validation",
            )["id"]
            print(f"  Launched: {iid}")
            return iid
        except (HTTPError, Exception) as exc:
            print(f"  Launch failed ({exc}), trying next...")
            time.sleep(2)

    raise RuntimeError("Could not launch any instance — all regions returned errors")


# ── upload / run helpers ────────────────────────────────────────────────────────

def upload_files(ip: str, key: str, n_mpnn_seqs: int) -> int:
    top_dir = ROOT / "pipeline_results" / "stage_a_selected" / "top_structures"
    pdbs = sorted(top_dir.glob("*.pdb"))
    if not pdbs:
        raise FileNotFoundError(f"No PDBs found in {top_dir}")

    run_script = ROOT / "scripts" / "gpu_setup" / "run_stage_a_validation.sh"

    # Remote directories
    ssh(ip, key, f"mkdir -p {REMOTE_PIPELINE}/top_structures {REMOTE_PIPELINE}/stage_a_validation")

    # Upload PDBs
    print(f"  Uploading {len(pdbs)} PDBs...")
    for pdb in pdbs:
        scp(str(pdb), ip, f"{REMOTE_PIPELINE}/top_structures/{pdb.name}", key)

    # Upload run script
    scp(str(run_script), ip, f"{REMOTE_PIPELINE}/run_stage_a_validation.sh", key)
    ssh(ip, key, f"chmod +x {REMOTE_PIPELINE}/run_stage_a_validation.sh")

    print(f"  Upload complete ({len(pdbs)} PDBs + run script)")
    return len(pdbs)


def run_validation(ip: str, key: str, n_mpnn_seqs: int) -> None:
    print(f"  Starting validation in tmux session '{SESSION}'...")
    ssh(ip, key,
        f"{REMOTE_TMUX} kill-session -t {SESSION} 2>/dev/null || true",
        check=False, quiet=True)
    env = f"N_MPNN_SEQS={n_mpnn_seqs}"
    ssh(ip, key,
        f"{REMOTE_TMUX} new-session -d -s {SESSION} -c {REMOTE_PIPELINE} "
        f"'env {env} bash {REMOTE_PIPELINE}/run_stage_a_validation.sh "
        f"2>&1 | tee {REMOTE_PIPELINE}/stage_a_val.log'")


def wait_validation(ip: str, key: str, n_pdbs: int, n_seqs: int,
                    timeout_s: int = 14400) -> bool:
    """Poll until Boltz-2 confidence JSON files reach n_pdbs * n_seqs or session ends."""
    total_expected = n_pdbs * n_seqs
    deadline = time.time() + timeout_s

    while time.time() < deadline:
        # Count completed Boltz-2 predictions
        done = int(subprocess.run(
            ["ssh", "-i", key,
             "-o", "StrictHostKeyChecking=no",
             "-o", "BatchMode=yes",
             "-o", "ConnectTimeout=15",
             f"{REMOTE_USER}@{ip}",
             f"ls {REMOTE_PIPELINE}/boltz_sc_outputs/*/predictions/confidence_*.json "
             f"2>/dev/null | wc -l"],
            capture_output=True, text=True,
        ).stdout.strip() or "0")

        # Also report MPNN status
        mpnn_done = int(subprocess.run(
            ["ssh", "-i", key,
             "-o", "StrictHostKeyChecking=no",
             "-o", "BatchMode=yes",
             "-o", "ConnectTimeout=15",
             f"{REMOTE_USER}@{ip}",
             f"ls {REMOTE_PIPELINE}/mpnn_outputs/seqs/*.fa 2>/dev/null | wc -l"],
            capture_output=True, text=True,
        ).stdout.strip() or "0")

        print(f"  {time.strftime('%H:%M:%S')}  MPNN: {mpnn_done}/{n_pdbs} designs | "
              f"Boltz-2 SC: {done}/{total_expected}", flush=True)

        if done >= total_expected:
            return True

        session_alive = ssh(ip, key,
                            f"{REMOTE_TMUX} has-session -t {SESSION}",
                            check=False, quiet=True) == 0
        if not session_alive:
            # Check if results CSV was written (script may have finished cleanly)
            csv_exists = ssh(ip, key,
                             f"test -f {REMOTE_PIPELINE}/stage_a_validation/stage_a_boltz2_sc.csv",
                             check=False, quiet=True) == 0
            if csv_exists:
                print("  Session ended; results CSV present — assuming success.")
                return True
            # Print last log lines for diagnosis
            ssh(ip, key, f"tail -30 {REMOTE_PIPELINE}/stage_a_val.log", check=False)
            raise RuntimeError("Validation session ended before all predictions completed")

        time.sleep(120)

    return False


def fetch_results(ip: str, key: str, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"  Fetching results → {out_dir}")

    # Pack and download
    ssh(ip, key,
        f"cd {REMOTE_PIPELINE}/stage_a_validation && "
        f"tar czf {REMOTE_PIPELINE}/stage_a_validation.tar.gz .")
    scp_from(ip, f"{REMOTE_PIPELINE}/stage_a_validation.tar.gz", str(out_dir), key)
    subprocess.run(["tar", "xzf", str(out_dir / "stage_a_validation.tar.gz"),
                    "-C", str(out_dir)], check=False)

    # Also grab full log
    scp_from(ip, f"{REMOTE_PIPELINE}/stage_a_val.log",
             str(out_dir / "stage_a_val.log"), key)


def summarise(out_dir: Path) -> None:
    csv_path = out_dir / "stage_a_boltz2_sc.csv"
    if not csv_path.exists():
        print("  No results CSV found — check stage_a_val.log for errors.")
        return

    import csv
    rows = list(csv.DictReader(open(csv_path)))
    valid = [r for r in rows if r["sc_rmsd_A"] and r["sc_rmsd_A"] != "nan"]
    valid.sort(key=lambda r: float(r["sc_rmsd_A"]))

    print(f"\n  Total predictions: {len(rows)}  Valid (structure found): {len(valid)}")
    print(f"\n{'Rank':<5} {'scRMSD':>7} {'ipTM_VL':>8} {'ipTM_CH1':>9} {'pLDDT':>7}  Design")
    print("  " + "-" * 70)
    for i, r in enumerate(valid[:20], 1):
        sc = float(r["sc_rmsd_A"])
        flag = " ★" if sc < 2.0 and float(r.get("mean_iptm", 0)) > 0.5 else ""
        print(f"{i:<5} {sc:>7.3f} {float(r['iptm_mb_vl']):>8.3f} "
              f"{float(r['iptm_mb_ch1']):>9.3f} {float(r['plddt']):>7.1f}  "
              f"{r['yaml_stem']}{flag}")

    passing = [r for r in valid if float(r["sc_rmsd_A"]) < 2.5]
    print(f"\n  Passing scRMSD < 2.5 Å: {len(passing)}/{len(valid)}")

    fasta = out_dir / "stage_a_validated.fasta"
    if fasta.exists():
        n = fasta.read_text().count(">")
        print(f"  Validated FASTA: {fasta}  ({n} sequences)")


# ── main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(description="Launch Stage A ProteinMPNN + Boltz-2 SC validation")
    p.add_argument("--results-dir", default="pipeline_results/stage_a_validation")
    p.add_argument("--n-mpnn-seqs", type=int, default=8,
                   help="ProteinMPNN sequences per design (default: 8)")
    p.add_argument("--no-terminate", action="store_true")
    p.add_argument("--ssh-key", default=os.environ.get("LAMBDA_SSH_KEY", ""))
    args = p.parse_args()

    client   = LambdaClient(os.environ["LAMBDA_API_KEY"])
    ssh_key, key_name = resolve_ssh_key(client, args.ssh_key, "cursor-agent")

    n_pdbs = len(list(
        (ROOT / "pipeline_results" / "stage_a_selected" / "top_structures").glob("*.pdb")
    ))
    total_predictions = n_pdbs * args.n_mpnn_seqs
    print(f"  Stage A validation: {n_pdbs} PDBs × {args.n_mpnn_seqs} MPNN seqs"
          f" = {total_predictions} Boltz-2 predictions")

    iid = launch_with_retry(client, key_name)

    ip = wait_until_reachable(client, iid, ssh_key, 1800)

    try:
        upload_files(ip, ssh_key, args.n_mpnn_seqs)
        run_validation(ip, ssh_key, args.n_mpnn_seqs)
        wait_validation(ip, ssh_key, n_pdbs, args.n_mpnn_seqs)

        out_dir = ROOT / args.results_dir
        fetch_results(ip, ssh_key, out_dir)
        summarise(out_dir)

    finally:
        if not args.no_terminate:
            try:
                client.terminate(iid)
                print(f"  Terminated: {iid}")
            except Exception as exc:
                print(f"  WARNING: terminate call failed ({exc}) — instance {iid} may still be running")


if __name__ == "__main__":
    main()
