# Lambda Cloud GPU Skill

Use this skill whenever you need to run GPU workloads on Lambda Cloud (launching instances, running scripts, fetching results, terminating).

---

## 1. Prerequisites — Secrets

Before anything works, **two secrets must be added in the Cursor Dashboard** (Settings → Cloud Agents → Secrets). They are injected as environment variables into every agent VM.

| Secret name | What it is |
|---|---|
| `LAMBDA_API_KEY` | Lambda Cloud API key — from https://cloud.lambdalabs.com/api-keys |
| `LAMBDA_SSH_KEY` | Full contents of an SSH private key (OpenSSH format, starts with `-----BEGIN OPENSSH PRIVATE KEY-----`). The corresponding public key must already be registered in your Lambda Cloud account under https://cloud.lambdalabs.com/ssh-keys |

If an agent says "cannot connect to Lambda Cloud", the most common cause is that these secrets are not yet added for that project. Add them in the Dashboard and re-run.

> **Tip — reuse an existing key**: if you already have `cursor-agent` registered in Lambda Cloud, export the matching private key and store it verbatim as the `LAMBDA_SSH_KEY` secret.

---

## 2. Self-Contained Lambda Client

Copy this file into your project as `lambda_client.py` (no external imports beyond `requests`):

```python
"""lambda_client.py — thin Lambda Cloud REST API wrapper."""
from __future__ import annotations
import os, time, base64, subprocess, tempfile
from pathlib import Path
import requests

BASE_URL = "https://cloud.lambda.ai/api/v1"
REMOTE_USER = "ubuntu"


class LambdaClient:
    def __init__(self, api_key: str | None = None):
        self.api_key = api_key or os.environ.get("LAMBDA_API_KEY", "")
        if not self.api_key:
            raise ValueError("Set LAMBDA_API_KEY env var or pass api_key=")
        self._s = requests.Session()
        self._s.headers["Authorization"] = f"Bearer {self.api_key}"

    def _get(self, path):
        r = self._s.get(f"{BASE_URL}{path}", timeout=30); r.raise_for_status()
        return r.json()

    def _post(self, path, payload):
        r = self._s.post(f"{BASE_URL}{path}", json=payload, timeout=30); r.raise_for_status()
        return r.json()

    def available_instance_types(self) -> list[dict]:
        """Return instance types that currently have capacity."""
        all_types = self._get("/instance-types")["data"]
        out = []
        for name, info in all_types.items():
            regions = [r["name"] for r in info.get("regions_with_capacity_available", [])]
            if not regions:
                continue
            it = info.get("instance_type", {}); specs = it.get("specs", {})
            out.append({
                "name": name,
                "price_per_hour": it.get("price_cents_per_hour", 0) / 100,
                "available_regions": regions,
                "gpus": specs.get("gpus", 0),
                "memory_gib": specs.get("memory_gibibytes", 0),
            })
        return out

    def list_ssh_keys(self) -> list[dict]:
        return self._get("/ssh-keys")["data"]

    def add_ssh_key(self, name: str, public_key: str) -> dict:
        return self._post("/ssh-keys", {"name": name, "public_key": public_key})["data"]

    def launch(self, instance_type: str, region: str,
               ssh_key_names: list[str], name: str = "agent") -> dict:
        payload = {"instance_type_name": instance_type, "region_name": region,
                   "ssh_key_names": ssh_key_names, "name": name, "quantity": 1}
        return self._post("/instance-operations/launch", payload)["data"]["instance_ids"][0] \
               if "instance_ids" in self._post("/instance-operations/launch", payload).get("data", {}) \
               else self._post("/instance-operations/launch", payload)["data"]

    def launch(self, instance_type: str, region: str,
               ssh_key_names: list[str], name: str = "agent") -> dict:
        r = self._post("/instance-operations/launch", {
            "instance_type_name": instance_type,
            "region_name": region,
            "ssh_key_names": ssh_key_names,
            "name": name,
            "quantity": 1,
        })
        ids = r["data"].get("instance_ids", [])
        return {"id": ids[0]} if ids else r["data"]

    def get_instance(self, instance_id: str) -> dict:
        return self._get(f"/instances/{instance_id}")["data"]

    def list_instances(self) -> list[dict]:
        return self._get("/instances")["data"]

    def terminate(self, instance_id: str) -> None:
        self._post("/instance-operations/terminate", {"instance_ids": [instance_id]})
```

---

## 3. SSH Helper Functions

```python
import subprocess, os, tempfile, base64
from pathlib import Path

REMOTE_USER = "ubuntu"

def ssh_opts(key_path: str) -> list[str]:
    return ["-i", key_path,
            "-o", "StrictHostKeyChecking=no",
            "-o", "ConnectTimeout=15",
            "-o", "ServerAliveInterval=30",
            "-o", "BatchMode=yes"]

def ssh(ip: str, key: str, cmd: str, check: bool = True, quiet: bool = False) -> int:
    kwargs = {"check": False}
    if quiet:
        kwargs["capture_output"] = True
        kwargs["stdin"] = subprocess.DEVNULL
    r = subprocess.run(["ssh"] + ssh_opts(key) + [f"{REMOTE_USER}@{ip}", cmd], **kwargs)
    if check and r.returncode != 0:
        raise RuntimeError(f"SSH failed ({r.returncode}): {cmd[:100]}")
    return r.returncode

def scp(local: str, ip: str, remote: str, key: str) -> None:
    subprocess.run(["scp"] + ssh_opts(key) + [local, f"{REMOTE_USER}@{ip}:{remote}"], check=True)

def scp_from(ip: str, remote: str, local: str, key: str) -> None:
    subprocess.run(["scp"] + ssh_opts(key) + [f"{REMOTE_USER}@{ip}:{remote}", local], check=True)
```

---

## 4. Resolve SSH Key from Secret

The `LAMBDA_SSH_KEY` secret holds the **raw private key text**. Write it to a temp file before passing to SSH:

```python
def resolve_ssh_key(client: LambdaClient, ssh_key_env: str = "",
                    key_name: str = "cursor-agent") -> tuple[str, str]:
    """
    Returns (path_to_private_key_file, registered_key_name).
    If LAMBDA_SSH_KEY contains raw key material, materialises it to a temp file
    and registers the public key with Lambda if not already present.
    """
    raw = ssh_key_env or os.environ.get("LAMBDA_SSH_KEY", "")

    if raw.startswith("-----BEGIN"):
        # Raw private key — write to temp file
        key_path = Path.home() / ".ssh" / "cursor_lambda_key"
        key_path.parent.mkdir(exist_ok=True)
        key_path.write_text(raw.strip() + "\n")
        key_path.chmod(0o600)

        # Register public key with Lambda if not already there
        pub = subprocess.run(["ssh-keygen", "-y", "-f", str(key_path)],
                             capture_output=True, text=True)
        if pub.returncode == 0:
            pub_key = pub.stdout.strip()
            existing = {k["name"] for k in client.list_ssh_keys()}
            if key_name not in existing:
                client.add_ssh_key(key_name, pub_key)
                print(f"  Registered SSH key '{key_name}' with Lambda")
            else:
                print(f"  Using inline SSH key material (registered as '{key_name}').")
        return str(key_path), key_name

    # Treat as a file path
    return os.path.expanduser(raw or f"~/.ssh/{key_name}"), key_name
```

---

## 5. Launch + Wait for SSH

```python
import time

# Preferred instance types (cheapest capable GPU first)
INSTANCE_PREFERENCE = [
    "gpu_1x_a100_sxm4", "gpu_1x_a100", "gpu_1x_a100_80gb_sxm4",
    "gpu_1x_h100_pcie", "gpu_1x_a6000", "gpu_1x_a10",
]
# ARM instances can't run x86 Docker images — skip them
ARM_INSTANCE_TYPES = {"gpu_1x_gh200"}

def launch_with_retry(client: LambdaClient, key_name: str,
                      instance_name: str = "agent") -> str:
    """Try every available (type, region) pair until one launches. Returns instance id."""
    from requests.exceptions import HTTPError
    avail = {t["name"]: t for t in client.available_instance_types()
             if t["name"] not in ARM_INSTANCE_TYPES}
    ordered = [n for n in INSTANCE_PREFERENCE if n in avail]
    ordered += sorted(n for n in avail if n not in INSTANCE_PREFERENCE)
    candidates = [(itype, region, avail[itype]["price_per_hour"])
                  for itype in ordered
                  for region in avail[itype]["available_regions"]]
    for itype, region, price in candidates:
        print(f"  Trying {itype} in {region} (${price:.2f}/h)...")
        try:
            iid = client.launch(itype, region,
                                ssh_key_names=[key_name],
                                name=instance_name)["id"]
            print(f"  Launched: {iid}")
            return iid
        except (HTTPError, Exception) as e:
            print(f"  Failed ({e}), trying next...")
            time.sleep(2)
    raise RuntimeError("No instance launched — all regions at capacity")


def wait_until_reachable(client: LambdaClient, instance_id: str,
                         key: str, timeout_s: int = 1800) -> str:
    """Block until SSH responds. Returns the instance IP."""
    print(f"  Waiting for SSH", end="", flush=True)
    deadline = time.time() + timeout_s
    ip = ""
    while time.time() < deadline:
        inst = client.get_instance(instance_id)
        status = inst.get("status", "unknown")
        if status in ("terminated", "terminating", "unhealthy"):
            raise RuntimeError(f"Instance entered state '{status}'")
        ip = inst.get("ip") or ""
        if ip and ssh(ip, key, "echo ok", check=False, quiet=True) == 0:
            print(f"\n  Reachable at {ip} (status: {status})")
            return ip
        print("." if ip else "_", end="", flush=True)
        time.sleep(15)
    raise TimeoutError(f"Instance {instance_id} did not accept SSH in {timeout_s}s")
```

---

## 6. Minimal End-to-End Template

```python
import os, subprocess, sys, time
from pathlib import Path

# -- paste lambda_client.py and helpers above, then: --

client   = LambdaClient(os.environ["LAMBDA_API_KEY"])
ssh_key, key_name = resolve_ssh_key(client)

iid = launch_with_retry(client, key_name, instance_name="my-job")
ip  = wait_until_reachable(client, iid, ssh_key)

try:
    # Upload files
    ssh(ip, ssh_key, "mkdir -p /home/ubuntu/pipeline")
    scp("my_script.sh", ip, "/home/ubuntu/pipeline/run.sh", ssh_key)
    ssh(ip, ssh_key, "chmod +x /home/ubuntu/pipeline/run.sh")

    # Run in tmux (survives SSH disconnects)
    ssh(ip, ssh_key,
        "tmux new-session -d -s job -c /home/ubuntu/pipeline "
        "'bash /home/ubuntu/pipeline/run.sh 2>&1 | tee /home/ubuntu/pipeline/run.log'")

    # Poll until done
    while True:
        done = ssh(ip, ssh_key,
                   "grep -c 'ALL DONE' /home/ubuntu/pipeline/run.log || true",
                   check=False, quiet=True)
        # ... your polling logic
        time.sleep(60)

    # Fetch results
    scp_from(ip, "/home/ubuntu/pipeline/results.tar.gz", ".", ssh_key)

finally:
    try:
        client.terminate(iid)
        print(f"Terminated {iid}")
    except Exception as e:
        print(f"WARNING: terminate failed — {iid} may still be running: {e}")
```

---

## 7. Installing Boltz-2 on Lambda (Known Issues)

Lambda instances have a system Python with many pre-installed packages that conflict with pip. **Always use a venv**:

```bash
# In your GPU run script (bash):
BOLTZ_VENV="$HOME/boltz_venv"
BOLTZ="$BOLTZ_VENV/bin/boltz"

if [ ! -f "$BOLTZ" ]; then
    python3 -m venv "$BOLTZ_VENV"
    PIP="$BOLTZ_VENV/bin/pip"
    "$PIP" install -q --upgrade pip
    # Install PyTorch for CUDA 12.4 — compatible with Lambda's NVIDIA driver
    "$PIP" install -q torch torchvision \
        --index-url https://download.pytorch.org/whl/cu124
    # Install boltz base (no [cuda] extras)
    "$PIP" install -q boltz
    # cuequivariance_ops_torch (compiled CUDA ext) is unavailable on Lambda →
    # remove it so boltz falls back to the pure-PyTorch triangular-mult path
    "$PIP" uninstall -y cuequivariance-torch cuequivariance 2>/dev/null || true
fi

# Run predictions (no --accelerator flag — Boltz-2 uses --devices only)
"$BOLTZ" predict /path/to/yaml_dir \
    --out_dir /path/to/out \
    --devices 1 \
    --num_workers 2 \
    --override
```

**Key gotchas:**
- Do NOT use `--accelerator gpu` — that is a Boltz-1 / PyTorch Lightning flag; Boltz-2 ignores it and may crash
- Do NOT use `boltz[cuda]` — it installs `cuequivariance-torch` whose compiled CUDA extension (`cuequivariance_ops_torch`) is not available on Lambda; install `boltz` then uninstall `cuequivariance*`
- `torchvision` from the system Python conflicts with pip-installed torch — the venv isolates this

---

## 8. Boltz-2 YAML Format

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: "RECEPTOR_SEQUENCE_HERE"
      msa: empty                        # skip MSA for speed
      templates:                        # optional: provide PDB as structural template
        - cif: /path/to/template.cif
          chain_id: A
  - protein:
      id: B
      sequence: "PEPTIDE_OR_BINDER"
      msa: empty
```

Confidence JSON output location: `<out_dir>/<input_name>/predictions/confidence_<name>_model_0.json`

Key fields in the JSON:
- `iptm`: dict of chain-pair scores e.g. `{"AB": 0.72, "BA": 0.72}` — use `AB` for interface quality
- `ptm`: overall predicted TM-score
- `mean_plddt`: average per-residue confidence

---

## 9. ProteinMPNN via Foundry Docker

For sequence design on a multi-chain PDB (design one chain, fix others), use the `rosettacommons/foundry` Docker image which has `atomworks` + `mpnn` pre-installed:

```bash
sudo chmod 666 /var/run/docker.sock
docker pull rosettacommons/foundry:latest

# Write a Python script to /home/ubuntu/pipeline/run_mpnn.py, then:
docker run --rm --gpus all \
    -v /home/ubuntu/pipeline:/workspace \
    rosettacommons/foundry:latest \
    python3 /workspace/run_mpnn.py
```

Inside `run_mpnn.py`:
```python
from atomworks.io.utils.io_utils import load_any
from mpnn.inference_engines.mpnn import MPNNInferenceEngine

engine = MPNNInferenceEngine(model_type="protein_mpnn", is_legacy_weights=True,
                             out_directory=None, write_structures=False, write_fasta=False)
raw = load_any("/workspace/top_structures/my_design.pdb")
aa  = raw[0] if hasattr(raw, "__getitem__") else raw

# Design chain M, fix A/B/C/D
chain_m_res = sorted(set(aa.res_id[aa.chain_id == "M"].tolist()))
designed    = [f"M{r}" for r in chain_m_res]

result = engine.run(atom_arrays=[aa], input_dicts=[{
    "batch_size": 8,
    "designed_residues": designed,
    "remove_waters": True,
}])
```

---

## 10. Available Instance Types (approximate, check API for current pricing)

| Type | GPU | VRAM | $/hr |
|---|---|---|---|
| `gpu_1x_a10` | A10 | 24 GB | ~$1.29 |
| `gpu_1x_a6000` | A6000 | 48 GB | ~$1.09 |
| `gpu_1x_a100_sxm4` | A100 SXM4 | 40 GB | ~$1.99 |
| `gpu_1x_a100_80gb_sxm4` | A100 80GB | 80 GB | ~$2.49 |
| `gpu_1x_h100_pcie` | H100 PCIe | 80 GB | ~$3.29 |
| `gpu_2x_a100` | 2× A100 | 80 GB | ~$3.99 |

Always use `launch_with_retry` — the API reports availability but launch can still fail (capacity race). The retry logic tries all (type, region) pairs in preference order.
