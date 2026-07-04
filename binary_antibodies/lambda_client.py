"""
lambda_client.py
----------------
Thin wrapper around the Lambda Cloud REST API for launching and managing
GPU instances to run the RFdiffusion / ProteinMPNN / AF2 pipeline.

Usage
-----
    from binary_antibodies.lambda_client import LambdaClient
    client = LambdaClient(api_key=os.environ["LAMBDA_API_KEY"])
    instance = client.launch("gpu_1x_a100_sxm4", "us-east-1", ssh_key_name="mactoby")
    client.wait_until_active(instance["id"])
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass
from typing import Any

import requests

BASE_URL = "https://cloud.lambda.ai/api/v1"
POLL_INTERVAL_S = 15


class LambdaAPIError(Exception):
    pass


class LambdaClient:
    """Lambda Cloud API client."""

    def __init__(self, api_key: str | None = None) -> None:
        self.api_key = api_key or os.environ.get("LAMBDA_API_KEY", "")
        if not self.api_key:
            raise ValueError(
                "Lambda Cloud API key required. "
                "Set LAMBDA_API_KEY environment variable or pass api_key=."
            )
        self._session = requests.Session()
        self._session.headers.update({"Authorization": f"Bearer {self.api_key}"})

    # ------------------------------------------------------------------
    # Low-level helpers
    # ------------------------------------------------------------------

    def _get(self, path: str) -> dict:
        resp = self._session.get(f"{BASE_URL}{path}", timeout=30)
        resp.raise_for_status()
        data = resp.json()
        if "error" in data:
            raise LambdaAPIError(data["error"]["message"])
        return data

    def _post(self, path: str, payload: dict) -> dict:
        resp = self._session.post(f"{BASE_URL}{path}", json=payload, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        if "error" in data:
            raise LambdaAPIError(data["error"]["message"])
        return data

    def _delete(self, path: str, payload: dict | None = None) -> dict:
        resp = self._session.delete(
            f"{BASE_URL}{path}", json=payload or {}, timeout=30
        )
        resp.raise_for_status()
        data = resp.json()
        if "error" in data:
            raise LambdaAPIError(data["error"]["message"])
        return data

    # ------------------------------------------------------------------
    # Instance types
    # ------------------------------------------------------------------

    def list_instance_types(self) -> dict[str, dict]:
        """Return all instance types with availability info."""
        return self._get("/instance-types")["data"]

    def available_instance_types(self) -> list[dict]:
        """Return only types with current capacity."""
        all_types = self.list_instance_types()
        available = []
        for name, info in sorted(all_types.items()):
            regions = [r["name"] for r in info.get("regions_with_capacity_available", [])]
            if regions:
                it = info.get("instance_type", {})
                specs = it.get("specs", {})
                available.append({
                    "name": name,
                    "price_per_hour": it.get("price_cents_per_hour", 0) / 100,
                    "gpus": specs.get("gpus", 0),
                    "vcpus": specs.get("vcpus", 0),
                    "memory_gib": specs.get("memory_gibibytes", 0),
                    "storage_gib": specs.get("storage_gibibytes", 0),
                    "available_regions": regions,
                })
        return available

    def print_available(self) -> None:
        """Print a table of available instance types."""
        types = self.available_instance_types()
        print(f"{'Name':<40} {'$/hr':>6}  {'GPUs':>4}  {'vCPUs':>5}  {'RAM(GiB)':>8}  Regions")
        print("-" * 90)
        for t in types:
            print(
                f"{t['name']:<40} {t['price_per_hour']:>6.2f}"
                f"  {t['gpus']:>4}  {t['vcpus']:>5}  {t['memory_gib']:>8}"
                f"  {', '.join(t['available_regions'])}"
            )

    # ------------------------------------------------------------------
    # Instances
    # ------------------------------------------------------------------

    def list_instances(self) -> list[dict]:
        """List all active instances."""
        return self._get("/instances")["data"]

    def get_instance(self, instance_id: str) -> dict:
        """Get details for a specific instance."""
        return self._get(f"/instances/{instance_id}")["data"]

    def launch(
        self,
        instance_type: str,
        region: str,
        ssh_key_names: list[str],
        name: str = "rfdiffusion-minibinder",
        file_system_names: list[str] | None = None,
        user_data: str = "",
    ) -> dict:
        """
        Launch a new instance.

        Parameters
        ----------
        instance_type : str
            e.g. 'gpu_1x_a100_sxm4'
        region : str
            e.g. 'us-east-1'
        ssh_key_names : list[str]
            SSH key names already registered in your Lambda account.
        name : str
            Human-readable instance name.
        file_system_names : list[str], optional
            Persistent file systems to attach.
        user_data : str
            Cloud-init script to run on first boot.

        Returns
        -------
        dict
            Instance info including 'id' and 'ip'.
        """
        payload: dict[str, Any] = {
            "region_name": region,
            "instance_type_name": instance_type,
            "ssh_key_names": ssh_key_names,
            "name": name,
        }
        if file_system_names:
            payload["file_system_names"] = file_system_names
        if user_data:
            payload["user_data"] = user_data

        resp = self._post("/instance-operations/launch", payload)
        instance_ids = resp.get("data", {}).get("instance_ids", [])
        if not instance_ids:
            raise LambdaAPIError(f"Launch failed: {resp}")
        instance_id = instance_ids[0]
        print(f"  Launched instance: {instance_id}")
        return {"id": instance_id}

    def terminate(self, instance_id: str) -> None:
        """Terminate an instance."""
        self._post("/instance-operations/terminate", {"instance_ids": [instance_id]})
        print(f"  Terminated: {instance_id}")

    def wait_until_active(
        self,
        instance_id: str,
        timeout_s: int = 600,
        poll_s: int = POLL_INTERVAL_S,
    ) -> dict:
        """
        Poll until instance is in 'active' state or timeout.

        Returns the instance dict (with 'ip' populated).
        """
        print(f"  Waiting for instance {instance_id} to become active", end="", flush=True)
        deadline = time.time() + timeout_s
        while time.time() < deadline:
            inst = self.get_instance(instance_id)
            status = inst.get("status", "unknown")
            ip = inst.get("ip")
            if status == "active" and ip:
                print(f"\n  Active! IP: {ip}")
                return inst
            if status in ("terminated", "terminating", "unhealthy"):
                raise LambdaAPIError(f"Instance entered unexpected status: {status}")
            print(".", end="", flush=True)
            time.sleep(poll_s)
        raise TimeoutError(f"Instance {instance_id} did not become active within {timeout_s}s")

    # ------------------------------------------------------------------
    # SSH keys
    # ------------------------------------------------------------------

    def list_ssh_keys(self) -> list[dict]:
        return self._get("/ssh-keys")["data"]

    def add_ssh_key(self, name: str, public_key: str) -> dict:
        return self._post("/ssh-keys", {"name": name, "public_key": public_key})["data"]

    # ------------------------------------------------------------------
    # File systems
    # ------------------------------------------------------------------

    def list_file_systems(self) -> list[dict]:
        return self._get("/file-systems")["data"]
