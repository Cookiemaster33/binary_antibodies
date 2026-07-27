"""
pisa_scoring.py
---------------
VH–VL interface energetics via PDBe-KB PISA (Docker or local binary).

Runs PISA in --as-is mode on Fv chains (A+B) extracted from holo Boltz structures
and reports interface area and solvation energy for ranking weakened designs.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any

DEFAULT_PISA_IMAGE = "pdbegroup/pisa:latest"
DEFAULT_PISA_CFG = "/usr/share/pisa/setup/pisa_cfg_tmp"


@dataclass(frozen=True)
class PisaConfig:
    enabled: bool = True
    docker_image: str = DEFAULT_PISA_IMAGE
    pisa_cfg: str = DEFAULT_PISA_CFG
    fv_chains: tuple[str, str] = ("A", "B")
    score_fv_only: bool = True
    use_docker: bool = True
    pisa_binary: str = "pisa"
    work_subdir: str = "pisa_work"


def load_pisa_config(config_path: Path | None) -> PisaConfig:
    if not config_path or not config_path.exists():
        return PisaConfig()
    cfg = __import__("json").loads(config_path.read_text())
    pisa = cfg.get("validation", {}).get("pisa", {})
    fv = pisa.get("fv_chains", ["A", "B"])
    return PisaConfig(
        enabled=bool(pisa.get("enabled", True)),
        docker_image=str(pisa.get("docker_image", DEFAULT_PISA_IMAGE)),
        pisa_cfg=str(pisa.get("pisa_cfg", DEFAULT_PISA_CFG)),
        fv_chains=(str(fv[0]), str(fv[1])),
        score_fv_only=bool(pisa.get("score_fv_only", True)),
        use_docker=bool(pisa.get("use_docker", True)),
        pisa_binary=str(pisa.get("pisa_binary", "pisa")),
        work_subdir=str(pisa.get("work_subdir", "pisa_work")),
    )


def _safe_float(text: str | None) -> float | None:
    if text is None:
        return None
    text = text.strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _safe_int(text: str | None) -> int | None:
    if text is None:
        return None
    text = text.strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def _parse_bond_count(section: ET.Element | None) -> int:
    if section is None:
        return 0
    return _safe_int(section.findtext("n_bonds")) or len(section.findall("bond"))


def parse_pisa_interfaces_xml(xml_text: str) -> list[dict[str, Any]]:
    """Parse PISA interfaces XML into a list of interface records."""
    root = ET.fromstring(xml_text)
    interfaces: list[dict[str, Any]] = []

    iface_nodes = root.findall(".//interface")
    if not iface_nodes and root.tag == "interface":
        iface_nodes = [root]

    for iface_el in iface_nodes:
        molecules = []
        chain_ids: set[str] = set()
        for mol in iface_el.findall("molecule"):
            chain_id = (mol.findtext("chain_id") or "").strip()
            if chain_id:
                chain_ids.add(chain_id)
            molecules.append(
                {
                    "chain_id": chain_id,
                    "int_nres": _safe_int(mol.findtext("int_nres")),
                    "int_area": _safe_float(mol.findtext("int_area")),
                    "int_solv_en": _safe_float(mol.findtext("int_solv_en")),
                }
            )

        solv_en = _safe_float(iface_el.findtext("int_solv_en"))
        if solv_en is None:
            solv_en = _safe_float(iface_el.findtext("Int_solv_en"))

        interfaces.append(
            {
                "id": _safe_int(iface_el.findtext("id")),
                "chain_ids": sorted(chain_ids),
                "int_area_A2": _safe_float(iface_el.findtext("int_area")),
                "int_solv_en_kcal": solv_en,
                "stab_en_kcal": _safe_float(iface_el.findtext("stab_en")),
                "pvalue": _safe_float(iface_el.findtext("pvalue")),
                "n_h_bonds": _parse_bond_count(iface_el.find("h-bonds")),
                "n_salt_bridges": _parse_bond_count(iface_el.find("salt-bridges")),
                "n_ss_bonds": _parse_bond_count(iface_el.find("ss-bonds")),
                "molecules": molecules,
            }
        )
    return interfaces


def select_chain_pair_interface(
    interfaces: list[dict[str, Any]],
    chain_a: str,
    chain_b: str,
) -> dict[str, Any] | None:
    """Return the interface between chain_a and chain_b (order-independent)."""
    want = {chain_a, chain_b}
    matches = [iface for iface in interfaces if set(iface.get("chain_ids", [])) == want]
    if not matches:
        # Fall back: interface whose chain set is a superset containing both chains.
        matches = [
            iface
            for iface in interfaces
            if want.issubset(set(iface.get("chain_ids", [])))
        ]
    if not matches:
        return None
    # Prefer the smallest-area match if multiple (e.g. duplicated crystal contacts).
    return min(matches, key=lambda iface: iface.get("int_area_A2") or float("inf"))


def extract_fv_structure(
    structure_path: Path,
    out_path: Path,
    chains: tuple[str, ...] = ("A", "B"),
) -> Path:
    """Write a PDB containing only the requested chains."""
    import numpy as np
    from biotite.structure.io.pdb import PDBFile
    from biotite.structure.io.pdbx import CIFFile, get_structure

    keep = set(chains)
    if structure_path.suffix.lower() in {".cif", ".mmcif"}:
        cif = CIFFile.read(str(structure_path))
        aa = get_structure(cif, model=1, include_bonds=False)
    else:
        pdb = PDBFile.read(str(structure_path))
        aa = pdb.get_structure(model=1)

    sub = aa[np.isin(aa.chain_id, list(keep))]
    if len(sub) == 0:
        raise ValueError(f"No atoms for chains {chains} in {structure_path}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_out = PDBFile()
    pdb_out.set_structure(sub)
    pdb_out.write(str(out_path))
    return out_path


def _sanitize_session_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)[:80]


def _run_pisa_local(
    session: str,
    structure_path: Path,
    work_dir: Path,
    cfg: str,
    binary: str,
) -> str:
    analyse_cmd = [binary, session, "-analyse", str(structure_path), "--as-is", cfg]
    xml_cmd = [binary, session, "-xml", "interfaces", "--as-is", cfg]
    subprocess.run(analyse_cmd, cwd=work_dir, check=True, capture_output=True, text=True)
    proc = subprocess.run(xml_cmd, cwd=work_dir, check=True, capture_output=True, text=True)
    return proc.stdout


def _run_pisa_docker(
    session: str,
    structure_path: Path,
    work_dir: Path,
    cfg: str,
    image: str,
) -> str:
    rel = structure_path.name
    mount = f"{work_dir.resolve()}:/data"
    in_container = f"/data/{rel}"
    analyse_cmd = [
        "docker",
        "run",
        "--rm",
        "-v",
        mount,
        image,
        "pisa",
        session,
        "-analyse",
        in_container,
        "--as-is",
        cfg,
    ]
    xml_cmd = [
        "docker",
        "run",
        "--rm",
        "-v",
        mount,
        image,
        "pisa",
        session,
        "-xml",
        "interfaces",
        "--as-is",
        cfg,
    ]
    subprocess.run(analyse_cmd, check=True, capture_output=True, text=True)
    proc = subprocess.run(xml_cmd, check=True, capture_output=True, text=True)
    return proc.stdout


def score_vh_vl_pisa(
    structure_path: Path,
    *,
    work_dir: Path,
    chain_a: str = "A",
    chain_b: str = "B",
    score_fv_only: bool = True,
    use_docker: bool = True,
    docker_image: str = DEFAULT_PISA_IMAGE,
    pisa_cfg: str = DEFAULT_PISA_CFG,
    pisa_binary: str = "pisa",
    session_name: str | None = None,
) -> dict[str, Any]:
    """
    Run PISA on a structure and return VH–VL (chain_a vs chain_b) interface metrics.

    Negative int_solv_en_kcal indicates a favorable interface; weaker designs are
    less negative (higher) than native WT.
    """
    import numpy as np  # noqa: F401 — biotite dependency

    work_dir.mkdir(parents=True, exist_ok=True)
    input_path = structure_path
    if score_fv_only:
        input_path = work_dir / f"{structure_path.stem}_fv_AB.pdb"
        extract_fv_structure(structure_path, input_path, chains=(chain_a, chain_b))

    session = _sanitize_session_name(session_name or structure_path.stem)

    if use_docker and shutil.which("docker"):
        xml_text = _run_pisa_docker(session, input_path, work_dir, pisa_cfg, docker_image)
    elif shutil.which(pisa_binary):
        xml_text = _run_pisa_local(session, input_path.resolve(), work_dir, pisa_cfg, pisa_binary)
    else:
        return {
            "pisa_status": "unavailable",
            "pisa_error": "docker and local pisa binary not found",
        }

    try:
        interfaces = parse_pisa_interfaces_xml(xml_text)
    except ET.ParseError as exc:
        return {"pisa_status": "parse_error", "pisa_error": str(exc)}

    iface = select_chain_pair_interface(interfaces, chain_a, chain_b)
    if iface is None:
        return {
            "pisa_status": "no_ab_interface",
            "pisa_n_interfaces": len(interfaces),
            "pisa_chain_pairs": [i.get("chain_ids") for i in interfaces],
        }

    return {
        "pisa_status": "ok",
        "pisa_int_area_A2": iface.get("int_area_A2"),
        "pisa_int_solv_en_kcal": iface.get("int_solv_en_kcal"),
        "pisa_stab_en_kcal": iface.get("stab_en_kcal"),
        "pisa_pvalue": iface.get("pvalue"),
        "pisa_n_h_bonds": iface.get("n_h_bonds"),
        "pisa_n_salt_bridges": iface.get("n_salt_bridges"),
        "pisa_n_ss_bonds": iface.get("n_ss_bonds"),
        "pisa_interface_id": iface.get("id"),
    }


def add_pisa_deltas(metrics: dict[str, Any], native: dict[str, Any]) -> dict[str, Any]:
    """Add delta metrics vs native reference PISA scores."""
    out = dict(metrics)
    native_solv = native.get("pisa_int_solv_en_kcal")
    design_solv = metrics.get("pisa_int_solv_en_kcal")
    if native_solv is not None and design_solv is not None:
        # Positive delta => design interface less favorable (weaker) than native.
        out["pisa_delta_int_solv_en_vs_native_kcal"] = round(design_solv - native_solv, 3)

    native_area = native.get("pisa_int_area_A2")
    design_area = metrics.get("pisa_int_area_A2")
    if native_area is not None and design_area is not None and native_area > 0:
        out["pisa_int_area_fraction_of_native"] = round(design_area / native_area, 4)
        out["pisa_delta_int_area_A2"] = round(design_area - native_area, 2)
    return out
