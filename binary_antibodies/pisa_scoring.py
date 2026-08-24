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

INTERFACE_SCOPE_FV = "fv"
INTERFACE_SCOPE_FULL_FAB = "full_fab"
INTERFACE_SCOPE_CHAIN_PAIRS: dict[str, tuple[tuple[str, str], ...]] = {
    INTERFACE_SCOPE_FV: (("A", "B"),),
    INTERFACE_SCOPE_FULL_FAB: (("A", "B"), ("C", "D")),
}
INTERFACE_PAIR_LABELS: dict[tuple[str, str], str] = {
    ("A", "B"): "vh_vl",
    ("C", "D"): "ch1_cl",
}


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
            residues = []
            for res in mol.findall(".//residue"):
                seq_num = _safe_int(res.findtext("seq_num"))
                if seq_num is None:
                    continue
                residues.append(
                    {
                        "seq_num": seq_num,
                        "name": (res.findtext("name") or "").strip(),
                        "bsa_A2": _safe_float(res.findtext("bsa")),
                        "asa_A2": _safe_float(res.findtext("asa")),
                        "solv_en_kcal": _safe_float(res.findtext("solv_en")),
                    }
                )
            molecules.append(
                {
                    "chain_id": chain_id,
                    "int_nres": _safe_int(mol.findtext("int_nres")),
                    "int_area": _safe_float(mol.findtext("int_area")),
                    "int_solv_en": _safe_float(mol.findtext("int_solv_en")),
                    "residues": residues,
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


def buried_interface_residues(
    iface: dict[str, Any],
    chain_a: str,
    chain_b: str,
    *,
    min_buried_sasa_A2: float = 0.0,
    max_resnum: dict[str, int] | None = None,
) -> dict[str, list[int]]:
    """
    Residue numbers per chain with PISA buried surface area (BSA) at the interface.

    Uses the molecule/residue lists from a PISA interface record.
    """
    out: dict[str, set[int]] = {chain_a: set(), chain_b: set()}
    want = {chain_a, chain_b}
    for mol in iface.get("molecules", []):
        cid = mol.get("chain_id")
        if cid not in want:
            continue
        for res in mol.get("residues", []):
            seq_num = res.get("seq_num")
            bsa = res.get("bsa_A2")
            if seq_num is None or bsa is None or bsa <= min_buried_sasa_A2:
                continue
            if max_resnum and cid in max_resnum and seq_num > max_resnum[cid]:
                continue
            out[cid].add(int(seq_num))
    return {chain: sorted(nums) for chain, nums in out.items()}


def residue_buried_sasa(
    iface: dict[str, Any],
    chain_id: str,
) -> dict[int, float]:
    """Map residue number → buried SASA (Å²) for one chain at the interface."""
    bsa: dict[int, float] = {}
    for mol in iface.get("molecules", []):
        if mol.get("chain_id") != chain_id:
            continue
        for res in mol.get("residues", []):
            seq_num = res.get("seq_num")
            val = res.get("bsa_A2")
            if seq_num is not None and val is not None and val > 0:
                bsa[int(seq_num)] = float(val)
    return bsa


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


def extract_chains_structure(
    structure_path: Path,
    out_path: Path,
    chains: tuple[str, ...] = ("A", "B"),
) -> Path:
    """Write a PDB containing only the requested chains."""
    return extract_fv_structure(structure_path, out_path, chains=chains)


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


def run_pisa_interfaces_xml(
    structure_path: Path,
    *,
    work_dir: Path,
    score_fv_only: bool = True,
    fv_chains: tuple[str, str] = ("A", "B"),
    analysis_chains: tuple[str, ...] | None = None,
    use_docker: bool = True,
    docker_image: str = DEFAULT_PISA_IMAGE,
    pisa_cfg: str = DEFAULT_PISA_CFG,
    pisa_binary: str = "pisa",
    session_name: str | None = None,
) -> str:
    """Run PISA --as-is and return interfaces XML text."""
    work_dir.mkdir(parents=True, exist_ok=True)
    input_path = structure_path
    chains = analysis_chains or (fv_chains if score_fv_only else ("A", "B", "C", "D"))
    if score_fv_only or analysis_chains is not None:
        input_path = work_dir / f"{structure_path.stem}_{'_'.join(chains)}.pdb"
        extract_chains_structure(structure_path, input_path, chains=chains)

    session = _sanitize_session_name(session_name or structure_path.stem)
    if use_docker and shutil.which("docker"):
        return _run_pisa_docker(session, input_path, work_dir, pisa_cfg, docker_image)
    if shutil.which(pisa_binary):
        return _run_pisa_local(session, input_path.resolve(), work_dir, pisa_cfg, pisa_binary)
    raise RuntimeError("PISA unavailable: docker and local pisa binary not found")


def identify_fab_interface_residues(
    structure_path: Path,
    *,
    work_dir: Path,
    interface_scope: str = INTERFACE_SCOPE_FV,
    chain_lengths: dict[str, int] | None = None,
    min_buried_sasa_A2: float = 0.0,
    exclude_cdr_resnums: dict[str, set[int]] | None = None,
    pisa_config: PisaConfig | None = None,
    session_name: str | None = None,
) -> dict[str, Any]:
    """
    Run PISA on WT Fab and return buried framework interface residues per chain.

    interface_scope:
        ``fv`` — VH–VL only (chains A+B).
        ``full_fab`` — VH–VL (A+B) and CH1–CL (C+D).
    """
    if interface_scope not in INTERFACE_SCOPE_CHAIN_PAIRS:
        return {
            "pisa_status": "invalid_scope",
            "pisa_error": f"Unknown interface_scope: {interface_scope}",
        }

    chain_pairs = INTERFACE_SCOPE_CHAIN_PAIRS[interface_scope]
    chain_lengths = chain_lengths or {"A": 113, "B": 107, "C": 101, "D": 113}
    analysis_chains = tuple(sorted({c for pair in chain_pairs for c in pair}))
    cfg = pisa_config or PisaConfig()

    try:
        xml_text = run_pisa_interfaces_xml(
            structure_path,
            work_dir=work_dir,
            score_fv_only=True,
            analysis_chains=analysis_chains,
            use_docker=cfg.use_docker,
            docker_image=cfg.docker_image,
            pisa_cfg=cfg.pisa_cfg,
            pisa_binary=cfg.pisa_binary,
            session_name=session_name or f"wt_{interface_scope}_interface",
        )
    except RuntimeError as exc:
        return {"pisa_status": "unavailable", "pisa_error": str(exc)}

    try:
        interfaces = parse_pisa_interfaces_xml(xml_text)
    except ET.ParseError as exc:
        return {"pisa_status": "parse_error", "pisa_error": str(exc)}

    framework: dict[str, list[int]] = {chain: [] for chain in analysis_chains}
    residue_bsa_by_chain: dict[str, dict[int, float]] = {chain: {} for chain in analysis_chains}
    interface_records: dict[str, dict[str, Any]] = {}
    buried_all: dict[str, list[int]] = {chain: [] for chain in analysis_chains}

    for chain_a, chain_b in chain_pairs:
        iface = select_chain_pair_interface(interfaces, chain_a, chain_b)
        label = INTERFACE_PAIR_LABELS.get((chain_a, chain_b), f"{chain_a}_{chain_b}")
        if iface is None:
            return {
                "pisa_status": f"no_{label}_interface",
                "pisa_n_interfaces": len(interfaces),
                "pisa_chain_pairs": [i.get("chain_ids") for i in interfaces],
            }

        max_resnum = {chain_a: chain_lengths[chain_a], chain_b: chain_lengths[chain_b]}
        buried = buried_interface_residues(
            iface,
            chain_a,
            chain_b,
            min_buried_sasa_A2=min_buried_sasa_A2,
            max_resnum=max_resnum,
        )
        for chain in (chain_a, chain_b):
            buried_all[chain] = sorted(set(buried_all.get(chain, [])) | set(buried.get(chain, [])))

        for chain in (chain_a, chain_b):
            nums = buried.get(chain, [])
            if exclude_cdr_resnums and chain in exclude_cdr_resnums:
                skip = exclude_cdr_resnums[chain]
                nums = [r for r in nums if r not in skip]
            framework[chain] = sorted(set(framework.get(chain, [])) | set(nums))
            residue_bsa_by_chain[chain].update(
                {
                    r: v
                    for r, v in residue_buried_sasa(iface, chain).items()
                    if r in framework[chain]
                }
            )

        interface_records[label] = {
            "chain_ids": [chain_a, chain_b],
            "pisa_int_area_A2": iface.get("int_area_A2"),
            "pisa_int_solv_en_kcal": iface.get("int_solv_en_kcal"),
            "pisa_interface_id": iface.get("id"),
            "buried_interface_residues": buried,
        }

    return {
        "pisa_status": "ok",
        "interface_scope": interface_scope,
        "analysis_chains": list(analysis_chains),
        "interfaces": interface_records,
        "buried_interface_residues": buried_all,
        "framework_interface_residues": framework,
        "residue_buried_sasa_A2": residue_bsa_by_chain,
        "pisa_xml_text": xml_text,
    }


def identify_fv_interface_residues(
    structure_path: Path,
    *,
    work_dir: Path,
    chain_a: str = "A",
    chain_b: str = "B",
    vh_len: int = 113,
    vl_len: int = 107,
    min_buried_sasa_A2: float = 0.0,
    exclude_cdr_resnums: dict[str, set[int]] | None = None,
    pisa_config: PisaConfig | None = None,
    session_name: str | None = None,
    interface_scope: str = INTERFACE_SCOPE_FV,
    ch1_len: int = 101,
    cl_len: int = 113,
) -> dict[str, Any]:
    """
    Run PISA on a WT Fab structure and return buried interface residues per chain.

    For ``interface_scope='full_fab'``, delegates to :func:`identify_fab_interface_residues`
    and returns CH1/CL interface residues in addition to VH/VL.
    """
    if interface_scope == INTERFACE_SCOPE_FULL_FAB:
        return identify_fab_interface_residues(
            structure_path,
            work_dir=work_dir,
            interface_scope=interface_scope,
            chain_lengths={"A": vh_len, "B": vl_len, "C": ch1_len, "D": cl_len},
            min_buried_sasa_A2=min_buried_sasa_A2,
            exclude_cdr_resnums=exclude_cdr_resnums,
            pisa_config=pisa_config,
            session_name=session_name,
        )

    result = identify_fab_interface_residues(
        structure_path,
        work_dir=work_dir,
        interface_scope=INTERFACE_SCOPE_FV,
        chain_lengths={chain_a: vh_len, chain_b: vl_len},
        min_buried_sasa_A2=min_buried_sasa_A2,
        exclude_cdr_resnums=exclude_cdr_resnums,
        pisa_config=pisa_config,
        session_name=session_name or "wt_fv_interface",
    )
    if result.get("pisa_status") != "ok":
        return result

    vh_vl = result.get("interfaces", {}).get("vh_vl", {})
    framework = result.get("framework_interface_residues", {})
    bsa = result.get("residue_buried_sasa_A2", {})
    return {
        "pisa_status": "ok",
        "pisa_int_area_A2": vh_vl.get("pisa_int_area_A2"),
        "pisa_int_solv_en_kcal": vh_vl.get("pisa_int_solv_en_kcal"),
        "pisa_interface_id": vh_vl.get("pisa_interface_id"),
        "buried_interface_residues": vh_vl.get("buried_interface_residues", {}),
        "framework_interface_residues": {
            chain_a: framework.get(chain_a, []),
            chain_b: framework.get(chain_b, []),
        },
        "residue_buried_sasa_A2": {
            chain_a: bsa.get(chain_a, {}),
            chain_b: bsa.get(chain_b, {}),
        },
        "pisa_xml_text": result.get("pisa_xml_text", ""),
    }


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
    try:
        xml_text = run_pisa_interfaces_xml(
            structure_path,
            work_dir=work_dir,
            score_fv_only=score_fv_only,
            fv_chains=(chain_a, chain_b),
            use_docker=use_docker,
            docker_image=docker_image,
            pisa_cfg=pisa_cfg,
            pisa_binary=pisa_binary,
            session_name=session_name,
        )
    except RuntimeError as exc:
        return {"pisa_status": "unavailable", "pisa_error": str(exc)}

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
