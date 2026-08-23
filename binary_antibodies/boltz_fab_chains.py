"""
boltz_fab_chains.py
-------------------
Fuse split Fab domains for Boltz holo validation.

MPNN / static scoring keep logical chains A=VH, B=VL, C=CH1, D=CL.
Boltz holo refold uses continuous immunoglobulin chains:
  H = VH + CH1
  L = VL + CL
  T = epitope stub (optional)
"""

from __future__ import annotations

from typing import Literal

import numpy as np

BOLTZ_HEAVY_ID = "H"
BOLTZ_LIGHT_ID = "L"
BOLTZ_EPITOPE_ID = "T"
BOLTZ_CHAIN_MODE_FUSED = "fused"
BOLTZ_CHAIN_MODE_SPLIT = "split"
BoltzChainMode = Literal["fused", "split"]


def default_boltz_chain_mode(interface_scope: str) -> BoltzChainMode:
    """full_fab validates with fused H/L; fv keeps legacy split chains."""
    return BOLTZ_CHAIN_MODE_FUSED if interface_scope == "full_fab" else BOLTZ_CHAIN_MODE_SPLIT


def load_boltz_chain_mode(config: dict | None, interface_scope: str) -> BoltzChainMode:
    if config:
        mode = config.get("validation", {}).get("boltz_chain_mode")
        if mode in (BOLTZ_CHAIN_MODE_FUSED, BOLTZ_CHAIN_MODE_SPLIT):
            return mode
    return default_boltz_chain_mode(interface_scope)


def fuse_fab_sequences(seq_a: str, seq_b: str, seq_c: str, seq_d: str) -> tuple[str, str]:
    """Concatenate VH+CH1 and VL+CL for Boltz."""
    return seq_a + seq_c, seq_b + seq_d


def boltz_holo_chain_sequences(
    chains: dict[str, str],
    *,
    mode: BoltzChainMode,
    epitope: str = "",
) -> list[tuple[str, str]]:
    """
    Return ordered (chain_id, sequence) entries for a Boltz YAML.

    fused: H (A+C), L (B+D), T
    split: A, B, C, D, T (legacy)
    """
    seq_a = chains.get("A", "")
    seq_b = chains.get("B", "")
    seq_c = chains.get("C", "")
    seq_d = chains.get("D", "")
    out: list[tuple[str, str]] = []

    if mode == BOLTZ_CHAIN_MODE_FUSED and seq_a and seq_b and seq_c and seq_d:
        heavy, light = fuse_fab_sequences(seq_a, seq_b, seq_c, seq_d)
        out.append((BOLTZ_HEAVY_ID, heavy))
        out.append((BOLTZ_LIGHT_ID, light))
    else:
        for ch in ("A", "B", "C", "D"):
            if chains.get(ch):
                out.append((ch, chains[ch]))

    if epitope:
        out.append((BOLTZ_EPITOPE_ID, epitope))
    return out


def format_boltz_yaml(chain_seqs: list[tuple[str, str]]) -> str:
    lines = ["sequences:"]
    for chain_id, sequence in chain_seqs:
        lines.append("  - protein:")
        lines.append(f"      id: {chain_id}")
        lines.append(f'      sequence: "{sequence}"')
        lines.append("      msa: empty")
    return "\n".join(lines) + "\n"


def is_fused_boltz_structure(aa, vh_len: int, vl_len: int) -> bool:
    """True when holo uses fused H/L (or long A/B) without separate C/D."""
    chains = set(aa.chain_id.tolist())
    if "C" in chains or "D" in chains:
        return False
    heavy_id = BOLTZ_HEAVY_ID if BOLTZ_HEAVY_ID in chains else "A"
    if heavy_id not in chains:
        return False
    heavy_res = {int(r) for r in aa.res_id[aa.chain_id == heavy_id]}
    return len(heavy_res) > vh_len


def _extract_relabel(
    aa,
    src_chain: str,
    res_lo: int,
    res_hi: int,
    dst_chain: str,
) -> np.ndarray:
    mask = (aa.chain_id == src_chain) & (aa.res_id >= res_lo) & (aa.res_id <= res_hi)
    sub = aa[mask]
    if len(sub) == 0:
        return sub
    sub = sub.copy()
    sub.chain_id = np.full(len(sub), dst_chain, dtype=sub.chain_id.dtype)
    sub.res_id = sub.res_id - res_lo + 1
    return sub


def expand_fused_fab_structure(aa, vh_len: int, vl_len: int):
    """
    Remap fused Boltz output (H/L[/T]) to logical chains A,B,C,D[,T] for scoring.

    VH residues on H → chain A (1..vh_len)
    CH1 residues on H → chain C (1..)
    VL residues on L → chain B (1..vl_len)
    CL residues on L → chain D (1..)
    """
    if not is_fused_boltz_structure(aa, vh_len, vl_len):
        return aa

    from biotite.structure import concatenate

    chains = set(aa.chain_id.tolist())
    heavy_id = BOLTZ_HEAVY_ID if BOLTZ_HEAVY_ID in chains else "A"
    light_id = BOLTZ_LIGHT_ID if BOLTZ_LIGHT_ID in chains else "B"

    heavy_res = sorted(int(r) for r in aa.res_id[aa.chain_id == heavy_id])
    light_res = sorted(int(r) for r in aa.res_id[aa.chain_id == light_id])
    if len(heavy_res) <= vh_len or len(light_res) <= vl_len:
        return aa

    parts = [
        _extract_relabel(aa, heavy_id, 1, vh_len, "A"),
        _extract_relabel(aa, heavy_id, vh_len + 1, heavy_res[-1], "C"),
        _extract_relabel(aa, light_id, 1, vl_len, "B"),
        _extract_relabel(aa, light_id, vl_len + 1, light_res[-1], "D"),
    ]
    if BOLTZ_EPITOPE_ID in chains:
        parts.append(aa[aa.chain_id == BOLTZ_EPITOPE_ID].copy())

    parts = [p for p in parts if len(p) > 0]
    if not parts:
        return aa
    return concatenate(parts)


def write_structure_pdb(aa, path: Path) -> None:
    """Write biotite AtomArray to PDB (used for PISA on expanded holo structures)."""
    from biotite.structure.io.pdb import PDBFile

    path.parent.mkdir(parents=True, exist_ok=True)
    pdb = PDBFile()
    pdb.set_structure(aa)
    pdb.write(str(path))
