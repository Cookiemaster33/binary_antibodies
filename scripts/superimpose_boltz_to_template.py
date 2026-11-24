#!/usr/bin/env python3
"""Kabsch-superimpose Boltz predicted CIFs onto CFTR template coordinates.

Predicted two-chain split numbering is 1-based within each chain. Template
auth_seq_id is UniProt P13569 numbering. We pair receptor Cα by UniProt
residue, compute the rigid transform, and apply it to EVERY atom so the
written file already overlaps the template in PyMOL/ChimeraX with no align.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

import numpy as np

# pred_res + offset = UniProt residue
TMD12_MAP = {
    # chain, pred_start, pred_end, offset
    "A": (1, 367, 27),    # CFTR 28-394
    "B": (1, 328, 844),   # CFTR 845-1172
}
FULLSPLIT_MAP = {
    "A": (1, 617, 27),    # CFTR 28-644
    "B": (1, 636, 844),   # CFTR 845-1480
}

PRESETS = {
    "tmd12": TMD12_MAP,
    "fullsplit": FULLSPLIT_MAP,
}


def tokenize_cif_line(line: str) -> list[str]:
    """Split a mmCIF data line, respecting quoted tokens."""
    tokens: list[str] = []
    i = 0
    n = len(line)
    while i < n:
        c = line[i]
        if c.isspace():
            i += 1
            continue
        if c in "'\"":
            quote = c
            i += 1
            start = i
            while i < n and line[i] != quote:
                i += 1
            tokens.append(line[start:i])
            i += 1
            continue
        start = i
        while i < n and not line[i].isspace():
            i += 1
        tokens.append(line[start:i])
    return tokens


def parse_atom_site(path: Path) -> tuple[list[str], list[str], list[list[str]]]:
    """Return (header_lines, column_names, atom_rows) for the first atom_site loop."""
    text = path.read_text()
    lines = text.splitlines(keepends=True)
    col_names: list[str] = []
    rows: list[list[str]] = []
    header: list[str] = []
    in_loop = False
    collecting_cols = False
    i = 0
    while i < len(lines):
        raw = lines[i]
        stripped = raw.strip()
        if not in_loop:
            if stripped == "loop_":
                # peek whether this loop is atom_site
                j = i + 1
                cols: list[str] = []
                while j < len(lines) and lines[j].lstrip().startswith("_"):
                    cols.append(lines[j].strip())
                    j += 1
                if cols and cols[0].startswith("_atom_site."):
                    in_loop = True
                    collecting_cols = True
                    col_names = []
                    i += 1
                    continue
            header.append(raw)
            i += 1
            continue
        if collecting_cols:
            if stripped.startswith("_atom_site."):
                col_names.append(stripped)
                i += 1
                continue
            collecting_cols = False
        if stripped.startswith("#") or stripped.startswith("loop_") or (
            stripped.startswith("_") and not stripped.startswith("_atom_site.")
        ):
            # remainder of file stays in header-after; we splice on write
            header.append(raw)
            i += 1
            while i < len(lines):
                header.append(lines[i])
                i += 1
            break
        if not stripped or stripped.startswith("#"):
            i += 1
            continue
        # data row (possibly the rest of a wrapped row — Boltz/PDB atom rows are one line)
        rows.append(tokenize_cif_line(stripped))
        i += 1
    if not col_names or not rows:
        raise ValueError(f"No atom_site loop in {path}")
    n = len(col_names)
    bad = [r for r in rows if len(r) != n]
    if bad:
        raise ValueError(
            f"{path}: {len(bad)} atom_site rows have {len(bad[0])} tokens, expected {n}"
        )
    return header, col_names, rows


def col_index(cols: list[str], *names: str) -> int:
    for name in names:
        if name in cols:
            return cols.index(name)
        bare = name if name.startswith("_") else f"_atom_site.{name}"
        if bare in cols:
            return cols.index(bare)
    raise KeyError(names)


def receptor_ca_by_uniprot(
    cols: list[str],
    rows: list[list[str]],
    *,
    is_prediction: bool,
    mapping: dict[str, tuple[int, int, int]],
    template_chain: str = "A",
) -> dict[int, np.ndarray]:
    """Map UniProt residue number -> Cα xyz for receptor chains only."""
    i_atom = col_index(cols, "_atom_site.label_atom_id")
    i_asym = col_index(cols, "_atom_site.label_asym_id")
    i_auth_asym = col_index(cols, "_atom_site.auth_asym_id") if "_atom_site.auth_asym_id" in cols else i_asym
    i_seq = col_index(cols, "_atom_site.label_seq_id")
    i_auth_seq = col_index(cols, "_atom_site.auth_seq_id") if "_atom_site.auth_seq_id" in cols else i_seq
    i_x = col_index(cols, "_atom_site.Cartn_x")
    i_y = col_index(cols, "_atom_site.Cartn_y")
    i_z = col_index(cols, "_atom_site.Cartn_z")
    i_group = col_index(cols, "_atom_site.group_PDB") if "_atom_site.group_PDB" in cols else None

    out: dict[int, np.ndarray] = {}
    for row in rows:
        if i_group is not None and row[i_group] not in ("ATOM", "HETATM"):
            continue
        if row[i_atom] not in ("CA", "CA"):
            continue
        xyz = np.array([float(row[i_x]), float(row[i_y]), float(row[i_z])], dtype=float)
        if is_prediction:
            chain = row[i_asym]
            if chain not in mapping:
                continue  # peptide / extra
            pred_start, pred_end, offset = mapping[chain]
            try:
                pred_res = int(row[i_auth_seq] if row[i_auth_seq] not in (".", "?") else row[i_seq])
            except ValueError:
                continue
            if pred_res < pred_start or pred_res > pred_end:
                continue
            uniprot = pred_res + offset
        else:
            chain = row[i_auth_asym]
            if chain != template_chain:
                continue
            try:
                uniprot = int(row[i_auth_seq] if row[i_auth_seq] not in (".", "?") else row[i_seq])
            except ValueError:
                continue
        # first CA wins (ignore altloc duplicates)
        out.setdefault(uniprot, xyz)
    return out


def all_atom_coords(cols: list[str], rows: list[list[str]]) -> np.ndarray:
    i_x = col_index(cols, "_atom_site.Cartn_x")
    i_y = col_index(cols, "_atom_site.Cartn_y")
    i_z = col_index(cols, "_atom_site.Cartn_z")
    return np.array([[float(r[i_x]), float(r[i_y]), float(r[i_z])] for r in rows], dtype=float)


def kabsch_transform(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Return R, t such that (X - P.mean) @ R.T + Q.mean superimposes P onto Q.

    P is mobile (predicted Cα), Q is target (template Cα). RMSD is after superposition.
    """
    if len(P) < 3:
        raise ValueError(f"Need >=3 paired Cα, got {len(P)}")
    Pc = P.mean(0)
    Qc = Q.mean(0)
    P0 = P - Pc
    Q0 = Q - Qc
    U, _S, Vt = np.linalg.svd(P0.T @ Q0)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    P_aligned = P0 @ R.T + Qc
    rmsd = float(np.sqrt(np.mean(np.sum((P_aligned - Q) ** 2, axis=1))))
    return R, Pc, Qc, rmsd


def apply_transform(X: np.ndarray, R: np.ndarray, Pc: np.ndarray, Qc: np.ndarray) -> np.ndarray:
    return (X - Pc) @ R.T + Qc


def write_transformed_cif(
    src: Path,
    dest: Path,
    header: list[str],
    cols: list[str],
    rows: list[list[str]],
    coords: np.ndarray,
) -> None:
    i_x = col_index(cols, "_atom_site.Cartn_x")
    i_y = col_index(cols, "_atom_site.Cartn_y")
    i_z = col_index(cols, "_atom_site.Cartn_z")
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Rebuild file: original header up to atom_site loop, then loop + transformed rows.
    # parse_atom_site stored pre-loop in header AND post-loop remainder appended after rows
    # were consumed. That's awkward — instead rewrite from original file, replacing xyz.
    original = src.read_text().splitlines(keepends=True)
    # Find atom_site data start
    loop_idx = None
    col_end = None
    for i, line in enumerate(original):
        if line.strip() == "loop_":
            j = i + 1
            if j < len(original) and original[j].lstrip().startswith("_atom_site."):
                loop_idx = i
                while j < len(original) and original[j].lstrip().startswith("_atom_site."):
                    j += 1
                col_end = j
                break
    if loop_idx is None:
        raise RuntimeError("atom_site loop not found on write")
    out_lines = original[:col_end]
    data_i = col_end
    row_i = 0
    while data_i < len(original) and row_i < len(rows):
        stripped = original[data_i].strip()
        if not stripped or stripped.startswith("#"):
            out_lines.append(original[data_i])
            data_i += 1
            continue
        if stripped.startswith("loop_") or (
            stripped.startswith("_") and not stripped.startswith("_atom_site.")
        ):
            break
        row = list(rows[row_i])
        x, y, z = coords[row_i]
        row[i_x] = f"{x:.5f}"
        row[i_y] = f"{y:.5f}"
        row[i_z] = f"{z:.5f}"
        out_lines.append(" ".join(row) + "\n")
        row_i += 1
        data_i += 1
    if row_i != len(rows):
        raise RuntimeError(f"Wrote {row_i} atom rows, expected {len(rows)}")
    out_lines.extend(original[data_i:])
    dest.write_text("".join(out_lines))


def superimpose_one(
    pred_path: Path,
    template_path: Path,
    mapping: dict[str, tuple[int, int, int]],
    out_path: Path,
    template_chain: str = "A",
) -> dict:
    _th, tcols, trows = parse_atom_site(template_path)
    ph, pcols, prows = parse_atom_site(pred_path)
    t_ca = receptor_ca_by_uniprot(
        tcols, trows, is_prediction=False, mapping=mapping, template_chain=template_chain
    )
    p_ca = receptor_ca_by_uniprot(
        pcols, prows, is_prediction=True, mapping=mapping, template_chain=template_chain
    )
    common = sorted(set(t_ca) & set(p_ca))
    P = np.stack([p_ca[k] for k in common])
    Q = np.stack([t_ca[k] for k in common])
    R, Pc, Qc, rmsd = kabsch_transform(P, Q)
    atoms = all_atom_coords(pcols, prows)
    aligned = apply_transform(atoms, R, Pc, Qc)
    write_transformed_cif(pred_path, out_path, ph, pcols, prows, aligned)

    # Verify the written file actually carries the transform
    _wh, wcols, wrows = parse_atom_site(out_path)
    w_ca = receptor_ca_by_uniprot(
        wcols, wrows, is_prediction=True, mapping=mapping, template_chain=template_chain
    )
    W = np.stack([w_ca[k] for k in common])
    written_rmsd = float(np.sqrt(np.mean(np.sum((W - Q) ** 2, axis=1))))
    centroid_dist = float(np.linalg.norm(W.mean(0) - Q.mean(0)))
    raw_centroid_dist = float(np.linalg.norm(P.mean(0) - Q.mean(0)))
    raw_rmsd = float(np.sqrt(np.mean(np.sum((P - Q) ** 2, axis=1))))
    per_res = np.sqrt(np.sum((W - Q) ** 2, axis=1))
    return {
        "n_paired_ca": len(common),
        "rmsd_after_A": rmsd,
        "rmsd_written_A": written_rmsd,
        "raw_rmsd_no_superpose_A": raw_rmsd,
        "centroid_dist_after_A": centroid_dist,
        "centroid_dist_before_A": raw_centroid_dist,
        "max_ca_dist_after_A": float(per_res.max()),
        "mean_ca_dist_after_A": float(per_res.mean()),
        "uniprot_min": common[0],
        "uniprot_max": common[-1],
    }


def run_dir(
    pred_dir: Path,
    templates: dict[str, Path],
    mapping: dict[str, tuple[int, int, int]],
    aligned_dir: Path,
) -> list[dict]:
    aligned_dir.mkdir(parents=True, exist_ok=True)
    for name, tpath in templates.items():
        shutil.copy2(tpath, aligned_dir / f"{name}_template.cif")
    models = sorted(pred_dir.glob("model_*.cif"), key=lambda p: int(p.stem.split("_")[1]))
    rows = []
    for model in models:
        idx = int(model.stem.split("_")[1])
        rec = {"model": idx, "file": model.name}
        for tname, tpath in templates.items():
            out = aligned_dir / f"model_{idx}_aligned_to_{tname}.cif"
            stats = superimpose_one(model, tpath, mapping, out)
            rec[f"n_paired_ca_{tname}"] = stats["n_paired_ca"]
            rec[f"rmsd_to_{tname}_A"] = round(stats["rmsd_after_A"], 3)
            rec[f"raw_rmsd_to_{tname}_A"] = round(stats["raw_rmsd_no_superpose_A"], 3)
            rec[f"centroid_before_{tname}_A"] = round(stats["centroid_dist_before_A"], 3)
            rec[f"centroid_after_{tname}_A"] = round(stats["centroid_dist_after_A"], 6)
            rec[f"max_ca_dist_to_{tname}_A"] = round(stats["max_ca_dist_after_A"], 3)
            rec[f"written_rmsd_to_{tname}_A"] = round(stats["rmsd_written_A"], 3)
            rec[f"uniprot_span_{tname}"] = f"{stats['uniprot_min']}-{stats['uniprot_max']}"
        rows.append(rec)
        print(
            f"model_{idx}: "
            + "  ".join(
                f"RMSD vs {t} = {rec[f'rmsd_to_{t}_A']:.3f} Å "
                f"(raw {rec[f'raw_rmsd_to_{t}_A']:.1f} Å, "
                f"centroid {rec[f'centroid_before_{t}_A']:.1f}→{rec[f'centroid_after_{t}_A']:.1e} Å)"
                for t in templates
            ),
            flush=True,
        )
    table = aligned_dir / "rmsd_table.csv"
    if rows:
        fieldnames = list(rows[0].keys())
        with table.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)
        md = aligned_dir / "rmsd_table.md"
        tnames = list(templates)
        lines = [
            "# Receptor Cα RMSD after Kabsch superposition",
            "",
            "RMSD is computed on paired receptor Cα (UniProt numbering) **after**",
            "the rigid transform that was applied to all atoms. `raw_rmsd` is the",
            "same pairs with **no** superposition (should be tens of Å). After",
            "superposition, template and model centroids coincide (~0 Å) so the",
            "structures overlap in PyMOL/ChimeraX without `align`.",
            "",
            "| model | " + " | ".join(f"RMSD vs {t} (Å)" for t in tnames) + " | paired Cα |",
            "|------:|" + "|".join(["----------------:"] * len(tnames)) + "|----------:|",
        ]
        for r in rows:
            t0 = tnames[0]
            lines.append(
                f"| {r['model']} | "
                + " | ".join(f"{r[f'rmsd_to_{t}_A']:.3f}" for t in tnames)
                + f" | {r[f'n_paired_ca_{t0}']} |"
            )
        md.write_text("\n".join(lines) + "\n")
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred-dir", type=Path, required=True)
    ap.add_argument("--aligned-dir", type=Path, required=True)
    ap.add_argument("--preset", choices=PRESETS, required=True)
    ap.add_argument("--templates-dir", type=Path, required=True)
    args = ap.parse_args()
    mapping = PRESETS[args.preset]
    templates = {
        "5UAK": args.templates_dir / "5UAK.cif",
        "6O2P": args.templates_dir / "6O2P.cif",
    }
    for p in templates.values():
        if not p.exists():
            sys.exit(f"missing template {p}")
    rows = run_dir(args.pred_dir, templates, mapping, args.aligned_dir)
    print(f"wrote {len(rows)} models → {args.aligned_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
