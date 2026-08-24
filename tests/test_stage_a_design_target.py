"""Tests for Stage A design target (full Fab, unlinked minibinder)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

from binary_antibodies.fab_hidden_switch import (  # noqa: E402
    DEFAULT_STAGE_A_SEPARATION_A,
    VH_END,
    VL_END,
    build_fab_context_pdb,
    build_holo_split_cif,
    build_stage_a_design_target_pdb,
    extract_chain_sequences,
    fab_has_split_chains,
    infer_already_split,
    stage_a_contig,
    stage_a_fixed_atoms,
    stage_a_rfd3_config,
    stage_a_unindex,
)


class TestStageADesignTarget(unittest.TestCase):
    def test_contig_unlinked_minibinder(self):
        contig = stage_a_contig(113, 107, 107, 107, "35-55")
        self.assertIn("/0,35-55,", contig)
        self.assertTrue(contig.startswith("B1-107"))
        self.assertTrue(contig.endswith("C1-107"))

    def test_unindex_fixed_context(self):
        unindex = stage_a_unindex(113, 107, 107, 12)
        self.assertIn("A1-113", unindex)
        self.assertIn("D1-107", unindex)
        self.assertIn("T1-12", unindex)

    def test_rfd3_config_uses_unindex(self):
        cfg = stage_a_rfd3_config(113, 107, 107, 107, "35-55", 12)
        self.assertIn("unindex", cfg)
        self.assertIn("A1-113", cfg["unindex"])
        self.assertNotIn("A1-113", cfg["contig"])

    def test_fixed_atoms_cover_full_fab(self):
        fixed = stage_a_fixed_atoms(113, 107, 107, 107, 12)
        self.assertEqual(set(fixed), {"A1-113", "B1-107", "C1-107", "D1-107", "T1-12"})

    def test_holo_split_cif_chains(self):
        fab_cif = (
            ROOT
            / "pipeline_results/stage_0_full_fab_fused_t025/structures/holo"
            / "rank079_s0_native_split_s296_model_0.cif"
        )
        if not fab_cif.exists():
            self.skipTest("rank079 holo CIF not in workspace")

        with tempfile.TemporaryDirectory() as tmp:
            split_cif = Path(tmp) / "rank079_split.cif"
            build_holo_split_cif(split_cif, fab_cif, separation_a=DEFAULT_STAGE_A_SEPARATION_A)
            self.assertTrue(fab_has_split_chains(split_cif))
            self.assertTrue(infer_already_split(split_cif))
            seqs = extract_chain_sequences(split_cif)
            self.assertEqual(set(seqs), {"A", "B", "C", "D", "T"})

    def test_split_cif_skips_retranslation(self):
        fab_cif = (
            ROOT
            / "pipeline_results/stage_0_full_fab_fused_t025/structures/holo"
            / "rank079_s0_native_split_s296_model_0.cif"
        )
        if not fab_cif.exists():
            self.skipTest("rank079 holo CIF not in workspace")

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            split_cif = tmp_path / "rank079_split.cif"
            stage_a_from_split = tmp_path / "from_split.pdb"
            stage_a_from_holo = tmp_path / "from_holo.pdb"
            build_holo_split_cif(split_cif, fab_cif, separation_a=DEFAULT_STAGE_A_SEPARATION_A)
            build_stage_a_design_target_pdb(stage_a_from_split, source_pdb=split_cif, already_split=True)
            build_stage_a_design_target_pdb(
                stage_a_from_holo,
                source_pdb=fab_cif,
                separation_a=DEFAULT_STAGE_A_SEPARATION_A,
            )

            from Bio.PDB import PDBParser

            parser = PDBParser(QUIET=True)
            s = parser.get_structure("s", str(stage_a_from_split))
            h = parser.get_structure("h", str(stage_a_from_holo))
            model_s = list(s.get_models())[0]
            model_h = list(h.get_models())[0]

            def ca_centroid(chain_id: str, struct):
                atoms = [a for a in struct[chain_id].get_atoms() if a.name == "CA"]
                import numpy as np

                return np.mean([a.coord for a in atoms], axis=0)

            dist_split = sum((ca_centroid("C", model_s) - ca_centroid("B", model_s)) ** 2) ** 0.5
            dist_holo = sum((ca_centroid("C", model_h) - ca_centroid("B", model_h)) ** 2) ** 0.5
            self.assertAlmostEqual(dist_split, dist_holo, delta=0.5)

    def test_stage_a_increases_vl_ch1_distance(self):
        fab_cif = (
            ROOT
            / "pipeline_results/stage_0_full_fab_fused_t025/structures/holo"
            / "rank079_s0_native_split_s296_model_0.cif"
        )
        if not fab_cif.exists():
            self.skipTest("rank079 holo CIF not in workspace")

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            native_pdb = tmp_path / "native.pdb"
            stage_a_pdb = tmp_path / "stage_a.pdb"
            build_fab_context_pdb(native_pdb, source_pdb=fab_cif)
            build_stage_a_design_target_pdb(
                stage_a_pdb,
                source_pdb=fab_cif,
                separation_a=DEFAULT_STAGE_A_SEPARATION_A,
            )

            from Bio.PDB import PDBParser

            parser = PDBParser(QUIET=True)
            native = parser.get_structure("n", str(native_pdb))
            opened = parser.get_structure("o", str(stage_a_pdb))
            model_n = list(native.get_models())[0]
            model_o = list(opened.get_models())[0]

            def ca_centroid(chain_id: str, struct):
                atoms = [a for a in struct[chain_id].get_atoms() if a.name == "CA"]
                import numpy as np

                return np.mean([a.coord for a in atoms], axis=0)

            dist_native = sum(
                (ca_centroid("C", model_n) - ca_centroid("B", model_n)) ** 2
            ) ** 0.5
            dist_opened = sum(
                (ca_centroid("C", model_o) - ca_centroid("B", model_o)) ** 2
            ) ** 0.5
            # Correct H/L assignment: VL–CH1 are already ~35 Å apart in the holo;
            # separation opens the VH–VL interface and modestly adjusts VL/CL position.
            self.assertGreater(dist_opened, 25.0)
            self.assertNotAlmostEqual(dist_opened, dist_native, delta=0.5)

    def test_build_script_config(self):
        fab_cif = (
            ROOT
            / "pipeline_results/stage_0_full_fab_fused_t025/structures/holo"
            / "rank079_s0_native_split_s296_model_0.cif"
        )
        if not fab_cif.exists():
            self.skipTest("rank079 holo CIF not in workspace")

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            out_json = tmp_path / "cfg.json"
            out_pdb = tmp_path / "target.pdb"
            import subprocess
            import sys

            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts/build_stage_a_design_target.py"),
                    "--fab-pdb",
                    str(fab_cif),
                    "--out-pdb",
                    str(out_pdb),
                    "--out-json",
                    str(out_json),
                ],
                check=True,
            )
            cfg = json.loads(out_json.read_text())
            self.assertIn("/0,35-55,", cfg["rfd3"]["contig"])
            self.assertIn("unindex", cfg["rfd3"])
            self.assertIn("A1-113", cfg["rfd3"]["select_fixed_atoms"])
            self.assertFalse(cfg["layout"]["minibinder_linked"])
            seqs = extract_chain_sequences(out_pdb)
            self.assertEqual(set(seqs), {"A", "B", "C", "D", "T"})
            self.assertEqual(len(seqs["A"]), VH_END)
            self.assertEqual(len(seqs["B"]), VL_END)


if __name__ == "__main__":
    unittest.main()
