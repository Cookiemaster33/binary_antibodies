"""Tests for fused Fab Boltz chain handling."""

from __future__ import annotations

import unittest

from binary_antibodies.boltz_fab_chains import (
    BOLTZ_CHAIN_MODE_FUSED,
    BOLTZ_CHAIN_MODE_SPLIT,
    boltz_holo_chain_sequences,
    default_boltz_chain_mode,
    expand_fused_fab_structure,
    format_boltz_yaml,
    fuse_fab_sequences,
    is_fused_boltz_structure,
    load_boltz_chain_mode,
)


class TestBoltzFabChains(unittest.TestCase):
    def test_default_mode_by_scope(self) -> None:
        self.assertEqual(default_boltz_chain_mode("full_fab"), BOLTZ_CHAIN_MODE_FUSED)
        self.assertEqual(default_boltz_chain_mode("fv"), BOLTZ_CHAIN_MODE_SPLIT)

    def test_load_mode_from_config(self) -> None:
        cfg = {"validation": {"boltz_chain_mode": "split"}}
        self.assertEqual(load_boltz_chain_mode(cfg, "full_fab"), BOLTZ_CHAIN_MODE_SPLIT)

    def test_fuse_sequences(self) -> None:
        heavy, light = fuse_fab_sequences("VH", "VL", "CH1", "CL")
        self.assertEqual(heavy, "VHCH1")
        self.assertEqual(light, "VLCL")

    def test_boltz_yaml_fused(self) -> None:
        chains = {"A": "VH", "B": "VL", "C": "CH1", "D": "CL"}
        entries = boltz_holo_chain_sequences(chains, mode=BOLTZ_CHAIN_MODE_FUSED, epitope="EP")
        self.assertEqual(entries[0], ("H", "VHCH1"))
        self.assertEqual(entries[1], ("L", "VLCL"))
        self.assertEqual(entries[2], ("T", "EP"))
        yaml = format_boltz_yaml(entries)
        self.assertIn('id: H', yaml)
        self.assertIn('sequence: "VHCH1"', yaml)
        self.assertNotIn('id: C', yaml)

    def test_boltz_yaml_split(self) -> None:
        chains = {"A": "VH", "B": "VL", "C": "CH1", "D": "CL"}
        entries = boltz_holo_chain_sequences(chains, mode=BOLTZ_CHAIN_MODE_SPLIT, epitope="EP")
        self.assertEqual([c for c, _ in entries], ["A", "B", "C", "D", "T"])

    def test_expand_fused_structure(self) -> None:
        try:
            import numpy as np
            from biotite.structure import AtomArray
        except ImportError:
            self.skipTest("biotite not installed")

        vh_len, vl_len = 3, 2
        n_heavy, n_light, n_t = 5, 4, 2
        n = n_heavy + n_light + n_t
        aa = AtomArray(n)
        aa.atom_name = np.array(["CA"] * n)
        aa.res_name = np.array(["ALA"] * n)
        aa.element = np.array(["C"] * n)
        aa.coord = np.zeros((n, 3))
        aa.chain_id = np.array(["H"] * n_heavy + ["L"] * n_light + ["T"] * n_t)
        aa.res_id = np.array([1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2])

        self.assertTrue(is_fused_boltz_structure(aa, vh_len, vl_len))
        expanded = expand_fused_fab_structure(aa, vh_len, vl_len)
        chains = set(expanded.chain_id.tolist())
        self.assertEqual(chains, {"A", "B", "C", "D", "T"})
        self.assertEqual(
            sorted(int(r) for r in expanded.res_id[expanded.chain_id == "C"]),
            [1, 2],
        )
        self.assertEqual(
            sorted(int(r) for r in expanded.res_id[expanded.chain_id == "A"]),
            [1, 2, 3],
        )


if __name__ == "__main__":
    unittest.main()
