"""Unit tests for PISA XML parsing (no Docker required)."""

from binary_antibodies.pisa_scoring import (
    parse_pisa_interfaces_xml,
    select_chain_pair_interface,
)

SAMPLE_XML = """<?xml version="1.0" encoding="UTF-8"?>
<pdb_entry>
  <interface>
    <id>1</id>
    <int_area>842.5</int_area>
    <int_solv_en>-8.3</int_solv_en>
    <stab_en>-4.1</stab_en>
    <pvalue>0.02</pvalue>
    <h-bonds><n_bonds>5</n_bonds></h-bonds>
    <salt-bridges><n_bonds>1</n_bonds></salt-bridges>
    <molecule>
      <chain_id>A</chain_id>
      <int_area>421.2</int_area>
      <int_solv_en>-8.3</int_solv_en>
      <int_nres>12</int_nres>
    </molecule>
    <molecule>
      <chain_id>B</chain_id>
      <int_area>421.3</int_area>
      <int_solv_en>-8.3</int_solv_en>
      <int_nres>11</int_nres>
    </molecule>
  </interface>
  <interface>
    <id>2</id>
    <int_area>120.0</int_area>
    <int_solv_en>-1.0</int_solv_en>
    <molecule><chain_id>A</chain_id></molecule>
    <molecule><chain_id>C</chain_id></molecule>
  </interface>
</pdb_entry>
"""


def test_parse_and_select_ab_interface():
    interfaces = parse_pisa_interfaces_xml(SAMPLE_XML)
    assert len(interfaces) == 2
    ab = select_chain_pair_interface(interfaces, "A", "B")
    assert ab is not None
    assert ab["int_area_A2"] == 842.5
    assert ab["int_solv_en_kcal"] == -8.3
    assert ab["n_h_bonds"] == 5
    assert ab["n_salt_bridges"] == 1
