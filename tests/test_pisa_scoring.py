"""Unit tests for PISA XML parsing (no Docker required)."""

from binary_antibodies.pisa_scoring import (
    buried_interface_residues,
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
      <residue>
        <seq_num>42</seq_num>
        <name>LEU</name>
        <bsa>18.5</bsa>
        <asa>12.0</asa>
      </residue>
      <residue>
        <seq_num>43</seq_num>
        <name>VAL</name>
        <bsa>0.0</bsa>
        <asa>40.0</asa>
      </residue>
      <residue>
        <seq_num>44</seq_num>
        <name>ALA</name>
        <bsa>6.2</bsa>
        <asa>8.0</asa>
      </residue>
    </molecule>
    <molecule>
      <chain_id>B</chain_id>
      <int_area>421.3</int_area>
      <int_solv_en>-8.3</int_solv_en>
      <int_nres>11</int_nres>
      <residue>
        <seq_num>45</seq_num>
        <name>PHE</name>
        <bsa>22.1</bsa>
        <asa>5.0</asa>
      </residue>
      <residue>
        <seq_num>46</seq_num>
        <name>TYR</name>
        <bsa>3.0</bsa>
        <asa>15.0</asa>
      </residue>
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


def test_buried_interface_residues_filters_bsa_and_cdr_truncation():
    interfaces = parse_pisa_interfaces_xml(SAMPLE_XML)
    ab = select_chain_pair_interface(interfaces, "A", "B")
    assert ab is not None

    buried = buried_interface_residues(ab, "A", "B", min_buried_sasa_A2=0.0)
    assert buried["A"] == [42, 44]
    assert buried["B"] == [45, 46]

    buried_min = buried_interface_residues(ab, "A", "B", min_buried_sasa_A2=7.0)
    assert buried_min["A"] == [42]
    assert buried_min["B"] == [45]

    truncated = buried_interface_residues(
        ab, "A", "B", min_buried_sasa_A2=0.0, max_resnum={"A": 43, "B": 107}
    )
    assert truncated["A"] == [42]
    assert truncated["B"] == [45, 46]
