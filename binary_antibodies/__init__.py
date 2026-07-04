"""
binary_antibodies
=================
Computational design tools for conditional proximity-gated nanobody constructs.

A nanobody is tethered to an antibody via a flexible linker. The nanobody
can only bind its target when the antibody is anchored to a specific antigen,
creating an AND-gate conditional binding behaviour.
"""

from .polymer import LinkerModel
from .design import ConditionalConstruct, SplitScFvConstruct
from .sequences import LinkerSequence
from .split_scfv import SplitScFvSwitch
from .minibinder import MinibinderDesignSpec
from .kicker import KickerGeometry, KickerConstruct
from .structures import (
    download_pdb, load_structure, extract_variable_domain,
    analyse_vh_vl_interface, interface_summary,
    EXAMPLE_STRUCTURES, print_design_overview,
)

__all__ = [
    "LinkerModel", "ConditionalConstruct", "SplitScFvConstruct",
    "LinkerSequence", "SplitScFvSwitch", "MinibinderDesignSpec",
    "KickerGeometry", "KickerConstruct",
    "download_pdb", "load_structure", "extract_variable_domain",
    "analyse_vh_vl_interface", "interface_summary",
    "EXAMPLE_STRUCTURES", "print_design_overview",
]

__all__ = [
    "LinkerModel",
    "ConditionalConstruct",
    "SplitScFvConstruct",
    "LinkerSequence",
    "SplitScFvSwitch",
    "MinibinderDesignSpec",
    "download_pdb",
    "load_structure",
    "extract_variable_domain",
    "analyse_vh_vl_interface",
    "interface_summary",
    "EXAMPLE_STRUCTURES",
    "print_design_overview",
]
__version__ = "0.1.0"
