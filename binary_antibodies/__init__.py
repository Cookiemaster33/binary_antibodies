"""
binary_antibodies
=================
Computational design tools for conditional proximity-gated nanobody constructs.

A nanobody is tethered to an antibody via a flexible linker. The nanobody
can only bind its target when the antibody is anchored to a specific antigen,
creating an AND-gate conditional binding behaviour.
"""

from .polymer import LinkerModel
from .design import ConditionalConstruct
from .sequences import LinkerSequence

__all__ = ["LinkerModel", "ConditionalConstruct", "LinkerSequence"]
__version__ = "0.1.0"
