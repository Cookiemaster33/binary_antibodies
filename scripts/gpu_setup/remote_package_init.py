"""Minimal ``binary_antibodies`` package marker for GPU hosts.

The real package ``__init__`` pulls in the whole design/scoring stack (scipy,
pandas, ...). GPU post-processing only needs ``fab_hidden_switch`` and
``backbone_ss``, so this stub is uploaded in its place to keep the remote
dependency surface to numpy + biopython.
"""

__version__ = "0.1.0"
