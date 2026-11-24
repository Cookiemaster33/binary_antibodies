# Receptor Cα RMSD after Kabsch superposition

RMSD is computed on paired receptor Cα (UniProt numbering) **after**
the rigid transform that was applied to all atoms. `raw_rmsd` is the
same pairs with **no** superposition (should be tens of Å). After
superposition, template and model centroids coincide (~0 Å) so the
structures overlap in PyMOL/ChimeraX without `align`.

| model | RMSD vs 5UAK (Å) | RMSD vs 6O2P (Å) | paired Cα |
|------:|----------------:|----------------:|----------:|
| 0 | 2.991 | 5.253 | 670 |
| 1 | 2.688 | 5.732 | 670 |
| 2 | 4.531 | 5.908 | 670 |
| 3 | 2.900 | 5.160 | 670 |
| 4 | 3.103 | 5.306 | 670 |
| 5 | 5.103 | 6.118 | 670 |
| 6 | 3.743 | 4.999 | 670 |
| 7 | 2.825 | 5.529 | 670 |
| 8 | 3.835 | 5.924 | 670 |
| 9 | 3.390 | 5.025 | 670 |
| 10 | 3.455 | 4.727 | 670 |
| 11 | 2.690 | 5.314 | 670 |
| 12 | 3.209 | 5.254 | 670 |
| 13 | 3.331 | 5.994 | 670 |
| 14 | 2.999 | 5.310 | 670 |
| 15 | 3.417 | 4.774 | 670 |
| 16 | 3.797 | 4.796 | 670 |
| 17 | 3.508 | 5.496 | 670 |
| 18 | 3.531 | 4.791 | 670 |
| 19 | 3.325 | 5.389 | 670 |
| 20 | 3.395 | 4.469 | 670 |
| 21 | 3.116 | 5.345 | 670 |
| 22 | 3.257 | 4.753 | 670 |
| 23 | 3.811 | 5.471 | 670 |
| 24 | 3.060 | 4.973 | 670 |
