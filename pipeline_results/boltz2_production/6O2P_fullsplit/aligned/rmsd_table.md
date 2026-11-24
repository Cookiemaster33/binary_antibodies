# Receptor Cα RMSD after Kabsch superposition

RMSD is computed on paired receptor Cα (UniProt numbering) **after**
the rigid transform that was applied to all atoms. `raw_rmsd` is the
same pairs with **no** superposition (should be tens of Å). After
superposition, template and model centroids coincide (~0 Å) so the
structures overlap in PyMOL/ChimeraX without `align`.

| model | RMSD vs 5UAK (Å) | RMSD vs 6O2P (Å) | paired Cα |
|------:|----------------:|----------------:|----------:|
| 0 | 7.647 | 3.725 | 1114 |
| 1 | 7.878 | 3.783 | 1114 |
| 2 | 7.839 | 3.685 | 1114 |
| 3 | 7.481 | 3.178 | 1114 |
| 4 | 7.355 | 3.764 | 1114 |
| 5 | 7.623 | 3.627 | 1114 |
| 6 | 7.002 | 3.624 | 1114 |
| 7 | 7.717 | 3.520 | 1114 |
| 8 | 7.587 | 3.863 | 1114 |
| 9 | 7.255 | 3.550 | 1114 |
| 10 | 7.563 | 3.486 | 1114 |
| 11 | 7.615 | 3.654 | 1114 |
| 12 | 7.670 | 3.654 | 1114 |
| 13 | 7.562 | 3.849 | 1114 |
| 14 | 7.395 | 3.858 | 1114 |
| 15 | 7.818 | 3.672 | 1114 |
| 16 | 7.466 | 3.408 | 1114 |
| 17 | 7.743 | 3.599 | 1114 |
| 18 | 7.523 | 3.483 | 1114 |
| 19 | 7.778 | 3.773 | 1114 |
| 20 | 7.770 | 3.750 | 1114 |
| 21 | 7.614 | 3.476 | 1114 |
| 22 | 7.703 | 3.759 | 1114 |
| 23 | 7.828 | 3.397 | 1114 |
| 24 | 7.517 | 3.403 | 1114 |
