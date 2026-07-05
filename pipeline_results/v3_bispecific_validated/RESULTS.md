# Boltz-2 Validation Results — Bispecific Bridging Minibinder

## Top-10 designs ranked by sqrt(ipTM_VH1 x ipTM_Nb)

| Rank | Backbone | sqrt_ipTM | ipTM↔VH1 | ipTM↔Nb | pLDDT% | scRMSD |
|------|----------|-----------|----------|---------|--------|--------|
| 1 | mb_b005_000 | 0.850 | 0.878 | 0.822 | 77% | 14.8 Å |
| 2 | mb_b018_005 | 0.841 | 0.833 | 0.850 | 78% | 13.1 Å |
| 3 | mb_b000_001 | 0.820 | 0.819 | 0.822 | 81% | 13.0 Å |
| 4 | mb_b001_005 | 0.806 | 0.826 | 0.786 | 83% | 13.7 Å |
| 5 | mb_b003_002 | 0.802 | 0.817 | 0.787 | 85% | 9.4 Å |
| 6 | mb_b019_000 | 0.784 | 0.820 | 0.750 | 78% | 15.5 Å |
| 7 | mb_b017_004 | 0.706 | 0.754 | 0.661 | 85% | 12.1 Å |
| 8 | mb_b010_001 | 0.695 | 0.694 | 0.697 | 84% | 14.1 Å |
| 9 | mb_b008_009 | 0.693 | 0.730 | 0.658 | 78% | 13.1 Å |
| 10 | mb_b015_004 | 0.677 | 0.646 | 0.709 | 82% | 12.8 Å |

## Structure files (in structures/)

Two CIF files per design:
- `top??_*_rfd3.cif` — RFdiffusion3 designed backbone (chains A=VH1, B=Minibinder, C=Nanobody)
- `top??_*_boltz2.cif` — Boltz-2 predicted complex

Open in PyMOL or ChimeraX. Superimpose on chains A+C to compare minibinder position.

## Note on scRMSD (9-16 Å)

This is expected for bispecific binders — not a failure. See full explanation in PR description.
The ipTM scores (0.65-0.88) are the correct quality metric for interface confidence.
