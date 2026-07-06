# Pipeline Results — Progress Across Runs

## scRMSD improvement trajectory

| Run | Designs | MB Length | Best scRMSD | Key change |
|-----|---------|-----------|-------------|------------|
| v2 | 200 | fixed 55 | 7.26 Å | Baseline (single-chain Boltz-2) |
| v3 | 500 | fixed 55 | 6.03 Å | 2.5× more designs |
| **v4** | **200** | **variable 35-70** | **5.18 Å** | **Variable length + select_fixed_atoms** |

## v4 — variable-length minibinder results

**Minibinder length distribution** (RFd3 sampled from 35-70):
- Short (35-40 res): 40 designs — best scRMSD candidates
- Medium (51-55 res): 100 designs — most common
- Long (58-70 res): 60 designs

**Top 5 designs (sorted by scRMSD):**

| Rank | Backbone | scRMSD | pLDDT% | pTM |
|------|----------|--------|--------|-----|
| 1 | mb_b003_000 | **5.18 Å** | 77.5% | 0.471 |
| 2 | mb_b006_005 | 5.87 Å | 80.3% | 0.512 |
| 3 | mb_b001_005 | 6.27 Å | 80.8% | 0.505 |
| 4 | mb_b004_002 | 6.30 Å | 76.4% | 0.532 |
| 5 | mb_b001_000 | 6.41 Å | 78.5% | 0.506 |

## Next steps to push below 2 Å

1. **Narrow the length range** to 35-45 residues (shorter = fewer DOF = lower scRMSD)
2. **Add Boltz-2 contact constraints** (intra-chain MB↔VH1 and MB↔Nb contacts)
3. **More designs** at the shorter length

## Structure files

`structures/top01_mb_b003_000_rfd3.cif`   — RFd3 backbone
`structures/top01_mb_b003_000_boltz2.cif` — Boltz-2 predicted fold
