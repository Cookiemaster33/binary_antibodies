# Stage 0 — holo + PISA run (2026-07-27)

**Instance:** `e0f037e22d4146c4b687ad9287005f74` (A100, us-east-1)  
**Config:** Whole interface (23 residues), 1000 MPNN → top 100 holo Boltz + PISA  
**Result:** **0/100** pass all filters

## Headline

**rank010 (s913)** is the clear Stage A candidate: passes every holo filter, PISA shows **+7.1 kcal/mol** weaker VH–VL vs native, excellent geometry (fw RMSD 1.89 Å, CDR RMSD 0.66 Å, 377 CDR–T contacts). Only fails the **static contact filter** at 47.6% WT (threshold 45%).

## Native reference (PISA)

| Metric | Value |
|--------|-------|
| Interface area | 663 Å² |
| Solvation energy | −10.5 kcal/mol |
| H-bonds | 7 |

## Filter breakdown

| Filter | Failures |
|--------|----------|
| Static ≤45% WT | **100/100** |
| Holo interface clashes | 20/100 |
| Framework RMSD ≤3.5 Å | 0/100 |
| CDR RMSD ≤6 Å | 0/100 |
| CDR–T ≥6 | 0/100 |

MPNN pool sent to Boltz: static contacts **51.1–58.9** (0/100 below 45% WT; 12/100 below 48%).

## Top candidates (relaxed static ≤48%)

| Design | Static (%WT) | PISA Δ | fw RMSD | CDR–T |
|--------|--------------|--------|---------|-------|
| **rank010 s913** | 53.8 (47.6%) | **+7.12** | 1.89 | 377 |
| rank006 s175 | 53.8 (47.6%) | +5.93 | 2.48 | 102 |
| rank008 s712 | 53.8 (47.6%) | +5.87 | 2.16 | 262 |
| rank007 s264 | 53.8 (47.6%) | +4.44 | 2.56 | 211 |
| rank009 s836 | 53.8 (47.6%) | +3.20 | 2.59 | 208 |
| rank004 s349 | 52.9 (46.8%) | +1.14 | 2.46 | 230 |

## PISA vs static

Correlation r = −0.12 (weak): lowest-static designs (s805/s209 at 45.2%) do **not** have the best PISA weakening. **rank010** and **rank005** are energetically much weaker per PISA despite similar static scores.

## Recommendation

1. **Advance rank010 (s913)** to Stage A with relaxed static filter (or replace static gate with PISA Δ ≥ +5 kcal/mol).
2. Relax `max_static_fraction_of_native_contacts` to **0.48** — would pass 6/100 including rank010.
3. Consider **partial interface** (14 rim residues) MPNN to push static below 45% while keeping holo fidelity.

## Pulled artifacts

| File | Description |
|------|-------------|
| `stage_0_results.json` | Full scoring for top 100 holo designs |
| `top_designs_stage_0.json` | MPNN-ranked top 100 sent to Boltz |
| `mpnn_all_sequences.json` | All 1000 MPNN sequences |
| `top5_candidates.json` | Metadata for top 5 candidates below |
| `structures/top5_holo/*.cif` | Holo Boltz model_0 CIFs (chains A+B+C+D+T) |

### Top 5 holo CIFs

| Rank | Design | CIF | PISA Δ |
|------|--------|-----|--------|
| 1 | rank010 s913 | `structures/top5_holo/rank010_s0_native_split_s913_model_0.cif` | +7.12 |
| 2 | rank006 s175 | `structures/top5_holo/rank006_s0_native_split_s175_model_0.cif` | +5.93 |
| 3 | rank008 s712 | `structures/top5_holo/rank008_s0_native_split_s712_model_0.cif` | +5.87 |
| 4 | rank007 s264 | `structures/top5_holo/rank007_s0_native_split_s264_model_0.cif` | +4.44 |
| 5 | rank009 s836 | `structures/top5_holo/rank009_s0_native_split_s836_model_0.cif` | +3.20 |
