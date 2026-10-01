# Stage B Candidate

## Design: sa_60to85_b004_001

Selected as rank #2 from Stage A scoring across 1,600 designs (8 runs).
Chosen for near-perfect helicity and the most balanced contact distribution of
any high-scoring design.

## Metrics

| Property           | Value        |
|--------------------|--------------|
| Source run         | stage_a_auto_split_30a_v1 (v1 hotspots, 30 Å gap) |
| Minibinder length  | 85 aa        |
| Helix fraction     | 98.8%        |
| CH1 contacts (C)   | 14           |
| VL contacts (B)    | 13           |
| Geometric mean     | 13.5         |
| Min Fab clearance  | 2.03 Å       |
| Radius of gyration | 12.12 Å      |
| Max span           | 36.81 Å      |
| Selection score    | 27.06        |

## Scoring formula

    score = sqrt(ch1_contacts × vl_contacts) × helix_frac × min_fab_clash_a

Hard filters applied before ranking:
- bridges_both = 1
- min_fab_clash_a ≥ 1.5 Å
- helix_frac ≥ 0.85
- ch1_contacts ≥ 3 AND vl_contacts ≥ 3

## Input file

`stage_b_input_sa_60to85_b004_001.pdb` — grafted structure with:
- Chain A: VH
- Chain B: VL
- Chain C: CH1
- Chain D: CL
- Chain M: designed minibinder (85 aa, fully helical)

## Next step (Stage B)

Design the target-binding arm (nanobody or VHH) against the HER2 epitope stub,
using this backbone as the fixed scaffold for the crosslinker (Chain M).
