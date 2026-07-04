# Binary Antibodies — Conditional Nanobody Design

## Concept

A nanobody that can **only** bind its target when a VH1 domain is bound to a
specific antigen on the same cell surface — AND-gate logic.

---

## Current Design: CH1-Kicker Mechanism

```
  ── INACTIVE ──────────────────────────────────────────────────────────

                [Nanobody CDRs blocked by Minibinder]
                         |
                    ~flexible linker~
                         |
               [CL]──[VL1]    ← CL/CH1 positioned away from nanobody
                    (free, unanchored)
         [VH1]
          (not yet bound to antigen)

         [Antigen] ··· membrane ···  [Target]


  ── ACTIVE ────────────────────────────────────────────────────────────

                [Nanobody]────────────────────────────╮
                     |                           binds [Target]
                ~flexible linker~
                     |
               [CL]  ← CH1/CL physically KICKS Minibinder off Nanobody CDRs
               [VL1]──[VH1]  ← VH1-VL1 pairing brings CL into kicking position
                    |
               [Antigen] ··· membrane ···  [Target]
```

### Components

| Component | Role | How to obtain |
|---|---|---|
| **VH1** | Binds antigen autonomously; anchors construct to membrane | Engineered from existing anti-antigen VH |
| **VL1** | Pairs with VH1 to form scFv; triggers the kicker mechanism | From same antibody as VH1 |
| **CL (CH1)** | "Kicker" domain — physically displaces minibinder from nanobody CDRs upon VH1-VL1 pairing | Natural constant light chain, attached to VL1 C-terminus |
| **Minibinder** | Blocks nanobody CDR face in free state | **Designed against nanobody paratope** (anti-idiotypic design) |
| **Flexible linker** | Connects VL1-CL block to nanobody; must be long enough for kicking geometry | (G₄S)ₙ, optimised |
| **Nanobody** | Anti-target single-domain antibody | Selected from library or existing |

### Mechanism (structural switch — NOT thermodynamic competition)

1. **Free state**: Minibinder occupies the nanobody's CDR face. Nanobody cannot bind target.
   VH1 and VL1 are not paired (or floating freely in solution).

2. **Anchoring**: VH1 binds the antigen on the cell surface.

3. **Pairing**: VH1 recruits VL1; VH1-VL1 pair into a functional scFv.
   The CL domain comes along with VL1 (it is covalently attached at VL1's C-terminus).

4. **Kicking**: The CL domain, now in its new position relative to the nanobody
   (set by the VH1-VL1 pairing geometry), **sterically clashes with the minibinder**
   and physically displaces it from the nanobody's CDR loops.

5. **Activation**: Nanobody CDRs are now free. The flexible linker allows the
   nanobody to extend and bind the target.

### Why this design is better than competitive displacement

| Property | Competitive displacement (old) | CH1-kicker (new) |
|---|---|---|
| Mechanism | Thermodynamic (VH1 surface conc. vs intramolecular minibinder) | Structural/mechanical (steric kick) |
| OFF-state fidelity | Good — intramolecular minibinder is always present | Good — minibinder blocks CDRs |
| ON-state trigger | Needs µM VH1 surface concentration | Needs VH1-VL1 pairing (binary event) |
| Precision needed | High — tight Kd windows, exact linker lengths | Low — CH1 kick is a physical event |
| Minibinders needed | 2 | **1** |
| Robustness | Sensitive to antigen expression level | **Robust** — switch is structural |

### Key design parameter: kicking geometry

The CL domain must be positioned such that, **after VH1-VL1 pairing**, it overlaps
with the minibinder's footprint on the nanobody CDRs. This is determined by:

- The Fab geometry (VH1-VL1-CL complex structure)
- The flexible linker length between CL and the nanobody
- The size and position of the minibinder on the nanobody CDR face

Use the structural analysis in `binary_antibodies/kicker.py` to assess geometry.

### Chain B assembly (N→C)

```
SP – [VL1] – [CL] – (G4S)n – [Minibinder] – (G4S)3 – [Nanobody VHH]
```

- `SP`: signal peptide
- `VL1`: variable light domain (pairs with VH1 to trigger kicking)
- `CL`: constant light domain (the "kicker" — repositions upon VL1-VH1 pairing)
- `(G4S)n`: flexible linker (length sets kicking geometry)
- `Minibinder`: blocks nanobody CDRs in free state
- `(G4S)3`: short linker
- `Nanobody VHH`: anti-target single-domain antibody

Chain A (separate): **VH1** (binds antigen; expressed as separate polypeptide)

---

## Example Biological System

| Component | Molecule | PDB |
|---|---|---|
| Antigen | EGFR | Cetuximab Fab 1YY9 (VH template) |
| Target | HER2 | — |
| VH1 + VL1 | Anti-HER2 scFv | Trastuzumab Fab 1N8Z |
| CL | Trastuzumab Cκ | 1N8Z chain B (C-terminal constant domain) |
| Nanobody | 2Rs15d anti-HER2 VHH | 5MY6 |
| Minibinder target | **Nanobody CDR face** (CDR1, CDR2, CDR3) | `structures/domains/her2_nanobody_VHH.pdb` |

---

## Repository Layout

```
binary_antibodies/
├── README.md
├── requirements.txt
├── binary_antibodies/
│   ├── polymer.py          # linker physics
│   ├── design.py           # ConditionalConstruct, SplitScFvConstruct
│   ├── kicker.py           # CH1-kicker geometry and displacement model ← NEW
│   ├── split_scfv.py       # competitive displacement model (reference)
│   ├── minibinder.py       # minibinder design guidance
│   ├── sequences.py        # linker sequences
│   └── structures.py       # PDB download/analysis
├── scripts/
│   ├── setup_example.py
│   ├── interface_analysis.py        # VH-VL interface (reference)
│   ├── nanobody_cdr_analysis.py     # Nanobody CDR face analysis ← NEW
│   └── gpu_setup/
│       ├── run_minibinder_design_rfd3.sh   # rfd3 pipeline (against nanobody CDRs)
│       └── ...
└── structures/
    ├── domains/
    │   ├── her2_nanobody_VHH.pdb    ← TARGET for minibinder design
    │   ├── trastuzumab_VH.pdb
    │   └── trastuzumab_VL.pdb
    └── interface/
```
