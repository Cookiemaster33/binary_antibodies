# Binary Antibodies — Conditional Nanobody Design

## Concept

A nanobody that can **only** bind its target when a separate antibody arm is bound to a specific antigen. The system acts as a molecular AND-gate on a cell surface.

---

## Refined Design: Split-scFv Conditional Switch

```
  ── INACTIVE (VH1 bound to antigen, VL1 locked by minibinders) ───────

         [Nanobody]
              |
       ~flexible linker~
              |
       [Minibinders]──[VL1]      ← VL1 locked; cannot pair with VH1
                                    (minibinders block VH1-pairing face)
         [VH1]
              |
          [Antigen] ··· membrane ···  [Target]


  ── ACTIVE (VH1 displaces minibinders, VH1+VL1 scFv formed) ──────────

         [Nanobody]────────────────────────────────╮
              |                               binds [Target]
       ~flexible linker~
              |
       [Minibinders]  (displaced)
       [VH1]──[VL1]      ← functional scFv formed at membrane surface
              |
          [Antigen] ··· membrane ···  [Target]
```

### Components

| Component | Description | How to obtain |
|---|---|---|
| **VH1** | VH domain of a known antibody, engineered to bind its antigen autonomously | Site-directed mutagenesis of a known VH/VL pair; select for antigen binding without VL |
| **VL1** | VL domain of the same antibody; together with VH1 forms a functional scFv that binds the **Target** | Taken directly from the known antibody |
| **Minibinders** | Designed small proteins that bind VL1 on its VH1-pairing interface | **Computational de novo design** (RFdiffusion + ProteinMPNN) against the VL1 framework |
| **Flexible linker** | (G₄S)ₙ peptide connecting minibinder–VL1 complex to the nanobody | Optimised by `binary_antibodies` using polymer physics |
| **Nanobody** | Single-domain antibody against the Target | Selected from llama/camel immune library or designed |

### Mechanism (competitive displacement switch)

1. **Free state** (no antigen): Minibinders occupy VL1's VH1-pairing interface (framework 2 and CDR-L2 region). VH1 and VL1 have been *engineered to have weak affinity in solution* (Kd,VH-VL ≈ 1–100 µM). The scFv is therefore non-functional. The nanobody, tethered to the locked VL1, cannot reach the target.

2. **Anchoring** (VH1 binds antigen): VH1 docks to the antigen on the cell membrane. This brings VH1 into proximity with VL1 (which is tethered via the minibinder–flexible-linker chain). The effective local concentration of VH1 near VL1 is now **µM range**.

3. **Displacement** (VH1 outcompetes minibinders for VL1): Because Kd(VH1-VL1) < Kd(minibinder-VL1) at the effective local concentration, VH1 displaces the minibinders from VL1's framework interface. VH1 and VL1 pair into a functional scFv.

4. **Activation** (nanobody reaches target): The minibinder displacement removes the steric constraint on the flexible linker. The nanobody, now free to diffuse in the hemisphere above the membrane, engages the Target.

### Thermodynamic AND-gate condition

For reliable activation, the following inequalities must hold:

```
Kd(VH1–VL1) < C_eff(VH1 | anchored)   →  VH1-VL1 pairing is driven by proximity
Kd(VH1–VL1) < Kd(minibinder–VL1)      →  VH1 outcompetes minibinders when anchored
Kd(minibinder–VL1) << 1/[construct]   →  VL1 is reliably locked in free state
```

The intermediate quantity C_eff is determined by linker length and antigen–target distance
(computed by `binary_antibodies.polymer.LinkerModel`).

---

## Repository Layout

```
binary_antibodies/
├── README.md
├── requirements.txt
├── binary_antibodies/
│   ├── __init__.py
│   ├── polymer.py          # FJC/WLC linker physics, effective concentration
│   ├── design.py           # ConditionalConstruct: end-to-end construct design
│   ├── split_scfv.py       # Split-scFv competitive displacement thermodynamics
│   ├── minibinder.py       # Minibinder design guidance and target interface analysis
│   └── sequences.py        # Linker sequence generation & composition
├── scripts/
│   ├── optimise_linker.py  # CLI: find optimal linker for a geometry
│   └── scan_geometry.py    # CLI: heatmap over d × N parameter space
└── notebooks/
    └── design_walkthrough.ipynb
```

---

## Quick Start

```bash
pip install -r requirements.txt

# Analyse the thermodynamic feasibility of a split-scFv switch
python -c "
from binary_antibodies.split_scfv import SplitScFvSwitch
s = SplitScFvSwitch(
    kd_vh_vl_M=50e-6,       # engineered weak VH-VL affinity: 50 µM
    kd_minibinder_vl_M=5e-6, # minibinder affinity for VL1: 5 µM
    kd_nanobody_M=10e-9,     # nanobody for target: 10 nM
)
print(s.summary(distance_nm=8.0, linker_n_residues=60))
"

# Optimise linker length
python scripts/optimise_linker.py --distance 8.0 --kd-nanobody 10e-9 --fold 100
```

---

## Design Workflow

### Step 1 — Choose the antigen/target pair
Select a membrane antigen (highly expressed on the target cell type) and a target
membrane protein whose engagement you wish to conditionalise.

### Step 2 — Select a known antibody for the Target
Take an existing high-affinity antibody (or scFv) against the Target. Split it into
VH1 and VL1. Engineer the VH1/VL1 interface to weaken their spontaneous association
(target: Kd ≈ 10–100 µM in solution).

### Step 3 — Design minibinders against VL1
Use RFdiffusion + ProteinMPNN to design small (40–80 residue) proteins that bind
VL1's VH1-pairing interface (FR2, CDR-L2). Target minibinder Kd ≈ 1–10 µM — tighter
than the weakened VH1-VL1 (so VL1 is locked in free state) but weaker than VH1's
EFFECTIVE affinity when anchored (Kd,VH-VL / C_eff).

### Step 4 — Select or design the nanobody
Choose or design a nanobody against the Target (or a different epitope on the same
target molecule).

### Step 5 — Optimise the (G₄S)ₙ linker
Use `binary_antibodies` to compute the minimum linker that achieves high C_eff at
the expected antigen–target distance.

### Step 6 — Assemble and express
```
[signal peptide] – VH1 – (G4S)3 – Minibinder – (G4S)3 – VL1 – (G4S)n – Nanobody VHH
```

VH1 is a separate polypeptide (or the same chain with a self-cleavage P2A site).

---

## Design Variants

### 1. Split-scFv with minibinder lock *(this repo — recommended)*
Minibinders occlude VL1; VH1 membrane-anchoring displaces them.
Highest OFF-state fidelity; requires minibinder design.

### 2. Simple proximity-gating
Nanobody tethered to anchored antibody via flexible linker only.
Simpler but weaker OFF-state (nanobody can still diffuse near target).

### 3. Steric-occlusion / Probody-like
Nanobody paratope masked by complementary peptide; antigen-binding drives unmasking.
Good OFF-state but harder to engineer structurally.
