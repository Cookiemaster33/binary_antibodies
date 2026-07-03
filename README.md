# Binary Antibodies — Conditional Nanobody Design

## Concept

A **conditional proximity-gated bispecific** construct in which a nanobody can **only** bind its target when a linked antibody has already bound a specific antigen.

```
         [Nanobody]
              |
       ~flexible linker~
              |
  [VH2─VL2 / VH1─VL1]          ←── antibody fragment
              |
          [Antigen]  ··· membrane ···  [Target]
```

### Mechanism (AND-gate logic)

| Antibody bound to antigen? | Nanobody near target? | Nanobody binds target? |
|---|---|---|
| No  | No  | **No**  |
| No  | Yes | **No** (linker not anchored, freely diffusing) |
| Yes | Yes | **Yes** (effective local concentration ≫ Kd) |

When the antibody is **not** bound to the antigen, the entire construct diffuses freely in solution — the nanobody never achieves a sufficiently high local concentration near the target. Once the antibody **anchors** to the antigen on the cell surface, the flexible linker constrains the nanobody within a small search volume above the membrane, dramatically increasing its effective local concentration near any surface target in that neighbourhood.

### Design Parameters

1. **Antigen–Target distance** (`d`): centre-to-centre distance on the membrane between the antigen and the target epitope.
2. **Linker length** (`N` residues): must be long enough for the nanobody to reach `d`, but short enough that the effective concentration is therapeutically meaningful.
3. **Nanobody Kd**: intrinsic affinity; the effective conditional Kd is modulated by the effective concentration.
4. **Antibody Kd for antigen**: determines when the construct is "anchored".

### Key Biophysical Insight

The **effective local concentration** of the nanobody near the target, given a flexible linker of `N` residues and an antigen–target distance `d`, is estimated via polymer physics (worm-like chain / freely-jointed chain model):

$$c_\text{eff}(d, N) = \left(\frac{3}{2\pi N l^2}\right)^{3/2} \exp\!\left(-\frac{3d^2}{2Nl^2}\right) \cdot N_A^{-1}$$

where `l ≈ 0.38 nm` is the Cα–Cα virtual bond length and `Nₐ` is Avogadro's number. The conditional apparent Kd becomes:

$$K_d^\text{app} = K_d^\text{nanobody} / c_\text{eff}$$

This framework lets us **optimise the linker length** for any given antigen–target geometry.

---

## Repository Layout

```
binary_antibodies/
├── README.md
├── requirements.txt
├── binary_antibodies/
│   ├── __init__.py
│   ├── polymer.py          # worm-like chain / FJC linker physics
│   ├── design.py           # end-to-end conditional construct design
│   └── sequences.py        # linker sequence generation & composition
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

# Find optimal linker for antigen–target distance of 8 nm,
# nanobody Kd of 10 nM, desired conditional activation fold > 100×
python scripts/optimise_linker.py --distance 8.0 --kd-nanobody 10e-9 --fold 100

# Generate a full parameter-space heatmap
python scripts/scan_geometry.py
```

---

## Design Variants

### 1. Proximity-gated (shown in figure, this repo)
Nanobody is tethered to the antibody C-terminus via a flexible (G4S)ₙ linker.
Activation requires antigen and target to be on the **same cell surface**.

### 2. Steric-occlusion / Probody-like
The nanobody paratope is masked by a complementary peptide that is displaced
upon antibody–antigen engagement. Requires structure-guided mask design.

### 3. Split-nanobody reconstitution
Nanobody is split into two non-functional halves (e.g. at an exposed loop).
One half is on the antibody, the other is soluble. Antigen binding drives
co-localisation and reconstitution.

**This repository focuses on variant 1** (proximity-gated), which is the most
straightforward to engineer and analyse computationally.
