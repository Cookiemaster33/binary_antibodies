# Pipeline Complete — All Instances Terminated

## v3 Bispecific Bridging Minibinder Results

Pipeline completed successfully. Instance terminated.

**Results in `pipeline_results/v3_bispecific/mpnn/`**

| File | Content |
|---|---|
| `minibinder_sequences.fasta` | Top 1600 sequences ranked by novelty |
| `all_sequences.json` | Full metadata (backbone, sequence, recovery score) |

### Stats
- 200 RFdiffusion3 backbones  
- 1600 ProteinMPNN sequences (8 per backbone)
- Sequence length: 56 residues (bispecific bridging)
- Sequence recovery range: 0.164–0.636 (lower = more novel)
- **Best design: 16.4% recovery** = 84% novel sequence

### Top 5 candidates (lowest sequence recovery)

```
rank_1  SDGSTGPPLSNCDPTNRGTTLNSNGNGVNGGLANSSNAGNTCYCENGVCMNETTSQ
rank_2  ADGTSPGTFCDRYTPGQPAVDTSLSPANYDASKLLQVQPFIDNYSKELLTQQDSSQ
rank_3  ADGTSPGTFCDCYVPGKPATDTKLERANYDPAKLLSVQPFKCLRSGEILTQEDSSQ
rank_4  ADGSTGPPLSNCDPTNRSTTLDANGNGVNGGLATASNAGNTGYCVNGVCDSETTSQ
rank_5  SDGSTGPPLSNCDPTNRSATLNSNGQGVNGGLATASNAKNSCLCSNGVCLNETTSQ
```

### Next steps

1. **ColabFold validation** — submit `minibinder:VH1:nanobody` as a 3-chain complex
   to https://colab.research.google.com/github/sokrypton/ColabFold
2. **Rosetta ΔΔG** — score binding to VH1-CH1-face and to nanobody CDRs
3. **Order top 3-5** as synthetic peptides for SPR/BLI assay

## Cost

~$1.50 total (38 min on A100 @ $1.99/hr)
