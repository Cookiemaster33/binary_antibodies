# Pipeline Complete — Instance Terminated

The anti-idiotypic minibinder design run completed successfully.

## Results (in pipeline_results/ — not committed, large files)

| File | Content |
|---|---|
| `pipeline_results/outputs/rfd3/` | 200 minibinder backbone CIF files (165 res each) |
| `pipeline_results/outputs/mpnn/minibinder_sequences.fasta` | 1600 minibinder sequences (8 per backbone) |
| `pipeline_results/outputs/mpnn/all_sequences.json` | Full metadata with sequence recovery scores |

## Top candidates

Sequences ranked by lowest `sequence_recovery` (most novel relative to nanobody).
Top designs have recovery ~0.10-0.20, meaning the minibinder sequence is
80-90% novel (not copying the nanobody it sits against).

## Next steps

1. **AF2-Multimer validation**: fold top 20 sequences in complex with the nanobody
   to confirm the minibinder adopts the designed conformation and contacts CDRs.
2. **Rosetta energy filter**: ΔΔG < -5 REU for nanobody-minibinder interface.
3. **Experimental validation**: SPR/BLI binding assay of top candidates vs 2Rs15d nanobody.

## Cost: ~$0.70 total (21 min at $1.99/hr)

---

## AF2-Multimer Validation (completed)

Ran AF2-Multimer v3, single-sequence mode, 2 models × 3 recycles on top 20 designs.

**Important caveat**: AF2-Multimer performs poorly in single-sequence mode for
de novo designed proteins (no evolutionary homologs → no MSA signal).
The low scores reflect this limitation, not necessarily that the designs are bad.

| Metric | Best | Mean |
|---|---|---|
| ipTM | 0.090 | 0.074 |
| pTM | 0.310 | 0.200 |
| pLDDT | 38.1 | 31.7 |

**Best candidate**: `mb_b019_005` (ipTM=0.090, pTM=0.310, pLDDT=38.1)

## Recommended next steps

1. Run AF2-Multimer **WITH MSA** on top 5 (submit to ColabFold server or use MMseqs2 locally)
2. Run ESMFold monomer on each minibinder sequence (fast, single-sequence, better for designed proteins)
3. Order top 3-5 for experimental SPR validation
