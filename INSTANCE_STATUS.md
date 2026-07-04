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
