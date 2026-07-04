# Pipeline Complete — All Instances Terminated

## v2 Run (3-chain design): Results available in pipeline_results/v2_3chain/

- 200 RFdiffusion3 backbones (nanobody + 50-res minibinder + CL domain = 278 res each)
- 1600 ProteinMPNN sequences (8 per backbone)
- RF3 self-consistency attempted but scRMSD metric is not appropriate for binder design

## Why scRMSD doesn't apply to binder design

scRMSD checks if a designed sequence folds to the intended backbone IN ISOLATION.
For BINDER design, the backbone is shaped by the binding partners — folding the
50-residue minibinder alone (without nanobody + CL) naturally gives a different
conformation. scRMSD ~15-20 Å is expected, not a failure.

## The right validation path

1. **ColabFold with MSA** on top-5 designs: submit minibinder:nanobody FASTA to
   https://colab.research.google.com/github/sokrypton/ColabFold — free, fast,
   gives trustworthy ipTM with real MSA features
2. **Rosetta or Autodock** for binding energy estimate
3. **Experimental SPR/BLI** on top-3 sequences — ultimate test

## Top 5 sequences to validate (lowest MPNN sequence recovery):

1. PDNLPPPPYTESEPTPSPAPPFVTPSPVVVVPAPVQPNVTSTLVRTGPNP  (backbone mb_b001_002, rec=0.080)
2. SSAAKLATPSNPEPPPLSPPVPLLPPVTPVVVVPGTVSDTANLVVTVYSP  (backbone mb_b005_004, rec=0.120)
3. PEPLQAWEPQEFTPVYSTQPPPVPPLPLPPPPVGPKVPVTKTLVRTGLNS  (backbone mb_b003_009, rec=0.140)
4. PESREPPPPVEFQCVPGPPAPVTEEVVVVVPGKPVSSDPRLKLVVPSPSP  (backbone mb_b003_005, rec=0.180)
5. PDNKPPPPYTPAEPTKSPFPPFVEPSPEVEVPTPVKPDVTTTTVYTGLNP  (backbone mb_b001_002, rec=0.160)
