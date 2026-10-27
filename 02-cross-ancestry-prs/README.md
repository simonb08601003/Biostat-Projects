# Cross-ancestry polygenic risk scores

How much accuracy does a European-trained PRS lose in East Asians, and how
much does a smaller East Asian GWAS win back? Simulated on real 1000 Genomes
genotypes so the true answer is known.

Report: [`report.html`](report.html)

## Results

- Same causal effects, causal SNPs genotyped: the European score keeps 98% of
  its accuracy in East Asians.
- When effects differ, the ratio follows ρ²: 0.63 at ρ = 0.8, 0.35 at ρ = 0.6.
- If the causal SNPs aren't genotyped, LD differences cost more: 0.85 at ρ = 1.
- A 20,000-person East Asian GWAS adds 34% R² at ρ = 0.8 and 114% at ρ = 0.6,
  and nothing at ρ = 1.
- LDpred2-auto beat C+T in every scenario and recovered h² ≈ 0.40 without tuning.

80 replicates: 3 values of ρ × 2 East Asian GWAS sizes × 10, plus 20 with the
causal SNPs removed.

## Setup

- 503 EUR and 504 EAS from 1000 Genomes phase 3; ~200k SNPs common in both, no MHC.
- h² = 0.4, 0.5% causal SNPs, effects correlated ρ across ancestries.
- GWAS simulated from LD (β̂ = Rβ + ε) in 2 Mb blocks; EUR N = 100k, EAS N = 5k or 20k.
- Scores: C+T, LDpred2-inf, LDpred2-auto, and a EUR + EAS combination tuned in East Asians.
- Accuracy: R² against the true genetic value in a held-out East Asian half.

LDpred2-inf is GBLUP on summary statistics, the same model as my thesis work on
genomic selection.

## Run

```bash
for c in $(seq 1 22); do bash data-raw/get_chr.sh $c; done
bash R/00_prepare_genotypes.sh
Rscript R/01_ld.R
bash R/run_simulations.sh     # ~2.5 h on 4 cores
Rscript R/03_summarise.R
quarto render report.qmd
```

Needs `plink2` and R with `bigsnpr` on the PATH. Genotype and LD files are a
few GB and not committed; per-replicate results are in `outputs/sim/`.
