# Lipids and coronary heart disease: Mendelian randomization

Which blood lipids cause coronary heart disease? Two-sample MR with public
summary statistics only.

Report: [`report.html`](report.html)

## Results

| Lipid (per 1 SD) | Univariable OR (95% CI) | Adjusted for ApoB, HDL, TG |
|---|---|---|
| LDL cholesterol | 1.46 (1.32–1.61) | 1.39 (1.32–1.47), LDL + HDL + TG model |
| Apolipoprotein B | 1.41 (1.28–1.54) | 1.30 (1.24–1.36) |
| Triglycerides | 1.23 (1.16–1.30) | 1.05 (1.00–1.11) |
| HDL cholesterol | 0.83 (0.79–0.87) | 0.89 (0.85–0.93) |

- The LDL / ApoB effect holds up under MR-Egger, weighted median, weighted
  mode, Steiger filtering, and dropping the LPA and 9p21 loci.
- Most of the triglyceride effect disappears once ApoB is in the model.
- The HDL result disagrees with the failed HDL-raising trials; heterogeneity
  is huge (Q p < 1e-90), so it's most likely pleiotropy.
- European instruments only partly carry over to East Asians: 61–71% have the
  same sign in Pan-UKB EAS, with about half the effect size.

## Data

| | Source | |
|---|---|---|
| Exposures | [Pan-UK Biobank](https://pan.ukbb.broadinstitute.org) biomarkers 30780, 30760, 30870, 30640 | EUR, ~400k, GRCh37 |
| Outcome | [FinnGen R12](https://www.finngen.fi/en/access_results) I9_CHD | ~500k, GRCh38, no overlap with UK Biobank |
| LD reference | 1000 Genomes, 503 Europeans | from `../02-cross-ancestry-prs` |

Instruments: p < 5e-8, clumped at r² < 0.001 within 10 Mb, mean F > 130.
Matched across builds by rsID; palindromic SNPs kept only when allele
frequency settles the strand.

## Run

```bash
bash R/00_download_sumstats.sh   # streams ~9 GB, keeps ~1.4 GB
Rscript R/01_instruments.R       # needs plink and tabix on the PATH
Rscript R/02_outcome_harmonise.R
Rscript R/03_mr.R
quarto render report.qmd
```
