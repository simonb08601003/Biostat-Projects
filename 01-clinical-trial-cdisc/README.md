# CDISC pilot: SDTM → ADaM → TLFs

The CDISC pilot study (xanomeline patch vs placebo in mild-to-moderate
Alzheimer's disease, 254 subjects) taken from raw SDTM to analysis datasets,
QC, and submission-style tables, using the pharmaverse packages.

Report: [`report.html`](report.html)

## Results

| | |
|---|---|
| ADaM | ADSL, ADQSADAS, ADAE, ADTTE built with `admiral` |
| QC vs reference ADaM | 29 variables, 13,400 values, 100% match |
| Primary analysis (ANCOVA, week 24, LOCF) | dose-response p = 0.245; high dose vs placebo −1.0 (−2.7, 0.7), p = 0.233. Same as Table 14-3.01 of the FDA R Submission Pilot |
| MMRM, observed data | high dose vs placebo −0.83 (−2.92, 1.27), p = 0.437 |
| Time to first dermatologic AE | Cox HR 5.28 (3.33, 8.36), high dose vs placebo; PH holds (p = 0.72) |

Getting to 100% meant working out rules the spec doesn't spell out: which
items make up ADAS-Cog(11), how proration works, visit windows, who gets LOCF,
partial-date imputation, and SAS vs R rounding. They're listed in the report.

## Run

```bash
Rscript run_all.R
quarto render report.qmd
```

| File | |
|---|---|
| `R/01–03_*.R` | ADaM derivations |
| `R/07_qc.R` | compare every derived variable to the reference ADaM |
| `R/04_efficacy.R` | ANCOVA and MMRM |
| `R/05_survival.R` | KM, log-rank, Cox, PH check, RMST |
| `R/06_tlf.R` | demographics and TEAE tables |
| `data-raw/` | SDTM `.xpt` |
| `data-reference/` | reference ADaM, used only for QC |

## Data

SDTM and reference ADaM from the CDISC pilot, distributed by PHUSE
(<https://github.com/phuse-org/phuse-scripts>). Published table from the R
Consortium submission pilot (<https://github.com/RConsortium/submissions-pilot1>).
