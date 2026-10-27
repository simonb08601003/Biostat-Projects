# Biostatistics portfolio

Three projects built on public data: a clinical trial analysis, a polygenic
risk score study, and a Mendelian randomization study. Each folder has its own
README, a rendered `report.html`, and the scripts to rebuild it from scratch.

| Project | Question | Main result |
|---|---|---|
| [Clinical trial (CDISC)](01-clinical-trial-cdisc) | Does xanomeline help cognition in Alzheimer's disease, and is it safe? | ADaM matches the sponsor's reference 100%; the primary table reproduces the FDA pilot exactly |
| [Cross-ancestry PRS](02-cross-ancestry-prs) | How well does a European PRS work in East Asians? | Accuracy loss follows ρ²; a 20k East Asian GWAS adds 34–114% R² |
| [Mendelian randomization](03-mendelian-randomization) | Which blood lipids cause heart disease? | LDL-C / ApoB do (OR 1.46 per SD); the triglyceride effect mostly goes away after adjusting for ApoB |

R 4.4, PLINK 1.9/2, Quarto.
