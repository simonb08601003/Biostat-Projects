# ADAS-Cog(11) change from baseline
# primary: ANCOVA at week 24 with LOCF (as in the SAP); sensitivity: MMRM on observed data
source("R/00_setup.R")
suppressPackageStartupMessages({ library(mmrm); library(emmeans) })

arms <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
adqs <- readRDS("adam/adqsadas.rds") |>
  filter(EFFFL == "Y", ANL01FL == "Y", AVISITN > 0) |>
  mutate(
    TRTP = factor(TRTP, levels = arms),
    AVISIT = factor(AVISIT, levels = c("Week 8", "Week 16", "Week 24")),
    SITEGR1 = factor(SITEGR1),
    USUBJID = factor(USUBJID)
  )

w24 <- filter(adqs, AVISIT == "Week 24")
fit_dose <- lm(CHG ~ TRTPN + SITEGR1 + BASE, data = w24)
fit_anc <- lm(CHG ~ TRTP + SITEGR1 + BASE, data = w24)
anc_emm <- emmeans(fit_anc, ~ TRTP)

obs <- filter(adqs, is.na(DTYPE))
fit_mmrm <- mmrm(CHG ~ TRTP * AVISIT + BASE * AVISIT + SITEGR1 + us(AVISIT | USUBJID),
                 data = obs, method = "Kenward-Roger")
mmrm_emm <- emmeans(fit_mmrm, ~ TRTP | AVISIT)

# unadjusted comparisons vs placebo, as in the study report
vs_pbo <- function(emm, label) {
  contrast(emm, method = "trt.vs.ctrl", ref = 1, adjust = "none") |>
    summary(infer = TRUE) |>
    as.data.frame() |>
    mutate(analysis = label)
}

res <- list(
  n_w24 = count(w24, TRTP),
  n_obs = count(obs, TRTP, AVISIT),
  dose_response = summary(fit_dose)$coefficients["TRTPN", , drop = FALSE],
  ancova_lsmeans = as.data.frame(summary(anc_emm)),
  ancova_contrasts = vs_pbo(anc_emm, "ANCOVA, Week 24 LOCF"),
  mmrm_lsmeans = as.data.frame(summary(mmrm_emm)),
  mmrm_contrasts = vs_pbo(mmrm_emm, "MMRM, observed"),
  mmrm_cov = VarCorr(fit_mmrm)
)
saveRDS(res, "outputs/efficacy_results.rds")

options(width = 120)
capture.output({
  cat("Dose-response (ANCOVA, week 24 LOCF, dose as continuous)\n"); print(res$dose_response)
  cat("\nANCOVA LS means\n"); print(res$ancova_lsmeans)
  cat("\nANCOVA vs placebo\n"); print(res$ancova_contrasts)
  cat("\nMMRM LS means (observed data, unstructured, Kenward-Roger)\n"); print(res$mmrm_lsmeans)
  cat("\nMMRM vs placebo\n"); print(res$mmrm_contrasts)
}, file = "outputs/efficacy_results.txt")
cat(readLines("outputs/efficacy_results.txt"), sep = "\n")
