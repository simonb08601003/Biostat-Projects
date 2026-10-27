# Time to first dermatologic AE: KM, log-rank, Cox, PH check, RMST
source("R/00_setup.R")
suppressPackageStartupMessages({ library(survival); library(ggsurvfit); library(ggplot2) })

adtte <- readRDS("adam/adtte.rds") |>
  filter(SAFFL == "Y", PARAMCD == "TTDE") |>
  mutate(TRTA = factor(TRTA, levels = c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")),
         EVENT = 1 - CNSR)

km <- survfit(Surv(AVAL, EVENT) ~ TRTA, data = adtte)
lr <- survdiff(Surv(AVAL, EVENT) ~ TRTA, data = adtte)
cox <- coxph(Surv(AVAL, EVENT) ~ TRTA + AGE + SEX, data = adtte)

p <- ggsurvfit(survfit2(Surv(AVAL, EVENT) ~ TRTA, data = adtte), linewidth = 0.8) +
  add_censor_mark(size = 1.5) +
  add_confidence_interval(alpha = 0.1) +
  add_risktable(risktable_stats = "n.risk") +
  scale_ggsurvfit(x_scales = list(breaks = seq(0, 210, 30))) +
  labs(x = "Days since first dose", y = "Proportion without a dermatologic event",
       title = "Time to first dermatologic adverse event (safety population)") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))
ggsave("outputs/km_ttde.png", p, width = 8, height = 6, dpi = 150)

res <- list(
  n = count(adtte, TRTA, EVENT),
  median = summary(km)$table[, c("records", "events", "median", "0.95LCL", "0.95UCL")],
  logrank_p = pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE),
  cox = summary(cox)$conf.int,
  cox_p = summary(cox)$coefficients[, "Pr(>|z|)"],
  zph = cox.zph(cox)$table,
  rmst180 = summary(km, rmean = 180)$table[, c("rmean", "se(rmean)")]
)
saveRDS(res, "outputs/survival_results.rds")

options(width = 120)
capture.output({
  print(res$median)
  cat("\nLog-rank p =", signif(res$logrank_p, 3), "\n\nCox (HR, 95% CI)\n")
  print(cbind(res$cox[, c(1, 3, 4)], p = res$cox_p))
  cat("\nPH test (cox.zph)\n"); print(res$zph)
  cat("\nRMST to day 180\n"); print(res$rmst180)
}, file = "outputs/survival_results.txt")
cat(readLines("outputs/survival_results.txt"), sep = "\n")
