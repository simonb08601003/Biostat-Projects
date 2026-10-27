# Two-sample MR: lipids -> coronary heart disease
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(MendelianRandomization); library(ggplot2)
})
harm <- fread("data/harmonised.tsv", colClasses = list(character = c("ea", "oa")))
labels <- c(ldl = "LDL cholesterol", hdl = "HDL cholesterol", tg = "Triglycerides",
            apob = "Apolipoprotein B")

ivw <- function(d) {
  mr_ivw(mr_input(bx = d$beta_x, bxse = d$se_x, by = d$beta_y, byse = d$se_y), model = "random")
}
est <- function(fit, method, n, lo = "CILower", hi = "CIUpper", p = "Pvalue") {
  tibble(method = method, estimate = fit@Estimate[1], lower = slot(fit, lo)[1],
         upper = slot(fit, hi)[1], p = slot(fit, p)[1], n_snps = n)
}

run_mr <- function(d) {
  inp <- mr_input(bx = d$beta_x, bxse = d$se_x, by = d$beta_y, byse = d$se_y, snps = d$rsid)
  fit_ivw <- mr_ivw(inp, model = "random")
  fit_egg <- mr_egger(inp)
  # Steiger: keep SNPs that explain more variance in the lipid than in CHD
  r2 <- function(b, se, n) (b / se)^2 / ((b / se)^2 + n)
  steiger <- r2(d$beta_x, d$se_x, 4e5) > r2(d$beta_y, d$se_y, 4.5e5)
  # LPA and 9p21 affect CHD outside the lipid pathway
  outlier <- (d$chr == 6 & between(d$pos, 160.5e6, 161.5e6)) |
             (d$chr == 9 & between(d$pos, 21.9e6, 22.2e6))
  n <- nrow(d)
  bind_rows(
    est(fit_ivw, "IVW (random effects)", n),
    est(fit_egg, "MR-Egger", n, "CILower.Est", "CIUpper.Est", "Pvalue.Est"),
    est(mr_median(inp, weighting = "weighted", iterations = 2000, seed = 2026), "Weighted median", n),
    est(mr_mbe(inp, weighting = "weighted", iterations = 2000, seed = 2026), "Weighted mode", n),
    est(ivw(d[steiger]), "IVW, Steiger-filtered", sum(steiger)),
    est(ivw(d[!outlier]), "IVW, without LPA and 9p21", sum(!outlier))
  ) |>
    mutate(egger_intercept = fit_egg@Intercept, egger_intercept_p = fit_egg@Pvalue.Int,
           cochran_q = fit_ivw@Heter.Stat[1], q_p = fit_ivw@Heter.Stat[2],
           mean_F = mean((d$beta_x / d$se_x)^2))
}

add_or <- function(df) mutate(df, OR = exp(estimate), OR_lower = exp(lower), OR_upper = exp(upper))

res <- harm |>
  group_by(exposure) |>
  group_modify(~ run_mr(as.data.table(.x))) |>
  ungroup() |>
  mutate(exposure_label = labels[exposure]) |>
  add_or()
fwrite(res, "outputs/mr_univariable.csv")

# --- Multivariable MR ---
# A: LDL + HDL + TG; B: ApoB + HDL + TG. SNPs = union of all instrument sets,
# with each lipid's effect read from its Pan-UKB file.
snps <- unique(harm[, .(rsid, chr, pos, ea, oa)])
lipid_effects <- function(nm) {
  ex <- fread(sprintf("data-raw/panukb_%s.tsv.gz", nm), colClasses = list(character = c("ref", "alt")))
  ex <- ex[chr %in% as.character(1:22)][, chr := as.integer(chr)]
  m <- merge(snps, ex[, .(chr, pos, ref, alt, beta = beta_eur, se = se_eur)], by = c("chr", "pos"))
  m <- m[(ea == alt & oa == ref) | (ea == ref & oa == alt)]
  m[ea == ref, beta := -beta]
  setnames(m[, .(rsid, beta, se)], c("rsid", paste0("b_", nm), paste0("se_", nm)))
}
mv <- Reduce(\(a, b) merge(a, b, by = "rsid"), lapply(names(labels), lipid_effects))
mv <- merge(mv, unique(harm[, .(rsid, beta_y, se_y)]), by = "rsid")

run_mv <- function(traits, model) {
  inp <- mr_mvinput(bx = as.matrix(mv[, paste0("b_", traits), with = FALSE]),
                    bxse = as.matrix(mv[, paste0("se_", traits), with = FALSE]),
                    by = mv$beta_y, byse = mv$se_y, exposure = labels[traits])
  ivw <- mr_mvivw(inp)
  egg <- mr_mvegger(inp)
  bind_rows(
    tibble(exposure_label = labels[traits], method = "Multivariable IVW",
           estimate = ivw@Estimate, lower = ivw@CILower, upper = ivw@CIUpper, p = ivw@Pvalue),
    tibble(exposure_label = labels[traits], method = "Multivariable MR-Egger",
           estimate = egg@Estimate, lower = egg@CILower.Est, upper = egg@CIUpper.Est, p = egg@Pvalue.Est)
  ) |>
    mutate(model = model)
}
mvres <- bind_rows(run_mv(c("ldl", "hdl", "tg"), "A: LDL + HDL + TG"),
                   run_mv(c("apob", "hdl", "tg"), "B: ApoB + HDL + TG")) |>
  add_or() |>
  mutate(n_snps = nrow(mv))
fwrite(mvres, "outputs/mr_multivariable.csv")

# --- Figures ---
method_order <- c("IVW (random effects)", "MR-Egger", "Weighted median", "Weighted mode",
                  "IVW, Steiger-filtered", "IVW, without LPA and 9p21", "MVMR A", "MVMR B")
fp <- bind_rows(
  res,
  mvres |> filter(method == "Multivariable IVW") |> mutate(method = paste("MVMR", substr(model, 1, 1)))
) |>
  mutate(method = factor(method, levels = rev(method_order)))
p1 <- ggplot(fp, aes(OR, method)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = OR_lower, xmax = OR_upper), width = 0.2, orientation = "y") +
  geom_point(size = 2) +
  facet_wrap(~ exposure_label, ncol = 2) +
  scale_x_log10() +
  labs(x = "Odds ratio for coronary heart disease per 1-SD higher lipid (log scale)", y = NULL) +
  theme_bw()
ggsave("outputs/mr_forest.png", p1, width = 7, height = 7, dpi = 150)

# orient each SNP so it raises the lipid
sc <- harm |> mutate(exposure_label = labels[exposure], s = sign(beta_x),
                     bx = beta_x * s, by = beta_y * s)
lines <- res |>
  filter(method %in% c("IVW (random effects)", "MR-Egger")) |>
  mutate(intercept = ifelse(method == "MR-Egger", egger_intercept, 0))
p2 <- ggplot(sc, aes(bx, by)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_errorbar(aes(ymin = by - 1.96 * se_y, ymax = by + 1.96 * se_y), colour = "grey70", width = 0) +
  geom_errorbar(aes(xmin = bx - 1.96 * se_x, xmax = bx + 1.96 * se_x), colour = "grey70",
                width = 0, orientation = "y") +
  geom_point(size = 1) +
  geom_abline(data = lines, aes(slope = estimate, intercept = intercept, colour = method)) +
  facet_wrap(~ exposure_label, scales = "free", ncol = 2) +
  labs(x = "SNP effect on lipid (SD)", y = "SNP effect on CHD (log-odds)", colour = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("outputs/mr_scatter.png", p2, width = 9, height = 7, dpi = 150)

# do the European instruments carry over to East Asians?
port <- harm |>
  filter(!is.na(beta_x_eas)) |>
  group_by(exposure) |>
  summarise(n = n(),
            slope_eas_on_eur = coef(lm(beta_x_eas ~ 0 + beta_x, weights = 1 / se_x_eas^2))[1],
            sign_concordance = mean(sign(beta_x_eas) == sign(beta_x)),
            r = cor(beta_x, beta_x_eas))
fwrite(port, "outputs/instrument_portability_eas.csv")

options(width = 160)
print(select(res, exposure, method, OR, OR_lower, OR_upper, p, n_snps, egger_intercept_p, q_p, mean_F))
print(select(mvres, model, exposure_label, method, OR, OR_lower, OR_upper, p, n_snps))
print(port)
