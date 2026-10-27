# Average the replicates and draw the main figure
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })

res <- bind_rows(lapply(list.files("outputs/sim", full.names = TRUE), readRDS))
hid <- filter(res, !causal_typed)
res <- filter(res, causal_typed)

methods <- c("C+T (EUR)", "LDpred2-inf (EUR)", "LDpred2-auto (EUR)",
             "LDpred2-auto (EAS)", "Multi-ancestry (EUR + EAS)")
mean_se <- function(df, x, ...) {
  df |>
    group_by(...) |>
    summarise(reps = n(), mean = mean({{ x }}), se = sd({{ x }}) / sqrt(n()), .groups = "drop")
}
# EAS / EUR accuracy of the EUR LDpred2-auto score
portability <- function(df, ...) {
  df |>
    filter(method == "LDpred2-auto (EUR)") |>
    select(rho, n_eas, rep, target, r2) |>
    pivot_wider(names_from = target, values_from = r2) |>
    mean_se(EAS / EUR, ...)
}

summ <- mean_se(res, r2, method, target, rho, n_eas) |>
  rename(r2_mean = mean, r2_se = se) |>
  mutate(method = factor(method, levels = methods))
write.csv(summ, "outputs/summary_r2.csv", row.names = FALSE)

port <- portability(res, rho) |> rename(ratio = mean, ratio_se = se)
write.csv(port, "outputs/portability.csv", row.names = FALSE)

gain <- res |>
  filter(target == "EAS", method %in% c("LDpred2-auto (EUR)", "Multi-ancestry (EUR + EAS)")) |>
  select(rho, n_eas, rep, method, r2) |>
  pivot_wider(names_from = method, values_from = r2) |>
  mean_se(`Multi-ancestry (EUR + EAS)` / `LDpred2-auto (EUR)` - 1, rho, n_eas) |>
  rename(rel_gain_mean = mean, rel_gain_se = se)
write.csv(gain, "outputs/multi_ancestry_gain.csv", row.names = FALSE)

pd <- summ |>
  filter(target == "EAS") |>
  mutate(panel = paste0("East Asian GWAS N = ", format(n_eas, big.mark = ",")),
         rho_lab = paste0("ρ = ", format(rho, nsmall = 1)))
p <- ggplot(pd, aes(r2_mean, method, colour = rho_lab)) +
  geom_pointrange(aes(xmin = r2_mean - 2 * r2_se, xmax = r2_mean + 2 * r2_se),
                  position = position_dodge(width = 0.6), size = 0.3) +
  facet_wrap(~ panel) +
  scale_y_discrete(limits = rev(methods)) +
  scale_colour_manual(values = c("ρ = 0.6" = "#b2182b", "ρ = 0.8" = "#ef8a62", "ρ = 1.0" = "#2166ac")) +
  labs(x = "Prediction R² in East Asian test samples (mean ± 2 SE)",
       y = NULL, colour = "EUR–EAS genetic correlation") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
ggsave("outputs/r2_by_method.png", p, width = 10, height = 5.5, dpi = 150)

# causal SNPs untyped vs typed, on the same seeds
paired <- filter(res, n_eas == 20000, rho %in% hid$rho, rep %in% hid$rep)
tag <- bind_rows(paired, hid) |>
  filter(method %in% c("LDpred2-auto (EUR)", "Multi-ancestry (EUR + EAS)")) |>
  mean_se(r2, method, target, rho, causal_typed) |>
  rename(r2_mean = mean, r2_se = se)
write.csv(tag, "outputs/tagging_scenario.csv", row.names = FALSE)
port_tag <- portability(hid, rho) |> rename(ratio = mean, ratio_se = se)
write.csv(port_tag, "outputs/portability_untyped_causal.csv", row.names = FALSE)

print(arrange(summ, target, rho, n_eas, method), n = Inf)
print(port); print(gain); print(tag, n = Inf); print(port_tag)
