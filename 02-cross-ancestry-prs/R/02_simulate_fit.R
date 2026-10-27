# One simulation replicate: Rscript R/02_simulate_fit.R <rho> <n_eas> <rep> [hidden]
#
# Effects are shared by EUR and EAS with genetic correlation rho. GWAS summary
# stats are simulated from each ancestry's LD (beta_hat = R beta + e), PRS are
# built from them, and accuracy is R2 against the true genetic value in a
# held-out half of the EAS samples. "hidden" drops the causal SNPs first, so
# scores have to rely on tagging SNPs whose LD differs between ancestries.
source("R/00_setup.R")
args <- commandArgs(trailingOnly = TRUE)
rho <- as.numeric(args[1])
n_eas <- as.numeric(args[2])
rep <- as.integer(args[3])
hidden <- isTRUE(args[4] == "hidden")
n_eur <- 100000
h2 <- 0.4
p_causal <- 0.005

out_file <- sprintf("outputs/sim/rho%.1f_neas%d_rep%02d%s.rds", rho, n_eas, rep,
                    if (hidden) "_hidden" else "")
if (file.exists(out_file)) quit(save = "no")
dir.create("outputs/sim", showWarnings = FALSE)
set.seed(1000 * rep + 100 * rho + n_eas / 1000)

d <- load_geno()
G <- d$G
map <- d$map
M <- ncol(G)
ld_eur <- readRDS("data/ld_eur.rds")
ld_eas <- readRDS("data/ld_eas.rds")
sd_eur <- ld_eur$stats$scale
sd_eas <- ld_eas$stats$scale

# per-allele effects, scaled so h2 = 0.4 in EUR
causal <- sort(sample(M, round(p_causal * M)))
z1 <- rnorm(length(causal))
z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(length(causal))
b_eur <- b_eas <- numeric(M)
b_eur[causal] <- z1
b_eas[causal] <- z2
k <- sqrt(h2 / sum((b_eur * sd_eur)^2))
b_eur <- b_eur * k
b_eas <- b_eas * k

# Summary stats, one 2 Mb block at a time. Using a genome-wide R from ~500
# people would add spurious long-range LD.
blocks <- split(seq_len(M), paste(map$chromosome, map$physical.pos %/% 2e6))
sim_sumstats <- function(b, ind, st, N) {
  n <- length(ind)
  beta_hat <- numeric(M)
  for (cols in blocks) {
    xb <- if (any(b[cols] != 0)) big_prodVec(G, b[cols], ind.row = ind, ind.col = cols) else rep(0, n)
    v <- (xb - mean(xb)) / (n - 1) + rnorm(n) / sqrt((n - 1) * N)
    beta_hat[cols] <- big_cprodVec(G, v, ind.row = ind, ind.col = cols,
                                   center = st$center[cols], scale = st$scale[cols])
  }
  data.frame(beta = beta_hat, beta_se = 1 / sqrt(N), n_eff = N)
}
ss_eur <- sim_sumstats(b_eur, d$ind_eur, ld_eur$stats, n_eur)
ss_eas <- sim_sumstats(b_eas, d$ind_eas, ld_eas$stats, n_eas)

typed <- seq_len(M)
if (hidden) {
  typed <- setdiff(typed, causal)
  tmp_ld <- c(tempfile(), tempfile())
  ld_eur <- c(build_ld(G, d$ind_eur, typed, map, tmp_ld[1]), list(stats = ld_eur$stats))
  ld_eas <- c(build_ld(G, d$ind_eas, typed, map, tmp_ld[2]), list(stats = ld_eas$stats))
  ss_eur <- ss_eur[typed, ]
  ss_eas <- ss_eas[typed, ]
}
expand <- function(w) replace(numeric(M), typed, w)

# --- PRS weights (per standardized genotype) ---
ldsc_h2 <- function(ss, ld) {
  snp_ldsc(ld$ld, length(ld$ld), chi2 = (ss$beta / ss$beta_se)^2,
           sample_size = ss$n_eff, blocks = NULL)[["h2"]]
}
h2_eur <- ldsc_h2(ss_eur, ld_eur)
h2_eas <- ldsc_h2(ss_eas, ld_eas)

ldpred_auto <- function(ss, ld, h2_init) {
  chains <- snp_ldpred2_auto(
    ld$corr, ss, h2_init = max(h2_init, 0.01),
    vec_p_init = seq_log(1e-4, 0.2, length.out = 8),
    allow_jump_sign = FALSE, shrink_corr = 0.95,
    burn_in = 300, num_iter = 200, ncores = NCORES
  )
  # drop chains that did not converge (standard LDpred2 filter)
  spread <- sapply(chains, \(ch) diff(range(ch$corr_est)))
  ok <- chains[spread > 0.95 * quantile(spread, 0.95, na.rm = TRUE)]
  list(beta = expand(rowMeans(sapply(ok, \(ch) ch$beta_est))),
       h2 = mean(sapply(ok, \(ch) ch$h2_est)),
       p = mean(sapply(ok, \(ch) ch$p_est)))
}
auto_eur <- ldpred_auto(ss_eur, ld_eur, h2_eur)
auto_eas <- ldpred_auto(ss_eas, ld_eas, h2_eas)
inf_eur <- expand(snp_ldpred2_inf(ld_eur$corr, ss_eur, h2 = max(h2_eur, 0.01)))

# C+T: clump on EUR LD (r2 < 0.1, 250 kb), then pick a p-value threshold
z_eur <- expand(ss_eur$beta / ss_eur$beta_se)
clumped <- snp_clumping(G, infos.chr = map$chromosome, infos.pos = map$physical.pos,
                        ind.row = d$ind_eur, exclude = setdiff(seq_len(M), typed),
                        S = abs(z_eur), thr.r2 = 0.1, size = 250, ncores = NCORES)
lp <- -pnorm(-abs(z_eur), log.p = TRUE) / log(10) - log10(2)
thr <- c(-log10(5e-8), 6, 5, 4, 3, 2, 1, 0.5, 0)
ct_weights <- function(t) {
  keep <- intersect(clumped, which(lp >= t))
  replace(numeric(M), keep, expand(ss_eur$beta)[keep])
}

# --- tune on half of EAS, test on the other half ---
half <- sample(rep(1:2, length.out = length(d$ind_eas)))
ind_tune <- d$ind_eas[half == 1]
ind_test <- d$ind_eas[half == 2]

prs <- function(w, sd, ind) {
  w <- w / sd
  w[!is.finite(w)] <- 0
  big_prodVec(G, w, ind.row = ind, ncores = NCORES)
}
g_val <- function(b, ind) { g <- big_prodVec(G, b, ind.row = ind); g - mean(g) }
r2 <- function(score, g, var_g) cor(score, g)^2 * var_g

g_tune <- g_val(b_eas, ind_tune)
g_test <- g_val(b_eas, ind_test)
g_eur <- g_val(b_eur, d$ind_eur)
var_g_eas <- sum((b_eas * sd_eas)^2)
y_tune <- g_tune + rnorm(length(g_tune), sd = sqrt(1 - var_g_eas))

ct_r2 <- sapply(thr, \(t) {
  w <- ct_weights(t)
  if (all(w == 0)) NA else cor(prs(w, sd_eur, ind_tune), y_tune)^2
})
best_thr <- thr[which.max(ct_r2)]
w_ct <- ct_weights(best_thr)

# multi-ancestry: regress EAS phenotype on the EUR and EAS scores
s_eur <- prs(auto_eur$beta, sd_eur, ind_tune)
s_eas <- prs(auto_eas$beta, sd_eas, ind_tune)
mix <- coef(lm(y_tune ~ s_eur + s_eas))[-1]
mix[is.na(mix)] <- 0

eur_scores <- list("C+T (EUR)" = w_ct, "LDpred2-inf (EUR)" = inf_eur,
                   "LDpred2-auto (EUR)" = auto_eur$beta)
test <- c(
  lapply(eur_scores, prs, sd = sd_eur, ind = ind_test),
  list("LDpred2-auto (EAS)" = prs(auto_eas$beta, sd_eas, ind_test))
)
test[["Multi-ancestry (EUR + EAS)"]] <- mix[1] * test[["LDpred2-auto (EUR)"]] +
  mix[2] * test[["LDpred2-auto (EAS)"]]

res <- bind_rows(
  data.frame(method = names(test), target = "EAS",
             r2 = sapply(test, r2, g = g_test, var_g = var_g_eas)),
  data.frame(method = names(eur_scores), target = "EUR",
             r2 = sapply(eur_scores, \(w) r2(prs(w, sd_eur, d$ind_eur), g_eur, h2)))
) |>
  mutate(rho = rho, n_eas = n_eas, rep = rep, causal_typed = !hidden, h2_eas_true = var_g_eas,
         h2_eur_ldsc = h2_eur, h2_eur_auto = auto_eur$h2, p_eur_auto = auto_eur$p,
         ct_threshold = best_thr, w_eur = unname(mix[1]), w_eas = unname(mix[2]))
rownames(res) <- NULL
saveRDS(res, out_file)
if (hidden) unlink(paste0(tmp_ld, ".sbk"))
print(res[, c("method", "target", "r2")])
