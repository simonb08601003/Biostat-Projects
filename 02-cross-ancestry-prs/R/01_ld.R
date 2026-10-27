# LD reference per ancestry, same 2 Mb blocks as the simulation
source("R/00_setup.R")
d <- load_geno()
for (anc in c("EUR", "EAS")) {
  ind <- if (anc == "EUR") d$ind_eur else d$ind_eas
  ld <- build_ld(d$G, ind, seq_len(ncol(d$G)), d$map,
                 backingfile = file.path("data", paste0("corr_", anc)))
  ld$stats <- pop_stats(d$G, ind)
  saveRDS(ld, file.path("data", paste0("ld_", tolower(anc), ".rds")))
  message(anc, ": LD for ", length(ld$ld), " SNPs, mean LD score ", round(mean(ld$ld), 2))
}
