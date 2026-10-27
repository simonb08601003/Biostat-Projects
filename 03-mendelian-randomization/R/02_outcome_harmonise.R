# Match instruments to FinnGen CHD by rsID and align alleles
suppressPackageStartupMessages({ library(data.table); library(dplyr) })

ins <- fread("data/instruments.tsv", colClasses = list(character = c("ref", "alt")))
chd <- fread("data-raw/finngen_r12_chd.tsv.gz", colClasses = list(character = c("ref", "alt")))
setnames(chd, "#chrom", "chrom")
# some rows list several rsIDs
chd <- chd[, .(rsid = unlist(strsplit(rsids, ","))), by = .(chrom, pos38 = pos, ref, alt, beta, sebeta, pval, af_alt)]
chd <- chd[rsid %in% ins$rsid]

palindromic <- function(a, b) paste0(a, b) %in% c("AT", "TA", "CG", "GC")
dat <- ins |>
  select(exposure, rsid, chr, pos, ea = alt, oa = ref, eaf = af_eur,
         beta_x = beta_eur, se_x = se_eur, mlogp_x = mlogp_eur,
         beta_x_eas = beta_eas, se_x_eas = se_eas) |>
  inner_join(chd |> select(rsid, o_alt = alt, o_ref = ref, beta_y = beta,
                           se_y = sebeta, p_y = pval, eaf_y = af_alt),
             by = "rsid", relationship = "many-to-many") |>
  mutate(
    same = ea == o_alt & oa == o_ref,
    swapped = ea == o_ref & oa == o_alt,
    beta_y = case_when(same ~ beta_y, swapped ~ -beta_y),
    eaf_y = case_when(same ~ eaf_y, swapped ~ 1 - eaf_y),
    pal = palindromic(ea, oa),
    # palindromic SNPs only if allele frequency settles the strand
    drop = is.na(beta_y) | (pal & (abs(eaf - 0.5) < 0.08 | abs(eaf - eaf_y) > 0.2))
  )

harm <- dat |> filter(!drop) |> distinct(exposure, rsid, .keep_all = TRUE)
fwrite(harm, "data/harmonised.tsv", sep = "\t")
dat |> group_by(exposure) |>
  summarise(instruments = n(), dropped_palindromic_or_mismatch = sum(drop), used = sum(!drop)) |>
  print()
