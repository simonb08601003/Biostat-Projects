# Lipid instruments: p < 5e-8 in Pan-UKB EUR, clumped r2 < 0.001 / 10 Mb
# against 1000 Genomes EUR (from project 02), then rsIDs for matching to FinnGen
suppressPackageStartupMessages({ library(data.table); library(dplyr) })
ref_bfile <- "../02-cross-ancestry-prs/data/all_raw"
exposures <- c("ldl", "hdl", "tg", "apob")

ref_ids <- fread(paste0(ref_bfile, ".bim"), header = FALSE, select = 2)[[1]]

clump <- function(name, p1 = 5e-8) {
  cand <- fread(sprintf("data-raw/panukb_%s.tsv.gz", name),
                colClasses = list(character = c("ref", "alt"))) |>
    filter(mlogp_eur >= -log10(p1)) |>
    mutate(
      id = paste(chr, pos, ref, alt, sep = ":"),
      id_flip = paste(chr, pos, alt, ref, sep = ":"),
      SNP = case_when(id %in% ref_ids ~ id, id_flip %in% ref_ids ~ id_flip),
      P = pmax(10^(-mlogp_eur), 1e-300)
    ) |>
    filter(!is.na(SNP))
  f <- tempfile()
  fwrite(select(cand, SNP, P), f, sep = "\t")
  system2("plink", c("--bfile", ref_bfile, "--keep", "data/eur_plink1.keep", "--clump", f,
                     "--clump-p1", p1, "--clump-r2", 0.001, "--clump-kb", 10000,
                     "--out", f, "--silent", "--threads", 4))
  kept <- fread(paste0(f, ".clumped"))$SNP
  message(name, ": ", nrow(cand), " genome-wide SNPs -> ", length(kept), " instruments")
  cand |> filter(SNP %in% kept) |> mutate(exposure = name)
}
instruments <- bind_rows(lapply(exposures, clump))

# rsIDs from the 1000 Genomes VCFs by remote tabix; cached because it is slow
kg <- "https://1000genomes.s3.amazonaws.com/release/20130502/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
rs_cache <- "data/instrument_rsids.tsv"
rs <- if (file.exists(rs_cache)) fread(rs_cache, colClasses = "character") else
  data.table(chr = character(), pos = character(), rsid = character(), ref = character(), alt = character())

todo <- instruments |> filter(!paste(chr, pos) %in% paste(rs$chr, rs$pos))
lookup <- function(d) {
  reg <- tempfile()
  fwrite(unique(d[, c("chr", "pos")]), reg, sep = "\t", col.names = FALSE)
  owd <- setwd(tempdir()); on.exit(setwd(owd))  # tabix drops index files in the cwd
  txt <- system2("tabix", c("-R", reg, sprintf(kg, d$chr[1])), stdout = TRUE)
  if (length(txt)) fread(text = txt, header = FALSE, select = 1:5, colClasses = "character",
                         col.names = c("chr", "pos", "rsid", "ref", "alt"))
}
rs <- rbind(rs, rbindlist(lapply(split(todo, todo$chr), lookup)))
fwrite(rs, rs_cache, sep = "\t")

rs <- rs |>
  mutate(chr = as.integer(chr), pos = as.integer(pos)) |>
  distinct(chr, pos, ref, alt, .keep_all = TRUE)
instruments <- instruments |>
  mutate(chr = suppressWarnings(as.integer(chr)), pos = as.integer(pos)) |>
  filter(!is.na(chr)) |>
  left_join(rs, by = c("chr", "pos", "ref", "alt")) |>
  filter(startsWith(rsid, "rs"))

fwrite(instruments, "data/instruments.tsv", sep = "\t")
print(count(instruments, exposure))
