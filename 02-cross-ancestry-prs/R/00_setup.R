suppressPackageStartupMessages({
  library(bigsnpr)
  library(bigstatsr)
  library(dplyr)
})
bigparallelr::set_blas_ncores(1)
NCORES <- as.integer(Sys.getenv("PRS_NCORES", min(4, nb_cores())))

load_geno <- function() {
  if (!file.exists("data/geno.rds")) snp_readBed("data/geno.bed")
  obj <- snp_attach("data/geno.rds")
  panel <- read.table("data-raw/panel.txt", header = TRUE, fill = TRUE)
  pop <- panel$super_pop[match(obj$fam$sample.ID, panel$sample)]
  list(G = obj$genotypes, map = obj$map, pop = pop,
       ind_eur = which(pop == "EUR"), ind_eas = which(pop == "EAS"))
}

# per-ancestry mean and SD of each SNP
pop_stats <- function(G, ind) {
  st <- big_colstats(G, ind.row = ind, ncores = NCORES)
  data.frame(center = st$sum / length(ind), scale = sqrt(st$var))
}

# block-diagonal LD (2 Mb blocks) for SNPs in `cols`
build_ld <- function(G, ind, cols, map, backingfile = tempfile(tmpdir = "data", pattern = "corr_")) {
  block <- paste(map$chromosome, map$physical.pos %/% 2e6)
  ld_all <- NULL
  for (chr in 1:22) {
    cols_chr <- cols[map$chromosome[cols] == chr]
    blk <- factor(block[cols_chr], levels = unique(block[cols_chr]))
    mats <- lapply(split(cols_chr, blk), function(cc) {
      snp_cor(G, ind.row = ind, ind.col = cc, size = length(cc))
    })
    corr0 <- as(Matrix::bdiag(mats), "generalMatrix")
    ld_all <- c(ld_all, Matrix::colSums(corr0^2))
    if (chr == 1) corr <- as_SFBM(corr0, backingfile, compact = TRUE)
    else corr$add_columns(corr0, nrow(corr))
  }
  list(corr = corr, ld = ld_all)
}
