#!/bin/bash
# Merge chromosomes and pick the SNP set: MAF >= 5% in both EUR and EAS,
# no MHC, thinned at random to 200k
set -euo pipefail
cd "$(dirname "$0")/.."
N_SNPS=${N_SNPS:-200000}

awk '$3=="EUR"{print $1}' data-raw/panel.txt > data/eur.keep
awk '$3=="EAS"{print $1}' data-raw/panel.txt > data/eas.keep

for c in $(seq 1 22); do echo "data-raw/chr$c"; done > data/merge_list.txt
plink2 --pmerge-list data/merge_list.txt bfile --make-bed --out data/all_raw \
  --threads 4 --memory 8000 > /dev/null

plink2 --bfile data/all_raw --keep data/eur.keep --freq --out data/eur > /dev/null
plink2 --bfile data/all_raw --keep data/eas.keep --freq --out data/eas > /dev/null
awk 'NR==FNR{if(NR>1 && $6>=0.05 && $6<=0.95) e[$2]=1; next} FNR>1 && ($2 in e) && $6>=0.05 && $6<=0.95 {print $2}' \
  data/eur.afreq data/eas.afreq > data/common_both.txt

plink2 --bfile data/all_raw --extract data/common_both.txt \
  --exclude bed1 <(printf '6\t25000000\t34000000\tMHC\n') \
  --thin-count "$N_SNPS" --seed 2026 --make-bed --out data/geno \
  --threads 4 --memory 8000 > /dev/null
echo "SNPs common in both: $(wc -l < data/common_both.txt); kept: $(wc -l < data/geno.bim)"

# drop the few SNPs with almost no variance within one ancestry
Rscript -e 'source("R/00_setup.R"); d <- load_geno();
  a <- pop_stats(d$G, d$ind_eur); b <- pop_stats(d$G, d$ind_eas);
  writeLines(d$map$marker.ID[a$scale < 0.3 | b$scale < 0.3], "data/low_var.txt")'
plink2 --bfile data/geno --exclude data/low_var.txt --make-bed --out data/geno_qc > /dev/null
for ext in bed bim fam; do mv data/geno_qc.$ext data/geno.$ext; done
rm -f data/geno.rds data/geno.bk
echo "after low-variance filter: $(wc -l < data/geno.bim)"
