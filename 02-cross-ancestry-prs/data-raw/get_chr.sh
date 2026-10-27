#!/bin/bash
# Download one 1000 Genomes chromosome, keep EUR + EAS biallelic SNPs
set -e
c=$1
cd "$(dirname "$0")"
[ -f chr$c.bed ] && exit 0
url=https://1000genomes.s3.amazonaws.com/release/20130502/ALL.chr$c.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
[ -f chr$c.vcf.gz ] || curl -sS --retry 4 -o chr$c.vcf.gz $url
plink2 --vcf chr$c.vcf.gz --keep keep_eur_eas.txt --snps-only just-acgt --max-alleles 2 \
  --maf 0.01 --rm-dup exclude-all --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1 missing \
  --make-bed --out chr$c --threads 1 --memory 3000 > chr$c.convert.log 2>&1
rm -f chr$c.vcf.gz
