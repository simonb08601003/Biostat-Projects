#!/bin/bash
# Stream the summary stats and keep only the columns we use
# Pan-UKB lipids (EUR + EAS, GRCh37); FinnGen R12 I9_CHD (GRCh38)
set -euo pipefail
cd "$(dirname "$0")/../data-raw"
PANUKB=https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files
for t in 30780:ldl 30760:hdl 30870:tg 30640:apob; do
  code=${t%%:*}; name=${t##*:}
  [ -s panukb_$name.tsv.gz ] && continue
  curl -sS --retry 4 "$PANUKB/biomarkers-$code-both_sexes-irnt.tsv.bgz" | gunzip -c |
    awk -F'\t' 'BEGIN{OFS="\t"; print "chr","pos","ref","alt","af_eur","beta_eur","se_eur","mlogp_eur","af_eas","beta_eas","se_eas","mlogp_eas"}
      NR>1 && $25!="NA" && $43=="false" && $19>=0.01 && $19<=0.99 {print $1,$2,$3,$4,$19,$25,$31,$37,$18,$24,$30,$36}' |
    gzip > panukb_$name.tsv.gz
  echo "done $name"
done
[ -s finngen_r12_chd.tsv.gz ] || {
  curl -sS --retry 4 https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_I9_CHD.gz |
    gunzip -c | awk -F'\t' 'BEGIN{OFS="\t"} NR==1 || ($11>=0.005 && $11<=0.995)' | gzip > finngen_r12_chd.tsv.gz
  echo "done chd"
}
