#!/bin/bash
# rho x EAS GWAS size x 10 replicates
cd "$(dirname "$0")/.."
export OPENBLAS_NUM_THREADS=1 PRS_NCORES=2
for rho in 1.0 0.8 0.6; do for n in 5000 20000; do for r in $(seq 1 10); do
  echo "$rho $n $r"; done; done; done |
  xargs -P 2 -L 1 bash -c 'Rscript R/02_simulate_fit.R $0 $1 $2 > /dev/null 2>&1 || echo "failed $0 $1 $2"'

# causal SNPs not genotyped
for rho in 1.0 0.8; do for r in $(seq 1 10); do echo "$rho 20000 $r hidden"; done; done |
  xargs -P 2 -L 1 bash -c 'Rscript R/02_simulate_fit.R $0 $1 $2 $3 > /dev/null 2>&1 || echo "failed $0 $1 $2 $3"'
