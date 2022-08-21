#! /usr/bin/env bash

#BSUB -J knit
#BSUB -o logs/knit_%J.out
#BSUB -e logs/knit_%J.err
#BSUB -q rna
#BSUB -R "rusage[mem=16] span[hosts=1]"

module load R/4.0.3

Rscript 'knit_rmd.R' \
    -i 'analysis.Rmd' \
    -o 'analysis.html'


