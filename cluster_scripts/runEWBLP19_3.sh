#!/bin/bash
#BSUB -n 22
#BSUB -W 10:00
#BSUB -R rusage[mem=2500]
#BSUB -R span[hosts=1]​
module load R/3.6.0
module load gcc/8.1.0
module load binutils/2.35

Rscript --vanilla EW_BLP_BMC_run.R 19 3 22 1
