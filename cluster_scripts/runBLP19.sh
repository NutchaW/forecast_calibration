#!/bin/bash
#BSUB -n 45
#BSUB -W 11:00
#BSUB -R rusage[mem=2500]
#BSUB -R span[hosts=1]​
module load R/3.6.0
module load gcc/8.1.0
module load binutils/2.35

Rscript --vanilla BLP_BMC_run.R 19 1 1
