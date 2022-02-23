#!/bin/bash
#BSUB -n 22
#BSUB -W 10:00
#BSUB -R rusage[mem=2500]
#BSUB -R span[hosts=1]​
module load R/3.6.0
module load gcc/8.1.0
module load binutils/2.35

Rscript --vanilla mix_trainRun.R 17 2 22 1
