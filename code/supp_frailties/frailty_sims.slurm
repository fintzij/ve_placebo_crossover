#!/bin/bash

#SBATCH --partition=multinode
#SBATCH --constraint=x2695
#SBATCH --qos=cv19
#SBATCH --nodes=3
#SBATCH --ntasks=84
#SBATCH --ntasks-per-node=28
#SBATCH --ntasks-per-core=1
#SBATCH --time=0-04:00:00
#SBATCH --exclusive
#SBATCH --mem=247g
#SBATCH --gres=lscratch:10
#SBATCH --array=1-16
#SBATCH --job-name="cr_sim"

module load R

mpiexec Rscript --vanilla frailty_sims.R ${SLURM_ARRAY_TASK_ID}
