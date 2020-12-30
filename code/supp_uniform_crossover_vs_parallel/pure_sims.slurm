#!/bin/bash

#SBATCH --partition=multinode
#SBATCH --constraint=x2680
#SBATCH --qos=cv19
#SBATCH --nodes=18
#SBATCH --ntasks=504
#SBATCH --ntasks-per-node=28
#SBATCH --ntasks-per-core=1
#SBATCH --time=0-04:00:00
#SBATCH --exclusive
#SBATCH --mem=247g
#SBATCH --gres=lscratch:10
#SBATCH --array=1-3
#SBATCH --job-name="cr_sim"

module load R

mpiexec Rscript --vanilla pure_sims.R ${SLURM_ARRAY_TASK_ID}
