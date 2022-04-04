#!/bin/sh

#SBATCH --nodes=1
#SBATCH --job-name=w_inputs
#SBATCH --error=input.slurm.log



#SBATCH --partition=short

source activate rmg_env
python-jl write_rmg_inputs.py
