#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=7-00:00:00
#SBATCH --job-name=fc_CH3CF3(1)_C2H5F(2)
#SBATCH --error=slurm/error/fc_error.slurm.log
#SBATCH --output=slurm/output/fc_output.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west

source activate cantera_env
python flame_speed_calc.py

