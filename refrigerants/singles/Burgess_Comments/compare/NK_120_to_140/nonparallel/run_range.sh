#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=4-00:00:00
#SBATCH --job-name=fc_nonpar
#SBATCH --error=fc_nonpar.error.slurm.log
#SBATCH --output=fc_nonpar.output.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west


source activate cantera_env
python flame_speed_calc_range_nonparallel.py
