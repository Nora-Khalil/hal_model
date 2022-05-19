#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=4-00:00:00
#SBATCH --job-name=fc_CH2F2_C2H5F_C2H3F3
#SBATCH --error=fc.slurm_1d1hr.log
#SBATCH --output=fc_output.slurm1d1hr.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west


source activate cantera_env
python flame_speed_calc1d1hr.py
