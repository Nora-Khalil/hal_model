#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --job-name=fc_C2H3F3_CH3F_1day
#SBATCH --error=fc.slurm1d.log
#SBATCH --output=fc_output_1day.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west


source activate cantera_env
python flame_speed_calc_1d.py
