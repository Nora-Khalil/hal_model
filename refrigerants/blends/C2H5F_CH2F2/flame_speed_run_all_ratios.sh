#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=FIXED_fcall_C2H5F_CH2F2
#SBATCH --error=FIXED_fcall.slurm.log
#SBATCH --output=FIXED_fc_outputall.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=short


source activate cantera_env
python flame_speed_calc_all_ratios.py


